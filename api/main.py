from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pathlib import Path
import sys
import traceback
import json

from api.services import structure_parser
from api.services.structure_parser import TubulinStructureParser

# --- Path Setup ---
# Add parent directory to path to allow imports from sibling packages
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

# --- Imports ---
# New strict schemas and mapper
from api.schemas import AlignmentRequest, AlignmentResponse
from api.services.alignment import TubulinAlignmentMapper

# Existing logic
from lib.tubulin_analyzer import SpatialGridGenerator, GridData, DebugData

# --- App Configuration ---
app = FastAPI(
    title="Tubulin Spatial Grid API",
    version="0.2.0", # Bumped version for refactor
    description="Generates idealized, aligned 2D grids and MSAs from 3D tubulin structures.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Service Initialization ---
grid_generator = SpatialGridGenerator()

# Paths
MUSCLE_BINARY = str(parent_dir / "muscle3.8.1")
MASTER_PROFILE = str(parent_dir / "data" / "alpha_tubulin" / "alpha_tubulin.afasta")

# Initialize Mapper (Fail fast if binary/profile missing)
try:
    alignment_mapper = TubulinAlignmentMapper(
        master_profile_path=MASTER_PROFILE, 
        muscle_binary=MUSCLE_BINARY
    )
except Exception as e:
    print(f"CRITICAL: Failed to initialize AlignmentMapper: {e}")
    # We don't exit here to allow other endpoints (like models) to work, 
    # but alignment endpoints will fail.
    alignment_mapper = None

# --- MSA / Alignment Endpoints ---

@app.post("/msaprofile/sequence", response_model=AlignmentResponse)
async def align_sequence(request: AlignmentRequest):
    """
    Align a tubulin sequence to the master alpha-tubulin profile.
    
    CRITICAL: For PDB structures, the frontend MUST provide 'auth_seq_ids' 
    that correspond 1:1 with the 'sequence'. This ensures the returned mapping
    uses strict PDB numbering rather than inferred indices.
    """
    if not alignment_mapper:
        raise HTTPException(status_code=503, detail="Alignment service unavailable (configuration error)")

    try:
        # 1. Logic: Align
        # Pass the strict auth_seq_ids to the mapper
        aligned_sequence, mapping = alignment_mapper.align_sequence(
            request.sequence_id, 
            request.sequence,
            request.auth_seq_ids
        )
        structure_parser = TubulinStructureParser()
    # 2. On-the-fly Verification
        # Only verify if we have a valid PDB ID structure (e.g., "1JFF_A")
        if request.sequence_id and "_" in request.sequence_id:
            pdb_id, chain_id = request.sequence_id.split("_")[:2]
            
            # Run verification (Fetch -> Parse -> Compare)
            print(f"Verifying {pdb_id} chain {chain_id} against RCSB...")
            verification = structure_parser.verify_integrity(
                pdb_id, chain_id, request.sequence, request.auth_seq_ids
            )
            
            if verification["status"] == "success":
                if verification["match"]:
                    print(f"✅ VERIFIED: Frontend data for {pdb_id}_{chain_id} matches PDB perfectly.")
                else:
                    print(f"⚠️ MISMATCH: Frontend vs PDB divergence for {pdb_id}_{chain_id}")
                    print(verification["details"])
                    # Optional: Raise HTTPException if strictness is required
            else:
                print(f"⚠️ Verification skipped: {verification.get('reason')}")
                    
        # 2. Logic: Map Annotations (if any)
        mapped_annotations = []
        if request.annotations:
            # The mapper expects Annotation objects, which match our Pydantic schema
            mapped_annotations = alignment_mapper.map_annotations(request.annotations, mapping)
        
        # 3. Logic: Statistics
        statistics = alignment_mapper.get_alignment_statistics(aligned_sequence, mapping)
        
        # 4. Response (Pydantic will validate this against AlignmentResponse)
        return {
            "sequence_id": request.sequence_id,
            "aligned_sequence": aligned_sequence,
            "mapping": mapping,
            "mapped_annotations": mapped_annotations,
            "statistics": statistics,
            "original_sequence": request.sequence
        }
        
    except ValueError as ve:
        # Handles validation errors (e.g. length mismatch)
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Alignment failed: {str(e)}")


@app.get("/msaprofile/master")
async def get_master_profile():
    """Get information about the master alignment profile including raw sequences."""
    try:
        profile_path = Path(MASTER_PROFILE)
        if not profile_path.exists():
            raise HTTPException(status_code=404, detail="Master profile not found")

        from Bio import AlignIO
        alignment = AlignIO.read(profile_path, "fasta")

        sequences_info = []
        for record in alignment:
            sequences_info.append({
                "id": record.id,
                "description": record.description,
                "length": len(record.seq),
                "gap_count": str(record.seq).count("-"),
                "sequence": str(record.seq),
            })

        # Generate full FASTA string
        full_alignment = ""
        for record in alignment:
            full_alignment += f">{record.description}\n{record.seq}\n"

        return {
            "profile_path": str(profile_path),
            "profile_exists": profile_path.exists(),
            "num_sequences": len(alignment),
            "alignment_length": alignment.get_alignment_length(),
            "sequences": sequences_info,
            "full_alignment": full_alignment,
            "muscle_binary": MUSCLE_BINARY,
        }

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Error reading master profile: {str(e)}")


# --- Grid Endpoints ---

@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """Generate idealized 2D grid from PDB structure with caching"""
    try:
        cache_dir = Path("cache")
        cache_dir.mkdir(exist_ok=True)
        cache_file = cache_dir / f"{pdb_id.lower()}_grid.json"

        if cache_file.exists():
            print(f"Loading cached grid data for {pdb_id.upper()}")
            with open(cache_file, "r") as f:
                cached_data = json.load(f)
                return GridData(**cached_data)

        print(f"Generating new grid data for {pdb_id.upper()}")
        grid_data = await grid_generator.generate_grid(pdb_id.lower())

        with open(cache_file, "w") as f:
            json.dump(grid_data.dict(), f, indent=2)

        return grid_data

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Internal error: {str(e)}")

# --- Root ---

@app.get("/")
async def root():
    return {
        "title": "Tubulin Spatial Grid API",
        "version": "0.2.0",
        "endpoints": {
            "/grid/{pdb_id}": "Generate 2D grid",
            "/msaprofile/sequence": "POST - Align sequence with explicit PDB numbering (auth_seq_ids)",
            "/msaprofile/master": "GET - Master alignment info",
        }
    }

if __name__ == "__main__":
    import uvicorn
    Path("debug_output").mkdir(exist_ok=True)
    uvicorn.run("api.main:app", host="127.0.0.1", port=8000, reload=True)
