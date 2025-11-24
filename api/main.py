from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pathlib import Path
import sys
import traceback
import json
from datetime import datetime

# --- Path Setup ---
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

# --- Imports ---
from api.schemas import AlignmentRequest, AlignmentResponse
from api.services.alignment import (
    TubulinAlignmentMapper, 
    TubulinStructureParser,
    TubulinIngestor  # Add this
)

from lib.tubulin_analyzer import SpatialGridGenerator, GridData, DebugData

# --- App Configuration ---
app = FastAPI(
    title="Tubulin Spatial Grid API",
    version="0.2.0",
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

# Create ingestion logs directory
INGESTION_LOGS_DIR = parent_dir / "ingestion_logs"
INGESTION_LOGS_DIR.mkdir(exist_ok=True)

# Initialize services
alignment_mapper = TubulinAlignmentMapper(
    master_profile_path=MASTER_PROFILE, 
    muscle_binary=MUSCLE_BINARY
)
mmcif_parser = TubulinStructureParser()
ingestor = TubulinIngestor(MASTER_PROFILE, MUSCLE_BINARY)

# --- MSA / Alignment Endpoints ---

@app.post("/msaprofile/sequence", response_model=AlignmentResponse)
async def align_sequence(request: AlignmentRequest):
    if not alignment_mapper:
        raise HTTPException(
            status_code=503, 
            detail="Alignment service unavailable (configuration error)"
        )

    try:
        # 1. Main alignment for frontend
        aligned_sequence, mapping = alignment_mapper.align_sequence_with_mapping(
            request.sequence_id, 
            request.sequence,
            request.auth_seq_ids
        )
        
        # 2. Verification
        if request.sequence_id and "_" in request.sequence_id:
            pdb_id, chain_id = request.sequence_id.split("_")[:2]
            print(f"Verifying {pdb_id} chain {chain_id} against RCSB...")
            
            verification = mmcif_parser.verify_integrity(
                pdb_id, chain_id, request.sequence, request.auth_seq_ids
            )
            
            if verification["status"] == "success":
                if verification["match"]:
                    print(f"✅ VERIFIED: {pdb_id}_{chain_id} matches PDB")
                else:
                    print(f"⚠️ MISMATCH: {pdb_id}_{chain_id}")
                    print(verification["details"])
            
            # 3. Run full ingestion pipeline (for logging/analysis)
            try:
                print(f"\n{'='*60}")
                print(f"RUNNING INGESTION PIPELINE FOR {pdb_id} CHAIN {chain_id}")
                print(f"{'='*60}")
                
                tubulin_class = "Alpha"  # Default, adjust as needed
                result = ingestor.process_chain(pdb_id, chain_id, tubulin_class)
                
                # Save to timestamped file
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_file = INGESTION_LOGS_DIR / f"{pdb_id}_{chain_id}_{timestamp}.json"
                
                with open(output_file, "w") as f:
                    from dataclasses import asdict
                    json.dump(asdict(result), f, indent=2)
                
                print(f"\n📊 INGESTION STATS:")
                print(f"  - Mutations detected: {result.stats['total_mutations']}")
                print(f"  - Insertions: {result.stats['insertions']}")
                print(f"  - MA Coverage: {result.stats['ma_coverage']}/{len(result.ma_to_auth_map)}")
                print(f"  - Results saved to: {output_file.name}")
                
                if result.mutations:
                    print(f"\n🧬 MUTATIONS:")
                    for mut in result.mutations[:5]:  # Show first 5
                        print(f"  - MA pos {mut.ma_position}: {mut.wild_type} → {mut.observed} (PDB ID: {mut.pdb_auth_id})")
                    if len(result.mutations) > 5:
                        print(f"  ... and {len(result.mutations) - 5} more")
                
                print(f"{'='*60}\n")
                
            except Exception as ing_error:
                # Don't fail the whole request if ingestion fails
                print(f"⚠️ INGESTION PIPELINE FAILED (non-critical): {ing_error}")
                traceback.print_exc()
        
        # 4. Return frontend response (unchanged)
        return {
            "sequence_id": request.sequence_id,
            "aligned_sequence": aligned_sequence,
            "mapping": mapping,
            "mapped_annotations": [],
            "statistics": {},
            "original_sequence": request.sequence
        }
        
    except ValueError as ve:
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
