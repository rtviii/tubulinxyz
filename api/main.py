import sys
sys.path.append('/Users/rtviii/dev/tubulinxyz')
from pathlib import Path
import traceback
import json
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

# Internal Modules
from api.config import settings
from api.schemas import AlignmentRequest, AlignmentResponse
from api.services.alignment import TubulinAlignmentMapper, Annotation

# External/Existing Libraries
from lib.tubulin_analyzer import SpatialGridGenerator, GridData

# --- App Init ---
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

mapper = TubulinAlignmentMapper(
    master_profile_path=settings.MASTER_PROFILE, 
    muscle_binary=settings.MUSCLE_BINARY
)

# --- Endpoints ---

@app.post("/msaprofile/sequence", response_model=AlignmentResponse)
async def align_sequence(request: AlignmentRequest): 
    """Align a tubulin sequence to the master alpha-tubulin profile."""
    try:
        # 1. Prepare ID
        seq_id = request.sequence_id or f"seq_{hash(request.sequence) % 10000:04d}"

        # 2. Logic: Align
        aligned_seq, mapping = mapper.align_sequence(
            seq_id, request.sequence, request.residue_numbers
        )

        # 3. Logic: Map Annotations
        mapped_anns = []
        if request.annotations:
            # Convert Schema -> Service Object
            service_anns = [
                Annotation(
                    start=a["start"], end=a["end"], 
                    label=a.get("label", "Feature"), metadata=a.get("metadata", {})
                ) for a in request.annotations
            ]
            mapped_anns = mapper.map_annotations(service_anns, mapping)

        # 4. Logic: Stats
        stats = mapper.get_alignment_statistics(aligned_seq, mapping)

        return {
            "sequence_id": seq_id,
            "aligned_sequence": aligned_seq,
            "mapping": mapping,
            "mapped_annotations": [
                {"start": a.start, "end": a.end, "label": a.label, "metadata": a.metadata}
                for a in mapped_anns
            ],
            "statistics": stats,
            "original_sequence": request.sequence,
        }

    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Alignment failed: {str(e)}")


@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """Generate idealized 2D grid from PDB structure with file-based caching."""
    grid_generator = SpatialGridGenerator() # Instantiating here, or could be global
    
    try:
        cache_file = settings.CACHE_DIR / f"{pdb_id.lower()}_grid.json"

        if cache_file.exists():
            print(f"Loading cached grid data for {pdb_id.upper()}")
            with open(cache_file, "r") as f:
                return GridData(**json.load(f))

        print(f"Generating new grid data for {pdb_id.upper()}")
        grid_data = await grid_generator.generate_grid(pdb_id.lower())

        with open(cache_file, "w") as f:
            json.dump(grid_data.dict(), f, indent=2)

        return grid_data

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Grid generation error: {str(e)}")


@app.get("/msaprofile/master")
async def get_master_profile():
    """Get information about the master alignment profile."""
    from Bio import AlignIO
    
    profile_path = Path(settings.MASTER_PROFILE)
    if not profile_path.exists():
        raise HTTPException(status_code=404, detail="Master profile not found")

    try:
        alignment = AlignIO.read(profile_path, "fasta")
        
        # --- RESTORED LOGIC START ---
        sequences_info = []
        for record in alignment:
            sequences_info.append({
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq),
                "length": len(record.seq),
                "gap_count": str(record.seq).count("-"),
            })
        return {
            "profile_path"    : str(profile_path),
            "num_sequences"   : len(alignment),
            "alignment_length": alignment.get_alignment_length(),
            "sequences"       : sequences_info,                     
            "full_alignment"  : format(alignment, "fasta"),
        }
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/")
async def root():
    return {
        "service": "Tubulin Spatial Grid API",
        "endpoints": ["/grid/{pdb_id}", "/msaprofile/sequence", "/msaprofile/master"]
    }

if __name__ == "__main__":
    import uvicorn
    # Note: reload=True works best when referencing the module path
    uvicorn.run("api.main:app", host="127.0.0.1", port=8000, reload=True)
