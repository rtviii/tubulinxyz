from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
from pathlib import Path
import sys
import traceback
import os

# --- Path Setup ---
# Add parent directory to path to allow imports from sibling packages
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

# --- Imports ---
# Import the updated mapper from your existing file (typo preserved: musle_alignment)
from api.musle_alignment import TubulinAlignmentMapper, Annotation
from tubulin_analyzer import SpatialGridGenerator, GridData, DebugData

# --- App Configuration ---
app = FastAPI(
    title="Tubulin Spatial Grid API",
    version="0.1.0",
    description="Generates an idealized, aligned 2D grid from 3D tubulin structures.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Service Initialization ---
grid_generator = SpatialGridGenerator()

# Paths configuration
MUSCLE_BINARY = str(parent_dir / "muscle3.8.1")
MASTER_PROFILE = str(parent_dir / "data" / "alpha_tubulin" / "alpha_tubulin.afasta")

# Initialize Alignment Mapper
alignment_mapper = TubulinAlignmentMapper(
    master_profile_path=MASTER_PROFILE, muscle_binary=MUSCLE_BINARY
)

# --- Data Models ---


class AlignmentRequest(BaseModel):
    sequence: str
    sequence_id: Optional[str] = None
    annotations: Optional[List[dict]] = []
    residue_numbers: Optional[List[int]] = None


# --- MSA / Alignment Endpoints ---


@app.post("/msaprofile/sequence")
async def align_sequence(request: AlignmentRequest):
    """
    Align a tubulin sequence to the master alpha-tubulin profile.

    The frontend should provide 'residue_numbers' (PDB Auth Seq IDs) corresponding
    exactly to the 'sequence'. The backend will generate a mapping that links
    MSA positions back to these specific PDB IDs.
    """
    try:
        # 1. Prepare ID
        sequence_id = request.sequence_id
        if sequence_id is None:
            # Generate a stable-ish ID if none provided
            sequence_id = f"seq_{hash(request.sequence) % 10000:04d}"

        # 2. Perform Alignment
        # We pass residue_numbers to ensure the mapping matches PDB coordinates
        aligned_sequence, mapping = alignment_mapper.align_sequence(
            sequence_id, request.sequence, request.residue_numbers
        )

        # 3. Process Annotations (Optional)
        # Input annotations must use the same coordinate system as residue_numbers
        mapped_annotations = []
        if request.annotations:
            annotation_objs = [
                Annotation(
                    start=ann["start"],
                    end=ann["end"],
                    label=ann.get("label", "Feature"),
                    metadata=ann.get("metadata", {}),
                )
                for ann in request.annotations
            ]
            mapped_annotations = alignment_mapper.map_annotations(
                annotation_objs, mapping
            )

        # 4. Generate Statistics
        statistics = alignment_mapper.get_alignment_statistics(
            aligned_sequence, mapping
        )

        return {
            "sequence_id": sequence_id,
            "aligned_sequence": aligned_sequence,
            "mapping": mapping,  # Array: index=MSA_pos, value=PDB_Auth_ID (or -1)
            "mapped_annotations": [
                {
                    "start": ann.start,
                    "end": ann.end,
                    "label": ann.label,
                    "metadata": ann.metadata,
                }
                for ann in mapped_annotations
            ],
            "statistics": statistics,
            "original_sequence": request.sequence,
        }

    except ValueError as ve:
        # Catch mismatch between sequence length and residue_number length
        print(f"Validation Error: {str(ve)}")
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Alignment failed: {str(e)}")


@app.get("/msaprofile/master")
async def get_master_profile():
    """Get information about the master alignment profile including sequences."""
    try:
        profile_path = Path(MASTER_PROFILE)
        if not profile_path.exists():
            raise HTTPException(status_code=404, detail="Master profile not found")

        from Bio import AlignIO

        alignment = AlignIO.read(profile_path, "fasta")

        sequences_info = []
        for record in alignment:
            sequences_info.append(
                {
                    "id": record.id,
                    "description": record.description,
                    "length": len(record.seq),
                    "gap_count": str(record.seq).count("-"),
                    "sequence": str(record.seq),
                }
            )

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
            "muscle_exists": Path(MUSCLE_BINARY).exists(),
        }

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(
            status_code=500, detail=f"Error reading master profile: {str(e)}"
        )


# --- Grid & Structure Endpoints ---


@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """Generate idealized 2D grid from PDB structure with caching"""
    try:
        cache_dir = Path("cache")
        cache_dir.mkdir(exist_ok=True)
        cache_file = cache_dir / f"{pdb_id.lower()}_grid.json"

        if cache_file.exists():
            print(f"Loading cached grid data for {pdb_id.upper()}")
            import json

            with open(cache_file, "r") as f:
                cached_data = json.load(f)
                return GridData(**cached_data)

        print(f"Generating new grid data for {pdb_id.upper()}")
        grid_data = await grid_generator.generate_grid(pdb_id.lower())

        with open(cache_file, "w") as f:
            import json

            json.dump(grid_data.dict(), f, indent=2)

        return grid_data

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(
            status_code=500, detail=f"An internal error occurred: {str(e)}"
        )


@app.post("/align/structure/{pdb_id}/{chain_id}")
async def align_structure_chain(pdb_id: str, chain_id: str, annotations: list = None):
    """
    Legacy/Placeholder: Align a specific chain from backend extraction.
    Note: Prefer /msaprofile/sequence with frontend-extracted data for visualization consistency.
    """
    try:
        # Placeholder logic
        return {"message": "Use /msaprofile/sequence for explicit frontend alignment"}
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(
            status_code=500, detail=f"Structure alignment failed: {str(e)}"
        )


# --- Debug Endpoints ---


# --- Root ---


@app.get("/")
async def root():
    return {
        "title": "Tubulin Spatial Grid API",
        "endpoints": {
            "/grid/{pdb_id}": "Generate 2D grid",
            "/msaprofile/sequence": "POST - Align sequence with explicit PDB numbering",
            "/msaprofile/master": "GET - Master alignment info",
        },
    }


if __name__ == "__main__":
    import uvicorn

    Path("debug_output").mkdir(exist_ok=True)
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)
