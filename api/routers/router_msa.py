# api/routers/router_msa.py
from fastapi import APIRouter, HTTPException
from pathlib import Path
import traceback

from api.config import settings
from lib.seq_aligner import TubulinAlignmentMapper


from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field, model_validator

class Annotation(BaseModel):
    start: int
    end: int
    label: str
    metadata: Dict[str, Any] = {}

class AlignmentRequest(BaseModel):
    sequence: str = Field(..., description="The observed single-letter amino acid sequence")
    sequence_id: Optional[str] = Field(None, description="Unique identifier for the sequence")
    annotations: Optional[List[dict]] = []
    
    auth_seq_ids: Optional[List[int]] = Field(
        None, 
        description="The PDB auth_seq_ids corresponding 1:1 to the sequence characters."
    )

    @model_validator(mode='after')
    def validate_lengths(self) -> 'AlignmentRequest':
        if self.auth_seq_ids is not None:
            if len(self.auth_seq_ids) != len(self.sequence):
                raise ValueError(
                    f"Length Mismatch: Sequence is {len(self.sequence)} chars, "
                    f"but {len(self.auth_seq_ids)} auth_seq_ids were provided."
                )
        return self

class AlignmentResponse(BaseModel):
    sequence_id: str
    aligned_sequence: str
    mapping: List[int]
    mapped_annotations: List[Annotation]
    statistics: Dict[str, Any]
    original_sequence: str

# Add this model:

class MasterProfileInfo(BaseModel):
    """Master alignment profile metadata and sequences."""
    profile_path: str
    profile_exists: bool
    num_sequences: int
    alignment_length: int
    sequences: List[Dict[str, Any]]
    full_alignment: str
    muscle_binary: str


# Endpoint signature changes:

router_msa = APIRouter()

# Initialize service
alignment_mapper = TubulinAlignmentMapper(
    master_profile_path=settings.MASTER_PROFILE,
    muscle_binary=settings.MUSCLE_BINARY
)



@router_msa.post("/sequence", response_model=AlignmentResponse, operation_id="align_sequence")
async def align_sequence(request: AlignmentRequest):
    """Align a sequence against the master profile and return mapping."""
    if not alignment_mapper:
        raise HTTPException(status_code=503, detail="Alignment service unavailable")

    try:
        aligned_sequence, mapping = alignment_mapper.align_sequence_with_mapping(
            request.sequence_id,
            request.sequence,
            request.auth_seq_ids
        )

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

@router_msa.get("/master", response_model=MasterProfileInfo, operation_id="get_master_profile")
async def get_master_profile():
    """Get information about the master alignment profile including raw sequences."""
    try:
        profile_path = Path(settings.MASTER_PROFILE)
        if not profile_path.exists():
            raise HTTPException(status_code=404, detail="Master profile not found")

        from Bio import AlignIO
        alignment = AlignIO.read(profile_path, "fasta")

        sequences_info = [
            {
                "id": record.id,
                "description": record.description,
                "length": len(record.seq),
                "gap_count": str(record.seq).count("-"),
                "sequence": str(record.seq),
            }
            for record in alignment
        ]

        full_alignment = "".join(
            f">{record.description}\n{record.seq}\n" for record in alignment
        )

        return {
            "profile_path": str(profile_path),
            "profile_exists": profile_path.exists(),
            "num_sequences": len(alignment),
            "alignment_length": alignment.get_alignment_length(),
            "sequences": sequences_info,
            "full_alignment": full_alignment,
            "muscle_binary": settings.MUSCLE_BINARY,
        }

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Error reading master profile: {str(e)}")