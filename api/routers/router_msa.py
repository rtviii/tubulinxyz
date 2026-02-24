# api/routers/router_msa.py
from fastapi import APIRouter, HTTPException, Query
from pathlib import Path
import traceback
from functools import lru_cache
from enum import Enum

from api.config import settings
from lib.seq_aligner import TubulinAlignmentMapper

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field, model_validator

# ----------------------------
# Family enum (canonical keys)
# ----------------------------
class TubulinFamily(str, Enum):
    ALPHA = "tubulin_alpha"
    BETA = "tubulin_beta"
    GAMMA = "tubulin_gamma"
    DELTA = "tubulin_delta"
    EPSILON = "tubulin_epsilon"

family_file_map = {
    TubulinFamily.ALPHA:   "/Users/rtviii/dev/tubulinxyz/data/alpha_tubulin/alpha_tubulin.afasta",
    TubulinFamily.BETA:    "/Users/rtviii/dev/tubulinxyz/data/beta_tubulin/beta_tubulin.afasta",
    TubulinFamily.GAMMA:   "/Users/rtviii/dev/tubulinxyz/data/gamma_tubulin/tubulin_gamma_clean.afasta",
    TubulinFamily.DELTA:   "/Users/rtviii/dev/tubulinxyz/data/delta_tubulin/tubulin_delta_clean.afasta",
    TubulinFamily.EPSILON: "/Users/rtviii/dev/tubulinxyz/data/epsilon_tubulin/tubulin_epsilon_clean.afasta",
}

@lru_cache(maxsize=16)
def get_alignment_mapper(family: TubulinFamily) -> TubulinAlignmentMapper:
    profile_path = Path(family_file_map[family])
    if not profile_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Master profile not found for {family}: {profile_path}"
        )
    return TubulinAlignmentMapper(
        master_profile_path=str(profile_path),
        muscle_binary=settings.MUSCLE_BINARY
    )

router_msa = APIRouter()

# ----------------------------
# Models (unchanged)
# ----------------------------
class Annotation(BaseModel):
    start: int
    end: int
    label: str
    metadata: Dict[str, Any] = {}

class AlignmentRequest(BaseModel):
    sequence: str = Field(..., description="The observed single-letter amino acid sequence")
    sequence_id: Optional[str] = Field(None, description="Unique identifier for the sequence")
    annotations: Optional[List[dict]] = []
    auth_seq_ids: Optional[List[int]] = Field(None, description="PDB auth_seq_ids, 1:1 with sequence")

    @model_validator(mode='after')
    def validate_lengths(self) -> 'AlignmentRequest':
        if self.auth_seq_ids is not None and len(self.auth_seq_ids) != len(self.sequence):
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

class MasterProfileInfo(BaseModel):
    profile_path: str
    profile_exists: bool
    num_sequences: int
    alignment_length: int
    sequences: List[Dict[str, Any]]
    full_alignment: str
    muscle_binary: str

# ----------------------------
# POST /sequence (FAMILY AWARE)
# ----------------------------
@router_msa.post("/sequence", response_model=AlignmentResponse, operation_id="align_sequence")
async def align_sequence(
    request: AlignmentRequest,
    family: TubulinFamily = Query(..., description="Which master alignment to align against"),
):
    try:
        mapper = get_alignment_mapper(family)

        aligned_sequence, mapping = mapper.align_sequence_with_mapping(
            request.sequence_id or f"seq_{family}",
            request.sequence,
            request.auth_seq_ids
        )

        return {
            "sequence_id": request.sequence_id or f"seq_{family}",
            "aligned_sequence": aligned_sequence,
            "mapping": mapping,
            "mapped_annotations": [],
            "statistics": {
                "family": family,
                "master_alignment_length": mapper.ref_len,
            },
            "original_sequence": request.sequence
        }

    except HTTPException:
        raise
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Alignment failed: {str(e)}")

# ----------------------------
# GET /master (FAMILY AWARE)
# ----------------------------
@router_msa.get("/master", response_model=MasterProfileInfo, operation_id="get_master_profile")
async def get_master_profile(
    family: TubulinFamily = Query(..., description="Which master alignment to return"),
):
    try:
        mapper = get_alignment_mapper(family)
        profile_path = Path(mapper.master_profile_path)

        from Bio import AlignIO
        alignment = AlignIO.read(profile_path, "fasta")

        sequences_info = [
            {
                "id": record.id,
                "description": record.description,
                "length": len(record.seq),
                "gap_count": str(record.seq).count("-"),
                "sequence": str(record.seq),
                "family": family,   # <- add this (nice for frontend)
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

    except HTTPException:
        raise
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Error reading master profile: {str(e)}")
