# api/routers/router_msa.py
from fastapi import APIRouter, HTTPException, Query
from pathlib import Path
import traceback
from functools import lru_cache

from api.config import settings

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field, model_validator

from lib.types import TubulinFamily
from lib.etl.sequence_alignment import SequenceAligner, get_aligner_for_family

router_msa = APIRouter()


@lru_cache(maxsize=16)
def get_cached_aligner(family: TubulinFamily) -> SequenceAligner:
    try:
        return get_aligner_for_family(family, Path(settings.PROJECT_ROOT))
    except Exception as e:
        raise HTTPException(status_code=404, detail=f"MSA not found for {family}: {e}")

# ----------------------------
# Models
# ----------------------------
class Annotation(BaseModel):
    start: int
    end: int
    label: str
    metadata: Dict[str, Any] = {}


class AlignmentRequest(BaseModel):
    sequence: str = Field(
        ..., description="The observed single-letter amino acid sequence"
    )
    sequence_id: Optional[str] = Field(
        None, description="Unique identifier for the sequence"
    )
    annotations: Optional[List[dict]] = []
    auth_seq_ids: Optional[List[int]] = Field(
        None, description="PDB auth_seq_ids, 1:1 with sequence"
    )

    @model_validator(mode="after")
    def validate_lengths(self) -> "AlignmentRequest":
        if self.auth_seq_ids is not None and len(self.auth_seq_ids) != len(
            self.sequence
        ):
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
# POST /sequence
# ----------------------------
@router_msa.post(
    "/sequence", response_model=AlignmentResponse, operation_id="align_sequence"
)
async def align_sequence(
    request: AlignmentRequest,
    family: TubulinFamily = Query(
        ..., description="Which master alignment to align against"
    ),
):
    try:
        aligner = get_cached_aligner(family)
        seq_id = request.sequence_id or f"seq_{family}"
        auth_seq_ids = request.auth_seq_ids or list(range(1, len(request.sequence) + 1))

        # Run MUSCLE alignment

        aln_target, is_original_col = aligner._run_muscle(
            seq_id, request.sequence
        )

        ref_len = aligner.consensus.length
        ma_to_auth = [-1] * ref_len

        master_idx = 0
        target_idx = 0

        for col, t_char in enumerate(aln_target):
            if col >= len(is_original_col):
                break
            if not is_original_col[col]:
                if t_char != "-":
                    target_idx += 1
            else:
                if (
                    t_char != "-"
                    and master_idx < ref_len
                    and target_idx < len(auth_seq_ids)
                ):
                    ma_to_auth[master_idx] = auth_seq_ids[target_idx]
                    target_idx += 1
                elif t_char == "-":
                    pass
                else:
                    target_idx += 1
                master_idx += 1

        return {
            "sequence_id": seq_id,
            "aligned_sequence": aln_target,
            "mapping": ma_to_auth,
            "mapped_annotations": [],
            "statistics": {
                "family": family,
                "master_alignment_length": ref_len,
            },
            "original_sequence": request.sequence,
        }

    except HTTPException:
        raise
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Alignment failed: {str(e)}")


# ----------------------------
# GET /master
# ----------------------------
@router_msa.get(
    "/master", response_model=MasterProfileInfo, operation_id="get_master_profile"
)
async def get_master_profile(
    family: TubulinFamily = Query(..., description="Which master alignment to return"),
):
    try:
        aligner = get_cached_aligner(family)
        profile_path = aligner.msa_path

        from Bio import AlignIO

        alignment = AlignIO.read(str(profile_path), "fasta")

        sequences_info = [
            {
                "id"         : record.id,
                "description": record.description,
                "length"     : len(record.seq),
                "gap_count"  : str(record.seq).count("-"),
                "sequence"   : str(record.seq),
                "family"     : family,
            }
            for record in alignment
        ]

        full_alignment = "".join(
            f">{record.description}\n{record.seq}\n" for record in alignment
        )

        return {
            "profile_path"    : str(profile_path),
            "profile_exists"  : profile_path.exists(),
            "num_sequences"   : len(alignment),
            "alignment_length": alignment.get_alignment_length(),
            "sequences"       : sequences_info,
            "full_alignment"  : full_alignment,
            "muscle_binary"   : str(aligner.muscle_binary),
        }

    except HTTPException:
        raise
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(
            status_code=500, detail=f"Error reading master profile: {str(e)}"
        )
