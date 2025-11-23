# api/schemas.py
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