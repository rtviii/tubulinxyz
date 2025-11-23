from typing import List, Optional, Dict, Any
from pydantic import BaseModel

# --- Alignment Models ---

class AlignmentRequest(BaseModel):
    sequence: str
    sequence_id: Optional[str] = None
    annotations: Optional[List[dict]] = []
    # residue_numbers corresponds to PDB Auth Seq IDs
    residue_numbers: Optional[List[int]] = None 

class MappedAnnotation(BaseModel):
    start: int
    end: int
    label: str
    metadata: Dict[str, Any] = {}

class AlignmentResponse(BaseModel):
    sequence_id: str
    aligned_sequence: str
    mapping: List[int]
    mapped_annotations: List[MappedAnnotation]
    statistics: Dict[str, Any]
    original_sequence: str

# --- Grid Models ---
# Assuming GridData is a Pydantic model imported from tubulin_analyzer.
# If it's just a class, you might want to wrap it here or import it directly in main.