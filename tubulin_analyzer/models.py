from pydantic import BaseModel
from typing import List, Dict, Any


class SubunitData(BaseModel):
    id: str
    auth_asym_id: str
    protofilament: int
    subunitIndex: int
    monomerType: str


class GridData(BaseModel):
    subunits: List[SubunitData]
    structure_type: str
    metadata: Dict[str, Any]


class DebugData(BaseModel):
    pdb_id: str
    num_tubulin_chains: int
    num_connections: int
    num_protofilaments: int
    protofilament_lengths: List[int]
    connection_success_rate: float
    debug_files: List[str]