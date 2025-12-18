from fastapi import APIRouter
from typing import List
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from neo4j_tubxz.db_lib_reader import dbqueries, PolymersFilterParams
from lib.types import TubulinFamily

router_polymers = APIRouter()

@router_polymers.post("/list")
def list_polymers(filters: PolymersFilterParams):
    proteins, total_structs, total_polys, next_cursor = dbqueries.list_polymers_filtered(filters)
    
    return {
        "proteins": proteins,
        "total_count": total_polys,
        "total_paren_structures_count": total_structs,
        "next_cursor": next_cursor,
        "has_more": next_cursor is not None
    }

@router_polymers.get("/families", response_model=List[str])
def get_tubulin_families():
    return [f.value for f in TubulinFamily]