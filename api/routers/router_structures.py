from fastapi import APIRouter, HTTPException
from typing import List, Dict, Any
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from neo4j_tubxz.db_lib_reader import dbqueries, StructureFilterParams

router_structures = APIRouter()

@router_structures.get("/all_ids", response_model=List[str])
def all_structure_ids():
    return dbqueries.all_ids()

@router_structures.get("/random")
def random_structure():
    return dbqueries.random_structure()

@router_structures.get("/{rcsb_id}")
def get_structure(rcsb_id: str):
    result = dbqueries.get_structure_by_id(rcsb_id.upper())
    if not result:
        raise HTTPException(status_code=404, detail=f"Structure {rcsb_id} not found")
    return result

@router_structures.post("/list")
def list_structures(filters: StructureFilterParams):
    structures, total_count, next_cursor = dbqueries.list_structs_filtered(filters)
    return {
        "structures": structures,
        "total_count": total_count,
        "next_cursor": next_cursor,
        "has_more": next_cursor is not None
    }

@router_structures.get("/taxa/{source_or_host}")
def get_taxa(source_or_host: str):
    """
    Returns simple list of IDs for logical ops.
    For the dropdown tree data, you might want a separate endpoint or formatting.
    """
    if source_or_host not in ["source", "host"]:
        raise HTTPException(status_code=400, detail="Must be 'source' or 'host'")
    return dbqueries.get_taxa(source_or_host)

@router_structures.get("/taxa_dict/all")
def get_taxa_dict():
    """Returns tree-compatible data for frontend dropdowns [[id, name], ...]"""
    raw_data = dbqueries.get_tax_dict()
    # Format for AntD TreeSelect: {value: id, title: name, key: id}
    return [{"value": item[0], "title": item[1], "key": item[0]} for item in raw_data]