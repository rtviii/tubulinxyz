# api/routers/router_ligands.py
from fastapi import APIRouter, Query
from typing import Optional, List

from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import LigandFilters, LigandListResponse

router_ligands = APIRouter()


@router_ligands.get("", response_model=LigandListResponse)
def list_ligands(
    cursor: Optional[str] = Query(None),
    limit: int = Query(25, ge=1, le=500),
    search: Optional[str] = Query(None),
    chemical_ids: Optional[List[str]] = Query(None, alias="ids"),
    has_drugbank: Optional[bool] = Query(None, alias="hasDrugbank"),
    in_structures: Optional[List[str]] = Query(None, alias="inStructures"),
):
    """List chemicals/ligands with filters and pagination."""
    filters = LigandFilters(
        cursor=cursor,
        limit=limit,
        search=search,
        chemical_ids=chemical_ids,
        has_drugbank=has_drugbank,
        in_structures=in_structures,
    )
    return db_reader.list_ligands(filters)

@router_ligands.get("/options", response_model=LigandListResponse)
def ligand_options(
    search: Optional[str] = Query(None, description="Search by ID or name"),
    limit: int = Query(30, ge=1, le=100),
):
    """
    Get ligands for filter dropdown/autocomplete.
    Lightweight endpoint optimized for UI selects.
    """
    filters = LigandFilters(search=search, limit=limit)
    return db_reader.list_ligands(filters)