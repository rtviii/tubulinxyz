# api/routers/router_ligands.py
from fastapi import APIRouter, Query, HTTPException
from typing import Optional, List
from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import (
    LigandFilters,
    LigandListResponse,
    PolymerNeighborhoodsResponse,
    CanonicalBindingSite,
)

router_ligands = APIRouter()


@router_ligands.get("", response_model=LigandListResponse, operation_id="list_ligands")
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


@router_ligands.get(
    "/options",
    response_model=LigandListResponse,
    operation_id="list_ligand_options",
)
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


@router_ligands.get(
    "/canonical-site/{chemical_id}/{family}",
    response_model=CanonicalBindingSite,
    operation_id="get_canonical_binding_site",
)
def get_canonical_binding_site(chemical_id: str, family: str):
    """
    Aggregate binding site for a ligand across all structures of a given family.

    Returns per-master-alignment-position frequencies: how often each residue position
    appears in the binding site across all structures in the database where this ligand
    contacts a polymer of the specified family.

    Useful for painting a "canonical" binding site heatmap onto any structure,
    even if that particular structure doesn't have the ligand bound.
    """
    result = db_reader.get_canonical_binding_site(chemical_id, family)
    if result is None:
        raise HTTPException(
            404,
            f"No binding data for {chemical_id.upper()} in {family}. "
            f"Either the ligand doesn't exist or it has never been observed "
            f"near a {family} polymer in any structure.",
        )
    return result


@router_ligands.get(
    "/neighborhoods/{rcsb_id}/{auth_asym_id}",
    response_model=PolymerNeighborhoodsResponse,
    operation_id="get_polymer_ligand_neighborhoods",
)
def get_ligand_neighborhoods(rcsb_id: str, auth_asym_id: str):
    """
    Get all ligand neighborhoods for a specific polymer chain.
    Returns ligands that are spatially near the polymer chain along with
    the residues that form the binding interface.
    """
    result = db_reader.get_ligand_neighborhoods_for_polymer(rcsb_id, auth_asym_id)
    if result.total_ligands == 0:
        # Check if the polymer exists at all
        exists_query = """
        MATCH (pi:PolypeptideInstance {parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id})
        RETURN count(pi) > 0 AS exists
        """
        with db_reader.adapter.driver.session() as session:
            exists = session.run(
                exists_query,
                {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id},
            ).single()["exists"]
        if not exists:
            raise HTTPException(
                404,
                f"Polymer chain {auth_asym_id} not found in structure {rcsb_id}",
            )
    return result
