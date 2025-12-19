# api/routers/router_structures.py
from fastapi import APIRouter, Query, HTTPException
from typing import Optional, List

from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import (
    FilterFacets,
    StructureFilters,
    StructureListResponse,
    ExpMethod,
    PolymerizationState,
)

router_structures = APIRouter()

# api/routers/router_structures.py - add this endpoint


@router_structures.get("/taxonomy-tree/{tax_type}")
def get_taxonomy_tree(tax_type: str = "source"):
    """Get taxonomy as tree structure for UI TreeSelect component."""
    if tax_type not in ("source", "host"):
        raise HTTPException(400, "tax_type must be 'source' or 'host'")
    return db_reader.get_taxonomy_tree_for_ui(tax_type)


def parse_list_param(value: Optional[List[str]]) -> Optional[List[str]]:
    """Handle both repeated params and comma-separated values."""
    if not value:
        return None
    # If it's a single string with commas, split it
    if len(value) == 1 and "," in value[0]:
        return value[0].split(",")
    return value


@router_structures.get("", response_model=StructureListResponse)
def list_structures(
    cursor: Optional[str] = Query(None),
    limit: int = Query(100, ge=1, le=500),
    search: Optional[str] = Query(None),
    rcsb_ids: Optional[List[str]] = Query(None, alias="ids"),
    resolution_min: Optional[float] = Query(None, alias="resMin"),
    resolution_max: Optional[float] = Query(None, alias="resMax"),
    year_min: Optional[int] = Query(None, alias="yearMin"),
    year_max: Optional[int] = Query(None, alias="yearMax"),
    exp_method: Optional[List[str]] = Query(None, alias="expMethod"),
    polymerization_state: Optional[List[str]] = Query(None, alias="polyState"),
    source_organism_ids: Optional[List[int]] = Query(None, alias="sourceTaxa"),
    host_organism_ids: Optional[List[int]] = Query(None, alias="hostTaxa"),
    has_ligand_ids: Optional[List[str]] = Query(None, alias="ligands"),
    has_polymer_family: Optional[List[str]] = Query(None, alias="family"),
    has_uniprot: Optional[List[str]] = Query(None, alias="uniprot"),
    has_mutations: Optional[bool] = Query(None, alias="hasMutations"),
    mutation_family: Optional[str] = Query(None, alias="mutationFamily"),
    mutation_position_min: Optional[int] = Query(None, alias="mutationPosMin"),
    mutation_position_max: Optional[int] = Query(None, alias="mutationPosMax"),
    mutation_from: Optional[str] = Query(None, alias="mutationFrom"),
    mutation_to: Optional[str] = Query(None, alias="mutationTo"),
    mutation_phenotype: Optional[str] = Query(None, alias="mutationPhenotype"),
):
    """List structures with cumulative filters and keyset pagination."""

    # Parse comma-separated list params
    parsed_families = parse_list_param(has_polymer_family)
    parsed_poly_state = parse_list_param(polymerization_state)
    parsed_exp_method = parse_list_param(exp_method)
    parsed_ligands = parse_list_param(has_ligand_ids)
    parsed_uniprot = parse_list_param(has_uniprot)
    parsed_ids = parse_list_param(rcsb_ids)

    filters = StructureFilters(
        cursor=cursor,
        limit=limit,
        search=search,
        rcsb_ids=parsed_ids,
        resolution_min=resolution_min,
        resolution_max=resolution_max,
        year_min=year_min,
        year_max=year_max,
        exp_method=[ExpMethod(m) for m in parsed_exp_method]
        if parsed_exp_method
        else None,
        polymerization_state=[PolymerizationState(s) for s in parsed_poly_state]
        if parsed_poly_state
        else None,
        source_organism_ids   = source_organism_ids,
        host_organism_ids     = host_organism_ids,
        has_ligand_ids        = parsed_ligands,
        has_polymer_family    = parsed_families,
        has_uniprot           = parsed_uniprot,
        has_mutations         = has_mutations,
        mutation_family       = mutation_family,
        mutation_position_min = mutation_position_min,
        mutation_position_max = mutation_position_max,
        mutation_from         = mutation_from,
        mutation_to           = mutation_to,
        mutation_phenotype    = mutation_phenotype,
    )

    return db_reader.list_structures(filters)


# Update the facets endpoint to have a response model
@router_structures.get("/facets", response_model=FilterFacets)
def get_facets():
    """Get available filter options for UI dropdowns."""
    return db_reader.get_filter_facets()


@router_structures.get("/taxonomy/{tax_type}")
def get_taxonomy(tax_type: str = "source"):
    """Get taxonomy options for filter dropdowns."""
    if tax_type not in ("source", "host"):
        raise HTTPException(400, "tax_type must be 'source' or 'host'")
    return db_reader.get_taxonomy_tree(tax_type)


@router_structures.get("/families")
def get_families():
    """Get tubulin family options with counts."""
    return db_reader.get_tubulin_families()


@router_structures.get("/{rcsb_id}")
def get_structure(rcsb_id: str):
    """Get full structure details."""
    result = db_reader.get_structure(rcsb_id)
    if not result:
        raise HTTPException(404, f"Structure {rcsb_id} not found")
    return result
