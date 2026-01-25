# api/routers/router_structures.py
from fastapi import APIRouter, Query, HTTPException
from typing import Any, Dict, Optional, List
import json
from pathlib import Path

from pydantic import BaseModel

from lib.types import TubulinStructure
from lib.etl.assets import TubulinStructureAssets
from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import (
    FilterFacets,
    StructureFilters,
    StructureListResponse,
    ExpMethod,
    PolymerizationState,
    VariantTypeFilter,
)

router_structures = APIRouter()

# Add these models near the top of the file:

class TaxonomyTreeNode(BaseModel):
    """Tree node for UI TreeSelect component."""
    value: int
    title: str
    children: Optional[List["TaxonomyTreeNode"]] = None

TaxonomyTreeNode.model_rebuild()  # For self-referential type


class TaxonomyFlatNode(BaseModel):
    """Flat taxonomy node with counts."""
    tax_id: int
    name: str
    rank: Optional[str] = None
    structure_count: int


class FamilyCount(BaseModel):
    """Family with structure count."""
    family: str
    count: int


class StructureDetail(BaseModel):
    """Full structure with related entities."""
    structure: Dict[str, Any]
    polypeptide_entities: List[Dict[str, Any]]
    ligand_entities: List[Dict[str, Any]]
    polypeptide_instances: List[Dict[str, Any]]
    ligand_instances: List[Dict[str, Any]]


# Endpoint signature changes:


@router_structures.get("/taxonomy-tree/{tax_type}", response_model=List[TaxonomyTreeNode], operation_id="get_taxonomy_tree")
def get_taxonomy_tree(tax_type: str = "source"):
    """Get taxonomy as tree structure for UI TreeSelect component."""
    if tax_type not in ("source", "host"):
        raise HTTPException(400, "tax_type must be 'source' or 'host'")
    return db_reader.get_taxonomy_tree_for_ui(tax_type)


def parse_list_param(value: Optional[List[str]]) -> Optional[List[str]]:
    if not value:
        return None
    result = []
    for item in value:
        if "," in item:
            result.extend([x.strip() for x in item.split(",") if x.strip()])
        else:
            result.append(item.strip())
    return result if result else None


def parse_int_list(value: Optional[List[str]]) -> Optional[List[int]]:
    """Handle both repeated params and comma-separated values for integers."""
    raw = parse_list_param(value)
    if raw:
        try:
            return [int(x) for x in raw if x.strip()]
        except ValueError:
            raise HTTPException(400, "Taxonomy IDs must be integers")
    return None


@router_structures.get("", response_model=StructureListResponse, operation_id="list_structures")
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
    source_organism_ids: Optional[List[str]] = Query(None, alias="sourceTaxa"),
    host_organism_ids: Optional[List[str]] = Query(None, alias="hostTaxa"),
    has_ligand_ids: Optional[List[str]] = Query(None, alias="ligands"),
    has_polymer_family: Optional[List[str]] = Query(None, alias="family"),
    has_uniprot: Optional[List[str]] = Query(None, alias="uniprot"),
    # Variant filters
    has_variants: Optional[bool] = Query(None, alias="hasVariants"),
    variant_family: Optional[str] = Query(None, alias="variantFamily"),
    variant_type: Optional[str] = Query(None, alias="variantType"),
    variant_position_min: Optional[int] = Query(None, alias="variantPosMin"),
    variant_position_max: Optional[int] = Query(None, alias="variantPosMax"),
    variant_wild_type: Optional[str] = Query(None, alias="variantWildType"),
    variant_observed: Optional[str] = Query(None, alias="variantObserved"),
    variant_source: Optional[str] = Query(None, alias="variantSource"),
    variant_phenotype: Optional[str] = Query(None, alias="variantPhenotype"),
):
    """List structures with cumulative filters and keyset pagination."""

    # Parse comma-separated list params
    parsed_families = parse_list_param(has_polymer_family)
    parsed_poly_state = parse_list_param(polymerization_state)
    parsed_exp_method = parse_list_param(exp_method)
    parsed_ligands = parse_list_param(has_ligand_ids)
    parsed_uniprot = parse_list_param(has_uniprot)
    parsed_ids = parse_list_param(rcsb_ids)

    parsed_source_taxa = parse_int_list(source_organism_ids)
    parsed_host_taxa = parse_int_list(host_organism_ids)

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
        source_organism_ids=parsed_source_taxa,
        host_organism_ids=parsed_host_taxa,
        has_ligand_ids=parsed_ligands,
        has_polymer_family=parsed_families,
        has_uniprot=parsed_uniprot,
        has_variants=has_variants,
        variant_family=variant_family,
        variant_type=VariantTypeFilter(variant_type) if variant_type else None,
        variant_position_min=variant_position_min,
        variant_position_max=variant_position_max,
        variant_wild_type=variant_wild_type,
        variant_observed=variant_observed,
        variant_source=variant_source,
        variant_phenotype=variant_phenotype,
    )

    return db_reader.list_structures(filters)


@router_structures.get("/facets", response_model=FilterFacets, operation_id="get_structure_facets")
def get_facets():
    """Get available filter options for UI dropdowns."""
    return db_reader.get_filter_facets()


@router_structures.get("/taxonomy/{tax_type}", response_model=List[TaxonomyFlatNode], operation_id="get_taxonomy_flat")
def get_taxonomy(tax_type: str = "source"):
    """Get taxonomy options for filter dropdowns."""
    if tax_type not in ("source", "host"):
        raise HTTPException(400, "tax_type must be 'source' or 'host'")
    return db_reader.get_taxonomy_tree(tax_type)


@router_structures.get("/families", response_model=List[FamilyCount], operation_id="list_families")
def get_families():
    """Get tubulin family options with counts."""
    return db_reader.get_tubulin_families()

@router_structures.get("/{rcsb_id}", response_model=StructureDetail, operation_id="get_structure")
def get_structure(rcsb_id: str):
    """Get full structure details."""
    result = db_reader.get_structure(rcsb_id)
    if not result:
        raise HTTPException(404, f"Structure {rcsb_id} not found")
    return result


@router_structures.get("/{rcsb_id}/profile", response_model=TubulinStructure, operation_id="get_structure_profile")
async def get_structure_profile(rcsb_id: str):
    """Fetches the pre-calculated TubulinStructure JSON profile from disk."""
    assets = TubulinStructureAssets(rcsb_id.upper())
    profile_path = Path(assets.paths.profile)

    if not profile_path.exists():
        raise HTTPException(
            status_code=404, detail=f"Profile for {rcsb_id} not collected yet."
        )

    with open(profile_path, "r") as f:
        return json.load(f)
