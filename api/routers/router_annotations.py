# api/routers/router_annotations.py
from fastapi import APIRouter, HTTPException, Query
from typing import List, Dict, Any, Optional
from pydantic import BaseModel
from neo4j import Transaction

from neo4j_tubxz.db_lib_reader import db_reader

router_annotations = APIRouter()


# =============================================================================
# Response Models
# =============================================================================


class VariantAnnotation(BaseModel):
    """Variant annotation from the database."""

    type          : str                   # substitution, insertion, deletion
    master_index  : Optional[int] = None
    observed_index: Optional[int] = None
    wild_type     : Optional[str] = None
    observed      : Optional[str] = None
    source        : str
    uniprot_id    : Optional[str] = None
    phenotype     : Optional[str] = None
    reference     : Optional[str] = None
    rcsb_id       : str
    entity_id     : str


class PositionAnnotationsResponse(BaseModel):
    """Annotations at a specific position."""

    position: int
    family: str
    variants: List[VariantAnnotation]
    total_count: int


class PolymerAnnotationsResponse(BaseModel):
    """All annotations for a polymer chain."""

    rcsb_id: str
    auth_asym_id: str
    entity_id: str
    family: Optional[str]
    variants: List[VariantAnnotation]
    total_count: int


# =============================================================================
# Query Functions
# =============================================================================


def get_variants_at_position(position: int, family: str) -> List[Dict[str, Any]]:
    """Get all variants at a specific master alignment position for a family."""
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE v.master_index = $position AND e.family = $family
    RETURN {
        type: v.type,
        master_index: v.master_index,
        observed_index: v.observed_index,
        wild_type: v.wild_type,
        observed: v.observed,
        source: v.source,
        uniprot_id: v.uniprot_id,
        phenotype: v.phenotype,
        reference: v.reference,
        rcsb_id: e.parent_rcsb_id,
        entity_id: e.entity_id
    } AS variant
    """
    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return [
                dict(r["variant"])
                for r in tx.run(query, {"position": position, "family": family})
            ]

        return session.execute_read(run)


def get_variants_for_polymer(
    rcsb_id: str, auth_asym_id: str
) -> tuple[List[Dict[str, Any]], Optional[str], Optional[str]]:
    """Get all variants for a specific polymer chain. Returns (variants, entity_id, family)."""
    query = """
    MATCH (i:PolypeptideInstance {parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id})
    MATCH (i)-[:INSTANCE_OF]->(e:PolypeptideEntity)
    OPTIONAL MATCH (e)-[:HAS_VARIANT]->(v:Variant)
    WITH e, collect(CASE WHEN v IS NOT NULL THEN {
        type: v.type,
        master_index: v.master_index,
        observed_index: v.observed_index,
        wild_type: v.wild_type,
        observed: v.observed,
        source: v.source,
        uniprot_id: v.uniprot_id,
        phenotype: v.phenotype,
        reference: v.reference,
        rcsb_id: e.parent_rcsb_id,
        entity_id: e.entity_id
    } END) AS variants
    RETURN variants, e.entity_id AS entity_id, e.family AS family
    """
    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            result = tx.run(
                query, {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id}
            ).single()
            if not result:
                return [], None, None
            variants = [v for v in result["variants"] if v is not None]
            return variants, result["entity_id"], result["family"]

        return session.execute_read(run)


# =============================================================================
# Position-based Endpoints
# =============================================================================


@router_annotations.get(
    "/variants/{family}/{position}", response_model=PositionAnnotationsResponse
)
async def get_variants_at_position_endpoint(
    family: str, position: int
) -> PositionAnnotationsResponse:
    """
    Get all variants at a specific master alignment position for a tubulin family.
    """
    variants = get_variants_at_position(position, family)
    return PositionAnnotationsResponse(
        position=position,
        family=family,
        variants=[VariantAnnotation(**v) for v in variants],
        total_count=len(variants),
    )


@router_annotations.get("/range/{family}")
async def get_variants_in_range(
    family: str,
    start: int = Query(..., description="Start position (inclusive)"),
    end: int = Query(..., description="End position (inclusive)"),
) -> Dict[str, Any]:
    """
    Get all variants within a position range for a family.
    """
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE v.master_index >= $start AND v.master_index <= $end AND e.family = $family
    WITH v.master_index AS position, collect({
        type: v.type,
        wild_type: v.wild_type,
        observed: v.observed,
        source: v.source,
        rcsb_id: e.parent_rcsb_id
    }) AS variants
    RETURN position, variants
    ORDER BY position
    """

    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return {
                r["position"]: r["variants"]
                for r in tx.run(query, {"start": start, "end": end, "family": family})
            }

        variants_by_pos = session.execute_read(run)

    return {
        "family": family,
        "range": {"start": start, "end": end},
        "positions_with_variants": len(variants_by_pos),
        "data": variants_by_pos,
    }


# =============================================================================
# Polymer-based Endpoints
# =============================================================================


@router_annotations.get(
    "/polymer/{rcsb_id}/{auth_asym_id}", response_model=PolymerAnnotationsResponse
)
async def get_polymer_annotations(
    rcsb_id: str, auth_asym_id: str
) -> PolymerAnnotationsResponse:
    """Get all variant annotations for a specific polymer chain."""
    variants, entity_id, family = get_variants_for_polymer(rcsb_id, auth_asym_id)

    if entity_id is None:
        raise HTTPException(
            404, f"Polymer chain {auth_asym_id} not found in structure {rcsb_id}"
        )

    return PolymerAnnotationsResponse(
        rcsb_id=rcsb_id.upper(),
        auth_asym_id=auth_asym_id,
        entity_id=entity_id,
        family=family,
        variants=[VariantAnnotation(**v) for v in variants],
        total_count=len(variants),
    )


# =============================================================================
# Summary/Stats Endpoints
# =============================================================================


@router_annotations.get("/stats/{family}")
async def get_variant_stats(family: str) -> Dict[str, Any]:
    """
    Get variant statistics for a tubulin family.
    """
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE e.family = $family
    WITH v.type AS variant_type, count(*) AS count
    RETURN variant_type, count
    ORDER BY count DESC
    """

    position_query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE e.family = $family AND v.master_index IS NOT NULL
    RETURN min(v.master_index) AS min_pos, max(v.master_index) AS max_pos, count(*) AS total
    """

    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            type_counts = {
                r["variant_type"]: r["count"] for r in tx.run(query, {"family": family})
            }
            pos_stats = tx.run(position_query, {"family": family}).single()
            return type_counts, pos_stats

        type_counts, pos_stats = session.execute_read(run)

    return {
        "family": family,
        "by_type": type_counts,
        "position_range": {
            "min": pos_stats["min_pos"] if pos_stats else None,
            "max": pos_stats["max_pos"] if pos_stats else None,
        },
        "total_variants": pos_stats["total"] if pos_stats else 0,
    }
