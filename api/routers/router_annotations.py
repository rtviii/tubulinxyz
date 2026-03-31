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

class VariantRangeSummary(BaseModel):
    """Variants grouped by position in a range."""
    family: str
    range: Dict[str, int]  # {"start": x, "end": y}
    positions_with_variants: int
    data: Dict[int, List[Dict[str, Any]]]  # position -> list of variant summaries


class VariantStats(BaseModel):
    """Variant statistics for a family."""
    family: str
    by_type: Dict[str, int]
    position_range: Dict[str, Optional[int]]  # {"min": x, "max": y}
    total_variants: int


# Endpoint signature changes:


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
    rcsb_id       : Optional[str] = None       # None for free-standing literature variants
    entity_id     : Optional[str] = None       # None for free-standing literature variants
    # Literature-specific fields (Morisette database, etc.)
    species       : Optional[str] = None
    tubulin_type  : Optional[str] = None
    family        : Optional[str] = None
    reference_link: Optional[str] = None
    keywords      : Optional[str] = None
    notes         : Optional[str] = None
    utn_position  : Optional[int] = None


class PositionAnnotationsResponse(BaseModel):
    """Annotations at a specific position."""

    position: int
    family: str
    variants: List[VariantAnnotation]
    total_count: int


class ModificationAnnotation(BaseModel):
    """PTM annotation from the Morisette database."""

    master_index      : int
    amino_acid        : str
    modification_type : str   # acetylation, phosphorylation, etc.
    uniprot_id        : Optional[str] = None
    species           : Optional[str] = None
    tubulin_type      : Optional[str] = None
    family            : Optional[str] = None
    phenotype         : Optional[str] = None
    utn_position      : Optional[int] = None
    database_source   : Optional[str] = None
    database_link     : Optional[str] = None
    keywords          : Optional[str] = None
    notes             : Optional[str] = None


class ModificationsResponse(BaseModel):
    """Modifications at a position or for a family."""

    family: str
    modifications: List[ModificationAnnotation]
    total_count: int


class PolymerAnnotationsResponse(BaseModel):
    """All annotations for a polymer chain."""

    rcsb_id: str
    auth_asym_id: str
    entity_id: str
    family: Optional[str]
    variants: List[VariantAnnotation]
    modifications: List[ModificationAnnotation] = []
    total_count: int


# =============================================================================
# Query Functions
# =============================================================================


_VARIANT_RETURN_FIELDS = """
    type: v.type,
    master_index: v.master_index,
    observed_index: v.observed_index,
    wild_type: v.wild_type,
    observed: v.observed,
    source: v.source,
    uniprot_id: v.uniprot_id,
    phenotype: v.phenotype,
    reference: v.reference,
    species: v.species,
    tubulin_type: v.tubulin_type,
    family: v.family,
    reference_link: v.reference_link,
    keywords: v.keywords,
    notes: v.notes,
    utn_position: v.utn_position
"""


def get_variants_at_position(position: int, family: str) -> List[Dict[str, Any]]:
    """Get all variants at a specific master alignment position for a family.
    Unions entity-linked (structural) and free-standing (literature) variants."""
    query = f"""
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE v.master_index = $position AND e.family = $family
    RETURN {{
        {_VARIANT_RETURN_FIELDS},
        rcsb_id: e.parent_rcsb_id,
        entity_id: e.entity_id
    }} AS variant
    UNION
    MATCH (v:Variant)
    WHERE v.master_index = $position AND v.family = $family
      AND NOT EXISTS {{ MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }}
    RETURN {{
        {_VARIANT_RETURN_FIELDS},
        rcsb_id: null,
        entity_id: null
    }} AS variant
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
    """Get all variants for a specific polymer chain.
    Returns structural (entity-linked) + literature (free-standing by family) variants,
    plus entity_id and family."""

    # Step 1: Get entity-linked variants and the entity metadata
    structural_query = f"""
    MATCH (i:PolypeptideInstance {{parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id}})
    MATCH (i)-[:INSTANCE_OF]->(e:PolypeptideEntity)
    OPTIONAL MATCH (e)-[:HAS_VARIANT]->(v:Variant)
    WITH e, collect(CASE WHEN v IS NOT NULL THEN {{
        {_VARIANT_RETURN_FIELDS},
        rcsb_id: e.parent_rcsb_id,
        entity_id: e.entity_id
    }} END) AS variants
    RETURN variants, e.entity_id AS entity_id, e.family AS family
    """

    # Step 2: Get literature variants by family (free-standing)
    literature_query = f"""
    MATCH (v:Variant)
    WHERE v.family = $family
      AND NOT EXISTS {{ MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }}
    RETURN {{
        {_VARIANT_RETURN_FIELDS},
        rcsb_id: null,
        entity_id: null
    }} AS variant
    """

    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            result = tx.run(
                structural_query, {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id}
            ).single()
            if not result:
                return [], None, None

            structural_variants = [v for v in result["variants"] if v is not None]
            entity_id = result["entity_id"]
            family = result["family"]

            # Only fetch literature variants if the entity belongs to a tubulin family
            literature_variants = []
            if family and family.startswith("tubulin_"):
                literature_variants = [
                    dict(r["variant"])
                    for r in tx.run(literature_query, {"family": family})
                ]

            return structural_variants + literature_variants, entity_id, family

        return session.execute_read(run)


# =============================================================================
# Position-based Endpoints
# =============================================================================


@router_annotations.get("/variants/{family}/{position}", response_model=PositionAnnotationsResponse, operation_id="get_variants_at_position")

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


@router_annotations.get("/range/{family}", response_model=VariantRangeSummary, operation_id="get_variants_in_range")
async def get_variants_in_range(
    family: str,
    start: int = Query(..., description="Start position (inclusive)"),
    end: int = Query(..., description="End position (inclusive)"),
) -> Dict[str, Any]:
    """
    Get all variants within a position range for a family.
    Unions entity-linked (structural) and free-standing (literature) variants.
    """
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE v.master_index >= $start AND v.master_index <= $end AND e.family = $family
    RETURN v.master_index AS position, {
        type: v.type,
        wild_type: v.wild_type,
        observed: v.observed,
        source: v.source,
        rcsb_id: e.parent_rcsb_id
    } AS variant
    UNION
    MATCH (v:Variant)
    WHERE v.master_index >= $start AND v.master_index <= $end AND v.family = $family
      AND NOT EXISTS { MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }
    RETURN v.master_index AS position, {
        type: v.type,
        wild_type: v.wild_type,
        observed: v.observed,
        source: v.source,
        rcsb_id: null
    } AS variant
    """

    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            variants_by_pos: Dict[int, List[Dict[str, Any]]] = {}
            for r in tx.run(query, {"start": start, "end": end, "family": family}):
                pos = r["position"]
                if pos not in variants_by_pos:
                    variants_by_pos[pos] = []
                variants_by_pos[pos].append(dict(r["variant"]))
            return variants_by_pos

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



@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}", response_model=PolymerAnnotationsResponse, operation_id="get_polymer_annotations")
async def get_polymer_annotations(
    rcsb_id: str, auth_asym_id: str
) -> PolymerAnnotationsResponse:
    """Get all variant and modification annotations for a specific polymer chain."""
    variants, entity_id, family = get_variants_for_polymer(rcsb_id, auth_asym_id)

    if entity_id is None:
        raise HTTPException(
            404, f"Polymer chain {auth_asym_id} not found in structure {rcsb_id}"
        )

    # Fetch modifications for this chain's family
    modifications = []
    if family and family.startswith("tubulin_"):
        modifications = get_modifications_for_family(family)

    return PolymerAnnotationsResponse(
        rcsb_id=rcsb_id.upper(),
        auth_asym_id=auth_asym_id,
        entity_id=entity_id,
        family=family,
        variants=[VariantAnnotation(**v) for v in variants],
        modifications=[ModificationAnnotation(**m) for m in modifications],
        total_count=len(variants),
    )


# =============================================================================
# Summary/Stats Endpoints
# =============================================================================


@router_annotations.get("/stats/{family}", response_model=VariantStats, operation_id="get_variant_stats")
async def get_variant_stats(family: str) -> Dict[str, Any]:
    """
    Get variant statistics for a tubulin family.
    Includes both entity-linked (structural) and free-standing (literature) variants.
    """
    # Entity-linked variants
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE e.family = $family
    RETURN v.type AS variant_type, v.source AS source, count(*) AS count
    UNION ALL
    MATCH (v:Variant)
    WHERE v.family = $family
      AND NOT EXISTS { MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }
    RETURN v.type AS variant_type, v.source AS source, count(*) AS count
    """

    position_query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE e.family = $family AND v.master_index IS NOT NULL
    RETURN min(v.master_index) AS min_pos, max(v.master_index) AS max_pos, count(*) AS total
    UNION ALL
    MATCH (v:Variant)
    WHERE v.family = $family AND v.master_index IS NOT NULL
      AND NOT EXISTS { MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }
    RETURN min(v.master_index) AS min_pos, max(v.master_index) AS max_pos, count(*) AS total
    """

    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            # Aggregate type counts across both sources
            type_counts: Dict[str, int] = {}
            for r in tx.run(query, {"family": family}):
                vtype = r["variant_type"]
                type_counts[vtype] = type_counts.get(vtype, 0) + r["count"]

            # Aggregate position stats across both sources
            overall_min = None
            overall_max = None
            total = 0
            for r in tx.run(position_query, {"family": family}):
                if r["min_pos"] is not None:
                    if overall_min is None or r["min_pos"] < overall_min:
                        overall_min = r["min_pos"]
                if r["max_pos"] is not None:
                    if overall_max is None or r["max_pos"] > overall_max:
                        overall_max = r["max_pos"]
                total += r["total"]

            return type_counts, overall_min, overall_max, total

        type_counts, min_pos, max_pos, total = session.execute_read(run)

    return {
        "family": family,
        "by_type": type_counts,
        "position_range": {
            "min": min_pos,
            "max": max_pos,
        },
        "total_variants": total,
    }


# =============================================================================
# Modification (PTM) Queries & Endpoints
# =============================================================================

_MODIFICATION_RETURN_FIELDS = """
    master_index: m.master_index,
    amino_acid: m.amino_acid,
    modification_type: m.modification_type,
    uniprot_id: m.uniprot_id,
    species: m.species,
    tubulin_type: m.tubulin_type,
    family: m.family,
    phenotype: m.phenotype,
    utn_position: m.utn_position,
    database_source: m.database_source,
    database_link: m.database_link,
    keywords: m.keywords,
    notes: m.notes
"""


def get_modifications_for_family(
    family: str,
    position_min: Optional[int] = None,
    position_max: Optional[int] = None,
    modification_type: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Get all modifications for a tubulin family, optionally filtered."""
    conditions = ["m.family = $family"]
    params: Dict[str, Any] = {"family": family}

    if position_min is not None:
        conditions.append("m.master_index >= $pos_min")
        params["pos_min"] = position_min
    if position_max is not None:
        conditions.append("m.master_index <= $pos_max")
        params["pos_max"] = position_max
    if modification_type:
        conditions.append("m.modification_type = $mod_type")
        params["mod_type"] = modification_type

    where_clause = " AND ".join(conditions)
    query = f"""
    MATCH (m:Modification)
    WHERE {where_clause}
    RETURN {{ {_MODIFICATION_RETURN_FIELDS} }} AS modification
    ORDER BY m.master_index
    """

    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return [dict(r["modification"]) for r in tx.run(query, params)]

        return session.execute_read(run)


@router_annotations.get("/modifications/{family}", response_model=ModificationsResponse, operation_id="get_modifications_for_family")
async def get_modifications_for_family_endpoint(
    family: str,
    position_min: Optional[int] = Query(None, alias="posMin"),
    position_max: Optional[int] = Query(None, alias="posMax"),
    modification_type: Optional[str] = Query(None, alias="modType"),
) -> ModificationsResponse:
    """Get all modifications (PTMs) for a tubulin family."""
    mods = get_modifications_for_family(family, position_min, position_max, modification_type)
    return ModificationsResponse(
        family=family,
        modifications=[ModificationAnnotation(**m) for m in mods],
        total_count=len(mods),
    )


@router_annotations.get("/modifications/{family}/{position}", response_model=ModificationsResponse, operation_id="get_modifications_at_position")
async def get_modifications_at_position_endpoint(
    family: str,
    position: int,
) -> ModificationsResponse:
    """Get all modifications at a specific master alignment position."""
    mods = get_modifications_for_family(family, position_min=position, position_max=position)
    return ModificationsResponse(
        family=family,
        modifications=[ModificationAnnotation(**m) for m in mods],
        total_count=len(mods),
    )
