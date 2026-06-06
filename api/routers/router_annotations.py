# api/routers/router_annotations.py
from fastapi import APIRouter, HTTPException, Query
from typing import List, Dict, Any, Optional, Literal, Tuple, Union, Annotated, Set
from pydantic import BaseModel, Field
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
    species           : Optional[str] = None      # source abbreviation, e.g. "H. sapiens"
    tax_id            : Optional[int] = None      # NCBI taxonomy id, e.g. 9606
    species_full_name : Optional[str] = None      # canonical scientific name
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
    chain_tax_id: Optional[int] = None              # NCBI taxonomy id of this chain's source organism (first entry of src_organism_ids)
    chain_species_full_name: Optional[str] = None   # canonical scientific name matching chain_tax_id
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
) -> tuple[List[Dict[str, Any]], Optional[str], Optional[str], Optional[int], Optional[str]]:
    """Get all variants for a specific polymer chain.
    Returns structural (entity-linked) + literature (free-standing by family) variants,
    plus entity_id, family, chain_tax_id and chain_species_full_name.

    chain_tax_id is **always coerced up to the 'species' rank** by walking the
    phylogeny graph (5..species)-[:descendant_of]->(508771..strain), so a strain
    taxid like 508771 (Toxoplasma gondii ME49) resolves to its species ancestor
    5811 (Toxoplasma gondii). Necessary because the Morisette PTM data is stored
    at species level -- without this coercion, structures whose PDB taxid is
    sub-species (strain/isolate/subspecies) get zero PTM matches even when their
    species has data. Falls back to the raw first src_organism_id if no species
    ancestor exists in the phylogeny graph."""

    # Step 1: Get entity-linked variants, entity metadata, AND coerce the
    # entity's first source organism to species rank.
    structural_query = f"""
    MATCH (i:PolypeptideInstance {{parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id}})
    MATCH (i)-[:INSTANCE_OF]->(e:PolypeptideEntity)
    OPTIONAL MATCH (e)-[:HAS_VARIANT]->(v:Variant)
    WITH e, collect(CASE WHEN v IS NOT NULL THEN {{
        {_VARIANT_RETURN_FIELDS},
        rcsb_id: e.parent_rcsb_id,
        entity_id: e.entity_id
    }} END) AS variants
    CALL {{
        WITH e
        OPTIONAL MATCH (start:PhylogenyNode {{ncbi_tax_id: e.src_organism_ids[0]}})
                       <-[:descendant_of*0..]-(sp:PhylogenyNode {{rank: 'species'}})
        RETURN sp.ncbi_tax_id     AS species_tax_id,
               sp.scientific_name AS species_full_name
        LIMIT 1
    }}
    RETURN variants,
           e.entity_id AS entity_id,
           e.family AS family,
           e.src_organism_ids AS src_organism_ids,
           e.src_organism_names AS src_organism_names,
           species_tax_id,
           species_full_name
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
                return [], None, None, None, None

            structural_variants = [v for v in result["variants"] if v is not None]
            entity_id = result["entity_id"]
            family = result["family"]

            org_ids = result["src_organism_ids"] or []
            org_names = result["src_organism_names"] or []
            # Prefer the species-rank ancestor (so strains resolve to their species).
            # Fall back to the raw first src_organism_id if the phylogeny graph
            # doesn't contain that taxid yet (rare; new strains not seeded).
            sp_id = result["species_tax_id"]
            sp_name = result["species_full_name"]
            chain_tax_id = (
                int(sp_id) if sp_id is not None
                else (int(org_ids[0]) if org_ids else None)
            )
            chain_species_full_name = (
                sp_name if sp_name is not None
                else (org_names[0] if org_names else None)
            )

            # Only fetch literature variants if the entity belongs to a tubulin family
            literature_variants = []
            if family and family.startswith("tubulin_"):
                literature_variants = [
                    dict(r["variant"])
                    for r in tx.run(literature_query, {"family": family})
                ]

            return (
                structural_variants + literature_variants,
                entity_id,
                family,
                chain_tax_id,
                chain_species_full_name,
            )

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
    variants, entity_id, family, chain_tax_id, chain_species_full_name = (
        get_variants_for_polymer(rcsb_id, auth_asym_id)
    )

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
        chain_tax_id=chain_tax_id,
        chain_species_full_name=chain_species_full_name,
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
    tax_id: m.tax_id,
    species_full_name: m.species_full_name,
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


# =============================================================================
# Annotation Tracks (filter -> master positions)
# =============================================================================
#
# A "track" is a chain-independent aux row whose paint targets are the set of
# master columns matching a typed FilterSpec. The same spec is consumed by:
#   - the frontend Add Track modal (preview + creation)
#   - the LLM viewer action AddAnnotationTrack (no position invention)
#   - URL deep-links (Phase 6)
#
# Composition: AND across fields, implicit OR within array-valued fields.
# No cross-field boolean OR in v1.
#
# Deferred to follow-up: aa_property_change, conservation_range,
# variability_range, min_substitution_count, chemical_classes, heatmap counts.


VariantSourceType = Literal['structural', 'literature']


class VariantFilterSpec(BaseModel):
    """Filter spec for variant-based annotation tracks (structural + literature)."""

    kind: Literal['variants'] = 'variants'
    family: str = Field(..., description="e.g. 'tubulin_alpha', 'tubulin_beta'")

    # --- AA-level (typed) ---
    wild_type_aas: Optional[List[str]] = Field(None, description="Wild-type AA in set (one-letter codes)")
    observed_aas: Optional[List[str]] = Field(None, description="Observed AA in set (one-letter codes)")
    substitution_pairs: Optional[List[Tuple[str, str]]] = Field(
        None, description="Specific (wt, obs) pairs; matches any pair in the list"
    )
    indel_present: Optional[bool] = Field(
        None,
        description="True: only insertion/deletion records. False: only substitutions. None: any.",
    )

    # --- Provenance / scoping (typed) ---
    sources: Optional[List[VariantSourceType]] = Field(
        None, description="Variant provenance; defaults to both. 'literature' maps to morisette internally."
    )
    uniprot_ids: Optional[List[str]] = Field(None, description="Restrict to specific UniProt ids (e.g. ['Q13509'] for TUBB3)")
    species_names: Optional[List[str]] = Field(
        None,
        description="Match variant.species text (e.g. 'H. sapiens'). Use species_tax_ids for structural variants.",
    )
    species_tax_ids: Optional[List[int]] = Field(
        None,
        description=(
            "NCBI tax ids; filters structural variants by the parent entity's src_organism_ids "
            "(e.g. [9606] for human). When set, literature variants are dropped from the result "
            "since they don't carry entity/tax info — use species_names for literature."
        ),
    )
    position_range: Optional[Tuple[int, int]] = Field(None, description="Master column [min, max], inclusive")
    positions: Optional[List[int]] = Field(
        None,
        description=(
            "Restrict to an explicit set of master indices. When set together with position_range, "
            "both are ANDed. Useful for site-scoped queries (e.g. residues of a binding site)."
        ),
    )
    co_occurs_with_mod_type: Optional[List[str]] = Field(
        None, description="Keep only positions where any of these PTM types is also recorded"
    )

    # --- Text overlay (secondary, opt-in) ---
    phenotype_contains: Optional[List[str]] = Field(
        None,
        description="Phenotype substring search (case-insensitive); OR within array. Only useful when sources includes 'literature'.",
    )


class ModificationFilterSpec(BaseModel):
    """Filter spec for modification (PTM) tracks."""

    kind: Literal['modifications'] = 'modifications'
    family: str

    modification_types: Optional[List[str]] = Field(None, description="e.g. ['phosphorylation', 'acetylation']")
    uniprot_ids: Optional[List[str]] = None
    species_tax_ids: Optional[List[int]] = Field(None, description="NCBI taxIds; modifications carry tax_id directly")
    species_names: Optional[List[str]] = Field(None, description="Fallback text match against species abbreviation")
    position_range: Optional[Tuple[int, int]] = None
    positions: Optional[List[int]] = Field(
        None,
        description=(
            "Restrict to an explicit set of master indices. ANDed with position_range when both set."
        ),
    )
    co_occurs_with_variant: Optional[bool] = Field(
        None, description="Keep only positions where any variant is also recorded"
    )
    evidence_source: Optional[List[str]] = Field(
        None, description="Morisette database_source column (e.g. 'SwissPalm', 'PhosphoSitePlus')"
    )
    phenotype_contains: Optional[List[str]] = None


class BindingContactFilterSpec(BaseModel):
    """Filter spec for ligand binding-contact tracks (chemical_ids only in v1)."""

    kind: Literal['binding_contacts'] = 'binding_contacts'
    family: str

    chemical_ids: List[str] = Field(..., description="Ligand chem comp ids (e.g. ['TA1', 'GTP', 'TXL'])")
    structure_ids: Optional[List[str]] = Field(None, description="Restrict to specific PDB ids")
    positions: Optional[List[int]] = Field(
        None,
        description=(
            "Restrict to an explicit set of master indices. Only contacts at these positions "
            "are returned (post-filter on parsed residues_json)."
        ),
    )


TrackFilterSpec = Annotated[
    Union[VariantFilterSpec, ModificationFilterSpec, BindingContactFilterSpec],
    Field(discriminator='kind'),
]


class ResolvedPosition(BaseModel):
    master_index: int
    match_count: int
    matched_records: List[Dict[str, Any]]


class ResolveResponse(BaseModel):
    family: str
    kind: str
    matched: int
    positions: List[ResolvedPosition]


class ResolveRequest(BaseModel):
    """Body for POST /annotations/track/resolve.

    The spec field is a discriminated union — Pydantic picks the variant by 'kind'.
    """
    spec: TrackFilterSpec


# -----------------------------------------------------------------------------
# Variant resolve
# -----------------------------------------------------------------------------

# Frontend-facing source label -> backend storage value
_VARIANT_SOURCE_MAP = {
    'structural': 'structural',
    'literature': 'morisette',
}


def _build_variant_where(spec: VariantFilterSpec) -> Tuple[str, Dict[str, Any]]:
    """Build the WHERE clause for a Variant match, plus the Cypher params.
    Returns ('TRUE' if no conditions, else 'cond1 AND cond2 ...')."""
    conditions: List[str] = ["v.master_index IS NOT NULL"]
    params: Dict[str, Any] = {}

    if spec.sources:
        backend_sources = [_VARIANT_SOURCE_MAP[s] for s in spec.sources]
        conditions.append("v.source IN $sources")
        params['sources'] = backend_sources

    if spec.wild_type_aas:
        conditions.append("v.wild_type IN $wt_aas")
        params['wt_aas'] = [a.upper() for a in spec.wild_type_aas]
    if spec.observed_aas:
        conditions.append("v.observed IN $obs_aas")
        params['obs_aas'] = [a.upper() for a in spec.observed_aas]
    if spec.substitution_pairs:
        pairs = [{"wt": wt.upper(), "obs": obs.upper()} for wt, obs in spec.substitution_pairs]
        conditions.append(
            "ANY(p IN $sub_pairs WHERE v.wild_type = p.wt AND v.observed = p.obs)"
        )
        params['sub_pairs'] = pairs
    if spec.indel_present is True:
        conditions.append("v.type IN ['insertion', 'deletion']")
    elif spec.indel_present is False:
        conditions.append("v.type = 'substitution'")

    if spec.uniprot_ids:
        conditions.append("v.uniprot_id IN $uniprot_ids")
        params['uniprot_ids'] = spec.uniprot_ids
    if spec.species_names:
        conditions.append("v.species IN $species_names")
        params['species_names'] = spec.species_names
    if spec.position_range:
        conditions.append("v.master_index >= $pos_min AND v.master_index <= $pos_max")
        params['pos_min'] = spec.position_range[0]
        params['pos_max'] = spec.position_range[1]
    if spec.positions:
        conditions.append("v.master_index IN $positions")
        params['positions'] = spec.positions

    if spec.phenotype_contains:
        conditions.append(
            "v.phenotype IS NOT NULL AND ANY(p IN $phen_patterns WHERE toLower(v.phenotype) CONTAINS toLower(p))"
        )
        params['phen_patterns'] = spec.phenotype_contains

    return " AND ".join(conditions), params


def _positions_with_mod_types(family: str, mod_types: List[str]) -> Set[int]:
    """Master positions where at least one Modification of the given types exists."""
    query = """
    MATCH (m:Modification)
    WHERE m.family = $family AND m.modification_type IN $mod_types
      AND m.master_index IS NOT NULL
    RETURN DISTINCT m.master_index AS pos
    """
    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            return {r["pos"] for r in tx.run(query, {"family": family, "mod_types": mod_types})}
        return session.execute_read(run)


def _positions_with_any_variant(family: str) -> Set[int]:
    """Master positions where at least one Variant exists for the family."""
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE e.family = $family AND v.master_index IS NOT NULL
    RETURN DISTINCT v.master_index AS pos
    UNION
    MATCH (v:Variant)
    WHERE v.family = $family AND v.master_index IS NOT NULL
      AND NOT EXISTS { MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }
    RETURN DISTINCT v.master_index AS pos
    """
    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            return {r["pos"] for r in tx.run(query, {"family": family})}
        return session.execute_read(run)


def resolve_variant_track(spec: VariantFilterSpec) -> List[Dict[str, Any]]:
    """Resolve a variant FilterSpec to master positions + their matched records.

    Two-arm query unions entity-linked (structural-typically) and free-standing
    (literature-typically) Variant nodes. Co-occurrence with PTM types is applied
    as a post-filter via _positions_with_mod_types.
    """
    where_clause, params = _build_variant_where(spec)
    params["family"] = spec.family

    # species_tax_ids is an entity-level filter (Variant nodes don't carry tax_id),
    # so it's injected into the structural arm only and drops the literature arm.
    structural_extra = ""
    if spec.species_tax_ids:
        structural_extra = " AND ANY(t IN $tax_ids WHERE t IN e.src_organism_ids)"
        params["tax_ids"] = spec.species_tax_ids

    structural_arm = f"""
    MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
    WHERE e.family = $family AND {where_clause}{structural_extra}
    RETURN v.master_index AS master_index, {{
        wild_type: v.wild_type,
        observed: v.observed,
        type: v.type,
        source: v.source,
        uniprot_id: v.uniprot_id,
        phenotype: v.phenotype,
        species: v.species,
        utn_position: v.utn_position,
        rcsb_id: e.parent_rcsb_id,
        entity_id: e.entity_id
    }} AS record
    """

    literature_arm = f"""
    MATCH (v:Variant)
    WHERE v.family = $family
      AND NOT EXISTS {{ MATCH (:PolypeptideEntity)-[:HAS_VARIANT]->(v) }}
      AND {where_clause}
    RETURN v.master_index AS master_index, {{
        wild_type: v.wild_type,
        observed: v.observed,
        type: v.type,
        source: v.source,
        uniprot_id: v.uniprot_id,
        phenotype: v.phenotype,
        species: v.species,
        utn_position: v.utn_position,
        rcsb_id: null,
        entity_id: null
    }} AS record
    """

    query = structural_arm if spec.species_tax_ids else f"{structural_arm}\nUNION ALL\n{literature_arm}"

    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            by_pos: Dict[int, List[Dict[str, Any]]] = {}
            for r in tx.run(query, params):
                pos = r["master_index"]
                by_pos.setdefault(pos, []).append(dict(r["record"]))
            return by_pos

        by_pos = session.execute_read(run)

    # Post-filter: co-occurrence with PTM types
    if spec.co_occurs_with_mod_type:
        mod_positions = _positions_with_mod_types(spec.family, spec.co_occurs_with_mod_type)
        by_pos = {p: recs for p, recs in by_pos.items() if p in mod_positions}

    return [
        {"master_index": p, "match_count": len(recs), "matched_records": recs}
        for p, recs in sorted(by_pos.items())
    ]


# -----------------------------------------------------------------------------
# Modification resolve
# -----------------------------------------------------------------------------


def _build_modification_where(spec: ModificationFilterSpec) -> Tuple[str, Dict[str, Any]]:
    conditions: List[str] = ["m.master_index IS NOT NULL"]
    params: Dict[str, Any] = {}

    if spec.modification_types:
        conditions.append("m.modification_type IN $mod_types")
        params['mod_types'] = spec.modification_types
    if spec.uniprot_ids:
        conditions.append("m.uniprot_id IN $uniprot_ids")
        params['uniprot_ids'] = spec.uniprot_ids
    if spec.species_tax_ids:
        conditions.append("m.tax_id IN $tax_ids")
        params['tax_ids'] = spec.species_tax_ids
    if spec.species_names:
        conditions.append("m.species IN $species_names")
        params['species_names'] = spec.species_names
    if spec.position_range:
        conditions.append("m.master_index >= $pos_min AND m.master_index <= $pos_max")
        params['pos_min'] = spec.position_range[0]
        params['pos_max'] = spec.position_range[1]
    if spec.positions:
        conditions.append("m.master_index IN $positions")
        params['positions'] = spec.positions
    if spec.evidence_source:
        conditions.append("m.database_source IN $ev_sources")
        params['ev_sources'] = spec.evidence_source
    if spec.phenotype_contains:
        conditions.append(
            "m.phenotype IS NOT NULL AND ANY(p IN $phen_patterns WHERE toLower(m.phenotype) CONTAINS toLower(p))"
        )
        params['phen_patterns'] = spec.phenotype_contains

    return " AND ".join(conditions), params


def resolve_modification_track(spec: ModificationFilterSpec) -> List[Dict[str, Any]]:
    where_clause, params = _build_modification_where(spec)
    params["family"] = spec.family

    query = f"""
    MATCH (m:Modification)
    WHERE m.family = $family AND {where_clause}
    RETURN m.master_index AS master_index, {{
        modification_type: m.modification_type,
        amino_acid: m.amino_acid,
        uniprot_id: m.uniprot_id,
        species: m.species,
        tax_id: m.tax_id,
        species_full_name: m.species_full_name,
        phenotype: m.phenotype,
        database_source: m.database_source,
        utn_position: m.utn_position
    }} AS record
    """

    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            by_pos: Dict[int, List[Dict[str, Any]]] = {}
            for r in tx.run(query, params):
                pos = r["master_index"]
                by_pos.setdefault(pos, []).append(dict(r["record"]))
            return by_pos

        by_pos = session.execute_read(run)

    if spec.co_occurs_with_variant:
        var_positions = _positions_with_any_variant(spec.family)
        by_pos = {p: recs for p, recs in by_pos.items() if p in var_positions}

    return [
        {"master_index": p, "match_count": len(recs), "matched_records": recs}
        for p, recs in sorted(by_pos.items())
    ]


# -----------------------------------------------------------------------------
# Binding contact resolve
# -----------------------------------------------------------------------------


def resolve_binding_contact_track(spec: BindingContactFilterSpec) -> List[Dict[str, Any]]:
    """Resolve ligand-contact positions across structures for a set of chemical ids.

    Mirrors get_canonical_binding_site: NonpolymerInstance -[INSTANCE_OF]-> NonpolymerEntity
    -[DEFINED_BY_CHEMICAL]-> Chemical, then NonpolymerInstance -[NEAR_POLYMER]-> PolypeptideInstance.
    Contact master indices are parsed from r.residues_json.
    """
    import json as _json

    structure_clause = ""
    params: Dict[str, Any] = {
        "family": spec.family,
        "chemical_ids": [c.upper() for c in spec.chemical_ids],
    }
    if spec.structure_ids:
        structure_clause = "AND li.parent_rcsb_id IN $structure_ids"
        params['structure_ids'] = [s.upper() for s in spec.structure_ids]

    positions_filter: Optional[Set[int]] = (
        set(spec.positions) if spec.positions else None
    )

    query = f"""
    MATCH (li:NonpolymerInstance)-[:INSTANCE_OF]->(ne:NonpolymerEntity)
          -[:DEFINED_BY_CHEMICAL]->(c:Chemical)
    WHERE c.chemical_id IN $chemical_ids
    MATCH (li)-[r:NEAR_POLYMER]->(pi:PolypeptideInstance)
          -[:INSTANCE_OF]->(pe:PolypeptideEntity)
    WHERE pe.family = $family
      {structure_clause}
    RETURN li.parent_rcsb_id AS rcsb_id,
           c.chemical_id AS chemical_id,
           pi.auth_asym_id AS auth_asym_id,
           r.residues_json AS residues_json
    """

    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            by_pos: Dict[int, List[Dict[str, Any]]] = {}
            for row in tx.run(query, params):
                residues_json = row["residues_json"]
                if not residues_json:
                    continue
                try:
                    residues = _json.loads(residues_json)
                except (ValueError, TypeError):
                    continue
                # Dedup master indices within this (ligand, chain) — multiple
                # contact atoms can hit the same residue.
                local_indices: Set[int] = set()
                for res in residues:
                    mi = res.get("master_index") if isinstance(res, dict) else None
                    if mi is not None:
                        local_indices.add(int(mi))
                if positions_filter is not None:
                    local_indices &= positions_filter
                for mi in local_indices:
                    by_pos.setdefault(mi, []).append({
                        "chemical_id": row["chemical_id"],
                        "rcsb_id": row["rcsb_id"],
                        "auth_asym_id": row["auth_asym_id"],
                    })
            return by_pos

        by_pos = session.execute_read(run)

    return [
        {"master_index": p, "match_count": len(recs), "matched_records": recs}
        for p, recs in sorted(by_pos.items())
    ]


# -----------------------------------------------------------------------------
# Resolve endpoint
# -----------------------------------------------------------------------------


@router_annotations.post(
    "/track/resolve",
    response_model=ResolveResponse,
    operation_id="resolve_annotation_track",
)
async def resolve_annotation_track(req: ResolveRequest) -> ResolveResponse:
    """Resolve a track FilterSpec to the set of master columns it paints.

    The spec is discriminated by `kind`:
      - 'variants': filter Variant nodes (structural + literature)
      - 'modifications': filter Modification nodes
      - 'binding_contacts': aggregate NEAR_POLYMER contacts for given chemical ids

    Returns master positions with their matched records as metadata payload.
    Composition: AND across fields, implicit OR within array-valued fields.
    """
    spec = req.spec

    if isinstance(spec, VariantFilterSpec):
        positions = resolve_variant_track(spec)
    elif isinstance(spec, ModificationFilterSpec):
        positions = resolve_modification_track(spec)
    elif isinstance(spec, BindingContactFilterSpec):
        positions = resolve_binding_contact_track(spec)
    else:
        raise HTTPException(400, f"Unknown FilterSpec kind: {getattr(spec, 'kind', '?')}")

    return ResolveResponse(
        family=spec.family,
        kind=spec.kind,
        matched=len(positions),
        positions=[ResolvedPosition(**p) for p in positions],
    )
