"""Read-only retrieval tools for the assistant orchestrator.

Every tool wraps an existing, vetted, parameterized query function over the
shared `db_reader`. The model SELECTS a tool and supplies typed arguments — it
never supplies Cypher. This is the hard guarantee: all DB access is through the
fixed methods below.

Returns are compact, JSON-serializable dicts (counts + small samples, not full
row dumps) so the orchestration loop's token budget stays bounded.

Argument models are kept FLAT (no deeply-nested unions) because that is what the
models call reliably; the richer typed FilterSpec objects are constructed
internally from the flat args before delegating to the existing query layer.
"""
from __future__ import annotations

from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Tuple, Type

from pydantic import BaseModel, Field

from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import LigandFilters, StructureFilters

from api.nl_translator.facets_loader import load_facet_context
from api.nl_translator.resolve import resolve_representative
from api.routers.router_annotations import (
    BindingContactFilterSpec,
    ModificationFilterSpec,
    VariantFilterSpec,
    resolve_binding_contact_track,
    resolve_modification_track,
    resolve_variant_track,
)

# Max rows we ever put in a `sample`/`positions` list returned to the model.
_SAMPLE = 8
_MAX_POSITIONS = 60


# ---------------------------------------------------------------------------
# Argument models (flat, minimal — only fields the model needs to express)
# ---------------------------------------------------------------------------

class FindStructuresArgs(BaseModel):
    """Filter the PDB structure catalogue. All fields optional; AND-combined."""
    search: Optional[str] = Field(None, description="Free-text search over title/keywords")
    has_ligand_ids: Optional[List[str]] = Field(None, description="PDB chem comp ids that must be bound, e.g. ['TA1','GTP']")
    has_polymer_family: Optional[List[str]] = Field(None, description="e.g. ['tubulin_alpha','tubulin_beta']")
    has_isotype: Optional[List[str]] = Field(None, description="Isotype codes, e.g. ['TUBB3']")
    source_organism_ids: Optional[List[int]] = Field(None, description="NCBI tax ids, e.g. [9606]")
    exp_method: Optional[List[str]] = Field(None, description="e.g. ['ELECTRON MICROSCOPY','X-RAY DIFFRACTION']")
    resolution_max: Optional[float] = Field(None, description="Best (highest) acceptable resolution in Angstroms; LOWER is better")
    year_min: Optional[int] = None
    has_variants: Optional[bool] = None


class GetStructureChainsArgs(BaseModel):
    """Per-chain map for one structure: which chain is which family and which
    ligands actually contact it. Use BEFORE focusing/highlighting a ligand to
    learn which chain it sits on."""
    rcsb_id: str = Field(..., description="PDB id, e.g. '1JFF'")


class ResolveStructureArgs(BaseModel):
    """Resolve intent (organism + family [+ ligand]) to a REAL (rcsb_id, chain)
    from the DB. Use instead of guessing a PDB id."""
    organism_id: Optional[int] = Field(None, description="NCBI tax id, e.g. 9606 human, 5811 Toxoplasma")
    family: Optional[str] = Field(None, description="e.g. 'tubulin_alpha'")
    ligand: Optional[str] = Field(None, description="Optional chem comp id that must bind, e.g. 'GTP'")


class GetBindingSiteArgs(BaseModel):
    """Canonical (cross-structure) binding-site residues for a ligand on a
    family, as master-alignment positions with frequency. This is the GROUNDED
    source of binding residues — never recall them from memory."""
    chemical_id: str = Field(..., description="PDB chem comp id, e.g. 'TA1', 'GTP', 'LOC'")
    family: str = Field(..., description="e.g. 'tubulin_beta'")


class CountModificationsArgs(BaseModel):
    """Count / summarize PTMs (modifications) for a family, with optional
    scoping. Returns distinct positions, total records, and a by-type breakdown."""
    family: str
    modification_types: Optional[List[str]] = Field(None, description="e.g. ['phosphorylation','acetylation','polyglutamylation']")
    species_tax_ids: Optional[List[int]] = Field(None, description="NCBI tax ids; modifications carry tax_id directly")
    position_range: Optional[Tuple[int, int]] = Field(None, description="[min,max] master positions, inclusive")
    positions: Optional[List[int]] = Field(None, description="Explicit master positions (e.g. a binding site)")
    phenotype_contains: Optional[List[str]] = Field(None, description="Case-insensitive phenotype substrings")


class CountVariantsArgs(BaseModel):
    """Count / summarize sequence variants for a family, with optional scoping."""
    family: str
    sources: Optional[List[str]] = Field(None, description="['structural','literature']; default both")
    species_tax_ids: Optional[List[int]] = Field(None, description="NCBI tax ids (structural variants only)")
    species_names: Optional[List[str]] = Field(None, description="Text species match for literature variants, e.g. ['H. sapiens']")
    position_range: Optional[Tuple[int, int]] = None
    positions: Optional[List[int]] = Field(None, description="Explicit master positions (e.g. a binding site); ANDed with position_range. Use to scope a count to a site.")
    wild_type_aas: Optional[List[str]] = None
    observed_aas: Optional[List[str]] = None
    phenotype_contains: Optional[List[str]] = None


class GetBindingContactsArgs(BaseModel):
    """Resolve ligand-contact master positions across a family (DB-wide), for a
    set of chemical ids. Use for "where does X bind on <family>" when no single
    loaded structure is in view."""
    family: str
    chemical_ids: List[str] = Field(..., description="e.g. ['TA1','GTP']")
    structure_ids: Optional[List[str]] = Field(None, description="Restrict to specific PDB ids")


class GetFacetsArgs(BaseModel):
    """Return the valid filter vocabulary (families, isotypes, top ligands,
    common organisms, ranges). Use to disambiguate names before filtering."""
    pass


# ---------------------------------------------------------------------------
# Tool implementations — each delegates to an existing query function.
# ---------------------------------------------------------------------------

def _summarize_positions(recs: List[Dict[str, Any]], type_key: Optional[str]) -> Dict[str, Any]:
    """Collapse a resolve_*_track result into a compact summary."""
    positions = [r["master_index"] for r in recs]
    total_records = sum(r.get("match_count", 0) for r in recs)
    by_type: Dict[str, int] = {}
    if type_key:
        counter: Counter = Counter()
        for r in recs:
            for rec in r.get("matched_records", []):
                v = rec.get(type_key)
                if v:
                    counter[v] += 1
        by_type = dict(counter.most_common())
    out: Dict[str, Any] = {
        "distinct_positions": len(positions),
        "total_records": total_records,
        "positions": sorted(positions)[:_MAX_POSITIONS],
        "truncated_positions": len(positions) > _MAX_POSITIONS,
    }
    if positions:
        out["position_range"] = [min(positions), max(positions)]
    if by_type:
        out["by_type"] = by_type
    return out


def find_structures(a: FindStructuresArgs) -> Dict[str, Any]:
    filters = StructureFilters(
        limit=_SAMPLE,
        search=a.search,
        has_ligand_ids=a.has_ligand_ids,
        has_polymer_family=a.has_polymer_family,
        has_isotype=a.has_isotype,
        source_organism_ids=a.source_organism_ids,
        exp_method=a.exp_method,
        resolution_max=a.resolution_max,
        year_min=a.year_min,
        has_variants=a.has_variants,
    )
    resp = db_reader.list_structures(filters)
    return {
        "total_count": resp.total_count,
        "has_more": resp.has_more,
        "sample": [
            {
                "rcsb_id": s.rcsb_id,
                "resolution": s.resolution,
                "exp_method": s.exp_method,
                "year": s.citation_year,
                "organisms": s.src_organism_names,
                "ligand_ids": s.ligand_ids,
                "families": s.polymer_families,
            }
            for s in resp.data[:_SAMPLE]
        ],
    }


def get_structure_chains(a: GetStructureChainsArgs) -> Dict[str, Any]:
    """Per-chain family + bound-ligand map. Fixed Cypher; no model input beyond rcsb_id."""
    query = """
    MATCH (s:Structure {rcsb_id: $rcsb})-[:HAS_INSTANCE]->(pi:PolypeptideInstance)
          -[:INSTANCE_OF]->(pe:PolypeptideEntity)
    OPTIONAL MATCH (ni:NonpolymerInstance)-[:NEAR_POLYMER]->(pi)
    OPTIONAL MATCH (ni)-[:INSTANCE_OF]->(:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
    WITH pi.auth_asym_id AS chain,
         pe.family AS family,
         pe.isotype AS isotype,
         pe.src_organism_ids AS tax_ids,
         collect(DISTINCT c.chemical_id) AS ligands
    RETURN chain, family, isotype, tax_ids,
           [l IN ligands WHERE l IS NOT NULL] AS bound_ligand_chem_ids
    ORDER BY chain
    """
    with db_reader.adapter.driver.session() as session:
        def run(tx):
            rows = []
            for r in tx.run(query, {"rcsb": a.rcsb_id.upper()}):
                tax = r["tax_ids"] or []
                rows.append({
                    "auth_asym_id": r["chain"],
                    "family": r["family"],
                    "isotype": r["isotype"],
                    "tax_id": tax[0] if tax else None,
                    "bound_ligand_chem_ids": list(r["bound_ligand_chem_ids"]),
                })
            return rows
        chains = session.execute_read(run)
    return {"rcsb_id": a.rcsb_id.upper(), "chains": chains}


def resolve_structure(a: ResolveStructureArgs) -> Dict[str, Any]:
    rep = resolve_representative(organism_id=a.organism_id, family=a.family, ligand=a.ligand)
    if rep is None:
        return {"found": False}
    return {"found": True, "rcsb_id": rep[0], "auth_asym_id": rep[1]}


def get_binding_site(a: GetBindingSiteArgs) -> Dict[str, Any]:
    site = db_reader.get_canonical_binding_site(a.chemical_id, a.family)
    if site is None:
        return {"found": False, "chemical_id": a.chemical_id.upper(), "family": a.family}
    residues = sorted(site.residues, key=lambda r: r.frequency, reverse=True)
    return {
        "found": True,
        "chemical_id": site.chemical_id,
        "chemical_name": site.chemical_name,
        "family": site.family,
        "structure_count": site.structure_count,
        "positions": sorted(r.master_index for r in site.residues)[:_MAX_POSITIONS],
        "top_positions": [
            {"master_index": r.master_index, "frequency": round(r.frequency, 3)}
            for r in residues[:_SAMPLE]
        ],
    }


def count_modifications(a: CountModificationsArgs) -> Dict[str, Any]:
    spec = ModificationFilterSpec(
        family=a.family,
        modification_types=a.modification_types,
        species_tax_ids=a.species_tax_ids,
        position_range=a.position_range,
        positions=a.positions,
        phenotype_contains=a.phenotype_contains,
    )
    recs = resolve_modification_track(spec)
    out = _summarize_positions(recs, type_key="modification_type")
    out["family"] = a.family
    return out


def count_variants(a: CountVariantsArgs) -> Dict[str, Any]:
    spec = VariantFilterSpec(
        family=a.family,
        sources=a.sources,
        species_tax_ids=a.species_tax_ids,
        species_names=a.species_names,
        position_range=a.position_range,
        positions=a.positions,
        wild_type_aas=a.wild_type_aas,
        observed_aas=a.observed_aas,
        phenotype_contains=a.phenotype_contains,
    )
    recs = resolve_variant_track(spec)
    out = _summarize_positions(recs, type_key="type")
    out["family"] = a.family
    return out


def get_binding_contacts(a: GetBindingContactsArgs) -> Dict[str, Any]:
    spec = BindingContactFilterSpec(
        family=a.family,
        chemical_ids=a.chemical_ids,
        structure_ids=a.structure_ids,
    )
    recs = resolve_binding_contact_track(spec)
    out = _summarize_positions(recs, type_key="chemical_id")
    out["family"] = a.family
    out["chemical_ids"] = [c.upper() for c in a.chemical_ids]
    return out


def get_facets(a: GetFacetsArgs) -> Dict[str, Any]:
    f = load_facet_context()
    return {
        "tubulin_families": f.tubulin_families,
        "isotypes": f.isotypes,
        "exp_methods": f.exp_methods,
        "top_ligands": f.top_ligands,
        "common_source_organisms": f.common_source_organisms,
        "year_range": f.year_range,
        "resolution_range": f.resolution_range,
    }


# ---------------------------------------------------------------------------
# Registry — drives tool-spec generation and dispatch in the orchestrator.
# ---------------------------------------------------------------------------

class RetrievalTool(BaseModel):
    name: str
    description: str
    args_model: Type[BaseModel]
    fn: Callable[[BaseModel], Dict[str, Any]]

    model_config = {"arbitrary_types_allowed": True}


RETRIEVAL_TOOLS: List[RetrievalTool] = [
    RetrievalTool(
        name="find_structures",
        description="Filter the PDB structure catalogue; returns total count + a small sample. Use for 'how many / which structures have X'.",
        args_model=FindStructuresArgs, fn=lambda x: find_structures(x),
    ),
    RetrievalTool(
        name="get_structure_chains",
        description="Per-chain family + bound-ligand map for one structure. Call BEFORE focusing a ligand to find which chain it actually sits on.",
        args_model=GetStructureChainsArgs, fn=lambda x: get_structure_chains(x),
    ),
    RetrievalTool(
        name="resolve_structure",
        description="Resolve organism+family[+ligand] intent to a real (rcsb_id, chain). Use instead of guessing PDB ids.",
        args_model=ResolveStructureArgs, fn=lambda x: resolve_structure(x),
    ),
    RetrievalTool(
        name="get_binding_site",
        description="Canonical cross-structure binding residues (master positions + frequency) for a ligand on a family. The grounded source of binding residues.",
        args_model=GetBindingSiteArgs, fn=lambda x: get_binding_site(x),
    ),
    RetrievalTool(
        name="count_modifications",
        description="Count/summarize PTMs for a family with optional scoping (types, species, positions, phenotype). Answers 'how many PTMs in ...'.",
        args_model=CountModificationsArgs, fn=lambda x: count_modifications(x),
    ),
    RetrievalTool(
        name="count_variants",
        description="Count/summarize sequence variants for a family with optional scoping.",
        args_model=CountVariantsArgs, fn=lambda x: count_variants(x),
    ),
    RetrievalTool(
        name="get_binding_contacts",
        description="DB-wide ligand-contact master positions for chemical ids on a family. Use for 'where does X bind on <family>' without a loaded structure.",
        args_model=GetBindingContactsArgs, fn=lambda x: get_binding_contacts(x),
    ),
    RetrievalTool(
        name="get_facets",
        description="Valid filter vocabulary (families, isotypes, top ligands, common organisms, ranges). Use to disambiguate names.",
        args_model=GetFacetsArgs, fn=lambda x: get_facets(x),
    ),
]

RETRIEVAL_TOOLS_BY_NAME: Dict[str, RetrievalTool] = {t.name: t for t in RETRIEVAL_TOOLS}


def run_retrieval_tool(name: str, raw_args: Dict[str, Any]) -> Dict[str, Any]:
    """Validate args against the tool's model and execute. Raises KeyError for
    an unknown tool; pydantic ValidationError for bad args (caller handles)."""
    tool = RETRIEVAL_TOOLS_BY_NAME[name]
    args = tool.args_model.model_validate(raw_args or {})
    return tool.fn(args)


__all__ = [
    "RETRIEVAL_TOOLS",
    "RETRIEVAL_TOOLS_BY_NAME",
    "RetrievalTool",
    "run_retrieval_tool",
]
