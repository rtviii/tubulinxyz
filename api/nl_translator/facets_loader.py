"""Adapts the existing Neo4jReader.get_filter_facets() into a compact
FacetContext suitable for injection into an LLM system prompt.
"""
from __future__ import annotations

import os
import time
from typing import List, Optional

from neo4j_tubxz.db_lib_reader import db_reader
from api.nl_translator.interface import FacetContext


# Cap per-field values to keep token usage bounded. Tune as needed.
_MAX_LIGANDS = 20
_MAX_ORGANISMS = 15
_MAX_FAMILIES = 15
_MAX_ISOTYPES = 15
_MAX_EXP_METHODS = 10

# TTL cache for the built FacetContext. Facets change only on ingestion,
# so a coarse TTL is fine. Override with NL_FACET_CACHE_TTL_SECONDS.
_DEFAULT_TTL_SECONDS = 3600
_cache: Optional[FacetContext] = None
_cache_ts: float = 0.0


def _top_n_values(items, n: int) -> List[str]:
    return [i.value for i in items[:n]]


def invalidate_facet_cache() -> None:
    """Drop the cached FacetContext. Call after ingestion runs."""
    global _cache, _cache_ts
    _cache = None
    _cache_ts = 0.0


def _ttl_seconds() -> float:
    raw = os.getenv("NL_FACET_CACHE_TTL_SECONDS")
    if raw is None:
        return float(_DEFAULT_TTL_SECONDS)
    try:
        return float(raw)
    except ValueError:
        return float(_DEFAULT_TTL_SECONDS)


def load_facet_context() -> FacetContext:
    """Query the DB facets endpoint and build a compact context object.

    Cached in-memory with a TTL (env NL_FACET_CACHE_TTL_SECONDS, default 1h).
    Call invalidate_facet_cache() after ingestion to force a rebuild.
    """
    global _cache, _cache_ts
    ttl = _ttl_seconds()
    now = time.monotonic()
    if _cache is not None and (now - _cache_ts) < ttl:
        return _cache

    f = db_reader.get_filter_facets()

    # Source-organism vocabulary: the facets model doesn't carry this today.
    # We pull a small set from the taxonomy reader so Claude can map
    # "human" -> 9606, etc.
    common_organisms = _load_common_organisms()

    ctx = FacetContext(
        exp_methods=_top_n_values(f.exp_methods, _MAX_EXP_METHODS),
        tubulin_families=_top_n_values(f.tubulin_families, _MAX_FAMILIES),
        isotypes=_top_n_values(f.isotypes, _MAX_ISOTYPES),
        top_ligands=[
            {"chemical_id": lf.chemical_id, "chemical_name": lf.chemical_name or ""}
            for lf in f.top_ligands[:_MAX_LIGANDS]
        ],
        common_source_organisms=common_organisms,
        year_range={
            "min": f.year_range.min if f.year_range else None,
            "max": f.year_range.max if f.year_range else None,
        },
        resolution_range={
            "min": f.resolution_range.min if f.resolution_range else None,
            "max": f.resolution_range.max if f.resolution_range else None,
        },
    )
    _cache = ctx
    _cache_ts = now
    return ctx


def _load_common_organisms() -> List[dict]:
    """Return top-N source organisms with tax_id + name + representative PDB IDs.

    Representative PDB IDs (up to 3 per organism) are passed to the LLM as
    grounding so it can name real structures instead of inventing IDs. If
    an organism isn't in this list, the LLM is instructed to fall back to
    an open_catalogue filter on source_organism_ids.
    """
    try:
        with db_reader.adapter.driver.session() as session:
            result = session.run(
                """
                MATCH (e:PolypeptideEntity)
                WHERE e.src_organism_names IS NOT NULL
                  AND e.src_organism_ids IS NOT NULL
                  AND size(e.src_organism_names) > 0
                  AND size(e.src_organism_ids) > 0
                WITH e.src_organism_ids[0] AS tax_id,
                     e.src_organism_names[0] AS name,
                     e.parent_rcsb_id AS pdb
                WITH tax_id, name, collect(DISTINCT pdb) AS pdbs
                RETURN tax_id, name, size(pdbs) AS cnt, pdbs[..3] AS rep_pdbs
                ORDER BY cnt DESC
                LIMIT $limit
                """,
                limit=_MAX_ORGANISMS,
            )
            return [
                {
                    "tax_id": r["tax_id"],
                    "name": r["name"],
                    "structure_count": r["cnt"],
                    "rep_pdbs": list(r["rep_pdbs"]) if r["rep_pdbs"] else [],
                }
                for r in result
            ]
    except Exception:
        # Facet loading must not break the NL route; fall back to empty.
        return []
