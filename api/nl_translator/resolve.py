"""Resolve LLM-expressed entity *intent* to REAL structure ids from the DB.

Why this exists
---------------
The front-page LLM was authoring PDB ids from its training memory. Two failure
modes followed:
  1. Invented ids that don't exist (e.g. "8FEH") — caught by hydration.
  2. Real ids cited with the WRONG organism/family (e.g. picking a pig
     structure and calling it human). Existence checks can't catch this.

The fix: for entity-bearing cards (open_structure / open_expert /
inspect_ligand) the model expresses *intent* — organism tax id, tubulin
family, optional ligand — and we resolve it HERE to a concrete
(rcsb_id, auth_asym_id) by querying Neo4j. The ids become real AND
semantically correct by construction.

Granularity note (important): organism is resolved at the ENTITY level
(`PolypeptideEntity.src_organism_ids`), NOT the parent structure's organism
union. A microtubule structure can mix organisms (pig tubulin + human MAP), so
a structure-level filter would mislabel. We filter the entity's own organism.

Representative selection: best (lowest) resolution wins. Positive results are
cached in-memory for 5 min — landing-page queries cluster on a few
(organism, family) pairs.
"""
from __future__ import annotations

import time
from typing import Dict, List, Optional, Tuple

from neo4j_tubxz.db_lib_reader import db_reader

from api.nl_translator.global_actions import ActionCard, AlignedRef, GlobalNLResponse


# (tax_id, family, ligand) -> (rcsb_id, auth_asym_id) | None. Negatives are
# NOT cached (data may land on the next ingest).
_REP_CACHE: Dict[Tuple[Optional[int], Optional[str], Optional[str]], Tuple[Tuple[str, str], float]] = {}
_TTL_SECONDS = 300.0


# Entity-level resolution. Ranked by resolution ascending (best first). The
# ligand clause mirrors PolypeptideEntityQueryBuilder's NEAR_POLYMER traversal
# so "structure with ligand X near this family" stays instance-accurate.
_RESOLVE_CYPHER = """
MATCH (s:Structure)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)
WHERE ($family IS NULL OR e.family = $family)
  AND ($tax IS NULL OR $tax IN e.src_organism_ids)
  AND size(coalesce(e.pdbx_strand_ids, [])) > 0
  AND ($lig IS NULL OR EXISTS {
        MATCH (pi:PolypeptideInstance)-[:INSTANCE_OF]->(e)
        MATCH (ni:NonpolymerInstance)-[:NEAR_POLYMER]->(pi)
        MATCH (ni)-[:INSTANCE_OF]->(:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
        WHERE c.chemical_id = $lig
      })
RETURN e.parent_rcsb_id AS rcsb, e.pdbx_strand_ids AS chains
ORDER BY coalesce(s.resolution, 9999.0) ASC
LIMIT 1
"""


def _cache_get(key: Tuple[Optional[int], Optional[str], Optional[str]]) -> Optional[Tuple[str, str]]:
    rec = _REP_CACHE.get(key)
    if rec is None:
        return None
    value, ts = rec
    if (time.monotonic() - ts) > _TTL_SECONDS:
        _REP_CACHE.pop(key, None)
        return None
    return value


def invalidate_resolve_cache() -> None:
    _REP_CACHE.clear()


def resolve_representative(
    *,
    organism_id: Optional[int] = None,
    family: Optional[str] = None,
    ligand: Optional[str] = None,
) -> Optional[Tuple[str, str]]:
    """Return (rcsb_id, auth_asym_id) for the best real polypeptide entity
    matching the intent, or None if nothing matches.

    At least one of organism_id / family / ligand should be set; calling with
    all-None returns None (we never pick a global "representative").
    """
    if organism_id is None and not family and not ligand:
        return None

    fam = family or None
    lig = ligand.upper() if ligand else None
    key = (organism_id, fam, lig)

    cached = _cache_get(key)
    if cached is not None:
        return cached

    try:
        with db_reader.adapter.driver.session() as session:
            rec = session.run(
                _RESOLVE_CYPHER, tax=organism_id, family=fam, lig=lig
            ).single()
    except Exception:
        return None

    if rec is None:
        return None
    rcsb = rec["rcsb"]
    chains = rec["chains"] or []
    if not rcsb or not chains:
        return None

    result = (str(rcsb).upper(), str(chains[0]))
    _REP_CACHE[key] = (result, time.monotonic())
    return result


# ---------------------------------------------------------------------------
# Card-level resolution
# ---------------------------------------------------------------------------

def _resolve_primary(card: ActionCard) -> bool:
    """Fill the card's primary structure+chain from intent. Returns False if
    resolution was required but failed (caller drops the card).

    Rule: if the model set `primary_organism_id`, resolve from intent and
    OVERWRITE any guessed rcsb_id (the selector is more trustworthy than a
    guessed id). If no organism selector and an rcsb_id is already present,
    treat it as a user-named structure and keep it.
    """
    a = card.action
    org = card.primary_organism_id
    fam = card.family
    lig = card.chemical_id if a == "inspect_ligand" else None

    needs_resolution = org is not None or (not card.rcsb_id and (bool(fam) or bool(lig)))

    if not needs_resolution:
        # User-named structure (rcsb_id present, no organism selector), or a
        # bare inspect_ligand that routes to the ligand catalogue.
        if card.rcsb_id:
            return True
        return a == "inspect_ligand" and bool(card.chemical_id)

    rep = resolve_representative(organism_id=org, family=fam, ligand=lig)
    if rep is None:
        if a == "inspect_ligand" and card.chemical_id:
            # No bound structure found, but the ligand is real enough to route
            # to the catalogue filtered by chemical id. Drop the structure
            # context so the dispatcher falls back to the ligand catalogue.
            card.rcsb_id = None
            card.suggested_chain = None
            return True
        return False

    rcsb, chain = rep
    card.rcsb_id = rcsb
    if a == "open_expert":
        card.primary_chain = chain
    elif a == "open_structure":
        card.focus_chains = [chain]
    elif a == "inspect_ligand":
        card.suggested_chain = chain
    return True


def _resolve_aligned(card: ActionCard) -> None:
    """Replace `aligned` with DB-resolved refs derived from
    `aligned_organism_ids` (each resolved against the card's family). If the
    model didn't supply organism selectors, leave `aligned` untouched (it may
    be None or already concrete from a user-named id)."""
    org_ids = card.aligned_organism_ids
    if not org_ids:
        return
    fam = card.family
    refs: List[AlignedRef] = []
    for org in org_ids:
        rep = resolve_representative(organism_id=org, family=fam)
        if rep:
            refs.append(AlignedRef(rcsb_id=rep[0], auth_asym_id=rep[1]))
    card.aligned = refs or None


def resolve_card(card: ActionCard) -> bool:
    """Resolve one entity card's ids from intent, in place. Returns False if
    the card should be dropped (its primary entity could not be resolved to a
    real structure). Non-entity cards (open_catalogue / view_variants /
    clarify) always return True untouched.
    """
    a = card.action
    if a in ("open_structure", "open_expert", "inspect_ligand"):
        if not _resolve_primary(card):
            return False
    if a == "open_expert":
        _resolve_aligned(card)
    return True


def resolve_response(resp: GlobalNLResponse) -> GlobalNLResponse:
    """Resolve all entity-bearing cards in place. Cards whose primary entity
    cannot be resolved are dropped — better to show fewer real cards than a
    card pointing at a guessed/wrong structure. Runs BEFORE hydration; the
    existence check in hydration.py is then a cheap second line of defense.
    """
    if not resp.cards:
        return resp
    resp.cards = [c for c in resp.cards if resolve_card(c)]
    return resp


__all__ = [
    "resolve_representative",
    "resolve_card",
    "resolve_response",
    "invalidate_resolve_cache",
]
