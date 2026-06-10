"""Server-side validation of EntityRefs / ActionCards before returning them
to the frontend.

The LLM may invent PDB ids, chain labels, or ligand chemical ids that aren't
in our database. Rather than ship those to the UI and let the user click into
dead pages, we run one batched Cypher pass per response, dedupe the entities
referenced by all cards, and stamp each card with `{ok, reason}`.

Positive existence results are cached in-memory for 5 min. Most landing-page
queries cluster on the same top ~50 structures and a small set of ligands,
so this cache amortizes well across the request stream.
"""
from __future__ import annotations

import os
import time
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from neo4j_tubxz.db_lib_reader import db_reader

from api.nl_translator.global_actions import (
    ActionCard,
    GlobalNLResponse,
    card_identity_key,
    card_stable_id,
)


# (kind, key) -> (exists, inserted_ts). Negative results not cached; the LLM
# may have a typo we want to catch on the next attempt after data ingest.
_EXIST_CACHE: Dict[Tuple[str, str], Tuple[bool, float]] = {}
_DEFAULT_TTL_SECONDS = 300.0


def _ttl_seconds() -> float:
    raw = os.getenv("NL_HYDRATION_CACHE_TTL_SECONDS")
    if raw is None:
        return _DEFAULT_TTL_SECONDS
    try:
        return float(raw)
    except ValueError:
        return _DEFAULT_TTL_SECONDS


def _cache_get(kind: str, key: str) -> Optional[bool]:
    rec = _EXIST_CACHE.get((kind, key))
    if rec is None:
        return None
    exists, ts = rec
    if (time.monotonic() - ts) > _ttl_seconds():
        _EXIST_CACHE.pop((kind, key), None)
        return None
    return exists if exists else None  # only positives are cached


def _cache_put_positives(kind: str, keys: Iterable[str]) -> None:
    now = time.monotonic()
    for k in keys:
        _EXIST_CACHE[(kind, k)] = (True, now)


def invalidate_hydration_cache() -> None:
    _EXIST_CACHE.clear()


# ---------------------------------------------------------------------------
# Reference extraction
# ---------------------------------------------------------------------------

def _refs_for_card(card: ActionCard) -> Dict[str, Set[str]]:
    """Return the entities this card depends on, grouped by kind.

    Keys per kind:
    - 'structure' -> set of "RCSB"
    - 'chain'     -> set of "RCSB:AUTH"
    - 'ligand'    -> set of "CHEM"
    - 'family'    -> set of "tubulin_xxx"
    """
    out: Dict[str, Set[str]] = {
        "structure": set(),
        "chain": set(),
        "ligand": set(),
        "family": set(),
    }
    a = card.action

    if a == "open_structure":
        if card.rcsb_id:
            out["structure"].add(card.rcsb_id.upper())
        for ch in (card.focus_chains or []):
            if card.rcsb_id and ch:
                out["chain"].add(f"{card.rcsb_id.upper()}:{ch}")
        for lig in (card.focus_ligands or []):
            if lig:
                out["ligand"].add(lig.upper())

    elif a == "open_expert":
        if card.rcsb_id:
            out["structure"].add(card.rcsb_id.upper())
        if card.rcsb_id and card.primary_chain:
            out["chain"].add(f"{card.rcsb_id.upper()}:{card.primary_chain}")
        for ar in (card.aligned or []):
            out["structure"].add(ar.rcsb_id.upper())
            out["chain"].add(f"{ar.rcsb_id.upper()}:{ar.auth_asym_id}")

    elif a == "inspect_ligand":
        if card.chemical_id:
            out["ligand"].add(card.chemical_id.upper())
        if card.rcsb_id:
            out["structure"].add(card.rcsb_id.upper())
        if card.rcsb_id and card.suggested_chain:
            out["chain"].add(f"{card.rcsb_id.upper()}:{card.suggested_chain}")

    elif a == "view_variants":
        if card.family:
            out["family"].add(card.family)

    # open_catalogue and clarify have no entity hydration requirements.
    return out


def _collect_all_refs(cards: List[ActionCard]) -> Dict[str, Set[str]]:
    agg: Dict[str, Set[str]] = {
        "structure": set(), "chain": set(), "ligand": set(), "family": set(),
    }
    for c in cards:
        for kind, keys in _refs_for_card(c).items():
            agg[kind].update(keys)
    return agg


# ---------------------------------------------------------------------------
# Batched Cypher existence checks
# ---------------------------------------------------------------------------

def _check_structures(ids: Set[str]) -> Set[str]:
    """Return the subset of `ids` that exist as :Structure(rcsb_id)."""
    unknown = {i for i in ids if _cache_get("structure", i) is None}
    if not unknown:
        return set(ids)
    cypher = "MATCH (s:Structure) WHERE s.rcsb_id IN $ids RETURN s.rcsb_id AS id"
    try:
        with db_reader.adapter.driver.session() as session:
            found = {r["id"] for r in session.run(cypher, ids=list(unknown))}
    except Exception:
        return set()
    _cache_put_positives("structure", found)
    cached_hits = {i for i in ids if _cache_get("structure", i) is True}
    return cached_hits | found


def _check_chains(pairs: Set[str]) -> Set[str]:
    """`pairs` are 'RCSB:AUTH' strings. Returns the subset that exist as
    :PolypeptideInstance.
    """
    unknown = {p for p in pairs if _cache_get("chain", p) is None}
    if not unknown:
        return set(pairs)
    rows = [{"rcsb": p.split(":", 1)[0], "chain": p.split(":", 1)[1]}
            for p in unknown if ":" in p]
    if not rows:
        return {p for p in pairs if _cache_get("chain", p) is True}
    cypher = """
    UNWIND $rows AS r
    MATCH (p:PolypeptideInstance {parent_rcsb_id: r.rcsb, auth_asym_id: r.chain})
    RETURN DISTINCT r.rcsb + ':' + r.chain AS key
    """
    try:
        with db_reader.adapter.driver.session() as session:
            found = {r["key"] for r in session.run(cypher, rows=rows)}
    except Exception:
        return set()
    _cache_put_positives("chain", found)
    cached_hits = {p for p in pairs if _cache_get("chain", p) is True}
    return cached_hits | found


def _check_ligands(cids: Set[str]) -> Set[str]:
    """Return the subset of chemical ids that exist as :NonpolymerEntity."""
    unknown = {c for c in cids if _cache_get("ligand", c) is None}
    if not unknown:
        return set(cids)
    cypher = (
        "MATCH (n:NonpolymerEntity) WHERE n.chemical_id IN $cids "
        "RETURN DISTINCT n.chemical_id AS id"
    )
    try:
        with db_reader.adapter.driver.session() as session:
            found = {r["id"] for r in session.run(cypher, cids=list(unknown))}
    except Exception:
        return set()
    _cache_put_positives("ligand", found)
    cached_hits = {c for c in cids if _cache_get("ligand", c) is True}
    return cached_hits | found


def _check_families(families: Set[str], known_families: List[str]) -> Set[str]:
    """Family is a controlled vocabulary — check against the facet list, no
    DB call needed.
    """
    known = set(known_families)
    return {f for f in families if f in known}


# ---------------------------------------------------------------------------
# Top-level: hydrate a GlobalNLResponse in place
# ---------------------------------------------------------------------------

def hydrate_response(
    resp: GlobalNLResponse,
    known_families: List[str],
) -> GlobalNLResponse:
    """Validate every card. Returns the same response with `validation` filled
    and any obviously-broken cards either downgraded or dropped.

    Downgrade rules (mirror plan.robustness):
    - open_expert / open_structure with a missing chain: drop the chain (still
      land on the structure page) instead of dropping the whole card.
    - aligned chains that don't exist: drop those entries from `aligned`.
    - Any card whose only point was a missing primary entity (e.g. structure
      doesn't exist at all) is dropped from `cards` entirely, so the UI never
      renders the hallucinated ids. `validation` is re-keyed to the kept
      cards' new contiguous indices.
    """
    if not resp.cards:
        return resp

    refs = _collect_all_refs(resp.cards)
    found_structures = _check_structures(refs["structure"])
    found_chains = _check_chains(refs["chain"])
    found_ligands = _check_ligands(refs["ligand"])
    found_families = _check_families(refs["family"], known_families)

    new_validation: Dict[str, Dict[str, Any]] = {}
    kept_cards: List[ActionCard] = []
    seen_keys: Set[str] = set()

    for card in resp.cards:
        ok = True
        reasons: List[str] = []
        a = card.action

        if a in ("open_structure", "open_expert", "inspect_ligand") and card.rcsb_id:
            if card.rcsb_id.upper() not in found_structures:
                ok = False
                reasons.append(f"structure {card.rcsb_id} not found")

        if ok and a == "open_structure" and card.focus_chains and card.rcsb_id:
            kept = [
                ch for ch in card.focus_chains
                if f"{card.rcsb_id.upper()}:{ch}" in found_chains
            ]
            if kept != list(card.focus_chains):
                reasons.append(f"dropped unknown chains from focus_chains")
            card.focus_chains = kept or None

        if ok and a == "open_structure" and card.focus_ligands:
            kept = [l for l in card.focus_ligands if l.upper() in found_ligands]
            if kept != list(card.focus_ligands):
                reasons.append(f"dropped unknown ligands from focus_ligands")
            card.focus_ligands = kept or None

        if ok and a == "open_expert" and card.rcsb_id and card.primary_chain:
            pair = f"{card.rcsb_id.upper()}:{card.primary_chain}"
            if pair not in found_chains:
                reasons.append(f"primary_chain {card.primary_chain} not in {card.rcsb_id}; dropping chain")
                card.primary_chain = None

        if ok and a == "open_expert" and card.aligned:
            originally = list(card.aligned)
            kept = [
                ar for ar in originally
                if f"{ar.rcsb_id.upper()}:{ar.auth_asym_id}" in found_chains
            ]
            if not kept:
                # All requested aligned chains were unknown. The card's
                # entire premise ("compare X with Y") collapses; surface that
                # as an invalid card instead of silently loading just X.
                missing = ", ".join(f"{a.rcsb_id}:{a.auth_asym_id}" for a in originally)
                ok = False
                reasons.append(f"aligned structure(s) not found: {missing}")
            elif len(kept) != len(originally):
                dropped = [
                    f"{a.rcsb_id}:{a.auth_asym_id}"
                    for a in originally
                    if a not in kept
                ]
                reasons.append(f"dropped unknown aligned chains: {', '.join(dropped)}")
            card.aligned = kept or None

        if ok and a == "inspect_ligand" and card.chemical_id:
            if card.chemical_id.upper() not in found_ligands:
                ok = False
                reasons.append(f"ligand {card.chemical_id} not found")

        if ok and a == "view_variants" and card.family:
            if card.family not in found_families:
                ok = False
                reasons.append(f"family {card.family} not in known set")

        # Drop broken cards entirely. A card that failed existence checks
        # (invented PDB id, missing chain that collapsed the card's premise,
        # unknown ligand/family) must NEVER reach the front page — a dimmed
        # card still leaks the hallucinated ids to the user. Partial
        # downgrades (ok=True with a reason, e.g. some aligned chains dropped)
        # are kept and stay clickable.
        if not ok:
            continue

        # Drop true duplicates and key validation by the card's stable id (not a
        # positional index) so the frontend can't mis-map after dedup/reorder.
        ident = card_identity_key(card)
        if ident in seen_keys:
            continue
        seen_keys.add(ident)
        card.id = card_stable_id(card)
        rec: Dict[str, Any] = {"ok": True}
        if reasons:
            rec["reason"] = "; ".join(reasons)
        new_validation[card.id] = rec
        kept_cards.append(card)

    resp.cards = kept_cards
    resp.validation = new_validation
    return resp


__all__ = ["hydrate_response", "invalidate_hydration_cache"]
