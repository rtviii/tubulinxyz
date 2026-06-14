"""Ground canonical (master) MSA positions to a real structure/chain's auth_seq_id.

Why this exists
---------------
Every position-returning read tool (get_binding_site, get_binding_contacts,
count_modifications, count_variants) reports MASTER positions — columns in the
family MSA consensus, NOT a structure's author residue numbering. Those numbers
cannot be handed to the 3D viewer directly: master_index 224 is not auth_seq_id
224. To highlight a residue on a *specific* loaded structure (e.g. the landing
demo 9MLF chain B = beta-tubulin) we must translate master -> auth_seq_id for
THAT chain.

The translation already exists, precomputed, in the structure profile:
`entities[*].chain_index_mappings[auth_asym_id].master_to_auth_seq_id`
(built during ETL, see lib/etl/collector.py). So grounding is a dict lookup, not
a live alignment.

Honesty contract: a master position with no residue in this chain (gap /
unmodeled density) maps to None and is dropped by callers — we never invent an
auth_seq_id. If the profile or the chain mapping is missing, we ground nothing.

Mirrors the cache conventions of resolve.py (5-min positive cache, keyed by the
chain identity).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from lib.etl.assets import TubulinStructureAssets


# (rcsb_id, auth_asym_id, family) -> (master_index -> auth_seq_id | None, ts).
# We cache the FULL chain map (cheap, a few hundred entries) so repeated
# groundings for the same chain within a turn are free.
_GROUNDING_CACHE: Dict[Tuple[str, str, Optional[str]], Tuple[Dict[int, Optional[int]], float]] = {}
_TTL_SECONDS = 300.0


def invalidate_grounding_cache() -> None:
    _GROUNDING_CACHE.clear()


def _load_chain_master_map(
    rcsb_id: str, auth_asym_id: str, family: Optional[str]
) -> Dict[int, Optional[int]]:
    """Return the chain's full master_index -> auth_seq_id map (ints), or {} if
    the profile / entity / chain mapping is absent. JSON object keys arrive as
    strings; we normalize both keys and values to ints (values may be None)."""
    key = (rcsb_id.upper(), auth_asym_id, family or None)

    rec = _GROUNDING_CACHE.get(key)
    if rec is not None:
        chain_map, ts = rec
        if (time.monotonic() - ts) <= _TTL_SECONDS:
            return chain_map
        _GROUNDING_CACHE.pop(key, None)

    try:
        assets = TubulinStructureAssets(rcsb_id.upper())
        profile_path = Path(assets.paths.profile)
        if not profile_path.exists():
            return {}
        with open(profile_path, "r") as f:
            profile = json.load(f)
    except Exception:
        return {}

    entities = profile.get("entities") or {}
    raw_map: Optional[Dict] = None
    for ent in entities.values():
        if not isinstance(ent, dict):
            continue
        if ent.get("type") != "polymer":
            continue
        mappings = ent.get("chain_index_mappings") or {}
        if auth_asym_id not in mappings:
            continue
        # Require the family to match when one was given — never ground onto a
        # chain of a different family (that would mislabel the residue).
        if family and ent.get("family") not in (None, family):
            continue
        raw_map = mappings[auth_asym_id].get("master_to_auth_seq_id")
        break

    if not isinstance(raw_map, dict):
        return {}

    chain_map: Dict[int, Optional[int]] = {}
    for mk, mv in raw_map.items():
        try:
            mi = int(mk)
        except (TypeError, ValueError):
            continue
        chain_map[mi] = int(mv) if isinstance(mv, (int, float)) else None

    _GROUNDING_CACHE[key] = (chain_map, time.monotonic())
    return chain_map


def ground_master_positions(
    rcsb_id: str,
    auth_asym_id: str,
    family: Optional[str],
    master_positions: List[int],
) -> Dict[int, Optional[int]]:
    """Map each master position to this chain's auth_seq_id.

    Returns {master_index: auth_seq_id | None}. None means the position has no
    residue in this chain (gap/unmodeled) OR the chain/profile couldn't be
    loaded — callers drop None entries rather than invent coordinates.
    """
    chain_map = _load_chain_master_map(rcsb_id, auth_asym_id, family)
    return {mp: chain_map.get(mp) for mp in master_positions}


def ground_auth_seq_ids(
    rcsb_id: str,
    auth_asym_id: str,
    family: Optional[str],
    master_positions: List[int],
) -> List[int]:
    """Convenience: just the resolved auth_seq_ids, ungroundable positions
    dropped, de-duplicated, ascending. Use when you only need the residues that
    actually exist on the chain (the common case for highlighting)."""
    grounded = ground_master_positions(rcsb_id, auth_asym_id, family, master_positions)
    return sorted({v for v in grounded.values() if v is not None})


__all__ = [
    "ground_master_positions",
    "ground_auth_seq_ids",
    "invalidate_grounding_cache",
]
