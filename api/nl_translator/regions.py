"""Runtime lookup over the per-family structural region tables.

The tables (data/genenames/regions_{family}.json) are built offline by
lib.etl.build_regions from UniProt curated features mapped into our master
alignment. Here we only READ them during harvest:

  - assign_region(family, master_index): which named region (if any) owns this
    master position, resolving overlap by precedence.

Honesty note: this never invents positions. A region's residue membership comes
straight from the UniProt-derived table; the assistant only ever surfaces the
subset of those positions it actually looked up and grounded.
"""

import json
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

from api.config import settings


def _table_path(family: str) -> str:
    return str(settings.PROJECT_ROOT / "data" / "genenames" / f"regions_{family}.json")


@lru_cache(maxsize=8)
def _load_table(family: str) -> Optional[dict]:
    try:
        with open(_table_path(family)) as f:
            return json.load(f)
    except FileNotFoundError:
        return None


@lru_cache(maxsize=8)
def _inverted_index(family: str) -> Dict[int, str]:
    """master_index -> winning region label. Highest precedence wins; ties go to
    the smaller (more specific) region."""
    table = _load_table(family)
    if not table:
        return {}
    # (master_index) -> (precedence, -size, label) for argmax comparison.
    best: Dict[int, Tuple[int, int, str]] = {}
    for r in table.get("regions", []):
        label = r["label"]
        prec = r.get("precedence", 0)
        mis = r.get("master_indices") or []
        key = (prec, -len(mis), label)
        for mi in mis:
            cur = best.get(mi)
            if cur is None or key > cur:
                best[mi] = key
    return {mi: v[2] for mi, v in best.items()}


def assign_region(family: Optional[str], master_index: int) -> Optional[str]:
    if not family:
        return None
    return _inverted_index(family).get(master_index)


def region_master_indices(family: str, label: str) -> List[int]:
    table = _load_table(family)
    if not table:
        return []
    for r in table.get("regions", []):
        if r["label"] == label:
            return list(r.get("master_indices") or [])
    return []
