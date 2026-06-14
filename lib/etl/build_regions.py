# lib/etl/build_regions.py
"""
Build per-family structural region tables in master-alignment coordinates.

Regions are sourced HONESTLY from UniProt curated features for each family's
human reference isotype (the same reference sequence that already lives in the
family master alignment .afasta), then mapped into our master_index space by
reading the reference's aligned row directly from the .afasta. No region
boundary is authored by hand; every range traces to a UniProt feature.

Why read the .afasta row directly instead of re-aligning with MUSCLE: the
master alignment IS the .afasta, and the human reference (TUBB / TUBA1A) is one
of its rows. master_index is the 1-based count of "original" alignment columns
(columns where any sequence has a residue) -- exactly how SequenceAligner builds
master indices for binding sites. So walking the reference row gives the same
coordinates with no MUSCLE dependency and no re-alignment drift.

v1 region kinds (Nogales named loops M-loop/T7 deferred to a cited follow-on;
each region carries a `source` so a later override can rename a generic loop):
  - C-terminal tail        (UniProt REGION 'Disordered' at the C-terminus)
  - GTP/nucleotide-binding (UniProt BINDING sites, pooled)
  - MREI motif             (UniProt MOTIF)
  - loop between X and Y   (gaps between consecutive secondary-structure elements)
  - helix H{n} / strand S{n} (UniProt secondary structure; low precedence)

Usage:
  python -m lib.etl.build_regions --family beta
  python -m lib.etl.build_regions --family alpha
  python -m lib.etl.build_regions --all
"""

import argparse
import json
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import AlignIO
from loguru import logger

from api.config import settings
from lib.types import TubulinFamily


# Human reference isotype per family (accession is the key into UniProt + the
# afasta header). These rows live verbatim in the family .afasta.
_REFERENCES: Dict[TubulinFamily, Dict[str, str]] = {
    TubulinFamily.ALPHA: {"accession": "Q71U36", "gene": "TUBA1A"},
    TubulinFamily.BETA: {"accession": "P07437", "gene": "TUBB"},
}

_AFASTA: Dict[TubulinFamily, str] = {
    TubulinFamily.ALPHA: "data/alpha_tubulin/alpha_tubulin.afasta",
    TubulinFamily.BETA: "data/beta_tubulin/beta_tubulin.afasta",
}

# Higher precedence wins when a master position belongs to several regions.
_PRECEDENCE = {"REGION": 30, "MOTIF": 30, "BINDING": 25, "LOOP": 15, "HELIX": 5, "STRAND": 5}

_GENENAMES = "data/genenames"
_UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"


# ---------------------------------------------------------------------------
# Master-alignment row -> master_index mapping
# ---------------------------------------------------------------------------

def _load_reference_mapping(
    family: TubulinFamily,
) -> Tuple[Dict[int, int], int]:
    """Read the family .afasta, find the human reference row, and return
    (uniprot_residue_number -> master_index, master_length).

    master_index = 1-based count of original columns (any-row-non-gap) up to and
    including the column, matching SequenceAligner. The reference row is non-gap
    only at original columns, so every residue maps to a real master_index.
    """
    accession = _REFERENCES[family]["accession"]
    path = settings.PROJECT_ROOT / _AFASTA[family]
    aln = AlignIO.read(str(path), "fasta")
    ncols = aln.get_alignment_length()

    # master_index per column (None for all-gap columns).
    master_of_col: List[Optional[int]] = []
    running = 0
    for c in range(ncols):
        column = aln[:, c]
        is_original = any(ch not in ("-", ".") for ch in column)
        if is_original:
            running += 1
            master_of_col.append(running)
        else:
            master_of_col.append(None)
    master_length = running

    ref_row = next((r for r in aln if accession in r.id or accession in r.description), None)
    if ref_row is None:
        raise ValueError(f"Reference {accession} not found in {path}")

    res_to_master: Dict[int, int] = {}
    res_no = 0
    for c, ch in enumerate(str(ref_row.seq)):
        if ch in ("-", "."):
            continue
        res_no += 1
        m = master_of_col[c]
        if m is not None:
            res_to_master[res_no] = m

    logger.info(
        f"{family.value}: reference {accession} has {res_no} residues; "
        f"master_length={master_length}"
    )
    return res_to_master, master_length


def _map_span(res_to_master: Dict[int, int], start: int, end: int) -> List[int]:
    return sorted({res_to_master[r] for r in range(start, end + 1) if r in res_to_master})


# ---------------------------------------------------------------------------
# UniProt features (cached)
# ---------------------------------------------------------------------------

def _fetch_features(accession: str) -> List[dict]:
    """UniProt curated features for an accession, cached under data/genenames."""
    cache = settings.PROJECT_ROOT / _GENENAMES / f"uniprot_features_{accession}.json"
    if cache.exists():
        data = json.loads(cache.read_text())
    else:
        url = _UNIPROT_URL.format(acc=accession)
        logger.info(f"Fetching UniProt {accession} from {url}")
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        cache.write_text(json.dumps(data))
    return data.get("features", [])


def _feat_span(f: dict) -> Optional[Tuple[int, int]]:
    loc = f.get("location") or {}
    try:
        s = int(loc["start"]["value"])
        e = int(loc["end"]["value"])
    except (KeyError, TypeError, ValueError):
        return None
    return (s, e)


# ---------------------------------------------------------------------------
# Region selection / labelling
# ---------------------------------------------------------------------------

def _select_regions(
    features: List[dict], res_to_master: Dict[int, int], ref_len: int
) -> List[dict]:
    regions: List[dict] = []

    def add(label: str, kind: str, source: str, span: Tuple[int, int]):
        master = _map_span(res_to_master, span[0], span[1])
        if len(master) < 1:
            return
        regions.append({
            "label": label,
            "kind": kind,
            "source": source,
            "ref_span": [span[0], span[1]],
            "master_indices": master,
            "master_span": [master[0], master[-1]],
            "precedence": _PRECEDENCE[kind],
        })

    # Pool all binding-site residues into one nucleotide region.
    binding_master: set = set()

    # NOTE: v1 deliberately does NOT emit secondary-structure-derived regions
    # (helix/strand/loop). UniProt lists ~24 short helices for tubulin, so a
    # sequential "H12" ordinal would NOT match the canonical Nogales H1-H12 / loop
    # nomenclature a structural biologist expects — labeling them that way would
    # imply a claim we aren't making. The honest contiguous-run fallback groups
    # those residues by master span instead; the deferred Nogales loop table
    # (M-loop, T7, ...) will name them with citations in a follow-on.

    # Named functional features (all correctly named from UniProt curation).
    for f in features:
        t = (f.get("type") or "").lower()
        span = _feat_span(f)
        if not span:
            continue
        desc = (f.get("description") or "")
        if t == "binding site":
            binding_master.update(_map_span(res_to_master, span[0], span[1]))
        elif t == "region":
            # C-terminal disordered tail.
            if "disordered" in desc.lower() or span[0] > 0.85 * ref_len:
                add("C-terminal tail", "REGION", "uniprot:REGION", span)
        elif t == "motif":
            label = "MREI motif" if "mrei" in desc.lower() else (f"{desc} motif" if desc else "motif")
            add(label, "MOTIF", "uniprot:MOTIF", span)

    if binding_master:
        bm = sorted(binding_master)
        regions.append({
            "label": "GTP/nucleotide-binding",
            "kind": "BINDING",
            "source": "uniprot:BINDING",
            "ref_span": None,
            "master_indices": bm,
            "master_span": [bm[0], bm[-1]],
            "precedence": _PRECEDENCE["BINDING"],
        })

    return regions


def build_family(family: TubulinFamily) -> dict:
    ref = _REFERENCES[family]
    res_to_master, master_length = _load_reference_mapping(family)
    ref_len = max(res_to_master) if res_to_master else 0
    features = _fetch_features(ref["accession"])
    regions = _select_regions(features, res_to_master, ref_len)
    table = {
        "family": family.value,
        "reference": {"accession": ref["accession"], "gene": ref["gene"]},
        "master_length": master_length,
        "regions": regions,
    }
    out = settings.PROJECT_ROOT / _GENENAMES / f"regions_{family.value}.json"
    out.write_text(json.dumps(table, indent=2))
    logger.info(f"{family.value}: wrote {len(regions)} regions -> {out}")
    return table


def main():
    p = argparse.ArgumentParser(description="Build per-family region tables in master coords.")
    p.add_argument("--family", choices=["alpha", "beta"])
    p.add_argument("--all", action="store_true")
    args = p.parse_args()

    fams = []
    if args.all:
        fams = [TubulinFamily.ALPHA, TubulinFamily.BETA]
    elif args.family == "alpha":
        fams = [TubulinFamily.ALPHA]
    elif args.family == "beta":
        fams = [TubulinFamily.BETA]
    else:
        p.error("pass --family alpha|beta or --all")

    for fam in fams:
        build_family(fam)


if __name__ == "__main__":
    main()
