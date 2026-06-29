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

# Cited secondary-structure loop names (M-loop, T7 loop, ...) are layered OVER the
# UniProt regions: precedence 40 beats REGION/MOTIF (30) and BINDING (25), so a
# pocket residue that is also e.g. a GTP-binding residue surfaces its loop name.
_NOGALES_PRECEDENCE = 40
_NOGALES_FILE = "data/genenames/nogales_loops.json"

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


def _reference_sequence(family: TubulinFamily) -> Dict[int, str]:
    """res_no (1-based UniProt) -> amino-acid letter, read from the .afasta
    reference row (ungapped). Used to ASSERT every cited loop anchor really is the
    residue the literature claims, so a published number can never silently point
    at the wrong residue in our coordinate frame."""
    accession = _REFERENCES[family]["accession"]
    path = settings.PROJECT_ROOT / _AFASTA[family]
    aln = AlignIO.read(str(path), "fasta")
    ref_row = next((r for r in aln if accession in r.id or accession in r.description), None)
    if ref_row is None:
        raise ValueError(f"Reference {accession} not found in {path}")
    seq: Dict[int, str] = {}
    res_no = 0
    for ch in str(ref_row.seq):
        if ch in ("-", "."):
            continue
        res_no += 1
        seq[res_no] = ch
    return seq


# ---------------------------------------------------------------------------
# UniProt features (cached)
# ---------------------------------------------------------------------------

def _load_uniprot(accession: str) -> dict:
    """Full UniProt entry (features + audit/version), cached under data/genenames."""
    cache = settings.PROJECT_ROOT / _GENENAMES / f"uniprot_features_{accession}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    url = _UNIPROT_URL.format(acc=accession)
    logger.info(f"Fetching UniProt {accession} from {url}")
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read().decode())
    cache.write_text(json.dumps(data))
    return data


def _feat_span(f: dict) -> Optional[Tuple[int, int]]:
    loc = f.get("location") or {}
    try:
        s = int(loc["start"]["value"])
        e = int(loc["end"]["value"])
    except (KeyError, TypeError, ValueError):
        return None
    return (s, e)


def _evidence(feats: List[dict]) -> List[dict]:
    """De-duplicated UniProt evidence (ECO code + source + id) across features."""
    seen, out = set(), []
    for f in feats:
        for ev in (f.get("evidences") or []):
            key = (ev.get("evidenceCode"), ev.get("source"), ev.get("id"))
            if key in seen:
                continue
            seen.add(key)
            out.append({"eco": ev.get("evidenceCode"), "source": ev.get("source"), "id": ev.get("id")})
    return out


# How regions are sourced + mapped. Embedded verbatim in every table so the JSON
# is self-documenting (the accountability the doc + ENTITY_TOOLING.md point at).
_PROVENANCE = {
    "region_source": "UniProt Knowledgebase (UniProtKB) curated sequence features",
    "fetch": "https://rest.uniprot.org/uniprotkb/{accession}.json (cached: data/genenames/uniprot_features_<acc>.json)",
    "master_mapping": (
        "A UniProt feature residue number is the 1-based position in the human "
        "reference sequence (TUBB/P07437, TUBA1A/Q71U36), which is itself a row in "
        "the family master alignment (.afasta). residue -> master_index is read by "
        "walking the reference's aligned row: master_index = running count of "
        "non-all-gap alignment columns. No MUSCLE re-run, no hand alignment, no "
        "invented positions. Same coordinate space binding sites already use "
        "(lib/etl/augmentation.py)."
    ),
    "eco_glossary": {
        "ECO:0000269": "experimental evidence used in manual assertion",
        "ECO:0000256": "sequence-model (computational) evidence, automatic assertion",
        "ECO:0000250": "sequence-similarity evidence (propagated from a related entry)",
    },
    "citations": [
        "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res.",
        "MREI motif (beta-tubulin autoregulation): PubMed 31727855.",
        "C-terminal disordered tail: MobiDB-lite (ECO:0000256) sequence-based disorder prediction.",
        "GTP / Mg(2+) binding residues: by similarity (ECO:0000250) from related tubulin entries (e.g. Q13509, P68363).",
    ],
}


# ---------------------------------------------------------------------------
# Region selection / labelling
# ---------------------------------------------------------------------------

def _select_regions(
    features: List[dict], res_to_master: Dict[int, int], ref_len: int
) -> List[dict]:
    regions: List[dict] = []

    def add(label: str, kind: str, source: str, span: Tuple[int, int], feats: List[dict]):
        master = _map_span(res_to_master, span[0], span[1])
        if len(master) < 1:
            return
        regions.append({
            "label": label,
            "kind": kind,
            "source": source,
            "ref_span": [span[0], span[1]],
            "ref_residues": [r for r in range(span[0], span[1] + 1) if r in res_to_master],
            "master_indices": master,
            "master_span": [master[0], master[-1]],
            "precedence": _PRECEDENCE[kind],
            "evidence": _evidence(feats),
        })

    # NOTE: v1 deliberately does NOT emit secondary-structure-derived regions
    # (helix/strand/loop). UniProt lists ~24 short helices for tubulin, so a
    # sequential "H12" ordinal would NOT match the canonical Nogales H1-H12 / loop
    # nomenclature a structural biologist expects — labeling them that way would
    # imply a claim we aren't making. The honest contiguous-run fallback groups
    # those residues by master span instead; the deferred Nogales loop table
    # (M-loop, T7, ...) will name them with citations in a follow-on.

    # Pool all binding-site residues into one nucleotide region.
    binding_feats: List[dict] = []
    binding_res: set = set()

    # Named functional features (all correctly named from UniProt curation).
    for f in features:
        t = (f.get("type") or "").lower()
        span = _feat_span(f)
        if not span:
            continue
        desc = (f.get("description") or "")
        if t == "binding site":
            binding_feats.append(f)
            binding_res.update(r for r in range(span[0], span[1] + 1) if r in res_to_master)
        elif t == "region":
            # C-terminal disordered tail.
            if "disordered" in desc.lower() or span[0] > 0.85 * ref_len:
                add("C-terminal tail", "REGION", "uniprot:REGION", span, [f])
        elif t == "motif":
            label = "MREI motif" if "mrei" in desc.lower() else (f"{desc} motif" if desc else "motif")
            add(label, "MOTIF", "uniprot:MOTIF", span, [f])

    if binding_res:
        ref_res = sorted(binding_res)
        master = sorted({res_to_master[r] for r in ref_res})
        ligands = sorted({(f.get("ligand") or {}).get("name") for f in binding_feats
                          if (f.get("ligand") or {}).get("name")})
        regions.append({
            "label": "GTP/nucleotide-binding",
            "kind": "BINDING",
            "source": "uniprot:BINDING",
            "ref_span": None,
            "ref_residues": ref_res,
            "ligands": ligands,
            "master_indices": master,
            "master_span": [master[0], master[-1]],
            "precedence": _PRECEDENCE["BINDING"],
            "evidence": _evidence(binding_feats),
        })

    return regions


# ---------------------------------------------------------------------------
# Nogales secondary-structure loops (cited, layered over the UniProt regions)
# ---------------------------------------------------------------------------

def _load_nogales_loops() -> dict:
    """Cited secondary-structure loop table (hand-curated; the literature analog of
    the cached UniProt snapshot). Returns {} if absent so builds still work without
    it. Each loop carries `anchors` (residue identities) that the selector asserts."""
    path = settings.PROJECT_ROOT / _NOGALES_FILE
    if not path.exists():
        logger.warning(f"no nogales loop table at {path}; skipping loop regions")
        return {}
    return json.loads(path.read_text())


def _select_nogales_regions(
    family: TubulinFamily,
    res_to_master: Dict[int, int],
    ref_seq: Dict[int, str],
) -> List[dict]:
    """Build LOOP regions from the cited Nogales table.

    The honesty gate: every loop's `anchors` (cited residue identities) are checked
    against the on-disk reference sequence. If a cited residue is NOT the residue we
    hold at that position, the build FAILS — the 'verifiable trace into our master
    alignment, or we don't use it' rule, enforced mechanically and re-checked on
    every regenerate so a number can never silently drift out of frame."""
    table = _load_nogales_loops()
    fam_loops = (table.get("families") or {}).get(family.value) or []
    out: List[dict] = []
    for loop in fam_loops:
        label = loop["label"]

        anchors = loop.get("anchors") or []
        if not anchors:
            raise ValueError(f"nogales loop {family.value} {label} has no anchors; "
                             f"refusing an unverifiable loop.")
        for a in anchors:
            pos, aa = a["pos"], a["aa"]
            have = ref_seq.get(pos)
            if have != aa:
                raise ValueError(
                    f"nogales anchor mismatch in {family.value} {label}: cited {aa}{pos} "
                    f"but reference {_REFERENCES[family]['accession']} has {have}{pos}. "
                    f"Refusing to build an untraceable loop."
                )

        # Map ref residues -> master: contiguous span, or an explicit residue set.
        if loop.get("ref_span"):
            s, e = loop["ref_span"][0], loop["ref_span"][1]
            ref_residues = [r for r in range(s, e + 1) if r in res_to_master]
            master = _map_span(res_to_master, s, e)
            ref_span = [s, e]
        else:
            ref_residues = sorted(r for r in (loop.get("ref_residues") or []) if r in res_to_master)
            master = sorted({res_to_master[r] for r in ref_residues})
            ref_span = None
        if not master:
            logger.warning(f"{family.value} {label}: no mapped master indices; skipping")
            continue

        region = {
            "label": label,
            "kind": loop.get("kind", "LOOP"),
            "source": f"nogales:{label}",
            "ref_span": ref_span,
            "ref_residues": ref_residues,
            "master_indices": master,
            "master_span": [master[0], master[-1]],
            "precedence": _NOGALES_PRECEDENCE,
            "evidence": loop.get("evidence") or [],
        }
        # Passthrough context fields if the table carries them (self-documenting JSON).
        for k in ("flanking", "drug", "contacts"):
            if loop.get(k):
                region[k] = loop[k]
        out.append(region)
    logger.info(f"{family.value}: {len(out)} nogales loop regions")
    return out


def build_family(family: TubulinFamily) -> dict:
    ref = _REFERENCES[family]
    res_to_master, master_length = _load_reference_mapping(family)
    ref_len = max(res_to_master) if res_to_master else 0
    doc = _load_uniprot(ref["accession"])
    audit = doc.get("entryAudit", {})
    regions = _select_regions(doc.get("features", []), res_to_master, ref_len)
    regions += _select_nogales_regions(family, res_to_master, _reference_sequence(family))
    table = {
        "family": family.value,
        "reference": {
            "accession": ref["accession"],
            "gene": ref["gene"],
            "uniprot_id": doc.get("uniProtkbId"),
            "entry_version": audit.get("entryVersion"),
            "sequence_version": audit.get("sequenceVersion"),
            "last_annotation_update": audit.get("lastAnnotationUpdateDate"),
        },
        "master_length": master_length,
        "provenance": _PROVENANCE,
        "regions": regions,
    }
    out = settings.PROJECT_ROOT / _GENENAMES / f"regions_{family.value}.json"
    out.write_text(json.dumps(table, indent=2))
    logger.info(f"{family.value}: wrote {len(regions)} regions -> {out}")
    return table


def _write_markdown(tables: List[dict]) -> None:
    """Human-readable accountability doc, generated from the same tables so it can
    never drift from the data. One place to see every region, its UniProt source +
    evidence, the original residue numbers, and the mapped master indices."""
    L: List[str] = []
    L.append("# Structural region tables — provenance & contents")
    L.append("")
    L.append("AUTO-GENERATED by `python -m lib.etl.build_regions --all`. Do not edit by hand.")
    L.append("")
    L.append("## What this is / where it lives")
    L.append("")
    L.append("Named structural regions the landing assistant groups grounded residues into")
    L.append("(hover a region pill -> the cluster lights in 3D). Recorded in:")
    L.append("")
    L.append("- `data/genenames/regions_tubulin_alpha.json`, `..._beta.json` — the tables (read at")
    L.append("  runtime by `api/nl_translator/regions.py::assign_region`).")
    L.append("- `lib/etl/build_regions.py` — the builder (this file's source).")
    L.append("- `data/genenames/uniprot_features_<acc>.json` — the cached raw UniProt snapshot.")
    L.append("")
    prov = tables[0]["provenance"]
    L.append("## How regions are sourced")
    L.append("")
    L.append(f"- **Source:** {prov['region_source']}")
    L.append(f"- **Fetch:** `{prov['fetch']}`")
    L.append(f"- **Master mapping:** {prov['master_mapping']}")
    L.append("")
    L.append("## Citations")
    L.append("")
    for c in prov["citations"]:
        L.append(f"- {c}")
    L.append("")
    L.append("### UniProt evidence codes (ECO)")
    L.append("")
    for code, desc in prov["eco_glossary"].items():
        L.append(f"- `{code}` — {desc}")
    L.append("")
    for t in tables:
        r = t["reference"]
        L.append(f"## {t['family']} — reference {r['gene']} ({r['accession']}, {r['uniprot_id']})")
        L.append("")
        L.append(f"UniProt entry version {r['entry_version']} (seq v{r['sequence_version']}), "
                 f"last annotation update {r['last_annotation_update']}. "
                 f"Master alignment length: {t['master_length']} columns.")
        L.append("")
        L.append("| Region | Source | Evidence | Original residues (ref #) | Master indices |")
        L.append("|--------|--------|----------|---------------------------|----------------|")
        for reg in t["regions"]:
            ev = "; ".join(
                f"{e['eco'].split(':')[-1] if e.get('eco') else '?'}"
                f"{'/' + e['source'] if e.get('source') else ''}"
                f"{':' + str(e['id']) if e.get('id') else ''}"
                for e in (reg.get("evidence") or [])
            ) or "—"
            src = reg["source"]
            if reg.get("ligands"):
                src += " (" + ", ".join(reg["ligands"]) + ")"
            rr = reg.get("ref_residues") or []
            rr_s = (f"{rr[0]}–{rr[-1]} ({len(rr)})" if reg.get("ref_span")
                    else (", ".join(map(str, rr)) if len(rr) <= 16 else f"{len(rr)} residues"))
            mi = reg["master_indices"]
            mi_s = f"{reg['master_span'][0]}–{reg['master_span'][1]} ({len(mi)})" if len(mi) > 4 else ", ".join(map(str, mi))
            L.append(f"| **{reg['label']}** | {src} | {ev} | {rr_s} | {mi_s} |")
        L.append("")
    out = settings.PROJECT_ROOT / _GENENAMES / "REGIONS.md"
    out.write_text("\n".join(L))
    logger.info(f"wrote provenance doc -> {out}")


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

    tables = [build_family(fam) for fam in fams]
    # The provenance doc covers all built families in one place; only (re)write it
    # on a full build so it never describes a partial set.
    if args.all:
        _write_markdown(tables)


if __name__ == "__main__":
    main()
