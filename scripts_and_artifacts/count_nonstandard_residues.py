#!/usr/bin/env python3
"""
Count non-standard residues across all tubulin chains in the ETL data.

Scans all structure profiles in TUBETL_DATA, identifies tubulin polypeptide
entities (alpha/beta/gamma/delta/epsilon), and catalogs every non-standard
residue comp_id found in those chains via the Molstar raw extraction.

Usage:
    TUBETL_DATA=/path/to/data python scripts_and_artifacts/count_nonstandard_residues.py
"""

import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

TUBETL_DATA = os.environ.get("TUBETL_DATA")
if not TUBETL_DATA:
    print("Error: TUBETL_DATA environment variable not set")
    sys.exit(1)

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

TUBULIN_FAMILIES = {"tubulin_alpha", "tubulin_beta", "tubulin_gamma", "tubulin_delta", "tubulin_epsilon"}


def scan_structure(rcsb_id: str, base_dir: str):
    """Returns list of (comp_id, auth_seq_id, chain, entity_id, family, rcsb_id) for non-standard residues."""
    profile_path = os.path.join(base_dir, f"{rcsb_id}.json")
    molstar_path = os.path.join(base_dir, f"{rcsb_id}_molstar_raw.json")

    if not os.path.exists(profile_path) or not os.path.exists(molstar_path):
        return []

    with open(profile_path) as f:
        profile = json.load(f)

    # Identify tubulin entities and their chain IDs
    tubulin_chains = {}  # auth_asym_id -> (entity_id, family)
    for eid, ent in profile.get("entities", {}).items():
        family = ent.get("family")
        if family and family in TUBULIN_FAMILIES:
            for chain_id in ent.get("pdbx_strand_ids", []):
                tubulin_chains[chain_id] = (eid, family)

    if not tubulin_chains:
        return []

    with open(molstar_path) as f:
        molstar = json.load(f)

    hits = []
    for seq in molstar.get("sequences", []):
        chain_id = seq["auth_asym_id"]
        if chain_id not in tubulin_chains:
            continue

        entity_id, family = tubulin_chains[chain_id]
        for r in seq["residues"]:
            if r["comp_id"] not in STANDARD_AA:
                hits.append((
                    r["comp_id"],
                    r["auth_seq_id"],
                    chain_id,
                    entity_id,
                    family,
                    rcsb_id,
                ))

    return hits


def main():
    data_dir = Path(TUBETL_DATA)
    structure_dirs = sorted([
        d.name for d in data_dir.iterdir()
        if d.is_dir() and (d / f"{d.name}.json").exists()
    ])

    print(f"Scanning {len(structure_dirs)} structures in {TUBETL_DATA}...")

    all_hits = []
    structures_with_nonstandard = set()
    structures_scanned = 0
    tubulin_chains_total = 0

    # Also count total tubulin residues for percentage
    total_tubulin_residues = 0

    for rcsb_id in structure_dirs:
        base = str(data_dir / rcsb_id)
        hits = scan_structure(rcsb_id, base)

        # Count total tubulin residues for this structure
        profile_path = os.path.join(base, f"{rcsb_id}.json")
        molstar_path = os.path.join(base, f"{rcsb_id}_molstar_raw.json")
        if os.path.exists(profile_path) and os.path.exists(molstar_path):
            with open(profile_path) as f:
                profile = json.load(f)
            tub_chain_ids = set()
            for eid, ent in profile.get("entities", {}).items():
                family = ent.get("family")
                if family and family in TUBULIN_FAMILIES:
                    for cid in ent.get("pdbx_strand_ids", []):
                        tub_chain_ids.add(cid)
            if tub_chain_ids:
                tubulin_chains_total += len(tub_chain_ids)
                with open(molstar_path) as f:
                    molstar = json.load(f)
                for seq in molstar.get("sequences", []):
                    if seq["auth_asym_id"] in tub_chain_ids:
                        total_tubulin_residues += len(seq["residues"])

        if hits:
            structures_with_nonstandard.add(rcsb_id)
            all_hits.extend(hits)

        structures_scanned += 1
        if structures_scanned % 100 == 0:
            print(f"  ...{structures_scanned}/{len(structure_dirs)}")

    # Summarize
    print(f"\n{'='*70}")
    print(f"RESULTS")
    print(f"{'='*70}")
    print(f"Structures scanned:              {structures_scanned}")
    print(f"Total tubulin chains:            {tubulin_chains_total}")
    print(f"Total tubulin residues:          {total_tubulin_residues}")
    print(f"Structures with non-standard:    {len(structures_with_nonstandard)}")
    print(f"Total non-standard residues:     {len(all_hits)}")
    if total_tubulin_residues > 0:
        print(f"Non-standard rate:               {100*len(all_hits)/total_tubulin_residues:.3f}%")

    if not all_hits:
        print("\nNo non-standard residues found in tubulin chains.")
        return

    # Count by comp_id
    comp_counts = Counter(h[0] for h in all_hits)
    print(f"\n--- Non-standard residue types (by count) ---")
    for comp_id, count in comp_counts.most_common():
        structs_with = len(set(h[5] for h in all_hits if h[0] == comp_id))
        print(f"  {comp_id:6s}  {count:5d} occurrences  in {structs_with:4d} structures")

    # Count by family
    family_counts = Counter(h[4] for h in all_hits)
    print(f"\n--- By tubulin family ---")
    for fam, count in family_counts.most_common():
        print(f"  {fam:20s}  {count:5d} non-standard residues")

    # Show a few example structures with the most non-standard residues
    per_struct = Counter(h[5] for h in all_hits)
    print(f"\n--- Top 10 structures by non-standard residue count ---")
    for rcsb_id, count in per_struct.most_common(10):
        comps = Counter(h[0] for h in all_hits if h[5] == rcsb_id)
        comp_str = ", ".join(f"{c}x{n}" for c, n in comps.most_common(5))
        print(f"  {rcsb_id}  {count:3d} residues  ({comp_str})")


if __name__ == "__main__":
    main()
