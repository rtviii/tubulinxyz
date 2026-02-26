# verify_ingestion.py
"""
Post-ingestion verification.
Checks that the reingested data is correct and complete.
"""
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

TUBETL = Path(os.environ["TUBETL_DATA"])

# Accumulators
total_structures = 0
structures_missing_profile = 0
structures_missing_binding_sites_file = 0
structures_missing_variants_file = 0

# Binding site checks
bs_total = 0
bs_any_master = 0
bs_all_master = 0
bs_no_master = 0
bs_uses_old_field_name = 0  # "observed_index" instead of "auth_seq_id"
bs_on_tubulin_chain_total = 0
bs_on_tubulin_chain_augmented = 0

# Variant checks
variant_total = 0
variant_subs = 0
variant_ins = 0
variant_del = 0
deletion_heavy_entities = []  # entities where deletions > 50% of canonical length

# Index mapping checks
entities_with_mapping = 0
entities_tubulin_without_mapping = 0
mapping_has_none_master = 0  # canonical positions that map to None (insertions)
mapping_coverage_values = []

# Chain mapping check (can we reconstruct per-chain from entity + molstar?)
chain_mapping_possible = 0
chain_mapping_impossible = 0

struct_dirs = sorted([
    d for d in TUBETL.iterdir()
    if d.is_dir() and len(d.name) == 4
])

print(f"Scanning {len(struct_dirs)} structure directories...\n")

for struct_dir in struct_dirs:
    rcsb_id = struct_dir.name
    profile_path = struct_dir / f"{rcsb_id}.json"
    bs_path = struct_dir / f"{rcsb_id}_ligand_binding_sites.json"
    variants_path = struct_dir / f"{rcsb_id}_variants.json"
    molstar_path = struct_dir / f"{rcsb_id}_molstar_raw.json"

    if not profile_path.exists():
        structures_missing_profile += 1
        continue

    total_structures += 1

    try:
        with open(profile_path) as f:
            profile = json.load(f)
    except Exception as e:
        print(f"  {rcsb_id}: FAILED to load profile: {e}")
        continue

    # ── Identify tubulin entities and their chains ──
    tubulin_entities = set()
    tubulin_chains = set()
    for eid, entity in profile.get("entities", {}).items():
        family = entity.get("family")
        if family and family.startswith("tubulin_"):
            tubulin_entities.add(eid)
            for strand in entity.get("pdbx_strand_ids", []):
                tubulin_chains.add(strand)

    # ── Check index mappings on tubulin entities ──
    for eid in tubulin_entities:
        entity = profile["entities"][eid]
        mapping = entity.get("index_mapping")
        if mapping and mapping.get("observed_to_master"):
            entities_with_mapping += 1
            o2m = mapping["observed_to_master"]
            total_positions = len(o2m)
            mapped_positions = sum(1 for v in o2m.values() if v is not None)
            none_positions = total_positions - mapped_positions
            mapping_has_none_master += none_positions
            if total_positions > 0:
                mapping_coverage_values.append(mapped_positions / total_positions)
        else:
            entities_tubulin_without_mapping += 1
            print(f"  {rcsb_id} entity {eid}: tubulin but NO index mapping")

    # ── Check variants ──
    if variants_path.exists():
        try:
            with open(variants_path) as f:
                var_data = json.load(f)
            for eid, variants in var_data.get("entities", {}).items():
                entity = profile["entities"].get(eid, {})
                can_len = len(entity.get("one_letter_code_can", ""))
                subs = sum(1 for v in variants if v["type"] == "substitution")
                ins = sum(1 for v in variants if v["type"] == "insertion")
                dels = sum(1 for v in variants if v["type"] == "deletion")
                variant_total += len(variants)
                variant_subs += subs
                variant_ins += ins
                variant_del += dels
                if can_len > 0 and dels > can_len * 0.5:
                    deletion_heavy_entities.append({
                        "rcsb_id": rcsb_id,
                        "entity_id": eid,
                        "deletions": dels,
                        "canonical_len": can_len,
                        "pct": round(100 * dels / can_len, 1),
                    })
        except Exception:
            pass
    else:
        if tubulin_entities:
            structures_missing_variants_file += 1

    # ── Check binding sites ──
    if bs_path.exists():
        try:
            with open(bs_path) as f:
                bs_data = json.load(f)

            for site in bs_data.get("binding_sites", []):
                bs_total += 1
                residues = site.get("residues", [])

                # Check field naming
                if residues:
                    keys = set(residues[0].keys())
                    if "observed_index" in keys and "auth_seq_id" not in keys:
                        bs_uses_old_field_name += 1

                has_any = any(r.get("master_index") is not None for r in residues)
                has_all = all(r.get("master_index") is not None for r in residues)
                if has_any:
                    bs_any_master += 1
                if has_all:
                    bs_all_master += 1
                if not has_any:
                    bs_no_master += 1

                # Check tubulin-chain residues specifically
                for r in residues:
                    chain = r.get("auth_asym_id")
                    if chain in tubulin_chains:
                        bs_on_tubulin_chain_total += 1
                        if r.get("master_index") is not None:
                            bs_on_tubulin_chain_augmented += 1

        except Exception as e:
            print(f"  {rcsb_id}: FAILED to load binding sites: {e}")
    else:
        if profile.get("ligand_binding_sites"):
            structures_missing_binding_sites_file += 1

# ═══════════════════════════════════════════════════════
# REPORT
# ═══════════════════════════════════════════════════════
print("=" * 70)
print("INGESTION VERIFICATION REPORT")
print("=" * 70)

print(f"\nStructures scanned:      {total_structures}")
print(f"Missing profiles:        {structures_missing_profile}")

print(f"\n--- INDEX MAPPINGS ---")
print(f"Tubulin entities with mapping:     {entities_with_mapping}")
print(f"Tubulin entities WITHOUT mapping:  {entities_tubulin_without_mapping}")
if mapping_coverage_values:
    avg_cov = sum(mapping_coverage_values) / len(mapping_coverage_values)
    min_cov = min(mapping_coverage_values)
    print(f"Mapping coverage (avg):            {avg_cov:.1%}")
    print(f"Mapping coverage (min):            {min_cov:.1%}")
print(f"Canonical positions -> None (insertions): {mapping_has_none_master}")

print(f"\n--- VARIANTS ---")
print(f"Total variants:       {variant_total}")
print(f"  Substitutions:      {variant_subs}")
print(f"  Insertions:         {variant_ins}")
print(f"  Deletions:          {variant_del}")
if variant_total > 0:
    print(f"  Deletion fraction:  {variant_del / variant_total:.1%}")
if deletion_heavy_entities:
    print(f"\nWARNING: {len(deletion_heavy_entities)} entities with deletions > 50% of canonical length:")
    for ex in deletion_heavy_entities[:10]:
        print(f"  {ex['rcsb_id']} entity {ex['entity_id']}: "
              f"{ex['deletions']} deletions / {ex['canonical_len']} canonical = {ex['pct']}%")

print(f"\n--- BINDING SITES ---")
print(f"Total binding sites:               {bs_total}")
print(f"With any master_index:             {bs_any_master}")
print(f"With ALL master_index:             {bs_all_master}")
print(f"With NO master_index:              {bs_no_master}")
if bs_total > 0:
    print(f"Site-level augmentation rate:       {100 * bs_any_master / bs_total:.1f}%")
print(f"Using old field name (observed_index): {bs_uses_old_field_name}")

print(f"\n--- TUBULIN-CHAIN RESIDUE AUGMENTATION ---")
print(f"Binding site residues on tubulin chains: {bs_on_tubulin_chain_total}")
print(f"Of those with master_index:              {bs_on_tubulin_chain_augmented}")
if bs_on_tubulin_chain_total > 0:
    rate = 100 * bs_on_tubulin_chain_augmented / bs_on_tubulin_chain_total
    print(f"Tubulin residue augmentation rate:       {rate:.1f}%")

# ── Pass/Fail summary ──
print(f"\n{'=' * 70}")
print("PASS/FAIL CRITERIA")
print("=" * 70)

checks = []

checks.append(("No tubulin entities without mapping",
                entities_tubulin_without_mapping == 0,
                f"{entities_tubulin_without_mapping} missing"))

checks.append(("No old field names in binding sites",
                bs_uses_old_field_name == 0,
                f"{bs_uses_old_field_name} using observed_index"))

checks.append(("Tubulin residue augmentation > 95%",
                bs_on_tubulin_chain_total == 0 or
                (bs_on_tubulin_chain_augmented / bs_on_tubulin_chain_total) > 0.95,
                f"{rate:.1f}%" if bs_on_tubulin_chain_total > 0 else "N/A"))

checks.append(("Deletion fraction < 20% of all variants",
                variant_total == 0 or (variant_del / variant_total) < 0.20,
                f"{100 * variant_del / variant_total:.1f}%" if variant_total > 0 else "N/A"))

checks.append(("No deletion-heavy entities (>50% canonical len)",
                len(deletion_heavy_entities) == 0,
                f"{len(deletion_heavy_entities)} found"))

for label, passed, detail in checks:
    status = "PASS" if passed else "FAIL"
    mark = "+" if passed else "!!!"
    print(f"  [{mark}] {status}: {label} ({detail})")