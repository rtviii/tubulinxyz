# check_deletion_heavy.py
import json, os
from pathlib import Path

TUBETL = Path(os.environ["TUBETL_DATA"])

for rcsb_id, eid in [("3RB8", "1"), ("4XCQ", "1")]:
    profile_path = TUBETL / rcsb_id / f"{rcsb_id}.json"
    with open(profile_path) as f:
        profile = json.load(f)

    entity = profile["entities"][eid]
    print(f"\n{rcsb_id} entity {eid}:")
    print(f"  family:          {entity.get('family')}")
    print(f"  description:     {entity.get('pdbx_description')}")
    print(f"  canonical len:   {len(entity.get('one_letter_code_can', ''))}")
    print(f"  organism:        {entity.get('src_organism_names')}")
    stats = entity.get("alignment_stats", {})
    print(f"  MA coverage:     {stats.get('ma_coverage')} / {stats.get('master_alignment_length')}")
    print(f"  subs/ins/del:    {stats.get('substitutions')}/{stats.get('insertions')}/{stats.get('deletions')}")
