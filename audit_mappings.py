import json
import os
from pathlib import Path

TUBETL_DATA = Path(os.environ.get("TUBETL_DATA", "/Users/rtviii/dev/TUBETL_DATA"))

def audit_all_mappings():
    print(f"Auditing all sequence mappings in {TUBETL_DATA}...\n")
    print(f"{'PDB':<6} | {'Entity':<6} | {'Family':<10} | {'Shift Detail':<25} | {'Augmentation Status'}")
    print("-" * 90)

    for struct_dir in TUBETL_DATA.iterdir():
        if not struct_dir.is_dir() or len(struct_dir.name) != 4: continue
        ingestion_path = struct_dir / "sequence_ingestion.json"
        if not ingestion_path.exists(): continue

        with open(ingestion_path, 'r') as f:
            ingestion_data = json.load(f)

        for entity_id, entry in ingestion_data.items():
            mapping = entry.get("data", {}).get("observed_to_ma_map", [])
            if not mapping: continue

            # Find the first real shift
            shifts = [i for i, ma in enumerate(mapping) if ma != (i + 1) and ma > 0]
            if shifts:
                idx = shifts[0]
                pdb_val, ma_val = idx + 1, mapping[idx]
                
                res = check_ligands(struct_dir, pdb_val, ma_val)
                print(f"{struct_dir.name:<6} | {entity_id:<6} | {entry.get('family', '??'):<10} | PDB {pdb_id:<3} -> MA {ma_val:<3} | {res}")

def check_ligands(struct_dir, pdb_val, ma_val):
    l_files = [f for f in struct_dir.glob(f"{struct_dir.name.upper()}_*_*.json") if f.name != f"{struct_dir.name.upper()}.json"]
    if not l_files: return "No Ligands"

    for lf in l_files:
        with open(lf, 'r') as f:
            raw = json.load(f)
            if not raw: continue # Handle empty list []
            data = raw[0] if isinstance(raw, list) else raw
        
        for res in data.get("neighborhood", []):
            if res[1] == pdb_val:
                if len(res) < 4: return "NOT AUGMENTED"
                return "SUCCESS" if res[3] == ma_val else f"BUG: Mirrored PDB {res[3]}"
    
    return "Not in site"

if __name__ == "__main__":
    audit_all_mappings()