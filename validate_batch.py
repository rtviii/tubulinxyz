import json
import os
from pathlib import Path
import random

TUBETL_DATA = Path("/Users/rtviii/dev/TUBETL_DATA")

def validate_structure(rcsb_dir: Path):
    rcsb_id_str = rcsb_dir.name.upper()
    profile_path = rcsb_dir / f"{rcsb_id_str}.json"
    ingestion_path = rcsb_dir / "sequence_ingestion.json"
    
    if not profile_path.exists() or not ingestion_path.exists():
        return None

    with open(profile_path, 'r') as f:
        profile = json.load(f)
    with open(ingestion_path, 'r') as f:
        ingestion = json.load(f)

    # Broader search for ligand JSONs: any JSON in the folder that isn't a known metadata file
    all_jsons = list(rcsb_dir.glob("*.json"))
    ligand_files = [
        f for f in all_jsons 
        if f.name != f"{rcsb_id_str}.json" 
        and f.name != "sequence_ingestion.json" 
        and "classification_report" not in f.name
    ]

    results = []

    for l_file in ligand_files:
        try:
            with open(l_file, 'r') as f:
                molstar_data = json.load(f)
            
            for entry in molstar_data:
                if not entry.get('interactions'): continue
                
                # We need to find a participant that is a polymer
                poly_part = None
                for ix in entry['interactions']:
                    for p in ix['participants']:
                        if p[4] is False: # isLigand == false
                            poly_part = p
                            break
                    if poly_part: break
                
                if not poly_part: continue

                chain_id = poly_part[0]
                auth_seq_id = poly_part[1]

                # Map Chain -> Entity
                entity_id = None
                for eid, ent in profile.get('entities', {}).items():
                    if chain_id in ent.get('pdbx_strand_ids', []):
                        entity_id = eid
                        break
                
                if not entity_id or entity_id not in ingestion:
                    continue

                # Run Round-trip logic
                o2m = ingestion[entity_id]['data']['observed_to_ma_map']
                m2a = ingestion[entity_id]['data']['ma_to_auth_map']
                
                # Check the first mapped residue in the ingestion data for this entity
                # This tells us the starting 'Auth ID' according to the Python ingestor
                first_mapped_auth = next((x for x in m2a if x != -1), None)

                idx = auth_seq_id - 1
                if 0 <= idx < len(o2m):
                    ma_pos = o2m[idx]
                    if ma_pos > 0:
                        reverse_id = m2a[ma_pos - 1]
                        success = (reverse_id == auth_seq_id)
                        results.append({
                            "struct": rcsb_id_str,
                            "chain": chain_id,
                            "auth_id": auth_seq_id,
                            "first_auth": first_mapped_auth,
                            "success": success,
                            "reason": "OK" if success else f"Mismatch: Got {reverse_id}"
                        })
                else:
                    results.append({
                        "struct": rcsb_id_str, 
                        "chain": chain_id, 
                        "auth_id": auth_seq_id,
                        "success": False, 
                        "reason": f"OOB (idx {idx} len {len(o2m)})"
                    })
        except Exception as e:
            continue
    
    return results

def run_batch_validation(sample_size=15):
    all_structs = [d for d in TUBETL_DATA.iterdir() if d.is_dir() and len(d.name) == 4]
    if not all_structs:
        print("No structure directories found.")
        return

    sample = random.sample(all_structs, min(len(all_structs), sample_size))
    
    print(f"{'STRUCT':<8} | {'CHAIN':<6} | {'START':<6} | {'TEST_ID':<8} | {'STATUS':<10} | {'NOTE'}")
    print("-" * 75)
    
    for s_dir in sample:
        v_results = validate_structure(s_dir)
        if not v_results: continue
        
        # Display results for each unique chain found in the ligand files
        seen_chains = set()
        for res in v_results:
            if res['chain'] in seen_chains: continue
            seen_chains.add(res['chain'])
            
            status = "✓" if res['success'] else "✗"
            print(f"{res['struct']:<8} | {res['chain']:<6} | {str(res['first_auth']):<6} | {res['auth_id']:<8} | {status:<10} | {res['reason']}")

if __name__ == "__main__":
    run_batch_validation(100)