import json
import os
from pathlib import Path

TUBETL_DATA = Path("/Users/rtviii/dev/TUBETL_DATA")

def verify_structure_mapping(rcsb_id: str):
    rcsb_id = rcsb_id.upper()
    struct_dir = TUBETL_DATA / rcsb_id
    
    # Files
    profile_path = struct_dir / f"{rcsb_id}.json"
    ingestion_path = struct_dir / "sequence_ingestion.json"
    all_jsons = list(struct_dir.glob(f"{rcsb_id}_*.json"))
    ligand_files = [f for f in all_jsons if not any(x in f.name for x in ["report", "ingestion"])]

    if not all([profile_path.exists(), ingestion_path.exists(), ligand_files]):
        print(f"Error: Missing required files in {struct_dir}")
        return

    # 1. Load Data
    with open(profile_path, 'r') as f:
        profile = json.load(f)
    with open(ingestion_path, 'r') as f:
        ingestion = json.load(f)
    with open(ligand_files[0], 'r') as f:
        molstar_data = json.load(f)

    # 2. Get a sample interaction from Chain D
    first_ligand_entry = molstar_data[0]
    sample_ix = first_ligand_entry['interactions'][0]
    poly_part = next(p for p in sample_ix['participants'] if not p[4])
    
    chain_id = poly_part[0]    # e.g., "D"
    auth_seq_id = poly_part[1] # e.g., 142
    res_name = poly_part[2]

    print(f"=== Verification Report for {rcsb_id} ===")
    print(f"Targeting Interaction: {res_name} {auth_seq_id} on Chain {chain_id}")

    # 3. Resolve Chain ID -> Entity ID
    entity_id = None
    for eid, ent_data in profile['entities'].items():
        # Check if our chain_id is in the pdbx_strand_ids list for this entity
        if chain_id in ent_data.get('pdbx_strand_ids', []):
            entity_id = eid
            break
    
    if not entity_id:
        print(f"Failure: Could not map Chain {chain_id} to any Entity in {rcsb_id}.json")
        return

    print(f"Resolved: Chain {chain_id} is an instance of Entity {entity_id}")

    # 4. Get Ingestion Data for that Entity
    chain_data = ingestion.get(entity_id)
    if not chain_data:
        print(f"Failure: Entity {entity_id} not found in sequence_ingestion.json")
        return

    o2m = chain_data['data']['observed_to_ma_map']
    m2a = chain_data['data']['ma_to_auth_map']

    # 5. VERACITY CHECK
    # Check first 5 residues for N-terminal alignment
    print(f"\n--- N-Terminal Diagnostic (Entity {entity_id}) ---")
    print(f"{'Index':<8} | {'Auth ID (Assumed)':<18} | {'MA Pos':<8}")
    for i in range(min(5, len(o2m))):
        print(f"{i:<8} | {i+1:<18} | {o2m[i]:<8}")

    idx = auth_seq_id - 1
    if 0 <= idx < len(o2m):
        ma_pos = o2m[idx]
        print(f"\nValidation:")
        print(f"1. Physical ID {auth_seq_id} maps to o2m index {idx}")
        print(f"2. MA Position: {ma_pos}")
        
        if ma_pos > 0:
            reverse_physical = m2a[ma_pos - 1]
            print(f"3. Reverse Check: MA {ma_pos} -> Physical ID {reverse_physical}")
            
            if reverse_physical == auth_seq_id:
                print("\nRESULT: [ MATCH SUCCESS ]")
            else:
                print(f"\nRESULT: [ OFFSET ERROR ] - Logic suggests ID {reverse_physical}")
        else:
            print("\nRESULT: [ GAP ]")
    else:
        print(f"\nRESULT: [ OUT OF BOUNDS ]")

if __name__ == "__main__":
    verify_structure_mapping("1FFX")