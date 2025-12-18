import json
import os
from pathlib import Path

TUBETL_DATA = Path(os.environ.get("TUBETL_DATA", "/Users/rtviii/dev/TUBETL_DATA"))

def scan_augmented_ligands():
    print(f"Searching for non-congruent ligand mappings in {TUBETL_DATA}...\n")
    print(f"{'PDB':<6} | {'Ligand':<10} | {'Residue':<8} | {'PDB ID':<6} -> {'MA Index':<8}")
    print("-" * 60)

    found_count = 0

    # Iterate through structure directories
    for struct_dir in TUBETL_DATA.iterdir():
        if not struct_dir.is_dir() or len(struct_dir.name) != 4:
            continue
        
        # Look for augmented ligand JSONs (e.g., 6WVR_GDP_B.json)
        # We exclude the main 6WVR.json profile
        ligand_files = [
            f for f in struct_dir.glob(f"{struct_dir.name}_*_*.json") 
            if f.name != f"{struct_dir.name}.json"
        ]

        for l_file in ligand_files:
            try:
                with open(l_file, 'r') as f:
                    data = json.load(f)
                
                # Check neighborhood residues
                neighborhood = data.get("neighborhood", [])
                for res in neighborhood:
                    # Tuple format: [auth_asym_id, auth_seq_id, auth_comp_id, master_index]
                    if len(res) < 4:
                        continue
                    
                    pdb_id = res[1]
                    ma_id = res[3]

                    if ma_id is not None and pdb_id != ma_id:
                        print(f"{struct_dir.name:<6} | {l_file.name.split('_')[1]:<10} | {res[2]:<8} | {pdb_id:<6} -> {ma_id:<8}")
                        found_count += 1

                # Check interaction participants
                interactions = data.get("interactions", [])
                for ix in interactions:
                    for part in ix.get("participants", []):
                        # Tuple format: [asym, seq, comp, atom, isLigand, master_index]
                        if len(part) < 6:
                            continue
                        
                        pdb_id = part[1]
                        ma_id = part[5]

                        if ma_id is not None and pdb_id != ma_id:
                            print(f"{struct_dir.name:<6} | {l_file.name.split('_')[1]:<10} | {part[2]:<8} | {pdb_id:<6} -> {ma_id:<8} (IX)")
                            found_count += 1
                            
            except Exception as e:
                # Silently skip malformed files during scan
                continue

    if found_count == 0:
        print("\nAll augmented residues currently have 1:1 PDB-to-MA mapping.")
    else:
        print(f"\nFound {found_count} discrepancies.")

if __name__ == "__main__":
    scan_augmented_ligands()