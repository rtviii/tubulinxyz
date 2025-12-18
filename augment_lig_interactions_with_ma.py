import json
import os
import glob
from pathlib import Path
from lib.etl.constants import TUBETL_DATA
from lib.types import (
    LigandNeighborhood, 
    NeighborResidue, 
    InteractionParticipant
)

def augment_all_data():
    print(f"\n{'='*80}")
    print(f"CANONICAL LIGAND AUGMENTATION")
    print(f"Target Directory: {TUBETL_DATA}")
    print(f"{'='*80}\n")

    data_dir = Path(TUBETL_DATA)
    # Get all structure directories
    struct_dirs = sorted([d for d in data_dir.iterdir() if d.is_dir() and len(d.name) == 4])
    
    grand_stats = {"files_updated": 0, "ix_mapped": 0, "nb_mapped": 0}

    for s_dir in struct_dirs:
        rcsb_id = s_dir.name.upper()
        profile_path = s_dir / f"{rcsb_id}.json"
        ingestion_path = s_dir / "sequence_ingestion.json"
        
        if not profile_path.exists() or not ingestion_path.exists():
            print(f"[-] Skipping {rcsb_id}: Missing metadata (Profile or Ingestion).")
            continue

        print(f"[*] Processing {rcsb_id}...")

        try:
            with open(profile_path, 'r') as f:
                profile_dict = json.load(f)
            with open(ingestion_path, 'r') as f:
                ingestion = json.load(f)

            # 1. Map Chain -> MA Map
            chain_to_ma = {}
            for eid, ent_data in profile_dict.get('entities', {}).items():
                if eid in ingestion:
                    ma_map = ingestion[eid]['data']['observed_to_ma_map']
                    for chain in ent_data.get('pdbx_strand_ids', []):
                        chain_to_ma[chain] = ma_map

            # 2. Find every potential ligand JSON in the folder
            # Pattern: 1FFX_GTP_A.json
            pattern = str(s_dir / f"{rcsb_id}_*_*[!.json].json")
            candidate_files = [Path(f) for f in glob.glob(pattern) 
                              if not any(x in f for x in ["report", "ingestion", f"{rcsb_id}.json"])]

            if not candidate_files:
                print(f"    (No ligand files found for {rcsb_id})")
                continue

            for f_path in candidate_files:
                print(f"    > {f_path.name}", end=" ")
                
                with open(f_path, 'r') as f:
                    ligand_data = json.load(f)

                file_modified = False
                current_file_ix = 0
                current_file_nb = 0

                for entry in ligand_data:
                    # Validate existing data using your updated Pydantic model
                    nb = LigandNeighborhood.from_raw(entry)

                    # A. Augment Interactions
                    for ix in nb.interactions:
                        for part in ix.participants:
                            if not part.is_ligand:
                                ma_map = chain_to_ma.get(part.auth_asym_id)
                                if ma_map:
                                    idx = part.auth_seq_id - 1
                                    if 0 <= idx < len(ma_map):
                                        new_ma = ma_map[idx]
                                        # Only update if different or new
                                        if part.master_index != new_ma:
                                            part.master_index = new_ma
                                            file_modified = True
                                            current_file_ix += 1

                    # B. Augment Neighborhoods
                    updated_nb_list = []
                    for raw_nb in entry.get('neighborhood', []):
                        n_res = NeighborResidue.from_tuple(raw_nb)
                        ma_map = chain_to_ma.get(n_res.auth_asym_id)
                        if ma_map:
                            idx = n_res.auth_seq_id - 1
                            if 0 <= idx < len(ma_map):
                                new_ma = ma_map[idx]
                                if n_res.master_index != new_ma:
                                    n_res.master_index = new_ma
                                    file_modified = True
                                    current_file_nb += 1
                        updated_nb_list.append(n_res)
                    
                    # C. Sync dictionaries back
                    if file_modified:
                        entry['interactions'] = [
                            {
                                "type": ix.type,
                                "participants": [p.to_tuple() for p in ix.participants]
                            } for ix in nb.interactions
                        ]
                        entry['neighborhood'] = [n.to_tuple() for n in updated_nb_list]

                if file_modified:
                    with open(f_path, 'w') as f:
                        json.dump(ligand_data, f, indent=2)
                    print(f"✓ (Mapped: {current_file_ix} ix, {current_file_nb} nb)")
                    grand_stats["files_updated"] += 1
                    grand_stats["ix_mapped"] += current_file_ix
                    grand_stats["nb_mapped"] += current_file_nb
                else:
                    print("- (Already augmented or no mapping available)")

        except Exception as e:
            print(f"\n    [!] ERROR processing {rcsb_id}: {e}")

    print(f"\n{'='*80}")
    print(f"AUGMENTATION SUMMARY")
    print(f"  Total Files Overwritten:   {grand_stats['files_updated']}")
    print(f"  Total Interactions Tagged: {grand_stats['ix_mapped']}")
    print(f"  Total Neighborhood Tagged: {grand_stats['nb_mapped']}")
    print(f"{'='*80}\n")

if __name__ == "__main__":
    augment_all_data()
