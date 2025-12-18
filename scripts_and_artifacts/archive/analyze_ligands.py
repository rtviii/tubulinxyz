import json
import csv
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, Any, List

# Configuration
DATA_DIR = Path("TUBETL_DATA")
OUTPUT_CSV = "tubulin_ligand_diversity.csv"

# Ligands to exclude from the "Interesting" list (Common biological buffers/ions)
# You can adjust this list later based on the output.
COMMON_NOISE = {
    "HOH", "DOD", # Water
    "SO4", "PO4", "CL", "NA", "MG", "K", "CA", "ZN", # Common Ions
    "GOL", "EDO", "PEG", "PGE", # Common Crystallization agents
    "GTP", "GDP", "G2P", # Endogenous Nucleotides (keep or remove depending on goal)
}

class LigandAnalyzer:
    def __init__(self):
        self.ligand_stats = defaultdict(lambda: {
            "name": "Unknown",
            "formula_weight": 0.0,
            "structure_count": 0,
            "total_instances": 0,
            "pdb_examples": set(),
            "is_drugbank_annotated": False,
            "drugbank_id": None
        })
        self.total_structures_processed = 0

    def process_file(self, file_path: Path):
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)

            rcsb_id = data.get("rcsb_id", "UNKNOWN")
            entities = data.get("entities", {})
            nonpolymers_instances = data.get("nonpolymers", [])

            # 1. Map Entity ID to Chemical Info
            # We need this because instances point to entities, and entities hold the metadata
            entity_lookup = {}

            for eid, entity in entities.items():
                if entity.get("type") == "non-polymer":
                    chem_id = entity.get("chemical_id")

                    # Check for DrugBank annotation in your schema
                    # Schema: nonpolymer_comp -> drugbank -> drugbank_container_identifiers -> drugbank_id
                    is_drug = False
                    db_id = None
                    if entity.get("nonpolymer_comp"):
                        try:
                            db_info = entity["nonpolymer_comp"].get("drugbank", {})
                            if db_info:
                                is_drug = True
                                # Try to grab ID if deeply nested
                                identifiers = db_info.get("drugbank_container_identifiers")
                                if identifiers:
                                    db_id = identifiers.get("drugbank_id")
                        except Exception:
                            pass

                    entity_lookup[eid] = {
                        "chem_id": chem_id,
                        "name": entity.get("chemical_name"),
                        "weight": entity.get("formula_weight"),
                        "is_drug": is_drug,
                        "db_id": db_id
                    }

            # 2. Count Instances
            # We track which ligands appear in THIS structure to increment structure_count only once per PDB
            ligands_in_this_structure = set()

            for instance in nonpolymers_instances:
                ent_id = instance.get("entity_id")
                info = entity_lookup.get(ent_id)

                if info:
                    chem_id = info["chem_id"]

                    # Update Stats
                    stat = self.ligand_stats[chem_id]
                    stat["name"] = info["name"]
                    stat["formula_weight"] = info["weight"] or 0.0
                    stat["total_instances"] += 1
                    stat["is_drugbank_annotated"] = info["is_drug"] or stat["is_drugbank_annotated"]
                    if info["db_id"]:
                         stat["drugbank_id"] = info["db_id"]

                    # Add example PDB (limit to first 5 to save memory in large scale)
                    if len(stat["pdb_examples"]) < 5:
                        stat["pdb_examples"].add(rcsb_id)

                    ligands_in_this_structure.add(chem_id)

            # Update structure counts
            for chem_id in ligands_in_this_structure:
                self.ligand_stats[chem_id]["structure_count"] += 1

            self.total_structures_processed += 1
            print(f"Processed {rcsb_id}...", end="\r")

        except Exception as e:
            print(f"\nError processing {file_path}: {e}")

    def report(self):
        print(f"\n\n--- Processing Complete ---")
        print(f"Total Structures: {self.total_structures_processed}")
        print(f"Total Unique Ligands Found: {len(self.ligand_stats)}")

        # Sort by structure count (descending)
        sorted_ligands = sorted(
            self.ligand_stats.items(),
            key=lambda item: item[1]['structure_count'],
            reverse=True
        )

        # Write to CSV
        with open(OUTPUT_CSV, 'w', newline='') as csvfile:
            fieldnames = [
                'Chemical_ID', 'Name', 'Category', 'Structure_Count',
                'Total_Instances', 'Avg_per_Structure', 'MW',
                'Is_DrugBank', 'DrugBank_ID', 'Example_PDBs'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for chem_id, data in sorted_ligands:
                # Basic Categorization Logic
                if chem_id in COMMON_NOISE:
                    category = "Common/Ion/Buffer"
                elif data['is_drugbank_annotated']:
                    category = "Drug/Lead"
                elif "GTP" in chem_id or "GDP" in chem_id:
                    category = "Nucleotide Analog"
                else:
                    category = "Unclassified"

                writer.writerow({
                    'Chemical_ID': chem_id,
                    'Name': data['name'],
                    'Category': category,
                    'Structure_Count': data['structure_count'],
                    'Total_Instances': data['total_instances'],
                    'Avg_per_Structure': round(data['total_instances'] / data['structure_count'], 2),
                    'MW': data['formula_weight'],
                    'Is_DrugBank': data['is_drugbank_annotated'],
                    'DrugBank_ID': data['drugbank_id'],
                    'Example_PDBs': ", ".join(list(data['pdb_examples']))
                })

        print(f"Results saved to {OUTPUT_CSV}")

def main():
    analyzer = LigandAnalyzer()

    # Recursively find all json files, excluding the sequence_ingestion ones
    all_files = list(DATA_DIR.glob("*/*.json"))
    structure_files = [f for f in all_files if "sequence_ingestion" not in f.name]

    if not structure_files:
        print("No JSON files found in TUBETL_DATA.")
        return

    for f in structure_files:
        analyzer.process_file(f)

    analyzer.report()

if __name__ == "__main__":
    main()
