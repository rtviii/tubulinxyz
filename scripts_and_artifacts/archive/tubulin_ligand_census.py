import json
from pprint import pprint
import sys
from pathlib import Path
from collections import defaultdict
import textwrap

# --- CONFIGURATION ---
DATA_DIR = Path("TUBETL_DATA")

# 1. THE IGNORE LIST
# These will be excluded from the main census to let actual drugs/leads shine.
IGNORE_IDS = {
    # Water
    "HOH", "DOD", 
    # Common Ions
    "MG", "CL", "NA", "K", "CA", "ZN", "MN", "SO4", "PO4", "NH4", "LI", "NI", "CO", 
    # Buffers & Crystallization Agents
    "MES", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "FMT", "ACT", "DMS", "IMD", "CIT", "ACY", "TRS", "BMA", "EPE",
    # Endogenous Nucleotides (Remove if you want to see them)
    "GTP", "GDP", "ATP", "ADP", 
    # High-abundance modified nucleotides often used as structural tools (Motor protein bound)
    "ANP", "ACP", "GSP", "G2P", "AGS", "GP2", "GCP"
}

class LigandCensus:
    def __init__(self):
        # Stats for the "Interesting" ligands
        self.stats = defaultdict(lambda: {
            "name": "Unknown",
            "mw": 0.0,
            "pdb_count": 0,
            "total_instances": 0,
            "description": "N/A",
            "db_id": None
        })
        # Stats for the "Excluded" ligands (just to capture their names)
        self.excluded_info = defaultdict(lambda: "Unknown")
        
        self.processed_count = 0
        self.skipped_count = 0

    def _update_name(self, current_name, new_name):
        """Helper to pick the most descriptive name found across files."""
        if not new_name:
            return current_name
        if current_name == "Unknown":
            return new_name
        # Heuristic: Longer names are usually more descriptive (e.g. "TAXOL" vs "PACLITAXEL")
        if len(new_name) > len(current_name):
            return new_name
        return current_name

    def process_file(self, file_path: Path):
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Use a set to count each ligand only once per structure (PDB Count)
            structure_ligands = set()
            
            # Map Entity -> Info
            entity_info = {}
            for eid, entity in data.get("entities", {}).items():
                if entity.get("type") == "non-polymer":
                    
                    # Extract DrugBank ID if available
                    db_id = None
                    desc = None
                    if entity.get("nonpolymer_comp"):
                        # Try DrugBank
                        db = entity["nonpolymer_comp"].get("drugbank", {})
                        if db:
                            # Safely navigate nested dicts
                            identifiers = db.get("drugbank_container_identifiers")
                            if identifiers:
                                db_id = identifiers.get("drugbank_id")
                            
                            info = db.get("drugbank_info")
                            if info:
                                desc = info.get("description")
                    
                    # Fallback description
                    if not desc:
                        desc = entity.get("pdbx_description")

                    entity_info[eid] = {
                        "chem_id": entity.get("chemical_id"),
                        "name": entity.get("chemical_name"),
                        "mw": entity.get("formula_weight", 0.0),
                        "desc": desc,
                        "db_id": db_id
                    }

            # Count Instances
            for instance in data.get("nonpolymers", []):
                eid = instance.get("entity_id")
                info = entity_info.get(eid)
                
                if info:
                    chem_id = info["chem_id"]
                    chem_name = info["name"]
                    
                    # --- PATH A: EXCLUDED LIGANDS ---
                    if chem_id in IGNORE_IDS:
                        self.excluded_info[chem_id] = self._update_name(self.excluded_info[chem_id], chem_name)
                        continue

                    # --- PATH B: CENSUS LIGANDS ---
                    entry = self.stats[chem_id]
                    
                    # Update Metadata
                    entry["name"] = self._update_name(entry["name"], chem_name)
                    
                    if info["mw"]: 
                        entry["mw"] = info["mw"]
                    
                    # Keep longest description found
                    if info["desc"] and (entry["description"] == "N/A" or len(str(info["desc"])) > len(str(entry["description"]))):
                        entry["description"] = info["desc"]
                    
                    if info["db_id"]:
                        entry["db_id"] = info["db_id"]
                    
                    entry["total_instances"] += 1
                    structure_ligands.add(chem_id)

            # Update Structure Counts (only for non-excluded)
            for chem_id in structure_ligands:
                self.stats[chem_id]["pdb_count"] += 1
            
            self.processed_count += 1
            print(f"Scanning... {self.processed_count} structures processed.", end="\r")

        except Exception as e:
            pass

    def report(self):
        print("\n" * 2)
        print("="*120)
        print(f"  TUBULIN PHARMACO-LIGAND CENSUS (Top 50)")
        print(f"  Scanned {self.processed_count} structures.")
        print("="*120)

        # Sort: Primary = PDB Count, Secondary = Total Instances
        sorted_ligands = sorted(
            self.stats.items(), 
            key=lambda x: (x[1]['pdb_count'], x[1]['total_instances']), 
            reverse=True
        )[:50]

        if not sorted_ligands:
            print("No relevant ligands found.")
            return

        # Header
        header = f"{'ID':<5} | {'PDBs':<5} | {'MW (Da)':<8} | {'Name':<22} | {'Annotation / Description / Link'}"
        print(header)
        print("-" * 120)

        for chem_id, data in sorted_ligands:
            pdbs = str(data['pdb_count'])
            mw = f"{data['mw']:.1f}" if data['mw'] else "N/A"
            
            # Smart Truncation for Name
            name_raw = data['name'] or "Unknown"
            name = name_raw[:20] + ".." if len(name_raw) > 22 else name_raw
            
            # Build Annotation Column
            annot_parts = []
            
            # 1. Add Link if available
            if data['db_id']:
                link = f"https://go.drugbank.com/drugs/{data['db_id']}"
                annot_parts.append(f"[{link}]")
            
            # 2. Add Description
            desc_text = str(data['description'] or "").replace("\n", " ")
            if desc_text and desc_text != "N/A":
                # Truncate description to fit remaining terminal width (approx 60 chars)
                max_desc_len = 65
                if len(desc_text) > max_desc_len:
                    desc_text = desc_text[:max_desc_len-3] + "..."
                annot_parts.append(desc_text)
            
            annotation = " ".join(annot_parts)
            if not annotation:
                annotation = "-"

            print(f"{chem_id:<5} | {pdbs:<5} | {mw:<8} | {name:<22} | {annotation}")

        print("-" * 120)
        
        # --- EXCLUDED LIST DISCLAIMER ---
        print("\nDisclaimer: The following common ions, buffers, and reagents were EXCLUDED from this census:")
        
        # Build list of "ID (Name)"
        excluded_display = []
        for eid in sorted(list(IGNORE_IDS)):
            name = self.excluded_info.get(eid)
            if name and name != "Unknown":
                excluded_display.append(f"{eid} ({name})")
            else:
                excluded_display.append(eid)
                
        # Print wrapped text
        # wrapper = textwrap.TextWrapper(width=120, initial_indent="  ", subsequent_indent="  ")
        pprint(excluded_display)
        print("="*120)

def main():
    census = LigandCensus()
    
    files = list(DATA_DIR.glob("*/*.json"))
    structure_files = [f for f in files if "sequence_ingestion" not in f.name]
    
    if not structure_files:
        print(f"No JSON files found in {DATA_DIR}")
        return

    for f in structure_files:
        census.process_file(f)
    
    census.report()

if __name__ == "__main__":
    main()