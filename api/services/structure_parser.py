import ssl
import urllib.request
import tempfile
import warnings
import os
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio import BiopythonWarning

# Suppress PDB construction warnings
warnings.simplefilter('ignore', BiopythonWarning)

# --- Configuration Data ---

AMINO_ACIDS = {
    "ALA": {"one_letter_code": "A", "charge": 0},
    "ARG": {"one_letter_code": "R", "charge": 1},
    "ASN": {"one_letter_code": "N", "charge": 0},
    "ASP": {"one_letter_code": "D", "charge": -1},
    "CYS": {"one_letter_code": "C", "charge": 0},
    "GLU": {"one_letter_code": "E", "charge": -1},
    "GLN": {"one_letter_code": "Q", "charge": 0},
    "GLY": {"one_letter_code": "G", "charge": 0},
    "HIS": {"one_letter_code": "H", "charge": 0},
    "ILE": {"one_letter_code": "I", "charge": 0},
    "LEU": {"one_letter_code": "L", "charge": 0},
    "LYS": {"one_letter_code": "K", "charge": 1},
    "MET": {"one_letter_code": "M", "charge": 0},
    "PHE": {"one_letter_code": "F", "charge": 0},
    "PRO": {"one_letter_code": "P", "charge": 0},
    "SER": {"one_letter_code": "S", "charge": 0},
    "THR": {"one_letter_code": "T", "charge": 0},
    "TRP": {"one_letter_code": "W", "charge": 0},
    "TYR": {"one_letter_code": "Y", "charge": 0},
    "VAL": {"one_letter_code": "V", "charge": 0},
}

# Common Tubulin/PDB modified residues mapping to parents
MOD_RES_MAP = {
    'MSE': 'M', 'CSX': 'C', 'CSO': 'C', 'SEP': 'S', 
    'TPO': 'T', 'PTR': 'Y', 'KCX': 'K', 'LLP': 'K', 
    'PCA': 'E', 'HIC': 'H', 
    # Nucleotides might appear as ligands, but usually we filter them out 
    # unless you want to track them. Leaving them out of AA map for now.
}

class TubulinStructureParser:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
        
        # Build a fast lookup table: 3-letter -> 1-letter
        # Merges your standard AAs with the modified residue map
        self.conversion_map = {}
        
        # Add Standard
        for code3, props in AMINO_ACIDS.items():
            self.conversion_map[code3] = props["one_letter_code"]
            
        # Add Modified
        for code3, code1 in MOD_RES_MAP.items():
            self.conversion_map[code3] = code1

    def fetch_cif_to_temp(self, pdb_id: str) -> Optional[Path]:
        """
        Downloads the .cif file from RCSB PDB to a temporary file.
        Returns path to the temp file.
        """
        pdb_id = pdb_id.lower()
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        
        try:
            # Create a named temp file that persists (delete=False)
            tf = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
            
            # Bypass SSL verification for reliability in dev environments
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE

            with urllib.request.urlopen(url, context=ctx) as response:
                tf.write(response.read())
            
            tf.close()
            return Path(tf.name)
            
        except Exception as e:
            print(f"Failed to download PDB {pdb_id}: {e}")
            if 'tf' in locals():
                tf.close()
                os.unlink(tf.name)
            return None

    def get_observed_data(self, file_path: Path, chain_id: str) -> Tuple[str, List[int]]:
        """
        Parses a CIF file path and returns the OBSERVED sequence and Auth Seq IDs.
        """
        try:
            structure = self.parser.get_structure('temp', str(file_path))
            
            # Robust model selection (handle cases where model is a list or dict)
            try:
                model = structure[0]
            except (KeyError, IndexError):
                model = next(iter(structure))

            if chain_id not in model:
                # If exact match fails, return empty (verification will fail safely)
                return "", []

            chain = model[chain_id]
            sequence_buffer = []
            auth_id_buffer = []

            for residue in chain:
                # residue.id = (hetero_flag, seq_num, insertion_code)
                hetero_flag, seq_num, ins_code = residue.id
                res_name = residue.resname.upper()

                # 1. Lookup: Is this a known AA (standard or modified)?
                if res_name not in self.conversion_map:
                    continue

                # 2. Filter: 
                # hetero_flag is ' ' for standard residues.
                # hetero_flag is 'H_XXX' for modified residues (e.g. H_MSE).
                # We want to exclude purely water ('W') or ligands that aren't in our map.
                # Since we checked `res_name in self.conversion_map`, we implicitly 
                # filtered out unknown ligands. We just need to be careful of conflicts.
                
                # Double check to ensure we don't pick up a ligand named "ALA" (unlikely but possible)
                if hetero_flag != ' ' and not res_name in MOD_RES_MAP:
                    # If it's hetero but not in our explicit modified map, it might be a ligand 
                    # that happens to share a 3-letter code, or just noise. 
                    # However, strictly adhering to your dictionary is safest.
                    # If it's in AMINO_ACIDS, we assume it's part of the polymer.
                    pass

                aa = self.conversion_map[res_name]
                
                sequence_buffer.append(aa)
                auth_id_buffer.append(seq_num)

            return "".join(sequence_buffer), auth_id_buffer
            
        except Exception as e:
            print(f"Error parsing structure {file_path}: {e}")
            # print(traceback.format_exc()) # Uncomment for deep debugging
            return "", []

    def verify_integrity(self, pdb_id: str, chain_id: str, 
                        fe_sequence: str, fe_auth_ids: List[int]) -> Dict:
        """
        Orchestrates the download, parse, and comparison.
        """
        temp_path = self.fetch_cif_to_temp(pdb_id)
        if not temp_path:
            return {"status": "skipped", "reason": "download_failed"}

        try:
            be_seq, be_ids = self.get_observed_data(temp_path, chain_id)
        finally:
            if temp_path.exists():
                os.unlink(temp_path)

        # Comparison Logic
        seq_match = fe_sequence == be_seq
        id_match = fe_auth_ids == be_ids
        
        mismatches = []
        if not seq_match or not id_match:
            limit = min(len(fe_sequence), len(be_seq))
            for i in range(limit):
                if fe_sequence[i] != be_seq[i] or fe_auth_ids[i] != be_ids[i]:
                    mismatches.append({
                        "pos": i,
                        "fe": f"{fe_sequence[i]}:{fe_auth_ids[i]}",
                        "be": f"{be_seq[i]}:{be_ids[i]}"
                    })
                    if len(mismatches) > 3: break 

        return {
            "status": "success",
            "match": seq_match and id_match,
            "details": {
                "seq_match": seq_match,
                "id_match": id_match,
                "fe_len": len(fe_sequence),
                "be_len": len(be_seq),
                "mismatches": mismatches
            }
        }