import ssl
import urllib.request
import tempfile
import warnings
import os
import json
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
from dataclasses import dataclass, field, asdict
from collections import Counter

# Biopython Imports
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio import BiopythonWarning, AlignIO

# Suppress PDB construction warnings
warnings.simplefilter('ignore', BiopythonWarning)

# ==========================================
# 1. CONSTANTS & DICTIONARIES
# ==========================================

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

MOD_RES_MAP = {
    'MSE': 'M', 'CSX': 'C', 'CSO': 'C', 'SEP': 'S', 
    'TPO': 'T', 'PTR': 'Y', 'KCX': 'K', 'LLP': 'K', 
    'PCA': 'E', 'HIC': 'H', 
}

# ==========================================
# 2. DATA CLASSES
# ==========================================

@dataclass
class MutationEntry:
    ma_position: int  # 1-based index in Master Alignment
    wild_type: str    # Consensus residue
    observed: str     # Actual residue in target
    pdb_auth_id: int  # The auth_seq_id of the mutation

@dataclass
class ProcessedChain:
    pdb_id: str
    chain_id: str
    tubulin_class: str
    sequence: str
    
    # Mappings
    # Index = MA Position (0-based), Value = Auth_Seq_ID or -1 (Missing)
    ma_to_auth_map: List[int]
    
    # Index = Observed Seq Index (0-based), Value = MA Position (1-based) or -2 (Insertion)
    observed_to_ma_map: List[int]

    mutations: List[MutationEntry]
    
    # Stats
    stats: Dict[str, Any]

# ==========================================
# 3. SERVICES
# ==========================================

class TubulinStructureParser:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
        self.conversion_map = {}
        
        # Add Standard
        for code3, props in AMINO_ACIDS.items():
            self.conversion_map[code3] = props["one_letter_code"]
            
        # Add Modified
        for code3, code1 in MOD_RES_MAP.items():
            self.conversion_map[code3] = code1

    def fetch_cif_to_temp(self, pdb_id: str) -> Optional[Path]:
        pdb_id = pdb_id.lower()
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        
        try:
            tf = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
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

        try:
            structure = self.parser.get_structure('temp', str(file_path))
            try:
                model = structure[0]
            except (KeyError, IndexError):
                model = next(iter(structure))

            if chain_id not in model:
                return "", []

            chain = model[chain_id]
            sequence_buffer = []
            auth_id_buffer = []

            for residue in chain:
                hetero_flag, seq_num, ins_code = residue.id
                res_name = residue.resname.upper()

                if res_name not in self.conversion_map:
                    continue

                if hetero_flag != ' ' and not res_name in MOD_RES_MAP:
                    pass

                aa = self.conversion_map[res_name]
                sequence_buffer.append(aa)
                auth_id_buffer.append(seq_num)

            return "".join(sequence_buffer), auth_id_buffer

        except Exception as e:
            print(f"Error parsing structure {file_path}: {e}")
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


class ConsensusCalculator:
    def __init__(self, master_profile_path: str):
        self.path = Path(master_profile_path)
        self.consensus_sequence: List[str] = []
        self._calculate()

    def _calculate(self):
        if not self.path.exists():
            raise FileNotFoundError(f"Master alignment not found at {self.path}")

        alignment = AlignIO.read(str(self.path), "fasta")
        length = alignment.get_alignment_length()
        
        consensus = []
        for i in range(length):
            column = alignment[:, i]
            residues = [r for r in column if r != '-' and r != '.']
            
            if not residues:
                consensus.append('-')
                continue

            most_common = Counter(residues).most_common(1)[0][0]
            consensus.append(most_common)
            
        self.consensus_sequence = consensus

    def get_residue_at(self, index_0_based: int) -> str:
        if 0 <= index_0_based < len(self.consensus_sequence):
            return self.consensus_sequence[index_0_based]
        return '?'


class TubulinAlignmentMapper:
    def __init__(self, master_profile_path: str, muscle_binary: str):
        self.master_profile_path = Path(master_profile_path)
        self.muscle_binary = muscle_binary
        
        if not Path(self.muscle_binary).exists():
            print(f"WARNING: MUSCLE binary not found at {self.muscle_binary}")

    def align_sequence(self, sequence_id: str, sequence: str) -> Tuple[str, str]:
        """
        Returns (aligned_master_sequence, aligned_target_sequence).
        We need BOTH to determine where gaps are in the master vs the target.
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as seq_file:
            seq_file.write(f">{sequence_id}\n{sequence}\n")
            seq_temp = seq_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.msa', delete=False) as output_file:
            try:
                cmd = [
                    self.muscle_binary, '-profile',
                    '-in1', str(self.master_profile_path),
                    '-in2', seq_temp,
                    '-out', output_file.name
                ]
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Parse to get both aligned sequences
                return self._extract_aligned_pair(output_file.name, sequence_id)
                
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"MUSCLE Error: {e.stderr.decode()}")
            finally:
                Path(seq_temp).unlink(missing_ok=True)
                Path(output_file.name).unlink(missing_ok=True)



    def align_sequence_with_mapping(
        self, 
        sequence_id: str, 
        sequence: str, 
        residue_numbers: Optional[List[int]] = None
    ) -> Tuple[str, List[int]]:
        """
        Legacy method for API compatibility.
        Returns (aligned_sequence, mapping) where mapping maps each aligned position 
        to auth_seq_id or -1 for gaps.
        """
        # Validate input
        if residue_numbers and len(sequence) != len(residue_numbers):
            raise ValueError(
                f"Length mismatch: Sequence len {len(sequence)} vs IDs len {len(residue_numbers)}."
            )

        # Use the new align_sequence to get both aligned strings
        aligned_master, aligned_target = self.align_sequence(sequence_id, sequence)
        
        # Build the mapping (same logic as old code)
        mapping = []
        residue_ptr = 0
        
        for char in aligned_target:
            if char == '-':
                mapping.append(-1)
            else:
                if residue_numbers:
                    mapping.append(residue_numbers[residue_ptr])
                else:
                    mapping.append(residue_ptr + 1)
                residue_ptr += 1
        
        return aligned_target, mapping

    def _extract_aligned_pair(self, output_file: str, sequence_id: str) -> Tuple[str, str]:
        alignment = AlignIO.read(output_file, "fasta")
        
        # Find the target record
        target_record = next((r for r in alignment if r.id == sequence_id), None)
        if not target_record:
            raise ValueError(f"Sequence {sequence_id} not found in alignment output")
            
        aligned_target = str(target_record.seq)
        
        # To find gaps in the Master, we take the first sequence of the profile.
        # (In MUSCLE profile mode, the preserved profile sequences usually come first)
        # This acts as our "Master Coordinate Ruler" in the alignment block.
        aligned_master = str(alignment[0].seq)
        
        return aligned_master, aligned_target



class TubulinIngestor:
    def __init__(self, master_profile: str, muscle_binary: str):
        self.parser = TubulinStructureParser()
        self.consensus = ConsensusCalculator(master_profile)
        self.mapper = TubulinAlignmentMapper(master_profile, muscle_binary)
        
        self.ref_len = len(self.consensus.consensus_sequence)

    def process_chain(self, pdb_id: str, chain_id: str, t_class: str) -> ProcessedChain:
        # 1. Fetch & Extract
        cif_path = self.parser.fetch_cif_to_temp(pdb_id)
        if not cif_path:
            raise ValueError(f"Could not fetch {pdb_id}")
        
        try:
            obs_seq, auth_ids = self.parser.get_observed_data(cif_path, chain_id)
        finally:
            if cif_path.exists():
                os.unlink(cif_path)

        if not obs_seq:
            raise ValueError(f"No sequence data for {pdb_id} chain {chain_id}")

        # 2. Align (Get both strings to see gaps)
        seq_id = f"{pdb_id}_{chain_id}"
        aln_master, aln_target = self.mapper.align_sequence(seq_id, obs_seq)

        # 3. Construct Mappings & Detect Mutations
        ma_to_auth = [-1] * self.ref_len  # Default -1 (Missing)
        observed_to_ma = []               # Will contain MA indices or -2 (Insertion)
        mutations = []
        
        # Pointers
        master_index = 0       # 0-based index into the Consensus Sequence
        target_res_index = 0   # Index into obs_seq / auth_ids

        # Iterate the alignment columns
        # Both strings are same length from MUSCLE
        for m_char, t_char in zip(aln_master, aln_target):
            
            # CASE 1: GAP IN MASTER (Insertion in Target)
            if m_char == '-':
                if t_char != '-':
                    # Target has residue, Master has gap -> Insertion
                    observed_to_ma.append(-2)
                    target_res_index += 1
                # If both are gaps, ignore (alignment artifact)
                
            # CASE 2: GAP IN TARGET (Deletion/Missing in Target)
            elif t_char == '-':
                # Master has residue, Target has gap -> Missing
                # ma_to_auth[master_index] remains -1
                master_index += 1
                
            # CASE 3: MATCH (Both have residues)
            else:
                # Map indices
                current_auth_id = auth_ids[target_res_index]
                
                # Store mapping
                if master_index < self.ref_len:
                    ma_to_auth[master_index] = current_auth_id
                    observed_to_ma.append(master_index + 1) # Store 1-based MA index
                    
                    # Check Mutation
                    wild_type = self.consensus.get_residue_at(master_index)
                    
                    # If wild_type is a gap in consensus (rare), skip mutation check
                    if wild_type != '-' and wild_type != '?' and t_char != wild_type:
                        mutations.append(MutationEntry(
                            ma_position=master_index + 1,
                            wild_type=wild_type,
                            observed=t_char,
                            pdb_auth_id=current_auth_id
                        ))
                
                master_index += 1
                target_res_index += 1

        # 4. Return Result
        stats = {
            "alignment_length": len(aln_target),
            "ma_coverage": len([x for x in ma_to_auth if x != -1]),
            "insertions": observed_to_ma.count(-2),
            "total_mutations": len(mutations)
        }

        return ProcessedChain(
            pdb_id=pdb_id,
            chain_id=chain_id,
            tubulin_class=t_class,
            sequence=obs_seq,
            ma_to_auth_map=ma_to_auth,
            observed_to_ma_map=observed_to_ma,
            mutations=mutations,
            stats=stats
        )


def main():
    # --- CONFIGURATION ---
    MUSCLE_BINARY = "./muscle3.8.1" # <--- ADJUST THIS PATH
    MASTER_PROFILE = "./data/alpha_tubulin/alpha_tubulin.afasta" # <--- ADJUST THIS PATH
    
    OUTPUT_FILE = "ingestion_results.json"

    targets = [
        {"pdb": "1JFF", "chain": "A", "class": "Alpha"},
        {"pdb": "1JFF", "chain": "B", "class": "Beta"},
        # {"pdb": "7SJ7", "chain": "A", "class": "Alpha"},
    ]
    
    print("Starting Ingestion Pipeline...")
    ingestor = TubulinIngestor(MASTER_PROFILE, MUSCLE_BINARY)
    
    results = []
    
    for target in targets:
        print(f"Processing {target['pdb']} Chain {target['chain']}...")
        try:
            res = ingestor.process_chain(target['pdb'], target['chain'], target['class'])
            
            # Convert to dict for JSON serialization
            results.append(asdict(res))
            print(f"  -> Success. Mutations: {res.stats['total_mutations']}, Insertions: {res.stats['insertions']}")
            
        except Exception as e:
            print(f"  -> FAILED: {e}")

    # Save to File
    with open(OUTPUT_FILE, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nDone. Results saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()