import os
from pprint import pprint
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import gemmi

from api.main import MUSCLE_BINARY
from etl.constants import TUBETL_DATA  # Better mmCIF parser than Biopython


@dataclass
class PolymerEntity:
    """Represents a unique polymer entity from mmCIF"""

    entity_id: str
    sequence : str  # pdbx_seq_one_letter_code_can
    
    
@dataclass
class AlignmentMappings:
    """Stores bidirectional mappings between sequences"""
    master_to_target: List[int]                   # MA index -> target index (-1 if gap in target)
    target_to_master: List[int]                   # target index -> MA index (-1 if insertion)
    mutations       : List[Tuple[int, str, str]]  # (MA_index, MA_residue, target_residue)

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def parse_mmcif_polymers(mmcif_path: str) -> List[PolymerEntity]:
    """
    Parse mmCIF file and extract unique polymer entities.
    Uses _entity_poly.pdbx_seq_one_letter_code_can as the canonical sequence.
    """
    mmcif_dict = MMCIF2Dict(mmcif_path)
    
    polymers = []
    
    # Get entity_poly data
    if '_entity_poly.entity_id' in mmcif_dict:
        entity_ids = mmcif_dict['_entity_poly.entity_id']
        sequences = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code_can']
        
        # Handle case where there's only one entity (not a list)
        if not isinstance(entity_ids, list):
            entity_ids = [entity_ids]
            sequences = [sequences]
        
        for entity_id, sequence in zip(entity_ids, sequences):
            # Clean up sequence (remove newlines and spaces)
            sequence = sequence.replace('\n', '').replace(' ', '')
            
            if sequence and sequence != '?':
                polymers.append(PolymerEntity(entity_id=entity_id, sequence=sequence))
    
    return polymers

def parse_master_alignment(msa_path: str) -> Tuple[str, str]:
    """
    Parse the master alignment FASTA file.
    Returns (sequence_id, aligned_sequence)
    Assumes single sequence in the file (the master alignment).
    """
    with open(msa_path, 'r') as f:
        lines = f.readlines()
    
    seq_id = lines[0].strip().lstrip('>')
    sequence = ''.join(line.strip() for line in lines[1:])
    
    return seq_id, sequence

def align_to_master_with_muscle(
    target_seq: str,
    master_aln_path: str,  # Changed: pass file path instead of sequence
    muscle_path: str = MUSCLE_BINARY
) -> Tuple[str, str]:
    """
    Align target sequence to master alignment using MUSCLE -profile.
    Returns (aligned_master, aligned_target).
    
    Args:
        target_seq: Unaligned target sequence
        master_aln_path: Path to master alignment FASTA file
        muscle_path: Path to muscle executable
    """
    target_fasta = f">target\n{target_seq}\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as target_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as output_file:
        
        target_file.write(target_fasta)
        target_file.flush()
        
        try:
            # Run MUSCLE with -profile option
            # -in1 is the master alignment (multiple sequences)
            # -in2 is the target sequence
            cmd = [
                muscle_path,
                '-profile',
                '-in1', master_aln_path,  # Use the file path directly
                '-in2', target_file.name,
                '-out', output_file.name
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse output alignment
            with open(output_file.name, 'r') as f:
                lines = f.readlines()
            
            # Extract sequences - the output will have all master sequences + target
            # We want the FIRST master sequence and the target
            sequences = {}
            current_id = None
            current_seq = []
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_id:
                sequences[current_id] = ''.join(current_seq)
            
            # Get the first sequence (master) and target
            seq_ids = list(sequences.keys())
            master_seq = sequences[seq_ids[0]]  # First sequence from master alignment
            target_seq = sequences['target']
            
            return master_seq, target_seq
            
        finally:
            Path(target_file.name).unlink(missing_ok=True)
            Path(output_file.name).unlink(missing_ok=True)


def parse_alignment_mappings(
    aligned_master: str,
    aligned_target: str,
    original_master_seq: Optional[str] = None
) -> AlignmentMappings:
    """
    Parse aligned sequences to create bidirectional mappings.
    
    Args:
        aligned_master: Master alignment with gaps
        aligned_target: Target sequence aligned to master (with gaps)
        original_master_seq: Original ungapped master (if None, derived from aligned_master)
    
    Returns:
        AlignmentMappings with master_to_target, target_to_master, and mutations
    """
    if original_master_seq is None:
        original_master_seq = aligned_master.replace('-', '')
    
    master_to_target = []
    target_to_master = []
    mutations = []
    
    master_pos = 0  # Position in original (ungapped) master
    target_pos = 0  # Position in original (ungapped) target
    
    for i, (m_char, t_char) in enumerate(zip(aligned_master, aligned_target)):
        m_is_gap = (m_char == '-')
        t_is_gap = (t_char == '-')
        
        if not m_is_gap and not t_is_gap:
            # Both have residues - they align
            master_to_target.append(target_pos)
            target_to_master.append(master_pos)
            
            # Check for mutation
            if m_char != t_char:
                mutations.append((master_pos, m_char, t_char))
            
            master_pos += 1
            target_pos += 1
            
        elif not m_is_gap and t_is_gap:
            # Master has residue, target has gap (deletion in target)
            master_to_target.append(-1)
            master_pos += 1
            
        elif m_is_gap and not t_is_gap:
            # Target has residue, master has gap (insertion in target)
            target_to_master.append(-2)
            target_pos += 1
            
        # Both gaps should not happen in a proper alignment, but skip if it does
    
    return AlignmentMappings(
        master_to_target=master_to_target,
        target_to_master=target_to_master,
        mutations=mutations
    )

def align_polymer_to_master(
  polymer                                   : PolymerEntity,
  master_aln_path                           : str,                # Changed to path
  muscle_path                               : str =MUSCLE_BINARY
) -> Tuple[PolymerEntity, AlignmentMappings]: 
    """
    Align a single polymer to the master alignment.
    Wrapper function for parallel processing.
    """
    aligned_master, aligned_target = align_to_master_with_muscle(
        polymer.sequence,
        master_aln_path,  # Pass path instead of sequence
        muscle_path
    )
    
    mappings = parse_alignment_mappings(aligned_master, aligned_target)
    
    return polymer, mappings

def process_polymers_parallel(
    polymers: List[PolymerEntity],
    master_aln_path: str,  # Changed to path
    muscle_path: str = MUSCLE_BINARY,
    max_workers: int = 4
) -> Dict[str, AlignmentMappings]:
    """
    Process multiple polymers in parallel.
    
    Returns:
        Dictionary mapping entity_id to AlignmentMappings
    """
    results = {}
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(align_polymer_to_master, polymer, master_aln_path, muscle_path): polymer
            for polymer in polymers
        }
        
        for future in as_completed(futures):
            polymer, mappings = future.result()
            results[polymer.entity_id] = mappings
            print(f"Completed alignment for entity {polymer.entity_id}")
    
    return results


from etl.assets import TubulinStructureAssets
from models.types_tubulin import TubulinFamily
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from typing import Dict

def get_observed_sequence_from_structure(cif_path: str, chain_id: str) -> str:
    """
    Extract the actual observed sequence from the structure file.
    Only includes residues that have at least one atom coordinate.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    
    # Get the first model (index 0)
    model = structure[0]
    
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in structure")
    
    chain = model[chain_id]
    sequence = []
    
    for residue in chain:
        # Skip heteroatoms (water, ligands, etc.)
        if residue.id[0] != ' ':
            continue
        
        resname = residue.resname
        try:
            # Convert 3-letter code to 1-letter
            sequence.append(protein_letters_3to1[resname])
        except KeyError:
            # Unknown residue, use X
            sequence.append('X')
    
    return ''.join(sequence)

def get_alpha_tubulins_from_structure(rcsb_id: str) -> List[PolymerEntity]:
    """
    Load profile to identify alpha tubulins, but get sequences from actual structure.
    """
    try:
        assets = TubulinStructureAssets(rcsb_id)
        profile = assets.profile()
        cif_path = assets.paths.cif
        
        if not os.path.exists(cif_path):
            print(f"CIF file not found: {cif_path}")
            return []
        
        alpha_tubulins = []
        for protein in profile.proteins:
            if protein.family == TubulinFamily.ALPHA:
                # Use auth_asym_id as the chain identifier
                chain_id = protein.auth_asym_id
                
                try:
                    observed_seq = get_observed_sequence_from_structure(cif_path, chain_id)
                    print(f"Entity {protein.entity_id}, Chain {chain_id}: {len(observed_seq)} residues observed")
                    
                    alpha_tubulins.append(
                        PolymerEntity(
                            entity_id=protein.entity_id,
                            sequence=observed_seq  # Actual coordinates
                        )
                    )
                except Exception as e:
                    print(f"Error getting observed sequence for chain {chain_id}: {e}")
        
        return alpha_tubulins
        
    except FileNotFoundError:
        print(f"Profile not found for {rcsb_id}")
        return []
    except Exception as e:
        print(f"Error loading profile for {rcsb_id}: {e}")
        return []

def get_alpha_tubulins_from_profile(rcsb_id: str) -> List[PolymerEntity]:
    """
    Load a structure profile and extract only alpha tubulin proteins.
    Uses the observed sequence (what's actually in the structure).
    """
    try:
        assets = TubulinStructureAssets(rcsb_id)
        profile = assets.profile()
        
        alpha_tubulins = []
        for protein in profile.proteins:
            if protein.family == TubulinFamily.ALPHA:
                alpha_tubulins.append(
                    PolymerEntity(
                        entity_id=protein.entity_id,
                        sequence=protein.entity_poly_seq_one_letter_code  # Observed, not canonical
                    )
                )
        
        return alpha_tubulins
        
    except FileNotFoundError:
        print(f"Profile not found for {rcsb_id}")
        return []
    except Exception as e:
        print(f"Error loading profile for {rcsb_id}: {e}")
        return []


def get_observed_sequence_and_numbering(cif_path: str, chain_id: str) -> Tuple[str, List[int]]:
    """
    Extracts the observed sequence and the corresponding auth_seq_ids.
    Strictly filters for protein residues only (skipping waters/ligands).
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    model = structure[0]
    
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found")
    
    chain = model[chain_id]
    sequence_chars = []
    auth_seq_ids = []
    
    for residue in chain:
        # 1. Filter Heteroatoms (Water, Ligands)
        # residue.id is a tuple: (hetero_flag, sequence_identifier, insertion_code)
        # hetero_flag is ' ' for standard residues, 'W' for water, 'H_x' for ligands
        if residue.id[0] != ' ':
            continue

        # 2. Extract Data
        resname = residue.resname
        # residue.id[1] is the auth_seq_id (integer)
        seq_id = residue.id[1] 
        
        try:
            one_letter = protein_letters_3to1.get(resname, 'X')
            sequence_chars.append(one_letter)
            auth_seq_ids.append(seq_id)
        except Exception:
            sequence_chars.append('X')
            auth_seq_ids.append(seq_id)
            
    return ''.join(sequence_chars), auth_seq_ids


if __name__ == "__main__":
    RCSB_ID = '1JFF'
    family_master_aln = "/Users/rtviii/dev/tubulinxyz/data/alpha_tubulin/alpha_tubulin.afasta"
    
    # Use the structure-based extraction (reads actual atoms from CIF)
    alpha_polymers = get_alpha_tubulins_from_structure(RCSB_ID)
    
    if not alpha_polymers:
        print(f"No alpha tubulins found in {RCSB_ID}")
    else:
        print(f"Found {len(alpha_polymers)} alpha tubulin(s) in {RCSB_ID}")
        mappings_dict = process_polymers_parallel(
            alpha_polymers, 
            family_master_aln,
            muscle_path="/Users/rtviii/dev/tubulinxyz/muscle3.8.1",
            max_workers=8
        )
        pprint(mappings_dict['1'])