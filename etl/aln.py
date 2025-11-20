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

def map_annotations_between_sequences(
    annotations_indices: List[int],
    source_mappings: AlignmentMappings,
    target_mappings: AlignmentMappings
) -> List[int]:
    """
    Map annotation indices from source sequence to target sequence via master alignment.
    
    Args:
        annotations_indices: List of indices in source sequence (0-based)
        source_mappings: AlignmentMappings for source sequence
        target_mappings: AlignmentMappings for target sequence
    
    Returns:
        List of corresponding indices in target sequence (-1 if deletion, -2 if insertion)
    """
    result = []
    
    for source_idx in annotations_indices:
        # Step 1: Map source index to master index
        if source_idx >= len(source_mappings.target_to_master):
            result.append(-1)
            continue
            
        master_idx = source_mappings.target_to_master[source_idx]
        
        # If source position is an insertion (-2), can't map it
        if master_idx == -2:  # Changed from -1
            result.append(-2)
            continue
        
        # Step 2: Map master index to target index
        if master_idx >= len(target_mappings.master_to_target):
            result.append(-1)
            continue
            
        target_idx = target_mappings.master_to_target[master_idx]
        result.append(target_idx)  # Could be -1 (deletion) or valid index
    
    return result


if __name__ == "__main__":
    RCSB_ID           = '1JFF'
    cif_path          = os.path.join(TUBETL_DATA, RCSB_ID, "{}.cif".format(RCSB_ID))
    family_master_aln = "/Users/rtviii/dev/tubulinxyz/data/alpha_tubulin/alpha_tubulin.afasta"
    
    polymers = parse_mmcif_polymers(cif_path)
    pprint(polymers)

    # Pass the file path, not the parsed sequence
    mappings_dict = process_polymers_parallel(polymers, family_master_aln, max_workers=8)
    pprint(mappings_dict)


    if len(polymers) >= 2:
        entity1_id = polymers[0].entity_id
        entity2_id = polymers[1].entity_id
        
        # Some PTM sites in entity 1
        ptm_sites = [10, 25, 47]
        
        mapped_sites = map_annotations_between_sequences(
            ptm_sites,
            mappings_dict[entity1_id],
            mappings_dict[entity2_id]
        )
        
        print(f"PTM sites in {entity1_id}: {ptm_sites}")
        print(f"Mapped to {entity2_id}: {mapped_sites}")