import tempfile
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import re

@dataclass
class Annotation:
    """Represents a residue range annotation"""
    start: int  # 1-based inclusive start position
    end: int    # 1-based inclusive end position  
    label: str
    metadata: Dict = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

class TubulinAlignmentMapper:
    def __init__(self, master_profile_path: str, muscle_binary: str = "muscle"):
        self.master_profile_path = Path(master_profile_path)
        self.muscle_binary = muscle_binary
        self.master_alignment = self._load_master_alignment()
        
    def _extract_new_sequence(self, output_file: str, sequence_id: str) -> str:
        """
        Extract the newly aligned sequence from the MUSCLE output file.
        
        Args:
            output_file: Path to the MUSCLE output file
            sequence_id: The ID of the sequence we're looking for
            
        Returns:
            The aligned sequence string
        """
        from Bio import AlignIO
        
        # Read the alignment
        alignment = AlignIO.read(output_file, "fasta")
        
        # Find our sequence (should be the last one in the profile alignment)
        for record in alignment:
            if record.id == sequence_id:
                return str(record.seq)
        
        # If not found by ID, return the last sequence (our custom sequence should be last)
        if alignment:
            return str(alignment[-1].seq)
        
        raise ValueError(f"Sequence {sequence_id} not found in alignment output")
    def _load_master_alignment(self) -> List[Tuple[str, str]]:
        """Load the master MSA profile"""
        try:
            from Bio import AlignIO
            
            if not self.master_profile_path.exists():
                print(f"Warning: Master profile not found at {self.master_profile_path}")
                return []
            
            # Read the alignment
            alignment = AlignIO.read(self.master_profile_path, "fasta")
            
            # Extract sequence IDs and aligned sequences
            master_alignment = []
            for record in alignment:
                master_alignment.append((record.id, str(record.seq)))
            
            print(f"Loaded master alignment with {len(master_alignment)} sequences, length {alignment.get_alignment_length()}")
            return master_alignment
            
        except Exception as e:
            print(f"Error loading master alignment: {e}")
            return []
            
    def align_sequence(self, sequence_id: str, sequence: str) -> Tuple[str, List[int]]:
        """
        Align a new sequence to the master profile
        
        Returns:
            aligned_sequence: The sequence with gaps inserted
            mapping: List where index i contains the original position for aligned position i
                    (or -1 if it's a gap)
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as seq_file:
            seq_file.write(f">{sequence_id}\n{sequence}\n")
            seq_temp = seq_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.msa', delete=False) as output_file:
            cmd = [
                self.muscle_binary,
                '-profile',
                '-in1', str(self.master_profile_path),
                '-in2', seq_temp,
                '-out', output_file.name
            ]
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Parse the result - new sequence is the last one
            aligned_seq = self._extract_new_sequence(output_file.name, sequence_id)
            
        # Clean up temp files
        Path(seq_temp).unlink()
        Path(output_file.name).unlink()
        
        # Create mapping: for each position in aligned sequence, what was the original position?
        mapping = self._create_alignment_mapping(sequence, aligned_seq)
        
        return aligned_seq, mapping
    
    def _create_alignment_mapping(self, original_seq: str, aligned_seq: str) -> List[int]:
        mapping = []
        orig_pos = 1  # Using 1-based indexing for residues
        
        for aligned_char in aligned_seq:
            if aligned_char == '-':
                mapping.append(-1)  # This position is a gap
            else:
                mapping.append(orig_pos)
                orig_pos += 1
                
        return mapping
    
    def map_annotations(self, annotations: List[Annotation], 
                       mapping: List[int]) -> List[Annotation]:
        mapped_annotations = []
        
        for ann in annotations:
            aligned_start = self._find_aligned_position(ann.start, mapping, 'start')
            aligned_end = self._find_aligned_position(ann.end, mapping, 'end')
            
            if aligned_start is not None and aligned_end is not None:
                mapped_ann = Annotation(
                    start=aligned_start,
                    end=aligned_end,
                    label=ann.label,
                    metadata=ann.metadata.copy()
                )
                mapped_ann.metadata['original_positions'] = (ann.start, ann.end)
                mapped_annotations.append(mapped_ann)
                
        return mapped_annotations
    
    def _find_aligned_position(self, original_pos: int, mapping: List[int], 
                             boundary: str) -> Optional[int]:
        """
        Find the aligned position corresponding to an original position.
        
        Args:
            original_pos: 1-based original sequence position
            mapping: The alignment mapping from _create_alignment_mapping
            boundary: 'start' or 'end' - affects how we handle gaps at boundaries
        """
        try:
            # Find all aligned positions that map to this original position
            aligned_positions = [i for i, orig in enumerate(mapping) 
                               if orig == original_pos]
            
            if aligned_positions:
                # Should only be one, but take the first if multiple
                return aligned_positions[0] + 1  # Convert to 1-based
                
            # If no exact match, we need to handle the edge case where the 
            # original position falls in a gapped region
            if boundary == 'start':
                # For start boundary, find the first aligned position AFTER the gap
                for i, orig in enumerate(mapping):
                    if orig != -1 and orig >= original_pos:
                        return i + 1  # 1-based
            else:  # 'end'
                # For end boundary, find the last aligned position BEFORE the gap
                for i in range(len(mapping)-1, -1, -1):
                    if mapping[i] != -1 and mapping[i] <= original_pos:
                        return i + 1  # 1-based
                        
        except (IndexError, ValueError):
            pass
            
        return None

    def get_alignment_statistics(self, aligned_seq: str, mapping: List[int]) -> Dict:
        """Get statistics about the alignment"""
        total_positions = len(aligned_seq)
        gap_count = aligned_seq.count('-')
        coverage = ((total_positions - gap_count) / total_positions) * 100
        
        return {
            'total_aligned_positions': total_positions,
            'gap_count': gap_count,
            'coverage_percent': round(coverage, 2),
            'original_length': len([x for x in mapping if x != -1])
        }