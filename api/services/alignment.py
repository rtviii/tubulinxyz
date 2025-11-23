import tempfile
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field
from Bio import AlignIO

@dataclass
class Annotation:
    """Represents a residue range annotation internally"""
    start: int
    end: int
    label: str
    metadata: Dict = field(default_factory=dict)

class TubulinAlignmentMapper:
    def __init__(self, master_profile_path: str, muscle_binary: str):
        self.master_profile_path = Path(master_profile_path)
        self.muscle_binary = muscle_binary
        
        if not Path(self.muscle_binary).exists():
            print(f"WARNING: MUSCLE binary not found at {self.muscle_binary}")

    def align_sequence(self, sequence_id: str, sequence: str, residue_numbers: Optional[List[int]] = None) -> Tuple[str, List[int]]:
        if residue_numbers and len(sequence) != len(residue_numbers):
            raise ValueError(f"Length mismatch: Sequence len {len(sequence)} vs IDs len {len(residue_numbers)}.")

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as seq_file:
            seq_file.write(f">{sequence_id}\n{sequence}\n")
            seq_temp = seq_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.msa', delete=False) as output_file:
            # Clean up temp files in finally block or use context managers properly in production
            try:
                cmd = [
                    self.muscle_binary, '-profile',
                    '-in1', str(self.master_profile_path),
                    '-in2', seq_temp,
                    '-out', output_file.name
                ]
                subprocess.run(cmd, check=True, capture_output=True)
                aligned_seq = self._extract_new_sequence(output_file.name, sequence_id)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"MUSCLE Error: {e.stderr.decode()}")
            finally:
                Path(seq_temp).unlink(missing_ok=True)
                # output_file is cleaned up by OS or should be unlinked here after read

        mapping = self._create_alignment_mapping(aligned_seq, residue_numbers)
        return aligned_seq, mapping

    def _extract_new_sequence(self, output_file: str, sequence_id: str) -> str:
        alignment = AlignIO.read(output_file, "fasta")
        for record in alignment:
            if record.id == sequence_id:
                return str(record.seq)
        # Fallback if ID got truncated or changed
        if alignment:
            return str(alignment[-1].seq)
        raise ValueError(f"Sequence {sequence_id} not found in alignment output")

    def _create_alignment_mapping(self, aligned_seq: str, residue_numbers: Optional[List[int]]) -> List[int]:
        mapping = []
        residue_ptr = 0 
        for char in aligned_seq:
            if char == '-':
                mapping.append(-1)
            else:
                if residue_numbers:
                    mapping.append(residue_numbers[residue_ptr])
                else:
                    mapping.append(residue_ptr + 1)
                residue_ptr += 1
        return mapping

    def map_annotations(self, annotations: List[Annotation], mapping: List[int]) -> List[Annotation]:
        mapped_annotations = []
        pdb_to_aligned = {}
        
        # Reverse lookup: PDB_ID -> Aligned_Index (1-based for frontend visualizers)
        for idx, pdb_id in enumerate(mapping):
            if pdb_id != -1:
                pdb_to_aligned[pdb_id] = idx + 1 

        for ann in annotations:
            aligned_start = pdb_to_aligned.get(ann.start)
            aligned_end = pdb_to_aligned.get(ann.end)
            
            if aligned_start is not None and aligned_end is not None:
                new_meta = ann.metadata.copy()
                new_meta['original_positions'] = (ann.start, ann.end)
                mapped_annotations.append(Annotation(
                    start=aligned_start, end=aligned_end, label=ann.label, metadata=new_meta
                ))
        return mapped_annotations

    def get_alignment_statistics(self, aligned_seq: str, mapping: List[int]) -> Dict:
        total = len(aligned_seq)
        gaps = aligned_seq.count('-')
        return {
            'total_aligned_positions': total,
            'gap_count': gaps,
            'coverage_percent': round(((total - gaps) / total) * 100, 2) if total > 0 else 0,
            'original_length': len([x for x in mapping if x != -1])
        }