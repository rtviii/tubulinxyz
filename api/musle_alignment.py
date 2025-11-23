import tempfile
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import re

from pydantic import BaseModel


# In musle_alignment.py
class AlignmentRequest(BaseModel):
    sequence: str
    sequence_id: Optional[str] = None
    annotations: Optional[List[dict]] = []
    residue_numbers: Optional[List[int]] = None  # ✨ NEW


import tempfile
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from Bio import AlignIO


@dataclass
class Annotation:
    """Represents a residue range annotation"""

    start: int  # 1-based inclusive start position
    end: int  # 1-based inclusive end position
    label: str
    metadata: Dict = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}


class TubulinAlignmentMapper:
    def __init__(self, master_profile_path: str, muscle_binary: str = "muscle"):
        self.master_profile_path = Path(master_profile_path)
        self.muscle_binary = muscle_binary

        # Verify binary exists
        if not Path(self.muscle_binary).exists():
            print(f"WARNING: MUSCLE binary not found at {self.muscle_binary}")

    def _extract_new_sequence(self, output_file: str, sequence_id: str) -> str:
        """Extract the newly aligned sequence from the MUSCLE output file."""
        alignment = AlignIO.read(output_file, "fasta")

        # Find our sequence
        for record in alignment:
            if record.id == sequence_id:
                return str(record.seq)

        # Fallback: return the last sequence
        if alignment:
            return str(alignment[-1].seq)

        raise ValueError(f"Sequence {sequence_id} not found in alignment output")

    def align_sequence(
        self,
        sequence_id: str,
        sequence: str,
        residue_numbers: Optional[List[int]] = None,
    ) -> Tuple[str, List[int]]:
        """
        Align a new sequence to the master profile and map residue numbers.

        Args:
            sequence_id: Identifier for the temp file
            sequence: The raw amino acid string (observed sequence)
            residue_numbers: List of PDB auth_seq_ids corresponding to the sequence

        Returns:
            aligned_seq: The gapped sequence string
            mapping: List of integers. -1 for gaps, auth_seq_id for residues.
        """

        # Validation: If numbers are provided, lengths MUST match
        if residue_numbers and len(sequence) != len(residue_numbers):
            raise ValueError(
                f"Length mismatch: Sequence has {len(sequence)} residues, but {len(residue_numbers)} IDs provided."
            )

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as seq_file:
            seq_file.write(f">{sequence_id}\n{sequence}\n")
            seq_temp = seq_file.name

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".msa", delete=False
        ) as output_file:
            cmd = [
                self.muscle_binary,
                "-profile",
                "-in1",
                str(self.master_profile_path),
                "-in2",
                seq_temp,
                "-out",
                output_file.name,
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
                aligned_seq = self._extract_new_sequence(output_file.name, sequence_id)
            except subprocess.CalledProcessError as e:
                print(f"MUSCLE Error: {e.stderr.decode()}")
                raise e
            finally:
                Path(seq_temp).unlink(missing_ok=True)
                Path(output_file.name).unlink(missing_ok=True)

        # Create the mapping using the provided PDB numbers
        mapping = self._create_alignment_mapping(aligned_seq, residue_numbers)

        return aligned_seq, mapping

    def _create_alignment_mapping(
        self, aligned_seq: str, residue_numbers: Optional[List[int]]
    ) -> List[int]:
        """
        Maps aligned positions (indices) to original PDB auth_seq_ids.
        """
        mapping = []
        residue_ptr = 0  # Pointer to the index of the original sequence/numbers

        for char in aligned_seq:
            if char == "-":
                mapping.append(-1)
            else:
                # It is a residue
                if residue_numbers:
                    # Use the explicit PDB number
                    mapping.append(residue_numbers[residue_ptr])
                else:
                    # Fallback to 1-based index if no numbers provided
                    mapping.append(residue_ptr + 1)

                residue_ptr += 1

        return mapping

    def map_annotations(
        self, annotations: List[Annotation], mapping: List[int]
    ) -> List[Annotation]:
        """
        Maps annotations from original coordinate space to aligned space.
        Uses the values in `mapping` to find indices.
        """
        mapped_annotations = []

        # Create a reverse lookup: PDB_ID -> Aligned_Index
        # We store the *first* occurrence if there are duplicates (unlikely in this context)
        pdb_to_aligned = {}
        for idx, pdb_id in enumerate(mapping):
            if pdb_id != -1:
                pdb_to_aligned[pdb_id] = idx + 1  # 1-based index for Nightingale

        for ann in annotations:
            # Look up the start and end PDB IDs in our reverse map
            aligned_start = pdb_to_aligned.get(ann.start)
            aligned_end = pdb_to_aligned.get(ann.end)

            if aligned_start is not None and aligned_end is not None:
                mapped_ann = Annotation(
                    start=aligned_start,
                    end=aligned_end,
                    label=ann.label,
                    metadata=ann.metadata.copy(),
                )
                mapped_ann.metadata["original_positions"] = (ann.start, ann.end)
                mapped_annotations.append(mapped_ann)

        return mapped_annotations

    def get_alignment_statistics(self, aligned_seq: str, mapping: List[int]) -> Dict:
        total_positions = len(aligned_seq)
        gap_count = aligned_seq.count("-")
        coverage = (
            ((total_positions - gap_count) / total_positions) * 100
            if total_positions > 0
            else 0
        )

        return {
            "total_aligned_positions": total_positions,
            "gap_count": gap_count,
            "coverage_percent": round(coverage, 2),
            "original_length": len([x for x in mapping if x != -1]),
        }
