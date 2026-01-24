"""
Sequence alignment against family Master Alignments.
"""

import json
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

from Bio import AlignIO
from loguru import logger

from lib.types import (
    TubulinFamily,
    ObservedSequenceData,
    IndexMappingData,
    MutationRecord,
    InsertionRecord,
    DeletionRecord,
    MutationsAndIndels,
)


@dataclass
class AlignmentResult:
    """Result of aligning an observed sequence to a family MSA."""

    family: TubulinFamily
    sequence: str
    auth_asym_id: str
    entity_id: str

    # Structured index mapping
    index_mapping: IndexMappingData

    # Mutations and indels
    mutations: List[MutationRecord]
    insertions: List[InsertionRecord]
    deletions: List[DeletionRecord]

    # Stats
    stats: Dict[str, Any]


class ConsensusCalculator:
    """Calculates consensus sequence from a family MSA."""

    def __init__(self, msa_path: Path):
        self.path = msa_path
        self.consensus: List[str] = []
        self._calculate()

    def _calculate(self):
        if not self.path.exists():
            raise FileNotFoundError(f"MSA not found: {self.path}")

        alignment = AlignIO.read(str(self.path), "fasta")
        length = alignment.get_alignment_length()

        for i in range(length):
            column = alignment[:, i]
            residues = [r for r in column if r not in ("-", ".")]

            if not residues:
                self.consensus.append("-")
            else:
                most_common = Counter(residues).most_common(1)[0][0]
                self.consensus.append(most_common)

    def get_residue(self, index: int) -> str:
        """Get consensus residue at 0-based index."""
        if 0 <= index < len(self.consensus):
            return self.consensus[index]
        return "?"

    @property
    def length(self) -> int:
        return len(self.consensus)


class SequenceAligner:
    """Aligns observed sequences to family MSAs using MUSCLE."""

    def __init__(self, msa_path: Path, muscle_binary: Path):
        self.msa_path = msa_path
        self.muscle_binary = muscle_binary
        self.consensus = ConsensusCalculator(msa_path)

    def align(
        self,
        observed: ObservedSequenceData,
        family: TubulinFamily,
        identifier: str,
    ) -> Optional[AlignmentResult]:
        """
        Align an observed sequence to the family MSA.
        """
        sequence = observed.sequence
        auth_seq_ids = observed.auth_seq_ids

        if len(sequence) < 10:
            logger.warning(
                f"{identifier}: Sequence too short ({len(sequence)} residues)"
            )
            return None

        try:
            aln_master, aln_target = self._run_muscle(identifier, sequence)
        except Exception as e:
            logger.error(f"{identifier}: MUSCLE alignment failed: {e}")
            return None

        return self._build_result(
            family=family,
            sequence=sequence,
            auth_seq_ids=auth_seq_ids,
            auth_asym_id=observed.auth_asym_id,
            entity_id=observed.entity_id,
            aln_master=aln_master,
            aln_target=aln_target,
            identifier=identifier,
        )

    def _run_muscle(self, seq_id: str, sequence: str) -> Tuple[str, str]:
        """Run MUSCLE profile alignment."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as seq_file:
            seq_file.write(f">{seq_id}\n{sequence}\n")
            seq_temp = Path(seq_file.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".msa", delete=False
        ) as out_file:
            out_temp = Path(out_file.name)

        try:
            cmd = [
                str(self.muscle_binary),
                "-profile",
                "-in1",
                str(self.msa_path),
                "-in2",
                str(seq_temp),
                "-out",
                str(out_temp),
            ]
            subprocess.run(cmd, check=True, capture_output=True)

            alignment = AlignIO.read(str(out_temp), "fasta")

            target_record = next((r for r in alignment if r.id == seq_id), None)
            if not target_record:
                raise ValueError(f"Sequence {seq_id} not found in alignment output")

            return str(alignment[0].seq), str(target_record.seq)

        finally:
            seq_temp.unlink(missing_ok=True)
            out_temp.unlink(missing_ok=True)

    def _build_result(
        self,
        family: TubulinFamily,
        sequence: str,
        auth_seq_ids: List[int],
        auth_asym_id: str,
        entity_id: str,
        aln_master: str,
        aln_target: str,
        identifier: str,
    ) -> AlignmentResult:
        """Build AlignmentResult from alignment strings."""
        ref_len = self.consensus.length

        # Build mappings
        observed_to_master: Dict[int, Optional[int]] = {}
        master_to_observed: Dict[int, Optional[int]] = {
            i + 1: None for i in range(ref_len)
        }

        mutations: List[MutationRecord] = []
        insertions: List[InsertionRecord] = []
        deletions: List[DeletionRecord] = []

        master_idx = 0  # 0-based index into consensus
        target_idx = 0  # index into sequence/auth_seq_ids

        for m_char, t_char in zip(aln_master, aln_target):
            if m_char == "-":
                # Gap in master = insertion in target
                if t_char != "-" and target_idx < len(auth_seq_ids):
                    auth_id = auth_seq_ids[target_idx]
                    observed_to_master[auth_id] = None  # No MA position
                    insertions.append(
                        InsertionRecord(
                            observed_index=auth_id,
                            residue=t_char,
                        )
                    )
                    target_idx += 1

            elif t_char == "-":
                # Gap in target = deletion/missing residue
                ma_pos = master_idx + 1
                expected = self.consensus.get_residue(master_idx)
                if expected not in ("-", "?"):
                    deletions.append(
                        DeletionRecord(
                            master_index=ma_pos,
                            expected_residue=expected,
                        )
                    )
                master_idx += 1

            else:
                # Match position
                if master_idx < ref_len and target_idx < len(auth_seq_ids):
                    auth_id = auth_seq_ids[target_idx]
                    ma_pos = master_idx + 1  # 1-based

                    observed_to_master[auth_id] = ma_pos
                    master_to_observed[ma_pos] = auth_id

                    # Check for mutation
                    wild_type = self.consensus.get_residue(master_idx)
                    if wild_type not in ("-", "?") and t_char != wild_type:
                        mutations.append(
                            MutationRecord(
                                master_index=ma_pos,
                                observed_index=auth_id,
                                wild_type=wild_type,
                                observed=t_char,
                            )
                        )

                master_idx += 1
                target_idx += 1

        # Build index mapping
        index_mapping = IndexMappingData(
            observed_to_master=observed_to_master,
            master_to_observed=master_to_observed,
        )

        coverage = sum(1 for v in master_to_observed.values() if v is not None)

        stats = {
            "alignment_length": len(sequence),
            "master_alignment_length": ref_len,
            "ma_coverage": coverage,
            "ma_coverage_pct": round(100 * coverage / ref_len, 1) if ref_len > 0 else 0,
            "insertions": len(insertions),
            "deletions": len(deletions),
            "mutations": len(mutations),
        }

        return AlignmentResult(
            family=family,
            sequence=sequence,
            auth_asym_id=auth_asym_id,
            entity_id=entity_id,
            index_mapping=index_mapping,
            mutations=mutations,
            insertions=insertions,
            deletions=deletions,
            stats=stats,
        )


def get_aligner_for_family(family: TubulinFamily, project_root: Path) -> SequenceAligner:
    family_file_map = {
        TubulinFamily.ALPHA: "tubulin_alpha",
        TubulinFamily.BETA: "tubulin_beta",
        TubulinFamily.GAMMA: "tubulin_gamma",
        TubulinFamily.DELTA: "tubulin_delta",
        TubulinFamily.EPSILON: "tubulin_epsilon",
    }
    
    file_name = family_file_map.get(family)
    if not file_name:
        raise ValueError(f"No MSA configured for family: {family}")
    
    # This path should match your actual files:
    msa_path = project_root / "data" / "sequences" / "tubulin" / f"{file_name}.afasta"
    muscle_path = project_root / "muscle3.8.1"
    
    if not msa_path.exists():
        raise FileNotFoundError(f"MSA not found: {msa_path}")
    
    return SequenceAligner(msa_path, muscle_path)