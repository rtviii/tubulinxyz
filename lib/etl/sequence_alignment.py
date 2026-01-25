# lib/etl/sequence_alignment.py
"""
Sequence alignment against family Master Alignments.
"""

import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

from Bio import AlignIO
from loguru import logger

from lib.types import (
    TubulinFamily,
    ObservedSequenceData,
    IndexMappingData,
    SequenceVariant,
    VariantType,
)


@dataclass
class AlignmentResult:
    """Result of aligning an observed sequence to a family MSA."""

    family: TubulinFamily
    sequence: str
    auth_asym_id: str
    entity_id: str

    index_mapping: IndexMappingData
    variants: List[SequenceVariant] = field(default_factory=list)
    stats: Dict[str, Any] = field(default_factory=dict)

    @property
    def substitutions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.SUBSTITUTION]

    @property
    def insertions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.INSERTION]

    @property
    def deletions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.DELETION]


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

        variants: List[SequenceVariant] = []

        master_idx = 0  # 0-based index into consensus
        target_idx = 0  # index into sequence/auth_seq_ids

        for m_char, t_char in zip(aln_master, aln_target):
            if m_char == "-":
                # Gap in master = insertion in target
                if t_char != "-" and target_idx < len(auth_seq_ids):
                    auth_id = auth_seq_ids[target_idx]
                    observed_to_master[auth_id] = None  # No MA position
                    variants.append(
                        SequenceVariant.insertion(
                            observed_index=auth_id,
                            residue=t_char,
                            source="structural",
                        )
                    )
                    target_idx += 1

            elif t_char == "-":
                # Gap in target = deletion/missing residue
                ma_pos = master_idx + 1
                expected = self.consensus.get_residue(master_idx)
                if expected not in ("-", "?"):
                    variants.append(
                        SequenceVariant.deletion(
                            master_index=ma_pos,
                            expected=expected,
                            source="structural",
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

                    # Check for substitution
                    wild_type = self.consensus.get_residue(master_idx)
                    if wild_type not in ("-", "?") and t_char != wild_type:
                        variants.append(
                            SequenceVariant.substitution(
                                master_index=ma_pos,
                                observed_index=auth_id,
                                wild_type=wild_type,
                                observed=t_char,
                                source="structural",
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
        substitution_count = sum(1 for v in variants if v.type == VariantType.SUBSTITUTION)
        insertion_count = sum(1 for v in variants if v.type == VariantType.INSERTION)
        deletion_count = sum(1 for v in variants if v.type == VariantType.DELETION)

        stats = {
            "alignment_length": len(sequence),
            "master_alignment_length": ref_len,
            "ma_coverage": coverage,
            "ma_coverage_pct": round(100 * coverage / ref_len, 1) if ref_len > 0 else 0,
            "insertions": insertion_count,
            "deletions": deletion_count,
            "substitutions": substitution_count,
        }

        return AlignmentResult(
            family=family,
            sequence=sequence,
            auth_asym_id=auth_asym_id,
            entity_id=entity_id,
            index_mapping=index_mapping,
            variants=variants,
            stats=stats,
        )


def get_aligner_for_family(family: TubulinFamily, project_root: Path) -> SequenceAligner:
    muscle_path = project_root / "muscle3.8.1"

    family_file_map = {
        TubulinFamily.ALPHA: "/Users/rtviii/dev/tubulinxyz/data/alpha_tubulin/alpha_tubulin.afasta",
        TubulinFamily.BETA: "/Users/rtviii/dev/tubulinxyz/data/beta_tubulin/beta_tubulin.afasta",
        TubulinFamily.GAMMA: "/Users/rtviii/dev/tubulinxyz/data/gamma_tubulin/tubulin_gamma_clean.afasta",
        TubulinFamily.DELTA: "/Users/rtviii/dev/tubulinxyz/data/delta_tubulin/tubulin_delta_clean.afasta",
        TubulinFamily.EPSILON: "/Users/rtviii/dev/tubulinxyz/data/epsilon_tubulin/tubulin_epsilon_clean.afasta",
    }

    msa_path = family_file_map.get(family)
    if not msa_path:
        raise ValueError(f"No MSA configured for family: {family}")

    return SequenceAligner(Path(msa_path), muscle_path)