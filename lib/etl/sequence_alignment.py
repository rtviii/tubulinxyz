# lib/etl/sequence_alignment.py
"""
Sequence alignment against family Master Alignments.

Architecture:
  1. Align CANONICAL sequence (per entity) to family MSA
     -> label_seq_id_to_master mapping + true genetic variants
  2. For each chain instance, use label_seq_id (from molstar) to build
     auth_seq_id -> master_index mapping per chain.

This means:
  - Variants (substitutions, insertions, deletions) are detected against
    the canonical sequence and represent true genetic differences, NOT
    unresolved density.
  - Per-chain mappings correctly handle the fact that different instances
    of the same entity can have different auth_seq_id numbering and
    different sets of resolved residues.
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
    EntityIndexMapping,
    ChainIndexMappingData,
    SequenceVariant,
    VariantType,
)


# ============================================================
# Result Types
# ============================================================


@dataclass
class EntityAlignmentResult:
    """Result of aligning an entity's CANONICAL sequence to a family MSA.

    Computed once per entity. Captures true genetic variants
    (canonical vs family consensus).

    canonical_to_master: label_seq_id (1-based) -> master_index (1-based) or None
    master_to_canonical: master_index (1-based) -> label_seq_id (1-based) or None
    """

    family: TubulinFamily
    entity_id: str
    canonical_sequence: str

    # label_seq_id (1-based) -> master_index (1-based) or None
    canonical_to_master: Dict[int, Optional[int]]
    # master_index (1-based) -> label_seq_id (1-based) or None
    master_to_canonical: Dict[int, Optional[int]]

    # True genetic variants (canonical vs consensus)
    variants: List[SequenceVariant] = field(default_factory=list)
    stats: Dict[str, Any] = field(default_factory=dict)

    def to_entity_index_mapping(self) -> EntityIndexMapping:
        """Convert to the pydantic model for storage/serialization."""
        return EntityIndexMapping(
            label_seq_id_to_master=self.canonical_to_master,
            master_to_label_seq_id=self.master_to_canonical,
        )

    @property
    def substitutions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.SUBSTITUTION]

    @property
    def insertions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.INSERTION]

    @property
    def deletions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.DELETION]


@dataclass
class ChainIndexMapping:
    """Per-chain mapping from auth_seq_id to master_index.

    Built cheaply from entity alignment + molstar label_seq_id data.

    auth_seq_id_to_master: author residue number -> master_index or None
    master_to_auth_seq_id: master_index -> author residue number or None
    """

    auth_asym_id: str
    entity_id: str

    # auth_seq_id -> master_index (1-based) or None
    auth_seq_id_to_master: Dict[int, Optional[int]]
    # master_index (1-based) -> auth_seq_id or None
    master_to_auth_seq_id: Dict[int, Optional[int]]

    # Canonical positions with a master mapping but not observed in this chain
    unresolved_positions: List[int] = field(default_factory=list)

    def to_chain_index_mapping_data(self) -> ChainIndexMappingData:
        """Convert to the pydantic model for storage/serialization."""
        return ChainIndexMappingData(
            auth_seq_id_to_master=self.auth_seq_id_to_master,
            master_to_auth_seq_id=self.master_to_auth_seq_id,
        )


# ============================================================
# Consensus
# ============================================================


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
                self.consensus.append(Counter(residues).most_common(1)[0][0])

    def get_residue(self, index: int) -> str:
        """Get consensus residue at 0-based index."""
        if 0 <= index < len(self.consensus):
            return self.consensus[index]
        return "?"

    @property
    def length(self) -> int:
        return len(self.consensus)


# ============================================================
# Aligner
# ============================================================


class SequenceAligner:
    """Aligns canonical sequences to family MSAs using MUSCLE."""

    def __init__(self, msa_path: Path, muscle_binary: Path):
        self.msa_path = msa_path
        self.muscle_binary = muscle_binary
        self.consensus = ConsensusCalculator(msa_path)

    def align_entity(
        self,
        canonical_sequence: str,
        entity_id: str,
        family: TubulinFamily,
        identifier: str,
    ) -> Optional[EntityAlignmentResult]:
        """
        Align an entity's canonical sequence to the family MSA.
        Returns label_seq_id-to-master mapping and true genetic variants.
        """
        if len(canonical_sequence) < 10:
            logger.warning(
                f"{identifier}: Canonical sequence too short ({len(canonical_sequence)})"
            )
            return None

        try:
            aln_target, is_original_col = self._run_muscle(
                identifier, canonical_sequence
            )
        except Exception as e:
            logger.error(f"{identifier}: MUSCLE alignment failed: {e}")
            return None

        return self._build_entity_result(
            family=family,
            entity_id=entity_id,
            canonical_sequence=canonical_sequence,
            aln_target=aln_target,
            is_original_col=is_original_col,
            identifier=identifier,
        )

    def build_chain_mapping(
        self,
        entity_result: EntityAlignmentResult,
        observed: ObservedSequenceData,
    ) -> ChainIndexMapping:
        """
        Build per-chain auth_seq_id -> master_index mapping.

        Uses label_seq_id (= canonical position) as the bridge:
            auth_seq_id  --(per residue, from molstar)--> label_seq_id
            label_seq_id = canonical_pos --(entity alignment)--> master_index
        """
        master_length = self.consensus.length

        auth_to_master: Dict[int, Optional[int]] = {}
        master_to_auth: Dict[int, Optional[int]] = {
            i: None for i in range(1, master_length + 1)
        }

        for residue in observed.residues:
            if residue.label_seq_id <= 0:
                continue
            canonical_pos = residue.label_seq_id
            master_idx = entity_result.canonical_to_master.get(canonical_pos)
            auth_to_master[residue.auth_seq_id] = master_idx
            if master_idx is not None:
                master_to_auth[master_idx] = residue.auth_seq_id

        # Unresolved: canonical positions that map to a master position
        # but have no observed residue in this chain
        observed_labels = set(
            r.label_seq_id for r in observed.residues if r.label_seq_id > 0
        )
        unresolved = sorted(
            can_pos
            for can_pos, ma_idx in entity_result.canonical_to_master.items()
            if ma_idx is not None and can_pos not in observed_labels
        )

        return ChainIndexMapping(
            auth_asym_id=observed.auth_asym_id,
            entity_id=observed.entity_id,
            auth_seq_id_to_master=auth_to_master,
            master_to_auth_seq_id=master_to_auth,
            unresolved_positions=unresolved,
        )

    # -- MUSCLE execution --

    def _run_muscle(self, seq_id: str, sequence: str) -> Tuple[str, List[bool]]:
        """
        Run MUSCLE profile alignment.

        Returns:
            (aligned_target, is_original_column_mask)

        The is_original_column mask distinguishes original MSA columns from
        new columns inserted by MUSCLE to accommodate target insertions.
        A column is "original" if ANY non-target sequence has a non-gap char.
        """
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

            n_cols = alignment.get_alignment_length()
            is_original_col = []
            for col_idx in range(n_cols):
                has_nongap = False
                for record in alignment:
                    if record.id == seq_id:
                        continue
                    if record.seq[col_idx] not in ("-", "."):
                        has_nongap = True
                        break
                is_original_col.append(has_nongap)

            return str(target_record.seq), is_original_col

        finally:
            seq_temp.unlink(missing_ok=True)
            out_temp.unlink(missing_ok=True)

    # -- Result construction --

    def _build_entity_result(
        self,
        family: TubulinFamily,
        entity_id: str,
        canonical_sequence: str,
        aln_target: str,
        is_original_col: List[bool],
        identifier: str,
    ) -> EntityAlignmentResult:
        ref_len = self.consensus.length

        canonical_to_master: Dict[int, Optional[int]] = {}
        master_to_canonical: Dict[int, Optional[int]] = {
            i + 1: None for i in range(ref_len)
        }
        variants: List[SequenceVariant] = []

        master_idx = 0  # 0-based into consensus columns
        target_idx = 0  # 0-based into canonical_sequence

        for col, t_char in enumerate(aln_target):
            if col >= len(is_original_col):
                logger.warning(
                    f"{identifier}: alignment column {col} exceeds is_original_col length "
                    f"{len(is_original_col)}. Stopping."
                )
                break

            if not is_original_col[col]:
                # Insertion column (not in original MSA)
                if t_char != "-":
                    can_pos = target_idx + 1
                    canonical_to_master[can_pos] = None
                    variants.append(
                        SequenceVariant.insertion(
                            observed_index=can_pos,
                            residue=t_char,
                            source="structural",
                        )
                    )
                    target_idx += 1
            else:
                # Original MSA column
                if t_char == "-":
                    # True genetic deletion
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
                else:
                    # Both present: match or substitution
                    can_pos = target_idx + 1

                    if master_idx < ref_len:
                        ma_pos = master_idx + 1
                        canonical_to_master[can_pos] = ma_pos
                        master_to_canonical[ma_pos] = can_pos

                        wild_type = self.consensus.get_residue(master_idx)
                        if wild_type not in ("-", "?") and t_char != wild_type:
                            variants.append(
                                SequenceVariant.substitution(
                                    master_index=ma_pos,
                                    observed_index=can_pos,
                                    wild_type=wild_type,
                                    observed=t_char,
                                    source="structural",
                                )
                            )
                    else:
                        logger.warning(
                            f"{identifier}: master_idx {master_idx} >= ref_len {ref_len} "
                            f"at column {col}, canonical_pos {can_pos}. Treating as unmapped."
                        )
                        canonical_to_master[can_pos] = None

                    target_idx += 1

                master_idx += 1

        coverage = sum(1 for v in master_to_canonical.values() if v is not None)

        stats = {
            "canonical_length": len(canonical_sequence),
            "master_alignment_length": ref_len,
            "ma_coverage": coverage,
            "ma_coverage_pct": (
                round(100 * coverage / ref_len, 1) if ref_len > 0 else 0
            ),
            "insertions": sum(1 for v in variants if v.type == VariantType.INSERTION),
            "deletions": sum(1 for v in variants if v.type == VariantType.DELETION),
            "substitutions": sum(
                1 for v in variants if v.type == VariantType.SUBSTITUTION
            ),
        }

        return EntityAlignmentResult(
            family=family,
            entity_id=entity_id,
            canonical_sequence=canonical_sequence,
            canonical_to_master=canonical_to_master,
            master_to_canonical=master_to_canonical,
            variants=variants,
            stats=stats,
        )


# ============================================================
# Factory
# ============================================================


def get_aligner_for_family(
    family: TubulinFamily, project_root: Path
) -> SequenceAligner:
    muscle_path = project_root / "muscle3.8.1"

    family_file_map = {
        TubulinFamily.ALPHA: project_root
        / "data"
        / "alpha_tubulin"
        / "alpha_tubulin.afasta",
        TubulinFamily.BETA: project_root
        / "data"
        / "beta_tubulin"
        / "beta_tubulin.afasta",
        TubulinFamily.GAMMA: project_root
        / "data"
        / "gamma_tubulin"
        / "tubulin_gamma_clean.afasta",
        TubulinFamily.DELTA: project_root
        / "data"
        / "delta_tubulin"
        / "tubulin_delta_clean.afasta",
        TubulinFamily.EPSILON: project_root
        / "data"
        / "epsilon_tubulin"
        / "tubulin_epsilon_clean.afasta",
    }

    msa_path = family_file_map.get(family)
    if not msa_path:
        raise ValueError(f"No MSA configured for family: {family}")

    return SequenceAligner(msa_path, muscle_path)
