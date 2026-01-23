"""
Sequence alignment against family Master Alignments.

Refactored from lib/seq_aligner.py to accept observed sequences with auth_seq_ids.
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

from lib.types import TubulinFamily, ObservedSequenceData, MutationEntryData


@dataclass
class AlignmentResult:
    """Result of aligning an observed sequence to a family MSA."""
    
    family: TubulinFamily
    sequence: str
    
    # ma_to_auth[ma_idx] = auth_seq_id or -1 if gap
    ma_to_auth_map: List[int]
    
    # auth_to_ma[auth_seq_id] = ma_index (1-based) - for augmentation
    auth_to_ma: Dict[int, int]
    
    mutations: List[MutationEntryData]
    
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
        
        Args:
            observed: Observed sequence data with auth_seq_ids
            family: The tubulin family
            identifier: Identifier for logging (e.g., "6WVR_1")
        
        Returns:
            AlignmentResult or None on failure
        """
        sequence = observed.sequence
        auth_seq_ids = observed.auth_seq_ids
        
        if len(sequence) < 10:
            logger.warning(f"{identifier}: Sequence too short ({len(sequence)} residues)")
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
            aln_master=aln_master,
            aln_target=aln_target,
            identifier=identifier,
        )
    
    def _run_muscle(self, seq_id: str, sequence: str) -> Tuple[str, str]:
        """Run MUSCLE profile alignment, return (aligned_master, aligned_target)."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as seq_file:
            seq_file.write(f">{seq_id}\n{sequence}\n")
            seq_temp = Path(seq_file.name)
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msa", delete=False) as out_file:
            out_temp = Path(out_file.name)
        
        try:
            cmd = [
                str(self.muscle_binary),
                "-profile",
                "-in1", str(self.msa_path),
                "-in2", str(seq_temp),
                "-out", str(out_temp),
            ]
            result = subprocess.run(cmd, check=True, capture_output=True)
            
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
        aln_master: str,
        aln_target: str,
        identifier: str,
    ) -> AlignmentResult:
        """Build AlignmentResult from alignment strings."""
        ref_len = self.consensus.length
        
        ma_to_auth: List[int] = [-1] * ref_len
        auth_to_ma: Dict[int, int] = {}
        mutations: List[MutationEntryData] = []
        
        master_idx = 0      # 0-based index into consensus
        target_idx = 0      # index into sequence/auth_seq_ids
        insertions = 0
        
        for m_char, t_char in zip(aln_master, aln_target):
            if m_char == "-":
                # Gap in master = insertion in target
                if t_char != "-":
                    insertions += 1
                    target_idx += 1
            elif t_char == "-":
                # Gap in target = missing residue
                master_idx += 1
            else:
                # Match position
                if master_idx < ref_len and target_idx < len(auth_seq_ids):
                    auth_id = auth_seq_ids[target_idx]
                    ma_to_auth[master_idx] = auth_id
                    auth_to_ma[auth_id] = master_idx + 1  # 1-based MA index
                    
                    # Check for mutation
                    wild_type = self.consensus.get_residue(master_idx)
                    if wild_type not in ("-", "?") and t_char != wild_type:
                        mutations.append(
                            MutationEntryData(
                                ma_position=master_idx + 1,
                                wild_type=wild_type,
                                observed=t_char,
                                pdb_auth_id=auth_id,
                            )
                        )
                
                master_idx += 1
                target_idx += 1
        
        coverage = sum(1 for x in ma_to_auth if x != -1)
        
        stats = {
            "alignment_length": len(aln_target.replace("-", "")),
            "ma_coverage": coverage,
            "ma_coverage_pct": round(100 * coverage / ref_len, 1),
            "insertions": insertions,
            "total_mutations": len(mutations),
        }
        
        return AlignmentResult(
            family=family,
            sequence=sequence,
            ma_to_auth_map=ma_to_auth,
            auth_to_ma=auth_to_ma,
            mutations=mutations,
            stats=stats,
        )


def get_aligner_for_family(family: TubulinFamily, project_root: Path) -> SequenceAligner:
    """Factory to get the appropriate aligner for a tubulin family."""
    family_dir_map = {
        TubulinFamily.ALPHA: "alpha_tubulin",
        TubulinFamily.BETA: "beta_tubulin",
        TubulinFamily.GAMMA: "gamma_tubulin",
        TubulinFamily.DELTA: "delta_tubulin",
        TubulinFamily.EPSILON: "epsilon_tubulin",
    }
    
    dir_name = family_dir_map.get(family)
    if not dir_name:
        raise ValueError(f"No MSA configured for family: {family}")
    
    msa_path = project_root / "data" / dir_name / f"{dir_name}.afasta"
    muscle_path = project_root / "muscle3.8.1"
    
    if not msa_path.exists():
        raise FileNotFoundError(f"MSA not found: {msa_path}")
    
    return SequenceAligner(msa_path, muscle_path)