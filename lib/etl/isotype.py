# lib/etl/isotype.py
"""
Isotype calling for alpha and beta tubulin entities.

Tier 1: UniProt accession lookup against known human isotype mappings.
Tier 2: Pairwise sequence alignment against human reference sequences.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Align import PairwiseAligner
from Bio import SeqIO
from loguru import logger

from lib.types import TubulinFamily


@dataclass
class IsotypeResult:
    isotype: Optional[str]
    method: Optional[str]
    confidence: Optional[float]
    details: str


MIN_ALIGNMENT_IDENTITY = 0.85

ISOTYPE_FAMILIES = {TubulinFamily.ALPHA, TubulinFamily.BETA}

_AFASTA_PATHS = {
    TubulinFamily.ALPHA: "data/alpha_tubulin/alpha_tubulin.afasta",
    TubulinFamily.BETA: "data/beta_tubulin/beta_tubulin.afasta",
}


class IsotypeCaller:

    def __init__(self, project_root: Path):
        self._project_root = project_root
        self._uniprot_map: Dict[str, str] = self._load_uniprot_map()
        self._reference_sequences: Dict[TubulinFamily, Dict[str, str]] = {}
        self._aligner = PairwiseAligner()
        self._aligner.mode = "global"
        self._aligner.match_score = 1
        self._aligner.mismatch_score = -1
        self._aligner.open_gap_score = -10
        self._aligner.extend_gap_score = -0.5

    def _load_uniprot_map(self) -> Dict[str, str]:
        path = self._project_root / "data" / "genenames" / "uniprot_to_isotype.json"
        with open(path) as f:
            return json.load(f)

    def _get_reference_sequences(self, family: TubulinFamily) -> Dict[str, str]:
        if family in self._reference_sequences:
            return self._reference_sequences[family]

        rel_path = _AFASTA_PATHS.get(family)
        if not rel_path:
            return {}

        path = self._project_root / rel_path
        refs = {}
        for record in SeqIO.parse(str(path), "fasta"):
            # Header format: "TUBA1A|Q71U36|TUBA1A_HUMAN ..."
            symbol = record.id.split("|")[0]
            raw_seq = str(record.seq).replace("-", "")
            refs[symbol] = raw_seq

        self._reference_sequences[family] = refs
        return refs

    def call_isotype(
        self,
        family: TubulinFamily,
        uniprot_accessions: List[str],
        canonical_sequence: str,
        entity_id: str,
        rcsb_id: str,
    ) -> IsotypeResult:
        if family not in ISOTYPE_FAMILIES:
            return IsotypeResult(None, None, None, f"Family {family} not supported")

        # Tier 1: UniProt lookup
        for acc in uniprot_accessions:
            if acc in self._uniprot_map:
                symbol = self._uniprot_map[acc]
                logger.debug(
                    f"  {rcsb_id} entity {entity_id}: "
                    f"UniProt {acc} -> {symbol}"
                )
                return IsotypeResult(
                    isotype=symbol,
                    method="uniprot_lookup",
                    confidence=1.0,
                    details=f"UniProt {acc} -> {symbol}",
                )

        # Tier 2: Sequence alignment fallback
        return self._align_to_references(
            family, canonical_sequence, entity_id, rcsb_id
        )

    def _align_to_references(
        self,
        family: TubulinFamily,
        sequence: str,
        entity_id: str,
        rcsb_id: str,
    ) -> IsotypeResult:
        refs = self._get_reference_sequences(family)
        if not refs:
            return IsotypeResult(None, None, None, "No reference sequences available")

        scores: List[Tuple[str, float]] = []
        for symbol, ref_seq in refs.items():
            identity = self._compute_identity(sequence, ref_seq)
            scores.append((symbol, identity))

        scores.sort(key=lambda x: x[1], reverse=True)
        best_symbol, best_identity = scores[0]
        top_str = ", ".join(f"{s}={v:.3f}" for s, v in scores[:3])

        if best_identity < MIN_ALIGNMENT_IDENTITY:
            logger.debug(
                f"  {rcsb_id} entity {entity_id}: below threshold "
                f"({best_identity:.3f} < {MIN_ALIGNMENT_IDENTITY}). Top: {top_str}"
            )
            return IsotypeResult(
                isotype=None,
                method="sequence_alignment",
                confidence=best_identity,
                details=f"Below threshold ({best_identity:.3f}). Top: {top_str}",
            )

        logger.debug(
            f"  {rcsb_id} entity {entity_id}: "
            f"aligned -> {best_symbol} ({best_identity:.3f}). Top: {top_str}"
        )
        return IsotypeResult(
            isotype=best_symbol,
            method="sequence_alignment",
            confidence=best_identity,
            details=f"Best match: {best_symbol} ({best_identity:.3f}). Top: {top_str}",
        )

    def _compute_identity(self, seq1: str, seq2: str) -> float:
        alignments = self._aligner.align(seq1, seq2)
        if not alignments:
            return 0.0

        best = alignments[0]
        # aligned is numpy array shape (2, N_blocks, 2) with [start, end] per block
        blocks_seq1 = best.aligned[0]
        blocks_seq2 = best.aligned[1]

        matches = 0
        for b1, b2 in zip(blocks_seq1, blocks_seq2):
            s1_start, s1_end = int(b1[0]), int(b1[1])
            s2_start, s2_end = int(b2[0]), int(b2[1])
            for i in range(s1_end - s1_start):
                if seq1[s1_start + i] == seq2[s2_start + i]:
                    matches += 1

        max_len = max(len(seq1), len(seq2))
        return matches / max_len if max_len > 0 else 0.0
