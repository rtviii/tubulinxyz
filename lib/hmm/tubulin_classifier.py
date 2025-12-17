"""
Tubulin and MAP family classifier using HMMs.
"""

import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict
from itertools import chain

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from pyhmmer.plan7 import HMM, Pipeline, TopHits
from loguru import logger

from lib.models.types_tubulin import TubulinFamily, MapFamily, HmmFamily
from lib.hmm import get_hmm_path, get_all_hmm_dirs


def get_classification_cache_dir() -> Path:
    # We can stick this in the first HMM dir or a dedicated cache dir
    # For now, let's put it in data/hmms/cache
    path = (
        Path(__file__).parent.parent.parent / "data" / "hmms" / "classification_cache"
    )
    path.mkdir(parents=True, exist_ok=True)
    return path


@dataclass
class DomainHit:
    score: float
    c_evalue: float
    i_evalue: float
    env_from: int
    env_to: int
    bias: float


@dataclass
class HMMHit:
    family: HmmFamily
    score: float
    evalue: float
    bias: float
    domains: list[DomainHit] = field(default_factory=list)


@dataclass
class ClassificationResult:
    rcsb_id: str
    auth_asym_id: str
    sequence_length: int
    hits: list[HMMHit] = field(default_factory=list)

    @property
    def best_hit(self) -> Optional[HMMHit]:
        if not self.hits:
            return None
        return max(self.hits, key=lambda h: h.score)

    @property
    def assigned_family(self) -> Optional[HmmFamily]:
        hit = self.best_hit
        return hit.family if hit else None

    def report(self) -> str:
        lines = [
            f"{'=' * 70}",
            f"Classification: {self.rcsb_id}.{self.auth_asym_id}",
            f"Length: {self.sequence_length}",
            f"{'=' * 70}",
            f"All hits:",
        ]
        sorted_hits = sorted(self.hits, key=lambda h: h.score, reverse=True)
        for i, hit in enumerate(sorted_hits):
            marker = "  <-- BEST" if i == 0 else ""
            lines.append(
                f"  {hit.family.value:<25} "
                f"score={hit.score:6.1f}  "
                f"evalue={hit.evalue:8.1e}  "
                f"{marker}"
            )

        if not self.hits:
            lines.append("  (no hits)")

        assigned = self.assigned_family
        lines.append(f"\nAssignment: {assigned.value.upper() if assigned else 'NONE'}")
        lines.append(f"{'=' * 70}")
        return "\n".join(lines)


class TubulinClassifier:
    """
    Classifies sequences against Tubulin and MAP HMMs.
    """

    def __init__(self, use_cache: bool = True):
        self.use_cache = use_cache
        self.alphabet = Alphabet.amino()
        self.hmms: Dict[HmmFamily, HMM] = {}
        self._load_hmms()

    def _load_hmms(self):
        """Load all available HMMs from disk."""

        # Iterate over all defined families in code
        all_families = list(chain(TubulinFamily, MapFamily))

        for family in all_families:
            hmm_path = get_hmm_path(family)

            if not hmm_path.exists():
                # It's common to not have HMMs for every single enum yet, so debug log is enough
                logger.debug(f"Skipping {family.value}, HMM not found at {hmm_path}")
                continue

            with open(hmm_path, "rb") as f:
                try:
                    hmm_file = pyhmmer.plan7.HMMFile(f)
                    hmm = hmm_file.read()
                    # Store mapping of Enum -> HMM Object
                    self.hmms[family] = hmm
                except Exception as e:
                    logger.error(f"Error loading {hmm_path}: {e}")

        logger.info(f"Classifier loaded {len(self.hmms)} HMM profiles")

    def classify(
        self,
        sequence: str,
        rcsb_id: str = "unknown",
        auth_asym_id: str = "X",
        force: bool = False,
    ) -> ClassificationResult:
        # Check Cache
        cache_path = get_classification_cache_dir() / f"{rcsb_id}_{auth_asym_id}.pkl"
        if self.use_cache and not force and cache_path.exists():
            try:
                with open(cache_path, "rb") as f:
                    return pickle.load(f)
            except Exception:
                pass  # Corrupt cache, re-compute

        # Prepare Sequence
        clean_seq = sequence.replace("-", "").replace(" ", "").replace("\n", "").upper()
        # Sanitize for PyHMMER (remove non-standard chars that might crash it, though Alphabet handles most)
        # For simplicity, we assume standard amino acids.

        result = ClassificationResult(
            rcsb_id=rcsb_id, auth_asym_id=auth_asym_id, sequence_length=len(clean_seq)
        )

        if len(clean_seq) < 10:
            return result  # Too short

        try:
            text_seq = TextSequence(name=b"query", sequence=clean_seq)
            digital_seq = text_seq.digitize(self.alphabet)
            seq_block = DigitalSequenceBlock(self.alphabet, [digital_seq])

            pipeline = Pipeline(self.alphabet)

            # Run search against all loaded HMMs
            for family, hmm in self.hmms.items():
                hits: TopHits = pipeline.search_hmm(hmm, seq_block)

                for hit in hits:
                    if hit.is_included():  # Only take significant hits
                        domains = [
                            DomainHit(
                                score=d.score,
                                c_evalue=d.c_evalue,
                                i_evalue=d.i_evalue,
                                env_from=d.env_from,
                                env_to=d.env_to,
                                bias=d.bias,
                            )
                            for d in hit.domains
                        ]

                        result.hits.append(
                            HMMHit(
                                family=family,
                                score=hit.score,
                                evalue=hit.evalue,
                                bias=hit.bias,
                                domains=domains,
                            )
                        )
        except Exception as e:
            logger.error(f"PyHMMER error on {rcsb_id}.{auth_asym_id}: {e}")

        # Save Cache
        if self.use_cache:
            with open(cache_path, "wb") as f:
                pickle.dump(result, f)

        return result
