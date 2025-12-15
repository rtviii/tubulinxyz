"""
Tubulin family classifier using HMMs.
"""

import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from pyhmmer.plan7 import HMM, Pipeline, TopHits
from loguru import logger

from lib.models.types_tubulin import TubulinFamily
from lib.hmm import get_hmm_dir, get_hmm_path


def get_classification_cache_dir() -> Path:

    return get_hmm_dir() / "classification_cache"



@dataclass
class DomainHit:
    """Single domain within an HMM hit."""
    score: float
    c_evalue: float
    i_evalue: float
    env_from: int    # envelope start (sequence coords)
    env_to: int      # envelope end
    bias: float

@dataclass
class HMMHit:
    """Single HMM hit result with domain details."""

    family: TubulinFamily
    score: float
    evalue: float
    bias: float
    domains: list[DomainHit] = field(default_factory=list)

    def __repr__(self):
        return f"HMMHit({self.family.value}: score={self.score:.1f}, evalue={self.evalue:.2e}, domains={len(self.domains)})"


@dataclass
class ClassificationResult:
    """Full classification result for a sequence."""

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
    def assigned_family(self) -> Optional[TubulinFamily]:
        hit = self.best_hit
        return hit.family if hit else None

    def report(self) -> str:
        """Generate a detailed report string."""
        lines = [
            f"{'=' * 70}",
            f"Classification: {self.rcsb_id}.{self.auth_asym_id}",
            f"Sequence length: {self.sequence_length}",
            f"{'=' * 70}",
            f"",
            f"All hits (sorted by bitscore):",
        ]
        sorted_hits = sorted(self.hits, key=lambda h: h.score, reverse=True)
        for i, hit in enumerate(sorted_hits):
            marker = "  <-- BEST" if i == 0 else ""
            lines.append(
                f"  {hit.family.value:10s}  "
                f"score={hit.score:8.2f}  "
                f"evalue={hit.evalue:12.2e}  "
                f"bias={hit.bias:6.2f}{marker}"
            )
            # Show domains
            for j, dom in enumerate(hit.domains):
                coverage = (dom.env_to - dom.env_from + 1) / self.sequence_length * 100
                lines.append(
                    f"      domain {j + 1}: score={dom.score:6.1f}  "
                    f"env[{dom.env_from}-{dom.env_to}] ({coverage:.1f}%)  "
                    f"evalue={dom.c_evalue:.2e}"
                )
        if not self.hits:
            lines.append("  (no hits above threshold)")
        lines.append("")
        assigned = self.assigned_family
        lines.append(f"Assignment: {assigned.value.upper() if assigned else 'NONE'}")
        lines.append(f"{'=' * 70}")
        return "\n".join(lines)

    def domain_report(self) -> str:
        """Detailed domain-level comparison between top hits."""
        if len(self.hits) < 2:
            return "Not enough hits for comparison"
        
        sorted_hits = sorted(self.hits, key=lambda h: h.score, reverse=True)
        top = sorted_hits[0]
        second = sorted_hits[1]
        
        lines = [
            f"Domain comparison: {top.family.value} vs {second.family.value}",
            f"Score difference: {top.score - second.score:.1f}",
            "",
            f"TOP ({top.family.value}, score={top.score:.1f}):",
        ]
        
        for i, dom in enumerate(top.domains):
            coverage = (dom.env_to - dom.env_from + 1) / self.sequence_length * 100
            lines.append(
                f"  dom{i+1}: score={dom.score:6.1f}  "
                f"env[{dom.env_from:4d}-{dom.env_to:4d}] ({coverage:4.1f}% of seq)  "
                f"evalue={dom.c_evalue:.2e}"
            )
        
        lines.append(f"\nSECOND ({second.family.value}, score={second.score:.1f}):")
        for i, dom in enumerate(second.domains):
            coverage = (dom.env_to - dom.env_from + 1) / self.sequence_length * 100
            lines.append(
                f"  dom{i+1}: score={dom.score:6.1f}  "
                f"env[{dom.env_from:4d}-{dom.env_to:4d}] ({coverage:4.1f}% of seq)  "
                f"evalue={dom.c_evalue:.2e}"
            )
        
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Serialize to dict for caching/JSON."""
        return {
            "rcsb_id": self.rcsb_id,
            "auth_asym_id": self.auth_asym_id,
            "sequence_length": self.sequence_length,
            "assigned_family": self.assigned_family.value if self.assigned_family else None,
            "hits": [
                {
                    "family": h.family.value,
                    "score": h.score,
                    "evalue": h.evalue,
                    "bias": h.bias,
                    "domains": [
                        {
                            "score": d.score,
                            "c_evalue": d.c_evalue,
                            "env_range": [d.env_from, d.env_to],
                        }
                        for d in h.domains
                    ]
                }
                for h in self.hits
            ]
        }


class TubulinClassifier:
    """
    Classifier that loads tubulin family HMMs and classifies sequences.
    """

    def __init__(
        self,
        hmm_dir: str | Path | None = None,
        bitscore_threshold: float = 50.0,
        use_cache: bool = True,
    ):
        self.hmm_dir = Path(hmm_dir) if hmm_dir else get_hmm_dir()
        self.bitscore_threshold = bitscore_threshold
        self.use_cache = use_cache
        self.alphabet = Alphabet.amino()
        self.hmms: dict[TubulinFamily, HMM] = {}

        self._load_hmms()

        if use_cache:
            get_classification_cache_dir().mkdir(exist_ok=True)

    def _load_hmms(self):
        """Load all family HMMs from disk."""
        for family in TubulinFamily:
            hmm_path = get_hmm_path(family.value)

            if not hmm_path.exists():
                logger.warning(f"HMM not found: {hmm_path}")
                continue

            with open(hmm_path, "rb") as f:
                hmm_file = pyhmmer.plan7.HMMFile(f)
                hmm = hmm_file.read()
                self.hmms[family] = hmm
                logger.debug(f"Loaded {family.value} HMM: M={hmm.M}")

        logger.info(f"Loaded {len(self.hmms)} HMMs")

    def _cache_key(self, rcsb_id: str, auth_asym_id: str) -> Path:
        return get_classification_cache_dir() / f"{rcsb_id}_{auth_asym_id}.pkl"

    def _load_cached(
        self, rcsb_id: str, auth_asym_id: str
    ) -> ClassificationResult | None:
        cache_path = self._cache_key(rcsb_id, auth_asym_id)
        if cache_path.exists():
            with open(cache_path, "rb") as f:
                logger.debug(f"Loaded cached result for {rcsb_id}.{auth_asym_id}")
                return pickle.load(f)
        return None

    def _save_cache(self, result: ClassificationResult):
        cache_path = self._cache_key(result.rcsb_id, result.auth_asym_id)
        with open(cache_path, "wb") as f:
            pickle.dump(result, f)
        logger.debug(f"Cached result for {result.rcsb_id}.{result.auth_asym_id}")

    def classify(
        self,
        sequence: str,
        rcsb_id: str = "unknown",
        auth_asym_id: str = "X",
        force: bool = False,
    ) -> ClassificationResult:
        """
        Classify a sequence against all tubulin family HMMs.
        """
        if self.use_cache and not force:
            cached = self._load_cached(rcsb_id, auth_asym_id)
            if cached:
                return cached

        clean_seq = sequence.replace("-", "").replace(" ", "").replace("\n", "").upper()

        result = ClassificationResult(
            rcsb_id=rcsb_id, auth_asym_id=auth_asym_id, sequence_length=len(clean_seq)
        )

        seq_name = f"{rcsb_id}_{auth_asym_id}".encode()
        text_seq = TextSequence(name=seq_name, sequence=clean_seq)
        digital_seq = text_seq.digitize(self.alphabet)
        seq_block = DigitalSequenceBlock(self.alphabet, [digital_seq])

        pipeline = Pipeline(self.alphabet)

        for family, hmm in self.hmms.items():
            hits: TopHits = pipeline.search_hmm(hmm, seq_block)
            
            for hit in hits:
                # Extract domain info
                domains = []
                for dom in hit.domains:
                    domains.append(DomainHit(
                        score=dom.score,
                        c_evalue=dom.c_evalue,
                        i_evalue=dom.i_evalue,
                        env_from=dom.env_from,
                        env_to=dom.env_to,
                        bias=dom.bias,
                    ))
                
                result.hits.append(HMMHit(
                    family=family,
                    score=hit.score,
                    evalue=hit.evalue,
                    bias=hit.bias,
                    domains=domains
                ))

        if result.best_hit:
            logger.info(
                f"{rcsb_id}.{auth_asym_id} -> {result.assigned_family.value} "
                f"(score={result.best_hit.score:.1f})"
            )
        else:
            logger.warning(f"{rcsb_id}.{auth_asym_id} -> no hits")

        if self.use_cache:
            self._save_cache(result)

        return result

    def classify_verbose(
        self,
        sequence: str,
        rcsb_id: str = "unknown",
        auth_asym_id: str = "X",
        force: bool = False,
    ) -> ClassificationResult:
        """Classify and print detailed report."""
        result = self.classify(sequence, rcsb_id, auth_asym_id, force)
        print(result.report())
        return result

    def hmm_info(self) -> dict:
        """Get info about loaded HMMs."""
        return {
            family.value: {
                "name": hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name,
                "M": hmm.M,
                "nseq": hmm.nseq,
            }
            for family, hmm in self.hmms.items()
        }
