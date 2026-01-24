"""
HMM-based classification of polymer sequences.
"""

from typing import Dict, List, Optional, Any
from loguru import logger

from lib.types import (
    PolymerClass,
    TubulinFamily,
    MapFamily,
    ObservedSequenceData,
    EntityClassificationResult,
    ClassificationReport,
)
from lib.hmm.classifier import TubulinClassifier, ClassificationResult


_classifier: Optional[TubulinClassifier] = None


def get_classifier() -> TubulinClassifier:
    """Get or create the singleton classifier."""
    global _classifier
    if _classifier is None:
        _classifier = TubulinClassifier(use_cache=True)
    return _classifier


def classify_sequence(
    observed: ObservedSequenceData,
    rcsb_id: str,
) -> tuple[Optional[PolymerClass], ClassificationResult]:
    """
    Classify an observed sequence using HMM.

    Returns:
        Tuple of (assigned family or None, full classification result)
    """
    classifier = get_classifier()

    result = classifier.classify(
        sequence=observed.sequence,
        rcsb_id=rcsb_id,
        auth_asym_id=observed.auth_asym_id,
    )

    return result.assigned_family, result


def is_tubulin_family(family: Optional[PolymerClass]) -> bool:
    """Check if a family is a tubulin (not MAP)."""
    return family in (
        TubulinFamily.ALPHA,
        TubulinFamily.BETA,
        TubulinFamily.GAMMA,
        TubulinFamily.DELTA,
        TubulinFamily.EPSILON,
    )


def is_map_family(family: Optional[PolymerClass]) -> bool:
    """Check if a family is a MAP."""
    return isinstance(family, MapFamily)


def build_entity_classification_result(
    entity_id: str,
    auth_asym_ids: List[str],
    sequence_length: int,
    hmm_result: ClassificationResult,
) -> EntityClassificationResult:
    """Build a structured classification result for an entity."""

    best_hit = hmm_result.best_hit

    all_hits = [
        {
            "family": hit.family.value,
            "score": hit.score,
            "evalue": hit.evalue,
            "bias": hit.bias,
        }
        for hit in sorted(hmm_result.hits, key=lambda h: h.score, reverse=True)
    ]

    return EntityClassificationResult(
        entity_id=entity_id,
        auth_asym_ids=auth_asym_ids,
        sequence_length=sequence_length,
        assigned_family=best_hit.family.value if best_hit else None,
        best_hit_score=best_hit.score if best_hit else None,
        best_hit_evalue=best_hit.evalue if best_hit else None,
        all_hits=all_hits,
    )
