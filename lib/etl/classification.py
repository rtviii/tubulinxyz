"""
HMM-based classification of polymer sequences.
"""

from typing import Dict, Optional
from loguru import logger

from lib.types import (
    PolymerClass,
    TubulinFamily,
    MapFamily,
    ObservedSequenceData,
)
from lib.hmm.classifier import TubulinClassifier, ClassificationResult


# Singleton classifier instance
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
) -> Optional[PolymerClass]:
    """
    Classify an observed sequence using HMM.
    
    Args:
        observed: Observed sequence data from Molstar
        rcsb_id: Structure ID for logging/caching
    
    Returns:
        The assigned PolymerClass (TubulinFamily or MapFamily), or None
    """
    classifier = get_classifier()
    
    result = classifier.classify(
        sequence=observed.sequence,
        rcsb_id=rcsb_id,
        auth_asym_id=observed.auth_asym_id,
    )
    
    return result.assigned_family


def classify_all_chains(
    sequences: Dict[str, ObservedSequenceData],
    rcsb_id: str,
) -> Dict[str, Optional[PolymerClass]]:
    """
    Classify all polymer chains in a structure.
    
    Args:
        sequences: Map of auth_asym_id -> ObservedSequenceData
        rcsb_id: Structure ID
    
    Returns:
        Map of auth_asym_id -> assigned family (or None)
    """
    results: Dict[str, Optional[PolymerClass]] = {}
    
    for auth_asym_id, observed in sequences.items():
        family = classify_sequence(observed, rcsb_id)
        results[auth_asym_id] = family
        
        status = family.value if family else "unclassified"
        marker = "+" if family else "-"
        logger.debug(f"  {auth_asym_id}: {status} {marker}")
    
    return results


def is_tubulin_family(family: Optional[PolymerClass]) -> bool:
    """Check if a family is a tubulin (not MAP)."""
    return family in (
        TubulinFamily.ALPHA,
        TubulinFamily.BETA,
        TubulinFamily.GAMMA,
        TubulinFamily.DELTA,
        TubulinFamily.EPSILON,
    )