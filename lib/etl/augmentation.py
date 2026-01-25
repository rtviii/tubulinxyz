# lib/etl/augmentation.py
"""
Post-alignment augmentation of extracted data.
"""

from typing import Dict, List, Optional
from loguru import logger

from lib.types import LigandBindingSite
from lib.etl.sequence_alignment import AlignmentResult


def augment_binding_sites(
    binding_sites: List[LigandBindingSite],
    chain_alignments: Dict[str, AlignmentResult],
) -> List[LigandBindingSite]:
    """
    Augment binding sites with master alignment indices.

    Args:
        binding_sites: Raw binding sites from Molstar extraction
        chain_alignments: Map of auth_asym_id -> AlignmentResult

    Returns:
        Augmented binding sites with master_index populated where applicable
    """
    augmented_count = 0

    for site in binding_sites:
        for residue in site.residues:
            alignment = chain_alignments.get(residue.auth_asym_id)
            if alignment:
                ma_idx = alignment.index_mapping.get_master_index(
                    residue.observed_index
                )
                if ma_idx is not None:
                    residue.master_index = ma_idx
                    augmented_count += 1

    if augmented_count > 0:
        logger.debug(
            f"  Augmented {augmented_count} binding site residues with master indices"
        )

    return binding_sites


# Stubs for future work


def extract_ptms(chain_alignments: Dict[str, AlignmentResult]) -> Dict[str, List[dict]]:
    """Stub for PTM extraction."""
    return {}


def extract_map_interfaces(
    chain_alignments: Dict[str, AlignmentResult],
    chain_families: Dict[str, Optional[str]],
) -> Dict[str, List[dict]]:
    """Stub for MAP interface extraction."""
    return {}
