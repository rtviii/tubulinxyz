"""
Post-alignment augmentation of extracted data.
"""

from typing import Dict, List, Optional
from loguru import logger

from lib.types import (
    SimplifiedLigandNeighborhood,
    IndexMappingData,
)
from lib.etl.sequence_alignment import AlignmentResult

from typing import Dict, List
from loguru import logger

from lib.types import LigandBindingSite
from lib.etl.sequence_alignment import AlignmentResult


def augment_ligand_neighborhoods(
    neighborhoods: List[SimplifiedLigandNeighborhood],
    chain_alignments: Dict[str, AlignmentResult],
) -> List[SimplifiedLigandNeighborhood]:
    """
    Augment ligand neighborhoods with master alignment indices.

    Args:
        neighborhoods: Raw ligand neighborhoods from Molstar
        chain_alignments: Map of auth_asym_id -> AlignmentResult

    Returns:
        Augmented neighborhoods with master_index populated where applicable
    """
    augmented_count = 0

    for neighborhood in neighborhoods:
        for residue in neighborhood.neighborhood_residues:
            alignment = chain_alignments.get(residue.auth_asym_id)
            if alignment:
                ma_idx = alignment.index_mapping.get_master_index(
                    residue.observed_index
                )
                if ma_idx is not None:
                    residue.master_index = ma_idx
                    augmented_count += 1

    if augmented_count > 0:
        logger.debug(f"  Augmented {augmented_count} residues with master indices")

    return neighborhoods


# Stubs for future


def extract_ptms(chain_alignments: Dict[str, AlignmentResult]) -> Dict[str, List[dict]]:
    """Stub for PTM extraction."""
    return {}


def extract_map_interfaces(
    chain_alignments: Dict[str, AlignmentResult],
    chain_families: Dict[str, Optional[str]],
) -> Dict[str, List[dict]]:
    """Stub for MAP interface extraction."""
    return {}



def augment_binding_sites(
    binding_sites: List[LigandBindingSite],
    chain_alignments: Dict[str, AlignmentResult],
) -> List[LigandBindingSite]:
    """
    Augment binding sites with master alignment indices.
    """
    augmented_count = 0
    
    for site in binding_sites:
        for residue in site.residues:
            alignment = chain_alignments.get(residue.auth_asym_id)
            if alignment:
                ma_idx = alignment.index_mapping.get_master_index(residue.observed_index)
                if ma_idx is not None:
                    residue.master_index = ma_idx
                    augmented_count += 1
    
    if augmented_count > 0:
        logger.debug(f"  Augmented {augmented_count} binding site residues with master indices")
    
    return binding_sites