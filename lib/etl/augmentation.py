# lib/etl/augmentation.py
"""
Post-alignment augmentation of extracted data.
"""

from typing import Dict, List
from loguru import logger
from lib.types import LigandBindingSite
from lib.etl.sequence_alignment import ChainIndexMapping


def augment_binding_sites(
    binding_sites: List[LigandBindingSite],
    chain_mappings: Dict[str, ChainIndexMapping],
) -> List[LigandBindingSite]:
    """
    Augment binding sites with master alignment indices.

    Args:
        binding_sites: Raw binding sites from Molstar extraction
        chain_mappings: Map of auth_asym_id -> ChainIndexMapping (all chains)

    Returns:
        Augmented binding sites with master_index populated where applicable
    """
    augmented_count = 0
    for site in binding_sites:
        for residue in site.residues:
            mapping = chain_mappings.get(residue.auth_asym_id)
            if mapping:
                ma_idx = mapping.auth_seq_id_to_master.get(residue.auth_seq_id)
                if ma_idx is not None:
                    residue.master_index = ma_idx
                    augmented_count += 1

    if augmented_count > 0:
        logger.debug(
            f"  Augmented {augmented_count} binding site residues with master indices"
        )

    return binding_sites
