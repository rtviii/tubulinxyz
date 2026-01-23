"""
Post-alignment augmentation of extracted data.

Adds master alignment indices to ligand neighborhoods, PTMs, etc.
"""

from typing import Dict, List, Optional
from loguru import logger

from lib.types import (
    LigandNeighborhood,
    PolymerClass,
    TubulinFamily,
)
from lib.etl.sequence_alignment import AlignmentResult


def augment_ligand_neighborhoods(
    neighborhoods: List[LigandNeighborhood],
    chain_alignments: Dict[str, AlignmentResult],
    chain_to_entity: Dict[str, str],
    entity_families: Dict[str, Optional[PolymerClass]],
) -> List[LigandNeighborhood]:
    """
    Augment ligand neighborhoods with master alignment indices.
    
    Args:
        neighborhoods: Raw ligand neighborhoods from Molstar
        chain_alignments: Map of auth_asym_id -> AlignmentResult (for tubulin chains)
        chain_to_entity: Map of auth_asym_id -> entity_id
        entity_families: Map of entity_id -> family
    
    Returns:
        Augmented neighborhoods with master_index populated where applicable
    """
    # Build unified auth_to_ma lookup across all aligned chains
    # Key: (auth_asym_id, auth_seq_id) -> ma_index
    unified_lookup: Dict[tuple, int] = {}
    
    for auth_asym_id, alignment in chain_alignments.items():
        for auth_seq_id, ma_idx in alignment.auth_to_ma.items():
            unified_lookup[(auth_asym_id, auth_seq_id)] = ma_idx
    
    augmented_count = 0
    
    for neighborhood in neighborhoods:
        # Augment interactions
        for interaction in neighborhood.interactions:
            for participant in interaction.participants:
                if not participant.is_ligand:
                    key = (participant.auth_asym_id, participant.auth_seq_id)
                    if key in unified_lookup:
                        participant.master_index = unified_lookup[key]
                        augmented_count += 1
        
        # Augment neighborhood residues
        for residue in neighborhood.neighborhood:
            key = (residue.auth_asym_id, residue.auth_seq_id)
            if key in unified_lookup:
                residue.master_index = unified_lookup[key]
                augmented_count += 1
    
    if augmented_count > 0:
        logger.debug(f"  Augmented {augmented_count} residues with MA indices")
    
    return neighborhoods


# ============================================================
# Stubs for future implementation
# ============================================================

def extract_ptms(
    chain_alignments: Dict[str, AlignmentResult],
) -> Dict[str, List[dict]]:
    """
    Extract PTMs from structure and map to master alignment.
    
    TODO: Parse pdbx_struct_mod_residue from mmCIF via Molstar.
    
    Returns:
        Map of entity_id -> list of PTM records with MA indices
    """
    # Stub - return empty
    return {}


def extract_map_interfaces(
    chain_alignments: Dict[str, AlignmentResult],
    chain_families: Dict[str, Optional[PolymerClass]],
) -> Dict[str, List[dict]]:
    """
    Extract tubulin-MAP interface residues.
    
    TODO: For each MAP chain, compute interactions with tubulin chains
    and record the tubulin-side residues with MA indices.
    
    Returns:
        Map of MAP auth_asym_id -> list of interface records
    """
    # Stub - return empty
    return {}