# neo4j_tubxz/node_partner_site.py
"""
MAP -> tubulin partner contact relationship creation.
Creates PARTNER_NEAR_POLYMER edges between a MAP PolypeptideInstance and a tubulin
PolypeptideInstance, with the tubulin-side contact residues stored on the edge.

This is the protein-protein analogue of node_binding_site.py (ligand NEAR_POLYMER).
It uses a DISTINCT edge label so the ligand queries -- which assume the left endpoint
is a NonpolymerInstance backed by a Chemical -- are never affected.
"""
import json
from typing import Callable, List, Dict
from neo4j import ManagedTransaction, Transaction

from lib.types import PartnerContact, BindingSiteResidue


def create_partner_site_relationships(
    tx: Transaction | ManagedTransaction,
    partner_contact: PartnerContact,
    parent_rcsb_id: str,
) -> int:
    # The contact residues are tubulin-side and may span several tubulin chains;
    # one PARTNER_NEAR_POLYMER edge is created per (MAP chain -> tubulin chain).
    residues_by_chain: Dict[str, List[BindingSiteResidue]] = {}
    for residue in partner_contact.residues:
        chain = residue.auth_asym_id
        if chain not in residues_by_chain:
            residues_by_chain[chain] = []
        residues_by_chain[chain].append(residue)

    count = 0
    for chain_id, residues in residues_by_chain.items():
        residues_json = json.dumps([r.to_dict() for r in residues])

        tx.run("""
            MATCH (mi:PolypeptideInstance {
                parent_rcsb_id: $rcsb_id,
                auth_asym_id:   $map_auth_id
            })

            MATCH (pi:PolypeptideInstance {
                parent_rcsb_id: $rcsb_id,
                auth_asym_id:   $tub_auth_id
            })

            MERGE (mi)-[r:PARTNER_NEAR_POLYMER]->(pi)
            SET r.residues_json = $residues_json,
                r.residue_count = $residue_count
        """, {
            "rcsb_id": parent_rcsb_id,
            "map_auth_id": partner_contact.partner_auth_asym_id,
            "tub_auth_id": chain_id,
            "residues_json": residues_json,
            "residue_count": len(residues),
        })
        count += 1
    return count


def process_all_partner_sites(
    rcsb_id: str,
    partner_contacts: List[PartnerContact],
) -> Callable[[Transaction | ManagedTransaction], int]:
    """
    Process all MAP->tubulin partner contacts for a structure.
    Returns a transaction function that creates all relationships.
    """
    def _(tx: Transaction | ManagedTransaction) -> int:
        total = 0
        for contact in partner_contacts:
            total += create_partner_site_relationships(tx, contact, rcsb_id)
        return total

    return _
