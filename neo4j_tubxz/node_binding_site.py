# neo4j_tubxz/node_binding_site.py
"""
Binding site relationship creation.
Creates NEAR_POLYMER edges between NonpolymerInstance and PolypeptideInstance
with residue data stored on the relationship.
"""
import json
from typing import Callable, List, Dict
from neo4j import ManagedTransaction, Transaction

from lib.types import LigandBindingSite, BindingSiteResidue


def create_binding_site_relationships(
    tx: Transaction | ManagedTransaction,
    binding_site: LigandBindingSite,
    parent_rcsb_id: str,
) -> int:
    """
    Create NEAR_POLYMER relationships for a ligand binding site.
    
    Groups residues by their parent chain (auth_asym_id) and creates
    one relationship per chain with all residues stored as a JSON array.
    
    Returns the number of relationships created.
    """
    # Group residues by chain
    residues_by_chain: Dict[str, List[BindingSiteResidue]] = {}
    for residue in binding_site.residues:
        chain = residue.auth_asym_id
        if chain not in residues_by_chain:
            residues_by_chain[chain] = []
        residues_by_chain[chain].append(residue)

    count = 0
    for chain_id, residues in residues_by_chain.items():
        residues_json = json.dumps([r.to_dict() for r in residues])

        tx.run("""
            // 1. Find the Ligand Instance. 
            // We use the Entity-link to be absolutely sure it's the right chemical
            MATCH (li:NonpolymerInstance {
                parent_rcsb_id: $rcsb_id,
                auth_asym_id:   $lig_auth_id
            })-[:INSTANCE_OF]->(e:NonpolymerEntity {chemical_id: $lig_comp_id})
            
            // Note: We relaxed auth_seq_id here because auth_asym_id is unique 
            // per structure for instances in your schema.

            // 2. Find the Protein Instance
            MATCH (pi:PolypeptideInstance {
                parent_rcsb_id: $rcsb_id,
                auth_asym_id:   $poly_auth_id
            })

            // 3. Merge the relationship
            MERGE (li)-[r:NEAR_POLYMER]->(pi)
            SET r.residues_json = $residues_json,
                r.residue_count = $residue_count
        """, {
            "rcsb_id": parent_rcsb_id,
            "lig_auth_id": binding_site.ligand_auth_asym_id,
            "lig_comp_id": binding_site.ligand_comp_id,
            "poly_auth_id": chain_id,
            "residues_json": residues_json,
            "residue_count": len(residues),
        })
        count += 1
    return count


def process_all_binding_sites(
    rcsb_id: str,
    binding_sites: List[LigandBindingSite],
) -> Callable[[Transaction | ManagedTransaction], int]:
    """
    Process all binding sites for a structure.
    Returns a transaction function that creates all relationships.
    """
    def _(tx: Transaction | ManagedTransaction) -> int:
        total = 0
        for site in binding_sites:
            total += create_binding_site_relationships(tx, site, rcsb_id)
        return total

    return _