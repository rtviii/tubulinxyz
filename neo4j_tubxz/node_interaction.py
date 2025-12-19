from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node
from lib.types import LigandNeighborhood

def node__ligand_interactions(
    rcsb_id: str,
    neighborhood: LigandNeighborhood
) -> Callable[[Transaction | ManagedTransaction], None]:
    
    def _(tx: Transaction | ManagedTransaction):
        # 1. Update Ligand Instance with sequence data
        tx.run("""
            MATCH (li:Instance:NonpolymerInstance {
                parent_rcsb_id: $rcsb_id, 
                auth_asym_id:   $lig_auth_id
            })
            SET li.auth_seq_id = $lig_seq_id
        """, {
            "rcsb_id": rcsb_id.upper(),
            "lig_auth_id": neighborhood.ligand_auth_asym_id,
            "lig_seq_id": neighborhood.ligand_auth_seq_id
        })

        # 2. Create Interaction Nodes
        for ix in neighborhood.interactions:
            # Participant format: [auth_asym_id, auth_seq_id, auth_comp_id, atom_id, is_ligand, master_index]
            # Get the polymer side of the interaction
            poly_part = ix.participants[1] if ix.participants[0].is_ligand else ix.participants[0]
            
            tx.run("""
                MATCH (li:Instance:NonpolymerInstance {
                    parent_rcsb_id: $rcsb_id, 
                    auth_asym_id:   $lig_auth_id
                })
                MATCH (pi:Instance:PolypeptideInstance {
                    parent_rcsb_id: $rcsb_id, 
                    auth_asym_id:   $poly_auth_id
                })
                CREATE (li)-[:HAS_INTERACTION]->(i:Interaction {
                    type:         $type,
                    atom_id:      $atom_id,
                    auth_seq_id:  $poly_seq_id,
                    auth_comp_id: $poly_comp_id,
                    master_index: $master_index
                })-[:INTERACTS_WITH]->(pi)
            """, {
                "rcsb_id":      rcsb_id.upper(),
                "lig_auth_id":  neighborhood.ligand_auth_asym_id,
                "poly_auth_id": poly_part.auth_asym_id,
                "type":         ix.type,
                "atom_id":      poly_part.atom_id,
                "poly_seq_id":  poly_part.auth_seq_id,
                "poly_comp_id": poly_part.auth_comp_id,
                "master_index": poly_part.master_index
            })

        # 3. Create proximity neighborhood links
        for res in neighborhood.neighborhood:
            tx.run("""
                MATCH (li:Instance:NonpolymerInstance {
                    parent_rcsb_id: $rcsb_id, 
                    auth_asym_id:   $lig_auth_id
                })
                MATCH (pi:Instance:PolypeptideInstance {
                    parent_rcsb_id: $rcsb_id, 
                    auth_asym_id:   $poly_auth_id
                })
                MERGE (li)-[r:NEAR_POLYMER]->(pi)
                ON CREATE SET r.residues = [$res_label]
                ON MATCH   SET r.residues = apoc.coll.toSet(r.residues + $res_label)
            """, {
                "rcsb_id":      rcsb_id.upper(),
                "lig_auth_id":  neighborhood.ligand_auth_asym_id,
                "poly_auth_id": res.auth_asym_id,
                "res_label":    f"{res.auth_comp_id}{res.auth_seq_id}"
            })
    return _