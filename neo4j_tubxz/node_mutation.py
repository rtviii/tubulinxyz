# tubexyz/neo4j_driver/node_mutation.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship

from lib.models.types_tubulin import Mutation

def node__mutation(
    mut: Mutation
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    mut_props = mut.model_dump()
    # Convert Enum to string value for DB
    mut_props["mutation_type"] = mut_props["mutation_type"].value

    def _(tx: Transaction | ManagedTransaction):
        # Mutations are instances, so we CREATE, not MERGE.
        # Your model includes parent_rcsb_id/auth_asym_id, so we set them.
        return tx.run("""
            CREATE (m:Mutation)
            SET m = $properties
            RETURN m
            """, {"properties": mut_props}).single(strict=True)['m']
    return _

def link__polymer_to_mutation(
    polymer_node: Node,
    mutation_node: Node
) -> Callable[[Transaction | ManagedTransaction], Relationship]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (p:Polymer) WHERE ELEMENTID(p) = $p_elem_id
            MATCH (m:Mutation) WHERE ELEMENTID(m) = $m_elem_id
            MERGE (p)-[r:HAS_MUTATION]->(m)
            RETURN r
            """, {
                "p_elem_id": polymer_node.element_id,
                "m_elem_id": mutation_node.element_id,
            }).single(strict=True)['r']
    return _