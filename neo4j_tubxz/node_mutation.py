# tubexyz/neo4j_driver/node_mutation.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from lib.models.types_tubulin import Mutation

def node__mutation(
    mut: Mutation
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    mut_props = mut.model_dump()
    
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            CREATE (m:Mutation)
            SET m = $properties
            RETURN m
            """, {"properties": mut_props}).single(strict=True)['m']
    return _

def link__mutation_to_master_alignment(
    mutation_node: Node,
    master_aln_node: Node
) -> Callable[[Transaction | ManagedTransaction], Relationship]:
    """Links a mutation to the master alignment it references"""
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (m:Mutation) WHERE ELEMENTID(m) = $m_elem_id
            MATCH (a:MasterAlignment) WHERE ELEMENTID(a) = $a_elem_id
            MERGE (m)-[r:ANNOTATES_POSITION_IN]->(a)
            RETURN r
            """, {
                "m_elem_id": mutation_node.element_id,
                "a_elem_id": master_aln_node.element_id,
            }).single(strict=True)['r']
    return _