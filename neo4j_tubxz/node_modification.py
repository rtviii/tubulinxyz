# tubexyz/neo4j_driver/node_modification.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from lib.models.types_tubulin import Modification

def node__modification(
    mod: Modification
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    mod_props = mod.model_dump()
    
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            CREATE (m:Modification)
            SET m = $properties
            RETURN m
            """, {"properties": mod_props}).single(strict=True)['m']
    return _

def link__modification_to_master_alignment(
    modification_node: Node,
    master_aln_node: Node
) -> Callable[[Transaction | ManagedTransaction], Relationship]:
    """Links a modification to the master alignment it references"""
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (m:Modification) WHERE ELEMENTID(m) = $m_elem_id
            MATCH (a:MasterAlignment) WHERE ELEMENTID(a) = $a_elem_id
            MERGE (m)-[r:ANNOTATES_POSITION_IN]->(a)
            RETURN r
            """, {
                "m_elem_id": modification_node.element_id,
                "a_elem_id": master_aln_node.element_id,
            }).single(strict=True)['r']
    return _