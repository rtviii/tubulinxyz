# tubexyz/neo4j_driver/node_modification.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship

from models.types_tubulin import Modification

def node__modification(
    mod: Modification  # noqa: F821
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    mod_props = mod.model_dump()
    mod_props["modification_type"] = mod_props["modification_type"].value

    def _(tx: Transaction | ManagedTransaction):
        # Modifications are instances, so we CREATE, not MERGE.
        return tx.run("""
            CREATE (m:Modification)
            SET m = $properties
            RETURN m
            """, {"properties": mod_props}).single(strict=True)['m']
    return _

def link__polymer_to_modification(
    polymer_node: Node,
    modification_node: Node
) -> Callable[[Transaction | ManagedTransaction], Relationship]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (p:Polymer) WHERE ELEMENTID(p) = $p_elem_id
            MATCH (m:Modification) WHERE ELEMENTID(m) = $m_elem_id
            MERGE (p)-[r:HAS_MODIFICATION]->(m)
            RETURN r
            """, {
                "p_elem_id": polymer_node.element_id,
                "m_elem_id": modification_node.element_id,
            }).single(strict=True)['r']
    return _