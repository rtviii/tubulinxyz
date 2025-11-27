# tubexyz/neo4j_driver/node_master_alignment.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship

from lib.models.types_tubulin import AlignmentMapping, MasterAlignment, TubulinFamily

def node__master_alignment(
    aln: MasterAlignment
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    aln_props = aln.model_dump()
    aln_props["family"] = aln_props["family"].value # Convert Enum

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MERGE (a:MasterAlignment { family: $family, version: $version })
            ON CREATE SET
                a.fasta_content = $fasta_content,
                a.created_date = $created_date,
                a.description = $description
            ON MATCH SET
                a.fasta_content = $fasta_content,
                a.created_date = $created_date,
                a.description = $description
            RETURN a
            """, aln_props).single(strict=True)['a']
    return _

def get_master_alignment(
    family: TubulinFamily, version: str
) -> Callable[[Transaction | ManagedTransaction], Node | None]:
    """Helper to fetch a specific master alignment node."""
    def _(tx: Transaction | ManagedTransaction):
        result = tx.run("""
            MATCH (a:MasterAlignment { family: $family, version: $version })
            RETURN a
            """, {"family": family.value, "version": version})
        
        record = result.single()
        return record['a'] if record else None
    return _

def link__polymer_to_master_alignment(
    polymer_node: Node,
    master_aln_node: Node,
    mapping: AlignmentMapping
) -> Callable[[Transaction | ManagedTransaction], Relationship]:
    
    mapping_props = mapping.model_dump()

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (p:Polymer) WHERE ELEMENTID(p) = $p_elem_id
            MATCH (a:MasterAlignment) WHERE ELEMENTID(a) = $a_elem_id
            MERGE (p)-[r:IS_MAPPED_TO]->(a)
            SET
                r.seqres_to_master = $seqres_to_master,
                r.master_to_seqres = $master_to_seqres
            RETURN r
            """, {
                "p_elem_id": polymer_node.element_id,
                "a_elem_id": master_aln_node.element_id,
                **mapping_props
            }).single(strict=True)['r']
    return _