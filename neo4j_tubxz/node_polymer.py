# neo4j_tubxz/node_polymer.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node
from lib.models.types_tubulin import PolypeptideEntity, PolynucleotideEntity, Polypeptide, Polynucleotide

def node__polypeptide_entity(
    entity: PolypeptideEntity,
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates/Merges a Polypeptide Entity. 
    This holds the SEQUENCE, Taxonomy, and Tubulin Family info.
    """
    props = entity.model_dump(include={
        'entity_id', 'pdbx_description', 'formula_weight', 
        'one_letter_code', 'one_letter_code_can', 'sequence_length',
        'src_organism_names', 'src_organism_ids', 
        'host_organism_names', 'host_organism_ids',
        'uniprot_accessions'
    })
    
    # Handle Enums manually for Neo4j compatibility
    family_val = entity.family.value if entity.family else None

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MERGE (e:Entity:PolypeptideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            ON CREATE SET 
                e += $props,
                e.family = $family
            ON MATCH SET
                e += $props,
                e.family = $family
            RETURN e
        """, {
            "rcsb_id": parent_rcsb_id,
            "entity_id": entity.entity_id,
            "props": props,
            "family": family_val
        }).single(strict=True)['e']
    return _

def node__polynucleotide_entity(
    entity: PolynucleotideEntity,
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates/Merges a DNA/RNA Entity.
    """
    props = entity.model_dump(include={
        'entity_id', 'pdbx_description', 'formula_weight', 
        'one_letter_code', 'one_letter_code_can', 'sequence_length',
        'src_organism_names', 'src_organism_ids', 'polymer_type'
    })

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MERGE (e:Entity:PolynucleotideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            ON CREATE SET e += $props
            ON MATCH SET e += $props
            RETURN e
        """, {
            "rcsb_id": parent_rcsb_id,
            "entity_id": entity.entity_id,
            "props": props
        }).single(strict=True)['e']
    return _

def node__polymer_instance(
    instance: Polypeptide | Polynucleotide
) -> Callable[[Transaction | ManagedTransaction], Node]:
    specific_label = "PolypeptideInstance" if isinstance(instance, Polypeptide) else "PolynucleotideInstance"
    def _(tx: Transaction | ManagedTransaction):
        return tx.run(f"""
            MERGE (i:Instance:{specific_label} {{parent_rcsb_id: $rcsb_id, asym_id: $asym_id}})
            ON CREATE SET 
                i.auth_asym_id = $auth_id,
                i.assembly_id = $assembly_id
            
            WITH i
            MATCH (e:Entity {{parent_rcsb_id: $rcsb_id, entity_id: $entity_id}})
            MERGE (i)-[:INSTANCE_OF]->(e)
            RETURN i
        """, {
            "rcsb_id"    : instance.parent_rcsb_id,
            "asym_id"    : instance.asym_id,         
            "auth_id"    : instance.auth_asym_id,     
            "assembly_id": instance.assembly_id,
            "entity_id"  : instance.entity_id
        }).single(strict=True)['i']
    return _

def link__entity_to_structure(
    entity_node: Node, 
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Links an Entity (Polymer or Nonpolymer) to the Structure root.
    """
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (e:Entity) WHERE ELEMENTID(e) = $e_id
            MATCH (s:Structure {rcsb_id: $rcsb_id})
            MERGE (s)-[:DEFINES_ENTITY]->(e)
            RETURN s
        """, {"e_id": entity_node.element_id, "rcsb_id": parent_rcsb_id}).single()
    return _

def link__instance_to_structure(
    instance_node: Node, 
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Links a physical Instance to the Structure root.
    """
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MATCH (i:Instance) WHERE ELEMENTID(i) = $i_id
            MATCH (s:Structure {rcsb_id: $rcsb_id})
            MERGE (s)-[:HAS_INSTANCE]->(i)
            RETURN s
        """, {"i_id": instance_node.element_id, "rcsb_id": parent_rcsb_id}).single()
    return _