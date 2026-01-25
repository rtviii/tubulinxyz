# neo4j_tubxz/node_polymer.py
"""
Polymer entity and instance node creation.
"""
import json
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node

from lib.types import (
    PolypeptideEntity,
    PolynucleotideEntity,
    Polypeptide,
    Polynucleotide,
)


def node__polypeptide_entity(
    entity: PolypeptideEntity,
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates/Merges a Polypeptide Entity.
    Stores sequence, taxonomy, family, index mappings, and alignment stats.
    """
    # Base properties
    props = {
        "entity_id": entity.entity_id,
        "parent_rcsb_id": parent_rcsb_id,
        "pdbx_description": entity.pdbx_description,
        "pdbx_strand_ids": entity.pdbx_strand_ids,
        "one_letter_code": entity.one_letter_code,
        "one_letter_code_can": entity.one_letter_code_can,
        "sequence_length": entity.sequence_length,
        "src_organism_names": entity.src_organism_names,
        "src_organism_ids": entity.src_organism_ids,
        "host_organism_names": entity.host_organism_names,
        "host_organism_ids": entity.host_organism_ids,
        "uniprot_accessions": entity.uniprot_accessions,
    }

    # Family (handle enum)
    family_val = entity.family.value if entity.family else None

    # Index mappings as JSON strings
    if entity.index_mapping:
        mapping_json = entity.index_mapping.to_json_dict()
        props["observed_to_master_json"] = mapping_json["observed_to_master_json"]
        props["master_to_observed_json"] = mapping_json["master_to_observed_json"]

    # Alignment stats as JSON
    if entity.alignment_stats:
        props["alignment_stats_json"] = json.dumps(entity.alignment_stats)

    # Clean None values
    props = {k: v for k, v in props.items() if v is not None}

    def _(tx: Transaction | ManagedTransaction) -> Node:
        return tx.run("""
            MERGE (e:Entity:PolypeptideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            SET e += $props,
                e.family = $family
            RETURN e
        """, {
            "rcsb_id": parent_rcsb_id,
            "entity_id": entity.entity_id,
            "props": props,
            "family": family_val,
        }).single(strict=True)['e']

    return _


def node__polynucleotide_entity(
    entity: PolynucleotideEntity,
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates/Merges a DNA/RNA Entity.
    """
    props = {
        "entity_id": entity.entity_id,
        "parent_rcsb_id": parent_rcsb_id,
        "pdbx_description": entity.pdbx_description,
        "formula_weight": entity.formula_weight,
        "pdbx_strand_ids": entity.pdbx_strand_ids,
        "polymer_type": entity.polymer_type,
        "one_letter_code": entity.one_letter_code,
        "one_letter_code_can": entity.one_letter_code_can,
        "sequence_length": entity.sequence_length,
        "src_organism_names": entity.src_organism_names,
        "src_organism_ids": entity.src_organism_ids,
    }
    props = {k: v for k, v in props.items() if v is not None}

    def _(tx: Transaction | ManagedTransaction) -> Node:
        return tx.run("""
            MERGE (e:Entity:PolynucleotideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            SET e += $props
            RETURN e
        """, {
            "rcsb_id": parent_rcsb_id,
            "entity_id": entity.entity_id,
            "props": props,
        }).single(strict=True)['e']

    return _


def node__polymer_instance(
    instance: Polypeptide | Polynucleotide
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates a polymer instance node using asym_id as unique key.
    Links to parent entity.
    """
    specific_label = "PolypeptideInstance" if isinstance(instance, Polypeptide) else "PolynucleotideInstance"

    def _(tx: Transaction | ManagedTransaction) -> Node:
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
            "rcsb_id": instance.parent_rcsb_id,
            "asym_id": instance.asym_id,
            "auth_id": instance.auth_asym_id,
            "assembly_id": instance.assembly_id,
            "entity_id": instance.entity_id,
        }).single(strict=True)['i']

    return _


def link__entity_to_structure(
    entity_node: Node,
    parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Links an Entity to the Structure root.
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