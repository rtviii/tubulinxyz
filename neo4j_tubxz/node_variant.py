# neo4j_tubxz/node_variant.py
"""
Variant node creation (Substitution, Insertion, Deletion).
"""
from typing import Callable, List
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node

from lib.types import SequenceVariant, VariantType


def node__variant(
    variant: SequenceVariant,
    parent_rcsb_id: str,
    entity_id: str,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates a Variant node with appropriate label (Substitution/Insertion/Deletion).
    Links it to the parent PolypeptideEntity.
    
    Node labels: Variant + one of (Substitution, Insertion, Deletion)
    """
    # Determine the specific label
    label_map = {
        VariantType.SUBSTITUTION: "Substitution",
        VariantType.INSERTION: "Insertion",
        VariantType.DELETION: "Deletion",
    }
    specific_label = label_map[variant.type]
    
    # Build properties dict, excluding None values
    props = {
        "type": variant.type.value,
        "source": variant.source,
        "master_index": variant.master_index,
        "observed_index": variant.observed_index,
        "wild_type": variant.wild_type,
        "observed": variant.observed,
        "uniprot_id": variant.uniprot_id,
        "phenotype": variant.phenotype,
        "reference": variant.reference,
    }
    props = {k: v for k, v in props.items() if v is not None}

    def _(tx: Transaction | ManagedTransaction) -> Node:
        # Create variant node and link to entity
        result = tx.run(f"""
            MATCH (e:Entity:PolypeptideEntity {{parent_rcsb_id: $rcsb_id, entity_id: $entity_id}})
            CREATE (v:Variant:{specific_label} $props)
            CREATE (e)-[:HAS_VARIANT]->(v)
            RETURN v
        """, {
            "rcsb_id": parent_rcsb_id,
            "entity_id": entity_id,
            "props": props,
        })
        return result.single(strict=True)['v']

    return _


def create_variants_for_entity(
    tx: Transaction | ManagedTransaction,
    variants: List[SequenceVariant],
    parent_rcsb_id: str,
    entity_id: str,
) -> int:
    """
    Batch create all variants for an entity.
    Returns the count of created variants.
    """
    if not variants:
        return 0

    # Group by type for batch creation
    substitutions = [v for v in variants if v.type == VariantType.SUBSTITUTION]
    insertions = [v for v in variants if v.type == VariantType.INSERTION]
    deletions = [v for v in variants if v.type == VariantType.DELETION]

    count = 0

    if substitutions:
        props_list = [_variant_to_props(v) for v in substitutions]
        tx.run("""
            MATCH (e:Entity:PolypeptideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            UNWIND $props_list AS props
            CREATE (v:Variant:Substitution)
            SET v = props
            CREATE (e)-[:HAS_VARIANT]->(v)
        """, {"rcsb_id": parent_rcsb_id, "entity_id": entity_id, "props_list": props_list})
        count += len(substitutions)

    if insertions:
        props_list = [_variant_to_props(v) for v in insertions]
        tx.run("""
            MATCH (e:Entity:PolypeptideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            UNWIND $props_list AS props
            CREATE (v:Variant:Insertion)
            SET v = props
            CREATE (e)-[:HAS_VARIANT]->(v)
        """, {"rcsb_id": parent_rcsb_id, "entity_id": entity_id, "props_list": props_list})
        count += len(insertions)

    if deletions:
        props_list = [_variant_to_props(v) for v in deletions]
        tx.run("""
            MATCH (e:Entity:PolypeptideEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            UNWIND $props_list AS props
            CREATE (v:Variant:Deletion)
            SET v = props
            CREATE (e)-[:HAS_VARIANT]->(v)
        """, {"rcsb_id": parent_rcsb_id, "entity_id": entity_id, "props_list": props_list})
        count += len(deletions)

    return count


def _variant_to_props(v: SequenceVariant) -> dict:
    """Convert variant to Neo4j properties dict."""
    props = {
        "type": v.type.value,
        "source": v.source,
        "master_index": v.master_index,
        "observed_index": v.observed_index,
        "wild_type": v.wild_type,
        "observed": v.observed,
        "uniprot_id": v.uniprot_id,
        "phenotype": v.phenotype,
        "reference": v.reference,
    }
    return {k: v for k, v in props.items() if v is not None}