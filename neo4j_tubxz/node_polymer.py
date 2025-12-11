# tubexyz/neo4j_driver/node_polymer.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from lib.models.types_tubulin import PolymerEntity, TubulinEntity, PolymerInstance

# --- ENTITY NODES ---


def node__entity(
    entity: PolymerEntity,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates or merges the PolymerEntity node.
    Contains Sequence, Taxonomy, and Descriptions.
    """
    E = entity.model_dump(
        exclude={
            "family",
            "pfam_accessions",
            "pfam_comments",
            "pfam_descriptions",
            "uniprot_accession",
        }
    )

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
    MERGE (e:Entity {
        parent_rcsb_id: $parent_rcsb_id,
        entity_id     : $entity_id
    })
    ON CREATE SET
        e.rcsb_pdbx_description                = $rcsb_pdbx_description,
        e.src_organism_names                   = $src_organism_names,
        e.host_organism_names                  = $host_organism_names,
        e.src_organism_ids                     = $src_organism_ids,
        e.host_organism_ids                    = $host_organism_ids,
        
        e.entity_poly_seq_one_letter_code      = $entity_poly_seq_one_letter_code,
        e.entity_poly_seq_one_letter_code_can  = $entity_poly_seq_one_letter_code_can,
        e.entity_poly_seq_length               = $entity_poly_seq_length,
        e.entity_poly_polymer_type             = $entity_poly_polymer_type,
        e.entity_poly_entity_type              = $entity_poly_entity_type
    RETURN e
    """,
            **E,
        ).single(strict=True)["e"]

    return _


def upsert__entity_to_tubulin(
    entity_node: Node, tubulin_entity: TubulinEntity
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Adds :Tubulin label and specific properties to the Entity.
    """
    props = tubulin_entity.model_dump(
        include={
            "family",
            "pfam_accessions",
            "pfam_comments",
            "pfam_descriptions",
            "uniprot_accession",
        }
    )
    if props["family"]:
        props["family"] = props["family"].value

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
        MATCH (e:Entity) WHERE ELEMENTID(e) = $elem_id
        SET
            e:Tubulin,
            e.family            = $family,
            e.pfam_accessions   = $pfam_accessions,
            e.pfam_comments     = $pfam_comments,
            e.pfam_descriptions = $pfam_descriptions,
            e.uniprot_accession = $uniprot_accession
        RETURN e
        """,
            {"elem_id": entity_node.element_id, **props},
        ).single(strict=True)["e"]

    return _


# --- INSTANCE NODES ---


def node__instance(
    instance: PolymerInstance,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Creates the lightweight Instance node.
    """
    I = instance.model_dump()

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
    MERGE (i:Instance {
        parent_rcsb_id: $parent_rcsb_id,
        auth_asym_id  : $auth_asym_id
    })
    ON CREATE SET
        i.assembly_id           = $assembly_id,
        i.asym_ids              = $asym_ids,
        i.entity_poly_strand_id = $entity_poly_strand_id
    RETURN i
    """,
            **I,
        ).single(strict=True)["i"]

    return _


# --- LINKS ---


def link__instance_to_structure(
    instance_node: Node, parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], list]:
    def _(tx: Transaction | ManagedTransaction):
        # Relationship: Structure HAS_INSTANCE Instance
        return tx.run(
            """//
            MATCH (i:Instance) WHERE ELEMENTID(i) = $ELEM_ID
            MATCH (s:Structure {rcsb_id: $PARENT})
            MERGE (s)-[rel:HAS_INSTANCE]->(i)
            RETURN s, rel, i
            """,
            {"ELEM_ID": instance_node.element_id, "PARENT": parent_rcsb_id},
        ).values("s", "rel", "i")

    return _


def link__instance_to_entity(
    instance_node: Node, entity_id: str, parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], list]:
    def _(tx: Transaction | ManagedTransaction):
        # Relationship: Instance IS_INSTANCE_OF Entity
        return tx.run(
            """//
            MATCH (i:Instance) WHERE ELEMENTID(i) = $I_ID
            MATCH (e:Entity {parent_rcsb_id: $PARENT, entity_id: $E_ID})
            MERGE (i)-[rel:IS_INSTANCE_OF]->(e)
            RETURN i, rel, e
            """,
            {
                "I_ID": instance_node.element_id,
                "PARENT": parent_rcsb_id,
                "E_ID": entity_id,
            },
        ).values("i", "rel", "e")

    return _
