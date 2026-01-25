# neo4j_tubxz/node_ligand.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from lib.types import NonpolymerEntity, Nonpolymer


def node__chemical(
    entity: NonpolymerEntity,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Merges the Global Chemical Node (shared across ALL PDB structures).
    Uses the chemical_id (e.g., 'GTP', 'TAX') as the unique key.
    """
    props = {
        "chemical_id": entity.chemical_id,
        "chemical_name": entity.chemical_name,
        "SMILES": entity.SMILES,
        "SMILES_stereo": entity.SMILES_stereo,
        "InChI": entity.InChI,
        "InChIKey": entity.InChIKey,
        # Flatten drugbank info if available
        "drugbank_id": (
            entity.nonpolymer_comp.drugbank.drugbank_container_identifiers.drugbank_id
            if entity.nonpolymer_comp
            and entity.nonpolymer_comp.drugbank
            and entity.nonpolymer_comp.drugbank.drugbank_container_identifiers
            else None
        ),
    }
    # Clean None values
    props = {k: v for k, v in props.items() if v is not None}

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """
            MERGE (c:Chemical {chemical_id: $chemical_id})
            ON CREATE SET c += $props
            ON MATCH SET 
                c.chemical_name = COALESCE(c.chemical_name, $props.chemical_name),
                c.SMILES = COALESCE(c.SMILES, $props.SMILES),
                c.InChIKey = COALESCE(c.InChIKey, $props.InChIKey)
            RETURN c
        """,
            {"chemical_id": entity.chemical_id, "props": props},
        ).single(strict=True)["c"]

    return _


def node__nonpolymer_entity(
    entity: NonpolymerEntity, parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Merges the Structure-Specific Entity Definition.
    Links to the Global Chemical node.
    """
    props = {
        "entity_id": entity.entity_id,
        "parent_rcsb_id": parent_rcsb_id,
        "pdbx_description": entity.pdbx_description,
        "formula_weight": entity.formula_weight,
        "chemical_id": entity.chemical_id,
    }

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """
            MERGE (e:Entity:NonpolymerEntity {parent_rcsb_id: $parent_rcsb_id, entity_id: $entity_id})
            ON CREATE SET e += $props
            ON MATCH SET e += $props
            
            WITH e
            MATCH (c:Chemical {chemical_id: $chemical_id})
            MERGE (e)-[:DEFINED_BY_CHEMICAL]->(c)
            RETURN e
        """,
            {
                "parent_rcsb_id": parent_rcsb_id,
                "entity_id": entity.entity_id,
                "chemical_id": entity.chemical_id,
                "props": props,
            },
        ).single(strict=True)["e"]

    return _


# node_ligand.py -> node__nonpolymer_instance
def node__nonpolymer_instance(instance: Nonpolymer):
    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """
            MERGE (i:Instance:NonpolymerInstance {parent_rcsb_id: $rcsb_id, asym_id: $asym_id})
            ON CREATE SET 
                i.auth_asym_id = $auth_id,
                i.auth_seq_id  = $seq_id,
                i.assembly_id  = $assembly_id
            
            WITH i
            MATCH (e:Entity:NonpolymerEntity {parent_rcsb_id: $rcsb_id, entity_id: $entity_id})
            MERGE (i)-[:INSTANCE_OF]->(e)
            RETURN i
        """,
            {
                "rcsb_id": instance.parent_rcsb_id,
                "asym_id": instance.asym_id,
                "auth_id": instance.auth_asym_id,
                "seq_id": instance.auth_seq_id,
                "assembly_id": instance.assembly_id,
                "entity_id": instance.entity_id,
            },
        ).single(strict=True)["i"]

    return _
