# neo4j_tubxz/node_modification.py
"""
Modification (PTM) node creation.
"""
from typing import Callable, List
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node

from lib.types import Modification


def node__modification(
    mod: Modification,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """Create or update a single Modification node."""

    mod_props = mod.model_dump(exclude_none=True)

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MERGE (m:Modification {
                master_index:      $master_index,
                uniprot_id:        $uniprot_id,
                modification_type: $modification_type
            })
            ON CREATE SET m += $properties
            ON MATCH SET m += $properties
            RETURN m
        """, {
            "master_index":      mod_props["master_index"],
            "uniprot_id":        mod_props["uniprot_id"],
            "modification_type": mod_props["modification_type"],
            "properties":        mod_props,
        }).single(strict=True)["m"]

    return _


def batch_create_modifications(
    tx: Transaction | ManagedTransaction,
    modifications: List[Modification],
) -> int:
    """
    Batch create Modification nodes for literature-sourced PTMs.
    Uses MERGE on (master_index, uniprot_id, modification_type) to avoid duplicates.
    """
    if not modifications:
        return 0

    props_list = [m.model_dump(exclude_none=True) for m in modifications]

    tx.run("""
        UNWIND $props_list AS props
        MERGE (m:Modification {
            master_index:      props.master_index,
            uniprot_id:        props.uniprot_id,
            modification_type: props.modification_type
        })
        SET m += props
    """, {"props_list": props_list})

    return len(props_list)
