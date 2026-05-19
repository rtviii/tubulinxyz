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
    """Create or update a single Modification node and its OCCURS_IN edge."""

    mod_props = mod.model_dump(exclude_none=True)

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MERGE (m:Modification {
                master_index:      $master_index,
                uniprot_id:        $uniprot_id,
                modification_type: $modification_type
            })
            ON CREATE SET m += $properties
            ON MATCH  SET m += $properties
            WITH m
            MERGE (p:PhylogenyNode {ncbi_tax_id: $tax_id})
              ON CREATE SET p.scientific_name = $species_full_name,
                            p.rank = 'species'
            MERGE (m)-[:OCCURS_IN]->(p)
            RETURN m
        """, {
            "master_index":      mod_props["master_index"],
            "uniprot_id":        mod_props["uniprot_id"],
            "modification_type": mod_props["modification_type"],
            "tax_id":            mod_props["tax_id"],
            "species_full_name": mod_props["species_full_name"],
            "properties":        mod_props,
        }).single(strict=True)["m"]

    return _


def batch_create_modifications(
    tx: Transaction | ManagedTransaction,
    modifications: List[Modification],
) -> int:
    """
    Batch create Modification nodes for literature-sourced PTMs.
    MERGE key: (master_index, uniprot_id, modification_type). Each Modification
    also gets an :OCCURS_IN edge to a :PhylogenyNode keyed by tax_id; the
    PhylogenyNode is MERGE'd on the spot (created with scientific_name + rank
    only on first sight, so existing taxonomy nodes from structure ingestion
    are reused unchanged).
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
        WITH m, props
        MERGE (p:PhylogenyNode {ncbi_tax_id: props.tax_id})
          ON CREATE SET p.scientific_name = props.species_full_name,
                        p.rank = 'species'
        MERGE (m)-[:OCCURS_IN]->(p)
    """, {"props_list": props_list})

    return len(props_list)
