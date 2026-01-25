# tubexyz/neo4j_driver/node_structure.py
from typing import Callable, Literal
from neo4j import ManagedTransaction, Record, Transaction
from neo4j.graph import Node

from lib.types import TubulinStructure

def struct_exists(rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], bool]:
    def _(tx: Transaction | ManagedTransaction):
        return (
            tx.run(
                """//
                MATCH (u:Structure {rcsb_id: $rcsb_id})
                return COUNT(u) > 0;
                """,
                parameters={"rcsb_id": rcsb_id},
            )
            .single()
            .value()
        )
    return _

def link__structure_to_lineage_member(
    rcsb_id: str,
    taxid,
    relationship: Literal["belongs_to_lineage_host", "belongs_to_lineage_source"],
) -> Callable[[Transaction | ManagedTransaction], list[Node]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
            match (struct:Structure {{rcsb_id: $rcsb_id}})
            match (phylo:PhylogenyNode {{ncbi_tax_id: $tax_id}}) 
            merge (struct)<-[organism: {} ]-(phylo)
            return  struct, phylo
            """.format(relationship),
            { "rcsb_id": rcsb_id, "tax_id": taxid, },
        ).values("struct", "phylo")
    return _


def link__structure_to_organism(
    rcsb_id: str, taxid, relationship: Literal["host", "source"]
) -> Callable[[Transaction | ManagedTransaction], list[Node]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
            match (struct:Structure {{rcsb_id: $rcsb_id}})
            match (phylo:PhylogenyNode {{ncbi_tax_id: $tax_id}}) 
            merge (struct)<-[organism: {} ]-(phylo)
            return  struct, phylo
            """.format(relationship),
            { "rcsb_id": rcsb_id, "tax_id": taxid, },
        ).values("struct", "phylo")
    return _


def node__structure(
    _struct: TubulinStructure,
) -> Callable[[Transaction | ManagedTransaction], Record | None]:
    # Only extract the fields needed for the Structure node
    params = {
        "rcsb_id"               : _struct.rcsb_id,
        "expMethod"             : _struct.expMethod,
        "resolution"            : _struct.resolution,
        "pdbx_keywords"         : _struct.pdbx_keywords,
        "pdbx_keywords_text"    : _struct.pdbx_keywords_text,
        "src_organism_ids"      : _struct.src_organism_ids,
        "src_organism_names"    : _struct.src_organism_names,
        "host_organism_ids"     : _struct.host_organism_ids,
        "host_organism_names"   : _struct.host_organism_names,
        "rcsb_external_ref_id"  : _struct.rcsb_external_ref_id,
        "rcsb_external_ref_type": _struct.rcsb_external_ref_type,
        "rcsb_external_ref_link": _struct.rcsb_external_ref_link,
        "citation_pdbx_doi"     : _struct.citation_pdbx_doi,
        "citation_year"         : _struct.citation_year,
        "citation_title"        : _struct.citation_title,
        "citation_rcsb_authors" : _struct.citation_rcsb_authors,
        "deposition_date"       : _struct.deposition_date,
        "polymerization_state"  : _struct.polymerization_state,
    }

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
            MERGE (struct:Structure { rcsb_id: $rcsb_id })
            SET
                struct.expMethod             = $expMethod,
                struct.resolution            = $resolution,
                struct.pdbx_keywords         = $pdbx_keywords,
                struct.pdbx_keywords_text    = $pdbx_keywords_text,
                
                // --- RESTORED PROPERTIES ---
                struct.src_organism_ids      = $src_organism_ids,
                struct.src_organism_names    = $src_organism_names,
                struct.host_organism_ids     = $host_organism_ids,
                struct.host_organism_names   = $host_organism_names,
                
                struct.rcsb_external_ref_id    = $rcsb_external_ref_id,
                struct.rcsb_external_ref_type = $rcsb_external_ref_type,
                struct.rcsb_external_ref_link = $rcsb_external_ref_link,
                
                struct.citation_pdbx_doi      = $citation_pdbx_doi,
                struct.citation_year          = $citation_year,
                struct.citation_title         = $citation_title,
                struct.citation_rcsb_authors  = $citation_rcsb_authors,
                
                struct.deposition_date        = $deposition_date,
                struct.polymerization_state   = $polymerization_state
            RETURN struct
            """,
            **params,
        ).single()

    return _