# tubexyz/neo4j_driver/node_polymer.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from lib.models.types_tubulin import Polymer, TubulinProtein

def node__polymer(poly:Polymer)->Callable[[Transaction | ManagedTransaction], Node ]:
    """
    Creates or merges the base Polymer node.
    This works for both base Polymers and TubulinProteins.
    """
    P = poly.model_dump()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
    MERGE (poly:Polymer {
        parent_rcsb_id: $parent_rcsb_id,
        auth_asym_id  : $auth_asym_id
    })
    ON CREATE SET
        poly.entity_id                           = $entity_id,
        poly.assembly_id                         = $assembly_id,
        poly.asym_ids                            = $asym_ids,
        poly.src_organism_names                  = $src_organism_names,
        poly.host_organism_names                 = $host_organism_names,
        poly.src_organism_ids                    = $src_organism_ids,
        poly.host_organism_ids                   = $host_organism_ids,
        poly.rcsb_pdbx_description               = $rcsb_pdbx_description,
        poly.entity_poly_strand_id               = $entity_poly_strand_id,
        poly.entity_poly_seq_one_letter_code     = $entity_poly_seq_one_letter_code,
        poly.entity_poly_seq_one_letter_code_can = $entity_poly_seq_one_letter_code_can,
        poly.entity_poly_seq_length              = $entity_poly_seq_length,
        poly.entity_poly_polymer_type            = $entity_poly_polymer_type,
        poly.entity_poly_entity_type             = $entity_poly_entity_type
    ON MATCH SET
        poly.entity_id = $entity_id,
        poly.entity_poly_seq_length = $entity_poly_seq_length
        // Add other ON MATCH properties as needed
    RETURN poly
    """, **P).single(strict=True)['poly']
    
    return _

def upsert__polymer_to_protein(
    polymer_node: Node, protein: TubulinProtein
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Adds the :Protein label and protein-specific properties
    to an existing Polymer node.
    """
    # Get only the properties specific to TubulinProtein
    prot_props = protein.model_dump(
        include={'family', 'pfam_accessions', 'pfam_comments', 'pfam_descriptions', 'uniprot_accession'}
    )
    # Convert Enum to string
    if prot_props['family']:
        prot_props['family'] = prot_props['family'].value

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        MATCH (p:Polymer) WHERE ELEMENTID(p) = $elem_id
        SET
            p:Protein,
            p.family = $family,
            p.pfam_accessions = $pfam_accessions,
            p.pfam_comments = $pfam_comments,
            p.pfam_descriptions = $pfam_descriptions,
            p.uniprot_accession = $uniprot_accession
        RETURN p
        """, {
            "elem_id": polymer_node.element_id,
            **prot_props
        }).single(strict=True)['p']
    return _


def link__polymer_to_structure(polymer_node: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    def _(tx: Transaction | ManagedTransaction):
        # : Relationship is HAS_POLYMER
        return tx.run("""//
            MATCH (polymer:Polymer) WHERE ELEMENTID(polymer) = $ELEM_ID
            MATCH (struct:Structure {rcsb_id: $PARENT})
            MERGE (struct)-[rel:HAS_POLYMER]->(polymer)
            RETURN polymer, rel, struct
            """,
            {
                "ELEM_ID": polymer_node.element_id,
                "PARENT": parent_rcsb_id
            }).values('polymer', 'rel', 'struct')
    return _