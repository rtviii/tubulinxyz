# tubexyz/neo4j_driver/node_ligand.py
from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship

from models.types_tubulin import NonpolymericLigand
# TUBE-UPDATE: Import new tubulin schema

def node__ligand(
    _ligand: NonpolymericLigand,
    parent_rcsb_id: str # TUBE-UPDATE: Added for bug fix
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    def _(tx: Transaction | ManagedTransaction):
        # Prepare the properties dictionary
        properties = {
            "chemicalId"         : _ligand.chemicalId,
            "chemicalName"       : _ligand.chemicalName,
            "formula_weight"     : _ligand.formula_weight,
            "pdbx_description"   : _ligand.pdbx_description,
            "number_of_instances": _ligand.number_of_instances,
            "drugbank_id"        : None,
            "drugbank_description": None,
            "SMILES"             : _ligand.SMILES,
            "SMILES_stereo"      : _ligand.SMILES_stereo,
            "InChI"              : _ligand.InChI,
            "InChIKey"           : _ligand.InChIKey,
        }

        # Safe dictionary navigation
        if _ligand.nonpolymer_comp and _ligand.nonpolymer_comp.drugbank:
            if _ligand.nonpolymer_comp.drugbank.drugbank_container_identifiers:
                properties["drugbank_id"] = _ligand.nonpolymer_comp.drugbank.drugbank_container_identifiers.drugbank_id
            if _ligand.nonpolymer_comp.drugbank.drugbank_info:
                properties["drugbank_description"] = _ligand.nonpolymer_comp.drugbank.drugbank_info.description

        properties = {k: v for k, v in properties.items() if v is not None}

        # TUBE-UPDATE: Implemented the robust MERGE query from your notes
        query = """
        MERGE (ligand:Ligand {chemicalId: $chemicalId})
        ON CREATE SET
            ligand += $properties,
            ligand.seen_in_structures = [$parent_rcsb_id]
        ON MATCH SET
            // Only update properties if they are currently null
            ligand.pdbx_description = COALESCE(ligand.pdbx_description, $properties.pdbx_description),
            ligand.SMILES = COALESCE(ligand.SMILES, $properties.SMILES),
            ligand.chemicalName = COALESCE(ligand.chemicalName, $properties.chemicalName),
            
            // Append to lists
            ligand.seen_in_structures = CASE
                WHEN $parent_rcsb_id IS NOT NULL AND NOT $parent_rcsb_id IN ligand.seen_in_structures
                THEN ligand.seen_in_structures + $parent_rcsb_id
                ELSE ligand.seen_in_structures
            END
        RETURN ligand
        """
        
        result = tx.run(query, {
            "chemicalId": _ligand.chemicalId,
            "properties": properties,
            "parent_rcsb_id": parent_rcsb_id
        })
        return result.single(strict=True)["ligand"]

    return _

def link__ligand_to_struct(
    ligand_node: Node, parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    parent_rcsb_id = parent_rcsb_id.upper()

    def _(tx: Transaction | ManagedTransaction):
        # TUBE-UPDATE: Renamed to Structure, relationship is CONTAINS
        return tx.run(
            """//
            MATCH (ligand:Ligand) WHERE ELEMENTID(ligand) = $ELEM_ID
            MATCH (struct:Structure {rcsb_id: $PARENT})
            MERGE (ligand)<-[contains:CONTAINS]-(struct)
            RETURN struct, ligand, contains
            """,
            {"ELEM_ID": ligand_node.element_id, "PARENT": parent_rcsb_id},
        ).values("struct", "ligand", "contains")

    return _