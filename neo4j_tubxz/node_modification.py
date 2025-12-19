from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from lib.types import Modification

def node__modification(
    mod: Modification
) -> Callable[[Transaction | ManagedTransaction], Node]:
    
    mod_props = mod.model_dump()
    
    def _(tx: Transaction | ManagedTransaction):
        # Assuming you want a unique constraint on master_index, uniprot_id, and type for modifications
        return tx.run("""
            MERGE (m:Modification {
                master_index: $master_index, 
                uniprot_id:   $uniprot_id, 
                modification_type: $modification_type
            })
            ON CREATE SET m += $properties
            ON MATCH SET m += $properties
            RETURN m
            """, {
                "master_index":      mod_props['master_index'],
                "uniprot_id":        mod_props['uniprot_id'],
                "modification_type": mod_props['modification_type'],
                "properties":        mod_props
            }).single(strict=True)['m']
    return _