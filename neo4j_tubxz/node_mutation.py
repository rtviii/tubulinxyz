from typing import Callable
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node
from lib.types import Mutation

def node__mutation(
    mut: Mutation
) -> Callable[[Transaction | ManagedTransaction], Node]:
    """
    Merges a mutation node based on the unique NODE KEY.
    The relationship to the MasterAlignment is omitted to avoid supernodes.
    """
    mut_props = mut.model_dump()
    
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""
            MERGE (m:Mutation {
                master_index: $master_index, 
                uniprot_id:   $uniprot_id, 
                from_residue: $from_residue, 
                to_residue:   $to_residue
            })
            ON CREATE SET m += $properties
            ON MATCH   SET m += $properties
            RETURN m
            """, {
                "master_index": mut_props['master_index'],
                "uniprot_id":   mut_props['uniprot_id'],
                "from_residue": mut_props['from_residue'],
                "to_residue":   mut_props['to_residue'],
                "properties":   mut_props
            }).single(strict=True)['m']
    return _