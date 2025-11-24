from typing import Callable
import typing
from pydantic import BaseModel
from neo4j import Driver, Transaction, ManagedTransaction
from neo4j.graph import Node

from lib.etl.libtax import PhylogenyNode

# Create constraint on phylogenyNode that ncbi_tax_id is unique
def node__phylogeny(phylogeny_obj:PhylogenyNode)->Callable[[Transaction | ManagedTransaction], Node ]:
    def _(tx: Transaction | ManagedTransaction):
        # TUBE-FIX: Robust MERGE.
        # MERGE only on the unique constraint (ncbi_tax_id).
        # SET other properties on create OR match to keep them updated.
        return tx.run("""//
            MERGE (phylonode:PhylogenyNode { ncbi_tax_id: $ncbi_tax_id })
            SET
                phylonode.scientific_name = $scientific_name,
                phylonode.rank = $rank
            RETURN phylonode
            """, **phylogeny_obj.model_dump()).single(strict=True)['phylonode']
    return _

def node__phylogeny_exists(taxid:int)->Callable[[Transaction | ManagedTransaction], bool]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
            match (p:PhylogenyNode {ncbi_tax_id:$TAXID}) return count(p) > 0;
        """, {"TAXID":taxid}).single().value()
    return _

def link__phylogeny(target:int, parent:int)->Callable[[Transaction | ManagedTransaction], list[Node] ]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
                match (target:PhylogenyNode {ncbi_tax_id:$target_taxid}) 
                match (parent:PhylogenyNode {ncbi_tax_id:$parent_taxid})
                merge (target)-[d:descendant_of]->(parent)
                return target, parent
        """, {"target_taxid": target, "parent_taxid": parent}).values()
    return _
