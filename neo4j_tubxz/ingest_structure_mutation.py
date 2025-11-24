# neo4j_tubxz/ingest_structure_mutations.py

"""
This will:

Find all sequence_ingestion.json files in your TUBETL_DATA directory
For each structure:

Create/update the IS_MAPPED_TO relationship with the mapping arrays
Create all mutation nodes
Link mutations to polymers via HAS_MUTATION
Link mutations to master alignment via ANNOTATES_POSITION_IN

The alignment mappings are now stored on the relationship, so you can query them like:
```
MATCH (p:Polymer {parent_rcsb_id: "1SA0", auth_asym_id: "A"})-[r:IS_MAPPED_TO]->(ma:MasterAlignment)
RETURN r.seqres_to_master, r.master_to_seqres
```
"""

import json
from pathlib import Path
from typing import Dict, Any, List
from neo4j import ManagedTransaction, Transaction
from neo4j.graph import Node

from neo4j_tubxz.db_driver import Neo4jAdapter
from lib.models.types_tubulin import AlignmentMapping


def load_sequence_ingestion_data(json_path: str) -> Dict[str, Any]:
    """Load the sequence_ingestion.json file"""
    with open(json_path, 'r') as f:
        return json.load(f)


def get_polymer_node(
    adapter: Neo4jAdapter,
    rcsb_id: str,
    auth_asym_id: str
) -> Node | None:
    """Fetch the polymer node from the database"""
    with adapter.driver.session() as session:
        def _(tx: Transaction | ManagedTransaction):
            result = tx.run("""
                MATCH (p:Polymer {parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id})
                RETURN p
                """, {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id}
            ).single()
            return result['p'] if result else None
        return session.execute_read(_)


def get_master_alignment_node(
    adapter: Neo4jAdapter,
    family: str,
    version: str = "v1.0"
) -> Node | None:
    """Fetch the master alignment node"""
    with adapter.driver.session() as session:
        def _(tx: Transaction | ManagedTransaction):
            result = tx.run("""
                MATCH (ma:MasterAlignment {family: $family, version: $version})
                RETURN ma
                """, {"family": family.lower(), "version": version}
            ).single()
            return result['ma'] if result else None
        return session.execute_read(_)


def create_or_update_alignment_mapping(
    adapter: Neo4jAdapter,
    polymer_node: Node,
    ma_node: Node,
    ma_to_auth_map: List[int],
    observed_to_ma_map: List[int]
) -> None:
    """
    Create or update the alignment mapping relationship between polymer and master alignment.
    Stores the bidirectional mapping arrays.
    """
    with adapter.driver.session() as session:
        def _(tx: Transaction | ManagedTransaction):
            # Convert lists to JSON strings
            mapping = AlignmentMapping(
                seqres_to_master=json.dumps(observed_to_ma_map),
                master_to_seqres=json.dumps(ma_to_auth_map)
            )
            
            tx.run("""
                MATCH (p:Polymer) WHERE ELEMENTID(p) = $p_elem_id
                MATCH (ma:MasterAlignment) WHERE ELEMENTID(ma) = $ma_elem_id
                MERGE (p)-[r:IS_MAPPED_TO]->(ma)
                SET r.seqres_to_master = $seqres_to_master,
                    r.master_to_seqres = $master_to_seqres
                RETURN r
                """, {
                    "p_elem_id": polymer_node.element_id,
                    "ma_elem_id": ma_node.element_id,
                    "seqres_to_master": mapping.seqres_to_master,
                    "master_to_seqres": mapping.master_to_seqres
                }
            )
        session.execute_write(_)


# neo4j_tubxz/ingest_structure_mutations.py

def create_structure_mutation_node(
    adapter: Neo4jAdapter,
    mutation_data: Dict[str, Any],
    rcsb_id: str,
    auth_asym_id: str,
    family: str
) -> Node:
    """Create a structure-derived mutation node"""
    with adapter.driver.session() as session:
        def _(tx: Transaction | ManagedTransaction):
            props = {
                "master_index": mutation_data["ma_position"],
                "from_residue": mutation_data["wild_type"],
                "to_residue": mutation_data["observed"],
                "pdb_auth_id": mutation_data["pdb_auth_id"],
                "parent_rcsb_id": rcsb_id.upper(),
                "auth_asym_id": auth_asym_id,
                "source_type": "structure_alignment",
                "tubulin_family": family.lower(),
                "uniprot_id": "",
                "species": "",
                "tubulin_type": family,
                "phenotype": f"Structural variant at position {mutation_data['ma_position']}",
                "database_source": "PDB",
                "reference_link": f"https://www.rcsb.org/structure/{rcsb_id}",
                "keywords": "structural_alignment",
                "notes": f"Detected via MUSCLE alignment to master {family} profile"
            }
            
            return tx.run("""
                CREATE (m:Mutation)
                SET m = $properties
                RETURN m
                """, {"properties": props}
            ).single(strict=True)['m']
        
        # FIX: Use execute_write instead of execute_read
        return session.execute_write(_)

def link_mutation_to_polymer(
    adapter: Neo4jAdapter,
    mutation_node: Node,
    polymer_node: Node
) -> None:
    """Link a mutation to its polymer"""
    with adapter.driver.session() as session:
        def _(tx: Transaction | ManagedTransaction):
            return tx.run("""
                MATCH (m:Mutation) WHERE ELEMENTID(m) = $m_elem_id
                MATCH (p:Polymer) WHERE ELEMENTID(p) = $p_elem_id
                MERGE (p)-[r:HAS_MUTATION]->(m)
                RETURN r
                """, {
                    "m_elem_id": mutation_node.element_id,
                    "p_elem_id": polymer_node.element_id
                }
            ).single()
        session.execute_write(_)  # <-- WRITE not READ


def link_mutation_to_master_alignment(
    adapter: Neo4jAdapter,
    mutation_node: Node,
    ma_node: Node
) -> None:
    """Link a mutation to the master alignment"""
    with adapter.driver.session() as session:
        def _(tx: Transaction | ManagedTransaction):
            return tx.run("""
                MATCH (m:Mutation) WHERE ELEMENTID(m) = $m_elem_id
                MATCH (ma:MasterAlignment) WHERE ELEMENTID(ma) = $ma_elem_id
                MERGE (m)-[r:ANNOTATES_POSITION_IN]->(ma)
                RETURN r
                """, {
                    "m_elem_id": mutation_node.element_id,
                    "ma_elem_id": ma_node.element_id
                }
            ).single()
        session.execute_write(_)  # <-- WRITE not READ


def ingest_structure_mutations_for_chain(
    adapter: Neo4jAdapter,
    rcsb_id: str,
    chain_data: Dict[str, Any],
    version: str = "v1.0"
) -> Dict[str, int]:
    """
    Ingest all mutations for a single chain AND create alignment mapping.
    Returns stats about what was ingested.
    """
    stats = {
        "mutations_created": 0,
        "mutations_skipped": 0,
        "alignment_mapped": 0,
        "errors": 0
    }
    
    data = chain_data.get("data", {})
    auth_asym_id = data.get("chain_id")
    family = chain_data.get("family", "").lower()
    mutations = data.get("mutations", [])
    ma_to_auth_map = data.get("ma_to_auth_map", [])
    observed_to_ma_map = data.get("observed_to_ma_map", [])
    
    if not auth_asym_id or not family:
        print(f"⚠️ Missing chain_id or family for {rcsb_id}")
        stats["errors"] += 1
        return stats
    
    # Get polymer node
    polymer_node = get_polymer_node(adapter, rcsb_id, auth_asym_id)
    if not polymer_node:
        print(f"⚠️ Polymer not found: {rcsb_id}.{auth_asym_id}")
        stats["mutations_skipped"] = len(mutations)
        return stats
    
    # Get master alignment node
    ma_node = get_master_alignment_node(adapter, family, version)
    if not ma_node:
        print(f"⚠️ Master alignment not found for family '{family}'")
        stats["mutations_skipped"] = len(mutations)
        return stats
    
    # Create/update alignment mapping
    if ma_to_auth_map and observed_to_ma_map:
        try:
            create_or_update_alignment_mapping(
                adapter, polymer_node, ma_node,
                ma_to_auth_map, observed_to_ma_map
            )
            stats["alignment_mapped"] = 1
            print(f"  ✓ Created alignment mapping")
        except Exception as e:
            print(f"  ⚠️ Failed to create alignment mapping: {e}")
            stats["errors"] += 1
    
    # Create and link each mutation
    print(f"  Processing {len(mutations)} mutations...")
    for mut_data in mutations:
        try:
            mut_node = create_structure_mutation_node(
                adapter, mut_data, rcsb_id, auth_asym_id, family
            )
            
            link_mutation_to_polymer(adapter, mut_node, polymer_node)
            link_mutation_to_master_alignment(adapter, mut_node, ma_node)
            
            stats["mutations_created"] += 1
            
        except Exception as e:
            print(f"  ⚠️ Failed to create mutation at MA pos {mut_data['ma_position']}: {e}")
            stats["errors"] += 1
    
    return stats


def ingest_structure_mutations_from_file(
    adapter: Neo4jAdapter,
    json_path: str,
    rcsb_id: str,
    version: str = "v1.0"
) -> Dict[str, Any]:
    """
    Ingest all structure mutations from a sequence_ingestion.json file.
    """
    print(f"\n{'='*60}")
    print(f"INGESTING STRUCTURE MUTATIONS: {rcsb_id}")
    print(f"Source: {json_path}")
    print(f"{'='*60}")
    
    data = load_sequence_ingestion_data(json_path)
    
    total_stats = {
        "chains_processed": 0,
        "mutations_created": 0,
        "mutations_skipped": 0,
        "alignments_mapped": 0,
        "errors": 0
    }
    
    for auth_asym_id, chain_data in data.items():
        print(f"\nChain {auth_asym_id}:")
        stats = ingest_structure_mutations_for_chain(
            adapter, rcsb_id, chain_data, version
        )
        
        total_stats["chains_processed"] += 1
        total_stats["mutations_created"] += stats["mutations_created"]
        total_stats["mutations_skipped"] += stats["mutations_skipped"]
        total_stats["alignments_mapped"] += stats["alignment_mapped"]
        total_stats["errors"] += stats["errors"]
        
        print(f"  Mutations: {stats['mutations_created']} created, {stats['mutations_skipped']} skipped")
        if stats["errors"] > 0:
            print(f"  Errors: {stats['errors']}")
    
    print(f"\n{'='*60}")
    print(f"SUMMARY:")
    print(f"  Chains processed: {total_stats['chains_processed']}")
    print(f"  Mutations created: {total_stats['mutations_created']}")
    print(f"  Mutations skipped: {total_stats['mutations_skipped']}")
    print(f"  Alignments mapped: {total_stats['alignments_mapped']}")
    print(f"  Errors: {total_stats['errors']}")
    print(f"{'='*60}\n")
    
    return total_stats


def batch_ingest_from_directory(
    adapter: Neo4jAdapter,
    base_dir: str,
    version: str = "v1.0"
) -> Dict[str, Any]:
    """
    Iterate over all structures in TUBETL_DATA and ingest their sequence_ingestion.json files.
    
    Expected structure:
    base_dir/
      1SA0/
        sequence_ingestion.json
      5JCO/
        sequence_ingestion.json
      ...
    """
    base_path = Path(base_dir)
    
    if not base_path.exists():
        raise ValueError(f"Base directory does not exist: {base_dir}")
    
    print(f"\n{'='*80}")
    print(f"BATCH INGESTION FROM: {base_dir}")
    print(f"{'='*80}\n")
    
    grand_total = {
        "structures_processed": 0,
        "structures_skipped": 0,
        "chains_processed": 0,
        "mutations_created": 0,
        "mutations_skipped": 0,
        "alignments_mapped": 0,
        "errors": 0
    }
    
    # Find all sequence_ingestion.json files
    ingestion_files = list(base_path.glob("*/sequence_ingestion.json"))
    
    print(f"Found {len(ingestion_files)} structures with sequence ingestion data\n")
    
    for ingestion_file in sorted(ingestion_files):
        rcsb_id = ingestion_file.parent.name.upper()
        
        try:
            stats = ingest_structure_mutations_from_file(
                adapter, str(ingestion_file), rcsb_id, version
            )
            
            grand_total["structures_processed"] += 1
            grand_total["chains_processed"] += stats["chains_processed"]
            grand_total["mutations_created"] += stats["mutations_created"]
            grand_total["mutations_skipped"] += stats["mutations_skipped"]
            grand_total["alignments_mapped"] += stats["alignments_mapped"]
            grand_total["errors"] += stats["errors"]
            
        except Exception as e:
            print(f"⚠️ FAILED TO PROCESS {rcsb_id}: {e}\n")
            grand_total["structures_skipped"] += 1
            grand_total["errors"] += 1
    
    print(f"\n{'='*80}")
    print(f"GRAND TOTAL:")
    print(f"  Structures processed: {grand_total['structures_processed']}")
    print(f"  Structures skipped: {grand_total['structures_skipped']}")
    print(f"  Chains processed: {grand_total['chains_processed']}")
    print(f"  Mutations created: {grand_total['mutations_created']}")
    print(f"  Alignments mapped: {grand_total['alignments_mapped']}")
    print(f"  Errors: {grand_total['errors']}")
    print(f"{'='*80}\n")
    
    return grand_total


def main():
    """Batch ingest all structure mutations"""
    from lib.etl.constants import NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD
    import os
    
    # Path to your TUBETL_DATA directory
    TUBETL_DATA = os.environ.get("TUBETL_DATA", "/path/to/TUBETL_DATA")
    
    adapter = Neo4jAdapter(
        uri=NEO4J_URI,
        user=NEO4J_USER,
        current_db=NEO4J_CURRENTDB,
        password=NEO4J_PASSWORD
    )
    
    try:
        batch_ingest_from_directory(adapter, TUBETL_DATA)
    finally:
        adapter.driver.close()


if __name__ == "__main__":
    main()