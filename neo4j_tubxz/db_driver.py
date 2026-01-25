# neo4j_tubxz/db_driver.py
import os
import sys
from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
from lib.etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, TUBETL_DATA
from lib.etl.assets import GlobalOps
from neo4j_tubxz.db_lib_builder import Neo4jAdapter

sys.dont_write_bytecode = True

def full_upload(overwrite: bool = True):
    """
    Recipe for initializing a new instance and ingesting all profiles on disk.
    """
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
    
    print("Initializing Neo4j instance (constraints and phylogeny)...")
    adapter.initialize_new_instance()

    all_structs_on_disk = GlobalOps.list_profiles()
    print(f"Found {len(all_structs_on_disk)} profiles on disk. Starting ingestion...")

    futures: list[Future] = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        for rcsb_id in sorted(all_structs_on_disk):
            fut = executor.submit(partial(adapter.add_structure, rcsb_id, overwrite))
            futures.append(fut)
    
    wait(futures, return_when=ALL_COMPLETED)
    print("\nFull upload complete.")
    adapter.close()

def rcsb_sync():
    """
    Compare local database against RCSB search results and sync missing structures.
    """
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
    
    # 1. Get existing IDs in DB
    db_ids = set(adapter.get_all_structure_ids())
    
    # 2. Get missing IDs (RCSB - DB)
    missing_ids = GlobalOps.missing_db_entries(list(db_ids))
    
    if not missing_ids:
        print("Database is already up to date with RCSB.")
        return

    print(f"Sync: {len(missing_ids)} structures found in RCSB but missing in DB.")
    
    # NOTE: This assumes ETL profiles have been generated already.
    # If not, you'd trigger the ETLCollector here.
    
    futures: list[Future] = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        for rcsb_id in sorted(missing_ids):
            fut = executor.submit(partial(adapter.add_structure, rcsb_id, False))
            futures.append(fut)
    
    wait(futures, return_when=ALL_COMPLETED)
    adapter.close()

if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1 and sys.argv[1] == "sync":
        rcsb_sync()
    else:
        full_upload(overwrite=True)
