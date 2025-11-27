# tubexyz/neo4j_driver/db_driver.py
from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
import os
import sys
from lib.etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_tubxz.db_lib_builder import Neo4jAdapter
sys.dont_write_bytecode = True

def full_upload():
    """
    Recipe for initializing a new instance from the data pool.
    """
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
    adapter.initialize_new_instance()
    futures: list[Future] = []

    # Get list of profiles on disk
    profile_dir = "tubet_data" # Or from etl.constants
    all_structs_on_disk = [
        d for d in os.listdir(profile_dir) 
        if os.path.isdir(os.path.join(profile_dir, d)) and len(d) == 4
    ]

    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in sorted(all_structs_on_disk):
            fut = executor.submit(partial(adapter.add_total_structure, rcsb_id, True))
            futures.append(fut)
    wait(futures, return_when=ALL_COMPLETED)

def rcsb_sync():
    """
    Add only the profiles that are missing versus the RCSB.
    """
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
    
    # This logic needs to be firmed up, but here's the idea:
    # 1. Get all structs from RCSB
    # all_rcsb_structs = get_all_tubulin_structs(YOUR_GQL_SEARCH_STRING)

    # 2. Get all structs from DB
    # all_db_structs = adapter.get_all_struct_ids()

    # 3. Find the difference
    # missing_structs = set(all_rcsb_structs) - set(all_db_structs)
    
    # Placeholder:
    missing_structs = [] # Replace with your diff logic
    print(f"Found {len(missing_structs)} missing structures to sync.")

    futures: list[Future] = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in sorted(missing_structs):
            # You'd need an ETL step here first!
            # 1. Run ETLCollector for rcsb_id
            # 2. Then submit the adapter job
            # For now, this assumes profiles are already downloaded
            fut = executor.submit(partial(adapter.add_total_structure, rcsb_id, True))
            futures.append(fut)
    wait(futures, return_when=ALL_COMPLETED)
