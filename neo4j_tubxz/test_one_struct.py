import asyncio
import os

from etl.collector import TubulinETLCollector
from etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_tubxz.db_lib_builder import Neo4jAdapter
import sys
sys.path.append("/Users/rtviii/dev/tubulinxyz")
TEST_RCSB_ID = "6O2T"

async def run_test():
    print(f"--- 1. Generating Profile for {TEST_RCSB_ID} ---")
    collector = TubulinETLCollector(TEST_RCSB_ID)
    try:
        await collector.generate_profile(overwrite=False)
        print(f"Profile for {TEST_RCSB_ID} is ready.")
    except Exception as e:
        print(f"CRITICAL: Failed to generate profile: {e}")
        return

    print(f"\n--- 2. Initializing Neo4j Adapter ---")
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
    try:
        print(f"Checking for and deleting any old version of {TEST_RCSB_ID}...")
        adapter.delete_structure(TEST_RCSB_ID)
    except ValueError as e:
        print(f"(Structure did not exist, which is fine)")
    except Exception as e:
        print(f"Error during deletion: {e}")

    print(f"\n--- 3. Calling add_total_structure for {TEST_RCSB_ID} ---")
    
    adapter.add_total_structure(TEST_RCSB_ID, disable_exists_check=True)
    print(f"\n--- 4. SUCCESS ---")
    print(f"Successfully added structure {TEST_RCSB_ID} to the database.")

if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    
    asyncio.run(run_test())