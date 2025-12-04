
import asyncio
import os
import sys
import time

from etl.collector import TubulinETLCollector
from etl.constants import TUBETL_DATA
from etl.assets import GlobalOps

# --- Helper to make sure we can find the project modules ---
def add_project_root_to_path():
    project_root = os.path.abspath(os.path.dirname(__file__))
    if project_root not in sys.path:
        sys.path.append(project_root)
        print(f"Added {project_root} to sys.path")

add_project_root_to_path()

SEMAPHORE = asyncio.Semaphore(10)

async def download_profile(rcsb_id: str):
    """
    Async worker to download and save a single structure profile.
    """
    async with SEMAPHORE:
        collector = TubulinETLCollector(rcsb_id)
        try:
            # Use overwrite=False to skip profiles that already exist
            await collector.generate_profile(overwrite=False)
            print(f"Successfully processed profile for {rcsb_id}")
            return rcsb_id, "Success"
        except Exception as e:
            print(f"!!! FAILED to process profile for {rcsb_id}: {e}")
            return rcsb_id, f"Failed: {e}"

async def main():
    """
    Fetches all tubulin structure IDs and runs the ETL collector for each one.
    """
    print("--- Step 1: Fetching all tubulin structure IDs from RCSB ---")
    
    # Ensure the base data directory exists
    os.makedirs(TUBETL_DATA, exist_ok=True)
    
    try:
        all_rcsb_ids = GlobalOps.current_rcsb_structs()
        print(f"Found {len(all_rcsb_ids)} total tubulin structures in RCSB.")
    except Exception as e:
        print(f"CRITICAL: Failed to fetch structure list from RCSB. {e}")
        return

    print(f"\n--- Step 2: Downloading all {len(all_rcsb_ids)} profiles ---")
    print(f"(This will take a while. Profiles will be saved in '{TUBETL_DATA}')\n")
    
    start_time = time.time()
    
    # Create a list of async tasks
    tasks = [download_profile(rcsb_id) for rcsb_id in all_rcsb_ids]
    
    # Run tasks concurrently
    results = await asyncio.gather(*tasks)
    
    end_time = time.time()
    
    success_count = len([res for res in results if res[1] == "Success"])
    fail_count = len(results) - success_count
    
    print("\n--- 3. COMPLETE ---")
    print(f"Processing finished in {end_time - start_time:.2f} seconds.")
    print(f"Successfully downloaded/verified: {success_count}")
    print(f"Failed: {fail_count}")
    print("\nYou can now run 'p3 init_database.py' to seed the database.")

if __name__ == "__main__":
    asyncio.run(main())