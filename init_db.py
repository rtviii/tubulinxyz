import os
import sys
from neo4j_tubxz.db_driver import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_tubxz.db_lib_builder import Neo4jAdapter

def add_project_root_to_path():
    project_root = os.path.abspath(os.path.dirname(__file__))
    if project_root not in sys.path:
        sys.path.append(project_root)
        print(f"Added {project_root} to sys.path")

if __name__ == "__main__":
    add_project_root_to_path()
    
    print("--- Initializing Neo4j Database ---")
    
    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
        adapter.initialize_new_instance()
        
        print("\n--- 4. SUCCESS ---")
        print("Database initialized successfully.")
        print("Constraints, master alignments, and phylogeny nodes are now in place.")
        print("You can now run 'test_one_structure.py'.")

    except Exception as e:
        print(f"An error occurred during database initialization: {e}")
        import traceback
        traceback.print_exc()