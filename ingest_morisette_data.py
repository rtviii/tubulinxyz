import json
import os
from pathlib import Path
from datetime import datetime
from api.config import PROJECT_ROOT
from lib.etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_tubxz.db_driver import Neo4jAdapter
from neo4j_tubxz.node_mutation import node__mutation, link__mutation_to_master_alignment
from neo4j_tubxz.node_modification import node__modification, link__modification_to_master_alignment
from neo4j_tubxz.node_master_alignment import node__master_alignment, get_master_alignment
from lib.models.types_tubulin import Modification, Mutation, MasterAlignment, TubulinFamily

def create_master_alignment(adapter: Neo4jAdapter, fasta_path: str, family: TubulinFamily, version: str):
    """Create or get the master alignment node"""

    with adapter.driver.session() as session:
        ma_node = session.execute_read(get_master_alignment(family, version))
        if ma_node:
            print(f"Master alignment {family.value} {version} already exists")
            return ma_node

    with open(fasta_path, 'r') as f:
        fasta_content = f.read()
    

    ma = MasterAlignment(
        family        = family,
        version       = version,
        fasta_content = fasta_content,
        created_date  = datetime.now().isoformat(),
        description   = f"Master alignment for {family.value} tubulin"
    )
    
    # Create node
    with adapter.driver.session() as session:
        ma_node = session.execute_write(node__master_alignment(ma))
    
    print(f"Created master alignment {family.value} {version}")
    return ma_node

def ingest_mutations(adapter: Neo4jAdapter, json_path: str, ma_node):
    """Load mutations from JSON and create nodes linked to MA"""
    with open(json_path, 'r') as f:
        mutations_data = json.load(f)
    
    print(f"Loading {len(mutations_data)} mutations...")
    created = 0
    skipped = 0
    
    with adapter.driver.session() as session:
        for i, mut_data in enumerate(mutations_data, 1):
            if i % 100 == 0:
                print(f"  Processed {i}/{len(mutations_data)} mutations (created: {created}, skipped: {skipped})...")
            
            # Skip if no MA mapping
            if mut_data['ma_status'] != 'mapped':
                skipped += 1
                continue
            
            mut = Mutation(
                master_index    = mut_data['ma_position'],
                utn_position    = mut_data.get('utn_position'),
                from_residue    = mut_data['wild_type'],
                to_residue      = mut_data['mutant'],
                uniprot_id      = mut_data['tubulin_id'],
                species         = mut_data['species'],
                tubulin_type    = mut_data['tubulin_type'],
                phenotype       = mut_data['phenotype'],
                database_source = mut_data.get('reference', 'Unknown'),
                reference_link  = mut_data.get('ref_link', ''),
                keywords        = mut_data.get('keywords', ''),
                notes           = mut_data.get('notes', '')
            )
            
            # Create mutation node
            mut_node = session.execute_write(node__mutation(mut))
            
            # Link to master alignment
            session.execute_write(link__mutation_to_master_alignment(mut_node, ma_node))
            
            created += 1
    
    print(f"Successfully loaded {created} mutations (skipped {skipped})")

def ingest_modifications(adapter: Neo4jAdapter, json_path: str, ma_node):
    """Load modifications from JSON and create nodes linked to MA"""
    with open(json_path, 'r') as f:
        mods_data = json.load(f)
    
    print(f"Loading {len(mods_data)} modifications...")
    
    created = 0
    skipped = 0
    
    with adapter.driver.session() as session:
        for i, mod_data in enumerate(mods_data, 1):
            if i % 100 == 0:
                print(f"  Processed {i}/{len(mods_data)} modifications (created: {created}, skipped: {skipped})...")
            
            # Skip if no MA mapping
            if mod_data['ma_status'] != 'mapped':
                skipped += 1
                continue
            
            mod = Modification(
                master_index      = mod_data['ma_position'],
                utn_position      = mod_data.get('utn_position'),
                amino_acid        = mod_data['amino_acid'],
                modification_type = mod_data['modification_type'],
                uniprot_id        = mod_data['tubulin_id'],
                species           = mod_data['species'],
                tubulin_type      = mod_data['tubulin_type'],
                phenotype         = mod_data['phenotype'],
                database_source   = mod_data.get('database', 'Unknown'),
                database_link     = mod_data.get('db_link', ''),
                keywords          = mod_data.get('keywords', ''),
                notes             = mod_data.get('notes', '')
            )
            
            mod_node = session.execute_write(node__modification(mod))
            session.execute_write(link__modification_to_master_alignment(mod_node, ma_node))
            
            created += 1
    
    print(f"Successfully loaded {created} modifications (skipped {skipped})")

def main():
    # Configuration
    MASTER_ALIGNMENT_FASTA = os.path.join(PROJECT_ROOT,"data/alpha_tubulin/alpha_tubulin.afasta" )
    MUTATIONS_JSON     = os.path.join(PROJECT_ROOT,'data/alpha_tubulin/alpha_mutations.json')
    MODIFICATIONS_JSON = os.path.join(PROJECT_ROOT,'data/alpha_tubulin/alpha_modifications.json')
    
    adapter = Neo4jAdapter(
        uri        = NEO4J_URI,
        user       = NEO4J_USER,
        current_db = NEO4J_CURRENTDB,
        password   = NEO4J_PASSWORD
    )
    
    
    try:
        ma_node = create_master_alignment(
            adapter,
            MASTER_ALIGNMENT_FASTA,
            TubulinFamily.ALPHA,
            "v1.0"
        )
        
        if Path(MUTATIONS_JSON).exists():
            ingest_mutations(adapter, MUTATIONS_JSON, ma_node)
        else:
            print(f"Warning: {MUTATIONS_JSON} not found")
        
        if Path(MODIFICATIONS_JSON).exists():
            ingest_modifications(adapter, MODIFICATIONS_JSON, ma_node)
        else:
            print(f"Warning: {MODIFICATIONS_JSON} not found")
    
    finally:
        adapter.driver.close()
    
    print("\nDone!")

if __name__ == "__main__":
    main()