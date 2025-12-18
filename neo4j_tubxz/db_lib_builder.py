# neo4j_tubxz/db_lib_builder.py
import sys, json
from neo4j import GraphDatabase, Driver, ManagedTransaction
from typing import Tuple, List, Dict
import datetime

from lib.etl.assets import GlobalOps, TubulinStructureAssets
from lib.etl.libtax import PhylogenyNode, Taxid

from lib.types import (
    TubulinStructure,
    PolypeptideEntity,
    PolynucleotideEntity,
    NonpolymerEntity,
    Polypeptide,
    Polynucleotide,
    Nonpolymer,
    TubulinFamily,
    Mutation,
)

# New Node Logic Imports
from neo4j_tubxz.node_ligand import (
    node__chemical,
    node__nonpolymer_entity,
    node__nonpolymer_instance,
)

from neo4j_tubxz.node_polymer import (
    node__polypeptide_entity,
    node__polynucleotide_entity,
    node__polymer_instance,
    link__entity_to_structure,
    link__instance_to_structure,
)

from neo4j_tubxz.node_master_alignment import (
    get_master_alignment,
    node__master_alignment,
)
from neo4j_tubxz.node_phylogeny import link__phylogeny, node__phylogeny
from neo4j_tubxz.node_mutation import node__mutation, link__mutation_to_master_alignment

from neo4j_tubxz.node_structure import (
    link__structure_to_lineage_member,
    link__structure_to_organism,
    node__structure,
    struct_exists,
)

NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT rcsb_id_unique IF NOT EXISTS FOR (s:Structure) REQUIRE s.rcsb_id IS UNIQUE;""",
    
    """CREATE CONSTRAINT entity_unique IF NOT EXISTS FOR (e:Entity) REQUIRE (e.parent_rcsb_id, e.entity_id) IS NODE KEY;""",
    
    # --- TUBE-FIX: Change constraint from auth_asym_id to asym_id ---
    """CREATE CONSTRAINT instance_unique IF NOT EXISTS FOR (i:Instance) REQUIRE (i.parent_rcsb_id, i.asym_id) IS NODE KEY;""",
    # -----------------------------------------------------------------
    
    """CREATE CONSTRAINT chemical_unique IF NOT EXISTS FOR (c:Chemical) REQUIRE c.chemical_id IS UNIQUE;""",
    """CREATE CONSTRAINT taxid_unique IF NOT EXISTS FOR (phylonode:PhylogenyNode) REQUIRE phylonode.ncbi_tax_id IS UNIQUE;""",
    """CREATE CONSTRAINT ma_unique IF NOT EXISTS FOR (a:MasterAlignment) REQUIRE (a.family, a.version) IS NODE KEY;""",
    """CREATE CONSTRAINT mutation_unique IF NOT EXISTS FOR (m:Mutation) REQUIRE (m.master_index, m.uniprot_id, m.from_residue, m.to_residue) IS NODE KEY;""",
]


class Neo4jAdapter:
    driver: Driver
    uri: str
    user: str
    databases: list[str]

    def __init__(
        self,
        uri: str,
        user: str,
        current_db: str,
        password: str | None = None,
    ) -> None:
        self.uri = uri
        self.user = user
        self.password = password
        try:
            self.driver = GraphDatabase.driver(
                uri, auth=(user, password), database=current_db
            )
        except Exception as ae:
            print(ae)
            raise ae

    def initialize_new_instance(self):
        self.init_constraints()
        self.init_phylogenies()

    def init_constraints(self) -> None:
        with self.driver.session() as session:
            for c in NODE_CONSTRAINTS:
                session.execute_write(lambda tx: tx.run(c))
                print("\nAdded constraint: ", c)

    def add_phylogeny_node(self, taxid: int):
        with self.driver.session() as session:
            node = session.execute_write(
                node__phylogeny(PhylogenyNode.from_taxid(taxid))
            )
            return node

    def init_phylogenies(self):
        taxa = GlobalOps.collect_all_taxa()
        for taxon in taxa:
            self._create_lineage(taxon.ncbi_tax_id)

    def _create_lineage(self, taxid: int) -> None:
        lin = Taxid.get_lineage(taxid)
        lin.reverse()
        previous_id: int | None = None
        with self.driver.session() as session:
            for taxid in lin:
                node = session.execute_write(
                    node__phylogeny(PhylogenyNode.from_taxid(taxid))
                )
                if previous_id is None:
                    previous_id = taxid
                    continue
                session.execute_write(link__phylogeny(taxid, previous_id))
                previous_id = taxid
        return

    def link_structure_to_phylogeny(
            self, rcsb_id: str, profile: TubulinStructure | None = None, verbose: bool = False
        ):
            rcsb_id = rcsb_id.upper()
            
            # Avoid variable shadowing
            _profile = profile
            if _profile is None:
                _profile = TubulinStructureAssets(rcsb_id).profile()

            # --- SIMPLIFIED LOGIC (Using Restored Root Fields) ---
            # Since we restored src_organism_ids/host_organism_ids to the TubulinStructure root,
            # we can access them directly. This bypasses the entity loop and the Pylance errors.
            
            host_ids = _profile.host_organism_ids
            src_ids = _profile.src_organism_ids
            
            with self.driver.session() as s:
                # Handle host organisms
                for organism_host in host_ids:
                    self._create_lineage(organism_host)
                    s.execute_write(
                        link__structure_to_organism(rcsb_id, organism_host, "host")
                    )
                    # Link full lineage
                    for org in Taxid.get_lineage(organism_host):
                        s.execute_write(
                            link__structure_to_lineage_member(
                                rcsb_id, org, "belongs_to_lineage_host"
                            )
                        )
                
                # Handle source organisms
                for organism_src in src_ids:
                    self._create_lineage(organism_src)
                    s.execute_write(
                        link__structure_to_organism(rcsb_id, organism_src, "source")
                    )
                    # Link full lineage
                    for org in Taxid.get_lineage(organism_src):
                        s.execute_write(
                            link__structure_to_lineage_member(
                                rcsb_id, org, "belongs_to_lineage_source"
                            )
                        )

    def check_structure_exists(self, rcsb_id: str) -> bool:
        rcsb_id = rcsb_id.upper()
        with self.driver.session() as session:
            return session.execute_read(struct_exists(rcsb_id))

    def process_entity_mutations(
            self, 
            entity: PolypeptideEntity, 
            entity_node_element_id: str,
            parent_rcsb_id: str
        ):
            """
            Process mutations stored on the Entity.
            Links: Entity -> HAS_MUTATION -> Mutation -> ANNOTATES -> MasterAlignment
            """
            if not entity.mutations or not entity.family:
                return
                
            with self.driver.session() as s:
                # TUBE-FIX: Pass family and version strings, not a Pydantic model
                ma_version = "v1.0"
                ma_node = s.execute_write(
                    node__master_alignment(entity.family, ma_version)
                )
                
                # Create mutation nodes and link them
                for mut in entity.mutations:
                    mut_node = s.execute_write(node__mutation(mut))
                    
                    # Link Mutation to MasterAlignment (Definition)
                    s.execute_write(
                        link__mutation_to_master_alignment(mut_node, ma_node)
                    )
                    
                    # Link Entity to Mutation (Occurrence)
                    s.execute_write(
                        lambda tx: tx.run("""
                            MATCH (m:Mutation) WHERE ELEMENTID(m) = $mut_id
                            MATCH (e:Entity) WHERE ELEMENTID(e) = $ent_id
                            MERGE (e)-[:HAS_MUTATION]->(m)
                            """, 
                            {
                                "mut_id": mut_node.element_id,
                                "ent_id": entity_node_element_id
                            }
                        )
                    )

    # --- Main ETL Function ---

    def add_total_structure(self, rcsb_id: str, disable_exists_check: bool = False):
        rcsb_id = rcsb_id.upper()

        print(f"Processing {rcsb_id}...")

        if not disable_exists_check:
            if self.check_structure_exists(rcsb_id):
                print(f"Structure {rcsb_id} already exists.")
                return

        # 1. EXTRACT - Load the profile
        S: TubulinStructure = TubulinStructureAssets(rcsb_id).profile()

        with self.driver.session() as s:
            # 2. CREATE STRUCTURE NODE
            print(f"  Creating Structure node...")
            structure_node = s.execute_write(node__structure(S))

            # 3. LINK TO PHYLOGENY
            print(f"  Linking to phylogeny...")
            self.link_structure_to_phylogeny(rcsb_id, S)

            # 4. PROCESS ENTITIES (The Blueprints)
            # We iterate the entities dictionary from the profile
            print(f"  Processing {len(S.entities)} entities...")

            # Keep track of created entity IDs to ensure successful creation
            created_entities = set()

            for entity_id, entity in S.entities.items():
                if isinstance(entity, NonpolymerEntity):
                    # 4a. Ligand: Ensure Global Chemical -> Create Local Entity
                    s.execute_write(node__chemical(entity))
                    e_node = s.execute_write(node__nonpolymer_entity(entity, S.rcsb_id))
                    s.execute_write(link__entity_to_structure(e_node, S.rcsb_id))
                    created_entities.add(entity_id)

                elif isinstance(entity, PolypeptideEntity):
                    # 4b. Polypeptide: Create Entity -> Process Mutations
                    e_node = s.execute_write(
                        node__polypeptide_entity(entity, S.rcsb_id)
                    )
                    s.execute_write(link__entity_to_structure(e_node, S.rcsb_id))

                    if entity.mutations:
                        self.process_entity_mutations(
                            entity, e_node.element_id, S.rcsb_id
                        )
                    created_entities.add(entity_id)

                elif isinstance(entity, PolynucleotideEntity):
                    # 4c. Polynucleotide: Create Entity
                    e_node = s.execute_write(
                        node__polynucleotide_entity(entity, S.rcsb_id)
                    )
                    s.execute_write(link__entity_to_structure(e_node, S.rcsb_id))
                    created_entities.add(entity_id)

            # 5. PROCESS INSTANCES (The Physical Copies)
            # These link to the Entities created above

            print(f"  Processing {len(S.polypeptides)} polypeptide instances...")
            for instance in S.polypeptides:
                if instance.entity_id in created_entities:
                    i_node = s.execute_write(node__polymer_instance(instance))
                    s.execute_write(link__instance_to_structure(i_node, S.rcsb_id))

            print(f"  Processing {len(S.polynucleotides)} polynucleotide instances...")
            for instance in S.polynucleotides:
                if instance.entity_id in created_entities:
                    i_node = s.execute_write(node__polymer_instance(instance))
                    s.execute_write(link__instance_to_structure(i_node, S.rcsb_id))

            print(f"  Processing {len(S.nonpolymers)} nonpolymer instances...")
            for instance in S.nonpolymers:
                if instance.entity_id in created_entities:
                    i_node = s.execute_write(node__nonpolymer_instance(instance))
                    s.execute_write(link__instance_to_structure(i_node, S.rcsb_id))

        print(f"✓ Successfully initialized structure {rcsb_id}\n")
        return structure_node

    def delete_structure(self, rcsb_id: str, dry_run: bool = False) -> dict[str, int]:
        rcsb_id = rcsb_id.upper()

        # TUBE-UPDATE: Updated delete logic for new schema
        query = """
        MATCH (s:Structure {rcsb_id: $rcsb_id})
        
        // 1. Get Instances linked to Structure
        OPTIONAL MATCH (s)-[r_inst_s]->(i:Instance)
        
        // 2. Get Entities linked to Structure
        OPTIONAL MATCH (s)-[r_ent_s]->(e:Entity)
        
        // 3. Get relationships between Instance and Entity
        OPTIONAL MATCH (i)-[r_inst_ent]->(e)
        
        // 4. Get Mutations linked to Entities (don't delete mutation nodes themselves if shared, 
        //    but here mutations are unique to the entity's context usually, or shared via merge.
        //    In this schema, mutations are unique nodes but linked to specific entities.
        //    We detach the entity.)
        
        // 5. Get relationships to Chemicals (don't delete Chemical node)
        OPTIONAL MATCH (e)-[r_def_chem]->(c:Chemical)

        WITH s, 
             collect(DISTINCT i) as instances,
             collect(DISTINCT e) as entities,
             collect(DISTINCT r_inst_s) as rels_inst_struct,
             collect(DISTINCT r_ent_s) as rels_ent_struct,
             collect(DISTINCT r_inst_ent) as rels_inst_ent,
             collect(DISTINCT r_def_chem) as rels_def_chem

        // Delete relationships
        FOREACH (r IN rels_inst_struct | DELETE r)
        FOREACH (r IN rels_ent_struct | DELETE r)
        FOREACH (r IN rels_inst_ent | DELETE r)
        FOREACH (r IN rels_def_chem | DELETE r)
        
        // Delete nodes (Entities and Instances are scoped to this structure, so safe to delete)
        FOREACH (i IN instances | DETACH DELETE i)
        FOREACH (e IN entities | DETACH DELETE e)
        
        // Finally delete structure
        DETACH DELETE s
        
        RETURN 1 as deleted
        """

        with self.driver.session() as session:
            if not self.check_structure_exists(rcsb_id):
                raise ValueError(f"Structure {rcsb_id} does not exist in database")
            try:
                result = session.execute_write(
                    lambda tx: tx.run(query, rcsb_id=rcsb_id).single()
                )
                print(f"Deleted structure {rcsb_id}")
                return {"deleted": result["deleted"]}
            except Exception as e:
                print(f"Error deleting structure {rcsb_id}: {str(e)}")
                raise
