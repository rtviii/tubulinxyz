# neo4j_tubxz/db_lib_builder.py
"""
Neo4j database builder - handles structure ingestion.
"""

import sys
from neo4j import GraphDatabase, Driver, ManagedTransaction, Transaction
from typing import List

from lib.etl.assets import GlobalOps, TubulinStructureAssets
from lib.etl.libtax import PhylogenyNode, Taxid

from lib.types import (
    TubulinStructure,
    PolypeptideEntity,
    PolynucleotideEntity,
    NonpolymerEntity,
)

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
from neo4j_tubxz.node_phylogeny import link__phylogeny, node__phylogeny
from neo4j_tubxz.node_structure import (
    link__structure_to_lineage_member,
    link__structure_to_organism,
    node__structure,
    struct_exists,
)
from neo4j_tubxz.node_variant import create_variants_for_entity
from neo4j_tubxz.node_binding_site import process_all_binding_sites

sys.dont_write_bytecode = True


NODE_CONSTRAINTS = [
    # Structure
    "CREATE CONSTRAINT rcsb_id_unique IF NOT EXISTS FOR (s:Structure) REQUIRE s.rcsb_id IS UNIQUE;",
    # Entity (unique per structure)
    "CREATE CONSTRAINT entity_unique IF NOT EXISTS FOR (e:Entity) REQUIRE (e.parent_rcsb_id, e.entity_id) IS NODE KEY;",
    # Instance (unique per structure by asym_id)
    "CREATE CONSTRAINT instance_unique IF NOT EXISTS FOR (i:Instance) REQUIRE (i.parent_rcsb_id, i.asym_id) IS NODE KEY;",
    # Chemical (global)
    "CREATE CONSTRAINT chemical_unique IF NOT EXISTS FOR (c:Chemical) REQUIRE c.chemical_id IS UNIQUE;",
    # Phylogeny
    "CREATE CONSTRAINT taxid_unique IF NOT EXISTS FOR (p:PhylogenyNode) REQUIRE p.ncbi_tax_id IS UNIQUE;",
    # Variant indexes for querying
    "CREATE INDEX variant_master_idx IF NOT EXISTS FOR (v:Variant) ON (v.master_index);",
    "CREATE INDEX variant_type_idx IF NOT EXISTS FOR (v:Variant) ON (v.type);",
    "CREATE INDEX variant_source_idx IF NOT EXISTS FOR (v:Variant) ON (v.source);",
    # Substitution-specific index for common queries
    "CREATE INDEX substitution_positions IF NOT EXISTS FOR (v:Substitution) ON (v.master_index, v.wild_type, v.observed);",
]


class Neo4jAdapter:
    driver: Driver
    uri: str
    user: str

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
        except Exception as e:
            print(f"Failed to connect to Neo4j: {e}")
            raise

    def close(self):
        self.driver.close()

    # =========================================================================
    # Initialization
    # =========================================================================

    def initialize_new_instance(self):
        """Set up constraints and seed phylogeny tree."""
        self.init_constraints()
        self.init_phylogenies()

    def init_constraints(self) -> None:
        with self.driver.session() as session:
            for constraint in NODE_CONSTRAINTS:
                try:
                    session.execute_write(lambda tx, c=constraint: tx.run(c))
                    print(f"Applied: {constraint[:60]}...")
                except Exception as e:
                    print(f"Skipped (may exist): {constraint[:40]}... ({e})")

    def init_phylogenies(self):
        """Seed phylogeny tree from all profiles on disk."""
        taxa = GlobalOps.collect_all_taxa()
        print(f"Seeding {len(taxa)} taxa into phylogeny tree...")
        for taxon in taxa:
            self._create_lineage(taxon.ncbi_tax_id)

    def _create_lineage(self, taxid: int) -> None:
        """Create full lineage path for a taxid."""
        lineage = Taxid.get_lineage(taxid)
        lineage.reverse()  # Root to leaf order

        previous_id: int | None = None
        with self.driver.session() as session:
            for tid in lineage:
                session.execute_write(node__phylogeny(PhylogenyNode.from_taxid(tid)))
                if previous_id is not None:
                    session.execute_write(link__phylogeny(tid, previous_id))
                previous_id = tid

    # =========================================================================
    # Structure Existence Check
    # =========================================================================

    def check_structure_exists(self, rcsb_id: str) -> bool:
        rcsb_id = rcsb_id.upper()
        with self.driver.session() as session:
            return session.execute_read(struct_exists(rcsb_id))

    # =========================================================================
    # Main Ingestion (backwards compatible name)
    # =========================================================================

    def add_total_structure(self, rcsb_id: str, disable_exists_check: bool = False):
        """
        Ingest a structure into the database.

        Args:
            rcsb_id: PDB ID
            disable_exists_check: If True, delete existing structure first (upsert behavior)
        """
        rcsb_id = rcsb_id.upper()
        print(f"\n{'=' * 60}")
        print(f"Ingesting {rcsb_id}")
        print(f"{'=' * 60}")

        # Check existence
        if self.check_structure_exists(rcsb_id):
            if disable_exists_check:
                print(f"  Deleting existing structure...")
                self.delete_structure(rcsb_id)
            else:
                print(f"  Structure already exists, skipping.")
                return

        # Load profile
        try:
            profile: TubulinStructure = TubulinStructureAssets(rcsb_id).profile()
        except FileNotFoundError:
            print(f"  ERROR: Profile not found for {rcsb_id}")
            return

        with self.driver.session() as s:
            # 1. Create Structure node
            print(f"  Creating Structure node...")
            s.execute_write(node__structure(profile))

            # 2. Link to phylogeny
            print(f"  Linking to phylogeny...")
            self._link_structure_to_phylogeny(s, rcsb_id, profile)

            # 3. Process entities
            print(f"  Processing {len(profile.entities)} entities...")
            entity_ids_created = set()

            for entity_id, entity in profile.entities.items():
                if isinstance(entity, NonpolymerEntity):
                    self._process_nonpolymer_entity(s, entity, rcsb_id)
                elif isinstance(entity, PolypeptideEntity):
                    self._process_polypeptide_entity(s, entity, rcsb_id)
                elif isinstance(entity, PolynucleotideEntity):
                    self._process_polynucleotide_entity(s, entity, rcsb_id)

                entity_ids_created.add(entity_id)

            # 4. Process instances
            print(f"  Processing instances...")
            self._process_instances(s, profile, entity_ids_created)

            # 5. Process binding sites
            if profile.ligand_binding_sites:
                print(
                    f"  Processing {len(profile.ligand_binding_sites)} binding sites..."
                )
                # process_all_binding_sites returns a transaction function
                tx_func = process_all_binding_sites(
                    rcsb_id, profile.ligand_binding_sites
                )

                # Execute the function within a write transaction
                count = s.execute_write(tx_func)
                print(f"    Created {count} NEAR_POLYMER relationships")
            else:
                print("  No ligand binding sites to ingest.")

        print(f"  Done: {rcsb_id}")

    def _link_structure_to_phylogeny(
        self, session, rcsb_id: str, profile: TubulinStructure
    ):
        """Link structure to source and host organisms."""
        # Source organisms
        for org_id in profile.src_organism_ids:
            self._create_lineage(org_id)
            session.execute_write(
                link__structure_to_organism(rcsb_id, org_id, "source")
            )
            for lineage_id in Taxid.get_lineage(org_id):
                session.execute_write(
                    link__structure_to_lineage_member(
                        rcsb_id, lineage_id, "belongs_to_lineage_source"
                    )
                )

        # Host organisms
        for org_id in profile.host_organism_ids:
            self._create_lineage(org_id)
            session.execute_write(link__structure_to_organism(rcsb_id, org_id, "host"))
            for lineage_id in Taxid.get_lineage(org_id):
                session.execute_write(
                    link__structure_to_lineage_member(
                        rcsb_id, lineage_id, "belongs_to_lineage_host"
                    )
                )

    def _process_nonpolymer_entity(
        self, session, entity: NonpolymerEntity, rcsb_id: str
    ):
        """Process a nonpolymer (ligand) entity."""
        # Ensure global Chemical node exists
        session.execute_write(node__chemical(entity))
        # Create entity node
        e_node = session.execute_write(node__nonpolymer_entity(entity, rcsb_id))
        # Link to structure
        session.execute_write(link__entity_to_structure(e_node, rcsb_id))

    def _process_polypeptide_entity(
        self, session, entity: PolypeptideEntity, rcsb_id: str
    ):
        """Process a polypeptide entity with variants."""
        # Create entity node (includes index mappings)
        e_node = session.execute_write(node__polypeptide_entity(entity, rcsb_id))
        # Link to structure
        session.execute_write(link__entity_to_structure(e_node, rcsb_id))

        # Create variant nodes
        if entity.variants:
            count = session.execute_write(
                lambda tx: create_variants_for_entity(
                    tx, entity.variants, rcsb_id, entity.entity_id
                )
            )
            if count > 0:
                print(f"    Entity {entity.entity_id}: {count} variants")

    def _process_polynucleotide_entity(
        self, session, entity: PolynucleotideEntity, rcsb_id: str
    ):
        """Process a polynucleotide entity."""
        e_node = session.execute_write(node__polynucleotide_entity(entity, rcsb_id))
        session.execute_write(link__entity_to_structure(e_node, rcsb_id))

    def _process_instances(self, session, profile: TubulinStructure, entity_ids: set):
        """Process all physical instances."""
        # Get ALL entity IDs that exist in the profile
        all_entity_ids = set(profile.entities.keys())
        
        # Create counters for debugging
        poly_count = 0
        nonpoly_count = 0
        
        # Polypeptides
        for instance in profile.polypeptides:
            if instance.entity_id in all_entity_ids:
                i_node = session.execute_write(node__polymer_instance(instance))
                session.execute_write(
                    link__instance_to_structure(i_node, profile.rcsb_id)
                )
                poly_count += 1

        # Polynucleotides  
        for instance in profile.polynucleotides:
            if instance.entity_id in all_entity_ids:
                i_node = session.execute_write(node__polymer_instance(instance))
                session.execute_write(
                    link__instance_to_structure(i_node, profile.rcsb_id)
                )
                poly_count += 1

        # Nonpolymers
        for instance in profile.nonpolymers:
            if instance.entity_id in all_entity_ids:
                i_node = session.execute_write(node__nonpolymer_instance(instance))
                session.execute_write(
                    link__instance_to_structure(i_node, profile.rcsb_id)
                )
                nonpoly_count += 1
            else:
                print(f"  Warning: Nonpolymer instance {instance.asym_id} references entity {instance.entity_id} not found in entities dict")

        # Print debug info
        if poly_count > 0 or nonpoly_count > 0:
            print(f"    Created {poly_count} polymer instances, {nonpoly_count} nonpolymer instances")

    # =========================================================================
    # Deletion
    # =========================================================================

    def delete_structure(self, rcsb_id: str) -> None:
        """
        Delete a structure and all its entities/instances/variants.
        Does NOT delete shared nodes (Chemical, PhylogenyNode).
        """
        rcsb_id = rcsb_id.upper()

        query = """
        MATCH (s:Structure {rcsb_id: $rcsb_id})

        // Get all entities
        OPTIONAL MATCH (s)-[:DEFINES_ENTITY]->(e:Entity)

        // Get all instances
        OPTIONAL MATCH (s)-[:HAS_INSTANCE]->(i:Instance)

        // Get all variants linked to polypeptide entities
        OPTIONAL MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)

        // Detach and delete variants
        DETACH DELETE v

        // Detach and delete instances
        DETACH DELETE i

        // Detach and delete entities
        DETACH DELETE e

        // Finally delete structure
        DETACH DELETE s
        """

        with self.driver.session() as session:
            session.execute_write(lambda tx: tx.run(query, {"rcsb_id": rcsb_id}))
            print(f"Deleted structure {rcsb_id}")

    # =========================================================================
    # Bulk Operations
    # =========================================================================

    def add_all_structures(self, overwrite: bool = False) -> None:
        """Ingest all structures from disk profiles."""
        profiles = GlobalOps.list_profiles()
        print(f"Found {len(profiles)} profiles to ingest")

        for rcsb_id in sorted(profiles):
            try:
                self.add_total_structure(rcsb_id, disable_exists_check=overwrite)
            except Exception as e:
                print(f"  ERROR processing {rcsb_id}: {e}")
                continue

    def get_all_structure_ids(self) -> List[str]:
        """Get all structure IDs in the database."""
        with self.driver.session() as session:
            result = session.run("MATCH (s:Structure) RETURN collect(s.rcsb_id) AS ids")
            return result.single()["ids"]
