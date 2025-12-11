import sys, json
from neo4j import GraphDatabase, Driver, ManagedTransaction
from typing import Tuple, List, Dict
import datetime

from lib.etl.assets import GlobalOps, TubulinStructureAssets
from lib.etl.libtax import PhylogenyNode, Taxid
from lib.models.types_tubulin import (
    TubulinStructure,
    NonpolymericLigand,
    MasterAlignment,
    TubulinFamily,
    AlignmentMapping,
    Mutation,
    Modification,
    TubulinEntity,
    PolymerInstance,
)
from neo4j_tubxz.node_ligand import link__ligand_to_struct, node__ligand
from neo4j_tubxz.node_master_alignment import (
    get_master_alignment,
    link__polymer_to_master_alignment,
    node__master_alignment,
)
from neo4j_tubxz.node_phylogeny import link__phylogeny, node__phylogeny

# TUBE-UPDATE: Split Polymer node logic
from neo4j_tubxz.node_polymer import (
    node__entity,
    upsert__entity_to_tubulin,
    node__instance,
    link__instance_to_structure,
    link__instance_to_entity,
)

from neo4j_tubxz.node_structure import (
    link__structure_to_lineage_member,
    link__structure_to_organism,
    node__structure,
    struct_exists,
)

NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT rcsb_id_unique IF NOT EXISTS FOR (s:Structure) REQUIRE s.rcsb_id IS UNIQUE;""",
    """CREATE CONSTRAINT entity_unique IF NOT EXISTS FOR (e:Entity) REQUIRE (e.parent_rcsb_id, e.entity_id) IS NODE KEY;""",
    """CREATE CONSTRAINT instance_unique IF NOT EXISTS FOR (i:Instance) REQUIRE (i.parent_rcsb_id, i.auth_asym_id) IS NODE KEY;""",
    """CREATE CONSTRAINT taxid_unique IF NOT EXISTS FOR (phylonode:PhylogenyNode) REQUIRE phylonode.ncbi_tax_id IS UNIQUE;""",
    """CREATE CONSTRAINT chemicalId IF NOT EXISTS FOR (ligand:Ligand) REQUIRE ligand.chemicalId IS UNIQUE;""",
    """CREATE CONSTRAINT ma_unique IF NOT EXISTS FOR (a:MasterAlignment) REQUIRE (a.family, a.version) IS NODE KEY;""",
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
        self, rcsb_id: str, profile: TubulinStructure, verbose: bool = False
    ):
        rcsb_id = rcsb_id.upper()
        if profile is None:
            profile: TubulinStructure = TubulinStructureAssets(rcsb_id).profile()

        with self.driver.session() as s:
            for organism_host in profile.host_organism_ids:
                self._create_lineage(organism_host)
                s.execute_write(
                    link__structure_to_organism(rcsb_id, organism_host, "host")
                )
                for org in Taxid.get_lineage(organism_host):
                    s.execute_write(
                        link__structure_to_lineage_member(
                            rcsb_id, org, "belongs_to_lineage_host"
                        )
                    )
            for organism_src in profile.src_organism_ids:
                self._create_lineage(organism_src)
                s.execute_write(
                    link__structure_to_organism(rcsb_id, organism_src, "source")
                )
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

    def upsert_ligand_node(self, ligand: NonpolymericLigand, parent_rcsb_id: str):
        with self.driver.session() as s:
            ligand_node = s.execute_write(node__ligand(ligand, parent_rcsb_id))
            s.execute_write(link__ligand_to_struct(ligand_node, parent_rcsb_id))

    # --- Main ETL Function ---

    def add_total_structure(self, rcsb_id: str, disable_exists_check: bool = False):
        rcsb_id = rcsb_id.upper()

        print("got ", rcsb_id)
        if not disable_exists_check:
            if self.check_structure_exists(rcsb_id):
                print(f"\nStruct node {rcsb_id} already exists.")
                return

        # 1. EXTRACT
        S: TubulinStructure = TubulinStructureAssets(rcsb_id).profile()

        with self.driver.session() as s:
            # 2. LOAD (A): Create Structure node
            structure_node = s.execute_write(node__structure(S))

            # 3. LINK: Link Structure to Phylogeny
            self.link_structure_to_phylogeny(rcsb_id, S)

            # 4. LOAD (B): Process Ligands
            if S.nonpolymeric_ligands is not None:
                for ligand in S.nonpolymeric_ligands:
                    self.upsert_ligand_node(ligand, S.rcsb_id)

            # 5. LOAD (C): Process Entities (Other Polymers)
            if S.other_entities is not None:
                for poly_entity in S.other_entities:
                    e_node = s.execute_write(node__entity(poly_entity))

            # 6. LOAD (D): Process Tubulin Entities
            if S.tubulin_entities is not None:
                for tub_entity in S.tubulin_entities:
                    # Create basic entity
                    e_node = s.execute_write(node__entity(tub_entity))
                    # Add Tubulin properties (Family, Uniprot, etc)
                    s.execute_write(upsert__entity_to_tubulin(e_node, tub_entity))

            # 7. LOAD (E): Process Instances (Chains) and Link
            if S.instances is not None:
                for instance in S.instances:
                    # Create Instance Node
                    i_node = s.execute_write(node__instance(instance))

                    # Link Instance -> Structure
                    s.execute_write(link__instance_to_structure(i_node, S.rcsb_id))

                    # Link Instance -> Entity (Parent)
                    s.execute_write(
                        link__instance_to_entity(i_node, instance.entity_id, S.rcsb_id)
                    )

            # TODO: Add logic here if you want to attach Mutation nodes to Specific Instances
            # or Alignment Mappings to Specific Entities.

        print(f"Successfully initialized structure {rcsb_id}.")
        return structure_node

    def delete_structure(self, rcsb_id: str, dry_run: bool = False) -> dict[str, int]:
        rcsb_id = rcsb_id.upper()

        query = """
        MATCH (s:Structure {rcsb_id: $rcsb_id})
        
        OPTIONAL MATCH (s)-[r1]-(i:Instance)
        OPTIONAL MATCH (i)-[r_ent]-(e:Entity) // We usually don't delete Entities as they might be shared, but in this specific schema, Entities are 'parent_rcsb_id' specific so we DO delete them.
        
        OPTIONAL MATCH (s)-[r2]-(l:Ligand)
        OPTIONAL MATCH (s)-[r3]-(t:PhylogenyNode)
        
        // Collect nodes to delete
        WITH s,
             collect(DISTINCT i) as instances,
             collect(DISTINCT e) as entities,
             collect(DISTINCT l) as ligands, // Careful, ligands are shared usually, relationships deleted
             collect(DISTINCT r1) as instRels,
             collect(DISTINCT r_ent) as entRels,
             collect(DISTINCT r2) as ligandRels,
             collect(DISTINCT r3) as taxaRels
        
        // Delete relationships
        FOREACH (r IN instRels | DELETE r)
        FOREACH (r IN entRels | DELETE r)
        FOREACH (r IN ligandRels | DELETE r)
        FOREACH (r IN taxaRels | DELETE r)
        
        // Delete nodes specific to this structure
        FOREACH (i IN instances | DETACH DELETE i)
        FOREACH (e IN entities | DETACH DELETE e) // Deleting entities because they are scoped by parent_rcsb_id in Pydantic model
        
        DETACH DELETE s
        
        RETURN count(s) as deleted
        """

        with self.driver.session() as session:
            if not self.check_structure_exists(rcsb_id):
                raise ValueError(f"Structure {rcsb_id} does not exist in database")
            try:
                # Note: Simplified delete for brevity.
                # In a shared-entity world, we would check if (e) has other connections.
                # Since your Entity model includes `parent_rcsb_id`, it is scoped to the structure, so safe to delete.
                result = session.execute_write(
                    lambda tx: tx.run(query, rcsb_id=rcsb_id).single()
                )
                print(f"Deleted structure {rcsb_id}")
                return {"deleted": 1}
            except Exception as e:
                print(f"Error deleting structure {rcsb_id}: {str(e)}")
                raise
