import sys, json
from neo4j import GraphDatabase, Driver, ManagedTransaction
from typing import Tuple, List, Dict
import datetime

# TUBE-UPDATE: Import all new node/link functions and schema
from etl.assets import GlobalOps
from etl.libtax import PhylogenyNode, Taxid
from models.types_tubulin import (
    TubulinStructure, NonpolymericLigand, MasterAlignment, TubulinFamily,
    AlignmentMapping, Mutation, Modification, TubulinProtein
)
from neo4j_tubxz.node_ligand import link__ligand_to_struct, node__ligand
from neo4j_tubxz.node_master_alignment import get_master_alignment, link__polymer_to_master_alignment, node__master_alignment
from neo4j_tubxz.node_modification import link__polymer_to_modification, node__modification
from neo4j_tubxz.node_mutation import link__polymer_to_mutation, node__mutation
from neo4j_tubxz.node_phylogeny import link__phylogeny, node__phylogeny
from neo4j_tubxz.node_polymer import link__polymer_to_structure, node__polymer, upsert__polymer_to_protein
from neo4j_tubxz.node_structure import link__structure_to_lineage_member, link__structure_to_organism, node__structure, struct_exists
from etl.assets import GlobalOps, TubulinStructureAssets


NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT rcsb_id_unique IF NOT EXISTS FOR (s:Structure) REQUIRE s.rcsb_id IS UNIQUE;""",
    """CREATE CONSTRAINT polymer_unique IF NOT EXISTS FOR (p:Polymer) REQUIRE (p.parent_rcsb_id, p.auth_asym_id) IS NODE KEY;""",
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
        # self.init_master_alignments()
        self.init_phylogenies()

    def init_constraints(self) -> None:
        with self.driver.session() as session:
            for c in NODE_CONSTRAINTS:
                session.execute_write(lambda tx: tx.run(c))
                print("\nAdded constraint: ", c)

    # def init_master_alignments(self):
    #     """
    #     Load canonical master alignments into the DB.
    #     """
    #     print("Initializing Master Alignments...")
    #     # Placeholder: In a real app, load this from files.
    #     today = datetime.date.today().isoformat()
    #     master_seqs = [
    #         # MasterAlignment(family=TubulinFamily.ALPHA, version="v1.0", fasta_content="MREVI...", created_date=today, description="Canonical Alpha Tubulin v1.0"),
    #         # MasterAlignment(family=TubulinFamily.BETA, version="v1.0", fasta_content="MREIV...", created_date=today, description="Canonical Beta Tubulin v1.0"),
    #         # MasterAlignment(family=TubulinFamily.GAMMA, version="v1.0", fasta_content="MRECI...", created_date=today, description="Canonical Gamma Tubulin v1.0"),
    #     ]
        
    #     with self.driver.session() as session:
    #         for aln in master_seqs:
    #             session.execute_write(node__master_alignment(aln))
    #     print(f"Initialized {len(master_seqs)} master alignments.")

    # ... (Phylogeny functions: add_phylogeny_node, init_phylogenies, _create_lineage - unchanged) ...
    def add_phylogeny_node(self, taxid: int): #-> Node:
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
                session.execute_write(
                    link__phylogeny(taxid, previous_id)
                )
                previous_id = taxid
        return
    
    def link_structure_to_phylogeny(
        self, rcsb_id: str, profile: TubulinStructure, verbose: bool = False
    ):
        rcsb_id = rcsb_id.upper()
        if profile is None:
            # TUBE-FIX: Use TubulinStructureAssets and call .profile()
            profile:TubulinStructure = TubulinStructureAssets(rcsb_id).profile()
        
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

    def upsert_ligand_node(
        self, ligand: NonpolymericLigand, parent_rcsb_id: str
    ):
        """Uglysert a ligand node and link it to its parent structure."""
        with self.driver.session() as s:
            # Pass parent_rcsb_id to node__ligand for the ON CREATE/ON MATCH logic
            ligand_node = s.execute_write(node__ligand(ligand, parent_rcsb_id))
            s.execute_write(link__ligand_to_struct(ligand_node, parent_rcsb_id))

    # --- ETL Placeholder Functions ---
    # These functions would be defined in your ETL/ops logic
    
    def _get_auth_to_seqres_map(self, protein: TubulinProtein) -> Dict[int, int]:
        """
        Placeholder: Gets map from auth_seq_id -> seqres_index.
        (From _pdbx_poly_seq_scheme)
        """
        # print(f"Stub: Getting auth_to_seqres map for {protein.auth_asym_id}")
        # This is just a dummy 1-to-1 mapping for illustration
        return {i+1: i for i in range(protein.entity_poly_seq_length)}

    def _run_muscle_profile(self, seqres_seq: str, master_fasta: str) -> str:
        """
        Placeholder: Runs MUSCLE alignment.
        Returns raw alignment text.
        """
        # print(f"Stub: Running MUSCLE for seq of length {len(seqres_seq)}")
        # Dummy alignment result
        return f">master\n{master_fasta}\n>seqres\n{seqres_seq}\n"

    def _parse_alignment(self, alignment_text: str, master_aln_node, parent_rcsb_id: str, auth_asym_id: str) -> Tuple[AlignmentMapping, List[Mutation]]:
        """
        Placeholder: Parses MUSCLE output.
        Returns mapping object and list of Mutation objects.
        """
        # print(f"Stub: Parsing alignment for {parent_rcsb_id}.{auth_asym_id}")
        # Dummy data
        map_list = [-1] * 200
        mapping = AlignmentMapping(
            seqres_to_master=json.dumps(map_list),
            master_to_seqres=json.dumps(map_list)
        )
        mutations = []
        return (mapping, mutations)

    def _find_ptms(self, protein: TubulinProtein) -> List[Dict]:
        """
        Placeholder: Finds PTMs for a polymer *in auth coordinates*.
        Returns a list of dictionaries.
        """
        # print(f"Stub: Finding PTMs for {protein.auth_asym_id}")
        # Dummy PTM
        # return [
        #     {"type": ModificationType.ACETYLATION, "auth_seq_id": 40, "evidence": "...", "pubmed": "..."}
        # ]
        return []

    # --- Main ETL Function ---

    def add_total_structure(self, rcsb_id: str, disable_exists_check: bool = False):
        rcsb_id = rcsb_id.upper()

        print("got ", rcsb_id)
        if not disable_exists_check:
            if self.check_structure_exists(rcsb_id):
                print(f"\nStruct node {rcsb_id} already exists.")
                return

        # 1. EXTRACT: Get base structure profile
        # TUBE-FIX: Use TubulinStructureAssets and call .profile()
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
            
            # 5. LOAD (C): Process "Other" Polymers
            if S.other_polymers is not None:
                for poly in S.other_polymers:
                    p_node = s.execute_write(node__polymer(poly))
                    s.execute_write(link__polymer_to_structure(p_node, S.rcsb_id))

            # 6. ETL Process for Tubulin Proteins
            for protein in S.proteins:
                
                # 6a. LOAD: Create base Polymer node
                p_node = s.execute_write(node__polymer(protein))
                
                # 6b. LOAD: Upsert to Protein, adding :Protein label and tubulin data
                s.execute_write(upsert__polymer_to_protein(p_node, protein))

                # 6c. LINK: Link Polymer to Structure
                s.execute_write(link__polymer_to_structure(p_node, S.rcsb_id))

                if protein.family is None:
                    print(f"WARNING: Skipping alignment for {rcsb_id}.{protein.auth_asym_id} (no family).")
                    continue

                # # 6d. GET MASTER: Fetch the master alignment node
                # master_aln_node = s.execute_read(
                #     get_master_alignment(protein.family, "v1.0") # Use "v1.0" as default
                # )

                # if master_aln_node is None:
                #     print(f"WARNING: No MasterAlignment for family {protein.family}. Skipping mapping for {rcsb_id}.{protein.auth_asym_id}")
                #     continue

                # # 6e. TRANSFORM (A): Get "Rosetta Stone"
                # auth_to_seqres_map = self._get_auth_to_seqres_map(protein)

                # # 6f. TRANSFORM (B): Run MUSCLE
                # alignment_text = self._run_muscle_profile(
                #     protein.entity_poly_seq_one_letter_code_can,
                #     master_aln_node["fasta_content"]
                # )

                # # 6g. TRANSFORM (C): Parse alignment for mappings & mutations
                # (mappings, mutations) = self._parse_alignment(
                #     alignment_text,
                #     master_aln_node,
                #     protein.parent_rcsb_id,
                #     protein.auth_asym_id
                # )

                # # 6h. LOAD (D): Store the mapping on the relationship
                # s.execute_write(link__polymer_to_master_alignment(
                #     p_node, master_aln_node, mappings
                # ))
                
                # # 6i. LOAD (E): Store all mutations
                # for mut in mutations:
                #     m_node = s.execute_write(node__mutation(mut))
                #     s.execute_write(link__polymer_to_mutation(p_node, m_node))

                # # 6j. TRANSFORM (D): Find PTMs
                # ptms = self._find_ptms(protein) # Gets PTMs in auth_seq_id coordinates
                
                # # Load seqres_to_master map
                # seqres_to_master = json.loads(mappings.seqres_to_master)

                # # 6k. LOAD (F): Store all PTMs
                # for ptm in ptms:
                #     seqres_idx = auth_to_seqres_map.get(ptm["auth_seq_id"])
                #     if seqres_idx is None:
                #         # print(f"Warning: Cannot map PTM at auth_seq_id {ptm['auth_seq_id']}")
                #         continue
                    
                #     master_idx = seqres_to_master[seqres_idx]
                #     if master_idx == -1:
                #         # print(f"Warning: PTM at auth_seq_id {ptm['auth_seq_id']} maps to a gap.")
                #         continue

                #     # TODO: Get master_residue from master_fasta at master_idx
                #     master_residue = "X" 

                #     mod_obj = Modification(
                #         modification_type = ptm["type"],
                #         master_index = master_idx,
                #         master_residue = master_residue,
                #         evidence_url = ptm.get("evidence", ""),
                #         pubmed_ids = ptm.get("pubmed", [])
                #     )
                    
                #     mod_node = s.execute_write(node__modification(mod_obj))
                #     s.execute_write(link__polymer_to_modification(p_node, mod_node))

        print(f"Successfully initialized structure {rcsb_id}.")
        return structure_node

    def upsert_structure_node(self, rcsb_id: str):
        rcsb_id = rcsb_id.upper()
        # TUBE-FIX: Use TubulinStructureAssets and call .profile()
        S: TubulinStructure = TubulinStructureAssets(rcsb_id).profile()
        with self.driver.session() as s:
            structure_node = s.execute_write(node__structure(S))
        print(f"Successfully merged structure {rcsb_id}.")
        return structure_node

    def delete_structure(self, rcsb_id: str, dry_run: bool = False) -> dict[str, int]:
        rcsb_id = rcsb_id.upper()
        
        query = """
        MATCH (s:Structure {rcsb_id: $rcsb_id})
        
        OPTIONAL MATCH (s)-[r1]-(p:Polymer)
        OPTIONAL MATCH (s)-[r2]-(l:Ligand)
        OPTIONAL MATCH (s)-[r3]-(t:PhylogenyNode)
        
        OPTIONAL MATCH (p)-[]-(mut:Mutation)
        OPTIONAL MATCH (p)-[]-(mod:Modification)

        WITH s,
             collect(DISTINCT p) as polymers,
             collect(DISTINCT mut) as mutations,
             collect(DISTINCT mod) as modifications,
             collect(DISTINCT r1) as polymerRels,
             collect(DISTINCT r2) as ligandRels,
             collect(DISTINCT r3) as taxaRels,
             count(DISTINCT p) as polymerCount,
             count(DISTINCT mut) as mutationCount,
             count(DISTINCT mod) as modificationCount,
             count(DISTINCT r1) + count(DISTINCT r2) + count(DISTINCT r3) as relCount
        
        FOREACH (r IN ligandRels | DELETE r)
        FOREACH (r IN taxaRels | DELETE r)
        
        // Detach and delete all polymers and their dependent nodes
        FOREACH (p IN polymers | DETACH DELETE p)
        
        // Ensure mutations/modifications are deleted (DETACH DELETE on polymer should cover this)
        // This is redundant if they are ONLY attached to polymers, but safe.
        FOREACH (m IN mutations | DETACH DELETE m)
        FOREACH (md IN modifications | DETACH DELETE md)
        
        DETACH DELETE s
        
        RETURN {
            structure: 1, polymers: polymerCount,
            mutations: mutationCount, modifications: modificationCount,
            relationships: relCount
        } as counts
        """
        
        dry_run_query = """
        MATCH (s:Structure {rcsb_id: $rcsb_id})
        OPTIONAL MATCH (s)-[r1]-(p:Polymer)
        OPTIONAL MATCH (s)-[r2]-(l:Ligand)
        OPTIONAL MATCH (s)-[r3]-(t:PhylogenyNode)
        OPTIONAL MATCH (p)-[]-(mut:Mutation)
        OPTIONAL MATCH (p)-[]-(mod:Modification)
        RETURN {
            structure: 1, polymers: count(DISTINCT p),
            mutations: count(DISTINCT mut), modifications: count(DISTINCT mod),
            relationships: count(DISTINCT r1) + count(DISTINCT r2) + count(DISTINCT r3)
        } as counts
        """
        
        with self.driver.session() as session:
            if not self.check_structure_exists(rcsb_id):
                raise ValueError(f"Structure {rcsb_id} does not exist in database")
            try:
                result = session.execute_write(
                    lambda tx: tx.run(
                        dry_run_query if dry_run else query,
                        rcsb_id=rcsb_id
                    ).single()
                )
                counts = result['counts'] if result else {}
                action = "Would delete" if dry_run else "Deleted"
                print(f"\n{action} for {rcsb_id}: {counts}")
                return counts
            except Exception as e:
                print(f"Error deleting structure {rcsb_id}: {str(e)}")
                raise