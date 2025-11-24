import sys
import typing
from typing import Optional, List, Tuple, Union, Any, Dict
from pydantic import BaseModel, Field
from neo4j import ManagedTransaction, Transaction

# Adjust imports to match your project structure
from neo4j_tubxz.db_lib_builder import Neo4jAdapter
from lib.etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from lib.models.types_tubulin import TubulinFamily

sys.dont_write_bytecode = True

# --- Pydantic Models for Filters ---

class StructureFilterParams(BaseModel):
    cursor              : Optional[str]                                         = None
    limit               : int                                                   = Field(default=20, ge=1, le=20000)
    year                : Optional[tuple[Optional[int], Optional[int]]]         = None
    search              : Optional[str]                                         = None
    resolution          : Optional[tuple[Optional[float], Optional[float]]]     = None
    source_taxa         : Optional[List[int]]                                   = None
    host_taxa           : Optional[List[int]]                                   = None
    polymerization_state: Optional[List[str]]                                   = None

class PolymersFilterParams(BaseModel):
    cursor           : Optional[Union[Tuple[Optional[str], Optional[str]], List[Optional[str]], str]] = None
    limit            : int = Field(default=20, ge=1, le=20000)
    
    # Structure-level filters
    year             : Optional[Tuple[Optional[int], Optional[int]]] = None
    search           : Optional[str] = None
    resolution       : Optional[Tuple[Optional[float], Optional[float]]] = None
    source_taxa      : Optional[List[int]] = None
    host_taxa        : Optional[List[int]] = None
    
    # Protein/Polymer-level filters
    family           : Optional[List[str]] = None
    uniprot_id       : Optional[str] = None
    has_motif        : Optional[str] = None

    def get_cursor(self) -> Optional[Tuple[Optional[str], Optional[str]]]:
        if self.cursor is None: return None
        if isinstance(self.cursor, str): return (self.cursor, None)
        if isinstance(self.cursor, (list, tuple)) and len(self.cursor) >= 2:
            return (self.cursor[0], self.cursor[1])
        return None

# --- Main Reader Class ---

class Neo4jReader:
    adapter: Neo4jAdapter

    def __init__(self, adapter: Neo4jAdapter | None = None) -> None:
        if adapter:
            self.adapter = adapter
            return
        self.adapter = Neo4jAdapter(
            NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD
        )

    # --- Basic Lookups ---

    def all_ids(self) -> List[str]:
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """MATCH (s:Structure) RETURN collect(s.rcsb_id)"""
                ).value()[0]
            return session.execute_read(_)

    def get_structure_by_id(self, rcsb_id: str) -> dict[str, Any] | None:
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = """
                MATCH (s:Structure {rcsb_id: $rcsb_id})
                RETURN s {
                    .*, 
                    ligands: [(s)-[:CONTAINS]-(l:Ligand) | properties(l)],
                    proteins: [(s)-[:HAS_POLYMER]-(p:Protein) | properties(p)],
                    other_polymers: [(s)-[:HAS_POLYMER]-(p:Polymer) WHERE NOT p:Protein | properties(p)]
                } as structure
                """
                result = tx.run(query, rcsb_id=rcsb_id.upper()).single()
                return result["structure"] if result else None
            return session.execute_read(_)

    def random_structure(self) -> dict[str, Any]:
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = """
                MATCH (s:Structure) 
                WITH s, rand() AS r 
                ORDER BY r LIMIT 1
                RETURN s {
                    .*, 
                    ligands: [(s)-[:CONTAINS]-(l:Ligand) | properties(l)],
                    proteins: [(s)-[:HAS_POLYMER]-(p:Protein) | properties(p)],
                    other_polymers: [(s)-[:HAS_POLYMER]-(p:Polymer) WHERE NOT p:Protein | properties(p)]
                } as structure
                """
                return tx.run(query).single()["structure"]
            return session.execute_read(_)

    def get_taxa(self, src_host: typing.Literal["source", "host"]) -> list[int]:
        rel_map = {
            "source": "belongs_to_lineage_source",
            "host": "belongs_to_lineage_host"
        }
        rel_name = rel_map.get(src_host, "belongs_to_lineage_source")
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = f"""
                MATCH (s:Structure)<-[:{rel_name}]-(p:PhylogenyNode)
                RETURN collect(DISTINCT p.ncbi_tax_id)
                """
                return tx.run(query).single().value()
            return session.execute_read(_)

    def get_tax_dict(self):
        """Returns [[id, name], ...] for the frontend dropdowns"""
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """MATCH (t:PhylogenyNode) RETURN collect([t.ncbi_tax_id, t.scientific_name])"""
                ).value()[0]
            return session.execute_read(_)

    # --- Filtering Logic (Ported) ---

    def list_structs_filtered(self, filters: StructureFilterParams):
        query_parts = ["MATCH (s:Structure)"]
        where_clauses = []
        params = {"limit": filters.limit}

        # 1. Search (Title, ID, Keywords)
        if filters.search:
            where_clauses.append(
                """
                (toLower(s.citation_title) CONTAINS $search OR 
                 toLower(s.rcsb_id) CONTAINS $search OR 
                 toLower(s.pdbx_keywords_text) CONTAINS $search OR
                 toLower(reduce(acc = '', str IN s.src_organism_names | acc + str)) CONTAINS $search)
                """
            )
            params["search"] = filters.search.lower()
        
        # 2. Year (Integer comparison)
        if filters.year:
            start, end = filters.year
            if start is not None:
                where_clauses.append("s.citation_year >= $year_start")
                params["year_start"] = start
            if end is not None:
                where_clauses.append("s.citation_year <= $year_end")
                params["year_end"] = end
        
        # 3. Resolution
        if filters.resolution:
            start, end = filters.resolution
            if start is not None:
                where_clauses.append("s.resolution >= $res_start")
                params["res_start"] = float(start)
            if end is not None:
                where_clauses.append("s.resolution <= $res_end")
                params["res_end"] = float(end)

        # 4. Taxa (Source/Host)
        if filters.source_taxa:
            where_clauses.append("EXISTS { MATCH (s)<-[:belongs_to_lineage_source]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $source_taxa }")
            params["source_taxa"] = filters.source_taxa
        
        if filters.host_taxa:
            where_clauses.append("EXISTS { MATCH (s)<-[:belongs_to_lineage_host]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $host_taxa }")
            params["host_taxa"] = filters.host_taxa
            
        # 5. Polymerization State (Tubulin specific)
        if filters.polymerization_state:
            where_clauses.append("s.polymerization_state IN $poly_state")
            params["poly_state"] = filters.polymerization_state

        # Apply Where
        if where_clauses:
            query_parts.append("WHERE " + " AND ".join(where_clauses))

        # Sort and Collect
        query_parts.extend([
            "WITH s ORDER BY s.rcsb_id DESC",
            "WITH collect(s) AS all_results",
            "WITH all_results, size(all_results) AS total_count"
        ])

        # Pagination Cursor
        if filters.cursor:
            query_parts.append("WITH total_count, [s IN all_results WHERE s.rcsb_id < $cursor] AS filtered_results")
            params["cursor"] = filters.cursor
        else:
            query_parts.append("WITH total_count, all_results AS filtered_results")

        # Page Slicing
        query_parts.extend([
            "WITH total_count, filtered_results[0..$limit] AS page, filtered_results",
            "WITH total_count, page, CASE WHEN size(filtered_results) > $limit THEN page[-1].rcsb_id ELSE null END AS next_cursor",
            "RETURN [x IN page | properties(x)] as structures, total_count, next_cursor"
        ])

        query = "\n".join(query_parts)

        with self.adapter.driver.session() as session:
            def _(tx):
                rec = tx.run(query, params).single()
                if not rec: return [], 0, None
                return rec["structures"], rec["total_count"], rec["next_cursor"]
            return session.execute_read(_)

    def list_polymers_filtered(self, filters: PolymersFilterParams):
        # 1. Filter Structures First
        query_parts = ["MATCH (s:Structure)"]
        struct_where = []
        params = {"limit": filters.limit}

        if filters.search:
            struct_where.append("(toLower(s.citation_title) CONTAINS $search OR toLower(s.rcsb_id) CONTAINS $search)")
            params["search"] = filters.search.lower()
        
        if filters.year:
            start, end = filters.year
            if start: 
                struct_where.append("s.citation_year >= $year_start")
                params["year_start"] = start
            if end:   
                struct_where.append("s.citation_year <= $year_end")
                params["year_end"] = end

        if filters.source_taxa:
            struct_where.append("EXISTS { MATCH (s)<-[:belongs_to_lineage_source]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $source_taxa }")
            params["source_taxa"] = filters.source_taxa

        if struct_where:
            query_parts.append("WHERE " + " AND ".join(struct_where))

        # 2. Collect filtered structures
        query_parts.extend([
            "WITH collect(s) as filtered_structs",
            "WITH filtered_structs, size(filtered_structs) as total_struct_count"
        ])

        # 3. Unwind and Match Proteins
        query_parts.extend([
            "UNWIND filtered_structs as s",
            "MATCH (s)-[:HAS_POLYMER]->(p:Protein)" # Ensure we only get proteins
        ])

        poly_where = []
        
        if filters.family:
            poly_where.append("p.family IN $family")
            params["family"] = filters.family
        
        if filters.uniprot_id:
            poly_where.append("$uniprot_id IN p.uniprot_accession")
            params["uniprot_id"] = filters.uniprot_id

        if filters.has_motif:
            poly_where.append("p.entity_poly_seq_one_letter_code_can CONTAINS $motif")
            params["motif"] = filters.has_motif

        if poly_where:
            query_parts.append("WHERE " + " AND ".join(poly_where))

        # 4. Sorting and grouping
        query_parts.extend([
            "WITH total_struct_count, s, p",
            "ORDER BY s.rcsb_id DESC, p.auth_asym_id ASC",
            "WITH total_struct_count, collect({struct: s.rcsb_id, protein: properties(p)}) AS all_rows",
            "WITH total_struct_count, all_rows, size(all_rows) AS total_poly_count"
        ])

        # 5. Cursor Logic (StructID + AuthAsymID)
        cursor = filters.get_cursor()
        if cursor and cursor[0]:
            query_parts.append("""
            WITH total_struct_count, total_poly_count, [r IN all_rows WHERE 
                r.struct < $cur_struct OR (r.struct = $cur_struct AND r.protein.auth_asym_id > $cur_auth)
            ] AS filtered_rows
            """)
            params["cur_struct"] = cursor[0]
            params["cur_auth"] = cursor[1] if cursor[1] else ""
        else:
             query_parts.append("WITH total_struct_count, total_poly_count, all_rows AS filtered_rows")

        # 6. Page Slicing
        query_parts.extend([
            "WITH total_struct_count, total_poly_count, filtered_rows[0..$limit] AS page, filtered_rows",
            "WITH total_struct_count, total_poly_count, page, CASE WHEN size(filtered_rows) > $limit THEN [page[-1].struct, page[-1].protein.auth_asym_id] ELSE null END AS next_cursor",
            "RETURN [x IN page | x.protein] AS proteins, total_struct_count, total_poly_count, next_cursor"
        ])

        query = "\n".join(query_parts)
        
        print("\n\033[96m" + query + "\033[0m\n")

        with self.adapter.driver.session() as session:
            def _(tx):
                rec = tx.run(query, params).single()
                if not rec: return [], 0, 0, None
                return rec["proteins"], rec["total_struct_count"], rec["total_poly_count"], rec["next_cursor"]
            return session.execute_read(_)
    # neo4j_tubxz/db_lib_reader.py (add these methods to Neo4jReader class)

    def get_mutations_at_position(
        self, 
        master_index: int, 
        family: str = "alpha",
        version: str = "v1.0"
    ) -> List[Dict[str, Any]]:
        """Get all mutations at a specific master alignment position"""
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = """
                MATCH (ma:MasterAlignment {family: $family, version: $version})
                MATCH (m:Mutation {master_index: $master_index})-[:ANNOTATES_POSITION_IN]->(ma)
                RETURN collect(properties(m)) as mutations
                """
                result = tx.run(query, {
                    "family": family,
                    "version": version,
                    "master_index": master_index
                }).single()
                
                return result["mutations"] if result else []
            return session.execute_read(_)

    def get_modifications_at_position(
        self,
        master_index: int,
        family: str = "alpha", 
        version: str = "v1.0"
    ) -> List[Dict[str, Any]]:
        """Get all modifications at a specific master alignment position"""
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = """
                MATCH (ma:MasterAlignment {family: $family, version: $version})
                MATCH (m:Modification {master_index: $master_index})-[:ANNOTATES_POSITION_IN]->(ma)
                RETURN collect(properties(m)) as modifications
                """
                result = tx.run(query, {
                    "family": family,
                    "version": version,
                    "master_index": master_index
                }).single()
                
                return result["modifications"] if result else []
            return session.execute_read(_)

    def get_all_annotations_at_position(
        self,
        master_index: int,
        family: str = "alpha",
        version: str = "v1.0"
    ) -> Dict[str, Any]:
        """Get both mutations and modifications at a specific position"""
        return {
            "position": master_index,
            "family": family,
            "version": version,
            "mutations": self.get_mutations_at_position(master_index, family, version),
            "modifications": self.get_modifications_at_position(master_index, family, version)
        }

    # neo4j_tubxz/db_lib_reader.py (add to Neo4jReader class)

    def get_mutations_for_polymer(
        self,
        rcsb_id: str,
        auth_asym_id: str
    ) -> List[Dict[str, Any]]:
        """Get all mutations for a specific polymer"""
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = """
                MATCH (p:Polymer {parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id})
                OPTIONAL MATCH (p)-[:HAS_MUTATION]->(m:Mutation)
                RETURN collect(properties(m)) as mutations
                """
                result = tx.run(query, {
                    "rcsb_id": rcsb_id.upper(),
                    "auth_asym_id": auth_asym_id
                }).single()
                
                # Filter out null entries
                muts = result["mutations"] if result else []
                return [m for m in muts if m is not None]
            return session.execute_read(_)

    def get_modifications_for_polymer(
        self,
        rcsb_id: str,
        auth_asym_id: str
    ) -> List[Dict[str, Any]]:
        """Get all modifications for a specific polymer"""
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                query = """
                MATCH (p:Polymer {parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id})
                OPTIONAL MATCH (p)-[:HAS_MODIFICATION]->(m:Modification)
                RETURN collect(properties(m)) as modifications
                """
                result = tx.run(query, {
                    "rcsb_id": rcsb_id.upper(),
                    "auth_asym_id": auth_asym_id
                }).single()
                
                # Filter out null entries from OPTIONAL MATCH
                mods = result["modifications"] if result else []
                return [m for m in mods if m is not None]
            return session.execute_read(_)

    def get_all_annotations_for_polymer(
        self,
        rcsb_id: str,
        auth_asym_id: str
    ) -> Dict[str, Any]:
        """Get both mutations and modifications for a specific polymer"""
        mutations = self.get_mutations_for_polymer(rcsb_id, auth_asym_id)
        modifications = self.get_modifications_for_polymer(rcsb_id, auth_asym_id)
        
        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "mutations": mutations,
            "modifications": modifications,
            "total_count": len(mutations) + len(modifications)
        }

dbqueries = Neo4jReader()