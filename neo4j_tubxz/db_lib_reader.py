import typing
import sys
from typing import Optional, List, Tuple, Union, Any

from pydantic import BaseModel, Field
from neo4j import ManagedTransaction, Transaction

# Adjust imports to match your project structure
from etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_tubxz.db_lib_builder import Neo4jAdapter

sys.dont_write_bytecode = True

class StructureFilterParams(BaseModel):

    cursor               : Optional[str]                                     = None
    limit                : int                                               = Field(default=20, ge=1, le=100)
    year                 : Optional[tuple[Optional[int], Optional[int]]]     = None
    search               : Optional[str]                                     = None
    resolution           : Optional[tuple[Optional[float], Optional[float]]] = None
    source_taxa          : Optional[List[int]]                               = None
    host_taxa            : Optional[List[int]]                               = None
    polymerization_state: Optional[List[str]]                                = None

class PolymersFilterParams(BaseModel):

    cursor        : Optional[ Union[Tuple[Optional[str], Optional[str]], List[Optional[str]], str] ] = None
    limit         : int = Field(default=20, ge=1, le=100)
    
    # Structure-level filters
    year          : Optional[Tuple[Optional[int], Optional[int]]] = None
    search        : Optional[str] = None
    resolution    : Optional[Tuple[Optional[float], Optional[float]]] = None
    source_taxa   : Optional[List[int]] = None
    host_taxa     : Optional[List[int]] = None
    
    # Protein/Polymer-level filters
    family        : Optional[List[str]] = None
    uniprot_id    : Optional[str] = None
    has_motif     : Optional[str] = None

    def get_cursor(self) -> Optional[Tuple[Optional[str], Optional[str]]]:
        if self.cursor is None: return None
        if isinstance(self.cursor, str): return (self.cursor, None) # Fallback
        if isinstance(self.cursor, (list, tuple)) and len(self.cursor) == 2:
            return (self.cursor[0], self.cursor[1])
        return None


class Neo4jReader:
    adapter: Neo4jAdapter

    def __init__(self, adapter: Neo4jAdapter | None = None) -> None:
        if adapter:
            self.adapter = adapter
            return
        self.adapter = Neo4jAdapter(
            NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD
        )

    def node_types(self) -> List[str]:
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("CALL db.labels() YIELD label RETURN collect(label)").value()[0]
            return session.execute_read(_)

    def all_ids(self) -> List[str]:
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """MATCH (s:Structure) RETURN collect(s.rcsb_id)"""
                ).value()[0]
            return session.execute_read(_)

    def get_structure_by_id(self, rcsb_id: str) -> dict[str, Any] | None:
        """Fetches a single structure with all components via map projection."""
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
    
    def list_structs_filtered(self, filters: StructureFilterParams):
        query_parts = ["MATCH (s:Structure)"]
        where_clauses = []
        params = {"limit": filters.limit}

        # --- Filters ---
        if filters.search:
            where_clauses.append("(toLower(s.citation_title) CONTAINS $search OR toLower(s.rcsb_id) CONTAINS $search OR toLower(s.pdbx_keywords_text) CONTAINS $search)")
            params["search"] = filters.search.lower()
        
        if filters.year:
            start, end = filters.year
            if start is not None:
                where_clauses.append("s.citation_year >= $year_start")
                params["year_start"] = start
            if end is not None:
                where_clauses.append("s.citation_year <= $year_end")
                params["year_end"] = end
        
        if filters.resolution:
            start, end = filters.resolution
            if start is not None:
                where_clauses.append("s.resolution >= $res_start")
                params["res_start"] = float(start)
            if end is not None:
                where_clauses.append("s.resolution <= $res_end")
                params["res_end"] = float(end)

        # For lineages, we use the relationships created in link_structure_to_phylogeny
        # Direction: (Structure)<-[belongs_to...]-(PhylogenyNode)
        if filters.source_taxa:
            where_clauses.append("EXISTS { MATCH (s)<-[:belongs_to_lineage_source]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $source_taxa }")
            params["source_taxa"] = filters.source_taxa
        
        if filters.host_taxa:
            where_clauses.append("EXISTS { MATCH (s)<-[:belongs_to_lineage_host]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $host_taxa }")
            params["host_taxa"] = filters.host_taxa
            
        if filters.polymerization_state:
            where_clauses.append("s.polymerization_state IN $poly_state")
            params["poly_state"] = filters.polymerization_state

        if where_clauses:
            query_parts.append("WHERE " + " AND ".join(where_clauses))

        # --- Pagination ---
        query_parts.extend([
            "WITH s ORDER BY s.rcsb_id DESC",
            "WITH collect(s) AS all_results",
            "WITH all_results, size(all_results) AS total_count"
        ])

        if filters.cursor:
            query_parts.append("WITH total_count, [s IN all_results WHERE s.rcsb_id < $cursor] AS filtered_results")
            params["cursor"] = filters.cursor
        else:
            query_parts.append("WITH total_count, all_results AS filtered_results")

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
        structure_query = ["MATCH (s:Structure)"]
        struct_where = []
        params = {"limit": filters.limit}

        # 1. Filter Structures first (Optimization)
        if filters.search:
            struct_where.append("(toLower(s.citation_title) CONTAINS $search OR toLower(s.rcsb_id) CONTAINS $search)")
            params["search"] = filters.search.lower()
        
        # Reuse tax logic
        if filters.source_taxa:
            struct_where.append("EXISTS { MATCH (s)<-[:belongs_to_lineage_source]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $source_taxa }")
            params["source_taxa"] = filters.source_taxa

        if struct_where:
            structure_query.append("WHERE " + " AND ".join(struct_where))

        # 2. Match into Proteins
        # Note: TubulinProtein nodes have labels :Polymer AND :Protein
        structure_query.append("MATCH (s)-[:HAS_POLYMER]-(p:Protein)") 
        
        poly_where = []
        if filters.family:
            poly_where.append("p.family IN $family")
            params["family"] = filters.family
        
        if filters.uniprot_id:
            # uniprot_accession is a list of strings in the model
            poly_where.append("$uniprot_id IN p.uniprot_accession")
            params["uniprot_id"] = filters.uniprot_id

        if filters.has_motif:
            poly_where.append("p.entity_poly_seq_one_letter_code_can CONTAINS $motif")
            params["motif"] = filters.has_motif

        if poly_where:
            structure_query.append("WHERE " + " AND ".join(poly_where))

        # 3. Sort and Collect
        structure_query.extend([
            "WITH s, p ORDER BY s.rcsb_id DESC, p.auth_asym_id ASC",
            "WITH collect({struct: s.rcsb_id, protein: properties(p)}) AS all_rows",
            "WITH all_rows, size(all_rows) AS total_count"
        ])

        # 4. Pagination (Complex Cursor: rcsb_id + auth_asym_id)
        cursor = filters.get_cursor()
        if cursor and cursor[0]:
            # Cypher doesn't handle tuple comparison well across mixed types, manual logic:
            # Row < Cursor means: (row.struct < c.struct) OR (row.struct == c.struct AND row.prot.auth < c.auth)
            structure_query.append("""
            WITH total_count, [r IN all_rows WHERE 
                r.struct < $cur_struct OR (r.struct = $cur_struct AND r.protein.auth_asym_id > $cur_auth)
            ] AS filtered_rows
            """)
            params["cur_struct"] = cursor[0]
            params["cur_auth"] = cursor[1] if cursor[1] else ""
        else:
             structure_query.append("WITH total_count, all_rows AS filtered_rows")

        structure_query.extend([
            "WITH total_count, filtered_rows[0..$limit] AS page, filtered_rows",
            "WITH total_count, page, CASE WHEN size(filtered_rows) > $limit THEN [page[-1].struct, page[-1].protein.auth_asym_id] ELSE null END AS next_cursor",
            "RETURN [x IN page | x.protein] AS proteins, total_count, next_cursor"
        ])

        query = "\n".join(structure_query)
        
        with self.adapter.driver.session() as session:
            def _(tx):
                rec = tx.run(query, params).single()
                if not rec: return [], 0, None
                return rec["proteins"], rec["total_count"], rec["next_cursor"]
            return session.execute_read(_)

    def get_taxa(self, src_host: typing.Literal["source", "host"]) -> list[int]:
        # Map input literal to relationship name
        rel_map = {
            "source": "belongs_to_lineage_source",
            "host": "belongs_to_lineage_host"
        }
        rel_name = rel_map.get(src_host, "belongs_to_lineage_source")

        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                # Direction: Structure <-[rel]- PhylogenyNode
                query = f"""
                MATCH (s:Structure)<-[:{rel_name}]-(p:PhylogenyNode)
                RETURN collect(DISTINCT p.ncbi_tax_id)
                """
                return tx.run(query).single().value()
            return session.execute_read(_)

dbqueries = Neo4jReader()