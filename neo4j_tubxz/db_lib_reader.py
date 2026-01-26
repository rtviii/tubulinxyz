# neo4j_tubxz/db_lib_reader.py
"""
Database reader with typed query methods.
Uses the query builders for filter logic.
"""

import json
import sys
from typing import Optional, List, Dict, Any
from neo4j import ManagedTransaction, Transaction
# Add this import at the top of db_lib_builder.py
from neo4j_tubxz.node_binding_site import process_all_binding_sites
from neo4j_tubxz.db_lib_builder import Neo4jAdapter
from neo4j_tubxz.models import (
    StructureFilters,
    PolypeptideEntityFilters,
    LigandFilters,
    StructureListResponse,
    PolypeptideListResponse,
    LigandListResponse,
    StructureSummary,
    PolypeptideEntitySummary,
    LigandSummary,
    FilterFacets,
    FacetValue,
    LigandFacet,
    RangeValue,
    VariantsByFamily,
    CommonVariant,
    VariantPositionRange,
    LigandNeighborhood,
    BindingSiteResidue,
    PolymerNeighborhoodsResponse,
)
from neo4j_tubxz.structure_query_builder import (
    StructureQueryBuilder,
    PolypeptideEntityQueryBuilder,
    LigandQueryBuilder,
)
from lib.etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER

sys.dont_write_bytecode = True


class Neo4jReader:
    """
    Read-only database operations with typed filters and responses.
    """

    def __init__(self, adapter: Optional[Neo4jAdapter] = None) -> None:
        self.adapter = adapter or Neo4jAdapter(
            NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD
        )

    # -------------------------------------------------------------------------
    # Structure Queries
    # -------------------------------------------------------------------------

    def list_structures(self, filters: StructureFilters) -> StructureListResponse:
        """
        List structures with filtering and keyset pagination.
        """
        builder = StructureQueryBuilder(filters)
        query, params = builder.build()

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                result = tx.run(query, params)
                records = list(result)

                if not records:
                    return StructureListResponse(
                        data=[], total_count=0, next_cursor=None, has_more=False
                    )

                total_count = records[0]["total_count"]
                next_cursor = records[0]["next_cursor"]

                structures = [
                    StructureSummary(
                        rcsb_id=r["rcsb_id"],
                        resolution=r["resolution"],
                        exp_method=r["exp_method"],
                        citation_title=r["citation_title"],
                        citation_year=r["citation_year"],
                        deposition_date=r["deposition_date"],
                        src_organism_names=r["src_organism_names"] or [],
                        pdbx_keywords=r["pdbx_keywords"],
                        entity_count=r["entity_count"],
                        ligand_count=r["ligand_count"],
                    )
                    for r in records
                ]

                return StructureListResponse(
                    data=structures,
                    total_count=total_count,
                    next_cursor=next_cursor,
                    has_more=next_cursor is not None,
                )

            return session.execute_read(run_query)

    def get_structure(self, rcsb_id: str) -> Optional[Dict[str, Any]]:
        """Get full structure details by ID"""
        query = """
        MATCH (s:Structure {rcsb_id: $rcsb_id})
        OPTIONAL MATCH (s)-[:DEFINES_ENTITY]->(pe:PolypeptideEntity)
        OPTIONAL MATCH (s)-[:DEFINES_ENTITY]->(ne:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
        OPTIONAL MATCH (s)-[:HAS_INSTANCE]->(pi:PolypeptideInstance)
        OPTIONAL MATCH (s)-[:HAS_INSTANCE]->(ni:NonpolymerInstance)
        WITH s,
             collect(DISTINCT properties(pe)) AS polypeptide_entities,
             collect(DISTINCT {entity: properties(ne), chemical: properties(c)}) AS ligand_entities,
             collect(DISTINCT properties(pi)) AS polypeptide_instances,
             collect(DISTINCT properties(ni)) AS ligand_instances
        RETURN properties(s) AS structure,
               polypeptide_entities,
               ligand_entities,
               polypeptide_instances,
               ligand_instances
        """

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                result = tx.run(query, {"rcsb_id": rcsb_id.upper()}).single()
                if not result:
                    return None
                return {
                    "structure": result["structure"],
                    "polypeptide_entities": [
                        e for e in result["polypeptide_entities"] if e
                    ],
                    "ligand_entities": [
                        e for e in result["ligand_entities"] if e.get("entity")
                    ],
                    "polypeptide_instances": [
                        i for i in result["polypeptide_instances"] if i
                    ],
                    "ligand_instances": [i for i in result["ligand_instances"] if i],
                }

            return session.execute_read(run_query)

    # -------------------------------------------------------------------------
    # Polypeptide Entity Queries
    # -------------------------------------------------------------------------

    def list_polypeptide_entities(
        self, filters: PolypeptideEntityFilters
    ) -> PolypeptideListResponse:
        """
        List polypeptide entities with filtering and pagination.
        """
        builder = PolypeptideEntityQueryBuilder(filters)
        query, params = builder.build()

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                result = tx.run(query, params)
                records = list(result)

                if not records:
                    return PolypeptideListResponse(
                        data=[], total_count=0, next_cursor=None, has_more=False
                    )

                total_count = records[0]["total_count"]
                next_cursor = records[0]["next_cursor"]

                entities = [
                    PolypeptideEntitySummary(
                        parent_rcsb_id=r["parent_rcsb_id"],
                        entity_id=r["entity_id"],
                        pdbx_description=r["pdbx_description"],
                        family=r["family"],
                        sequence_length=r["sequence_length"],
                        src_organism_names=r["src_organism_names"] or [],
                        uniprot_accessions=r["uniprot_accessions"] or [],
                        variant_count=r["variant_count"],
                    )
                    for r in records
                ]

                return PolypeptideListResponse(
                    data=entities,
                    total_count=total_count,
                    next_cursor=next_cursor,
                    has_more=next_cursor is not None,
                )

            return session.execute_read(run_query)

    # -------------------------------------------------------------------------
    # Ligand/Chemical Queries
    # -------------------------------------------------------------------------

    def list_ligands(self, filters: LigandFilters) -> LigandListResponse:
        """
        List chemicals/ligands with filtering and pagination.
        """
        builder = LigandQueryBuilder(filters)
        query, params = builder.build()

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                result = tx.run(query, params)
                records = list(result)

                if not records:
                    return LigandListResponse(
                        data=[], total_count=0, next_cursor=None, has_more=False
                    )

                total_count = records[0]["total_count"]
                next_cursor = records[0]["next_cursor"]

                ligands = [
                    LigandSummary(
                        chemical_id=r["chemical_id"],
                        chemical_name=r["chemical_name"],
                        drugbank_id=r["drugbank_id"],
                        formula_weight=r["formula_weight"],
                        structure_count=r["structure_count"],
                    )
                    for r in records
                ]

                return LigandListResponse(
                    data=ligands,
                    total_count=total_count,
                    next_cursor=next_cursor,
                    has_more=next_cursor is not None,
                )

            return session.execute_read(run_query)

    # -------------------------------------------------------------------------
    # Ligand Neighborhoods
    # -------------------------------------------------------------------------

    def get_ligand_neighborhoods_for_polymer(
        self, rcsb_id: str, auth_asym_id: str
    ) -> PolymerNeighborhoodsResponse:
        """
        Get all ligands in the neighborhood of a polymer chain with their nearby residues.
        """
        query = """
        MATCH (ni:NonpolymerInstance)-[r:NEAR_POLYMER]->(pi:PolypeptideInstance)
        WHERE pi.parent_rcsb_id = $rcsb_id AND pi.auth_asym_id = $auth_asym_id
        MATCH (ni)-[:INSTANCE_OF]->(ne:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
        RETURN 
            c.chemical_id AS ligand_id,
            c.chemical_name AS ligand_name,
            ni.auth_asym_id AS ligand_auth_asym_id,
            r.residues_json AS residues_json,
            r.residue_count AS residue_count,
            c.drugbank_id AS drugbank_id
        ORDER BY c.chemical_id
        """

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                records = list(
                    tx.run(
                        query,
                        {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id},
                    )
                )

                neighborhoods = []
                total_residues = 0

                for r in records:
                    residues_json = r["residues_json"]
                    residues = []
                    if residues_json:
                        raw_residues = json.loads(residues_json)
                        residues = [BindingSiteResidue(**res) for res in raw_residues]

                    neighborhoods.append(
                        LigandNeighborhood(
                            ligand_id=r["ligand_id"],
                            ligand_name=r["ligand_name"],
                            ligand_auth_asym_id=r["ligand_auth_asym_id"],
                            residues=residues,
                            residue_count=r["residue_count"] or 0,
                            drugbank_id=r["drugbank_id"],
                        )
                    )
                    total_residues += r["residue_count"] or 0

                return PolymerNeighborhoodsResponse(
                    rcsb_id=rcsb_id.upper(),
                    auth_asym_id=auth_asym_id,
                    neighborhoods=neighborhoods,
                    total_ligands=len(neighborhoods),
                    total_residues=total_residues,
                )

            return session.execute_read(run_query)

    # -------------------------------------------------------------------------
    # Facet/Aggregation Queries (for filter UI dropdowns)
    # -------------------------------------------------------------------------

    def get_taxonomy_tree(self, tax_type: str = "source") -> List[Dict[str, Any]]:
        """
        Get taxonomy nodes linked to structures for filter dropdowns.
        """
        rel_type = f"belongs_to_lineage_{tax_type}"
        query = f"""
        MATCH (s:Structure)<-[:{rel_type}]-(p:PhylogenyNode)
        WITH p, count(DISTINCT s) AS structure_count
        RETURN p.ncbi_tax_id AS tax_id,
               p.scientific_name AS name,
               p.rank AS rank,
               structure_count
        ORDER BY structure_count DESC
        LIMIT 100
        """

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                return [dict(r) for r in tx.run(query)]

            return session.execute_read(run_query)

    def get_tubulin_families(self) -> List[Dict[str, Any]]:
        """Get tubulin family options with counts"""
        query = """
        MATCH (e:PolypeptideEntity)
        WHERE e.family IS NOT NULL
        WITH e.family AS family, count(*) AS count
        RETURN family, count
        ORDER BY count DESC
        """

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                return [dict(r) for r in tx.run(query)]

            return session.execute_read(run_query)

    # -------------------------------------------------------------------------
    # Simple lookups
    # -------------------------------------------------------------------------

    def all_structure_ids(self) -> List[str]:
        """Get all structure IDs"""
        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                return tx.run("MATCH (s:Structure) RETURN collect(s.rcsb_id)").value()[
                    0
                ]

            return session.execute_read(run_query)

    def structure_exists(self, rcsb_id: str) -> bool:
        """Check if structure exists"""
        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                result = tx.run(
                    "MATCH (s:Structure {rcsb_id: $id}) RETURN count(s) > 0 AS exists",
                    {"id": rcsb_id.upper()},
                ).single()
                return result["exists"] if result else False

            return session.execute_read(run_query)

    def get_taxonomy_tree_for_ui(
        self, tax_type: str = "source"
    ) -> List[Dict[str, Any]]:
        """
        Get taxonomy as a tree structure for antd TreeSelect.
        Returns format: { value, title, children }
        """
        rel_type = f"belongs_to_lineage_{tax_type}"

        query = f"""
        MATCH (s:Structure)<-[:{rel_type}]-(p:PhylogenyNode)
        WITH p, count(DISTINCT s) AS structure_count
        OPTIONAL MATCH (p)-[:descendant_of]->(parent:PhylogenyNode)
        RETURN
            p.ncbi_tax_id AS tax_id,
            p.scientific_name AS name,
            p.rank AS rank,
            structure_count,
            parent.ncbi_tax_id AS parent_id
        ORDER BY structure_count DESC
        """

        with self.adapter.driver.session() as session:

            def run_query(tx):
                records = list(tx.run(query))

                nodes = {}
                for r in records:
                    nodes[r["tax_id"]] = {
                        "value": r["tax_id"],
                        "title": f"{r['name']} ({r['structure_count']})",
                        "rank": r["rank"],
                        "parent_id": r["parent_id"],
                        "children": [],
                    }

                roots = []
                for tax_id, node in nodes.items():
                    parent_id = node.pop("parent_id")
                    node.pop("rank")
                    if parent_id and parent_id in nodes:
                        nodes[parent_id]["children"].append(node)
                    else:
                        roots.append(node)

                def clean_children(node):
                    if not node["children"]:
                        del node["children"]
                    else:
                        for child in node["children"]:
                            clean_children(child)
                    return node

                return [clean_children(r) for r in roots]

            return session.execute_read(run_query)

    def get_ligand_options(
        self, search: Optional[str] = None, limit: int = 50
    ) -> List[Dict[str, Any]]:
        """
        Get ligands for filter dropdown/autocomplete.
        Returns lightweight list: chemical_id, name, structure_count.
        """
        query = """
        MATCH (c:Chemical)
        """
        params: Dict[str, Any] = {"limit": limit}

        if search:
            query += """
            WHERE toLower(c.chemical_id) CONTAINS $search
            OR toLower(c.chemical_name) CONTAINS $search
            """
            params["search"] = search.lower()

        query += """
        OPTIONAL MATCH (c)<-[:DEFINED_BY_CHEMICAL]-(ne:NonpolymerEntity)<-[:DEFINES_ENTITY]-(s:Structure)
        WITH c, count(DISTINCT s) AS structure_count
        WHERE structure_count > 0
        RETURN
            c.chemical_id AS chemical_id,
            c.chemical_name AS chemical_name,
            structure_count
        ORDER BY structure_count DESC
        LIMIT $limit
        """

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                return [dict(r) for r in tx.run(query, params)]

            return session.execute_read(run_query)

    def get_filter_facets(self) -> FilterFacets:
        """
        Get available filter options for the UI.
        Returns counts for each facet value.
        """
        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                # Base counts
                base_result = tx.run("""
                    MATCH (s:Structure)
                    RETURN count(s) AS total_structures,
                           min(s.citation_year) AS min_year,
                           max(s.citation_year) AS max_year,
                           min(s.resolution) AS min_res,
                           max(s.resolution) AS max_res
                """).single()

                total_structures = base_result["total_structures"]
                year_range = RangeValue(
                    min=base_result["min_year"], max=base_result["max_year"]
                )
                resolution_range = RangeValue(
                    min=base_result["min_res"], max=base_result["max_res"]
                )

                # Experimental methods
                exp_methods = [
                    FacetValue(value=r["value"], count=r["count"])
                    for r in tx.run("""
                        MATCH (s:Structure)
                        WHERE s.expMethod IS NOT NULL
                        RETURN s.expMethod AS value, count(*) AS count
                        ORDER BY count DESC
                    """)
                ]

                # Tubulin families
                tubulin_families = [
                    FacetValue(value=r["value"], count=r["count"])
                    for r in tx.run("""
                        MATCH (e:PolypeptideEntity)
                        WHERE e.family IS NOT NULL
                        RETURN e.family AS value, count(DISTINCT e.parent_rcsb_id) AS count
                        ORDER BY count DESC
                    """)
                ]

                # Top ligands
                top_ligands = [
                    LigandFacet(
                        chemical_id=r["chemical_id"],
                        chemical_name=r["chemical_name"],
                        count=r["count"],
                    )
                    for r in tx.run("""
                        MATCH (c:Chemical)<-[:DEFINED_BY_CHEMICAL]-(ne:NonpolymerEntity)<-[:DEFINES_ENTITY]-(s:Structure)
                        WITH c, count(DISTINCT s) AS cnt
                        ORDER BY cnt DESC
                        LIMIT 50
                        RETURN c.chemical_id AS chemical_id, c.chemical_name AS chemical_name, cnt AS count
                    """)
                ]

                # Variant stats by family
                variants_by_family = [
                    VariantsByFamily(
                        family=r["family"],
                        variant_count=r["variant_count"],
                        structure_count=r["structure_count"],
                    )
                    for r in tx.run("""
                        MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
                        WHERE e.family IS NOT NULL
                        WITH e.family AS family, count(DISTINCT v) AS variant_count, count(DISTINCT e.parent_rcsb_id) AS structure_count
                        RETURN family, variant_count, structure_count
                        ORDER BY variant_count DESC
                    """)
                ]

                # Common variants (top 30)
                common_variants = [
                    CommonVariant(
                        family=r["family"],
                        variant_type=r["variant_type"],
                        position=r["position"],
                        wild_type=r["wild_type"],
                        observed=r["observed"],
                        count=r["cnt"],
                    )
                    for r in tx.run("""
                        MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
                        WITH e.family AS family, v.type AS variant_type, v.master_index AS position, 
                             v.wild_type AS wild_type, v.observed AS observed, count(*) AS cnt
                        ORDER BY cnt DESC
                        LIMIT 30
                        RETURN family, variant_type, position, wild_type, observed, cnt
                    """)
                ]

                # Variant position ranges by family
                variant_position_ranges = [
                    VariantPositionRange(
                        family=r["family"],
                        min_position=r["min_pos"],
                        max_position=r["max_pos"],
                    )
                    for r in tx.run("""
                        MATCH (e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
                        WHERE e.family IS NOT NULL AND v.master_index IS NOT NULL
                        WITH e.family AS family, min(v.master_index) AS min_pos, max(v.master_index) AS max_pos
                        RETURN family, min_pos, max_pos
                    """)
                ]

                return FilterFacets(
                    total_structures=total_structures,
                    exp_methods=exp_methods,
                    tubulin_families=tubulin_families,
                    year_range=year_range,
                    resolution_range=resolution_range,
                    top_ligands=top_ligands,
                    variants_by_family=variants_by_family,
                    common_variants=common_variants,
                    variant_position_ranges=variant_position_ranges,
                )

            return session.execute_read(run_query)


# Singleton instance for convenience
db_reader = Neo4jReader()
