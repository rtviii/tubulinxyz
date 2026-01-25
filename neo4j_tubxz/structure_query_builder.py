# neo4j_tubxz/structure_query_builder.py
"""
Cypher query builders for filtering structures and entities.
"""

from typing import Dict, Any, Tuple, List
from neo4j_tubxz.models import StructureFilters, PolypeptideEntityFilters, LigandFilters


class StructureQueryBuilder:
    """
    Builds Cypher queries for filtering Structures with keyset pagination.
    """

    def __init__(self, filters: StructureFilters):
        self.filters = filters
        self._where_clauses: List[str] = []
        self._params: Dict[str, Any] = {"limit": filters.limit}

    def build(self) -> Tuple[str, Dict[str, Any]]:
        self._process_filters()
        return self._assemble(), self._params

    def _process_filters(self):
        f = self.filters

        if f.search:
            self._add_search_filter(f.search)

        if f.rcsb_ids:
            self._where_clauses.append("s.rcsb_id IN $rcsb_ids")
            self._params["rcsb_ids"] = [x.upper() for x in f.rcsb_ids]

        if f.resolution_min is not None:
            self._where_clauses.append("s.resolution >= $res_min")
            self._params["res_min"] = f.resolution_min

        if f.resolution_max is not None:
            self._where_clauses.append("s.resolution <= $res_max")
            self._params["res_max"] = f.resolution_max

        if f.year_min is not None:
            self._where_clauses.append("s.citation_year >= $year_min")
            self._params["year_min"] = f.year_min

        if f.year_max is not None:
            self._where_clauses.append("s.citation_year <= $year_max")
            self._params["year_max"] = f.year_max

        if f.exp_method:
            self._where_clauses.append("s.expMethod IN $exp_methods")
            self._params["exp_methods"] = [m.value for m in f.exp_method]

        if f.polymerization_state:
            self._where_clauses.append("s.polymerization_state IN $poly_states")
            self._params["poly_states"] = [p.value for p in f.polymerization_state]

        if f.source_organism_ids:
            self._add_taxonomy_filter(f.source_organism_ids, "source")

        if f.host_organism_ids:
            self._add_taxonomy_filter(f.host_organism_ids, "host")

        if f.has_ligand_ids:
            self._add_ligand_filter(f.has_ligand_ids)

        if f.has_polymer_family:
            self._add_polymer_family_filter(f.has_polymer_family)

        if f.has_uniprot:
            self._add_uniprot_filter(f.has_uniprot)

        self._add_variant_filters()

        if f.cursor:
            self._where_clauses.append("s.rcsb_id < $cursor")
            self._params["cursor"] = f.cursor.upper()

    def _add_search_filter(self, search: str):
        self._where_clauses.append("""(
            toLower(s.rcsb_id) CONTAINS $search OR
            toLower(s.citation_title) CONTAINS $search OR
            toLower(s.pdbx_keywords_text) CONTAINS $search OR
            ANY(name IN s.src_organism_names WHERE toLower(name) CONTAINS $search)
        )""")
        self._params["search"] = search.lower()

    def _add_taxonomy_filter(self, tax_ids: List[int], tax_type: str):
        rel_type = f"belongs_to_lineage_{tax_type}"
        param_name = f"{tax_type}_taxa"
        self._where_clauses.append(f"""
            EXISTS {{
                MATCH (s)<-[:{rel_type}]-(p:PhylogenyNode)
                WHERE p.ncbi_tax_id IN ${param_name}
            }}
        """)
        self._params[param_name] = tax_ids

    def _add_ligand_filter(self, chemical_ids: List[str]):
        for i, chem_id in enumerate(chemical_ids):
            param_name = f"ligand_id_{i}"
            self._where_clauses.append(f"""
                EXISTS {{
                    MATCH (s)-[:DEFINES_ENTITY]->(e:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
                    WHERE c.chemical_id = ${param_name}
                }}
            """)
            self._params[param_name] = chem_id.upper()

    def _add_polymer_family_filter(self, families: List[str]):
        for i, family in enumerate(families):
            param_name = f"polymer_family_{i}"
            self._where_clauses.append(f"""
                EXISTS {{
                    MATCH (s)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)
                    WHERE e.family = ${param_name}
                }}
            """)
            self._params[param_name] = family

    def _add_uniprot_filter(self, accessions: List[str]):
        self._where_clauses.append("""
            EXISTS {
                MATCH (s)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)
                WHERE ANY(acc IN e.uniprot_accessions WHERE acc IN $uniprot_ids)
            }
        """)
        self._params["uniprot_ids"] = accessions

    def _add_variant_filters(self):
        f = self.filters
        variant_conditions = []

        if f.variant_family:
            variant_conditions.append("e.family = $variant_family")
            self._params["variant_family"] = f.variant_family

        if f.variant_type:
            variant_conditions.append("v.type = $variant_type")
            self._params["variant_type"] = f.variant_type.value

        if f.variant_position_min is not None:
            variant_conditions.append("v.master_index >= $var_pos_min")
            self._params["var_pos_min"] = f.variant_position_min

        if f.variant_position_max is not None:
            variant_conditions.append("v.master_index <= $var_pos_max")
            self._params["var_pos_max"] = f.variant_position_max

        if f.variant_wild_type:
            variant_conditions.append("v.wild_type = $var_wild_type")
            self._params["var_wild_type"] = f.variant_wild_type.upper()

        if f.variant_observed:
            variant_conditions.append("v.observed = $var_observed")
            self._params["var_observed"] = f.variant_observed.upper()

        if f.variant_source:
            variant_conditions.append("v.source = $var_source")
            self._params["var_source"] = f.variant_source

        if f.variant_phenotype:
            variant_conditions.append("toLower(v.phenotype) CONTAINS $var_phenotype")
            self._params["var_phenotype"] = f.variant_phenotype.lower()

        if variant_conditions:
            where_clause = " AND ".join(variant_conditions)
            self._where_clauses.append(f"""
                EXISTS {{
                    MATCH (s)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)-[:HAS_VARIANT]->(v:Variant)
                    WHERE {where_clause}
                }}
            """)

        if f.has_variants is not None:
            if f.has_variants:
                if f.variant_family:
                    self._where_clauses.append("""
                        EXISTS {
                            MATCH (s)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)-[:HAS_VARIANT]->(:Variant)
                            WHERE e.family = $has_var_family
                        }
                    """)
                    self._params["has_var_family"] = f.variant_family
                else:
                    self._where_clauses.append("""
                        EXISTS {
                            MATCH (s)-[:DEFINES_ENTITY]->(:PolypeptideEntity)-[:HAS_VARIANT]->(:Variant)
                        }
                    """)
            else:
                if f.variant_family:
                    self._where_clauses.append("""
                        NOT EXISTS {
                            MATCH (s)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)-[:HAS_VARIANT]->(:Variant)
                            WHERE e.family = $no_var_family
                        }
                    """)
                    self._params["no_var_family"] = f.variant_family
                else:
                    self._where_clauses.append("""
                        NOT EXISTS {
                            MATCH (s)-[:DEFINES_ENTITY]->(:PolypeptideEntity)-[:HAS_VARIANT]->(:Variant)
                        }
                    """)

    def _assemble(self) -> str:
        parts = ["MATCH (s:Structure)"]

        if self._where_clauses:
            parts.append("WHERE " + "\n  AND ".join(self._where_clauses))

        parts.append("WITH s ORDER BY s.rcsb_id DESC")
        parts.append("""
WITH collect(s) AS all_results
WITH all_results, size(all_results) AS total_count
WITH total_count, all_results[0..$limit] AS page
WITH total_count, page,
     CASE WHEN size(page) = $limit AND size(page) > 0
          THEN page[-1].rcsb_id
          ELSE null
     END AS next_cursor
UNWIND page AS s
OPTIONAL MATCH (s)-[:DEFINES_ENTITY]->(e:Entity)
WITH total_count, next_cursor, s, count(e) AS entity_count
OPTIONAL MATCH (s)-[:DEFINES_ENTITY]->(ne:NonpolymerEntity)
WITH total_count, next_cursor, s, entity_count, count(ne) AS ligand_count
RETURN
    s.rcsb_id AS rcsb_id,
    s.resolution AS resolution,
    s.expMethod AS exp_method,
    s.citation_title AS citation_title,
    s.citation_year AS citation_year,
    s.deposition_date AS deposition_date,
    s.src_organism_names AS src_organism_names,
    s.pdbx_keywords AS pdbx_keywords,
    entity_count,
    ligand_count,
    total_count,
    next_cursor
        """)

        return "\n".join(parts)


class PolypeptideEntityQueryBuilder:
    """Builds Cypher queries for filtering PolypeptideEntity nodes."""

    def __init__(self, filters: PolypeptideEntityFilters):
        self.filters = filters
        self._where_clauses: List[str] = []
        self._struct_where_clauses: List[str] = []
        self._params: Dict[str, Any] = {"limit": filters.limit}

    def build(self) -> Tuple[str, Dict[str, Any]]:
        self._process_filters()
        return self._assemble(), self._params

    def _process_filters(self):
        f = self.filters

        if f.structure_filters:
            self._apply_structure_filters(f.structure_filters)

        if f.family:
            self._where_clauses.append("e.family IN $families")
            self._params["families"] = f.family

        if f.uniprot_accession:
            self._where_clauses.append("$uniprot IN e.uniprot_accessions")
            self._params["uniprot"] = f.uniprot_accession

        if f.sequence_contains:
            self._where_clauses.append("e.one_letter_code_can CONTAINS $motif")
            self._params["motif"] = f.sequence_contains.upper()

        if f.sequence_length_min is not None:
            self._where_clauses.append("e.sequence_length >= $seq_len_min")
            self._params["seq_len_min"] = f.sequence_length_min

        if f.sequence_length_max is not None:
            self._where_clauses.append("e.sequence_length <= $seq_len_max")
            self._params["seq_len_max"] = f.sequence_length_max

        if f.has_variants is not None:
            if f.has_variants:
                self._where_clauses.append(
                    "EXISTS { MATCH (e)-[:HAS_VARIANT]->(:Variant) }"
                )
            else:
                self._where_clauses.append(
                    "NOT EXISTS { MATCH (e)-[:HAS_VARIANT]->(:Variant) }"
                )

        if f.cursor:
            parts = f.cursor.split(":")
            if len(parts) == 2:
                rcsb_id, entity_id = parts
                self._where_clauses.append("""
                    (e.parent_rcsb_id < $cursor_rcsb OR
                     (e.parent_rcsb_id = $cursor_rcsb AND e.entity_id > $cursor_entity))
                """)
                self._params["cursor_rcsb"] = rcsb_id.upper()
                self._params["cursor_entity"] = entity_id

    def _apply_structure_filters(self, sf: StructureFilters):
        if sf.resolution_min is not None:
            self._struct_where_clauses.append("s.resolution >= $res_min")
            self._params["res_min"] = sf.resolution_min

        if sf.resolution_max is not None:
            self._struct_where_clauses.append("s.resolution <= $res_max")
            self._params["res_max"] = sf.resolution_max

        if sf.year_min is not None:
            self._struct_where_clauses.append("s.citation_year >= $year_min")
            self._params["year_min"] = sf.year_min

        if sf.year_max is not None:
            self._struct_where_clauses.append("s.citation_year <= $year_max")
            self._params["year_max"] = sf.year_max

        if sf.source_organism_ids:
            self._struct_where_clauses.append("""
                EXISTS {
                    MATCH (s)<-[:belongs_to_lineage_source]-(p:PhylogenyNode)
                    WHERE p.ncbi_tax_id IN $source_taxa
                }
            """)
            self._params["source_taxa"] = sf.source_organism_ids

    def _assemble(self) -> str:
        parts = ["MATCH (s:Structure)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)"]

        all_where = self._struct_where_clauses + self._where_clauses
        if all_where:
            parts.append("WHERE " + "\n  AND ".join(all_where))

        parts.append("""
WITH e, s
ORDER BY e.parent_rcsb_id DESC, e.entity_id ASC
WITH collect({entity: e, struct: s}) AS all_results
WITH all_results, size(all_results) AS total_count
WITH total_count, all_results[0..$limit] AS page
WITH total_count, page,
     CASE WHEN size(page) = $limit AND size(page) > 0
          THEN page[-1].entity.parent_rcsb_id + ':' + page[-1].entity.entity_id
          ELSE null
     END AS next_cursor
UNWIND page AS row
WITH total_count, next_cursor, row.entity AS e, row.struct AS s
OPTIONAL MATCH (e)-[:HAS_VARIANT]->(v:Variant)
WITH total_count, next_cursor, e, count(v) AS variant_count
RETURN
    e.parent_rcsb_id AS parent_rcsb_id,
    e.entity_id AS entity_id,
    e.pdbx_description AS pdbx_description,
    e.family AS family,
    e.sequence_length AS sequence_length,
    e.src_organism_names AS src_organism_names,
    e.uniprot_accessions AS uniprot_accessions,
    variant_count,
    total_count,
    next_cursor
        """)

        return "\n".join(parts)


class LigandQueryBuilder:
    """Builds Cypher queries for Chemical/Ligand filtering."""

    def __init__(self, filters: LigandFilters):
        self.filters = filters
        self._where_clauses: List[str] = []
        self._params: Dict[str, Any] = {"limit": filters.limit}

    def build(self) -> Tuple[str, Dict[str, Any]]:
        self._process_filters()
        return self._assemble(), self._params

    def _process_filters(self):
        f = self.filters

        if f.search:
            self._where_clauses.append("""(
                toLower(c.chemical_id) CONTAINS $search OR
                toLower(c.chemical_name) CONTAINS $search
            )""")
            self._params["search"] = f.search.lower()

        if f.chemical_ids:
            self._where_clauses.append("c.chemical_id IN $chem_ids")
            self._params["chem_ids"] = [x.upper() for x in f.chemical_ids]

        if f.has_drugbank is not None:
            if f.has_drugbank:
                self._where_clauses.append("c.drugbank_id IS NOT NULL")
            else:
                self._where_clauses.append("c.drugbank_id IS NULL")

        if f.in_structures:
            self._where_clauses.append("""
                EXISTS {
                    MATCH (c)<-[:DEFINED_BY_CHEMICAL]-(:NonpolymerEntity)<-[:DEFINES_ENTITY]-(s:Structure)
                    WHERE s.rcsb_id IN $struct_ids
                }
            """)
            self._params["struct_ids"] = [x.upper() for x in f.in_structures]

        if f.cursor:
            self._where_clauses.append("c.chemical_id > $cursor")
            self._params["cursor"] = f.cursor.upper()

    def _assemble(self) -> str:
        parts = ["MATCH (c:Chemical)"]

        if self._where_clauses:
            parts.append("WHERE " + "\n  AND ".join(self._where_clauses))

        parts.append("""
WITH c ORDER BY c.chemical_id ASC
WITH collect(c) AS all_results
WITH all_results, size(all_results) AS total_count
WITH total_count, all_results[0..$limit] AS page
WITH total_count, page,
     CASE WHEN size(page) = $limit AND size(page) > 0
          THEN page[-1].chemical_id
          ELSE null
     END AS next_cursor
UNWIND page AS c
OPTIONAL MATCH (c)<-[:DEFINED_BY_CHEMICAL]-(ne:NonpolymerEntity)<-[:DEFINES_ENTITY]-(s:Structure)
WITH total_count, next_cursor, c, count(DISTINCT s) AS structure_count
RETURN
    c.chemical_id AS chemical_id,
    c.chemical_name AS chemical_name,
    c.drugbank_id AS drugbank_id,
    c.formula_weight AS formula_weight,
    structure_count,
    total_count,
    next_cursor
        """)

        return "\n".join(parts)
