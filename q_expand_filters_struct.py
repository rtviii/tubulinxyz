So this is a database application for tubulin structures, their polymers, ligands, mutations, annotations etc. We have a neo4j database, fastapi in the backend and a react/nextjs in the frontend. 

Currently we are building out the infrastructure for the backend to be able to pass info from the db to the frontend based on the UI's queries and filters.

-----------

Our immediate objective is to move to a structured **API Communication Layer**. This layer will serve as a typed contract between the Neo4j graph data and the React/Redux frontend, utilizing FastAPI’s ability to generate OpenAPI schemas for automatic client-side code generation.

Let's begin this by prototyping a very good set of frontend filters for the structures of tubulin and their facets (all the data nad connections that we can associate with them ).

### The Goal: Dynamic, Type-Safe Filtering

We need a **Cumulative Filtering System** designed to handle complex, real-time refetching.

* **Cumulative Query Logic**: Every filter parameter (e.g., tubulin family, resolution, specific ligands) must be additive. The backend must dynamically append `MATCH` and `WHERE` clauses to the Cypher query based on the active state of the frontend UI.

* **Keyset Pagination (The Cursor)**: To ensure performance at scale, we will avoid `SKIP` and use the `rcsb_id` as a sortable "anchor" or cursor. The query will always fetch the next N results where the ID is greater than the provided cursor.

* **Real-time Refetching**: The system must be optimized for "as-you-type" or "as-you-toggle" updates, relying on indexed database lookups to maintain sub-second response times.

---

### Enforcing the Data Contract (FastAPI & RTK-Query)

To ensure that our Redux Toolkit (RTK) Query hooks are perfectly typed on the frontend, we will implement strict **Pydantic models** for every input and output

* **Request Models (`POST` / `GET` params)**: Instead of loose query parameters, we will define `FilterRequest` models. This allows RTK-Query to auto-generate interfaces that prevent you from passing a string where a number is expected.

* **Response Models (`PaginatedResponse`)**: Every response will follow a standard wrapper that includes the data list, total count, and the `next_cursor`.


----

Let me show you my database schema now and perhaps we can start playing with this.
```
[
  {
    "label": "descendant_of",
    "property": "PhylogenyNode",
    "count": 389,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["PhylogenyNode"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "source",
    "property": "PhylogenyNode",
    "count": 85,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 20,
    "right": 0,
    "other": ["Structure"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "belongs_to_lineage_source",
    "property": "PhylogenyNode",
    "count": 463,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 68,
    "right": 0,
    "other": ["Structure"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "DEFINES_ENTITY",
    "property": "Structure",
    "count": 855,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 9,
    "right": 0,
    "other": [
      "Entity",
      "NonpolymerEntity",
      "PolypeptideEntity",
      "PolynucleotideEntity"
    ],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "ANNOTATES_POSITION_IN",
    "property": "Mutation",
    "count": 9,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["MasterAlignment"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "HAS_MUTATION",
    "property": "Entity",
    "count": 174,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 16,
    "right": 0,
    "other": ["Mutation"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "HAS_MUTATION",
    "property": "PolypeptideEntity",
    "count": 465,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 19,
    "right": 0,
    "other": ["Mutation"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "DEFINED_BY_CHEMICAL",
    "property": "Entity",
    "count": 554,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Chemical"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "DEFINED_BY_CHEMICAL",
    "property": "NonpolymerEntity",
    "count": 1903,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Chemical"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "INSTANCE_OF",
    "property": "Instance",
    "count": 1160,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "PolypeptideEntity", "NonpolymerEntity"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "INSTANCE_OF",
    "property": "PolypeptideInstance",
    "count": 1019,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "PolypeptideEntity"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "INSTANCE_OF",
    "property": "NonpolymerInstance",
    "count": 978,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "NonpolymerEntity"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "INSTANCE_OF",
    "property": "PolynucleotideInstance",
    "count": 10,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "PolynucleotideEntity"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "HAS_INSTANCE",
    "property": "Structure",
    "count": 854,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 50,
    "right": 0,
    "other": [
      "Instance",
      "NonpolymerInstance",
      "PolypeptideInstance",
      "PolynucleotideInstance"
    ],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "HAS_INTERACTION",
    "property": "Instance",
    "count": 533,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 26,
    "right": 0,
    "other": ["Interaction"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "HAS_INTERACTION",
    "property": "NonpolymerInstance",
    "count": 902,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 26,
    "right": 0,
    "other": ["Interaction"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "INTERACTS_WITH",
    "property": "Interaction",
    "count": 1087,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Instance", "PolypeptideInstance"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "NEAR_POLYMER",
    "property": "Instance",
    "count": 533,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Instance", "PolypeptideInstance"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "NEAR_POLYMER",
    "property": "residues",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "NEAR_POLYMER",
    "property": "NonpolymerInstance",
    "count": 902,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Instance", "PolypeptideInstance"],
    "otherLabels": [],
    "elementType": "relationship"
  },
  {
    "label": "Structure",
    "property": "HAS_INSTANCE",
    "count": 854,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 50,
    "right": 0,
    "other": [
      "Instance",
      "NonpolymerInstance",
      "PolypeptideInstance",
      "PolynucleotideInstance"
    ],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "DEFINES_ENTITY",
    "count": 855,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 9,
    "right": 0,
    "other": [
      "Entity",
      "NonpolymerEntity",
      "PolypeptideEntity",
      "PolynucleotideEntity"
    ],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "deposition_date",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "citation_title",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "citation_rcsb_authors",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "citation_pdbx_doi",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "citation_year",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "host_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "rcsb_external_ref_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "rcsb_external_ref_type",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "rcsb_external_ref_link",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "pdbx_keywords_text",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "pdbx_keywords",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "host_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "resolution",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "FLOAT",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "src_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "src_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "rcsb_id",
    "count": 0,
    "unique": true,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Structure",
    "property": "expMethod",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "one_letter_code_can",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "src_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "entity_id",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "pdbx_description",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "sequence_length",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "one_letter_code",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "src_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "host_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "host_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "uniprot_accessions",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "family",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "DEFINED_BY_CHEMICAL",
    "count": 554,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Chemical"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "chemical_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Entity",
    "property": "HAS_MUTATION",
    "count": 174,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 16,
    "right": 0,
    "other": ["Mutation"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "INSTANCE_OF",
    "count": 1160,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "PolypeptideEntity", "NonpolymerEntity"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "asym_id",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "assembly_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "auth_asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "HAS_INTERACTION",
    "count": 533,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 26,
    "right": 0,
    "other": ["Interaction"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "NEAR_POLYMER",
    "count": 533,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Instance", "PolypeptideInstance"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Instance",
    "property": "auth_seq_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "SMILES_stereo",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "chemical_name",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "InChI",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "SMILES",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "drugbank_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "chemical_id",
    "count": 0,
    "unique": true,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Chemical",
    "property": "InChIKey",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PhylogenyNode",
    "property": "source",
    "count": 85,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 20,
    "right": 0,
    "other": ["Structure"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PhylogenyNode",
    "property": "belongs_to_lineage_source",
    "count": 463,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 68,
    "right": 0,
    "other": ["Structure"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PhylogenyNode",
    "property": "rank",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PhylogenyNode",
    "property": "ncbi_tax_id",
    "count": 0,
    "unique": true,
    "index": true,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PhylogenyNode",
    "property": "scientific_name",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PhylogenyNode",
    "property": "descendant_of",
    "count": 389,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["PhylogenyNode"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "MasterAlignment",
    "property": "family",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "MasterAlignment",
    "property": "version",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "MasterAlignment",
    "property": "description",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "ANNOTATES_POSITION_IN",
    "count": 9,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["MasterAlignment"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "tubulin_type",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "database_source",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "uniprot_id",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "species",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "keywords",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "to_residue",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "phenotype",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "master_index",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "from_residue",
    "count": 0,
    "unique": false,
    "index": true,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Mutation",
    "property": "reference_link",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "one_letter_code_can",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "src_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "entity_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "pdbx_description",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "sequence_length",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "one_letter_code",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "src_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "host_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "host_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "uniprot_accessions",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "HAS_MUTATION",
    "count": 465,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 19,
    "right": 0,
    "other": ["Mutation"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideEntity",
    "property": "family",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerEntity",
    "property": "DEFINED_BY_CHEMICAL",
    "count": 1903,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Chemical"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerEntity",
    "property": "entity_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerEntity",
    "property": "pdbx_description",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerEntity",
    "property": "chemical_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerEntity",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideInstance",
    "property": "INSTANCE_OF",
    "count": 1019,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "PolypeptideEntity"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideInstance",
    "property": "asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideInstance",
    "property": "assembly_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideInstance",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolypeptideInstance",
    "property": "auth_asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "INSTANCE_OF",
    "count": 978,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "NonpolymerEntity"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "assembly_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "auth_asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "HAS_INTERACTION",
    "count": 902,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 26,
    "right": 0,
    "other": ["Interaction"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "NEAR_POLYMER",
    "count": 902,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": true,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Instance", "PolypeptideInstance"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "NonpolymerInstance",
    "property": "auth_seq_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Interaction",
    "property": "INTERACTS_WITH",
    "count": 1087,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Instance", "PolypeptideInstance"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Interaction",
    "property": "type",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Interaction",
    "property": "auth_seq_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Interaction",
    "property": "auth_comp_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Interaction",
    "property": "atom_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "Interaction",
    "property": "master_index",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "entity_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "polymer_type",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "sequence_length",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "one_letter_code",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "src_organism_ids",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "src_organism_names",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "LIST",
    "array": true,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideEntity",
    "property": "one_letter_code_can",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideInstance",
    "property": "INSTANCE_OF",
    "count": 10,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "RELATIONSHIP",
    "array": false,
    "sample": null,
    "left": 1,
    "right": 0,
    "other": ["Entity", "PolynucleotideEntity"],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideInstance",
    "property": "asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideInstance",
    "property": "assembly_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "INTEGER",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideInstance",
    "property": "parent_rcsb_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  },
  {
    "label": "PolynucleotideInstance",
    "property": "auth_asym_id",
    "count": 0,
    "unique": false,
    "index": false,
    "existence": false,
    "type": "STRING",
    "array": false,
    "sample": null,
    "left": 0,
    "right": 0,
    "other": [],
    "otherLabels": [],
    "elementType": "node"
  }
]
```



I will further attach the shape of my backend application and, crucially, some type definition files for the models in which i operate. If we can reuse those across the application -- that'd be awesome (don't feel too constrained by this, but just know that indeed all data in the neo4j are these files and we can mix and match them as we please.).
```
25 directories, 124 files
(venv) ᢹ saeta.rtviii[ dev/tubulinxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.json|*.npy|*.ply|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta' .
.
├── audit_mappings.py
├── augment_lig_interactions_with_ma.py
├── check_ligand_aug.py
├── cli.py
├── data
│   ├── alpha_tubulin
│   ├── beta_tubulin
│   ├── hmms
│   │   ├── classification_cache
│   │   ├── maps
│   │   └── tubulin
│   ├── maxim_data
│   │   ├── 7sj7_tails.pdb
│   │   ├── 7sj7_with_metadata.cif
│   │   ├── add_computed_res_annotation.py
│   │   ├── curved_hum_wt_GDP_6S8L.pdb
│   │   ├── fold_htuba1a_model_0.cif
│   │   ├── fold_htuba1a.zip
│   │   ├── mini_mt_patch_hum_wt_GDP_7SJ7.pdb
│   │   ├── README.txt
│   │   └── straight_hum_wt_GDP_7SJ7.pdb
│   └── sequences
│       ├── maps
│       └── tubulin
├── dirtocontext.py
├── docs
├── docs.md
├── ingest_morisette_data.py
├── ingestion_logs
├── init_all_profiles.py
├── init_db.py
├── lib
│   ├── etl
│   │   ├── assets.py
│   │   ├── collector.py
│   │   ├── constants.py
│   │   ├── libtax.py
│   │   └── ligand_extraction.py
│   ├── ingestion_logs
│   ├── seq_aligner.py
│   ├── tubulin_analyzer
│   │   ├── __init__.py
│   │   ├── geometric_analyzer.py
│   │   ├── interface_analyzer.py
│   │   ├── models.py
│   │   ├── structural_analyzer.py
│   │   └── visualization.py
│   └── types.py
├── map_sequences
├── muscle3.8.1
├── neo4j_tubxz
│   ├── db_driver.py
│   ├── db_lib_builder.py
│   ├── db_lib_reader.py
│   ├── ingest_structure_mutation.py
│   ├── node_interaction.py
│   ├── node_ligand.py
│   ├── node_master_alignment.py
│   ├── node_modification.py
│   ├── node_mutation.py
│   ├── node_phylogeny.py
│   ├── node_polymer.py
│   ├── node_structure.py
│   ├── query_builder.py
│   └── test_one_struct.py
├── notes
│   ├── 0_general_context.md
│   ├── hmm_building.md
│   ├── hmm_classifier_eval_report.txt
│   ├── lig_seq_id_ingestion.md
│   ├── ligand_census.md
│   ├── MAP_report.md
│   ├── minutes_november24th.md
│   ├── notes.md
│   └── tubulin_classes_seq_clustering.md
├── plan.md
├── q_api_design.py
├── requirements.txt
├── scripts_and_artifacts
│   ├── archive
│   │   ├── analyze_ligands.py
│   │   ├── tubulin_ligand_census.py
│   │   └── visualize_ligands.py
│   ├── extract_ixs.tsx
│   ├── hmm_building
│   │   ├── analyze_cluster_composition.py
│   │   ├── eval_classifier.py
│   │   ├── fetch_mipmaps.py
│   │   ├── fetch_tubulin.py
│   │   ├── process_mipmaps.py
│   │   ├── process_tubulin.py
│   │   └── validate_pdb_structs.py
│   └── morisette_stuff
│       ├── morisette_alpha_beta_gamma_uniprot.md
│       ├── morisette_alpha.py
│       ├── mset_consensus.py
│       └── mset_parser.py
├── taxdump.tar.gz
├── validate_batch.py
├── validate_ligand_augmentation.py
├── verify_mapping.py
└── yarn.lock

25 directories, 78 files

```

```
# types_tubulin.py
from typing import Dict, Optional, List, Literal, Tuple, Union, Any
from enum import Enum
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, Field

class TubulinFamily(str, Enum):
    ALPHA   = "tubulin_alpha"
    BETA    = "tubulin_beta"
    GAMMA   = "tubulin_gamma"
    DELTA   = "tubulin_delta"
    EPSILON = "tubulin_epsilon"

class MapFamily(str, Enum):
    """
    Microtubule Associated Proteins (MAPs) and related enzymes.
    Values correspond to the filename base.
    """
    ATAT1           = "map_atat1"
    CAMSAP1         = "map_camsap1"
    CAMSAP2         = "map_camsap2"
    CAMSAP3         = "map_camsap3"
    CCP             = "map_ccp_deglutamylase"
    CFAP53          = "map_cfap53"
    CKAP5           = "map_ckap5_chtog"
    CLASP           = "map_clasp"
    CLIP115         = "map_clip115"
    CLIP170         = "map_clip170"
    DOUBLECORTIN    = "map_doublecortin"
    EB_FAMILY       = "map_eb_family"
    FAP20           = "map_fap20_cfap20"
    GCP2_3          = "map_gcp2_3"
    GCP4            = "map_gcp4"
    GCP5_6          = "map_gcp5_6"
    KATANIN         = "map_katanin_p60"
    KINESIN13       = "map_kinesin13"
    MAP1_HEAVY      = "map_map1_heavy"
    MAP1S           = "map_map1s"
    MAP2            = "map_map2"
    MAP4            = "map_map4"
    MAP7            = "map_map7"
    NME7            = "map_nme7"
    NME8            = "map_nme8"
    NUMA            = "map_numa"
    PACRG           = "map_pacrg"
    PRC1            = "map_prc1"
    RIB72           = "map_rib72_efhc"
    SPAG6           = "map_spag6"
    SPASTIN         = "map_spastin"
    STATHMIN        = "map_stathmin"
    TACC            = "map_tacc"
    TAU             = "map_tau"
    TPX2            = "map_tpx2"
    TTLL_LONG       = "map_ttll_glutamylase_long"
    TTLL_SHORT      = "map_ttll_glutamylase_short"
    VASH            = "map_vash_detyrosinase"

# Type alias for any family
PolymerClass = Union[TubulinFamily, MapFamily]

# --- Ligand Interaction Models ---

class InteractionType(str, Enum):

    UNKNOWN            = "Unknown"
    IONIC              = "Ionic"
    CATION_PI          = "Cation-Pi Interaction"
    PI_STACKING        = "Pi Stacking"
    HYDROGEN_BOND      = "Hydrogen Bond"
    HALOGEN_BOND       = "Halogen Bond"
    HYDROPHOBIC        = "Hydrophobic Contact"
    METAL_COORDINATION = "Metal Coordination"
    WEAK_HYDROGEN_BOND = "Weak Hydrogen Bond"

class InteractionParticipant(BaseModel):
    """An atom participating in an interaction."""
    auth_asym_id: str
    auth_seq_id : int
    auth_comp_id: str
    atom_id     : str
    is_ligand   : bool
    master_index: Optional[int] = None  # Added field, made optional

    @classmethod
    def from_tuple(cls, t: list) -> "InteractionParticipant":
        return cls(
            auth_asym_id=t[0],
            auth_seq_id=t[1],
            auth_comp_id=t[2],
            atom_id=t[3],
            is_ligand=t[4],
            # Safely handle both 5-tuples and 6-tuples
            master_index=t[5] if len(t) > 5 else None 
        )
    
    def to_tuple(self) -> list:
        """Helper to convert back to the list format for JSON storage."""
        base = [
            self.auth_asym_id,
            self.auth_seq_id,
            self.auth_comp_id,
            self.atom_id,
            self.is_ligand
        ]
        if self.master_index is not None:
            base.append(self.master_index)
        return base

class LigandInteraction(BaseModel):
    """A single interaction between ligand and polymer."""
    type: str
    participants: Tuple[InteractionParticipant, InteractionParticipant]

    @classmethod
    def from_raw(cls, raw: dict) -> "LigandInteraction":
        return cls(
            type=raw["type"],
            participants=(
                InteractionParticipant.from_tuple(raw["participants"][0]),
                InteractionParticipant.from_tuple(raw["participants"][1]),
            ),
        )

class NeighborResidue(BaseModel):
    """A residue in the ligand's neighborhood."""
    auth_asym_id: str
    auth_seq_id: int
    auth_comp_id: str
    master_index: Optional[int] = None # Added

    @classmethod
    def from_tuple(cls, t: list) -> "NeighborResidue":
        return cls(
            auth_asym_id=t[0], 
            auth_seq_id=t[1], 
            auth_comp_id=t[2],
            master_index=t[3] if len(t) > 3 else None
        )
    
    def to_tuple(self) -> list:
        base = [self.auth_asym_id, self.auth_seq_id, self.auth_comp_id]
        if self.master_index is not None:
            base.append(self.master_index)
        return base

class LigandNeighborhood(BaseModel):
    ligand_auth_asym_id: str
    ligand_auth_seq_id: int
    ligand_comp_id: str
    interactions: List[LigandInteraction]
    neighborhood: List[NeighborResidue]

    @classmethod
    def from_raw(cls, raw: Any) -> "LigandNeighborhood":
        # If it's the raw list from Molstar (TSX), take the first entry
        if isinstance(raw, list):
            if not raw:
                raise ValueError("Empty ligand data list.")
            raw = raw[0]
            
        # The schema uses the "ligand" key which is [auth_asym_id, auth_seq_id, comp_id]
        ligand_info = raw["ligand"]
        return cls(
            ligand_auth_asym_id = ligand_info[0],
            ligand_auth_seq_id  = ligand_info[1],
            ligand_comp_id      = ligand_info[2],
            interactions        = [LigandInteraction.from_raw(i) for i in raw["interactions"]],
            neighborhood        = [NeighborResidue.from_tuple(n) for n in raw["neighborhood"]],
        )


# --- Enums ---
class MasterAlignment(BaseModel):
    """Versioned canonical reference MSA for a tubulin family"""
    version      : str
    family       : TubulinFamily
    fasta_content: str
    created_date : str
    description  : Optional[str] = None

class AlignmentMapping(BaseModel):
    seqres_to_master: str  # JSON array: list[int] (seqres_idx -> master_idx | -1)
    master_to_seqres: str  # JSON array: list[int] (master_idx -> seqres_idx | -1)



class MutationType(str, Enum):
    SUBSTITUTION = "substitution"
    INSERTION    = "insertion"
    DELETION     = "deletion"


class ModificationType(str, Enum):
    ACETYLATION     = "acetylation"
    PHOSPHORYLATION = "phosphorylation"
    METHYLATION     = "methylation"
    UBIQUITINATION  = "ubiquitination"
    SUMOYLATION     = "sumoylation"
    PALMITOYLATION  = "palmitoylation"
    NITROSYLATION   = "nitrosylation"
    GLUTAMYLATION   = "glutamylation"
    GLYCYLATION     = "glycylation"
    TYROSINATION    = "tyrosination"
    DETYROSINATION  = "detyrosination"


# --- Support Models ---

class Modification(BaseModel):
    """Post-translational modification from literature/databases"""
    
    master_index: int
    utn_position: Optional[int] = None
    
    amino_acid       : str
    modification_type: str
    
    uniprot_id  : str
    species     : str
    tubulin_type: str
    
    phenotype      : str
    database_source: str
    database_link  : str
    keywords       : str
    notes          : Optional[str] = None

class NonpolymericLigand(BaseModel):
    """Ligand model - unchanged from riboxyz"""
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    
    class NonpolymerComp(BaseModel):
        class Drugbank(BaseModel):
            class DrugbankInfo(BaseModel):
                cas_number: Optional[str] = None
                description: Optional[str] = None

class Mutation(BaseModel):
    master_index: int
    utn_position: Optional[int] = None
    from_residue: str
    to_residue: str

    uniprot_id     : str
    species        : str
    tubulin_type   : TubulinFamily
    phenotype      : str
    database_source: str
    reference_link : str
    keywords       : str
    notes          : Optional[str] = None


class NonpolymerComp(BaseModel):
    class Drugbank(BaseModel):
        class DrugbankInfo(BaseModel):
            cas_number: Optional[str] = None
            description: Optional[str] = None

        class DrugbankContainerIdentifiers(BaseModel):
            drugbank_id: str

        drugbank_container_identifiers: Optional[DrugbankContainerIdentifiers] = None
        drugbank_info: Optional[DrugbankInfo] = None

    class RcsbChemCompTarget(BaseModel):
        interaction_type: Optional[str] = None
        name: Optional[str] = None
        provenance_source: Optional[str] = None
        reference_database_accession_code: Optional[str] = None
        reference_database_name: Optional[str] = None

    drugbank: Optional[Drugbank] = None
    rcsb_chem_comp_target: Optional[List[RcsbChemCompTarget]] = None

# This one is for db only.
class ChemicalCompound(BaseModel):
    """
    Global chemical identity - shared across all structures.
    One instance per unique chemical_id in the entire PDB.
    """
    chemical_id  : str  # Primary key: "TAX", "PMM", "GTP"
    chemical_name: str
    
    SMILES        : Optional[str]   = None
    SMILES_stereo : Optional[str]   = None
    InChI         : Optional[str]   = None
    InChIKey      : Optional[str]   = None
    formula_weight: Optional[float] = None
    
    nonpolymer_comp: Optional[NonpolymerComp] = None



class BaseInstance(BaseModel):
    parent_rcsb_id: str
    auth_asym_id  : str
    asym_id       : str
    entity_id     : str
    assembly_id   : int

    def __hash__(self):
        return hash(self.asym_id + self.parent_rcsb_id)


class BaseEntity(BaseModel):
    entity_id       : str
    type            : Literal["polymer", "non-polymer", "water", "branched"]
    pdbx_description: Optional[str] = None
    formula_weight  : Optional[float] = None
    pdbx_strand_ids : List[str] = []




class NonpolymerEntity(BaseEntity):
    type: Literal["non-polymer"] = "non-polymer"
    
    # Reference to the shared chemical
    chemical_id: str
    chemical_name: str
    
    pdbx_description: Optional[str] = None
    formula_weight: Optional[float] = None
    
    # --- MISSING FIELDS ADDED BELOW ---
    # These are needed so node_ligand.py can create the Global Chemical Node
    nonpolymer_comp: Optional[NonpolymerComp] = None
    
    SMILES        : Optional[str] = None
    SMILES_stereo : Optional[str] = None
    InChI         : Optional[str] = None
    InChIKey      : Optional[str] = None

    num_instances : int = 0

class PolypeptideEntity(BaseEntity):
    type: Literal["polymer"] = "polymer"
    polymer_type: Literal["Protein"] = "Protein"

    one_letter_code: str
    one_letter_code_can: str
    sequence_length: int

    src_organism_names: List[str] = []
    host_organism_names: List[str] = []
    src_organism_ids: List[int] = []
    host_organism_ids: List[int] = []

    family: Optional[PolymerClass] = None  # Changed from Optional[TubulinFamily]
    uniprot_accessions: List[str] = []

    mutations: List[Mutation] = []
    alignment_stats: Dict[str, Any] = {}

    def to_SeqRecord(self, rcsb_id: str) -> SeqRecord:
        return SeqRecord(
            seq=Seq(self.one_letter_code_can),
            id=f"{rcsb_id}_{self.entity_id}",
            description=self.pdbx_description or "",
            name=f"{rcsb_id}_entity_{self.entity_id}",
        )

class PolynucleotideEntity(BaseEntity):

    type               : Literal["polymer"] = "polymer"
    polymer_type       : str                             
    one_letter_code    : str
    one_letter_code_can: str
    sequence_length    : int
    src_organism_names: List[str] = []
    src_organism_ids  : List[int] = []




class Polypeptide(BaseInstance):
    pass


class Polynucleotide(BaseInstance):
    pass


class Nonpolymer(BaseInstance):
    """
    Represents a single nonpolymer instance (one molecule copy) in the structure.
    
    Inherits from BaseInstance:

        - parent_rcsb_id: str - The PDB ID (e.g., "6WVR")
        - auth_asym_id  : str - Author chain ID (e.g., "E")
        - asym_id       : str - Internal chain ID (e.g., "E")
        - entity_id     : str - References the NonpolymerEntity (e.g., "4")
        - assembly_id   : int - Which biological assembly (e.g., 1)
    
    Example:
        Three Taxol molecules in 6WVR would create three Nonpolymer instances,
        all pointing to the same NonpolymerEntity (which points to the same
        ChemicalCompound).
    """
    
    # Currently no additional fields beyond BaseInstance
    # Could add instance-specific data later if needed:
    # occupancy: Optional[float] = None
    # b_factor: Optional[float] = None




# --- 3. STRUCTURE ROOT ---


class AssemblyInstancesMap(BaseModel):
    class InstanceIdentifier(BaseModel):
        entity_id: str
        auth_asym_id: str
        asym_id: Optional[str] = None

    rcsb_id: str
    nonpolymer_entity_instances: Optional[List[Dict[str, InstanceIdentifier]]] = None
    polymer_entity_instances: List[Dict[str, InstanceIdentifier]]


class RCSBStructureMetadata(BaseModel):
    model_config = {"json_encoders": {Enum: lambda v: v.value}}
    rcsb_id: str
    expMethod: str
    resolution: float
    deposition_date: Optional[str] = None

    pdbx_keywords: Optional[str] = None
    pdbx_keywords_text: Optional[str] = None

    rcsb_external_ref_id: List[str]
    rcsb_external_ref_type: List[str]
    rcsb_external_ref_link: List[str]

    citation_year: Optional[int] = None
    citation_rcsb_authors: Optional[List[str]] = None
    citation_title: Optional[str] = None
    citation_pdbx_doi: Optional[str] = None
    
    # --- RESTORED FIELDS ---
    src_organism_ids   : List[int] = []
    src_organism_names : List[str] = []
    host_organism_ids  : List[int] = []
    host_organism_names: List[str] = []


class TubulinStructure(RCSBStructureMetadata):
    entities: Dict[
        str, Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]
    ]

    polypeptides   : List[Polypeptide]
    polynucleotides: List[Polynucleotide]
    nonpolymers    : List[Nonpolymer]

    assembly_map: Optional[List[AssemblyInstancesMap]] = None
    polymerization_state: Optional[
        Literal["monomer", "dimer", "oligomer", "filament", "unknown"]
    ] = None




# ------
# Add to lib/types.py

class MutationEntryData(BaseModel):
    """Mutation entry as stored in sequence ingestion results."""
    ma_position: int
    wild_type: str
    observed: str
    pdb_auth_id: int


class ProcessedChainData(BaseModel):
    """Result of sequence alignment/ingestion for a single entity."""
    pdb_id: str
    chain_id: str
    tubulin_class: str
    sequence: str
    
    # ma_to_auth_map[ma_idx] = auth_seq_id (or -1 if missing)
    ma_to_auth_map: List[int]
    
    # observed_to_ma_map[obs_idx] = MA position (1-based) or -2 for insertions
    observed_to_ma_map: List[int]
    
    mutations: List[MutationEntryData]
    stats: Dict[str, Any]


class SequenceIngestionEntry(BaseModel):
    """A single entity's ingestion record."""
    processed_at: str
    family: str
    data: ProcessedChainData
    
    def build_auth_to_ma_map(self) -> Dict[int, int]:
        """
        Build reverse lookup: auth_seq_id -> MA index (1-based).
        
        The ma_to_auth_map stores: ma_to_auth_map[ma_idx] = auth_seq_id
        We invert this to get: auth_seq_id -> ma_idx + 1
        """
        auth_to_ma: Dict[int, int] = {}
        for ma_idx, auth_id in enumerate(self.data.ma_to_auth_map):
            if auth_id != -1:
                auth_to_ma[auth_id] = ma_idx + 1  # MA positions are 1-based
        return auth_to_ma

```
```
What i want to do now is to expand my structure querying to also be able to filter by the ligands that belong to them.
For that we need to first be able to fetch a list of ligands/chemicals that exist in the db in general and then add a filter 
to structures based on that nonpolymer belonging there.

Let me show you the code.


neo4j_tubxz/db_lib_reader.py
```py
# neo4j_tubxz/db_lib_reader.py
"""
Database reader with typed query methods.
Uses the query builders for filter logic.
"""

import sys
from typing import Optional, List, Dict, Any
from neo4j import ManagedTransaction, Transaction

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

                # All records have the same total_count and next_cursor
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
                        mutation_count=r["mutation_count"],
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
    # Facet/Aggregation Queries (for filter UI dropdowns)
    # -------------------------------------------------------------------------

    def get_filter_facets(self) -> Dict[str, Any]:
        """
        Get available filter options for the UI.
        Returns counts for each facet value.
        """
        query = """
        MATCH (s:Structure)
        WITH collect(s) AS all_structs
        
        // Experimental methods
        UNWIND all_structs AS s
        WITH all_structs, s.expMethod AS method
        WHERE method IS NOT NULL
        WITH all_structs, method, count(*) AS method_count
        WITH all_structs, collect({value: method, count: method_count}) AS exp_methods
        
        // Year range
        UNWIND all_structs AS s
        WITH all_structs, exp_methods, min(s.citation_year) AS min_year, max(s.citation_year) AS max_year
        
        // Resolution range  
        UNWIND all_structs AS s
        WITH exp_methods, min_year, max_year, 
             min(s.resolution) AS min_res, max(s.resolution) AS max_res,
             count(s) AS total_structures
        
        RETURN {
            total_structures: total_structures,
            exp_methods: exp_methods,
            year_range: {min: min_year, max: max_year},
            resolution_range: {min: min_res, max: max_res}
        } AS facets
        """

        with self.adapter.driver.session() as session:

            def run_query(tx: Transaction):
                result = tx.run(query).single()
                return result["facets"] if result else {}

            return session.execute_read(run_query)

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

    # neo4j_tubxz/db_lib_reader.py - add this method to Neo4jReader class

    def get_taxonomy_tree_for_ui(self, tax_type: str = "source") -> List[Dict[str, Any]]:
        """
        Get taxonomy as a tree structure for antd TreeSelect.
        Returns format: { value, title, children }
        """
        rel_type = f"belongs_to_lineage_{tax_type}"
        
        # Get all taxa linked to structures with their parent info
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
                
                # Build flat list with parent refs
                nodes = {}
                for r in records:
                    nodes[r["tax_id"]] = {
                        "value": r["tax_id"],
                        "title": f"{r['name']} ({r['structure_count']})",
                        "rank": r["rank"],
                        "parent_id": r["parent_id"],
                        "children": []
                    }
                
                # Build tree
                roots = []
                for tax_id, node in nodes.items():
                    parent_id = node.pop("parent_id")
                    node.pop("rank")  # Clean up
                    if parent_id and parent_id in nodes:
                        nodes[parent_id]["children"].append(node)
                    else:
                        roots.append(node)
                
                # Remove empty children arrays
                def clean_children(node):
                    if not node["children"]:
                        del node["children"]
                    else:
                        for child in node["children"]:
                            clean_children(child)
                    return node
                
                return [clean_children(r) for r in roots]
            
            return session.execute_read(run_query)

# Singleton instance for convenience
db_reader = Neo4jReader()

```

neo4j_tubxz/structure_query_builder.py
```py
# neo4j_tubxz/structure_query_builder.py
"""
Cypher query builder for Structure filtering.
Designed to be used standalone or via the API layer.
"""
from typing import Dict, Any, Tuple, List, Optional
from neo4j_tubxz.models import StructureFilters


class StructureQueryBuilder:
    """
    Builds Cypher queries for filtering Structures with keyset pagination.
    
    Usage:
        builder = StructureQueryBuilder(filters)
        query, params = builder.build()
        # Execute query with params
    """
    
    def __init__(self, filters: StructureFilters):
        self.filters = filters
        self._where_clauses: List[str] = []
        self._params: Dict[str, Any] = {"limit": filters.limit}
        self._needs_entity_match = False
        self._needs_ligand_match = False
        
    def build(self) -> Tuple[str, Dict[str, Any]]:
        """
        Build the complete query and return (query_string, params).
        """
        self._process_filters()
        return self._assemble(), self._params
    
    def _process_filters(self):
        """Process all filter conditions"""
        f = self.filters
        
        # Text search
        if f.search:
            self._add_search_filter(f.search)
        
        # Exact ID match
        if f.rcsb_ids:
            self._where_clauses.append("s.rcsb_id IN $rcsb_ids")
            self._params["rcsb_ids"] = [x.upper() for x in f.rcsb_ids]
        
        # Resolution range
        if f.resolution_min is not None:
            self._where_clauses.append("s.resolution >= $res_min")
            self._params["res_min"] = f.resolution_min
        if f.resolution_max is not None:
            self._where_clauses.append("s.resolution <= $res_max")
            self._params["res_max"] = f.resolution_max
        
        # Year range
        if f.year_min is not None:
            self._where_clauses.append("s.citation_year >= $year_min")
            self._params["year_min"] = f.year_min
        if f.year_max is not None:
            self._where_clauses.append("s.citation_year <= $year_max")
            self._params["year_max"] = f.year_max
        
        # Experimental method
        if f.exp_method:
            self._where_clauses.append("s.expMethod IN $exp_methods")
            self._params["exp_methods"] = [m.value for m in f.exp_method]
        
        # Polymerization state
        if f.polymerization_state:
            self._where_clauses.append("s.polymerization_state IN $poly_states")
            self._params["poly_states"] = [p.value for p in f.polymerization_state]
        
        # Taxonomy filters (using EXISTS subqueries for performance)
        if f.source_organism_ids:
            self._add_taxonomy_filter(f.source_organism_ids, "source")
        if f.host_organism_ids:
            self._add_taxonomy_filter(f.host_organism_ids, "host")
        
        # Related entity filters
        if f.has_ligand_ids:
            self._add_ligand_filter(f.has_ligand_ids)
        if f.has_polymer_family:
            self._add_polymer_family_filter(f.has_polymer_family)
        if f.has_uniprot:
            self._add_uniprot_filter(f.has_uniprot)
        
        # Cursor (keyset pagination)
        if f.cursor:
            self._where_clauses.append("s.rcsb_id < $cursor")
            self._params["cursor"] = f.cursor.upper()
    
    def _add_search_filter(self, search: str):
        """Full-text search across multiple fields"""
        self._where_clauses.append("""(
            toLower(s.rcsb_id) CONTAINS $search OR
            toLower(s.citation_title) CONTAINS $search OR
            toLower(s.pdbx_keywords_text) CONTAINS $search OR
            ANY(name IN s.src_organism_names WHERE toLower(name) CONTAINS $search)
        )""")
        self._params["search"] = search.lower()
    
    def _add_taxonomy_filter(self, tax_ids: List[int], tax_type: str):
        """
        Filter by taxonomy. Uses belongs_to_lineage relationship which
        includes all descendants of a taxon.
        """
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
        """Filter to structures containing specific ligands"""
        self._where_clauses.append("""
            EXISTS {
                MATCH (s)-[:DEFINES_ENTITY]->(e:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
                WHERE c.chemical_id IN $ligand_ids
            }
        """)
        self._params["ligand_ids"] = [x.upper() for x in chemical_ids]
    

    def _add_polymer_family_filter(self, families: List[str]):
        """
        Filter to structures with ALL specified tubulin families.
        E.g., selecting [alpha, beta] returns structures that have BOTH.
        """
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
        """Filter to structures with specific UniProt accessions"""
        self._where_clauses.append("""
            EXISTS {
                MATCH (s)-[:DEFINES_ENTITY]->(e:PolypeptideEntity)
                WHERE ANY(acc IN e.uniprot_accessions WHERE acc IN $uniprot_ids)
            }
        """)
        self._params["uniprot_ids"] = accessions
    
    def _assemble(self) -> str:
        """Assemble the final query"""
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
    """
    Builds Cypher queries for filtering PolypeptideEntity nodes.
    """
    
    def __init__(self, filters):
        from neo4j_tubxz.models import PolypeptideEntityFilters
        self.filters: PolypeptideEntityFilters = filters
        self._where_clauses: List[str] = []
        self._struct_where_clauses: List[str] = []
        self._params: Dict[str, Any] = {"limit": filters.limit}
    
    def build(self) -> Tuple[str, Dict[str, Any]]:
        self._process_filters()
        return self._assemble(), self._params
    
    def _process_filters(self):
        f = self.filters
        
        # Apply structure-level filters if provided
        if f.structure_filters:
            self._apply_structure_filters(f.structure_filters)
        
        # Entity-specific filters
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
        
        if f.has_mutations is not None:
            if f.has_mutations:
                self._where_clauses.append("EXISTS { MATCH (e)-[:HAS_MUTATION]->(:Mutation) }")
            else:
                self._where_clauses.append("NOT EXISTS { MATCH (e)-[:HAS_MUTATION]->(:Mutation) }")
        
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
    
    def _apply_structure_filters(self, sf):
        """Apply structure-level filters as EXISTS subqueries"""
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
OPTIONAL MATCH (e)-[:HAS_MUTATION]->(m:Mutation)
WITH total_count, next_cursor, e, count(m) AS mutation_count
RETURN 
    e.parent_rcsb_id AS parent_rcsb_id,
    e.entity_id AS entity_id,
    e.pdbx_description AS pdbx_description,
    e.family AS family,
    e.sequence_length AS sequence_length,
    e.src_organism_names AS src_organism_names,
    e.uniprot_accessions AS uniprot_accessions,
    mutation_count,
    total_count,
    next_cursor
        """)
        
        return "\n".join(parts)


class LigandQueryBuilder:
    """Builds Cypher queries for Chemical/Ligand filtering"""
    
    def __init__(self, filters):
        from neo4j_tubxz.models import LigandFilters
        self.filters: LigandFilters = filters
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
```




Can you write that shit or do u want me to show you more???