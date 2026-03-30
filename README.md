[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/rtviii/tubulinxyz)

# TubXYZ Backend

Structural biology database for tubulin structures. Neo4j graph database, FastAPI backend, ETL pipeline that processes every tubulin-containing structure from the PDB.

## Setup

### Prerequisites

- Python 3.11+
- Node.js (for Molstar-based structure extraction via `tsx`)
- Neo4j 5.x instance (local or remote)
- MUSCLE 3.8.1 binary at the project root (`./muscle3.8.1`)

### Installation

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
npm install          # installs molstar, tsx, jsdom (for CIF parsing)
```

### Environment

Create a `.env` file (or export these variables). All are required unless noted:

```
TUBETL_DATA=/path/to/tubulinxyz/data     # root of the data/ directory where profiles are stored
NEO4J_URI=bolt://localhost:7687
NEO4J_USER=neo4j
NEO4J_PASSWORD=<password>
NEO4J_CURRENTDB=neo4j                    # neo4j database name
NCBI_TAXA_SQLITE=/path/to/taxa.sqlite    # ete3 NCBITaxa sqlite file (for taxonomy resolution)
```

Source it before running anything:

```bash
source .env   # or use python-dotenv which is in requirements
```

## CLI reference

All commands go through `cli.py`:

```bash
python cli.py --help
```

Key commands:

| Command | What it does |
|---|---|
| `collect-one <PDB_ID> [-o]` | Run full ETL for one structure, save profile JSON to disk |
| `collect-missing` | Collect only structures present on RCSB but missing locally |
| `collect-all [-o] [-w 8]` | Collect all tubulin structures from RCSB (`-w` = parallel workers) |
| `upload-one <PDB_ID> [-o]` | Upload one profile from disk into Neo4j |
| `upload-all [-o] [-w 8]` | Upload all local profiles into Neo4j |

The `-o` flag forces overwrite (re-collect or re-upload even if already exists).

### Full data rebuild

```bash
# 1. Regenerate all profiles on disk (fetches from RCSB, runs extraction + alignment)
python cli.py collect-all -o -w 8

# 2. Upload all profiles to Neo4j
python cli.py upload-all -o -w 8

# 3. Ingest Morisette literature annotations (separate from structural ETL)
python -m lib.etl.ingest_morisette --family tubulin_alpha
python -m lib.etl.ingest_morisette --family tubulin_beta
```

## How structures are discovered

The pipeline queries the RCSB Search API for all experimental structures annotated with InterPro `IPR000217` (Tubulin family). This is the parent family encompassing alpha, beta, gamma, delta, epsilon, and zeta tubulins. The query uses `annotation_lineage.id` so it matches any child family annotation. This excludes FtsZ (bacterial tubulin homolog that shares the GTPase superfamily but is not tubulin). The query lives in `lib/etl/assets.py` as `TUBULIN_SEARCH_QUERY` and currently returns ~809 structures.

## ETL pipeline

`TubulinETLCollector.generate_profile()` in `lib/etl/collector.py` runs the following phases for each structure:

### Phase 1: Acquire raw data

- Downloads the `.cif` file from RCSB
- Fetches polymer/nonpolymer entity metadata via RCSB GraphQL (sequences, UniProt accessions, organism info, ligand entities)
- Runs Molstar CIF extraction (`scripts_and_artifacts/extract_structure_data.tsx` via `tsx`) to get observed sequences per chain and ligand binding site neighborhoods from the 3D coordinates

### Phase 2: Classification

Each polypeptide entity's observed sequence is classified into a tubulin family (alpha, beta, gamma, delta, epsilon) or a MAP family using HMM profiles built from curated seed alignments. The classifier lives in `lib/etl/classification.py`, profiles in `data/hmms/`.

### Phase 2.5: Isotype calling

For alpha and beta tubulin entities, resolves the specific human isotype (HGNC symbol like TUBA1A, TUBB3):

- **Tier 1**: Direct lookup of the entity's UniProt accessions against a verified mapping (`data/genenames/uniprot_to_isotype.json`)
- **Tier 2 fallback**: Pairwise sequence alignment against human isotype reference sequences from the master alignment FASTA files, picking the closest match above 85% identity

The isotype caller is in `lib/etl/isotype.py`.

### Phase 3: Sequence alignment

Each tubulin entity's canonical sequence is profile-aligned against the family's master alignment (MSA) using MUSCLE. The MSA files live in `data/{alpha,beta,gamma,delta,epsilon}_tubulin/`. This produces:

- **Entity index mapping**: bidirectional map between the entity's `label_seq_id` (canonical residue positions) and master alignment indices
- **Chain index mappings**: for each chain instance, maps `auth_seq_id` (observed/crystallographic residue numbers) to master alignment indices, accounting for unresolved residues
- **Variants**: substitutions, insertions, and deletions relative to the MSA consensus

This mapping system is what allows the frontend to show residue-level annotations (mutations, modifications, binding contacts) in a unified coordinate system across all structures.

### Phase 4: Binding site augmentation

Takes the raw ligand neighborhood data from Molstar (which residues contact each ligand) and augments it with master alignment indices using the chain mappings from Phase 3. This means binding site residues can be compared across structures in alignment coordinates.

### Phase 5: Assemble and persist

Assembles all the above into a `TubulinStructure` profile and writes it as JSON to `$TUBETL_DATA/<PDB_ID>/profile.json`.

## Uploading to Neo4j

`upload-one` / `upload-all` reads profile JSONs from disk and creates the graph in Neo4j. `Neo4jAdapter.add_total_structure()` in `neo4j_tubxz/db_lib_builder.py` creates:

- `Structure` node with metadata
- `PolypeptideEntity` / `PolynucleotideEntity` / `NonpolymerEntity` nodes linked via `DEFINES_ENTITY`
- Instance nodes (`PolypeptideInstance`, `NonpolymerInstance`) linked via `INSTANCE_OF`
- `Chemical` nodes for ligands, linked via `DEFINED_BY_CHEMICAL`
- `Variant` nodes for sequence variants, linked via `HAS_VARIANT`
- `PhylogenyNode` links for taxonomy
- `NEAR_POLYMER` relationships between ligand and polymer instances

## Morisette annotations

The Morisette et al. database (PLOS ONE 2023) catalogs known tubulin mutations and post-translational modifications from the literature. This is ingested separately from the structural ETL since it's not per-structure data.

The Morisette data uses "Universal Tubulin Numbering" (UTN) -- a 1-based numbering on a gapless consensus of ~90 alpha/beta sequences. Our pipeline needs to translate UTN positions into our master alignment indices.

`lib/etl/utn_mapper.py` handles this: it profile-aligns the UTN consensus sequence (hardcoded from Morisette Supplementary Table S1) against our family MSA using MUSCLE, producing a bidirectional UTN position <-> master alignment index map.

`lib/etl/ingest_morisette.py` then:

1. Reads the mutation and modification CSVs from `data/{alpha,beta}_tubulin/`
2. Parses the AA1 notation (e.g. `M1L` for mutations, `C4-PLM` for modifications)
3. Translates UTN positions to master alignment indices via the mapper
4. Batch-inserts `Variant` (source="morisette") and `Modification` nodes into Neo4j

```bash
python -m lib.etl.ingest_morisette --family tubulin_alpha
python -m lib.etl.ingest_morisette --family tubulin_beta
python -m lib.etl.ingest_morisette --family tubulin_alpha --dry-run   # validate without inserting
```

## API

FastAPI application in `api/main.py`. Start with:

```bash
python -m api.main
# or
uvicorn api.main:app --reload --port 8000
```

OpenAPI docs at `http://localhost:8000/docs`. The frontend (`../fend_tubulinxyz`) auto-generates its API client from this spec via RTK Query codegen.

Routers:

- `/structures` -- list, filter, get structure details
- `/polymers` -- list/filter polypeptide entities (by family, isotype, organism, variants, ligands)
- `/ligands` -- list/filter chemicals
- `/msa` -- master alignment data for visualization
- `/annotations` -- Morisette literature mutations and modifications

## Directory structure

```
lib/
  etl/
    collector.py          # main ETL orchestrator
    classification.py     # HMM-based family classification
    isotype.py            # isotype calling (alpha/beta)
    sequence_alignment.py # MUSCLE profile alignment, variant calling
    molstar_bridge.py     # bridge to tsx extraction script
    augmentation.py       # binding site augmentation
    ingest_morisette.py   # literature annotation ingestion
    utn_mapper.py         # UTN <-> master alignment mapping
    assets.py             # file paths, RCSB search query, GlobalOps
    constants.py          # env vars
  types.py                # all data models (Pydantic)

neo4j_tubxz/
  db_lib_builder.py       # Neo4j write operations (ingestion)
  db_lib_reader.py        # Neo4j read operations (API queries)
  models.py               # API response/filter models
  structure_query_builder.py  # Cypher query builders
  node_*.py               # per-node-type creation functions

api/
  main.py                 # FastAPI app
  routers/                # route handlers
  config.py               # project root, paths

data/
  alpha_tubulin/          # alpha MSA, individual isotype FASTAs, mutation CSVs
  beta_tubulin/           # beta MSA, mutation CSVs
  gamma_tubulin/          # gamma MSA
  delta_tubulin/          # delta MSA
  epsilon_tubulin/        # epsilon MSA
  hmms/                   # HMM profiles for classification
  genenames/              # HGNC gene tables, UniProt-to-isotype mapping
  sequences/              # cross-species multi-sequence alignments
```
