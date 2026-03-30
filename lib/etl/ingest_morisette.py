# lib/etl/ingest_morisette.py
"""
Standalone ingestion script for the Morisette et al. tubulin mutation/modification
database into Neo4j. Reads raw CSVs, translates UTN positions to our master alignment
indices via UTNMapper, and batch-inserts Variant and Modification nodes.

This is separate from the structural ETL pipeline -- Morisette data is not per-structure.

Usage:
    python -m lib.etl.ingest_morisette --family tubulin_alpha
    python -m lib.etl.ingest_morisette --family tubulin_beta
    python -m lib.etl.ingest_morisette --family tubulin_alpha --dry-run
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from loguru import logger

from api.config import settings
from lib.types import SequenceVariant, VariantType, Modification
from lib.etl.utn_mapper import get_mapper, UTNMapper


# ---------------------------------------------------------------------------
# Morisette modification type abbreviations -> our naming
# ---------------------------------------------------------------------------

MORISETTE_MOD_MAP: Dict[str, str] = {
    "ACE": "acetylation",
    "PHO": "phosphorylation",
    "PLM": "palmitoylation",
    "UBI": "ubiquitination",
    "MTL": "methylation",
    "NTR": "nitrosylation",
    "SMO": "sumoylation",
    "GLU": "glutamylation",
    "GLY": "glycylation",
    "TYR": "tyrosination",
}

# ---------------------------------------------------------------------------
# CSV column name normalization (different CSVs use slightly different names)
# ---------------------------------------------------------------------------

# Mutations CSV columns (after skiprows=1):
# UTN, AA1, Tubulin, TUBULIN Link, Species, Phenotype, AA3, Reference, Ref Link Info, Keywords, Notes
#
# Modifications CSV columns (after skiprows=1):
# UTN, AA1, Tubulin, TubLink, Species, Phenotype, AA3, Database, DBLink:, Key words:, Notes:

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names across mutation and modification CSVs."""
    rename = {}
    for col in df.columns:
        low = col.strip().lower().replace(" ", "_").rstrip(":")
        if low in ("tubulin_link", "tublink"):
            rename[col] = "tubulin_link"
        elif low in ("ref_link_info", "dblink"):
            rename[col] = "ref_link"
        elif low in ("keywords", "key_words"):
            rename[col] = "keywords"
        elif low in ("notes"):
            rename[col] = "notes"
        elif low in ("reference", "database"):
            rename[col] = "reference"
        elif low == "tubulin":
            rename[col] = "uniprot_id"
        elif low == "utn":
            rename[col] = "utn_isotype"
        elif low == "aa1":
            rename[col] = "aa1"
        elif low == "species":
            rename[col] = "species"
        elif low == "phenotype":
            rename[col] = "phenotype"
    return df.rename(columns=rename)


# ---------------------------------------------------------------------------
# Parsers for AA1 notation
# ---------------------------------------------------------------------------

def parse_mutation_aa1(aa1: str) -> Optional[Tuple[str, int, str]]:
    """Parse 'M1L' -> (wild_type, utn_position, mutant). Returns None on failure."""
    m = re.match(r"([A-Z])(\d+)([A-Z*])", aa1.strip())
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None


def parse_modification_aa1(aa1: str) -> Optional[Tuple[str, int, str]]:
    """Parse 'C4-PLM' -> (amino_acid, utn_position, mod_type_abbrev). Returns None on failure."""
    m = re.match(r"([A-Za-z])(\d+)-([A-Z]+)", aa1.strip())
    if m:
        return m.group(1).upper(), int(m.group(2)), m.group(3)
    return None


# ---------------------------------------------------------------------------
# CSV -> typed objects
# ---------------------------------------------------------------------------

def parse_mutations_csv(
    csv_path: Path,
    family: str,
    mapper: UTNMapper,
    family_prefix: str,
) -> Tuple[List[SequenceVariant], int, int]:
    """
    Parse a Morisette mutations CSV into SequenceVariant objects.
    Returns (variants, total_rows, skipped_count).
    """
    df = pd.read_csv(csv_path, skiprows=1)
    # Drop fully empty columns (the CSVs have leading empty columns)
    df = df.dropna(axis=1, how="all")
    df = _normalize_columns(df)

    # Filter to rows that have the family prefix in utn_isotype (e.g. "alpha" rows contain "α")
    if "utn_isotype" in df.columns:
        df = df[df["utn_isotype"].str.contains(family_prefix, na=False)]

    variants: List[SequenceVariant] = []
    skipped = 0

    for _, row in df.iterrows():
        aa1 = str(row.get("aa1", "")).strip()
        parsed = parse_mutation_aa1(aa1)
        if not parsed:
            logger.warning(f"Could not parse AA1: '{aa1}'")
            skipped += 1
            continue

        wild_type, utn_pos, mutant = parsed
        master_idx = mapper.utn_to_master(utn_pos)

        if master_idx is None:
            logger.debug(f"UTN position {utn_pos} unmapped to MA (AA1: {aa1})")
            skipped += 1
            continue

        v = SequenceVariant(
            type=VariantType.SUBSTITUTION,
            source="morisette",
            master_index=master_idx,
            wild_type=wild_type,
            observed=mutant,
            uniprot_id=_clean(row.get("uniprot_id")),
            species=_clean(row.get("species")),
            tubulin_type=_clean(row.get("utn_isotype")),
            family=family,
            phenotype=_clean(row.get("phenotype")),
            reference=_clean(row.get("reference")),
            reference_link=_clean(row.get("ref_link")),
            keywords=_clean(row.get("keywords")),
            notes=_clean(row.get("notes")),
            utn_position=utn_pos,
        )
        variants.append(v)

    return variants, len(df), skipped


def parse_modifications_csv(
    csv_path: Path,
    family: str,
    mapper: UTNMapper,
    family_prefix: str,
) -> Tuple[List[Modification], int, int]:
    """
    Parse a Morisette modifications CSV into Modification objects.
    Returns (modifications, total_rows, skipped_count).
    """
    df = pd.read_csv(csv_path, skiprows=1)
    df = df.dropna(axis=1, how="all")
    df = _normalize_columns(df)

    if "utn_isotype" in df.columns:
        df = df[df["utn_isotype"].str.contains(family_prefix, na=False)]

    modifications: List[Modification] = []
    skipped = 0

    for _, row in df.iterrows():
        aa1 = str(row.get("aa1", "")).strip()
        parsed = parse_modification_aa1(aa1)
        if not parsed:
            logger.warning(f"Could not parse modification AA1: '{aa1}'")
            skipped += 1
            continue

        amino_acid, utn_pos, mod_abbrev = parsed
        master_idx = mapper.utn_to_master(utn_pos)

        if master_idx is None:
            logger.debug(f"UTN position {utn_pos} unmapped to MA (AA1: {aa1})")
            skipped += 1
            continue

        mod_type = MORISETTE_MOD_MAP.get(mod_abbrev, mod_abbrev.lower())

        m = Modification(
            master_index=master_idx,
            amino_acid=amino_acid,
            modification_type=mod_type,
            uniprot_id=_clean(row.get("uniprot_id")) or "",
            species=_clean(row.get("species")) or "",
            tubulin_type=_clean(row.get("utn_isotype")) or "",
            family=family,
            utn_position=utn_pos,
            phenotype=_clean(row.get("phenotype")) or "",
            database_source=_clean(row.get("reference")) or "",
            database_link=_clean(row.get("ref_link")) or "",
            keywords=_clean(row.get("keywords")) or "",
            notes=_clean(row.get("notes")),
        )
        modifications.append(m)

    return modifications, len(df), skipped


def _clean(val) -> Optional[str]:
    """Clean a CSV cell value: strip whitespace, convert nan to None."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return None
    s = str(val).strip()
    return s if s and s.lower() != "nan" else None


# ---------------------------------------------------------------------------
# Family config: where the CSVs live
# ---------------------------------------------------------------------------

_FAMILY_CSV_CONFIG = {
    "tubulin_alpha": {
        "mutations_csv": "alpha_tubulin/alpha_tubulin_mutations.csv",
        "modifications_csv": "alpha_tubulin/alpha_tubulin_modifications.csv",
        "prefix": "\u03b1",  # α
    },
    "tubulin_beta": {
        "mutations_csv": "beta_tubulin/beta_tubulin_mutations.csv",
        "modifications_csv": "beta_tubulin/beta_tubulin_modifications.csv",
        "prefix": "\u03b2",  # β
    },
}


# ---------------------------------------------------------------------------
# Main ingestion
# ---------------------------------------------------------------------------

def ingest_family(family: str, dry_run: bool = False):
    """Run full Morisette ingestion for a family."""

    if family not in _FAMILY_CSV_CONFIG:
        logger.error(f"No CSV config for family '{family}'. Available: {list(_FAMILY_CSV_CONFIG.keys())}")
        sys.exit(1)

    cfg = _FAMILY_CSV_CONFIG[family]
    data_dir = settings.PROJECT_ROOT / "data"

    mutations_csv = data_dir / cfg["mutations_csv"]
    modifications_csv = data_dir / cfg["modifications_csv"]

    # Build UTN -> MA mapping
    logger.info(f"Building UTN mapper for {family}...")
    mapper = get_mapper(family)

    # Parse mutations
    variants: List[SequenceVariant] = []
    if mutations_csv.exists():
        logger.info(f"Parsing mutations from {mutations_csv.name}...")
        v, total, skip = parse_mutations_csv(mutations_csv, family, mapper, cfg["prefix"])
        variants = v
        logger.info(f"  {len(v)} variants from {total} rows ({skip} skipped)")
    else:
        logger.warning(f"Mutations CSV not found: {mutations_csv}")

    # Parse modifications
    modifications: List[Modification] = []
    if modifications_csv.exists():
        logger.info(f"Parsing modifications from {modifications_csv.name}...")
        m, total, skip = parse_modifications_csv(modifications_csv, family, mapper, cfg["prefix"])
        modifications = m
        logger.info(f"  {len(m)} modifications from {total} rows ({skip} skipped)")
    else:
        logger.warning(f"Modifications CSV not found: {modifications_csv}")

    if dry_run:
        logger.info("[DRY RUN] Would insert:")
        logger.info(f"  {len(variants)} literature variants")
        logger.info(f"  {len(modifications)} modifications")
        return

    # Insert into Neo4j
    from neo4j import GraphDatabase
    from neo4j_tubxz.node_variant import batch_create_literature_variants
    from neo4j_tubxz.node_modification import batch_create_modifications

    # Use same connection params as the main adapter
    import os
    uri = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
    user = os.environ.get("NEO4J_USER", "neo4j")
    password = os.environ.get("NEO4J_PASSWORD", "neo4j")
    db = os.environ.get("NEO4J_DATABASE", "neo4j")

    driver = GraphDatabase.driver(uri, auth=(user, password), database=db)

    try:
        with driver.session() as session:
            if variants:
                logger.info(f"Inserting {len(variants)} literature variants...")
                count = session.execute_write(
                    lambda tx: batch_create_literature_variants(tx, variants)
                )
                logger.info(f"  Created/merged {count} Variant nodes")

            if modifications:
                logger.info(f"Inserting {len(modifications)} modifications...")
                count = session.execute_write(
                    lambda tx: batch_create_modifications(tx, modifications)
                )
                logger.info(f"  Created/merged {count} Modification nodes")

        logger.info("Done.")
    finally:
        driver.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Ingest Morisette mutation/modification data into Neo4j")
    parser.add_argument("--family", required=True, choices=list(_FAMILY_CSV_CONFIG.keys()),
                        help="Tubulin family to ingest")
    parser.add_argument("--dry-run", action="store_true",
                        help="Parse and validate without inserting into Neo4j")
    args = parser.parse_args()

    ingest_family(args.family, dry_run=args.dry_run)
