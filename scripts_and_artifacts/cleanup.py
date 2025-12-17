#!/usr/bin/env python3
"""
Cleanup Script:
1. Deletes clustering artifacts (*_nr.fasta, *.clstr) from data/sequences.
2. Renames *_final.fasta files to match the strict Enum values (removing '_final').
   e.g. 'map_atat1_final.fasta' -> 'map_atat1.fasta'
"""

import os
from pathlib import Path
from loguru import logger

# Import your Enums to validate filenames against the schema
from lib.models.types_tubulin import TubulinFamily, MapFamily

# Define base paths
PROJECT_ROOT = Path(__file__).parent.parent
SEQ_DIR = PROJECT_ROOT / "data" / "sequences"
MAPS_DIR = SEQ_DIR / "maps"
TUBULIN_DIR = SEQ_DIR / "tubulin"

def clean_and_rename(directory: Path, family_enum):
    """
    Process a specific directory: remove artifacts, rename finals, validate against Enum.
    """
    if not directory.exists():
        logger.warning(f"Directory not found: {directory}")
        return

    logger.info(f"Processing directory: {directory}")

    # 1. Remove artifacts
    # We look for files ending in _nr.fasta (clusters) or .clstr (CD-HIT output)
    artifacts = list(directory.glob("*_nr.fasta")) + list(directory.glob("*.clstr"))
    
    for artifact in artifacts:
        try:
            os.remove(artifact)
            logger.debug(f"Deleted artifact: {artifact.name}")
        except OSError as e:
            logger.error(f"Error deleting {artifact.name}: {e}")

    # 2. Rename _final.fasta files
    final_files = list(directory.glob("*_final.fasta"))
    
    for source_path in final_files:
        # Determine strict name: map_atat1_final.fasta -> map_atat1.fasta
        old_name = source_path.name
        new_name = old_name.replace("_final.fasta", ".fasta")
        target_path = directory / new_name
        
        # Validation: Ensure the resulting stem matches an Enum value
        stem = new_name.replace(".fasta", "")
        
        # specific check because your TubulinFamily enum values are "tubulin_alpha" 
        # but files might be "alpha_final.fasta" depending on how previous script ran.
        # Based on your tree, they are "tubulin_alpha_final.fasta", so strictly matching enum is fine.
        
        is_valid = any(member.value == stem for member in family_enum)
        
        if not is_valid:
            logger.warning(f"⚠️  Filename {new_name} does not match any member in {family_enum.__name__}. Skipping rename.")
            continue

        try:
            source_path.rename(target_path)
            logger.success(f"Renamed: {old_name} -> {new_name}")
        except OSError as e:
            logger.error(f"Error renaming {old_name}: {e}")

def main():
    logger.info("Starting sequence directory cleanup...")
    
    # Process MAPs
    clean_and_rename(MAPS_DIR, MapFamily)
    
    # Process Tubulins
    clean_and_rename(TUBULIN_DIR, TubulinFamily)
    
    logger.info("Cleanup complete. Sequence folder is now synced with Enums.")

if __name__ == "__main__":
    main()