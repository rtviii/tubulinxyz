"""
Bridge to Molstar TypeScript extraction scripts.
"""

import json
import subprocess
from pathlib import Path
from typing import Optional
from loguru import logger

from lib.types import (
    MolstarExtractionResult,
    ObservedSequenceData,
    ObservedResidue,
    LigandBindingSite,
    BindingSiteResidue,
)


def run_molstar_extraction(
    cif_path: Path,
    rcsb_id: str,
    output_path: Path,
    script_path: Path,
    project_root: Path,
    timeout: int = 600,
) -> Optional[MolstarExtractionResult]:
    local_tsx = project_root / "node_modules" / ".bin" / "tsx"

    if local_tsx.exists():
        cmd = [str(local_tsx), str(script_path)]
    else:
        logger.warning("Local tsx not found, falling back to npx")
        cmd = ["npx", "tsx", str(script_path)]

    cmd.extend([str(cif_path), rcsb_id.upper(), str(output_path)])

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=project_root
        )

        if result.returncode != 0:
            logger.error(
                f"Molstar extraction failed for {rcsb_id}: {result.stderr[:500]}"
            )
            return None

        return parse_extraction_result(output_path)

    except subprocess.TimeoutExpired:
        logger.error(f"Molstar extraction timed out for {rcsb_id}")
        return None
    except Exception as e:
        logger.error(f"Molstar extraction error for {rcsb_id}: {e}")
        return None


def parse_extraction_result(output_path: Path) -> Optional[MolstarExtractionResult]:
    if not output_path.exists():
        return None

    try:
        with open(output_path) as f:
            raw = json.load(f)

        sequences = [
            ObservedSequenceData(
                auth_asym_id=seq["auth_asym_id"],
                entity_id=seq["entity_id"],
                residues=[
                    ObservedResidue(
                        auth_seq_id=r["auth_seq_id"],
                        label_seq_id=r["label_seq_id"],
                        comp_id=r["comp_id"],
                        one_letter=r["one_letter"],
                    )
                    for r in seq["residues"]
                ],
            )
            for seq in raw.get("sequences", [])
        ]

        # Parse as LigandBindingSite (renamed from SimplifiedLigandNeighborhood)
        binding_sites = [
            LigandBindingSite(
                ligand_comp_id=lig["ligand_comp_id"],
                ligand_auth_asym_id=lig["ligand_auth_asym_id"],
                ligand_auth_seq_id=lig["ligand_auth_seq_id"],
                residues=[
                    BindingSiteResidue(
                        auth_asym_id=r["auth_asym_id"],
                        observed_index=r["observed_index"],
                        comp_id=r["comp_id"],
                    )
                    for r in lig.get("neighborhood_residues", [])
                ],
            )
            for lig in raw.get("ligand_neighborhoods", [])
        ]

        return MolstarExtractionResult(
            rcsb_id=raw["rcsb_id"],
            sequences=sequences,
            ligand_neighborhoods=binding_sites,
        )

    except Exception as e:
        logger.error(f"Failed to parse extraction result from {output_path}: {e}")
        return None
