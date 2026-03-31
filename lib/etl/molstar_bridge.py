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


RENDER_SCRIPT = Path(__file__).resolve().parent.parent.parent / "scripts_and_artifacts" / "render_thumbnail.tsx"


def run_thumbnail_render(
    cif_path: Path,
    rcsb_id: str,
    profile_path: Path,
    output_path: Path,
    project_root: Path,
    timeout: int = 120,
) -> bool:
    """Render a structure thumbnail using headless Molstar. Returns True on success."""
    local_tsx = project_root / "node_modules" / ".bin" / "tsx"

    if local_tsx.exists():
        cmd = [str(local_tsx), str(RENDER_SCRIPT)]
    else:
        logger.warning("Local tsx not found, falling back to npx")
        cmd = ["npx", "tsx", str(RENDER_SCRIPT)]

    cmd.extend([str(cif_path), rcsb_id.upper(), str(profile_path), str(output_path)])

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=project_root
        )

        if result.returncode != 0:
            logger.error(f"Thumbnail render failed for {rcsb_id}: {result.stderr[:500]}")
            return False

        logger.debug(f"Thumbnail rendered: {output_path}")
        return True

    except subprocess.TimeoutExpired:
        logger.error(f"Thumbnail render timed out for {rcsb_id}")
        return False
    except Exception as e:
        logger.error(f"Thumbnail render error for {rcsb_id}: {e}")
        return False


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

        # Parse ligand neighborhoods into LigandBindingSite objects.
        # Handles both old format (observed_index) and new format (auth_seq_id)
        # via the BindingSiteResidue validation alias.
        binding_sites = []
        for lig in raw.get("ligand_neighborhoods", []):
            residues = []
            for r in lig.get("neighborhood_residues", []):
                # Support both field names during transition
                seq_id = r.get("auth_seq_id", r.get("observed_index"))
                residues.append(
                    BindingSiteResidue(
                        auth_asym_id=r["auth_asym_id"],
                        auth_seq_id=seq_id,
                        comp_id=r["comp_id"],
                    )
                )
            binding_sites.append(
                LigandBindingSite(
                    ligand_comp_id=lig["ligand_comp_id"],
                    ligand_auth_asym_id=lig["ligand_auth_asym_id"],
                    ligand_auth_seq_id=lig["ligand_auth_seq_id"],
                    residues=residues,
                )
            )

        return MolstarExtractionResult(
            rcsb_id=raw["rcsb_id"],
            sequences=sequences,
            ligand_neighborhoods=binding_sites,
        )

    except Exception as e:
        logger.error(f"Failed to parse extraction result from {output_path}: {e}")
        return None
