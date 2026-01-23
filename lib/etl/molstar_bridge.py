"""
Bridge to Molstar TypeScript extraction scripts.

Handles subprocess calls and result parsing.
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
    LigandNeighborhood,
    LigandInteraction,
    NeighborResidue,
    InteractionParticipant,
)


def run_molstar_extraction(
    cif_path: Path,
    rcsb_id: str,
    output_path: Path,
    script_path: Path,
    project_root: Path,
    timeout: int = 600,
) -> Optional[MolstarExtractionResult]:
    """
    Run the unified Molstar extraction script.
    
    Returns MolstarExtractionResult on success, None on failure.
    """
    # Use local tsx binary to ensure it sees local node_modules
    local_tsx = project_root / "node_modules" / ".bin" / "tsx"
    
    if local_tsx.exists():
        cmd = [str(local_tsx), str(script_path)]
    else:
        # Fallback to npx, but this may fail if molstar isn't globally installed
        logger.warning("Local tsx not found, falling back to npx (may fail)")
        cmd = ["npx", "tsx", str(script_path)]
    
    cmd.extend([
        str(cif_path),
        rcsb_id.upper(),
        str(output_path),
    ])
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=project_root,
        )
        
        if result.returncode != 0:
            logger.error(f"Molstar extraction failed for {rcsb_id}: {result.stderr[:500]}")
            return None
        
        return parse_extraction_result(output_path)
        
    except subprocess.TimeoutExpired:
        logger.error(f"Molstar extraction timed out for {rcsb_id}")
        return None
    except Exception as e:
        logger.error(f"Molstar extraction error for {rcsb_id}: {e}")
        return None


def parse_extraction_result(output_path: Path) -> Optional[MolstarExtractionResult]:
    """Parse the JSON output from Molstar extraction."""
    if not output_path.exists():
        return None
    
    try:
        with open(output_path) as f:
            raw = json.load(f)
        
        # Parse sequences
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
                ]
            )
            for seq in raw.get("sequences", [])
        ]
        
        # Parse ligand neighborhoods
        neighborhoods = []
        for lig in raw.get("ligand_neighborhoods", []):
            ligand_info = lig["ligand"]
            
            interactions = []
            for ix in lig.get("interactions", []):
                interactions.append(
                    LigandInteraction(
                        type=ix["type"],
                        participants=(
                            InteractionParticipant.from_tuple(ix["participants"][0]),
                            InteractionParticipant.from_tuple(ix["participants"][1]),
                        )
                    )
                )
            
            neighborhood_residues = [
                NeighborResidue(
                    auth_asym_id=n[0],
                    auth_seq_id=n[1],
                    auth_comp_id=n[2],
                )
                for n in lig.get("neighborhood", [])
            ]
            
            neighborhoods.append(
                LigandNeighborhood(
                    ligand_auth_asym_id=ligand_info[0],
                    ligand_auth_seq_id=ligand_info[1],
                    ligand_comp_id=ligand_info[2],
                    interactions=interactions,
                    neighborhood=neighborhood_residues,
                )
            )
        
        return MolstarExtractionResult(
            rcsb_id=raw["rcsb_id"],
            sequences=sequences,
            ligand_neighborhoods=neighborhoods,
        )
        
    except Exception as e:
        logger.error(f"Failed to parse extraction result from {output_path}: {e}")
        return None