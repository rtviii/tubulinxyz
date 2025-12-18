"""
Ligand interaction extraction module.

Runs the TypeScript molstar script to extract binding site interactions
for each ligand instance in a structure and augments them with Master Alignment indices.
"""

import json
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional, Dict, List

from lib.types import LigandNeighborhood, SequenceIngestionEntry

# Ligands to skip (waters, common ions, buffers)
SKIP_LIGANDS = {
    "HOH",
    "DOD",
    "WAT",
    "NA",
    "CL",
    "K",
    "MG",
    "CA",
    "ZN",
    "FE",
    "MN",
    "CO",
    "NI",
    "CU",
    "SO4",
    "PO4",
    "NO3",
    "CO3",
    "GOL",
    "EDO",
    "PEG",
    "PGE",
    "ACT",
    "FMT",
    "MES",
    "TRS",
    "HEP",
    "EPE",
    "CIT",
    "TAR",
    "DMS",
    "DMF",
    "BME",
    "DTT",
}


@dataclass
class ExtractionTarget:
    rcsb_id: str
    comp_id: str
    auth_asym_id: str
    cif_path: Path
    output_path: Path


@dataclass
class ExtractionResult:
    target: ExtractionTarget
    success: bool
    augmented: bool = False
    error: Optional[str] = None


def load_ingestion_data(ingestion_path: Path) -> Dict[str, SequenceIngestionEntry]:
    """Load and validate the sequence ingestion file."""
    if not ingestion_path.exists():
        return {}

    with open(ingestion_path) as f:
        raw = json.load(f)

    return {
        entity_id: SequenceIngestionEntry.model_validate(entry)
        for entity_id, entry in raw.items()
    }


def augment_with_canonical_indices(
    output_path: Path,
    asym_to_entity: Dict[str, str],
    ingestion_data: Dict[str, SequenceIngestionEntry],
) -> bool:
    """
    Reads the generated ligand JSON and adds Master Alignment indices.

    For each polymeric residue in interactions/neighborhood:
      1. Map auth_asym_id (chain) -> entity_id
      2. Get that entity's SequenceIngestionEntry
      3. Use the reverse map (auth_seq_id -> MA index) to tag the residue
    """
    if not output_path.exists():
        return False

    # Pre-build auth_to_ma maps for each entity
    entity_auth_to_ma: Dict[str, Dict[int, int]] = {}
    for entity_id, entry in ingestion_data.items():
        entity_auth_to_ma[entity_id] = entry.build_auth_to_ma_map()

    try:
        with open(output_path) as f:
            raw_data = json.load(f)

        # Handle list wrapper from TypeScript batch output
        if isinstance(raw_data, list):
            if not raw_data:
                return False
            raw_data = raw_data[0]

        neighborhood = LigandNeighborhood.from_raw(raw_data)
        updated_count = 0

        def get_ma_index(auth_asym_id: str, auth_seq_id: int) -> Optional[int]:
            entity_id = asym_to_entity.get(auth_asym_id)
            if entity_id is None:
                return None

            auth_to_ma = entity_auth_to_ma.get(entity_id)
            if auth_to_ma is None:
                return None

            return auth_to_ma.get(auth_seq_id)

        # Augment interactions
        for interaction in neighborhood.interactions:
            for participant in interaction.participants:
                if not participant.is_ligand:
                    ma_idx = get_ma_index(
                        participant.auth_asym_id, participant.auth_seq_id
                    )
                    if ma_idx is not None:
                        participant.master_index = ma_idx
                        updated_count += 1

        # Augment neighborhood residues
        for res in neighborhood.neighborhood:
            ma_idx = get_ma_index(res.auth_asym_id, res.auth_seq_id)
            if ma_idx is not None:
                res.master_index = ma_idx
                updated_count += 1

        if updated_count > 0:
            output = {
                "ligand": [
                    neighborhood.ligand_auth_asym_id,
                    neighborhood.ligand_auth_seq_id,
                    neighborhood.ligand_comp_id,
                ],
                "interactions": [
                    {
                        "type": ix.type,
                        "participants": [p.to_tuple() for p in ix.participants],
                    }
                    for ix in neighborhood.interactions
                ],
                "neighborhood": [res.to_tuple() for res in neighborhood.neighborhood],
            }
            with open(output_path, "w") as f:
                json.dump(output, f, indent=2)
            return True

        return False

    except Exception as e:
        print(f"      Augmentation error for {output_path.name}: {e}")
        return False


def run_single_extraction(
    target: ExtractionTarget,
    script_path: Path,
    project_root: Path,
    asym_to_entity: Optional[Dict[str, str]] = None,
    ingestion_data: Optional[Dict[str, SequenceIngestionEntry]] = None,
) -> ExtractionResult:
    """Run extraction for a single ligand instance and augment with MA indices."""
    cmd = [
        "npx",
        "tsx",
        str(script_path),
        str(target.cif_path),
        target.comp_id,
        target.auth_asym_id,
        str(target.output_path),
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300, cwd=project_root
        )

        if result.returncode != 0:
            return ExtractionResult(
                target=target,
                success=False,
                error=result.stderr[:500] if result.stderr else "Unknown error",
            )

        augmented = False
        if asym_to_entity and ingestion_data:
            augmented = augment_with_canonical_indices(
                target.output_path, asym_to_entity, ingestion_data
            )

        return ExtractionResult(target=target, success=True, augmented=augmented)

    except subprocess.TimeoutExpired:
        return ExtractionResult(target=target, success=False, error="Timeout")
    except Exception as e:
        return ExtractionResult(target=target, success=False, error=str(e))


def extract_ligands_parallel(
    rcsb_id: str,
    cif_path: Path,
    ligand_instances: List[tuple[str, str]],
    output_dir: Path,
    script_path: Path,
    project_root: Path,
    max_workers: int = 8,
    overwrite: bool = False,
    skip_ligands: Optional[set[str]] = None,
    asym_to_entity: Optional[Dict[str, str]] = None,
    ingestion_path: Optional[Path] = None,
) -> List[ExtractionResult]:
    """
    Extract interactions for multiple ligand instances in parallel with augmentation.

    Args:
        ingestion_path: Path to sequence_ingestion.json (replaces raw dict)
    """
    skip = skip_ligands if skip_ligands is not None else SKIP_LIGANDS
    rcsb_id = rcsb_id.upper()

    # Load ingestion data with proper types
    ingestion_data = load_ingestion_data(ingestion_path) if ingestion_path else {}

    targets = []
    for comp_id, auth_asym_id in ligand_instances:
        if comp_id in skip:
            continue
        output_path = output_dir / f"{rcsb_id}_{comp_id}_{auth_asym_id}.json"
        if output_path.exists() and not overwrite:
            continue

        targets.append(
            ExtractionTarget(
                rcsb_id=rcsb_id,
                comp_id=comp_id,
                auth_asym_id=auth_asym_id,
                cif_path=cif_path,
                output_path=output_path,
            )
        )

    if not targets:
        return []

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_single_extraction,
                t,
                script_path,
                project_root,
                asym_to_entity,
                ingestion_data,
            ): t
            for t in targets
        }
        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            status = "ok" if res.success else "FAIL"
            aug_status = "+MA" if res.augmented else ""
            print(
                f"      {res.target.comp_id}_{res.target.auth_asym_id}: {status} {aug_status}"
            )

    return results
