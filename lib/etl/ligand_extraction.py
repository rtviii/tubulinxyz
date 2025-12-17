"""
Ligand interaction extraction module.

Runs the TypeScript molstar script to extract binding site interactions
for each ligand instance in a structure.
"""

import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional

from lib.etl.constants import TUBETL_DATA

# Ligands to skip (waters, common ions, buffers)
SKIP_LIGANDS = {
    "HOH", "DOD",  # water
    "NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN", "CO", "NI", "CU",  # ions
    "SO4", "PO4", "NO3",  # common anions
    "GOL", "EDO", "PEG", "PGE",  # polyols
    "ACT", "FMT",  # acetate, formate
    "MES", "TRS", "HEP", "EPE",  # buffers
    "CIT", "TAR",  # acids
    "DMS", "DMF",  # solvents
    "BME", "DTT",  # reducing agents
}


@dataclass
class ExtractionTarget:
    """A single ligand instance to extract."""
    rcsb_id: str
    comp_id: str
    auth_asym_id: str
    cif_path: Path
    output_path: Path


@dataclass 
class ExtractionResult:
    """Result of a single extraction."""
    target: ExtractionTarget
    success: bool
    error: Optional[str] = None


def run_single_extraction(
    target: ExtractionTarget,
    script_path: Path,
    project_root: Path,
) -> ExtractionResult:
    """Run extraction for a single ligand instance."""
    cmd = [
        "npx", "tsx",
        str(script_path),
        str(target.cif_path),
        target.comp_id,
        target.auth_asym_id,
        str(target.output_path),
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=project_root,
        )
        
        if result.returncode != 0:
            return ExtractionResult(
                target=target,
                success=False,
                error=result.stderr[:500] if result.stderr else "Unknown error",
            )
        return ExtractionResult(target=target, success=True)
        
    except subprocess.TimeoutExpired:
        return ExtractionResult(target=target, success=False, error="Timeout")
    except Exception as e:
        return ExtractionResult(target=target, success=False, error=str(e))


def extract_ligands_parallel(
    rcsb_id: str,
    cif_path: Path,
    ligand_instances: list[tuple[str, str]],  # [(comp_id, auth_asym_id), ...]
    output_dir: Path,
    script_path: Path,
    project_root: Path,
    max_workers: int = 4,
    overwrite: bool = False,
    skip_ligands: Optional[set[str]] = None,
) -> list[ExtractionResult]:
    """
    Extract interactions for multiple ligand instances in parallel.
    
    Args:
        rcsb_id: PDB ID
        cif_path: Path to the CIF file
        ligand_instances: List of (comp_id, auth_asym_id) tuples
        output_dir: Directory to write output files
        script_path: Path to extract_ixs.tsx
        project_root: Project root for subprocess cwd
        max_workers: Number of parallel workers
        overwrite: Re-extract even if output exists
        skip_ligands: Set of comp_ids to skip
    
    Returns:
        List of ExtractionResult objects
    """
    skip = skip_ligands if skip_ligands is not None else SKIP_LIGANDS
    rcsb_id = rcsb_id.upper()
    
    # Build targets
    targets = []
    for comp_id, auth_asym_id in ligand_instances:
        if comp_id in skip:
            continue
        
        output_path = output_dir / f"{rcsb_id}_{comp_id}_{auth_asym_id}.json"
        
        if output_path.exists() and not overwrite:
            continue
        
        targets.append(ExtractionTarget(
            rcsb_id=rcsb_id,
            comp_id=comp_id,
            auth_asym_id=auth_asym_id,
            cif_path=cif_path,
            output_path=output_path,
        ))
    
    if not targets:
        return []
    
    print(f"  Extracting {len(targets)} ligand(s) with {max_workers} workers...")
    
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_single_extraction, target, script_path, project_root
            ): target
            for target in targets
        }
        
        for future in as_completed(futures):
            result = future.result()
            status = "ok" if result.success else f"FAILED: {result.error}"
            print(f"    {result.target.comp_id}/{result.target.auth_asym_id}: {status}")
            results.append(result)
    
    return results


# Convenience function for standalone use
def extract_for_structure(
    rcsb_id: str,
    max_workers: int = 4,
    overwrite: bool = False,
) -> list[ExtractionResult]:
    """
    Extract all ligand interactions for a structure.
    
    Loads the profile to get ligand instances, then runs parallel extraction.
    """
    from api.config import PROJECT_ROOT
    from lib.etl.assets import TubulinStructureAssets
    
    rcsb_id = rcsb_id.upper()
    assets = TubulinStructureAssets(rcsb_id)
    
    cif_path = Path(assets.paths.cif)
    if not cif_path.exists():
        raise FileNotFoundError(f"CIF not found: {cif_path}")
    
    profile = assets.profile()
    
    # Gather (comp_id, auth_asym_id) for each nonpolymer instance
    ligand_instances = []
    for nonpoly in profile.nonpolymers:
        entity = profile.entities.get(nonpoly.entity_id)
        if entity is None:
            continue
        ligand_instances.append((entity.chemical_id, nonpoly.auth_asym_id))
    
    return extract_ligands_parallel(
        rcsb_id=rcsb_id,
        cif_path=cif_path,
        ligand_instances=ligand_instances,
        output_dir=Path(assets.paths.base_dir),
        script_path=PROJECT_ROOT / "scripts_and_artifacts" / "extract_ixs.tsx",
        project_root=PROJECT_ROOT,
        max_workers=max_workers,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python -m lib.etl.ligand_extraction <RCSB_ID> [--overwrite] [--workers=N]")
        sys.exit(1)
    
    rcsb_id = sys.argv[1]
    overwrite = "--overwrite" in sys.argv
    
    workers = 4
    for arg in sys.argv:
        if arg.startswith("--workers="):
            workers = int(arg.split("=")[1])
    
    results = extract_for_structure(rcsb_id, max_workers=workers, overwrite=overwrite)
    
    succeeded = sum(1 for r in results if r.success)
    print(f"\nDone: {succeeded}/{len(results)} succeeded")