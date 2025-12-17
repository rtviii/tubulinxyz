"""
Ligand interaction extraction module.

Downloads CIF structures and runs the TypeScript extraction script
to generate interaction reports for each ligand instance.
"""

import json
import subprocess
import requests
from pathlib import Path
from typing import Optional

from lib.etl.constants import TUBETL_DATA
from lib.etl.assets import TubulinStructureAssets
from lib.models.types_tubulin import TubulinStructure, LigandNeighborhood


# Ligands to skip (waters, common ions, etc.)
SKIP_LIGANDS = {"HOH", "DOD", "NA", "CL", "K", "MG", "CA", "ZN", "SO4", "PO4", "GOL", "EDO", "ACT"}

# Path to the extraction script (adjust as needed)
SCRIPT_PATH = Path(__file__).parent.parent.parent / "scripts_and_artifacts" / "extract_ixs.tsx"


def download_cif(rcsb_id: str, output_path: Path) -> bool:
    """Download CIF file from RCSB PDB."""
    url = f"https://files.rcsb.org/download/{rcsb_id.upper()}.cif"
    
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        output_path.write_text(resp.text)
        print(f"Downloaded: {output_path}")
        return True
    except requests.RequestException as e:
        print(f"Failed to download {rcsb_id}: {e}")
        return False


def run_extraction(
    cif_path: Path,
    ligand_comp_id: str,
    auth_asym_id: str,
    output_path: Path,
) -> bool:
    """Run the TypeScript extraction script."""
    cmd = [
        "npx", "tsx",
        str(SCRIPT_PATH),
        str(cif_path),
        ligand_comp_id,
        auth_asym_id,
        str(output_path),
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )
        
        if result.returncode != 0:
            print(f"Extraction failed for {ligand_comp_id}/{auth_asym_id}:")
            print(result.stderr)
            return False
            
        return True
        
    except subprocess.TimeoutExpired:
        print(f"Extraction timed out for {ligand_comp_id}/{auth_asym_id}")
        return False
    except Exception as e:
        print(f"Extraction error: {e}")
        return False


def extract_ligand_interactions(
    rcsb_id: str,
    overwrite: bool = False,
    skip_ligands: Optional[set] = None,
) -> list[Path]:
    """
    Extract interactions for all ligand instances in a structure.
    
    Args:
        rcsb_id: PDB ID
        overwrite: Re-run extraction even if output exists
        skip_ligands: Set of comp_ids to skip (defaults to SKIP_LIGANDS)
    
    Returns:
        List of paths to generated report files
    """
    rcsb_id = rcsb_id.upper()
    skip = skip_ligands if skip_ligands is not None else SKIP_LIGANDS
    
    assets = TubulinStructureAssets(rcsb_id)
    assets._verify_dir_exists()
    
    # Ensure CIF exists
    cif_path = Path(assets.paths.cif)
    if not cif_path.exists():
        if not download_cif(rcsb_id, cif_path):
            raise RuntimeError(f"Could not obtain CIF for {rcsb_id}")
    
    # Load profile to get ligand instances
    profile = assets.profile()
    
    generated_files = []
    
    # Group instances by (comp_id, auth_asym_id)
    for nonpoly in profile.nonpolymers:
        entity = profile.entities.get(nonpoly.entity_id)
        if entity is None:
            continue
            
        comp_id = entity.chemical_id
        
        if comp_id in skip:
            continue
        
        auth_asym_id = nonpoly.auth_asym_id
        output_name = f"{rcsb_id}_{comp_id}_{auth_asym_id}.json"
        output_path = Path(assets.paths.base_dir) / output_name
        
        if output_path.exists() and not overwrite:
            print(f"Skipping (exists): {output_name}")
            generated_files.append(output_path)
            continue
        
        print(f"Extracting: {rcsb_id} / {comp_id} / {auth_asym_id}")
        
        if run_extraction(cif_path, comp_id, auth_asym_id, output_path):
            generated_files.append(output_path)
    
    return generated_files


def load_neighborhood_report(path: Path) -> list[LigandNeighborhood]:
    """Load and parse a neighborhood report JSON file."""
    with open(path) as f:
        raw = json.load(f)
    return [LigandNeighborhood.from_raw(item) for item in raw]


# --- CLI Entry Point ---

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python -m lib.etl.ligand_extraction <RCSB_ID> [--overwrite]")
        sys.exit(1)
    
    rcsb_id = sys.argv[1]
    overwrite = "--overwrite" in sys.argv
    
    try:
        files = extract_ligand_interactions(rcsb_id, overwrite=overwrite)
        print(f"\nGenerated {len(files)} report(s):")
        for f in files:
            print(f"  {f}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)