# api/routers/router_grid.py
from fastapi import APIRouter, HTTPException
from pathlib import Path
import json
import traceback

from lib.tubulin_analyzer import SpatialGridGenerator, GridData
from api.config import settings

router_grid = APIRouter()

grid_generator = SpatialGridGenerator()


@router_grid.get("/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """Generate idealized 2D grid from PDB structure with caching."""
    try:
        cache_file = settings.CACHE_DIR / f"{pdb_id.lower()}_grid.json"

        if cache_file.exists():
            print(f"Loading cached grid data for {pdb_id.upper()}")
            with open(cache_file, "r") as f:
                cached_data = json.load(f)
                return GridData(**cached_data)

        print(f"Generating new grid data for {pdb_id.upper()}")
        grid_data = await grid_generator.generate_grid(pdb_id.lower())

        with open(cache_file, "w") as f:
            json.dump(grid_data.dict(), f, indent=2)

        return grid_data

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Internal error: {str(e)}")