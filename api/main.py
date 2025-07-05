from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
import sys
import traceback
import os

# Add the parent directory to the path to import tubulin_analyzer
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

from tubulin_analyzer import SpatialGridGenerator, GridData, DebugData

# --- FastAPI App ---
app = FastAPI(
    title="Tubulin Spatial Grid API", 
    version="0.1.0",
    description="Generates an idealized, aligned 2D grid from 3D tubulin structures using N-terminus connectivity analysis."
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- API Endpoints ---
grid_generator = SpatialGridGenerator()

@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """Generate idealized 2D grid from PDB structure with caching"""
    try:
        # Check cache first
        cache_dir = Path("cache")
        cache_dir.mkdir(exist_ok=True)
        cache_file = cache_dir / f"{pdb_id.lower()}_grid.json"
        
        if cache_file.exists():
            print(f"Loading cached grid data for {pdb_id.upper()}")
            with open(cache_file, 'r') as f:
                import json
                cached_data = json.load(f)
                # Reconstruct GridData object from cached JSON
                grid_data = GridData(**cached_data)
                return grid_data
        
        # Generate new grid data
        print(f"Generating new grid data for {pdb_id.upper()}")
        grid_data = await grid_generator.generate_grid(pdb_id.lower())
        
        # Save to cache
        with open(cache_file, 'w') as f:
            import json
            json.dump(grid_data.dict(), f, indent=2)
        
        print(f"Grid data cached to: {cache_file}")
        return grid_data
        
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"An internal error occurred: {str(e)}")

@app.get("/debug/{pdb_id}", response_model=DebugData)
async def debug_analysis(pdb_id: str):
    """Debug N-terminus connectivity and protofilament tracing"""
    try:
        # Run the grid generation which includes debugging
        grid_data = await grid_generator.generate_grid(pdb_id.lower())
        
        # Extract debug info from generated files
        debug_dir = Path("debug_output")
        debug_files = []
        
        # Check for generated debug files
        potential_files = [
            "nterm_debug.json",
            f"{pdb_id}_visualize.pml", 
            f"view_{pdb_id}.sh",
            f"protofilaments_{pdb_id}.png",
            f"grid_{pdb_id}_{grid_data.structure_type}.png"
        ]
        
        for file in potential_files:
            if (debug_dir / file).exists():
                debug_files.append(file)
        
        # Calculate debug metrics
        pf_lengths = []
        for subunit in grid_data.subunits:
            pf_idx = subunit.protofilament
            while len(pf_lengths) <= pf_idx:
                pf_lengths.append(0)
            pf_lengths[pf_idx] += 1
        
        # Remove empty protofilaments
        pf_lengths = [length for length in pf_lengths if length > 0]
        
        return DebugData(
            pdb_id=pdb_id.upper(),
            num_tubulin_chains=grid_data.metadata.get("num_tubulin_chains", 0),
            num_connections=grid_data.metadata.get("nterm_connections", 0),
            num_protofilaments=grid_data.metadata.get("num_protofilaments", 0),
            protofilament_lengths=pf_lengths,
            connection_success_rate=(
                grid_data.metadata.get("nterm_connections", 0) / 
                max(grid_data.metadata.get("num_tubulin_chains", 1), 1)
            ),
            debug_files=debug_files
        )
        
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Debug analysis failed: {str(e)}")

@app.get("/debug/{pdb_id}/pymol")
async def get_pymol_instructions(pdb_id: str):
    """Get PyMOL visualization instructions"""
    debug_dir = Path("debug_output")
    
    instructions = {
        "pdb_id": pdb_id.upper(),
        "steps": [
            f"1. Download PDB file: wget https://files.rcsb.org/download/{pdb_id.upper()}.cif",
            f"2. Copy the .cif file to the debug_output/ directory",
            f"3. Run the batch script: ./debug_output/view_{pdb_id}.sh",
            f"   OR manually: pymol {pdb_id.upper()}.cif -d \"@debug_output/{pdb_id}_visualize.pml\""
        ],
        "files_needed": [
            f"{pdb_id.upper()}.cif (structure file)",
            f"{pdb_id}_visualize.pml (PyMOL script)"
        ],
        "visualization_shows": [
            "N-terminus residues colored green",
            "Connected chain pairs colored red", 
            "Distance measurements between connected N-termini",
            "Overall tubulin structure in cartoon representation"
        ]
    }
    
    # Check if files exist
    script_file = debug_dir / f"{pdb_id}_visualize.pml"
    batch_file = debug_dir / f"view_{pdb_id}.sh"
    
    if not script_file.exists():
        instructions["error"] = f"PyMOL script not found. Run /debug/{pdb_id} first to generate visualization files."
    
    instructions["files_status"] = {
        "pymol_script": script_file.exists(),
        "batch_script": batch_file.exists()
    }
    
    return instructions

@app.get("/")
async def root():
    """API information and usage"""
    return {
        "title": "Tubulin Spatial Grid API",
        "version": "0.1.0", 
        "description": "Converts microtubule PDB structures into idealized 2D grids",
        "endpoints": {
            "/grid/{pdb_id}": "Generate 2D grid layout (cached automatically)",
            "/debug/{pdb_id}": "Debug N-terminus connectivity analysis",
            "/debug/{pdb_id}/pymol": "Get PyMOL visualization instructions"
        },
        "algorithm": {
            "method": "N-terminus connectivity tracing",
            "cutoff": "7.0 Ångström for N-terminus neighbors",
            "min_contacts": "3 atoms for valid connection",
            "output": "Protofilament traces → 2D grid coordinates"
        },
        "example_usage": [
            "curl http://localhost:8000/grid/6o2t",
            "curl http://localhost:8000/debug/6o2t"
        ]
    }

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    debug_dir = Path("debug_output")
    return {
        "status": "healthy",
        "version": "0.1.0",
        "debug_directory": str(debug_dir),
        "debug_directory_exists": debug_dir.exists()
    }

if __name__ == "__main__":
    import uvicorn
    
    # Ensure debug output directory exists
    Path("debug_output").mkdir(exist_ok=True)
    
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)