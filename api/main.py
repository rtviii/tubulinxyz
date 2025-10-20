from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
import sys
import traceback
import os

from fastapi.responses import PlainTextResponse

# Add the parent directory to the path to import tubulin_analyzer
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

from tubulin_analyzer import SpatialGridGenerator, GridData, DebugData

# Import the alignment mapper
from api_muscle_alignment import TubulinAlignmentMapper, Annotation

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

# --- Initialize Services ---
grid_generator = SpatialGridGenerator()

# Initialize alignment mapper with your specific paths
MUSCLE_BINARY = str(parent_dir / "muscle3.8.1")
MASTER_PROFILE = str(parent_dir / "data" / "alpha_tubulin"/"alpha_tubulin.afasta")  # or alpha_tubulin.fasta

# Create the mapper instance
alignment_mapper = TubulinAlignmentMapper(
    master_profile_path=MASTER_PROFILE,
    muscle_binary=MUSCLE_BINARY
)

# --- API Endpoints ---

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

# --- NEW ALIGNMENT ENDPOINTS ---

@app.post("/align/sequence")
async def align_sequence(sequence: str, sequence_id: str = None, annotations: list = None):
    """
    Align a tubulin sequence to the master alpha-tubulin profile
    
    Args:
        sequence: Amino acid sequence to align
        sequence_id: Optional identifier for the sequence
        annotations: Optional list of annotations in format [{"start": 1, "end": 10, "label": "PTM"}]
    """
    try:
        if sequence_id is None:
            sequence_id = f"user_sequence_{hash(sequence) % 10000:04d}"
        
        if annotations is None:
            annotations = []
        
        # Convert to Annotation objects
        annotation_objs = [
            Annotation(
                start=ann['start'],
                end=ann['end'],
                label=ann['label'],
                metadata=ann.get('metadata', {})
            ) for ann in annotations
        ]
        
        # Perform alignment
        aligned_sequence, mapping = alignment_mapper.align_sequence(sequence_id, sequence)
        mapped_annotations = alignment_mapper.map_annotations(annotation_objs, mapping)
        statistics = alignment_mapper.get_alignment_statistics(aligned_sequence, mapping)
        
        return {
            "sequence_id": sequence_id,
            "aligned_sequence": aligned_sequence,
            "mapped_annotations": [
                {
                    "start": ann.start,
                    "end": ann.end,
                    "label": ann.label,
                    "metadata": ann.metadata
                } for ann in mapped_annotations
            ],
            "statistics": statistics,
            "mapping": mapping,  # For debugging
            "original_sequence": sequence
        }
        
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Alignment failed: {str(e)}")

@app.get("/align/master-profile")
async def get_master_profile():
    """Get information about the master alignment profile INCLUDING THE ACTUAL SEQUENCES"""
    try:
        profile_path = Path(MASTER_PROFILE)
        if not profile_path.exists():
            raise HTTPException(status_code=404, detail="Master profile not found")
        
        # Read and parse the master profile
        from Bio import AlignIO
        alignment = AlignIO.read(profile_path, "fasta")
        
        sequences_info = []
        for record in alignment:
            sequences_info.append({
                "id": record.id,
                "description": record.description,
                "length": len(record.seq),
                "gap_count": str(record.seq).count('-'),
                "sequence": str(record.seq)  # <- HERE'S THE FUCKING SEQUENCE!
            })
        
        # Also return the full alignment as a string for good measure
        full_alignment = ""
        for record in alignment:
            full_alignment += f">{record.description}\n{record.seq}\n"
        
        return {
            "profile_path": str(profile_path),
            "profile_exists": profile_path.exists(),
            "num_sequences": len(alignment),
            "alignment_length": alignment.get_alignment_length(),
            "sequences": sequences_info,
            "full_alignment": full_alignment,  # <- AND THE WHOLE DAMN THING!
            "muscle_binary": MUSCLE_BINARY,
            "muscle_exists": Path(MUSCLE_BINARY).exists()
        }
        
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Error reading master profile: {str(e)}")

@app.post("/align/structure/{pdb_id}/{chain_id}")
async def align_structure_chain(pdb_id: str, chain_id: str, annotations: list = None):
    """
    Extract and align a specific chain from a PDB structure
    """
    try:
        # This would need your PDB parsing logic - here's a skeleton
        pdb_sequence = await extract_chain_sequence(pdb_id, chain_id)
        
        if pdb_sequence is None:
            raise HTTPException(status_code=404, detail=f"Chain {chain_id} not found in {pdb_id}")
        
        # Use the alignment endpoint
        return await align_sequence(pdb_sequence, f"{pdb_id}_{chain_id}", annotations)
        
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Structure alignment failed: {str(e)}")

async def extract_chain_sequence(pdb_id: str, chain_id: str) -> str:
    """
    Extract sequence from PDB chain - you'll need to implement this
    based on your existing PDB processing code
    """
    # TODO: Integrate with your existing PDB processing code
    # This is a placeholder - you'll want to use your actual PDB parsing
    try:
        # Example using your existing structure processing
        # You might need to adapt this based on your actual code
        from tubulin_analyzer.structural_analyzer import load_structure
        structure = await load_structure(pdb_id)
        
        # Extract chain sequence logic here
        # This is pseudocode - adapt to your actual implementation
        for chain in structure.get_chains():
            if chain.id == chain_id:
                return extract_aa_sequence(chain)
        
        return None
    except Exception as e:
        print(f"Error extracting sequence from {pdb_id} chain {chain_id}: {e}")
        return None

# --- Existing endpoints (keep all your current functionality) ---

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

@app.get("/models/{filename}")
async def get_model_file(filename: str):
    """Serve mmCIF model files with computed metadata"""
    try:
        print(f"Current working directory: {os.getcwd()}")
        
        # Try different path strategies
        paths_to_try = [
            f"/Users/rtviii/dev/tubulinxyz/api/maxim_data/{filename}",
            f"maxim_data/{filename}",
            f"./maxim_data/{filename}",
            f"../maxim_data/{filename}",
        ]
        
        print(f"Looking for file: {filename}")
        for path in paths_to_try:
            print(f"Trying path: {path}")
            if os.path.exists(path):
                print(f"✅ Found file at: {path}")
                with open(path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                return PlainTextResponse(
                    content=content,
                    media_type="chemical/x-mmcif",
                    headers={
                        "Content-Disposition": f"inline; filename={filename}",
                        "Access-Control-Allow-Origin": "*"
                    }
                )
            else:
                print(f"❌ Not found at: {path}")
        
        # List what's actually in the directory
        try:
            maxim_dir = "/Users/rtviii/dev/tubulinxyz/api/maxim_data"
            if os.path.exists(maxim_dir):
                files = os.listdir(maxim_dir)
                print(f"Files in {maxim_dir}: {files}")
            else:
                print(f"Directory {maxim_dir} doesn't exist")
        except Exception as e:
            print(f"Error listing directory: {e}")
        
        raise HTTPException(status_code=404, detail=f"File '{filename}' not found in any location")
        
    except HTTPException:
        raise
    except Exception as e:
        print(f"Error serving model file {filename}: {e}")

@app.get("/models")
async def list_models():
    """List available model files"""
    try:
        models_dir = Path("maxim_data")
        if not models_dir.exists():
            models_dir.mkdir(exist_ok=True)
            return {"models": [], "message": "maxim_data directory created but empty"}
        
        # Find all .cif files
        cif_files = list(models_dir.glob("*.cif"))
        
        models = []
        for file_path in cif_files:
            stat = file_path.stat()
            models.append({
                "filename": file_path.name,
                "size_bytes": stat.st_size,
                "modified": stat.st_mtime
            })
        
        return {
            "models": models,
            "count": len(models),
            "directory": str(models_dir.absolute())
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error listing models: {str(e)}")

# Update root endpoint to include alignment endpoints
@app.get("/")
async def root():
    """API information and usage"""
    return {
        "title": "Tubulin Spatial Grid API",
        "version": "0.1.0", 
        "description": "Converts microtubule PDB structures into idealized 2D grids with sequence alignment",
        "endpoints": {
            "/grid/{pdb_id}": "Generate 2D grid layout (cached automatically)",
            "/debug/{pdb_id}": "Debug N-terminus connectivity analysis", 
            "/debug/{pdb_id}/pymol": "Get PyMOL visualization instructions",
            "/align/sequence": "POST - Align a tubulin sequence to master profile",
            "/align/master-profile": "GET - Info about master alignment",
            "/align/structure/{pdb_id}/{chain_id}": "POST - Align a PDB chain",
            "/models": "List available model files in maxim_data/",
            "/models/{filename}": "Serve mmCIF model files from maxim_data/"
        },
        "alignment": {
            "master_profile": MASTER_PROFILE,
            "method": "MUSCLE profile alignment",
            "version": "3.8.1"
        },
        "example_usage": [
            "curl http://localhost:8000/grid/6o2t",
            "curl http://localhost:8000/align/master-profile",
            'curl -X POST http://localhost:8000/align/sequence -H "Content-Type: application/json" -d \'{"sequence": "MREIVHIQ...", "annotations": [{"start": 40, "end": 40, "label": "K40_acetylation"}]]}\'',
            "curl http://localhost:8000/models"
        ]
    }

@app.get("/health")
async def health_check():
    """Health check endpoint including alignment service"""
    debug_dir = Path("debug_output")
    models_dir = Path("maxim_data")
    profile_path = Path(MASTER_PROFILE)
    muscle_path = Path(MUSCLE_BINARY)
    
    return {
        "status": "healthy",
        "version": "0.1.0",
        "debug_directory": str(debug_dir),
        "debug_directory_exists": debug_dir.exists(),
        "models_directory": str(models_dir),
        "models_directory_exists": models_dir.exists(),
        "alignment_service": {
            "master_profile": str(profile_path),
            "master_profile_exists": profile_path.exists(),
            "muscle_binary": str(muscle_path),
            "muscle_exists": muscle_path.exists(),
            "service_ready": profile_path.exists() and muscle_path.exists()
        }
    }

if __name__ == "__main__":
    import uvicorn
    
    # Ensure debug output directory exists
    Path("debug_output").mkdir(exist_ok=True)
    
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)