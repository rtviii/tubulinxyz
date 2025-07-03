from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Any
import uvicorn
from structural_analyzer import SpatialGridGenerator

# --- Pydantic Models ---
class SubunitData(BaseModel):
    id: str
    auth_asym_id: str
    protofilament: int
    subunitIndex: int
    monomerType: str

class GridData(BaseModel):
    subunits: List[SubunitData]
    structure_type: str
    metadata: Dict[str, Any]

# --- FastAPI App ---
app = FastAPI(title="Tubulin Spatial Grid API", version="4.2.0",
              description="Generates an idealized, aligned 2D grid from 3D tubulin structures using interface analysis.")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- API Endpoint ---
grid_generator = SpatialGridGenerator()

@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    try:
        return await grid_generator.generate_grid(pdb_id.lower())
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"An internal error occurred: {str(e)}")

if __name__ == "__main__":
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)
