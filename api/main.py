# api/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import RedirectResponse
from pathlib import Path

from api.routers import (
    router_structures,
    router_polymers,
    router_ligands,
    router_msa,
    router_grid,
)

from api.routers.router_annotations import router_annotations
app = FastAPI(
    title="TubXYZ API",
    version="0.2.0",
    description="API for tubulin structure data, spatial grids, and MSA alignment.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:5173", "*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Routers ---
app.include_router(router_structures, prefix="/structures", tags=["Structures"])
app.include_router(router_polymers, prefix="/polymers", tags=["Polymers"])
app.include_router(router_ligands, prefix="/ligands", tags=["Ligands"])
app.include_router(router_msa, prefix="/msa", tags=["MSA Alignment"])
app.include_router(router_grid, prefix="/grid", tags=["Grid"])
app.include_router(router_annotations, prefix="/annotations", tags=["Annotations"])

@app.get("/", include_in_schema=False)
def root():
    return RedirectResponse(url="/docs")


@app.get("/health", tags=["Health"])
def health():
    return {"status": "ok", "version": "0.2.0"}


if __name__ == "__main__":
    import uvicorn
    Path("debug_output").mkdir(exist_ok=True)
    uvicorn.run("api.main:app", host="127.0.0.1", port=8000, reload=True)