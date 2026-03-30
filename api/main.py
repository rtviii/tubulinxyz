# api/main.py
import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import RedirectResponse, JSONResponse
from pathlib import Path

from api.routers import (
    router_structures,
    router_polymers,
    router_ligands,
    router_msa,
)

from api.routers.router_annotations import router_annotations
app = FastAPI(
    title="TubXYZ API",
    version="0.2.0",
    description="API for tubulin structure data, spatial grids, and MSA alignment.",
)

cors_origins = os.environ.get("CORS_ALLOWED_ORIGINS", "http://localhost:3000,http://localhost:5173").split(",")
app.add_middleware(
    CORSMiddleware,
    allow_origins=[o.strip() for o in cors_origins],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Routers ---
app.include_router(router_structures, prefix="/structures", tags=["Structures"])
app.include_router(router_polymers, prefix="/polymers", tags=["Polymers"])
app.include_router(router_ligands, prefix="/ligands", tags=["Ligands"])
app.include_router(router_msa, prefix="/msa", tags=["MSA Alignment"])
app.include_router(router_annotations, prefix="/annotations", tags=["Annotations"])
for route in app.routes:
    print(f"  {getattr(route, 'methods', '')} {route.path}")

@app.get("/", include_in_schema=False)
def root():
    return RedirectResponse(url="/docs")


@app.get("/health", tags=["Health"])
def health():
    from neo4j import GraphDatabase
    neo4j_uri  = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
    neo4j_user = os.environ.get("NEO4J_USER", "neo4j")
    neo4j_pw   = os.environ.get("NEO4J_PASSWORD", "")
    neo4j_db   = os.environ.get("NEO4J_CURRENTDB", "neo4j")
    try:
        driver = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_pw))
        with driver.session(database=neo4j_db) as session:
            session.run("RETURN 1").single()
        driver.close()
        return {"status": "ok", "neo4j": "connected", "version": "0.2.0"}
    except Exception as e:
        return JSONResponse(
            status_code=503,
            content={"status": "degraded", "neo4j": str(e), "version": "0.2.0"},
        )


@app.get("/ingest-status", tags=["Health"])
def ingest_status():
    """Returns the last ingestion run status (written by the scheduler service)."""
    import json
    status_file = Path("/var/log/tubxz/last_ingest_status.json")
    if not status_file.exists():
        return {"status": "no ingestion has run yet"}
    try:
        return json.loads(status_file.read_text())
    except Exception:
        return {"status": "error reading status file"}


if __name__ == "__main__":
    import uvicorn
    Path("debug_output").mkdir(exist_ok=True)
    uvicorn.run("api.main:app", host="127.0.0.1", port=8000, reload=True)