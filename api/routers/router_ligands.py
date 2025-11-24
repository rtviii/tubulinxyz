# api/routers/router_ligands.py
from fastapi import APIRouter

router_ligands = APIRouter()

@router_ligands.get("/all")
def list_all_ligands():
    """Get all ligands"""
    # Add this query to db_lib_reader if needed
    return {"message": "TODO: implement ligand listing"}