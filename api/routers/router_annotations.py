# api/routers/router_annotations.py (NEW FILE)
from fastapi import APIRouter, HTTPException, Query
from typing import List, Dict, Any, Optional
from neo4j_tubxz.db_lib_reader import dbqueries

router_annotations = APIRouter()

@router_annotations.get("/mutations/{family}/{version}/{position}")
async def get_mutations_at_position(
    family: str,
    version: str,
    position: int
) -> Dict[str, Any]:
    """
    Get all mutations at a specific master alignment position.
    
    - **family**: Tubulin family (e.g., "alpha", "beta", "gamma")
    - **version**: Master alignment version (e.g., "v1.0")
    - **position**: Master alignment position (1-based index)
    """
    try:
        mutations = dbqueries.get_mutations_at_position(position, family, version)
        
        return {
            "position": position,
            "family": family,
            "version": version,
            "count": len(mutations),
            "mutations": mutations
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/modifications/{family}/{version}/{position}")
async def get_modifications_at_position(
    family: str,
    version: str, 
    position: int
) -> Dict[str, Any]:
    """
    Get all modifications at a specific master alignment position.
    
    - **family**: Tubulin family (e.g., "alpha", "beta", "gamma")
    - **version**: Master alignment version (e.g., "v1.0")
    - **position**: Master alignment position (1-based index)
    """
    try:
        modifications = dbqueries.get_modifications_at_position(position, family, version)
        
        return {
            "position": position,
            "family": family,
            "version": version,
            "count": len(modifications),
            "modifications": modifications
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/all/{family}/{version}/{position}")
async def get_all_annotations_at_position(
    family: str,
    version: str,
    position: int
) -> Dict[str, Any]:
    """
    Get both mutations and modifications at a specific master alignment position.
    
    - **family**: Tubulin family (e.g., "alpha", "beta", "gamma")
    - **version**: Master alignment version (e.g., "v1.0")
    - **position**: Master alignment position (1-based index)
    """
    try:
        data = dbqueries.get_all_annotations_at_position(position, family, version)
        
        return {
            **data,
            "total_count": len(data["mutations"]) + len(data["modifications"])
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/range/{family}/{version}")
async def get_annotations_in_range(
    family: str,
    version: str,
    start: int = Query(..., description="Start position (inclusive)"),
    end: int = Query(..., description="End position (inclusive)")
) -> Dict[str, Any]:
    """
    Get all annotations within a position range.
    """
    try:
        annotations_by_position = {}
        
        for pos in range(start, end + 1):
            data = dbqueries.get_all_annotations_at_position(pos, family, version)
            if data["mutations"] or data["modifications"]:
                annotations_by_position[pos] = data
        
        return {
            "family": family,
            "version": version,
            "range": {"start": start, "end": end},
            "positions_with_annotations": len(annotations_by_position),
            "data": annotations_by_position
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/mutations")
async def get_polymer_mutations(
    rcsb_id: str,
    auth_asym_id: str
) -> Dict[str, Any]:
    """
    Get all mutations for a specific polymer.
    
    - **rcsb_id**: PDB ID (e.g., "1SA0")
    - **auth_asym_id**: Chain ID (e.g., "A")
    """
    try:
        mutations = dbqueries.get_mutations_for_polymer(rcsb_id, auth_asym_id)
        
        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "count": len(mutations),
            "mutations": mutations
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/modifications")
async def get_polymer_modifications(
    rcsb_id: str,
    auth_asym_id: str
) -> Dict[str, Any]:
    """
    Get all modifications for a specific polymer.
    
    - **rcsb_id**: PDB ID (e.g., "1SA0")
    - **auth_asym_id**: Chain ID (e.g., "A")
    """
    try:
        modifications = dbqueries.get_modifications_for_polymer(rcsb_id, auth_asym_id)
        
        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "count": len(modifications),
            "modifications": modifications
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/all")
async def get_polymer_all_annotations(
    rcsb_id: str,
    auth_asym_id: str
) -> Dict[str, Any]:
    """
    Get both mutations and modifications for a specific polymer.
    
    - **rcsb_id**: PDB ID (e.g., "1SA0")
    - **auth_asym_id**: Chain ID (e.g., "A")
    """
    try:
        data = dbqueries.get_all_annotations_for_polymer(rcsb_id, auth_asym_id)
        return data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))