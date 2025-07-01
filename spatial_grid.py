from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel
from typing import List, Tuple, Optional, Dict, Any
import numpy as np
import asyncio
import httpx
from pathlib import Path
import tempfile
import logging
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure
import math
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)

# Pydantic models for the response
class SubunitData(BaseModel):
    id: str
    auth_asym_id: str
    protofilament: int
    subunitIndex: int
    monomerType: str  # 'α' or 'β'

class GridData(BaseModel):
    subunits: List[SubunitData]
    structure_type: str  # 'dimer', 'linear', 'sheet', 'cylindrical'
    metadata: Dict[str, Any]

router = APIRouter()

class SpatialGridGenerator:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
    
    async def download_mmcif(self, pdb_id: str) -> Path:
        """Download mmCIF file from PDB"""
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        
        async with httpx.AsyncClient() as client:
            response = await client.get(url)
            response.raise_for_status()
            
            # Save to temporary file
            temp_file = Path(tempfile.mktemp(suffix=".cif"))
            temp_file.write_bytes(response.content)
            return temp_file
    
    def extract_chain_positions(self, structure: Structure) -> Dict[str, Tuple[float, float, float]]:
        """Extract first residue position for each chain"""
        chain_positions = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                
                # Get first residue
                residues = list(chain.get_residues())
                if not residues:
                    continue
                    
                first_residue = residues[0]
                
                # Get first atom of first residue (usually backbone atom)
                atoms = list(first_residue.get_atoms())
                if not atoms:
                    continue
                    
                first_atom = atoms[0]
                coord = first_atom.get_coord()
                chain_positions[chain_id] = tuple(coord)
                
        return chain_positions
    
    def detect_structure_type(self, positions: Dict[str, Tuple[float, float, float]]) -> str:
        """Detect if structure is dimer, linear, sheet, or cylindrical"""
        coords = np.array(list(positions.values()))
        n_chains = len(coords)
        
        if n_chains <= 2:
            return "dimer"
        elif n_chains <= 4:
            return "linear"
        
        # Calculate distances and analyze geometry
        distances = pdist(coords)
        mean_dist = np.mean(distances)
        
        # Check for cylindrical arrangement
        if self._is_cylindrical(coords):
            return "cylindrical"
        elif self._is_sheet_like(coords):
            return "sheet"
        else:
            return "linear"
    
    def _is_cylindrical(self, coords: np.ndarray) -> bool:
        """Check if points are arranged in a roughly cylindrical pattern"""
        if len(coords) < 8:  # Need enough points to detect cylinder
            return False
            
        # Find the main axis using PCA
        pca = PCA(n_components=3)
        pca.fit(coords)
        
        # Project onto the first two principal components (radial plane)
        coords_2d = pca.transform(coords)[:, :2]
        
        # Calculate distances from origin in 2D
        radial_distances = np.linalg.norm(coords_2d, axis=1)
        
        # Check if radial distances are roughly consistent (cylindrical)
        radial_std = np.std(radial_distances)
        radial_mean = np.mean(radial_distances)
        
        # If coefficient of variation is low, likely cylindrical
        cv = radial_std / radial_mean if radial_mean > 0 else float('inf')
        
        return cv < 0.3 and len(coords) > 12  # Threshold for cylindrical detection
    
    def _is_sheet_like(self, coords: np.ndarray) -> bool:
        """Check if points are arranged in a roughly planar sheet"""
        if len(coords) < 4:
            return False
            
        # Use PCA to find the best-fit plane
        pca = PCA(n_components=3)
        pca.fit(coords)
        
        # Check if most variance is in first two components (planar)
        explained_variance = pca.explained_variance_ratio_
        
        # If third component explains little variance, it's roughly planar
        return explained_variance[2] < 0.2
    
    def map_dimer(self, positions: Dict[str, Tuple[float, float, float]]) -> List[SubunitData]:
        """Map dimer structure to grid"""
        chain_ids = list(positions.keys())
        subunits = []
        
        for i, chain_id in enumerate(chain_ids):
            subunits.append(SubunitData(
                id=f"pf{i}-{chain_id}",
                auth_asym_id=chain_id,
                protofilament=i,
                subunitIndex=0,
                monomerType='α' if i % 2 == 0 else 'β'
            ))
            
        return subunits
    
    def map_linear(self, positions: Dict[str, Tuple[float, float, float]]) -> List[SubunitData]:
        """Map linear filament to grid"""
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        
        # Find the main axis using PCA
        pca = PCA(n_components=1)
        pca.fit(coords)
        
        # Project coordinates onto main axis
        projected = pca.transform(coords).flatten()
        
        # Sort by position along axis
        sorted_indices = np.argsort(projected)
        
        subunits = []
        for i, idx in enumerate(sorted_indices):
            chain_id = chain_ids[idx]
            subunits.append(SubunitData(
                id=f"pf0-{chain_id}",
                auth_asym_id=chain_id,
                protofilament=0,
                subunitIndex=i,
                monomerType='α' if i % 2 == 0 else 'β'
            ))
            
        return subunits
    
    def map_sheet(self, positions: Dict[str, Tuple[float, float, float]]) -> List[SubunitData]:
        """Map sheet structure to grid"""
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        
        # Use PCA to find the sheet plane
        pca = PCA(n_components=2)
        pca.fit(coords)
        coords_2d = pca.transform(coords)
        
        # Cluster into protofilaments using DBSCAN or simple grid detection
        # For now, use a simple approach: group by X coordinate ranges
        x_coords = coords_2d[:, 0]
        y_coords = coords_2d[:, 1]
        
        # Find approximate grid structure
        x_sorted_indices = np.argsort(x_coords)
        
        # Group into protofilaments by clustering X coordinates
        x_clusters = self._cluster_coordinates(x_coords)
        
        subunits = []
        for chain_idx, chain_id in enumerate(chain_ids):
            pf = x_clusters[chain_idx]
            
            # Within each protofilament, sort by Y coordinate
            pf_indices = [i for i, cluster in enumerate(x_clusters) if cluster == pf]
            pf_y_coords = [y_coords[i] for i in pf_indices]
            pf_chain_ids = [chain_ids[i] for i in pf_indices]
            
            # Sort by Y within protofilament
            y_sorted = sorted(enumerate(pf_y_coords), key=lambda x: x[1])
            
            # Find position of current chain in this protofilament
            chain_pos_in_pf = next(i for i, (orig_i, _) in enumerate(y_sorted) 
                                 if pf_chain_ids[orig_i] == chain_id)
            
            subunits.append(SubunitData(
                id=f"pf{pf}-{chain_id}",
                auth_asym_id=chain_id,
                protofilament=pf,
                subunitIndex=chain_pos_in_pf,
                monomerType='α' if (pf + chain_pos_in_pf) % 2 == 0 else 'β'
            ))
            
        return subunits
    
    def map_cylindrical(self, positions: Dict[str, Tuple[float, float, float]]) -> List[SubunitData]:
        """Map cylindrical structure (microtubule) to grid"""
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        
        # Find cylinder axis using PCA
        pca = PCA(n_components=3)
        pca.fit(coords)
        
        # Assume main axis is the first principal component
        axis = pca.components_[0]
        
        # Project coordinates onto plane perpendicular to axis
        coords_centered = coords - np.mean(coords, axis=0)
        
        # Project onto the plane (remove component along axis)
        projected = coords_centered - np.outer(np.dot(coords_centered, axis), axis)
        
        # Convert to cylindrical coordinates
        x, y = projected[:, 0], projected[:, 1]
        angles = np.arctan2(y, x)
        
        # Normalize angles to [0, 2π]
        angles = (angles + 2 * np.pi) % (2 * np.pi)
        
        # Z coordinates (along cylinder axis)
        z_coords = np.dot(coords_centered, axis)
        
        # Group into protofilaments by angle
        angle_clusters = self._cluster_angles(angles)
        
        subunits = []
        for chain_idx, chain_id in enumerate(chain_ids):
            pf = angle_clusters[chain_idx]
            
            # Within each protofilament, sort by Z coordinate
            pf_indices = [i for i, cluster in enumerate(angle_clusters) if cluster == pf]
            pf_z_coords = [z_coords[i] for i in pf_indices]
            pf_chain_ids = [chain_ids[i] for i in pf_indices]
            
            # Sort by Z within protofilament
            z_sorted = sorted(enumerate(pf_z_coords), key=lambda x: x[1])
            
            # Find position of current chain in this protofilament
            chain_pos_in_pf = next(i for i, (orig_i, _) in enumerate(z_sorted) 
                                 if pf_chain_ids[orig_i] == chain_id)
            
            subunits.append(SubunitData(
                id=f"pf{pf}-{chain_id}",
                auth_asym_id=chain_id,
                protofilament=pf,
                subunitIndex=chain_pos_in_pf,
                monomerType='α' if (pf + chain_pos_in_pf) % 2 == 0 else 'β'
            ))
            
        return subunits
    
    def _cluster_coordinates(self, coords: np.ndarray, threshold: float = None) -> List[int]:
        """Cluster 1D coordinates into groups"""
        if threshold is None:
            # Auto-determine threshold based on data spread
            threshold = np.std(coords) / 2
            
        coords_reshaped = coords.reshape(-1, 1)
        clustering = DBSCAN(eps=threshold, min_samples=1).fit(coords_reshaped)
        return clustering.labels_.tolist()
    
    def _cluster_angles(self, angles: np.ndarray) -> List[int]:
        """Cluster angles into protofilaments"""
        # Sort angles and detect gaps to determine protofilament boundaries
        sorted_indices = np.argsort(angles)
        sorted_angles = angles[sorted_indices]
        
        # Find gaps in angle distribution
        angle_diffs = np.diff(sorted_angles)
        
        # Add wraparound difference
        wraparound_diff = (sorted_angles[0] + 2*np.pi) - sorted_angles[-1]
        angle_diffs = np.append(angle_diffs, wraparound_diff)
        
        # Determine cluster threshold
        threshold = np.mean(angle_diffs) + np.std(angle_diffs)
        
        # Assign clusters
        clusters = np.zeros(len(angles), dtype=int)
        current_cluster = 0
        
        for i in range(len(sorted_indices)):
            clusters[sorted_indices[i]] = current_cluster
            if i < len(angle_diffs) and angle_diffs[i] > threshold:
                current_cluster += 1
                
        return clusters.tolist()
    
    async def generate_grid(self, pdb_id: str) -> GridData:
        """Main method to generate grid data from PDB structure"""
        try:
            # Download and parse mmCIF
            mmcif_path = await self.download_mmcif(pdb_id)
            structure = self.parser.get_structure(pdb_id, mmcif_path)
            
            # Extract chain positions
            positions = self.extract_chain_positions(structure)
            
            if not positions:
                raise HTTPException(status_code=400, detail="No chains found in structure")
            
            # Detect structure type
            structure_type = self.detect_structure_type(positions)
            
            # Map to grid based on structure type
            if structure_type == "dimer":
                subunits = self.map_dimer(positions)
            elif structure_type == "linear":
                subunits = self.map_linear(positions)
            elif structure_type == "sheet":
                subunits = self.map_sheet(positions)
            elif structure_type == "cylindrical":
                subunits = self.map_cylindrical(positions)
            else:
                raise HTTPException(status_code=400, detail=f"Unknown structure type: {structure_type}")
            
            # Clean up temporary file
            mmcif_path.unlink()
            
            return GridData(
                subunits=subunits,
                structure_type=structure_type,
                metadata={
                    "pdb_id": pdb_id,
                    "num_chains": len(positions),
                    "chain_ids": list(positions.keys())
                }
            )
            
        except Exception as e:
            logger.error(f"Error generating grid for {pdb_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Failed to generate grid: {str(e)}")

# Initialize generator
grid_generator = SpatialGridGenerator()

@router.post("/structures/{pdb_id}/generate_grid", response_model=GridData)
async def generate_spatial_grid(pdb_id: str):
    """
    Generate 2D grid layout from 3D molecular structure
    
    Args:
        pdb_id: PDB identifier (e.g., "6dpv")
        
    Returns:
        GridData with spatial grid layout for frontend
    """
    return await grid_generator.generate_grid(pdb_id)

@router.get("/structures/{pdb_id}/chain_info")
async def get_chain_info(pdb_id: str):
    """Debug endpoint to see raw chain positions"""
    try:
        mmcif_path = await grid_generator.download_mmcif(pdb_id)
        structure = grid_generator.parser.get_structure(pdb_id, mmcif_path)
        positions = grid_generator.extract_chain_positions(structure)
        structure_type = grid_generator.detect_structure_type(positions)
        
        mmcif_path.unlink()
        
        return {
            "pdb_id": pdb_id,
            "structure_type": structure_type,
            "chain_positions": positions,
            "num_chains": len(positions)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))