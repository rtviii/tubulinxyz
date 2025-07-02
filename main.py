from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Any, Tuple
import numpy as np
import asyncio
import httpx
from pathlib import Path
import tempfile
import uvicorn
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
import math
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
app = FastAPI(title="Tubulin Spatial Grid API", version="2.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Core Logic: SpatialGridGenerator ---
class SpatialGridGenerator:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
        self.profile_base_path = Path(os.getenv("PROFILE_BASE_PATH", "profiles"))
        self.profile_base_path.mkdir(exist_ok=True)

    async def get_profile(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch PDB profile from RCSB API or load from local cache."""
        profile_path = self.profile_base_path / f"{pdb_id.upper()}_profile.json"
        if profile_path.exists():
            print(f"✅ Loading profile from cache: {profile_path}")
            with open(profile_path, 'r') as f:
                return json.load(f)

        print(f"📥 Fetching profile from RCSB for {pdb_id.upper()}")
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        async with httpx.AsyncClient() as client:
            try:
                response = await client.get(url)
                response.raise_for_status()
                profile_data = response.json()
                with open(profile_path, 'w') as f:
                    json.dump(profile_data, f)
                return profile_data
            except httpx.HTTPStatusError as e:
                print(f"🔥 Failed to fetch profile for {pdb_id}: {e}")
                raise HTTPException(status_code=404, detail=f"PDB profile for {pdb_id} not found.")

    def filter_tubulin_chains(self, profile: Dict[str, Any]) -> Dict[str, str]:
        """Extract alpha/beta tubulin chain mappings from the profile."""
        chain_types = {}
        entry_data = profile.get("entry", {})
        if not entry_data or "polymer_entities" not in entry_data:
            print("⚠️ 'polymer_entities' key not found in profile['entry'].")
            return chain_types

        for entity in entry_data.get("polymer_entities", []):
            description = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "").lower()
            is_alpha = "alpha" in description and "tubulin" in description
            is_beta = "beta" in description and "tubulin" in description

            if is_alpha or is_beta:
                monomer_type = "alpha" if is_alpha else "beta"
                container = entity.get("rcsb_polymer_entity_container_identifiers", {})
                auth_asym_ids = container.get("auth_asym_ids", [])
                for chain_id in auth_asym_ids:
                    chain_types[chain_id] = monomer_type

        print(f"🧬 Found tubulin chains: {len(chain_types)} total")
        if chain_types:
            alpha_count = sum(1 for t in chain_types.values() if t == "alpha")
            beta_count = sum(1 for t in chain_types.values() if t == "beta")
            print(f"   α-tubulin: {alpha_count}, β-tubulin: {beta_count}")
        return chain_types


    async def download_mmcif(self, pdb_id: str) -> Path:
        """Download mmCIF file from RCSB."""
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        async with httpx.AsyncClient(timeout=60.0) as client:
            try:
                response = await client.get(url)
                response.raise_for_status()
                temp_file = Path(tempfile.mktemp(suffix=".cif"))
                temp_file.write_bytes(response.content)
                return temp_file
            except httpx.HTTPStatusError as e:
                print(f"🔥 Failed to download CIF for {pdb_id}: {e}")
                raise HTTPException(status_code=404, detail=f"mmCIF file for {pdb_id} not found.")

    def extract_chain_positions(self, structure, tubulin_chains: Dict[str, str]) -> Dict[str, np.ndarray]:
        """Extract the centroid of each tubulin chain."""
        chain_positions = {}
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                if chain_id in tubulin_chains:
                    atom_coords = [atom.get_coord() for atom in chain.get_atoms()]
                    if atom_coords:
                        chain_positions[chain_id] = np.mean(atom_coords, axis=0)
        print(f"📍 Extracted centroids for {len(chain_positions)} tubulin chains.")
        return chain_positions

    def visualize_point_cloud(self, positions: Dict[str, np.ndarray],
                              tubulin_chains: Dict[str, str], pdb_id: str, stage: str = "raw") -> str:
        """Create 3D scatter plot of chain positions colored by alpha/beta."""
        if not positions:
            print("⚠️ Cannot visualize point cloud: No positions provided.")
            return ""
            
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        colors = ['#1f77b4' if tubulin_chains.get(cid) == 'alpha' else '#ff7f0e' for cid in chain_ids]
        
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
        
        ax.set_xlabel('X (Å)', fontweight='bold')
        ax.set_ylabel('Y (Å)', fontweight='bold')
        ax.set_zlabel('Z (Å)', fontweight='bold')
        ax.set_title(f'3D Point Cloud - {pdb_id.upper()} ({stage.title()})\n{len(coords)} chains (Blue=α, Orange=β)',
                     fontsize=14, fontweight='bold')
        
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.grid(True, linestyle='--', alpha=0.6)
        
        alpha_patch = patches.Patch(color='#1f77b4', label='α-Tubulin')
        beta_patch = patches.Patch(color='#ff7f0e', label='β-Tubulin')
        ax.legend(handles=[alpha_patch, beta_patch])

        plt.tight_layout()
        viz_filename = f"pointcloud_{pdb_id}_{stage}.png"
        plt.savefig(viz_filename, dpi=150)
        plt.close()
        print(f"🎨 Saved 3D point cloud visualization: {viz_filename}")
        return viz_filename

    def detect_structure_type(self, coords: np.ndarray) -> str:
        """Classify structure type based on 3D geometry of chain centroids."""
        n_chains = len(coords)
        if n_chains <= 2:
            return "dimer"

        pca = PCA(n_components=3)
        pca.fit(coords)
        variances = pca.explained_variance_ratio_
        print(f"📊 PCA Explained Variance Ratios: {variances}")

        is_elongated = variances[0] > 0.7
        is_circular_cross_section = (variances[1] / variances[2] > 0.5) if variances[2] > 1e-6 else False

        if n_chains > 10 and is_elongated and is_circular_cross_section:
            return "cylindrical"
        elif n_chains > 4 and variances[2] < 0.05:
            return "sheet"
        elif variances[0] > 0.95:
            return "linear"
        else:
            return "sheet" if n_chains > 2 else "linear"

    def map_linear(self, positions: Dict[str, np.ndarray], tubulin_chains: Dict[str, str]) -> List[SubunitData]:
        """Map a linear arrangement of subunits."""
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        
        pca = PCA(n_components=1)
        projected = pca.fit_transform(coords).flatten()
        sorted_indices = np.argsort(projected)
        
        subunits = []
        for i, idx in enumerate(sorted_indices):
            chain_id = chain_ids[idx]
            subunits.append(SubunitData(
                id=f"pf0-{chain_id}",
                auth_asym_id=chain_id,
                protofilament=0,
                subunitIndex=i,
                monomerType=tubulin_chains[chain_id]
            ))
        return subunits

    def map_sheet(self, positions: Dict[str, np.ndarray], tubulin_chains: Dict[str, str]) -> List[SubunitData]:
        """Map a 2D sheet of subunits."""
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())

        pca = PCA(n_components=2)
        coords_2d = pca.fit_transform(coords)

        pf_clustering = DBSCAN(eps=25, min_samples=1).fit(coords_2d[:, 0].reshape(-1, 1))
        pf_labels = pf_clustering.labels_

        unique_labels = sorted(list(set(pf_labels)), key=lambda l: np.mean(coords_2d[pf_labels == l, 0]))
        pf_map = {label: i for i, label in enumerate(unique_labels)}
        
        subunits = []
        for pf_label in unique_labels:
            pf_idx = pf_map[pf_label]
            indices_in_pf = np.where(pf_labels == pf_label)[0]
            
            y_coords_in_pf = coords_2d[indices_in_pf, 1]
            sorted_pf_indices = indices_in_pf[np.argsort(y_coords_in_pf)]
            
            for i, chain_idx in enumerate(sorted_pf_indices):
                chain_id = chain_ids[chain_idx]
                subunits.append(SubunitData(
                    id=f"pf{pf_idx}-{chain_id}",
                    auth_asym_id=chain_id,
                    protofilament=pf_idx,
                    subunitIndex=i,
                    monomerType=tubulin_chains[chain_id]
                ))
        return subunits

    def map_cylindrical(self, positions: Dict[str, np.ndarray], tubulin_chains: Dict[str, str]) -> List[SubunitData]:
        """Unwrap a cylindrical structure into a 2D grid."""
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())

        pca = PCA(n_components=3)
        pca.fit(coords)
        axis_vector = pca.components_[0]

        coords_centered = coords - np.mean(coords, axis=0)
        axial_positions = np.dot(coords_centered, axis_vector)

        radial_vectors = coords_centered - np.outer(axial_positions, axis_vector)
        basis_v1 = pca.components_[1]
        basis_v2 = pca.components_[2]
        x_plane = np.dot(radial_vectors, basis_v1)
        y_plane = np.dot(radial_vectors, basis_v2)
        angles = np.arctan2(y_plane, x_plane)

        angle_coords = np.array([[math.cos(a), math.sin(a)] for a in angles])
        pf_clustering = DBSCAN(eps=0.3, min_samples=2).fit(angle_coords)
        pf_labels = pf_clustering.labels_

        unique_labels = sorted([l for l in set(pf_labels) if l != -1], 
                               key=lambda l: np.arctan2(
                                   np.mean(y_plane[pf_labels == l]),
                                   np.mean(x_plane[pf_labels == l])
                               ))
        pf_map = {label: i for i, label in enumerate(unique_labels)}
        
        subunits = []
        for pf_label in unique_labels:
            pf_idx = pf_map[pf_label]
            indices_in_pf = np.where(pf_labels == pf_label)[0]
            
            axial_pos_in_pf = axial_positions[indices_in_pf]
            sorted_pf_indices = indices_in_pf[np.argsort(axial_pos_in_pf)]

            for i, chain_idx in enumerate(sorted_pf_indices):
                chain_id = chain_ids[chain_idx]
                subunits.append(SubunitData(
                    id=f"pf{pf_idx}-{chain_id}",
                    auth_asym_id=chain_id,
                    protofilament=pf_idx,
                    subunitIndex=i,
                    monomerType=tubulin_chains[chain_id]
                ))
        return subunits

    def visualize_grid(self, grid_data: GridData, pdb_id: str) -> str:
        """Create a 2D visualization of the final grid layout."""
        subunits = grid_data.subunits
        if not subunits: return ""

        max_pf = max((s.protofilament for s in subunits), default=0)
        max_su = max((s.subunitIndex for s in subunits), default=0)
        
        fig_width = max(8, max_pf + 2)
        fig_height = max(6, (max_su + 2) * 0.8)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        for subunit in subunits:
            x, y = subunit.protofilament, max_su - subunit.subunitIndex
            color = 'lightblue' if subunit.monomerType == 'alpha' else '#FFC107' # Amber
            
            circle = patches.Circle((x, y), 0.45, facecolor=color, edgecolor='black', linewidth=1)
            ax.add_patch(circle)
            
            monomer_symbol = 'α' if subunit.monomerType == 'alpha' else 'β'
            ax.text(x, y, monomer_symbol, ha='center', va='center', fontsize=12, fontweight='bold', color='black')

        ax.set_xlim(-0.8, max_pf + 0.8)
        ax.set_ylim(-0.8, max_su + 0.8)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xticks(range(max_pf + 1))
        ax.set_yticks(range(max_su + 1))
        ax.set_yticklabels(reversed(range(max_su + 1)))
        ax.grid(True, linestyle='--', alpha=0.5)

        ax.set_xlabel('Protofilament Index', fontsize=12)
        ax.set_ylabel('Subunit Index', fontsize=12)
        ax.set_title(f'Grid Layout - {pdb_id.upper()} ({grid_data.structure_type})', fontsize=14, fontweight='bold')
        
        alpha_patch = patches.Patch(color='lightblue', label='α-Tubulin')
        beta_patch = patches.Patch(color='#FFC107', label='β-Tubulin')
        ax.legend(handles=[alpha_patch, beta_patch], loc='upper right', bbox_to_anchor=(1.1, 1.1))
        
        plt.tight_layout()
        viz_filename = f"grid_{pdb_id}_{grid_data.structure_type}.png"
        plt.savefig(viz_filename, dpi=150)
        plt.close()
        return viz_filename

    async def generate_grid(self, pdb_id: str) -> GridData:
        """The main pipeline to generate the 2D grid from a PDB ID."""
        profile = await self.get_profile(pdb_id)
        tubulin_chains = self.filter_tubulin_chains(profile)
        if not tubulin_chains:
            raise HTTPException(status_code=400, detail=f"No alpha or beta tubulin chains found in PDB profile for {pdb_id}.")

        mmcif_path = await self.download_mmcif(pdb_id)
        structure = self.parser.get_structure(pdb_id, mmcif_path)
        
        positions = self.extract_chain_positions(structure, tubulin_chains)
        if not positions:
            raise HTTPException(status_code=400, detail="Could not extract coordinates for tubulin chains.")
        
        self.visualize_point_cloud(positions, tubulin_chains, pdb_id, stage="raw")

        coords_array = np.array(list(positions.values()))
        structure_type = self.detect_structure_type(coords_array)
        print(f"🔬 Detected structure type: {structure_type.upper()}")

        mapping_functions = {
            "dimer": self.map_linear,
            "linear": self.map_linear,
            "sheet": self.map_sheet,
            "cylindrical": self.map_cylindrical,
        }
        
        if structure_type in mapping_functions:
            subunits = mapping_functions[structure_type](positions, tubulin_chains)
        else:
            raise HTTPException(status_code=500, detail=f"No mapping function for unknown structure type: {structure_type}")

        mmcif_path.unlink()
        
        grid_data = GridData(
            subunits=subunits,
            structure_type=structure_type,
            metadata={
                "pdb_id": pdb_id,
                "num_tubulin_chains": len(positions),
                "chain_ids": list(positions.keys())
            }
        )
        
        json_filename = f"grid_{pdb_id}.json"
        with open(json_filename, 'w') as f:
            f.write(grid_data.model_dump_json(indent=2))
        print(f"💾 Saved grid data: {json_filename}")

        viz_filename = self.visualize_grid(grid_data, pdb_id)
        print(f"🎨 Saved 2D grid visualization: {viz_filename}")

        return grid_data

# --- API Endpoint ---
grid_generator = SpatialGridGenerator()

@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """
    Generates a 2D spatial grid for tubulin subunits from a given PDB ID.
    - **Fetches** protein metadata to identify α and β tubulin chains.
    - **Downloads** the 3D structure from RCSB PDB.
    - **Plots** the 3D point cloud of subunits for debugging.
    - **Analyzes** the 3D coordinates to classify the structure.
    - **Maps** the arrangement to a 2D grid.
    - **Returns** the grid data as JSON and saves visualizations.
    """
    print(f"🚀 Received request for PDB ID: {pdb_id}")
    try:
        return await grid_generator.generate_grid(pdb_id)
    except HTTPException as e:
        raise e
    except Exception as e:
        print(f"🔥 An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"An internal error occurred: {str(e)}")


if __name__ == "__main__":
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)
