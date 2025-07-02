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
from scipy.spatial.distance import pdist
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA


# Pydantic models
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


# FastAPI app
app = FastAPI(title="Spatial Grid API", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


class SpatialGridGenerator:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
        self.profile_base_path = (
            "/Users/rtviii/dev/tubulinxyz"  # Hardcoded base path for profiles
        )

    def load_profile(self, pdb_id: str) -> Dict[str, Any]:
        """Load PDB profile JSON to get entity information"""
        profile_path = os.path.join(
            self.profile_base_path, f"{pdb_id.upper()}_profile.json"
        )
        try:
            print(f"🧬 Loading profile: {profile_path}")
            with open(profile_path, "r") as f:
                jj= json.load(f)
                print(f"🧬 Loaded profile: {jj}")
                return jj

        except FileNotFoundError:
            print(f"⚠️  Profile not found: {profile_path}")
            return {}

    def filter_tubulin_chains(self, profile: Dict[str, Any]) -> Dict[str, str]:
        """Extract alpha/beta tubulin chain mappings from profile
        Returns: {chain_id: 'alpha'/'beta'}
        """
        chain_types = {}

        if not profile or "entry" not in profile:
            return chain_types

        polymer_entities = profile["entry"].get("polymer_entities", [])

        for entity in polymer_entities:
            entity_id = entity["rcsb_polymer_entity_container_identifiers"]["entity_id"]
            description = entity["rcsb_polymer_entity"]["pdbx_description"].lower()
            auth_asym_ids = entity["rcsb_polymer_entity_container_identifiers"][
                "auth_asym_ids"
            ]

            # Determine if this entity is alpha or beta tubulin
            if "tubulin alpha" in description or "alpha tubulin" in description:
                for chain_id in auth_asym_ids:
                    chain_types[chain_id] = "alpha"
            elif "tubulin beta" in description or "beta tubulin" in description:
                for chain_id in auth_asym_ids:
                    chain_types[chain_id] = "beta"

        print(f"🧬 Found tubulin chains: {len(chain_types)} total")
        alpha_count = sum(1 for t in chain_types.values() if t == "alpha")
        beta_count = sum(1 for t in chain_types.values() if t == "beta")
        print(f"   α-tubulin: {alpha_count}, β-tubulin: {beta_count}")

        return chain_types

    async def download_mmcif(self, pdb_id: str) -> Path:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        async with httpx.AsyncClient() as client:
            response = await client.get(url)
            response.raise_for_status()
            temp_file = Path(tempfile.mktemp(suffix=".cif"))
            temp_file.write_bytes(response.content)
            return temp_file

    def extract_chain_positions(self, structure) -> Dict[str, tuple]:
        chain_positions = {}
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                residues = list(chain.get_residues())
                if not residues:
                    continue
                first_residue = residues[0]
                atoms = list(first_residue.get_atoms())
                if not atoms:
                    continue
                coord = atoms[0].get_coord()
                chain_positions[chain_id] = tuple(coord)
        return chain_positions

    def detect_structure_type(self, positions: Dict[str, tuple]) -> str:
        coords = np.array(list(positions.values()))
        n_chains = len(coords)

        if n_chains <= 2:
            return "dimer"
        elif n_chains <= 4:
            return "linear"

        if self._is_cylindrical(coords):
            return "cylindrical"
        elif self._is_sheet_like(coords):
            return "sheet"
        else:
            return "linear"

    def _is_cylindrical(self, coords: np.ndarray) -> bool:
        if len(coords) < 8:
            return False
        pca = PCA(n_components=3)
        pca.fit(coords)
        coords_2d = pca.transform(coords)[:, :2]
        radial_distances = np.linalg.norm(coords_2d, axis=1)
        radial_std = np.std(radial_distances)
        radial_mean = np.mean(radial_distances)
        cv = radial_std / radial_mean if radial_mean > 0 else float("inf")
        return cv < 0.3 and len(coords) > 12

    def _is_sheet_like(self, coords: np.ndarray) -> bool:
        if len(coords) < 4:
            return False
        pca = PCA(n_components=3)
        pca.fit(coords)
        explained_variance = pca.explained_variance_ratio_
        return explained_variance[2] < 0.2

    def map_dimer(self, positions: Dict[str, tuple]) -> List[SubunitData]:
        chain_ids = list(positions.keys())
        subunits = []
        for i, chain_id in enumerate(chain_ids):
            subunits.append(
                SubunitData(
                    id=f"pf{i}-{chain_id}",
                    auth_asym_id=chain_id,
                    protofilament=i,
                    subunitIndex=0,
                    monomerType="α" if i % 2 == 0 else "β",
                )
            )
        return subunits

    def map_linear(self, positions: Dict[str, tuple]) -> List[SubunitData]:
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        pca = PCA(n_components=1)
        pca.fit(coords)
        projected = pca.transform(coords).flatten()
        sorted_indices = np.argsort(projected)

        subunits = []
        for i, idx in enumerate(sorted_indices):
            chain_id = chain_ids[idx]
            subunits.append(
                SubunitData(
                    id=f"pf0-{chain_id}",
                    auth_asym_id=chain_id,
                    protofilament=0,
                    subunitIndex=i,
                    monomerType="α" if i % 2 == 0 else "β",
                )
            )
        return subunits

    def map_sheet(self, positions: Dict[str, tuple]) -> List[SubunitData]:
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())
        pca = PCA(n_components=2)
        pca.fit(coords)
        coords_2d = pca.transform(coords)

        x_coords = coords_2d[:, 0]
        y_coords = coords_2d[:, 1]
        x_clusters = self._cluster_coordinates(x_coords)

        subunits = []
        for chain_idx, chain_id in enumerate(chain_ids):
            pf = x_clusters[chain_idx]
            pf_indices = [i for i, cluster in enumerate(x_clusters) if cluster == pf]
            pf_y_coords = [y_coords[i] for i in pf_indices]
            pf_chain_ids = [chain_ids[i] for i in pf_indices]

            y_sorted = sorted(enumerate(pf_y_coords), key=lambda x: x[1])
            chain_pos_in_pf = next(
                i
                for i, (orig_i, _) in enumerate(y_sorted)
                if pf_chain_ids[orig_i] == chain_id
            )

            subunits.append(
                SubunitData(
                    id=f"pf{pf}-{chain_id}",
                    auth_asym_id=chain_id,
                    protofilament=pf,
                    subunitIndex=chain_pos_in_pf,
                    monomerType="α" if (pf + chain_pos_in_pf) % 2 == 0 else "β",
                )
            )
        return subunits

    def map_cylindrical(self, positions: Dict[str, tuple]) -> List[SubunitData]:
        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())

        pca = PCA(n_components=3)
        pca.fit(coords)
        axis = pca.components_[0]
        coords_centered = coords - np.mean(coords, axis=0)
        projected = coords_centered - np.outer(np.dot(coords_centered, axis), axis)

        x, y = projected[:, 0], projected[:, 1]
        angles = np.arctan2(y, x)
        angles = (angles + 2 * np.pi) % (2 * np.pi)
        z_coords = np.dot(coords_centered, axis)

        angle_clusters = self._cluster_angles(angles)

        subunits = []
        for chain_idx, chain_id in enumerate(chain_ids):
            pf = angle_clusters[chain_idx]
            pf_indices = [
                i for i, cluster in enumerate(angle_clusters) if cluster == pf
            ]
            pf_z_coords = [z_coords[i] for i in pf_indices]
            pf_chain_ids = [chain_ids[i] for i in pf_indices]

            z_sorted = sorted(enumerate(pf_z_coords), key=lambda x: x[1])
            chain_pos_in_pf = next(
                i
                for i, (orig_i, _) in enumerate(z_sorted)
                if pf_chain_ids[orig_i] == chain_id
            )

            subunits.append(
                SubunitData(
                    id=f"pf{pf}-{chain_id}",
                    auth_asym_id=chain_id,
                    protofilament=pf,
                    subunitIndex=chain_pos_in_pf,
                    monomerType="α" if (pf + chain_pos_in_pf) % 2 == 0 else "β",
                )
            )
        return subunits

    def _cluster_coordinates(
        self, coords: np.ndarray, threshold: float = None
    ) -> List[int]:
        if threshold is None:
            threshold = np.std(coords) / 2
        coords_reshaped = coords.reshape(-1, 1)
        clustering = DBSCAN(eps=threshold, min_samples=1).fit(coords_reshaped)
        return clustering.labels_.tolist()

    def visualize_point_cloud(
        self,
        positions: Dict[str, Tuple[float, float, float]],
        tubulin_chains: Dict[str, str],
        pdb_id: str,
        stage: str = "raw",
    ) -> str:
        """Create 3D scatter plot of chain positions colored by alpha/beta"""
        from mpl_toolkits.mplot3d import Axes3D

        coords = np.array(list(positions.values()))
        chain_ids = list(positions.keys())

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="3d")

        # Color by alpha/beta type
        colors = []
        for chain_id in chain_ids:
            chain_type = tubulin_chains.get(chain_id, "unknown")
            if chain_type == "alpha":
                colors.append("blue")
            elif chain_type == "beta":
                colors.append("orange")
            else:
                colors.append("gray")

        # Create scatter plot
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=50, alpha=0.7
        )

        # Add labels for some points
        for i, (chain_id, coord) in enumerate(zip(chain_ids, coords)):
            if i % 10 == 0:  # Label every 10th point to avoid clutter
                ax.text(coord[0], coord[1], coord[2], chain_id, fontsize=8)

        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
        ax.set_title(
            f"{pdb_id.upper()} - {stage.title()} Point Cloud\n"
            f"{len(coords)} chains (Blue=α, Orange=β)"
        )

        # Add legend
        alpha_proxy = plt.Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="blue", markersize=8
        )
        beta_proxy = plt.Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="orange", markersize=8
        )
        ax.legend([alpha_proxy, beta_proxy], ["α-tubulin", "β-tubulin"])

        plt.tight_layout()

        # Save visualization
        viz_filename = f"pointcloud_{pdb_id}_{stage}.png"
        plt.savefig(viz_filename, dpi=150, bbox_inches="tight")
        plt.close()

        return viz_filename

    def visualize_grid(self, grid_data: GridData, pdb_id: str) -> str:
        """Create visualization of the grid layout"""
        subunits = grid_data.subunits

        if not subunits:
            return ""

        # Get grid dimensions
        max_pf = max(s.protofilament for s in subunits)
        max_su = max(s.subunitIndex for s in subunits)

        # Create figure
        fig, ax = plt.subplots(figsize=(max(8, max_pf + 2), max(6, max_su + 2)))

        # Plot each subunit
        for subunit in subunits:
            x = subunit.protofilament
            y = max_su - subunit.subunitIndex  # Flip Y so 0 is at top

            # Color based on monomer type
            color = "lightblue" if subunit.monomerType == "alpha" else "orange"

            # Draw circle
            circle = patches.Circle(
                (x, y), 0.4, facecolor=color, edgecolor="black", linewidth=1.5
            )
            ax.add_patch(circle)

            # Add chain ID label
            ax.text(
                x,
                y,
                subunit.auth_asym_id,
                ha="center",
                va="center",
                fontsize=12,
                fontweight="bold",
            )

        # Customize plot
        ax.set_xlim(-0.8, max_pf + 0.8)
        ax.set_ylim(-0.8, max_su + 0.8)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Protofilament", fontsize=12)
        ax.set_ylabel("Subunit Index", fontsize=12)
        ax.set_title(
            f"Grid Layout - {pdb_id.upper()} ({grid_data.structure_type})",
            fontsize=14,
            fontweight="bold",
        )

        # Add protofilament labels
        for pf in range(max_pf + 1):
            ax.text(
                pf,
                -0.6,
                f"PF{pf}",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
            )

        # Add legend
        alpha_patch = patches.Patch(color="lightblue", label="α monomer")
        beta_patch = patches.Patch(color="orange", label="β monomer")
        ax.legend(handles=[alpha_patch, beta_patch], loc="upper right")

        plt.tight_layout()

        # Save visualization
        viz_filename = f"grid_{pdb_id}_{grid_data.structure_type}.png"
        plt.savefig(viz_filename, dpi=150, bbox_inches="tight")
        plt.close()

        return viz_filename
        """Create visualization of the grid layout"""
        subunits = grid_data.subunits

        if not subunits:
            return ""

        # Get grid dimensions
        max_pf = max(s.protofilament for s in subunits)
        max_su = max(s.subunitIndex for s in subunits)

        # Create figure
        fig, ax = plt.subplots(figsize=(max(8, max_pf + 2), max(6, max_su + 2)))

        # Plot each subunit
        for subunit in subunits:
            x = subunit.protofilament
            y = max_su - subunit.subunitIndex  # Flip Y so 0 is at top

            # Color based on monomer type
            color = "lightblue" if subunit.monomerType == "α" else "orange"

            # Draw circle
            circle = patches.Circle(
                (x, y), 0.4, facecolor=color, edgecolor="black", linewidth=1.5
            )
            ax.add_patch(circle)

            # Add chain ID label
            ax.text(
                x,
                y,
                subunit.auth_asym_id,
                ha="center",
                va="center",
                fontsize=12,
                fontweight="bold",
            )

        # Customize plot
        ax.set_xlim(-0.8, max_pf + 0.8)
        ax.set_ylim(-0.8, max_su + 0.8)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Protofilament", fontsize=12)
        ax.set_ylabel("Subunit Index", fontsize=12)
        ax.set_title(
            f"Grid Layout - {pdb_id.upper()} ({grid_data.structure_type})",
            fontsize=14,
            fontweight="bold",
        )

        # Add protofilament labels
        for pf in range(max_pf + 1):
            ax.text(
                pf,
                -0.6,
                f"PF{pf}",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
            )

        # Add legend
        alpha_patch = patches.Patch(color="lightblue", label="α monomer")
        beta_patch = patches.Patch(color="orange", label="β monomer")
        ax.legend(handles=[alpha_patch, beta_patch], loc="upper right")

        plt.tight_layout()

        # Save visualization
        viz_filename = f"grid_{pdb_id}_{grid_data.structure_type}.png"
        plt.savefig(viz_filename, dpi=150, bbox_inches="tight")
        plt.close()

        return viz_filename

    def _cluster_angles(self, angles: np.ndarray) -> List[int]:
        sorted_indices = np.argsort(angles)
        sorted_angles = angles[sorted_indices]
        angle_diffs = np.diff(sorted_angles)
        wraparound_diff = (sorted_angles[0] + 2 * np.pi) - sorted_angles[-1]
        angle_diffs = np.append(angle_diffs, wraparound_diff)
        threshold = np.mean(angle_diffs) + np.std(angle_diffs)

        clusters = np.zeros(len(angles), dtype=int)
        current_cluster = 0

        for i in range(len(sorted_indices)):
            clusters[sorted_indices[i]] = current_cluster
            if i < len(angle_diffs) and angle_diffs[i] > threshold:
                current_cluster += 1
        return clusters.tolist()

    async def generate_grid(self, pdb_id: str) -> GridData:
        mmcif_path = await self.download_mmcif(pdb_id)
        structure = self.parser.get_structure(pdb_id, mmcif_path)
        positions = self.extract_chain_positions(structure)

        if not positions:
            raise HTTPException(status_code=400, detail="No chains found")

        structure_type = self.detect_structure_type(positions)

        if structure_type == "dimer":
            subunits = self.map_dimer(positions)
        elif structure_type == "linear":
            subunits = self.map_linear(positions)
        elif structure_type == "sheet":
            subunits = self.map_sheet(positions)
        elif structure_type == "cylindrical":
            subunits = self.map_cylindrical(positions)
        else:
            raise HTTPException(
                status_code=400, detail=f"Unknown structure type: {structure_type}"
            )

        mmcif_path.unlink()

        grid_data = GridData(
            subunits=subunits,
            structure_type=structure_type,
            metadata={
                "pdb_id": pdb_id,
                "num_chains": len(positions),
                "chain_ids": list(positions.keys()),
            },
        )

        # Save as JSON
        json_filename = f"grid_{pdb_id}_{structure_type}.json"

        print("dumping grid data")
        with open(json_filename, "w") as f:
            json.dump(grid_data.model_dump(), f, indent=2)
        print("dumped grid data")
        print(f"💾 Saved grid data: {json_filename}")

        # Create visualization
        viz_filename = self.visualize_grid(grid_data, pdb_id)
        print(f"🎨 Saved visualization: {viz_filename}")

        return grid_data


# Single endpoint
grid_generator = SpatialGridGenerator()


@app.get("/grid/{pdb_id}", response_model=GridData)
async def get_grid(pdb_id: str):
    """Generate 2D grid from PDB structure"""
    print(f"🚀 NEW VERSION RUNNING - Processing {pdb_id}")  # Debug line
    return await grid_generator.generate_grid(pdb_id)


if __name__ == "__main__":
    uvicorn.run(
        "dev_server:app", host="127.0.0.1", port=8000, reload=True, log_level="info"
    )
