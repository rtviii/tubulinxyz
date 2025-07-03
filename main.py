from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Any, Tuple, Set
import numpy as np
import asyncio
import httpx
from pathlib import Path
import tempfile
import uvicorn
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Residue import Residue
from scipy.spatial._ckdtree import cKDTree
from scipy.spatial.transform import Rotation as R
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
import math

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
app = FastAPI(title="Tubulin Spatial Grid API", version="4.1.0",
              description="Generates an idealized, aligned 2D grid from 3D tubulin structures.")

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
        self.NEIGHBOR_CUTOFF_Å = 75.0
        self.TERMINUS_RESIDUE_COUNT = 20
        self.MIN_ANGLE_DEG = 150.0

    async def get_profile(self, pdb_id: str) -> Dict[str, Any]:
        profile_path = self.profile_base_path / f"{pdb_id.upper()}_profile.json"
        if profile_path.exists():
            with open(profile_path, 'r') as f: return json.load(f)
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        async with httpx.AsyncClient() as client:
            response = await client.get(url)
            response.raise_for_status()
            profile_data = response.json()
            with open(profile_path, 'w') as f: json.dump(profile_data, f)
            return profile_data

    def filter_tubulin_chains(self, profile: Dict[str, Any]) -> Dict[str, str]:
        chain_types: Dict[str, str] = {}
        entry_data = profile.get("entry", {})
        for entity in entry_data.get("polymer_entities", []):
            desc = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "").lower()
            if "tubulin" in desc and ("alpha" in desc or "beta" in desc):
                m_type = "alpha" if "alpha" in desc else "beta"
                ids = entity.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", [])
                for chain_id in ids: chain_types[chain_id] = m_type
        print(f"🧬 Found {len(chain_types)} tubulin chains.")
        return chain_types

    async def download_mmcif(self, pdb_id: str) -> Path:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.get(url)
            response.raise_for_status()
            temp_file = Path(tempfile.mktemp(suffix=".cif"))
            temp_file.write_bytes(response.content)
            return temp_file

    def _extract_chain_features(self, structure: Any, tubulin_chains: Dict[str, str]
                                ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        positions, orientations = {}, {}
        for model in structure:
            for chain in model:
                if chain.id not in tubulin_chains: continue
                residues: List[Residue.Residue] = [res for res in chain.get_residues() if "CA" in res]
                if len(residues) < self.TERMINUS_RESIDUE_COUNT * 2: continue
                all_atoms = [atom.get_coord() for res in residues for atom in res.get_atoms()]
                if not all_atoms: continue
                positions[chain.id] = np.mean(all_atoms, axis=0)
                residues.sort(key=lambda r: r.get_id()[1])
                n_term_atoms = [atom.get_coord() for res in residues[:self.TERMINUS_RESIDUE_COUNT] for atom in res.get_atoms()]
                c_term_atoms = [atom.get_coord() for res in residues[-self.TERMINUS_RESIDUE_COUNT:] for atom in res.get_atoms()]
                orientation_vec = np.mean(c_term_atoms, axis=0) - np.mean(n_term_atoms, axis=0)
                norm = np.linalg.norm(orientation_vec)
                if norm > 1e-6: orientations[chain.id] = orientation_vec / norm
        print(f"📍 Extracted features for {len(positions)} chains.")
        return positions, orientations

    def _build_neighborhood_graph(self, positions: Dict[str, np.ndarray]) -> Dict[str, List[str]]:
        chain_ids = list(positions.keys())
        coords = np.array([positions[cid] for cid in chain_ids])
        kdtree = cKDTree(coords)
        pairs = kdtree.query_pairs(r=self.NEIGHBOR_CUTOFF_Å)
        graph = {cid: [] for cid in chain_ids}
        for i, j in pairs:
            id1, id2 = chain_ids[i], chain_ids[j]
            graph[id1].append(id2)
            graph[id2].append(id1)
        print(f"📊 Built neighborhood graph with {len(pairs)} connections.")
        return graph

    def _find_dimer_pairs_by_orientation(self, graph: Dict[str, List[str]], positions: Dict[str, np.ndarray],
                                         orientations: Dict[str, np.ndarray], tubulin_chains: Dict[str, str]) -> Dict[str, str]:
        dimer_pairs = {}
        alpha_chains = {cid for cid, mtype in tubulin_chains.items() if mtype == 'alpha' and cid in positions}
        for alpha_id in alpha_chains:
            beta_neighbors = [nid for nid in graph.get(alpha_id, []) if tubulin_chains.get(nid) == 'beta' and nid in orientations]
            if not beta_neighbors: continue
            min_score, best_partner = float('inf'), None
            v_alpha, pos_alpha = orientations[alpha_id], positions[alpha_id]
            for beta_id in beta_neighbors:
                d_vec = positions[beta_id] - pos_alpha
                d_vec /= np.linalg.norm(d_vec)
                score = abs(np.dot(v_alpha, d_vec)) + abs(np.dot(orientations[beta_id], d_vec))
                if score < min_score:
                    min_score, best_partner = score, beta_id
            if best_partner:
                dimer_pairs[alpha_id] = best_partner
                dimer_pairs[best_partner] = alpha_id
        print(f"🤝 Identified {len(dimer_pairs)//2} α-β dimer pairs.")
        return dimer_pairs

    def _trace_protofilaments(self, graph: Dict[str, List[str]], dimer_pairs: Dict[str, str],
                              tubulin_chains: Dict[str, str], positions: Dict[str, np.ndarray]) -> List[List[str]]:
        dimer_graph: Dict[str, str] = {}
        alpha_chains = {cid for cid, mtype in tubulin_chains.items() if mtype == 'alpha' and cid in dimer_pairs}
        for alpha_id in alpha_chains:
            beta_partner = dimer_pairs[alpha_id]
            alpha_neighbors = [nid for nid in graph.get(beta_partner, []) 
                               if tubulin_chains.get(nid) == 'alpha' and nid != alpha_id]
            if not alpha_neighbors: continue
            best_candidate, max_angle = None, -1.0
            for candidate_id in alpha_neighbors:
                vec1 = positions[alpha_id] - positions[beta_partner]
                vec2 = positions[candidate_id] - positions[beta_partner]
                cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
                if angle > max_angle:
                    max_angle, best_candidate = angle, candidate_id
            if best_candidate and max_angle > np.deg2rad(self.MIN_ANGLE_DEG):
                dimer_graph[alpha_id] = best_candidate
        
        all_successors = set(dimer_graph.values())
        starters = [alpha_id for alpha_id in dimer_graph.keys() if alpha_id not in all_successors]
        
        protofilaments = []
        visited: Set[str] = set()
        for starter_alpha in starters:
            if starter_alpha in visited: continue
            current_pf, current_alpha = [], starter_alpha
            while current_alpha and current_alpha not in visited:
                visited.add(current_alpha)
                beta_partner = dimer_pairs.get(current_alpha)
                if beta_partner:
                    current_pf.extend([current_alpha, beta_partner])
                    visited.add(beta_partner)
                current_alpha = dimer_graph.get(current_alpha)
            if current_pf: protofilaments.append(current_pf)
        print(f"🧵 Traced {len(protofilaments)} protofilaments.")
        return protofilaments
    
    def _finalize_grid_idealized(self, protofilaments: List[List[str]], positions: Dict[str, np.ndarray],
                                 tubulin_chains: Dict[str, str]) -> Tuple[List[SubunitData], str]:
        if not protofilaments: return [], "unknown"

        all_chains_in_pfs = [cid for pf in protofilaments for cid in pf]
        all_pf_coords = np.array([positions[cid] for cid in all_chains_in_pfs])
        pca_global = PCA(n_components=3).fit(all_pf_coords)
        
        # --- BUGFIX IS HERE ---
        # Use the robust R.align_vectors method to find the rotation
        target_z_axis = np.array([0, 0, 1])
        source_axis = pca_global.components_[0]
        rot, _ = R.align_vectors([target_z_axis], [source_axis])
        rotated_positions = {cid: rot.apply(pos) for cid, pos in positions.items()}

        pf_centroids = [np.mean([rotated_positions[cid] for cid in pf], axis=0) for pf in protofilaments]
        angles = np.arctan2([c[1] for c in pf_centroids], [c[0] for c in pf_centroids])
        ordered_pfs = [pf for _, pf in sorted(zip(angles, protofilaments))]
        structure_type = "cylindrical" if len(ordered_pfs) > 4 else "sheet"

        z_coords_map = {cid: rotated_positions[cid][2] for cid in all_chains_in_pfs}
        chain_order = list(z_coords_map.keys())
        z_coords_array = np.array([z_coords_map[cid] for cid in chain_order]).reshape(-1, 1)
        z_clusters = DBSCAN(eps=20, min_samples=1).fit(z_coords_array)
        
        chain_to_z_label = {chain: label for chain, label in zip(chain_order, z_clusters.labels_)}
        unique_z_labels = sorted(list(set(z_clusters.labels_)), 
                                 key=lambda l: np.mean(z_coords_array[z_clusters.labels_ == l]))
        z_level_map = {label: i for i, label in enumerate(unique_z_labels)}
        
        all_subunits = []
        for pf_idx, pf_chains in enumerate(ordered_pfs):
            for chain_id in pf_chains:
                if chain_id not in chain_to_z_label: continue
                cluster_label = chain_to_z_label[chain_id]
                subunit_idx = z_level_map.get(cluster_label, -1)
                all_subunits.append(SubunitData(id=f"pf{pf_idx}-{chain_id}", auth_asym_id=chain_id,
                                                protofilament=pf_idx, subunitIndex=subunit_idx,
                                                monomerType=tubulin_chains[chain_id]))

        print(f"✅ Finalized idealized grid. Structure type determined as: {structure_type.upper()}")
        return all_subunits, structure_type

    async def generate_grid(self, pdb_id: str) -> GridData:
        profile = await self.get_profile(pdb_id)
        tubulin_chains = self.filter_tubulin_chains(profile)
        if not tubulin_chains: raise HTTPException(400, f"No tubulin chains for {pdb_id}.")
        mmcif_path = await self.download_mmcif(pdb_id)
        structure = self.parser.get_structure(pdb_id, mmcif_path)
        positions, orientations = self._extract_chain_features(structure, tubulin_chains)
        if not positions: raise HTTPException(400, "Could not extract features.")
        graph = self._build_neighborhood_graph(positions)
        dimer_pairs = self._find_dimer_pairs_by_orientation(graph, positions, orientations, tubulin_chains)
        protofilaments = self._trace_protofilaments(graph, dimer_pairs, tubulin_chains, positions)
        subunits, structure_type = self._finalize_grid_idealized(protofilaments, positions, tubulin_chains)
        mmcif_path.unlink()
        grid_data = GridData(subunits=subunits, structure_type=structure_type, 
                             metadata={"pdb_id": pdb_id, "num_tubulin_chains": len(positions)})
        self.visualize_grid(grid_data, pdb_id)
        return grid_data

    def visualize_grid(self, grid_data: GridData, pdb_id: str):
        if not grid_data.subunits: return
        max_pf = max((s.protofilament for s in grid_data.subunits), default=0)
        max_su = max((s.subunitIndex for s in grid_data.subunits), default=0)
        fig_width = max(8, (max_pf + 1) * 1.0)
        fig_height = max(6, (max_su + 1) * 1.0)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        for s in grid_data.subunits:
            x, y = s.protofilament, s.subunitIndex
            color = '#1f77b4' if s.monomerType == 'alpha' else '#ff7f0e'
            ax.add_patch(patches.Circle((x, y), 0.45, fc=color, ec='black', lw=0.5))
            ax.text(x, y, s.auth_asym_id, ha='center', va='center',
                    fontsize=6, weight='normal', color='white')
        ax.set_aspect('equal')
        ax.set_xlim(-0.7, max_pf + 0.7)
        ax.set_ylim(-0.7, max_su + 0.7)
        ax.set_xticks(range(max_pf + 1))
        ax.set_xticklabels([f"PF{i}" for i in range(max_pf + 1)])
        ax.set_yticks(range(max_su + 1))
        ax.set_xlabel("Protofilament Index", fontsize=12)
        ax.set_ylabel("Subunit Index", fontsize=12)
        ax.invert_yaxis()
        ax.set_title(f'Grid Layout - {pdb_id.upper()} ({grid_data.structure_type})', fontsize=14, weight='bold')
        alpha_patch = patches.Patch(color='#1f77b4', label='α-Tubulin')
        beta_patch = patches.Patch(color='#ff7f0e', label='β-Tubulin')
        ax.legend(handles=[alpha_patch, beta_patch], bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout(rect=[0, 0, 0.9, 1])
        viz_filename = f"grid_{pdb_id}_{grid_data.structure_type}.png"
        plt.savefig(viz_filename, dpi=150)
        plt.close()
        print(f"🎨 Saved final 2D grid visualization: {viz_filename}")

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
