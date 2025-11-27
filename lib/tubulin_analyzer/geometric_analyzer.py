import numpy as np
from typing import List, Dict, Tuple
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from scipy.spatial.transform import Rotation as R
import json
from pathlib import Path


class GeometricAnalyzer:
    """Handles geometric operations and transformations"""

    def __init__(self):
        self.Z_CLUSTERING_EPS = 20.0
        self.MIN_SAMPLES = 1

    def align_structure_to_z_axis(
        self, positions: Dict[str, np.ndarray], protofilaments: List[List[str]]
    ) -> Dict[str, np.ndarray]:
        """Align the entire structure so main axis is along Z"""
        all_chains_in_pfs = [cid for pf in protofilaments for cid in pf]
        if not all_chains_in_pfs:
            return positions

        all_pf_coords = np.array(
            [positions[cid] for cid in all_chains_in_pfs if cid in positions]
        )
        if len(all_pf_coords) < 3:
            return positions

        # Find main axis using PCA
        pca_global = PCA(n_components=3).fit(all_pf_coords)

        # Align main axis to Z-axis
        target_z_axis = np.array([0, 0, 1])
        source_axis = pca_global.components_[0]

        try:
            rot, _ = R.align_vectors([target_z_axis], [source_axis])
            rotated_positions = {cid: rot.apply(pos) for cid, pos in positions.items()}
            return rotated_positions
        except Exception as e:
            print(f"Warning: Failed to align structure: {e}")
            return positions

    def order_protofilaments_by_angle(
        self, protofilaments: List[List[str]], positions: Dict[str, np.ndarray]
    ) -> List[List[str]]:
        """Order protofilaments by their angular position around Z-axis"""
        if not protofilaments:
            return []

        # Calculate centroid of each protofilament
        pf_centroids = []
        valid_pfs = []

        for pf in protofilaments:
            pf_coords = [positions[cid] for cid in pf if cid in positions]
            if pf_coords:
                centroid = np.mean(pf_coords, axis=0)
                pf_centroids.append(centroid)
                valid_pfs.append(pf)

        if not pf_centroids:
            return []

        # Calculate angles around Z-axis
        angles = []
        for centroid in pf_centroids:
            angle = np.arctan2(centroid[1], centroid[0])
            angles.append(angle)

        # Sort protofilaments by angle
        ordered_pfs = [pf for _, pf in sorted(zip(angles, valid_pfs))]
        return ordered_pfs

    def cluster_z_coordinates_per_protofilament(
        self, positions: Dict[str, np.ndarray], protofilaments: List[List[str]]
    ) -> Dict[str, int]:
        """Cluster chains by Z-coordinate within each protofilament to create subunit indices"""
        chain_to_subunit_idx = {}

        for pf_idx, pf_chains in enumerate(protofilaments):
            if not pf_chains:
                continue

            # Get Z-coordinates for this protofilament
            chain_z_data = []
            for chain_id in pf_chains:
                if chain_id in positions:
                    z = positions[chain_id][2]
                    chain_z_data.append((chain_id, z))

            if len(chain_z_data) < 2:
                # Single chain or no chains
                for chain_id, _ in chain_z_data:
                    chain_to_subunit_idx[chain_id] = 0
                continue

            # Sort chains by Z-coordinate
            chain_z_data.sort(key=lambda x: x[1])

            # Group chains by similar Z-coordinates (same layer)
            Z_LAYER_TOLERANCE = 5.0
            layers = []  # List of [chain_ids] at each layer
            current_layer = [chain_z_data[0]]

            for i in range(1, len(chain_z_data)):
                chain_id, z = chain_z_data[i]
                prev_z = current_layer[-1][1]

                if abs(z - prev_z) <= Z_LAYER_TOLERANCE:
                    current_layer.append((chain_id, z))
                else:
                    layers.append(current_layer)
                    current_layer = [(chain_id, z)]

            layers.append(current_layer)

            # Assign subunit indices based on layer
            for layer_idx, layer_chains in enumerate(layers):
                for chain_id, z in layer_chains:
                    chain_to_subunit_idx[chain_id] = layer_idx

        return chain_to_subunit_idx

    def determine_structure_type(self, num_protofilaments: int) -> str:
        """Determine if structure is cylindrical or sheet-like"""
        if num_protofilaments > 6:
            return "cylindrical"
        else:
            return "sheet"

    def save_geometric_debug(
        self,
        pdb_id: str,
        protofilaments: List[List[str]],
        positions: Dict[str, np.ndarray],
        tubulin_chains: Dict[str, str],
    ):
        """Save geometric analysis debug data"""
        debug_dir = Path("debug_output")
        debug_dir.mkdir(exist_ok=True)

        geometric_data = {
            "pdb_id": pdb_id,
            "num_protofilaments": len(protofilaments),
            "structure_type": self.determine_structure_type(len(protofilaments)),
            "protofilament_details": [],
        }

        # Align structure for analysis
        aligned_positions = self.align_structure_to_z_axis(positions, protofilaments)
        ordered_pfs = self.order_protofilaments_by_angle(
            protofilaments, aligned_positions
        )

        for i, pf_chains in enumerate(ordered_pfs):
            if not pf_chains:
                continue

            # Calculate protofilament centroid
            pf_coords = [
                aligned_positions[cid] for cid in pf_chains if cid in aligned_positions
            ]
            if pf_coords:
                centroid = np.mean(pf_coords, axis=0)
                angle = np.arctan2(centroid[1], centroid[0])

                pf_info = {
                    "index": i,
                    "chain_count": len(pf_chains),
                    "chains": pf_chains,
                    "centroid": [
                        float(centroid[0]),
                        float(centroid[1]),
                        float(centroid[2]),
                    ],
                    "angle_rad": float(angle),
                    "angle_deg": float(np.degrees(angle)),
                    "chain_types": [
                        tubulin_chains.get(c, "unknown") for c in pf_chains
                    ],
                }

                geometric_data["protofilament_details"].append(pf_info)

        with open(debug_dir / f"{pdb_id}_geometric_debug.json", "w") as f:
            json.dump(geometric_data, f, indent=2)

    def process_protofilaments(
        self,
        protofilaments: List[List[str]],
        positions: Dict[str, np.ndarray],
        tubulin_chains: Dict[str, str],
    ) -> Tuple[List, str]:
        """Process protofilaments into final grid by unrolling cylinder"""
        if not protofilaments:
            return [], "unknown"

        # Import here to avoid circular imports
        from lib.tubulin_analyzer.models import SubunitData

        # Align structure to Z-axis for consistent orientation
        aligned_positions = self.align_structure_to_z_axis(positions, protofilaments)

        # Order protofilaments by angular position around Z-axis
        ordered_pfs = self.order_protofilaments_by_angle(
            protofilaments, aligned_positions
        )

        # Determine structure type
        structure_type = self.determine_structure_type(len(ordered_pfs))

        # Create subunit indices by clustering Z-coordinates within each protofilament
        chain_to_subunit_idx = self.cluster_z_coordinates_per_protofilament(
            aligned_positions, ordered_pfs
        )

        # Create subunit data for the unrolled grid
        subunits = []
        seen_chains = set()

        for pf_idx, pf_chains in enumerate(ordered_pfs):
            for chain_id in pf_chains:
                if (
                    chain_id in chain_to_subunit_idx
                    and chain_id in tubulin_chains
                    and chain_id not in seen_chains
                ):
                    subunit_idx = chain_to_subunit_idx[chain_id]
                    subunit = SubunitData(
                        id=f"pf{pf_idx}-su{subunit_idx}-{chain_id}",
                        auth_asym_id=chain_id,
                        protofilament=pf_idx,
                        subunitIndex=subunit_idx,
                        monomerType=tubulin_chains[chain_id],
                    )
                    subunits.append(subunit)
                    seen_chains.add(chain_id)

        return subunits, structure_type

    def calculate_interface_center(self, residues: List) -> np.ndarray:
        """Calculate the center of mass of interface residues"""
        if not residues:
            return np.array([0, 0, 0])

        all_coords = []
        for residue in residues:
            for atom in residue:
                all_coords.append(atom.get_coord())

        if all_coords:
            return np.mean(all_coords, axis=0)
        else:
            return np.array([0, 0, 0])

    def calculate_vector_between_points(
        self, point1: np.ndarray, point2: np.ndarray
    ) -> np.ndarray:
        """Calculate normalized vector between two points"""
        vector = point2 - point1
        norm = np.linalg.norm(vector)
        if norm > 1e-6:
            return vector / norm
        else:
            return np.array([0, 0, 0])

    def calculate_angle_between_vectors(
        self, vec1: np.ndarray, vec2: np.ndarray
    ) -> float:
        """Calculate angle between two vectors in degrees"""
        dot_product = np.dot(vec1, vec2)
        dot_product = np.clip(dot_product, -1.0, 1.0)
        angle_rad = np.arccos(dot_product)
        return np.degrees(angle_rad)