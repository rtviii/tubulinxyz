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
        # Get all positions in protofilaments
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

        # Use robust alignment
        try:
            rot, _ = R.align_vectors([target_z_axis], [source_axis])
            rotated_positions = {cid: rot.apply(pos) for cid, pos in positions.items()}
            print(f"🔄 Aligned structure to Z-axis")
            return rotated_positions
        except Exception as e:
            print(f"⚠️ Failed to align structure: {e}")
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
        print(f"📐 Ordered {len(ordered_pfs)} protofilaments by angular position")
        return ordered_pfs

    def cluster_z_coordinates(
        self, positions: Dict[str, np.ndarray], chain_list: List[str]
    ) -> Dict[str, int]:
        """Cluster chains by their Z-coordinate to create subunit indices"""
        if not chain_list:
            return {}

        # Get Z-coordinates
        z_coords = []
        valid_chains = []

        for chain_id in chain_list:
            if chain_id in positions:
                z_coords.append(positions[chain_id][2])
                valid_chains.append(chain_id)

        if not z_coords:
            return {}

        # Cluster Z-coordinates
        z_array = np.array(z_coords).reshape(-1, 1)
        clustering = DBSCAN(
            eps=self.Z_CLUSTERING_EPS, min_samples=self.MIN_SAMPLES
        ).fit(z_array)

        # Map cluster labels to subunit indices
        unique_labels = sorted(list(set(clustering.labels_)))
        if -1 in unique_labels:  # Remove noise label
            unique_labels.remove(-1)

        # Sort labels by average Z-coordinate
        label_to_avg_z = {}
        for label in unique_labels:
            mask = clustering.labels_ == label
            avg_z = np.mean(z_array[mask])
            label_to_avg_z[label] = avg_z

        sorted_labels = sorted(unique_labels, key=lambda l: label_to_avg_z[l])
        label_to_subunit_idx = {label: i for i, label in enumerate(sorted_labels)}

        # Create chain to subunit index mapping
        chain_to_subunit_idx = {}
        for i, chain_id in enumerate(valid_chains):
            label = clustering.labels_[i]
            if label in label_to_subunit_idx:
                chain_to_subunit_idx[chain_id] = label_to_subunit_idx[label]
            else:
                chain_to_subunit_idx[chain_id] = -1  # Noise

        print(f"🎯 Clustered chains into {len(unique_labels)} Z-levels")
        return chain_to_subunit_idx

    def determine_structure_type(self, num_protofilaments: int) -> str:
        """Determine if structure is cylindrical or sheet-like"""
        if num_protofilaments > 4:
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

        print(
            f"🔍 Saved geometric debug data to debug_output/{pdb_id}_geometric_debug.json"
        )

    def process_protofilaments(
        self,
        protofilaments: List[List[str]],
        positions: Dict[str, np.ndarray],
        tubulin_chains: Dict[str, str],
    ) -> Tuple[List, str]:
        """Process protofilaments into final grid"""
        if not protofilaments:
            return [], "unknown"

        # Import here to avoid circular imports
        from main import SubunitData

        # Save debug data first
        pdb_id = (
            "current"  # This should be passed from caller, but keeping simple for now
        )
        self.save_geometric_debug(pdb_id, protofilaments, positions, tubulin_chains)

        # Align structure
        aligned_positions = self.align_structure_to_z_axis(positions, protofilaments)

        # Order protofilaments
        ordered_pfs = self.order_protofilaments_by_angle(
            protofilaments, aligned_positions
        )

        # Determine structure type
        structure_type = self.determine_structure_type(len(ordered_pfs))

        # Cluster all chains by Z-coordinate
        all_chains = [cid for pf in ordered_pfs for cid in pf]
        chain_to_subunit_idx = self.cluster_z_coordinates(aligned_positions, all_chains)

        # Create subunit data
        subunits = []
        seen_chains = set()  # Track to avoid duplicates

        for pf_idx, pf_chains in enumerate(ordered_pfs):
            for chain_id in pf_chains:
                if (
                    chain_id in chain_to_subunit_idx
                    and chain_id in tubulin_chains
                    and chain_id not in seen_chains
                ):  # Avoid duplicates

                    subunit_idx = chain_to_subunit_idx[chain_id]
                    if subunit_idx >= 0:  # Skip noise points
                        subunit = SubunitData(
                            id=f"pf{pf_idx}-{chain_id}",
                            auth_asym_id=chain_id,
                            protofilament=pf_idx,
                            subunitIndex=subunit_idx,
                            monomerType=tubulin_chains[chain_id],
                        )
                        subunits.append(subunit)
                        seen_chains.add(chain_id)

        print(
            f"📊 Created {len(subunits)} unique subunits from {len(seen_chains)} chains"
        )
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
