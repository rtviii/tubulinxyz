import numpy as np
import asyncio
import httpx
from pathlib import Path
import tempfile
import json
import os
from typing import Dict, List, Any, Tuple, Set
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from fastapi import HTTPException
from geometric_operations import GeometricAnalyzer
from visualization_utils import VisualizationUtils


class InterfaceAnalyzer:
    """Analyzes interfaces between protein chains"""

    def __init__(self, interface_cutoff: float = 10.0):  # Increased to 10Å
        self.INTERFACE_CUTOFF = interface_cutoff

    def find_chain_neighbors(
        self, chain_positions: Dict[str, np.ndarray], neighbor_cutoff: float = 75.0
    ) -> Dict[str, List[str]]:
        """Find neighboring chains based on distance cutoff"""
        neighbors = {chain_id: [] for chain_id in chain_positions.keys()}

        chain_ids = list(chain_positions.keys())
        for i, chain_id1 in enumerate(chain_ids):
            for j, chain_id2 in enumerate(chain_ids[i + 1 :], i + 1):
                dist = np.linalg.norm(
                    chain_positions[chain_id1] - chain_positions[chain_id2]
                )
                if dist <= neighbor_cutoff:
                    neighbors[chain_id1].append(chain_id2)
                    neighbors[chain_id2].append(chain_id1)

        print(f"🔗 Found chain neighbors for {len(neighbors)} chains")
        return neighbors

    def find_nterm_connections(
        self, structure, tubulin_chains: Dict[str, str]
    ) -> Dict[str, str]:
        """Find protofilament connections using N-terminus to next-monomer interfaces"""
        connections = {}  # chain_id -> next_chain_id in protofilament

        # Get all atoms for neighbor search
        all_atoms = []
        chain_to_atoms = {}

        for model in structure:
            for chain in model:
                if chain.id in tubulin_chains:
                    chain_atoms = []
                    for residue in chain:
                        for atom in residue:
                            all_atoms.append(atom)
                            chain_atoms.append(atom)
                    chain_to_atoms[chain.id] = chain_atoms

        # Create neighbor search object
        neighbor_search = NeighborSearch(all_atoms)

        # For each chain, find its N-terminus and what it connects to
        for model in structure:
            for chain in model:
                if chain.id not in tubulin_chains:
                    continue

                # Get N-terminus residue (first by sequence number)
                residues = list(chain.get_residues())
                if not residues:
                    continue

                residues.sort(key=lambda r: r.get_id()[1])
                n_terminus = residues[0]

                # Find atoms from OTHER chains within 5Å of N-terminus
                nearby_chains = {}  # chain_id -> count of close atoms

                for atom in n_terminus:
                    nearby_atoms = neighbor_search.search(atom.coord, 5.0)

                    for nearby_atom in nearby_atoms:
                        nearby_chain_id = nearby_atom.get_parent().get_parent().id

                        # Skip same chain and non-tubulin chains
                        if (
                            nearby_chain_id != chain.id
                            and nearby_chain_id in tubulin_chains
                        ):

                            if nearby_chain_id not in nearby_chains:
                                nearby_chains[nearby_chain_id] = 0
                            nearby_chains[nearby_chain_id] += 1

                # The chain with the most contacts to our N-terminus is the "next" one
                if nearby_chains:
                    next_chain = max(nearby_chains, key=nearby_chains.get)
                    contact_count = nearby_chains[next_chain]

                    if contact_count >= 3:  # Minimum meaningful contact
                        connections[chain.id] = next_chain
                        print(
                            f"🔗 N-term connection: {chain.id} → {next_chain} ({contact_count} contacts)"
                        )

        print(f"🔗 Found {len(connections)} N-terminus connections")
        return connections

    def get_n_terminus_residue(self, structure, chain_id: str) -> Residue:
        """Get the N-terminus residue of a chain"""
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    residues = list(chain.get_residues())
                    if residues:
                        # Sort by residue number to ensure we get the first one
                        residues.sort(key=lambda r: r.get_id()[1])
                        return residues[0]
        return None


class SpatialGridGenerator:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
        self.profile_base_path = Path(os.getenv("PROFILE_BASE_PATH", "profiles"))
        self.profile_base_path.mkdir(exist_ok=True)
        self.NEIGHBOR_CUTOFF_Å = 75.0
        self.interface_analyzer = InterfaceAnalyzer()
        self.geometric_analyzer = GeometricAnalyzer()
        self.visualizer = VisualizationUtils()

    async def get_profile(self, pdb_id: str) -> Dict[str, Any]:
        """Get or download PDB profile data"""
        profile_path = self.profile_base_path / f"{pdb_id.upper()}_profile.json"
        if profile_path.exists():
            with open(profile_path, "r") as f:
                return json.load(f)

        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        async with httpx.AsyncClient() as client:
            response = await client.get(url)
            response.raise_for_status()
            profile_data = response.json()
            with open(profile_path, "w") as f:
                json.dump(profile_data, f)
            return profile_data

    def filter_tubulin_chains(self, profile: Dict[str, Any]) -> Dict[str, str]:
        """Identify alpha and beta tubulin chains from profile"""
        chain_types: Dict[str, str] = {}
        entry_data = profile.get("entry", {})

        for entity in entry_data.get("polymer_entities", []):
            desc = (
                entity.get("rcsb_polymer_entity", {})
                .get("pdbx_description", "")
                .lower()
            )
            if "tubulin" in desc and ("alpha" in desc or "beta" in desc):
                m_type = "alpha" if "alpha" in desc else "beta"
                ids = entity.get("rcsb_polymer_entity_container_identifiers", {}).get(
                    "auth_asym_ids", []
                )
                for chain_id in ids:
                    chain_types[chain_id] = m_type

        print(
            f"🧬 Found {len(chain_types)} tubulin chains: {sum(1 for t in chain_types.values() if t == 'alpha')} alpha, {sum(1 for t in chain_types.values() if t == 'beta')} beta"
        )
        return chain_types

    async def download_mmcif(self, pdb_id: str) -> Path:
        """Download mmCIF file"""
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.get(url)
            response.raise_for_status()
            temp_file = Path(tempfile.mktemp(suffix=".cif"))
            temp_file.write_bytes(response.content)
            return temp_file

    def extract_chain_positions(
        self, structure, tubulin_chains: Dict[str, str]
    ) -> Dict[str, np.ndarray]:
        """Extract centroid positions for each chain"""
        positions = {}

        for model in structure:
            for chain in model:
                if chain.id not in tubulin_chains:
                    continue

                # Get all atoms in chain
                all_atoms = []
                for residue in chain:
                    for atom in residue:
                        all_atoms.append(atom.get_coord())

                if all_atoms:
                    positions[chain.id] = np.mean(all_atoms, axis=0)

        print(f"📍 Extracted positions for {len(positions)} chains")
        return positions

    def determine_chain_directionality(
        self,
        structure,
        interfaces: Dict[str, Dict[str, List[Residue]]],
        tubulin_chains: Dict[str, str],
    ) -> Dict[str, str]:
        """Determine directionality using N-terminus interface analysis"""
        chain_directions = {}

        for chain_id in tubulin_chains.keys():
            if chain_id not in interfaces:
                continue

            # Get N-terminus residue
            n_terminus = self.interface_analyzer.get_n_terminus_residue(
                structure, chain_id
            )
            if not n_terminus:
                continue

            # Check which neighbor chains have the N-terminus in their interface
            n_terminus_neighbors = []
            for neighbor_id, interface_residues in interfaces[chain_id].items():
                if n_terminus in interface_residues:
                    n_terminus_neighbors.append(neighbor_id)

            # Store directionality information
            chain_directions[chain_id] = {
                "n_terminus_residue": n_terminus,
                "n_terminus_neighbors": n_terminus_neighbors,
            }

        print(f"🧭 Determined directionality for {len(chain_directions)} chains")
        return chain_directions

    def build_dimer_pairs(
        self,
        interfaces: Dict[str, Dict[str, List[Residue]]],
        tubulin_chains: Dict[str, str],
        positions: Dict[str, np.ndarray],
    ) -> Dict[str, str]:
        """Build dimer pairs using interface analysis"""
        dimer_pairs = {}

        # For each alpha chain, find the best beta partner
        alpha_chains = [
            cid for cid, mtype in tubulin_chains.items() if mtype == "alpha"
        ]

        print(
            f"🔍 Debug: Processing {len(alpha_chains)} alpha chains for dimer pairing"
        )

        for alpha_id in alpha_chains:
            if alpha_id not in interfaces:
                print(f"⚠️ Alpha {alpha_id} has no interfaces")
                continue

            best_beta = None
            max_interface_size = 0

            print(
                f"🔍 Alpha {alpha_id} has interfaces with: {list(interfaces[alpha_id].keys())}"
            )

            # Check all beta neighbors
            for neighbor_id, interface_residues in interfaces[alpha_id].items():
                neighbor_type = tubulin_chains.get(neighbor_id, "unknown")
                interface_size = len(interface_residues)

                print(
                    f"  - Neighbor {neighbor_id} ({neighbor_type}): {interface_size} interface residues"
                )

                if neighbor_type == "beta":
                    if interface_size > max_interface_size:
                        max_interface_size = interface_size
                        best_beta = neighbor_id
                        print(
                            f"    ✓ New best beta partner: {best_beta} with {interface_size} residues"
                        )

            if best_beta and max_interface_size > 5:  # Minimum interface size threshold
                dimer_pairs[alpha_id] = best_beta
                dimer_pairs[best_beta] = alpha_id
                print(
                    f"✅ Paired α{alpha_id} ↔ β{best_beta} ({max_interface_size} interface residues)"
                )
            else:
                print(
                    f"❌ No suitable beta partner for α{alpha_id} (best had {max_interface_size} residues)"
                )

        print(f"🤝 Built {len(dimer_pairs)//2} dimer pairs")
        return dimer_pairs

    def build_dimer_pairs_simple(
        self, positions: Dict[str, np.ndarray], tubulin_chains: Dict[str, str]
    ) -> Dict[str, str]:
        """Build dimer pairs using simple distance-based approach"""
        dimer_pairs = {}

        alpha_chains = [
            cid
            for cid, mtype in tubulin_chains.items()
            if mtype == "alpha" and cid in positions
        ]
        beta_chains = [
            cid
            for cid, mtype in tubulin_chains.items()
            if mtype == "beta" and cid in positions
        ]

        print(
            f"🔍 Simple pairing: {len(alpha_chains)} alphas, {len(beta_chains)} betas"
        )

        for alpha_id in alpha_chains:
            alpha_pos = positions[alpha_id]

            # Find closest beta
            min_dist = float("inf")
            closest_beta = None

            for beta_id in beta_chains:
                if beta_id in dimer_pairs:  # Already paired
                    continue

                beta_pos = positions[beta_id]
                dist = np.linalg.norm(alpha_pos - beta_pos)

                if dist < min_dist:
                    min_dist = dist
                    closest_beta = beta_id

            # Pair if distance is reasonable (< 50Å for dimers in this large structure)
            if closest_beta and min_dist < 50.0:
                dimer_pairs[alpha_id] = closest_beta
                dimer_pairs[closest_beta] = alpha_id
                print(
                    f"✅ Simple pair: α{alpha_id} ↔ β{closest_beta} ({min_dist:.1f}Å)"
                )

        print(f"🤝 Built {len(dimer_pairs)//2} simple dimer pairs")
        return dimer_pairs

    def trace_protofilaments_simple(
        self, nterm_connections: Dict[str, str], tubulin_chains: Dict[str, str]
    ) -> List[List[str]]:
        """Trace protofilaments using simple linked list approach"""

        # Build linked lists for each protofilament
        protofilaments = []  # List of chain lists
        chain_to_pf = {}  # chain_id -> protofilament_index

        print(f"🔍 Building protofilaments from {len(nterm_connections)} connections")

        # Process each N-term connection
        for source_chain, target_chain in nterm_connections.items():
            print(f"🔗 Processing connection: {source_chain} → {target_chain}")

            source_pf = chain_to_pf.get(source_chain)
            target_pf = chain_to_pf.get(target_chain)

            if source_pf is None and target_pf is None:
                # Create new protofilament
                new_pf_idx = len(protofilaments)
                protofilaments.append([source_chain, target_chain])
                chain_to_pf[source_chain] = new_pf_idx
                chain_to_pf[target_chain] = new_pf_idx
                print(
                    f"  ✨ Created new protofilament {new_pf_idx}: [{source_chain}, {target_chain}]"
                )

            elif source_pf is None:
                # Add source to existing target protofilament
                protofilaments[target_pf].insert(0, source_chain)  # Add to front
                chain_to_pf[source_chain] = target_pf
                print(
                    f"  ➕ Added {source_chain} to front of protofilament {target_pf}"
                )

            elif target_pf is None:
                # Add target to existing source protofilament
                protofilaments[source_pf].append(target_chain)  # Add to end
                chain_to_pf[target_chain] = source_pf
                print(f"  ➕ Added {target_chain} to end of protofilament {source_pf}")

            elif source_pf != target_pf:
                # Merge two protofilaments
                source_pf_chains = protofilaments[source_pf]
                target_pf_chains = protofilaments[target_pf]

                # Merge target into source
                merged_chains = source_pf_chains + target_pf_chains
                protofilaments[source_pf] = merged_chains

                # Update chain mappings
                for chain in target_pf_chains:
                    chain_to_pf[chain] = source_pf

                # Remove the target protofilament (mark as empty)
                protofilaments[target_pf] = []
                print(f"  🔗 Merged protofilament {target_pf} into {source_pf}")

            # else: both already in same protofilament, nothing to do

        # Filter out empty protofilaments
        final_protofilaments = []
        for pf_chains in protofilaments:
            if len(pf_chains) >= 2:  # At least 2 chains
                final_protofilaments.append(pf_chains)
                print(
                    f"✅ Final protofilament with {len(pf_chains)} chains: {pf_chains[:10]}..."
                )

        print(f"🧵 Built {len(final_protofilaments)} protofilaments")
        return final_protofilaments

    def finalize_grid_from_interfaces(
        self,
        protofilaments: List[List[str]],
        positions: Dict[str, np.ndarray],
        tubulin_chains: Dict[str, str],
    ) -> Tuple[List, str]:
        """Create final grid using geometric analysis"""
        if not protofilaments:
            return [], "unknown"

        # Use geometric analyzer to process the protofilaments
        subunits, structure_type = self.geometric_analyzer.process_protofilaments(
            protofilaments, positions, tubulin_chains
        )

        print(
            f"✅ Finalized grid with {len(subunits)} subunits. Structure type: {structure_type.upper()}"
        )
        return subunits, structure_type

    async def generate_grid(self, pdb_id: str):
        """Main method to generate grid from PDB structure"""
        print(f"🚀 Starting grid generation for {pdb_id.upper()}")

        # Get profile and identify tubulin chains
        profile = await self.get_profile(pdb_id)
        tubulin_chains = self.filter_tubulin_chains(profile)
        if not tubulin_chains:
            raise HTTPException(400, f"No tubulin chains found for {pdb_id}")

        # Download and parse structure
        mmcif_path = await self.download_mmcif(pdb_id)
        structure = self.parser.get_structure(pdb_id, mmcif_path)

        # Extract positions
        positions = self.extract_chain_positions(structure, tubulin_chains)
        if not positions:
            raise HTTPException(400, "Could not extract chain positions")

        # Find chain neighbors
        chain_neighbors = self.interface_analyzer.find_chain_neighbors(
            positions, self.NEIGHBOR_CUTOFF_Å
        )

        # Find N-terminus connections (simpler approach)
        nterm_connections = self.interface_analyzer.find_nterm_connections(
            structure, tubulin_chains
        )

        # Trace protofilaments using N-terminus connections
        protofilaments = self.trace_protofilaments_simple(
            nterm_connections, tubulin_chains
        )

        # Create final grid
        subunits, structure_type = self.finalize_grid_from_interfaces(
            protofilaments, positions, tubulin_chains
        )

        # Debug output
        alpha_count = sum(1 for subunit in subunits if subunit.monomerType == "alpha")
        beta_count = sum(1 for subunit in subunits if subunit.monomerType == "beta")
        print(
            f"🔍 Debug: {alpha_count} alpha subunits, {beta_count} beta subunits in final grid"
        )

        # Create visualization
        from main import GridData

        grid_data = GridData(
            subunits=subunits,
            structure_type=structure_type,
            metadata={
                "pdb_id": pdb_id,
                "num_tubulin_chains": len(positions),
                "num_protofilaments": len(protofilaments),
                "nterm_connections": len(nterm_connections),
            },
        )

        # Generate visualizations
        # Create simple debug data for visualization
        debug_data = {
            "nterm_connections": nterm_connections,
            "protofilaments": protofilaments,
        }

        self.visualizer.create_debug_visualization(
            debug_data, pdb_id, "N-term Connections"
        )
        self.visualizer.create_protofilament_tracing_plot(
            positions, protofilaments, tubulin_chains, pdb_id
        )
        self.visualizer.create_grid_visualization(grid_data, pdb_id)

        # Cleanup
        mmcif_path.unlink()

        print(f"🎉 Grid generation completed for {pdb_id.upper()}")
        return grid_data
