import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Any
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from collections import defaultdict


class InterfaceAnalyzer:
    """Analyzes interfaces between protein chains"""

    def __init__(self, interface_cutoff: float = 7.0):
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

        return neighbors

    def find_nterm_connections(
        self, structure, tubulin_chains: Dict[str, str]
    ) -> Dict[str, str]:
        """Find protofilament connections using N-terminus to next-monomer interfaces"""
        connections = {}  # chain_id -> next_chain_id in protofilament
        debug_data = {}

        # Get all atoms for neighbor search
        all_atoms = []
        for model in structure:
            for chain in model:
                if chain.id in tubulin_chains:
                    for residue in chain:
                        for atom in residue:
                            all_atoms.append(atom)

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

                # Find atoms from OTHER chains within cutoff of N-terminus
                nearby_chains = {}  # chain_id -> count of close atoms
                nearby_residues = defaultdict(set)  # chain_id -> set of residues

                for atom in n_terminus:
                    nearby_atoms = neighbor_search.search(atom.coord, self.INTERFACE_CUTOFF)

                    for nearby_atom in nearby_atoms:
                        nearby_residue = nearby_atom.get_parent()
                        nearby_chain_id = nearby_residue.get_parent().id

                        # Skip same chain and non-tubulin chains
                        if (nearby_chain_id != chain.id and 
                            nearby_chain_id in tubulin_chains):
                            
                            if nearby_chain_id not in nearby_chains:
                                nearby_chains[nearby_chain_id] = 0
                            nearby_chains[nearby_chain_id] += 1
                            nearby_residues[nearby_chain_id].add(nearby_residue)

                # Store debug data
                debug_data[chain.id] = {
                    "n_terminus_resnum": n_terminus.get_id()[1],
                    "neighbor_counts": dict(nearby_chains),
                    "neighbor_residues": {k: len(v) for k, v in nearby_residues.items()}
                }

                # The chain with the most contacts to our N-terminus is the "next" one
                if nearby_chains:
                    next_chain = max(nearby_chains, key=nearby_chains.get)
                    contact_count = nearby_chains[next_chain]

                    if contact_count >= 3:  # Minimum meaningful contact
                        connections[chain.id] = next_chain

        print(f"Found {len(connections)} N-terminus connections")
        
        # Save debug data
        self._save_nterm_debug(debug_data, connections)
        
        return connections

    def _save_nterm_debug(self, debug_data: Dict, connections: Dict):
        """Save N-terminus debug data for PyMOL visualization"""
        debug_dir = Path("debug_output")
        debug_dir.mkdir(exist_ok=True)
        
        pymol_data = {
            "n_terminus_data": debug_data,
            "connections": connections,
            "total_chains": len(debug_data),
            "successful_connections": len(connections)
        }
        
        with open(debug_dir / "nterm_debug.json", 'w') as f:
            json.dump(pymol_data, f, indent=2)

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