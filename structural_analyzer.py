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
from collections import defaultdict


class InterfaceAnalyzer:
    """Analyzes interfaces between protein chains"""

    def __init__(self, interface_cutoff: float = 7.0):  # Increased to 7Å for better N-terminus detection
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
                        print(f"🔗 N-term connection: {chain.id} → {next_chain} ({contact_count} contacts)")

        print(f"🔗 Found {len(connections)} N-terminus connections")
        
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
        
        print(f"🔍 Saved N-terminus debug data to debug_output/nterm_debug.json")

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
        """Get or download PDB profile data using GraphQL"""
        profile_path = self.profile_base_path / f"{pdb_id.upper()}_profile.json"
        if profile_path.exists():
            with open(profile_path, "r") as f:
                return json.load(f)

        # Use GraphQL query for comprehensive metadata
        profile_data = await self.generate_profile_graphql(pdb_id.upper())
        
        # Save to cache
        with open(profile_path, "w") as f:
            json.dump(profile_data, f, indent=4)
            
        return profile_data

    async def generate_profile_graphql(self, rcsb_id: str) -> Dict[str, Any]:
        """Generate raw structure profile using WholeStructureTemplate GraphQL query"""
        
        # GraphQL query template
        whole_structure_template = """{
          entry(entry_id: "$RCSB_ID") {
            assemblies{
                rcsb_id 
               	nonpolymer_entity_instances{
              
                  rcsb_nonpolymer_entity_instance_container_identifiers{
                    auth_asym_id
                    auth_seq_id
                    entity_id
                  }
                }
                polymer_entity_instances{
                   rcsb_polymer_entity_instance_container_identifiers {        
                    auth_asym_id
                    entity_id
                  }
                }
              }
            rcsb_id
              rcsb_accession_info{
                deposit_date
              }
              struct_keywords {
                pdbx_keywords
                text
              }
              rcsb_entry_info {
                resolution_combined
              }
              rcsb_external_references {
                link
                type
                id
              }
              exptl {
                method
              }
              citation {
                rcsb_authors
                year
                title
                pdbx_database_id_DOI
              }
              struct_keywords {
                pdbx_keywords
                text
              }
            polymer_entities {
              entry{
                rcsb_id
              }
              rcsb_polymer_entity_container_identifiers {
                asym_ids
                auth_asym_ids
                entry_id
                entity_id
              }
              pfams {
                rcsb_pfam_accession
                rcsb_pfam_comment
                rcsb_pfam_description
              }
              rcsb_entity_source_organism {
                ncbi_taxonomy_id
                scientific_name
              }
              rcsb_entity_host_organism {
                ncbi_taxonomy_id
                scientific_name
              }
              uniprots {
                rcsb_id
              }
              rcsb_polymer_entity {
                pdbx_description
              }
              entity_poly {
                pdbx_seq_one_letter_code
                pdbx_seq_one_letter_code_can
                pdbx_strand_id
                rcsb_entity_polymer_type
                rcsb_sample_sequence_length
                type
              }
            }
            nonpolymer_entities {
              pdbx_entity_nonpoly {
                comp_id
                name
                entity_id
              }
              rcsb_nonpolymer_entity {
                details
                formula_weight
                pdbx_description
                pdbx_number_of_molecules
              }
              nonpolymer_comp {
                chem_comp {
                  id
                  name
                  three_letter_code
                }
                rcsb_chem_comp_target {
                  interaction_type
                  name
                  provenance_source
                  reference_database_accession_code
                  reference_database_name
                }
                drugbank {
                  drugbank_container_identifiers {
                    drugbank_id
                  }
                  drugbank_info {
                    cas_number
                    description
                    indication
                    mechanism_of_action
                    name
                    pharmacology
                  }
                }
              }
            }
          }
        }"""
        
        # Replace template variable
        query = whole_structure_template.replace("$RCSB_ID", rcsb_id)
        
        # Execute GraphQL query
        raw_data = await self.query_rcsb_api(query)
        
        print(f"📊 Retrieved comprehensive profile for {rcsb_id} via GraphQL")
        return raw_data

    async def query_rcsb_api(self, query: str) -> Dict[str, Any]:
        """Execute GraphQL query against RCSB API"""
        url = "https://data.rcsb.org/graphql"
        
        payload = {
            "query": query
        }
        
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.post(
                url, 
                json=payload,
                headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()
            
            result = response.json()
            
            # Check for GraphQL errors
            if "errors" in result:
                error_msg = "; ".join(err.get("message", "Unknown error") for err in result["errors"])
                raise HTTPException(500, f"GraphQL query failed: {error_msg}")
            
            return result.get("data", {})

    def filter_tubulin_chains(self, profile: Dict[str, Any]) -> Dict[str, str]:
        """Identify alpha and beta tubulin chains from GraphQL profile"""
        chain_types: Dict[str, str] = {}
        
        # Handle GraphQL response structure
        entry_data = profile.get("entry", {})
        if not entry_data:
            print("⚠️ No entry data found in profile")
            return chain_types

        polymer_entities = entry_data.get("polymer_entities", [])
        if not polymer_entities:
            print("⚠️ No polymer entities found in profile")
            return chain_types

        print(f"🔍 Processing {len(polymer_entities)} polymer entities...")

        for entity in polymer_entities:
            # Get description from the GraphQL structure
            polymer_entity_info = entity.get("rcsb_polymer_entity", {})
            desc = polymer_entity_info.get("pdbx_description", "").lower()
            
            if "tubulin" in desc and ("alpha" in desc or "beta" in desc):
                m_type = "alpha" if "alpha" in desc else "beta"
                
                # Get chain IDs from GraphQL structure
                container_ids = entity.get("rcsb_polymer_entity_container_identifiers", {})
                auth_asym_ids = container_ids.get("auth_asym_ids", [])
                
                for chain_id in auth_asym_ids:
                    chain_types[chain_id] = m_type
                    
                print(f"🧬 Found {m_type}-tubulin entity: {desc}")
                print(f"   Chains: {auth_asym_ids}")

        print(
            f"🧬 Total tubulin chains found: {len(chain_types)} "
            f"({sum(1 for t in chain_types.values() if t == 'alpha')} alpha, "
            f"{sum(1 for t in chain_types.values() if t == 'beta')} beta)"
        )
        
        if not chain_types:
            print("❌ No tubulin chains found! Check if this structure contains tubulin.")
            print("🔍 Available polymer entity descriptions:")
            for i, entity in enumerate(polymer_entities[:5]):  # Show first 5
                polymer_info = entity.get("rcsb_polymer_entity", {})
                desc = polymer_info.get("pdbx_description", "No description")
                print(f"   Entity {i}: {desc}")
        
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

    def trace_protofilaments_from_connections(
        self, nterm_connections: Dict[str, str], tubulin_chains: Dict[str, str]
    ) -> List[List[str]]:
        """Build full protofilaments by following N-terminus connections"""
        
        if not nterm_connections:
            print("❌ No N-term connections to trace")
            return []
        
        print(f"🔍 Building protofilaments from {len(nterm_connections)} connections")
        
        # Build reverse mapping: target -> source
        incoming = {target: source for source, target in nterm_connections.items()}
        
        visited = set()
        protofilaments = []
        
        # Find TRUE start chains (chains that are not targets of any connection)
        all_chains_in_connections = set(nterm_connections.keys()) | set(nterm_connections.values())
        true_start_chains = []
        
        # A true start is a chain that:
        # 1. Has outgoing connections (is a source) AND
        # 2. Is NOT a target of any connection
        for chain in nterm_connections.keys():  # All sources
            if chain not in incoming:  # Not a target of anyone
                true_start_chains.append(chain)
        
        print(f"🎯 Found {len(true_start_chains)} TRUE protofilament starts: {true_start_chains[:10]}...")
        
        # If we have very few true starts, we might have cycles or missing connections
        # In that case, also consider high-numbered chains as starts
        expected_pfs = max(8, len(nterm_connections) // 4)  # Estimate based on connections
        if len(true_start_chains) < expected_pfs // 2:  # Less than half expected
            print(f"⚠️ Too few true starts ({len(true_start_chains)}), looking for backup starts...")
            backup_starts = []
            for chain in all_chains_in_connections:
                # Look for chains that could be starts (higher numbered layers or specific patterns)
                if any(chain.startswith(prefix) for prefix in ['3', '4']) and chain not in visited:
                    backup_starts.append(chain)
            
            print(f"🔄 Adding {len(backup_starts)} backup starts: {backup_starts[:10]}...")
            true_start_chains.extend(backup_starts)
        
        # Trace each protofilament from its start
        for start_chain in true_start_chains:
            if start_chain in visited:
                continue
                
            # Follow the chain of connections
            protofilament = []
            current = start_chain
            pf_visited = set()  # Local visited set to detect cycles
            
            while current and current not in visited:
                if current in pf_visited:  # Cycle detection
                    print(f"  ⚠️ Cycle detected at {current}, breaking")
                    break
                    
                protofilament.append(current)
                visited.add(current)
                pf_visited.add(current)
                print(f"  → Adding {current} to PF{len(protofilaments)}")
                
                # Move to next chain via N-terminus connection
                current = nterm_connections.get(current)
            
            if len(protofilament) >= 2:  # Only keep substantial protofilaments
                protofilaments.append(protofilament)
                print(f"✅ Completed PF{len(protofilaments)-1}: {len(protofilament)} chains")
            elif len(protofilament) == 1:
                # Single chains still count as protofilaments
                protofilaments.append(protofilament) 
                print(f"⚠️ Single-chain PF{len(protofilaments)-1}: {protofilament}")
        
        # Handle remaining connected chains (might be in cycles)
        remaining_chains = all_chains_in_connections - visited
        if remaining_chains:
            print(f"🔍 Processing {len(remaining_chains)} remaining chains...")
            
            for chain in remaining_chains:
                if chain in visited:
                    continue
                    
                # Try to build a protofilament from this chain
                protofilament = []
                current = chain
                seen_in_this_pf = set()
                
                while current and current not in visited and current not in seen_in_this_pf:
                    protofilament.append(current)
                    visited.add(current)
                    seen_in_this_pf.add(current)
                    current = nterm_connections.get(current)
                
                if len(protofilament) >= 1:
                    protofilaments.append(protofilament)
                    print(f"✅ Found remaining PF: {len(protofilament)} chains")
        
        # Add isolated chains (chains not in any connections) - but group them smarter
        all_tubulin_chains = set(tubulin_chains.keys())
        isolated_chains = all_tubulin_chains - all_chains_in_connections
        
        if isolated_chains:
            print(f"🏝️ Found {len(isolated_chains)} isolated chains - grouping by layer/position...")
            
            # Group isolated chains by their layer number (1A, 1B, 1C vs 2A, 2B, 2C, etc.)
            isolated_by_layer = {}
            for chain in isolated_chains:
                layer = chain[0] if chain else "unknown"  # First character (1,2,3,4)
                if layer not in isolated_by_layer:
                    isolated_by_layer[layer] = []
                isolated_by_layer[layer].append(chain)
            
            # Add each layer as separate protofilaments, but don't create too many tiny ones
            for layer, chains in isolated_by_layer.items():
                if len(chains) <= 3:  # Small groups get individual protofilaments
                    for chain in chains:
                        protofilaments.append([chain])
                else:  # Large groups get chunked
                    chunk_size = max(2, len(chains) // 3)  # Try to make ~3 protofilaments per layer
                    for i in range(0, len(chains), chunk_size):
                        chunk = chains[i:i+chunk_size]
                        protofilaments.append(chunk)
                        print(f"✅ Isolated layer {layer} PF: {len(chunk)} chains")
        
        print(f"🧵 Built {len(protofilaments)} protofilaments:")
        for i, pf in enumerate(protofilaments):
            alpha_count = sum(1 for c in pf if tubulin_chains.get(c) == 'alpha')
            beta_count = sum(1 for c in pf if tubulin_chains.get(c) == 'beta')
            print(f"  PF{i}: {len(pf)} chains ({alpha_count}α, {beta_count}β)")
        
        return protofilaments

    def finalize_grid_from_protofilaments(
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

        # Find N-terminus connections
        nterm_connections = self.interface_analyzer.find_nterm_connections(
            structure, tubulin_chains
        )

        # Trace complete protofilaments by following connections
        protofilaments = self.trace_protofilaments_from_connections(
            nterm_connections, tubulin_chains
        )

        # Create final grid
        subunits, structure_type = self.finalize_grid_from_protofilaments(
            protofilaments, positions, tubulin_chains
        )

        # Debug output
        alpha_count = sum(1 for subunit in subunits if subunit.monomerType == "alpha")
        beta_count = sum(1 for subunit in subunits if subunit.monomerType == "beta")
        print(
            f"🔍 Final grid: {alpha_count} alpha subunits, {beta_count} beta subunits"
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
        self.visualizer.create_protofilament_tracing_plot(
            positions, protofilaments, tubulin_chains, pdb_id
        )
        self.visualizer.create_grid_visualization(grid_data, pdb_id)

        # Cleanup
        mmcif_path.unlink()

        print(f"🎉 Grid generation completed for {pdb_id.upper()}")
        return grid_data
