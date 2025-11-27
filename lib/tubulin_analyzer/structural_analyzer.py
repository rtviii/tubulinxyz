import numpy as np
import asyncio
import httpx
from pathlib import Path
import tempfile
import json
import os
from typing import Dict, List, Any, Tuple, Set
from Bio.PDB.MMCIFParser import MMCIFParser
from fastapi import HTTPException

from lib.tubulin_analyzer.interface_analyzer import InterfaceAnalyzer
from lib.tubulin_analyzer.geometric_analyzer import GeometricAnalyzer
from lib.tubulin_analyzer.visualization import VisualizationUtils
from lib.tubulin_analyzer.models import GridData, SubunitData


class SpatialGridGenerator:
    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)
        self.profile_base_path = Path(os.getenv("PROFILE_BASE_PATH", "profiles"))
        self.profile_base_path.mkdir(exist_ok=True)
        self.NEIGHBOR_CUTOFF_Å  = 75.0
        self.interface_analyzer = InterfaceAnalyzer()
        self.geometric_analyzer = GeometricAnalyzer()
        self.visualizer         = VisualizationUtils()

    async def get_profile(self, pdb_id: str) -> Dict[str, Any]:

        """Get or download PDB profile data using GraphQL"""
        profile_path = self.profile_base_path / f"{pdb_id.upper()}.json"
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
            print("Warning: No entry data found in profile")
            return chain_types

        polymer_entities = entry_data.get("polymer_entities", [])
        if not polymer_entities:
            print("Warning: No polymer entities found in profile")
            return chain_types

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
        
        if not chain_types:
            print("Error: No tubulin chains found! Check if this structure contains tubulin.")
        
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

        return positions

    def trace_protofilaments_from_connections(
        self, nterm_connections: Dict[str, str], tubulin_chains: Dict[str, str]
    ) -> List[List[str]]:
        """Build full protofilaments by following N-terminus connections"""
        
        if not nterm_connections:
            print("Error: No N-term connections to trace")
            return []
        
        # Build reverse mapping: target -> source
        incoming = {target: source for source, target in nterm_connections.items()}
        
        visited = set()
        protofilaments = []
        
        # Find TRUE start chains (chains that are not targets of any connection)
        all_chains_in_connections = set(nterm_connections.keys()) | set(nterm_connections.values())
        true_start_chains = []
        
        # A true start is a chain that has outgoing connections AND is NOT a target
        for chain in nterm_connections.keys():
            if chain not in incoming:
                true_start_chains.append(chain)
        
        # If we have very few true starts, add backup starts
        expected_pfs = max(8, len(nterm_connections) // 4)
        if len(true_start_chains) < expected_pfs // 2:
            backup_starts = []
            for chain in all_chains_in_connections:
                if any(chain.startswith(prefix) for prefix in ['3', '4']) and chain not in visited:
                    backup_starts.append(chain)
            true_start_chains.extend(backup_starts)
        
        # Trace each protofilament from its start
        for start_chain in true_start_chains:
            if start_chain in visited:
                continue
                
            # Follow the chain of connections
            protofilament = []
            current = start_chain
            pf_visited = set()
            
            while current and current not in visited:
                if current in pf_visited:  # Cycle detection
                    break
                    
                protofilament.append(current)
                visited.add(current)
                pf_visited.add(current)
                
                # Move to next chain via N-terminus connection
                current = nterm_connections.get(current)
            
            if len(protofilament) >= 1:
                protofilaments.append(protofilament)
        
        # Handle remaining connected chains (might be in cycles)
        remaining_chains = all_chains_in_connections - visited
        if remaining_chains:
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
        
        # Add isolated chains (chains not in any connections)
        all_tubulin_chains = set(tubulin_chains.keys())
        isolated_chains = all_tubulin_chains - all_chains_in_connections
        
        if isolated_chains:
            # Group isolated chains by their layer number
            isolated_by_layer = {}
            for chain in isolated_chains:
                layer = chain[0] if chain else "unknown"
                if layer not in isolated_by_layer:
                    isolated_by_layer[layer] = []
                isolated_by_layer[layer].append(chain)
            
            # Add each layer as separate protofilaments
            for layer, chains in isolated_by_layer.items():
                if len(chains) <= 3:
                    for chain in chains:
                        protofilaments.append([chain])
                else:
                    chunk_size = max(2, len(chains) // 3)
                    for i in range(0, len(chains), chunk_size):
                        chunk = chains[i:i+chunk_size]
                        protofilaments.append(chunk)
        
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

        return subunits, structure_type

    async def generate_grid(self, pdb_id: str):
        """Main method to generate grid from PDB structure"""
        print(f"Starting grid generation for {pdb_id.upper()}")

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

        # Create visualization
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

        print(f"Grid generation completed for {pdb_id.upper()}")
        return grid_data