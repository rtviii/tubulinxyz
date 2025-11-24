import asyncio
import json
import os
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
from datetime import datetime
from pydantic import ValidationError

from api.config import PROJECT_ROOT
from api.main import MASTER_PROFILE, MUSCLE_BINARY
from lib.etl.assets import AssemblyIdentificationString, EntryInfoString, LigandsChemInfo, NonpolymerEntitiesString, PolymerEntitiesString, TubulinStructureAssets, query_rcsb_api
from lib.models.types_tubulin import (
    TubulinStructure,
    TubulinProtein,
    Polymer,
    NonpolymericLigand,
    AssemblyInstancesMap,
    TubulinFamily
)

from lib.etl.constants import TUBETL_DATA

class TubulinETLCollector:
    """
    - Fetches data from RCSB using multiple GQL queries.
    - Parses and validates data against tubexyz Pydantic models.
    - Classifies tubulin proteins by family (Alpha, Beta, etc.)
    - Saves the final TubulinStructure profile to disk.
    """

    rcsb_id : str
    assets  : TubulinStructureAssets
    asm_maps: List[AssemblyInstancesMap] = []

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id.upper()
        self.assets  = TubulinStructureAssets(self.rcsb_id)

    def _classify_tubulin_family(self, polymer_data: dict) -> Optional[TubulinFamily]:
            """
            Simple classifier to determine TubulinFamily based on pdbx_description,
            TODO: This is a stand-in for a more robust classifier.
            """
            
            # Get the description string, make it lowercase, and handle None
            desc = (polymer_data
                    .get("rcsb_polymer_entity", {})
                    .get("pdbx_description", "") 
                    or "").lower()

            if "alpha" in desc or "tuba" in desc:
                return TubulinFamily.ALPHA
            if "beta" in desc or "tubb" in desc:
                return TubulinFamily.BETA
            if "gamma" in desc or "tubg" in desc:
                return TubulinFamily.GAMMA

            return None

    def _raw_to_polymer_base(self, 
                                poly_entity: dict, 
                                auth_asym_id: str, 
                                assembly_id: int) -> Polymer:
            """
            Helper to create a single, instance-specific Polymer model.
            """
            src_orgs = poly_entity.get("rcsb_entity_source_organism", []) or []
            host_orgs = poly_entity.get("rcsb_entity_host_organism", []) or []
            identifiers = poly_entity["rcsb_polymer_entity_container_identifiers"]

            return Polymer(
                entity_id      = identifiers["entity_id"],
                parent_rcsb_id = self.rcsb_id,
                auth_asym_id   = auth_asym_id,
                assembly_id    = assembly_id,
                asym_ids=identifiers["asym_ids"], # Full list of non-auth IDs
                
                src_organism_names                  = [org.get("scientific_name") for org in src_orgs if org.get("scientific_name")],
                host_organism_names                 = [org.get("scientific_name") for org in host_orgs if org.get("scientific_name")],
                src_organism_ids                    = [org.get("ncbi_taxonomy_id") for org in src_orgs if org.get("ncbi_taxonomy_id")],
                host_organism_ids                   = [org.get("ncbi_taxonomy_id") for org in host_orgs if org.get("ncbi_taxonomy_id")],
                rcsb_pdbx_description               = poly_entity.get("rcsb_polymer_entity", {}).get("pdbx_description"),
                entity_poly_strand_id               = poly_entity["entity_poly"]["pdbx_strand_id"],
                entity_poly_seq_one_letter_code     = poly_entity["entity_poly"]["pdbx_seq_one_letter_code"],
                entity_poly_seq_one_letter_code_can = poly_entity["entity_poly"]["pdbx_seq_one_letter_code_can"],
                entity_poly_seq_length              = poly_entity["entity_poly"]["rcsb_sample_sequence_length"],
                entity_poly_polymer_type            = poly_entity["entity_poly"]["rcsb_entity_polymer_type"],
                entity_poly_entity_type             = poly_entity["entity_poly"]["type"],
            )

    def _process_polymers(self, polymers_data: dict) -> Tuple[List[TubulinProtein], List[Polymer]]:
            """
            Process polymer entities, creating one object per INSTANCE.
            """
            tubulin_proteins: List[TubulinProtein] = []
            other_polymers: List[Polymer] = []
            
            print(f"    [Debug] _process_polymers received polymers_data: {'Yes' if polymers_data else 'No'}")
            
            if not polymers_data:
                return [], []

            # Safely get the list of polymer entities, defaulting to []
            polymer_entities_list = polymers_data.get("polymer_entities") or []
            
            print(f"    [Debug] Found {len(polymer_entities_list)} polymer entities.")

            if not polymer_entities_list:
                return [], []

            # Get *all* auth_asym_ids from *all* entities first
            all_auth_asym_ids = []
            for poly in polymer_entities_list:
                identifiers = poly.get("rcsb_polymer_entity_container_identifiers") or {}
                # Safely get auth_asym_ids, defaulting to []
                auth_ids = identifiers.get("auth_asym_ids") or []
                all_auth_asym_ids.extend(auth_ids)
            
            # Build a single, efficient map of instance -> assembly
            assembly_map = self._get_assembly_mappings(all_auth_asym_ids)

            for i, poly_entity in enumerate(polymer_entities_list):
                
                entity_desc_dict = poly_entity.get("rcsb_polymer_entity") or {}
                entity_desc = entity_desc_dict.get("pdbx_description", "N/A")
                print(f"    [Debug] Processing entity {i+1}/{len(polymer_entities_list)}: {entity_desc[:50]}...")
                
                # Get shared entity data
                family = self._classify_tubulin_family(poly_entity)
                
                # --- TUBE-FIX: REMOVED ALL PFAM LOGIC ---
                uniprots_list = poly_entity.get("uniprots") or []
                uniprot_accessions = [up.get("rcsb_id") for up in uniprots_list if up.get("rcsb_id")]
                # --- END FIX ---

                # Safely get the list of auth_asym_ids, defaulting to an empty list
                identifiers = poly_entity.get("rcsb_polymer_entity_container_identifiers") or {}
                auth_asym_ids = identifiers.get("auth_asym_ids") or []
                
                if not auth_asym_ids:
                    print(f"[Debug] ---> SKIPPING entity: No 'auth_asym_ids' (instances) found.")
                    continue

                # Now, loop over INSTANCES (the auth_asym_ids)
                for auth_asym_id in auth_asym_ids:
                    
                    assembly_id = assembly_map.get(auth_asym_id, 0)
                    
                    # 1. Create the base Polymer model for this specific instance
                    base_polymer = self._raw_to_polymer_base(poly_entity, auth_asym_id, assembly_id)
                    
                    if family and base_polymer.entity_poly_polymer_type == "Protein":
                        # 2. If tubulin, create a TubulinProtein
                        try:
                            tub_protein = TubulinProtein(
                                **base_polymer.model_dump(),
                                family=family,
                                # --- TUBE-FIX: REMOVED PFAM, ADDED UNIPROT ---
                                uniprot_accession=uniprot_accessions,
                                pfam_accessions=[],   # Set to empty list as per model
                                pfam_comments=[],     # Set to empty list as per model
                                pfam_descriptions=[]  # Set to empty list as per model
                            )
                            tubulin_proteins.append(tub_protein)
                        except ValidationError as e:
                            print(f"    [Debug] ---> ERROR: Validation error for TubulinProtein {auth_asym_id}: {e}")
                    else:
                        # 3. Otherwise, add as a generic Polymer
                        other_polymers.append(base_polymer)

            print(f"    [Debug] Finished processing. Found {len(tubulin_proteins)} tub-proteins, {len(other_polymers)} other.")
            return tubulin_proteins, other_polymers
            
    def _process_nonpolymers(self, nonpolymers_data: dict) -> List[NonpolymericLigand]:
            """Process non-polymer entities into Ligand models."""
            
            ligands: List[NonpolymericLigand] = []
            if not nonpolymers_data: 
                return []

            ligand_entities = nonpolymers_data.get("nonpolymer_entities") or []
            
            if not ligand_entities:
                return []
                
            lig_chem_ids = [
                np["pdbx_entity_nonpoly"]["comp_id"] 
                for np in ligand_entities 
                if "pdbx_entity_nonpoly" in np
            ]
            
            if not lig_chem_ids:
                return []
                
            chem_info_map = {}
            try:
                chem_info_query = LigandsChemInfo.replace(
                    "$COMP_IDS", 
                    json.dumps(list(set(lig_chem_ids)))
                )
                chem_data = query_rcsb_api(chem_info_query)
                
                if chem_data and "chem_comps" in chem_data:
                    for comp in chem_data["chem_comps"]:
                        chem_id = comp["chem_comp"]["id"]
                        chem_info_map[chem_id] = comp.get("rcsb_chem_comp_descriptor", {})
            except Exception as e:
                print(f"Warning: Could not fetch chemical info for ligands. {e}")

            # Build ligand models
            for np_entity in ligand_entities:
                try:
                    chem_id = np_entity["pdbx_entity_nonpoly"]["comp_id"]
                    chem_info = chem_info_map.get(chem_id, {})
                    
                    ligand = NonpolymericLigand(
                        chemicalId          = chem_id,
                        chemicalName        = np_entity["pdbx_entity_nonpoly"]["name"],
                        formula_weight      = np_entity["rcsb_nonpolymer_entity"].get("formula_weight"),
                        pdbx_description    = np_entity["rcsb_nonpolymer_entity"]["pdbx_description"],
                        number_of_instances = np_entity["rcsb_nonpolymer_entity"]["pdbx_number_of_molecules"],
                        nonpolymer_comp     = np_entity.get("nonpolymer_comp"),                                  # This carries drugbank, etc.
                        SMILES              = chem_info.get("SMILES"),
                        SMILES_stereo       = chem_info.get("SMILES_stereo"),
                        InChI               = chem_info.get("InChI"),
                        InChIKey            = chem_info.get("InChIKey"),
                    )
                    ligands.append(ligand)
                except (ValidationError, KeyError) as e:
                    print(f"Skipping ligand: Validation/Key error {e} for entity {np_entity}")
            
            return ligands

    def _infer_organisms(self, polymers: List[Polymer]) -> dict[str, Any]:
        """
        Infer parent organisms from all polymers, same logic as riboxyz.
        """
        all_src_ids: List[int] = []
        all_host_ids: List[int] = []
        all_src_names: List[str] = []
        all_host_names: List[str] = []

        for poly in polymers:
            all_src_ids.extend(poly.src_organism_ids)
            all_host_ids.extend(poly.host_organism_ids)
            all_src_names.extend(poly.src_organism_names)
            all_host_names.extend(poly.host_organism_names)
            
        # This just returns all unique IDs/names found.
        # The "tally" logic can be added here if you prefer a single ID.
        return {
            "src_organism_ids": sorted(list(set(all_src_ids))),
            "src_organism_names": sorted(list(set(all_src_names))),
            "host_organism_ids": sorted(list(set(all_host_ids))),
            "host_organism_names": sorted(list(set(all_host_names))),
        }

    def _get_assembly_mappings(self, auth_asym_ids: List[str]) -> Dict[str, int]:
            """
            Returns a DICTIONARY mapping auth_asym_id -> assembly_id.
            """
            if len(self.asm_maps) <= 1:
                return {auth_id: 0 for auth_id in auth_asym_ids}  # Default to assembly 0

            mapping = {}
            all_ids = set(auth_asym_ids)
            
            # First, build a map of all polymer instances
            for asm in self.asm_maps:
                assembly_idx = int(asm.rcsb_id.split("-")[1]) - 1
                for poly_inst in asm.polymer_entity_instances:
                    auth_id = poly_inst.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id
                    if auth_id in all_ids:
                        mapping[auth_id] = assembly_idx
            
            # (Optional) Map non-polymers too
            for asm in self.asm_maps:
                if asm.nonpolymer_entity_instances:
                    assembly_idx = int(asm.rcsb_id.split("-")[1]) - 1
                    for nonpoly_inst in asm.nonpolymer_entity_instances:
                        auth_id = nonpoly_inst.rcsb_nonpolymer_entity_instance_container_identifiers.auth_asym_id
                        if auth_id in all_ids and auth_id not in mapping:
                            mapping[auth_id] = assembly_idx

            # Assign 0 to any unmapped
            for auth_id in auth_asym_ids:
                if auth_id not in mapping:
                    mapping[auth_id] = 0
                    
            return mapping


    def sequence_ingestion(self, protein: TubulinProtein) -> None:
        """
        Runs the full alignment/mutation detection pipeline for a single tubulin protein.
        Only processes Alpha and Beta tubulins.
        Results are stored per-structure in sequence_ingestion.json, keyed by auth_asym_id.
        """
        # DEBUG OUTPUT
        print(f"DEBUG: TUBETL_DATA = {TUBETL_DATA}")
        print(f"DEBUG: self.assets.paths.base_dir = {self.assets.paths.base_dir}")
        print(f"DEBUG: self.assets.paths.sequence_ingestion = {self.assets.paths.sequence_ingestion}")
        print(f"DEBUG: Does base_dir exist? {os.path.exists(self.assets.paths.base_dir)}")
        
        # Only process alpha and beta
        if protein.family not in [TubulinFamily.ALPHA, TubulinFamily.BETA]:
            print(f"    [Sequence Ingestion] Skipping {protein.auth_asym_id} - not Alpha/Beta (is {protein.family.value})")
            return

        try:
            # Import here to avoid circular dependencies
            from api.services.alignment import TubulinIngestor
            from dataclasses import asdict
            
            MUSCLE_BINARY = str(PROJECT_ROOT / "muscle3.8.1")
            
            # Select correct master profile based on family
            if protein.family == TubulinFamily.ALPHA:
                MASTER_PROFILE = str(PROJECT_ROOT / "data" / "alpha_tubulin" / "alpha_tubulin.afasta")
                profile_type = "Alpha"
            elif protein.family == TubulinFamily.BETA:
                MASTER_PROFILE = str(PROJECT_ROOT / "data" / "beta_tubulin" / "beta_tubulin.afasta")
                profile_type = "Beta"
            else:
                print(f"    [Sequence Ingestion] Unsupported family: {protein.family.value}")
                return
            
            print(f"\n{'='*60}")
            print(f"SEQUENCE INGESTION: {self.rcsb_id} CHAIN {protein.auth_asym_id}")
            print(f"Family: {protein.family.value} (using {profile_type} profile)")
            print(f"{'='*60}")
            
            # Initialize ingestor with the correct profile
            ingestor = TubulinIngestor(MASTER_PROFILE, MUSCLE_BINARY)
            
            # Run ingestion
            result = ingestor.process_chain(
                self.rcsb_id, 
                protein.auth_asym_id, 
                protein.family.value
            )
            
            # ENSURE DIRECTORY EXISTS
            output_file = Path(self.assets.paths.sequence_ingestion)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            print(f"DEBUG: About to save to: {output_file}")
            print(f"DEBUG: Output file parent exists? {output_file.parent.exists()}")
            
            # Load existing results or create new dict
            if output_file.exists():
                with open(output_file, "r") as f:
                    all_results = json.load(f)
            else:
                all_results = {}
            
            # Add/update this chain's results
            timestamp = datetime.now().isoformat()
            all_results[protein.auth_asym_id] = {
                "processed_at": timestamp,
                "family": profile_type,
                "data": asdict(result)
            }
            
            # Save back
            with open(output_file, "w") as f:
                json.dump(all_results, f, indent=2)
            
            print(f"DEBUG: File written. Exists now? {output_file.exists()}")
            
            # Log stats
            print(f"\n📊 INGESTION STATS:")
            print(f"  - Mutations detected: {result.stats['total_mutations']}")
            print(f"  - Insertions: {result.stats['insertions']}")
            print(f"  - MA Coverage: {result.stats['ma_coverage']}/{len(result.ma_to_auth_map)}")
            print(f"  - Results saved to: {output_file}")
            
            if result.mutations:
                print(f"\n🧬 MUTATIONS:")
                for mut in result.mutations[:5]:
                    print(f"  - MA pos {mut.ma_position}: {mut.wild_type} → {mut.observed} (PDB ID: {mut.pdb_auth_id})")
                if len(result.mutations) > 5:
                    print(f"  ... and {len(result.mutations) - 5} more")
            
            print(f"{'='*60}\n")
            
        except Exception as e:
            print(f"⚠️ SEQUENCE INGESTION FAILED for {protein.auth_asym_id}: {e}")
            import traceback
            traceback.print_exc()

    async def generate_profile(self, overwrite: bool = False, run_sequence_ingestion: bool = True) -> TubulinStructure:
        """
        Main ETL method to generate and save the TubulinStructure profile.
        
        Args:
            overwrite: Whether to overwrite existing profile
            run_sequence_ingestion: Whether to run alignment/mutation detection for Alpha/Beta chains
        """
        
        profile_path = self.assets.paths.profile
        if os.path.exists(profile_path) and not overwrite:
            print(f"Profile for {self.rcsb_id} already exists. Loading from disk.")
            return self.assets.profile()
            
        print(f"Generating new profile for {self.rcsb_id}...")
        
        # 1. Fetch Assembly Info (needed first for mapping)
        asm_data = query_rcsb_api(
            AssemblyIdentificationString.replace("$RCSB_ID", self.rcsb_id)
        )
        self.asm_maps = [
            AssemblyInstancesMap.model_validate(asm) 
            for asm in asm_data.get("entry", {}).get("assemblies", [])
        ]

        # 2. Fetch and Process Polymers
        polymers_data = query_rcsb_api(
            PolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id)
        )["entry"]
        proteins, other_polymers = self._process_polymers(polymers_data)

        # 3. Fetch and Process Non-Polymers (Ligands)
        nonpolymers_data = query_rcsb_api(
            NonpolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id)
        )["entry"]
        ligands = self._process_nonpolymers(nonpolymers_data)

        # 4. Fetch Structure-level Metadata
        entry_data = query_rcsb_api(
            EntryInfoString.replace("$RCSB_ID", self.rcsb_id)
        )["entry"]
        
        # Extract citation info safely
        citation = (entry_data.get("citation") or [{}])[0]
        
        # Extract external refs
        external_refs = entry_data.get("rcsb_external_references", []) or []

        # 5. Infer Organisms
        organism_info = self._infer_organisms(proteins + other_polymers)

        # 6. Assemble the Final TubulinStructure Profile
        try:
            structure_profile = TubulinStructure(
                rcsb_id         = entry_data["rcsb_id"],
                expMethod       = entry_data["exptl"][0]["method"],
                resolution      = (entry_data["rcsb_entry_info"]["resolution_combined"] or [None])[0],
                deposition_date = entry_data["rcsb_accession_info"].get("deposit_date"),
                
                pdbx_keywords      = entry_data.get("struct_keywords", {}).get("pdbx_keywords"),
                pdbx_keywords_text = entry_data.get("struct_keywords", {}).get("text"),

                rcsb_external_ref_id   = [ref.get("id") for ref in external_refs],
                rcsb_external_ref_type = [ref.get("type") for ref in external_refs],
                rcsb_external_ref_link = [ref.get("link") for ref in external_refs],

                citation_year         = citation.get("year"),
                citation_rcsb_authors = citation.get("rcsb_authors"),
                citation_title        = citation.get("title"),
                citation_pdbx_doi     = citation.get("pdbx_database_id_DOI"),

                **organism_info,

                proteins             = proteins,
                other_polymers       = other_polymers,
                nonpolymeric_ligands = ligands,
                assembly_map         = self.asm_maps,
            )

            # 7. Save Profile to Disk
            self.assets._verify_dir_exists()
            with open(profile_path, "w") as f:
                f.write(structure_profile.model_dump_json(indent=4))
            
            print(f"Successfully generated and saved profile to: {profile_path}")
            
            # 8. Run Sequence Ingestion for Alpha/Beta chains
            if run_sequence_ingestion and proteins:
                print(f"\n{'='*60}")
                print(f"RUNNING SEQUENCE INGESTION FOR ALPHA/BETA CHAINS")
                print(f"{'='*60}\n")
                
                for protein in proteins:
                    self.sequence_ingestion(protein)
            
            return structure_profile

        except ValidationError as e:
            print(f"CRITICAL: Failed to validate final TubulinStructure for {self.rcsb_id}.")
            print(e)
            raise e


async def main(rcsb_id: str):
    collector = TubulinETLCollector(rcsb_id)
    try:
        profile = await collector.generate_profile(overwrite=True, run_sequence_ingestion=True)
        print(f"\n--- Profile Summary for {rcsb_id} ---")
        print(f"Resolution: {profile.resolution}")
        print(f"Tubulin Proteins: {len(profile.proteins)}")
        
    except Exception as e:
        print(f"An error occurred while processing {rcsb_id}: {e}")

if __name__ == "__main__":
    os.makedirs(TUBETL_DATA, exist_ok=True)
    TEST_RCSB_ID = "1JFF" 
    asyncio.run(main(TEST_RCSB_ID))