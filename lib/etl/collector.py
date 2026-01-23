"""
Tubulin structure ETL collector.

Orchestrates the complete pipeline:
  1. Acquire raw data (RCSB API + Molstar extraction)
  2. Classify polymer chains
  3. Align tubulin sequences to family MSAs
  4. Augment extracted data with MA indices
  5. Assemble and persist the profile
"""

import asyncio
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

from loguru import logger

from api.config import PROJECT_ROOT
from lib.etl.assets import (
    AssemblyIdentificationString,
    EntryInfoString,
    LigandsChemInfo,
    NonpolymerEntitiesString,
    PolymerEntitiesString,
    TubulinStructureAssets,
    query_rcsb_api,
)
from lib.etl.molstar_bridge import run_molstar_extraction
from lib.etl.classification import classify_sequence, is_tubulin_family
from lib.etl.sequence_alignment import get_aligner_for_family, AlignmentResult
from lib.etl.augmentation import (
    augment_ligand_neighborhoods,
    extract_ptms,
    extract_map_interfaces,
)
from lib.types import (
    TubulinStructure,
    PolypeptideEntity,
    PolynucleotideEntity,
    NonpolymerEntity,
    Polypeptide,
    Polynucleotide,
    Nonpolymer,
    TubulinFamily,
    AssemblyInstancesMap,
    Mutation,
    MolstarExtractionResult,
    ObservedSequenceData,
    ProcessedChainData,
    SequenceIngestionEntry,
    MutationEntryData,
    PolymerClass,
)


EXTRACT_SCRIPT = PROJECT_ROOT / "scripts_and_artifacts" / "extract_structure_data.tsx"


class TubulinETLCollector:
    """
    Collects and processes tubulin structure data.
    """

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id.upper()
        self.assets = TubulinStructureAssets(self.rcsb_id)
        self.asm_maps: List[AssemblyInstancesMap] = []

    async def generate_profile(
        self,
        overwrite: bool = False,
    ) -> TubulinStructure:
        """
        Generate a complete structure profile.

        Phases:
          1. Acquire raw data
          2. Classify chains
          3. Align tubulin sequences
          4. Augment with MA indices
          5. Assemble and persist
        """
        # Check for existing profile
        if (
            self.assets.paths.profile
            and Path(self.assets.paths.profile).exists()
            and not overwrite
        ):
            return self.assets.profile()

        logger.info(f"{'─' * 60}")
        logger.info(f"Processing {self.rcsb_id}")
        logger.info(f"{'─' * 60}")

        self.assets._verify_dir_exists()

        # ========================================
        # Phase 1: Acquire Raw Data
        # ========================================
        logger.info("Phase 1: Acquiring raw data")

        # 1a. Download structure file
        cif_path = self._download_cif()

        # 1b. Fetch RCSB metadata
        logger.debug("  Fetching RCSB metadata...")
        (
            entry_data,
            poly_entities,
            ligand_entities,
            polypeptides,
            polynucleotides,
            nonpolymers,
        ) = self._fetch_rcsb_metadata()

        # 1c. Run Molstar extraction
        logger.debug("  Running Molstar extraction...")
        extraction_output = Path(self.assets.paths.base_dir) / "molstar_extraction.json"
        molstar_result = run_molstar_extraction(
            cif_path=cif_path,
            rcsb_id=self.rcsb_id,
            output_path=extraction_output,
            script_path=EXTRACT_SCRIPT,
            project_root=PROJECT_ROOT,
        )

        if not molstar_result:
            raise RuntimeError(f"Molstar extraction failed for {self.rcsb_id}")

        logger.debug(
            f"    Extracted {len(molstar_result.sequences)} chains, "
            f"{len(molstar_result.ligand_neighborhoods)} ligand neighborhoods"
        )

        # Build lookup maps
        chain_to_entity = self._build_chain_to_entity_map(polypeptides)
        observed_by_entity = self._map_observed_to_entities(
            molstar_result.sequences, chain_to_entity
        )

        # ========================================
        # Phase 2: Classification
        # ========================================
        logger.info("Phase 2: Classifying chains")

        entity_families: Dict[str, Optional[PolymerClass]] = {}

        for entity_id, entity in poly_entities.items():
            if not isinstance(entity, PolypeptideEntity):
                continue

            observed = observed_by_entity.get(entity_id)
            if not observed:
                logger.warning(f"  {entity_id}: No observed sequence found")
                continue

            family = classify_sequence(observed, self.rcsb_id)
            entity_families[entity_id] = family
            entity.family = family

            status = family.value if family else "unclassified"
            marker = "+" if family else "-"
            logger.debug(f"  Entity {entity_id}: {status} {marker}")

        # ========================================
        # Phase 3: Sequence Alignment
        # ========================================
        logger.info("Phase 3: Aligning tubulin sequences")

        chain_alignments: Dict[str, AlignmentResult] = {}
        ingestion_data: Dict[str, SequenceIngestionEntry] = {}

        for entity_id, family in entity_families.items():
            if not is_tubulin_family(family):
                continue

            observed = observed_by_entity.get(entity_id)
            if not observed:
                continue

            try:
                aligner = get_aligner_for_family(family, PROJECT_ROOT)
                identifier = f"{self.rcsb_id}_{entity_id}"

                result = aligner.align(observed, family, identifier)

                if result:
                    chain_alignments[observed.auth_asym_id] = result

                    # Store for entity
                    entity = poly_entities[entity_id]
                    if isinstance(entity, PolypeptideEntity):
                        entity.mutations = [
                            Mutation(
                                master_index=m.ma_position,
                                from_residue=m.wild_type,
                                to_residue=m.observed,
                                uniprot_id=entity.uniprot_accessions[0]
                                if entity.uniprot_accessions
                                else "UNK",
                                species=entity.src_organism_names[0]
                                if entity.src_organism_names
                                else "UNK",
                                tubulin_type=family,
                                phenotype="Canonical",
                                database_source="PDB",
                                reference_link=f"https://rcsb.org/structure/{self.rcsb_id}",
                                keywords="canonical_mutation",
                            )
                            for m in result.mutations
                        ]
                        entity.alignment_stats = result.stats

                    # Build ingestion entry
                    ingestion_data[entity_id] = SequenceIngestionEntry(
                        processed_at=datetime.now().isoformat(),
                        family=family.value,
                        data=ProcessedChainData(
                            pdb_id=self.rcsb_id,
                            chain_id=observed.auth_asym_id,
                            tubulin_class=family.value,
                            sequence=result.sequence,
                            ma_to_auth_map=result.ma_to_auth_map,
                            auth_to_ma_json=json.dumps(result.auth_to_ma),
                            mutations=[
                                MutationEntryData(
                                    ma_position=m.ma_position,
                                    wild_type=m.wild_type,
                                    observed=m.observed,
                                    pdb_auth_id=m.pdb_auth_id,
                                )
                                for m in result.mutations
                            ],
                            stats=result.stats,
                        ),
                    )

                    coverage = result.stats.get("ma_coverage", 0)
                    logger.debug(
                        f"  Entity {entity_id} ({observed.auth_asym_id}): "
                        f"aligned, {coverage} residues mapped"
                    )

            except Exception as e:
                logger.warning(f"  Entity {entity_id}: Alignment failed - {e}")
                continue

        # Persist ingestion data
        if ingestion_data:
            ingestion_path = Path(self.assets.paths.sequence_ingestion)
            with open(ingestion_path, "w") as f:
                json.dump(
                    {k: v.model_dump() for k, v in ingestion_data.items()},
                    f,
                    indent=2,
                )

        # ========================================
        # Phase 4: Augmentation
        # ========================================
        logger.info("Phase 4: Augmenting extracted data")

        # Augment ligand neighborhoods
        augmented_neighborhoods = augment_ligand_neighborhoods(
            neighborhoods=molstar_result.ligand_neighborhoods,
            chain_alignments=chain_alignments,
            chain_to_entity=chain_to_entity,
            entity_families=entity_families,
        )

        # Write individual ligand files (for compatibility)
        for neighborhood in augmented_neighborhoods:
            output_path = Path(
                self.assets.paths.ligand_neighborhood(
                    neighborhood.ligand_comp_id,
                    neighborhood.ligand_auth_asym_id,
                )
            )
            with open(output_path, "w") as f:
                json.dump(
                    {
                        "ligand": [
                            neighborhood.ligand_auth_asym_id,
                            neighborhood.ligand_auth_seq_id,
                            neighborhood.ligand_comp_id,
                        ],
                        "interactions": [
                            {
                                "type": ix.type,
                                "participants": [p.to_tuple() for p in ix.participants],
                            }
                            for ix in neighborhood.interactions
                        ],
                        "neighborhood": [
                            r.to_tuple() for r in neighborhood.neighborhood
                        ],
                    },
                    f,
                    indent=2,
                )

        logger.debug(
            f"  Wrote {len(augmented_neighborhoods)} ligand neighborhood files"
        )

        # Future: PTMs and MAP interfaces (stubs)
        _ = extract_ptms(chain_alignments)
        _ = extract_map_interfaces(
            chain_alignments,
            {k: entity_families.get(chain_to_entity.get(k)) for k in chain_alignments},
        )

        # ========================================
        # Phase 5: Assemble and Persist
        # ========================================
        logger.info("Phase 5: Assembling profile")

        all_entities = {**poly_entities, **ligand_entities}
        org_info = self._infer_organisms(list(all_entities.values()))
        citation = (entry_data.get("citation") or [{}])[0]

        structure = TubulinStructure(
            rcsb_id=self.rcsb_id,
            expMethod=entry_data["exptl"][0]["method"],
            resolution=(entry_data["rcsb_entry_info"]["resolution_combined"] or [0.0])[
                0
            ],
            deposition_date=entry_data["rcsb_accession_info"].get("deposit_date"),
            pdbx_keywords=(entry_data.get("struct_keywords") or {}).get(
                "pdbx_keywords"
            ),
            pdbx_keywords_text=(entry_data.get("struct_keywords") or {}).get("text"),
            rcsb_external_ref_id=[
                r.get("id") for r in (entry_data.get("rcsb_external_references") or [])
            ],
            rcsb_external_ref_type=[
                r.get("type")
                for r in (entry_data.get("rcsb_external_references") or [])
            ],
            rcsb_external_ref_link=[
                r.get("link")
                for r in (entry_data.get("rcsb_external_references") or [])
            ],
            citation_title=citation.get("title"),
            citation_year=citation.get("year"),
            citation_rcsb_authors=citation.get("rcsb_authors"),
            citation_pdbx_doi=citation.get("pdbx_database_id_DOI"),
            **org_info,
            entities=all_entities,
            polypeptides=polypeptides,
            polynucleotides=polynucleotides,
            nonpolymers=nonpolymers,
            assembly_map=self.asm_maps,
        )

        with open(self.assets.paths.profile, "w") as f:
            f.write(structure.model_dump_json(indent=4))

        logger.info(f"Profile written to {self.assets.paths.profile}")

        return structure

    # ========================================
    # Helper Methods
    # ========================================

    def _download_cif(self) -> Path:
        """Download structure file if not present."""
        import requests

        cif_path = Path(self.assets.paths.cif)
        if cif_path.exists():
            return cif_path

        cif_path.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://files.rcsb.org/download/{self.rcsb_id}.cif"
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        cif_path.write_text(resp.text)

        return cif_path

    def _fetch_rcsb_metadata(self):
        """Fetch all RCSB metadata via GraphQL."""
        # Entry info
        entry_data = query_rcsb_api(EntryInfoString.replace("$RCSB_ID", self.rcsb_id))[
            "entry"
        ]

        # Assembly maps
        asm_data = query_rcsb_api(
            AssemblyIdentificationString.replace("$RCSB_ID", self.rcsb_id)
        )
        self.asm_maps = [
            AssemblyInstancesMap.model_validate(a)
            for a in asm_data.get("entry", {}).get("assemblies", [])
        ]

        # Polymers
        poly_data = query_rcsb_api(
            PolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id)
        )["entry"]
        poly_entities, polypeptides, polynucleotides = self._process_polymers(poly_data)

        # Non-polymers
        nonpoly_data = query_rcsb_api(
            NonpolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id)
        )["entry"]
        ligand_entities, nonpolymers = self._process_nonpolymers(nonpoly_data)

        return (
            entry_data,
            poly_entities,
            ligand_entities,
            polypeptides,
            polynucleotides,
            nonpolymers,
        )

    def _process_polymers(self, polymers_data: dict):
        """Process polymer entities from RCSB data."""
        entity_map = {}
        polypeptides = []
        polynucleotides = []

        raw_entities = polymers_data.get("polymer_entities") or []
        all_auth_ids = [
            aid
            for raw in raw_entities
            for aid in raw.get("rcsb_polymer_entity_container_identifiers", {}).get(
                "auth_asym_ids", []
            )
        ]
        assembly_lookup = self._get_assembly_mappings(all_auth_ids)

        for raw in raw_entities:
            ids = raw["rcsb_polymer_entity_container_identifiers"]
            entity_id = ids["entity_id"]
            poly_type = raw["entity_poly"]["rcsb_entity_polymer_type"]
            strand_ids = (
                [
                    s.strip()
                    for s in raw["entity_poly"].get("pdbx_strand_id", "").split(",")
                ]
                if raw["entity_poly"].get("pdbx_strand_id")
                else []
            )

            if poly_type == "Protein":
                ent = PolypeptideEntity(
                    entity_id=entity_id,
                    pdbx_description=raw.get("rcsb_polymer_entity", {}).get(
                        "pdbx_description"
                    ),
                    pdbx_strand_ids=strand_ids,
                    one_letter_code=raw["entity_poly"]["pdbx_seq_one_letter_code"],
                    one_letter_code_can=raw["entity_poly"][
                        "pdbx_seq_one_letter_code_can"
                    ],
                    sequence_length=raw["entity_poly"]["rcsb_sample_sequence_length"],
                    src_organism_ids=[
                        o.get("ncbi_taxonomy_id")
                        for o in raw.get("rcsb_entity_source_organism", []) or []
                    ],
                    src_organism_names=[
                        o.get("scientific_name")
                        for o in raw.get("rcsb_entity_source_organism", []) or []
                    ],
                    uniprot_accessions=[
                        u.get("rcsb_id")
                        for u in (raw.get("uniprots") or [])
                        if u.get("rcsb_id")
                    ],
                )
            else:
                ent = PolynucleotideEntity(
                    entity_id=entity_id,
                    polymer_type=poly_type,
                    one_letter_code=raw["entity_poly"]["pdbx_seq_one_letter_code"],
                    one_letter_code_can=raw["entity_poly"][
                        "pdbx_seq_one_letter_code_can"
                    ],
                    sequence_length=raw["entity_poly"]["rcsb_sample_sequence_length"],
                )

            entity_map[entity_id] = ent

            for auth_id, asym_id in zip(
                ids.get("auth_asym_ids", []), ids.get("asym_ids", [])
            ):
                asm_id = assembly_lookup.get(auth_id, 0)
                if poly_type == "Protein":
                    polypeptides.append(
                        Polypeptide(
                            parent_rcsb_id=self.rcsb_id,
                            auth_asym_id=auth_id,
                            asym_id=asym_id,
                            entity_id=entity_id,
                            assembly_id=asm_id,
                        )
                    )
                else:
                    polynucleotides.append(
                        Polynucleotide(
                            parent_rcsb_id=self.rcsb_id,
                            auth_asym_id=auth_id,
                            asym_id=asym_id,
                            entity_id=entity_id,
                            assembly_id=asm_id,
                        )
                    )

        return entity_map, polypeptides, polynucleotides

    def _process_nonpolymers(self, nonpoly_data: dict):
        """Process non-polymer entities from RCSB data."""
        entity_map = {}
        instances = []

        raw_list = nonpoly_data.get("nonpolymer_entities") or []
        chem_info_map = self._fetch_chem_info(
            list(set(r["pdbx_entity_nonpoly"]["comp_id"] for r in raw_list))
        )
        all_auth_ids = [
            aid
            for raw in raw_list
            for aid in raw.get("rcsb_nonpolymer_entity_container_identifiers", {}).get(
                "auth_asym_ids", []
            )
        ]
        assembly_lookup = self._get_assembly_mappings(all_auth_ids)

        for raw in raw_list:
            ids = raw["rcsb_nonpolymer_entity_container_identifiers"]
            comp_id = raw["pdbx_entity_nonpoly"]["comp_id"]
            chem = chem_info_map.get(comp_id, {})

            ent = NonpolymerEntity(
                entity_id=ids["entity_id"],
                chemical_id=comp_id,
                chemical_name=raw["pdbx_entity_nonpoly"]["name"],
                pdbx_description=raw["rcsb_nonpolymer_entity"]["pdbx_description"],
                SMILES=chem.get("SMILES"),
                InChIKey=chem.get("InChIKey"),
                nonpolymer_comp=raw.get("nonpolymer_comp"),
            )
            entity_map[ids["entity_id"]] = ent

            for auth_id, asym_id in zip(
                ids.get("auth_asym_ids", []), ids.get("asym_ids", [])
            ):
                instances.append(
                    Nonpolymer(
                        parent_rcsb_id=self.rcsb_id,
                        auth_asym_id=auth_id,
                        asym_id=asym_id,
                        entity_id=ids["entity_id"],
                        assembly_id=assembly_lookup.get(auth_id, 0),
                    )
                )

        return entity_map, instances

    def _fetch_chem_info(self, comp_ids: List[str]) -> dict:
        """Fetch chemical descriptors for ligands."""
        if not comp_ids:
            return {}
        try:
            data = query_rcsb_api(
                LigandsChemInfo.replace("$COMP_IDS", json.dumps(comp_ids))
            )
            return {
                c["chem_comp"]["id"]: c.get("rcsb_chem_comp_descriptor", {})
                for c in data.get("chem_comps", [])
            }
        except Exception:
            return {}

    def _get_assembly_mappings(self, auth_asym_ids: List[str]) -> Dict[str, int]:
        """Map chain IDs to assembly indices."""
        mapping = {}
        target_ids = set(auth_asym_ids)

        for asm in self.asm_maps:
            try:
                aidx = int(asm.rcsb_id.split("-")[1])
            except Exception:
                aidx = 0

            for wrapper in asm.polymer_entity_instances or []:
                ident = wrapper.get(
                    "rcsb_polymer_entity_instance_container_identifiers"
                )
                if ident and ident.auth_asym_id in target_ids:
                    mapping[ident.auth_asym_id] = aidx

            for wrapper in asm.nonpolymer_entity_instances or []:
                ident = wrapper.get(
                    "rcsb_nonpolymer_entity_instance_container_identifiers"
                )
                if ident and ident.auth_asym_id in target_ids:
                    mapping[ident.auth_asym_id] = aidx

        return mapping

    def _build_chain_to_entity_map(
        self, polypeptides: List[Polypeptide]
    ) -> Dict[str, str]:
        """Build auth_asym_id -> entity_id map."""
        return {p.auth_asym_id: p.entity_id for p in polypeptides}

    def _map_observed_to_entities(
        self,
        sequences: List[ObservedSequenceData],
        chain_to_entity: Dict[str, str],
    ) -> Dict[str, ObservedSequenceData]:
        """Map observed sequences to entity IDs."""
        result = {}
        for seq in sequences:
            entity_id = chain_to_entity.get(seq.auth_asym_id)
            if entity_id and entity_id not in result:
                result[entity_id] = seq
        return result

    def _infer_organisms(self, entities) -> dict:
        """Aggregate organism info from entities."""
        all_src_ids, all_host_ids, all_src_names, all_host_names = [], [], [], []

        for ent in entities:
            if hasattr(ent, "src_organism_ids") and ent.src_organism_ids:
                all_src_ids.extend(ent.src_organism_ids)
            if hasattr(ent, "host_organism_ids") and ent.host_organism_ids:
                all_host_ids.extend(ent.host_organism_ids)
            if hasattr(ent, "src_organism_names") and ent.src_organism_names:
                all_src_names.extend(ent.src_organism_names)
            if hasattr(ent, "host_organism_names") and ent.host_organism_names:
                all_host_names.extend(ent.host_organism_names)

        return {
            "src_organism_ids": sorted(list(set(all_src_ids))),
            "src_organism_names": sorted(list(set(all_src_names))),
            "host_organism_ids": sorted(list(set(all_host_ids))),
            "host_organism_names": sorted(list(set(all_host_names))),
        }


async def main(rcsb_id: str):
    """CLI entry point."""
    await TubulinETLCollector(rcsb_id).generate_profile(overwrite=True)


if __name__ == "__main__":
    import sys

    asyncio.run(main(sys.argv[1] if len(sys.argv) > 1 else "6WVR"))
