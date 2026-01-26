"""
Tubulin structure ETL collector.
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


from lib.etl.classification import (
    classify_sequence,
    build_entity_classification_result,
    is_map_family,
    is_tubulin_family,
)
from lib.types import ClassificationReport, EntityClassificationResult, SequenceVariant
from lib.etl.molstar_bridge import run_molstar_extraction
from lib.etl.sequence_alignment import get_aligner_for_family, AlignmentResult
from lib.etl.augmentation import augment_binding_sites
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
    MolstarExtractionResult,
    ObservedSequenceData,
    PolymerClass,
    VariantsFile,
    LigandBindingSitesFile,
)


EXTRACT_SCRIPT = PROJECT_ROOT / "scripts_and_artifacts" / "extract_structure_data.tsx"


class TubulinETLCollector:
    """Collects and processes tubulin structure data."""

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id.upper()
        self.assets = TubulinStructureAssets(self.rcsb_id)
        self.asm_maps: List[AssemblyInstancesMap] = []

    async def generate_profile(self, overwrite: bool = False) -> TubulinStructure:
        """Generate a complete structure profile."""

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

        cif_path = self._download_cif()

        logger.debug("  Fetching RCSB metadata...")
        (
            entry_data,
            poly_entities,
            ligand_entities,
            polypeptides,
            polynucleotides,
            nonpolymers,
        ) = self._fetch_rcsb_metadata()
        logger.debug(
            f"  RCSB nonpolymer: {len(ligand_entities)} entities, {len(nonpolymers)} instances"
        )
        if len(nonpolymers) == 0 and len(ligand_entities) > 0:
            logger.warning(
                "  Nonpolymer entities exist but instances are empty. "
                "Likely NonpolymerEntitiesString query missing nonpolymer_entity_instances."
            )


        logger.debug("  Running Molstar extraction...")
        molstar_output = Path(self.assets.paths.molstar_raw_extraction)
        molstar_result = run_molstar_extraction(
            cif_path=cif_path,
            rcsb_id=self.rcsb_id,
            output_path=molstar_output,
            script_path=EXTRACT_SCRIPT,
            project_root=PROJECT_ROOT,
        )

        if not molstar_result:
            raise RuntimeError(f"Molstar extraction failed for {self.rcsb_id}")

        logger.debug(
            f"    {len(molstar_result.sequences)} chains, "
            f"{len(molstar_result.ligand_neighborhoods)} ligand neighborhoods"
        )

        chain_to_entity = self._build_chain_to_entity_map(polypeptides)
        observed_by_entity = self._map_observed_to_entities(
            molstar_result.sequences, chain_to_entity
        )

        # ========================================
        # Phase 2: Classification
        # ========================================
        logger.info("Phase 2: Classifying chains")

        entity_families: Dict[str, Optional[PolymerClass]] = {}
        classification_results: Dict[str, EntityClassificationResult] = {}

        tubulin_count = 0
        map_count = 0
        classified_count = 0
        total_count = 0

        for entity_id, entity in poly_entities.items():
            if not isinstance(entity, PolypeptideEntity):
                continue

            total_count += 1

            observed = observed_by_entity.get(entity_id)
            if not observed:
                logger.warning(f"  {entity_id}: No observed sequence")
                classification_results[entity_id] = EntityClassificationResult(
                    entity_id=entity_id,
                    auth_asym_ids=entity.pdbx_strand_ids,
                    sequence_length=entity.sequence_length,
                    assigned_family=None,
                )
                continue
            logger.debug(
                f"  Entity {entity_id}: sequence length = {len(observed.sequence)}, first 50 chars: {observed.sequence[:50]}"
            )

            family, hmm_result = classify_sequence(observed, self.rcsb_id)
            entity_families[entity_id] = family
            entity.family = family

            classification_results[entity_id] = build_entity_classification_result(
                entity_id=entity_id,
                auth_asym_ids=entity.pdbx_strand_ids,
                sequence_length=len(observed.sequence),
                hmm_result=hmm_result,
            )

            if family:
                classified_count += 1
                if is_tubulin_family(family):
                    tubulin_count += 1
                elif is_map_family(family):
                    map_count += 1

            status = family.value if family else "unclassified"
            marker = "+" if family else "-"
            score_str = (
                f" (score={hmm_result.best_hit.score:.1f})"
                if hmm_result.best_hit
                else ""
            )
            logger.debug(f"  Entity {entity_id}: {status}{score_str} {marker}")

        classification_report = ClassificationReport(
            rcsb_id=self.rcsb_id,
            generated_at=datetime.now().isoformat(),
            summary={
                "total": total_count,
                "classified": classified_count,
                "unclassified": total_count - classified_count,
                "tubulin": tubulin_count,
                "map": map_count,
            },
            entities=classification_results,
        )

        with open(self.assets.paths.classification_report, "w") as f:
            f.write(classification_report.model_dump_json(indent=2))

        logger.debug(
            f"  Classification: {classified_count}/{total_count} entities "
            f"({tubulin_count} tubulin, {map_count} MAP)"
        )

        # ========================================
        # Phase 3: Sequence Alignment
        # ========================================
        logger.info("Phase 3: Aligning tubulin sequences")

        chain_alignments: Dict[str, AlignmentResult] = {}
        all_variants: Dict[str, List[SequenceVariant]] = {}

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

                    entity = poly_entities[entity_id]
                    if isinstance(entity, PolypeptideEntity):
                        entity.index_mapping = result.index_mapping
                        entity.variants = result.variants
                        entity.alignment_stats = result.stats

                    all_variants[entity_id] = result.variants

                    logger.debug(
                        f"  Entity {entity_id}: {result.stats['ma_coverage']} residues mapped, "
                        f"{len(result.substitutions)} sub, {len(result.insertions)} ins, {len(result.deletions)} del"
                    )

            except Exception as e:
                logger.warning(f"  Entity {entity_id}: Alignment failed - {e}")
                continue

        if all_variants:
            variants_file = VariantsFile(
                rcsb_id=self.rcsb_id,
                generated_at=datetime.now().isoformat(),
                entities=all_variants,
            )
            with open(self.assets.paths.variants_file, "w") as f:
                f.write(variants_file.model_dump_json(indent=2))
            logger.debug(f"  Wrote variants to {self.assets.paths.variants_file}")

        # ========================================
        # Phase 4: Augmenting extracted data
        # ========================================
        logger.info("Phase 4: Augmenting extracted data")

        augmented_binding_sites = augment_binding_sites(
            binding_sites=molstar_result.ligand_neighborhoods,
            chain_alignments=chain_alignments,
        )

        if augmented_binding_sites:
            sites_file = LigandBindingSitesFile(
                rcsb_id=self.rcsb_id,
                generated_at=datetime.now().isoformat(),
                binding_sites=augmented_binding_sites,
            )
            with open(self.assets.paths.binding_sites_file, "w") as f:
                f.write(sites_file.model_dump_json(indent=2))
            logger.debug(f"  Wrote {len(augmented_binding_sites)} binding sites")

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
            ligand_binding_sites=augmented_binding_sites,
            assembly_map=self.asm_maps,
        )

        with open(self.assets.paths.profile, "w") as f:
            f.write(structure.model_dump_json(indent=2))

        logger.info(f"Profile written to {self.assets.paths.profile}")

        return structure

    # ========================================
    # Helper Methods
    # ========================================

    def _download_cif(self) -> Path:
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
        entry_data = query_rcsb_api(EntryInfoString.replace("$RCSB_ID", self.rcsb_id))[
            "entry"
        ]

        asm_data = query_rcsb_api(
            AssemblyIdentificationString.replace("$RCSB_ID", self.rcsb_id)
        )
        self.asm_maps = [
            AssemblyInstancesMap.model_validate(a)
            for a in asm_data.get("entry", {}).get("assemblies", [])
        ]

        poly_data = query_rcsb_api(
            PolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id)
        )["entry"]
        poly_entities, polypeptides, polynucleotides = self._process_polymers(poly_data)

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
                src_organisms = raw.get("rcsb_entity_source_organism") or []
                host_organisms = raw.get("rcsb_entity_host_organism") or []

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
                        for o in src_organisms
                        if o.get("ncbi_taxonomy_id") is not None
                    ],
                    src_organism_names=[
                        o.get("scientific_name")
                        for o in src_organisms
                        if o.get("scientific_name") is not None
                    ],
                    host_organism_ids=[
                        o.get("ncbi_taxonomy_id")
                        for o in host_organisms
                        if o.get("ncbi_taxonomy_id") is not None
                    ],
                    host_organism_names=[
                        o.get("scientific_name")
                        for o in host_organisms
                        if o.get("scientific_name") is not None
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

    # lib/etl/collector.py

    def _process_nonpolymers(self, nonpoly_data: dict):
        entity_map = {}
        instances = []

        # 1) Blueprints
        raw_entities = nonpoly_data.get("nonpolymer_entities") or []

        # Batch fetch chemical info
        comp_ids = list(set(r["pdbx_entity_nonpoly"]["comp_id"] for r in raw_entities))
        chem_info_map = self._fetch_chem_info(comp_ids)

        for raw_ent in raw_entities:
            ids = raw_ent["pdbx_entity_nonpoly"]
            ent_id = ids["entity_id"]
            comp_id = ids["comp_id"]
            chem = chem_info_map.get(comp_id, {})

            ent = NonpolymerEntity(
                entity_id        = ent_id,
                chemical_id      = comp_id,
                chemical_name    = ids["name"],
                pdbx_description = raw_ent["rcsb_nonpolymer_entity"]["pdbx_description"],
                formula_weight   = raw_ent["rcsb_nonpolymer_entity"].get("formula_weight"),
                SMILES           = chem.get("SMILES"),
                InChIKey         = chem.get("InChIKey"),
                nonpolymer_comp  = raw_ent.get("nonpolymer_comp"),
            )
            entity_map[ent_id] = ent

        # 2) Instances
        # Your schema puts instances UNDER each nonpolymer entity.
        raw_instances = []
        for raw_ent in raw_entities:
            raw_instances.extend(raw_ent.get("nonpolymer_entity_instances") or [])

        # Extract auth_asym_ids for assembly mapping
        all_auth_ids = []
        for inst in raw_instances:
            ident = inst.get("rcsb_nonpolymer_entity_instance_container_identifiers") or {}
            if ident.get("auth_asym_id"):
                all_auth_ids.append(ident["auth_asym_id"])

        assembly_lookup = self._get_assembly_mappings(all_auth_ids)

        for raw_inst in raw_instances:
            ident = raw_inst["rcsb_nonpolymer_entity_instance_container_identifiers"]

            auth_id = ident["auth_asym_id"]
            instances.append(
                Nonpolymer(
                    parent_rcsb_id=self.rcsb_id,
                    auth_asym_id=auth_id,
                    asym_id=ident["asym_id"],
                    auth_seq_id=ident["auth_seq_id"],
                    entity_id=ident["entity_id"],
                    assembly_id=assembly_lookup.get(auth_id, 0),
                )
            )

        return entity_map, instances


    def _fetch_chem_info(self, comp_ids: List[str]) -> dict:
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
        return {p.auth_asym_id: p.entity_id for p in polypeptides}

    def _map_observed_to_entities(
        self,
        sequences: List[ObservedSequenceData],
        chain_to_entity: Dict[str, str],
    ) -> Dict[str, ObservedSequenceData]:
        result = {}
        for seq in sequences:
            entity_id = chain_to_entity.get(seq.auth_asym_id)
            if entity_id and entity_id not in result:
                result[entity_id] = seq
        return result

    def _infer_organisms(self, entities) -> dict:
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
    await TubulinETLCollector(rcsb_id).generate_profile(overwrite=True)


if __name__ == "__main__":
    import sys

    asyncio.run(main(sys.argv[1] if len(sys.argv) > 1 else "6WVR"))
