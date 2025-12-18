import asyncio
import json
import os
import requests
from typing import Dict, List, Optional, Tuple, Union, Any
from pathlib import Path
from datetime import datetime
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

from lib.etl.ligand_extraction import extract_ligands_parallel, SKIP_LIGANDS
from lib.hmm.classifier import TubulinClassifier
from lib.types import (
    PolymerClass,
    MapFamily,
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
)

from lib.etl.constants import TUBETL_DATA

EXTRACT_SCRIPT = PROJECT_ROOT / "scripts_and_artifacts" / "extract_ixs.tsx"


class TubulinETLCollector:
    rcsb_id: str
    assets: TubulinStructureAssets
    asm_maps: List[AssemblyInstancesMap] = []
    _classifier: Optional[TubulinClassifier] = None

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id.upper()
        self.assets = TubulinStructureAssets(self.rcsb_id)
        self._classifier = None

    @property
    def classifier(self) -> TubulinClassifier:
        if self._classifier is None:
            self._classifier = TubulinClassifier(use_cache=True)
        return self._classifier

    def _parse_family_enum(self, family_value: Optional[str]) -> Optional[PolymerClass]:
        if not family_value:
            return None
        for fam in TubulinFamily:
            if fam.value == family_value:
                return fam
        for fam in MapFamily:
            if fam.value == family_value:
                return fam
        return None

    def run_classification(
        self,
        entities: Dict[
            str, Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]
        ],
    ) -> None:
        MIN_SCORE_THRESHOLD = 100.0
        MIN_DELTA_THRESHOLD = 50.0
        report = {
            "rcsb_id": self.rcsb_id,
            "summary": {"total": 0, "classified": 0},
            "entities": {},
        }

        for entity_id, entity in entities.items():
            if not isinstance(entity, PolypeptideEntity):
                continue
            report["summary"]["total"] += 1
            result = self.classifier.classify(
                entity.one_letter_code_can, self.rcsb_id, entity_id
            )
            entity.family = self._parse_family_enum(
                result.assigned_family.value if result.assigned_family else None
            )

            status = (
                result.assigned_family.value
                if result.assigned_family
                else "unclassified"
            )
            marker = "✓" if entity.family else "⚠"
            print(f"      {entity_id}: {status} {marker}")
            if entity.family:
                report["summary"]["classified"] += 1

        with open(self.assets.paths.classification_report, "w") as f:
            json.dump(report, f, indent=2)

    def download_cif(self) -> Path:
        cif_path = Path(self.assets.paths.cif)
        if cif_path.exists():
            return cif_path
        cif_path.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://files.rcsb.org/download/{self.rcsb_id}.cif"
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        cif_path.write_text(resp.text)
        return cif_path

    def _infer_organisms(self, entities: List[Any]) -> dict:
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

    def _process_polymers(self, polymers_data: dict):
        entity_map, polypeptides, polynucleotides = {}, [], []
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
        entity_map, instances = {}, []
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
        except:
            return {}

    def _get_assembly_mappings(self, auth_asym_ids: List[str]) -> Dict[str, int]:
        mapping = {}
        target_ids = set(auth_asym_ids)
        for asm in self.asm_maps:
            try:
                aidx = int(asm.rcsb_id.split("-")[1])
            except:
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

    def _sequence_ingestion_quiet(self, entity: PolypeptideEntity) -> None:
        if not entity.family or entity.family not in [
            TubulinFamily.ALPHA,
            TubulinFamily.BETA,
        ]:
            return
        try:
            from lib.seq_aligner import TubulinIngestor
            from dataclasses import asdict

            MBIN = str(PROJECT_ROOT / "muscle3.8.1")
            MDIR = (
                "alpha_tubulin"
                if entity.family == TubulinFamily.ALPHA
                else "beta_tubulin"
            )
            MPROF = str(PROJECT_ROOT / "data" / MDIR / f"{MDIR}.afasta")
            ingestor = TubulinIngestor(MPROF, MBIN)
            result = ingestor.process_sequence(
                entity.one_letter_code_can,
                f"{self.rcsb_id}_{entity.entity_id}",
                entity.family.value,
            )

            # Safely handle mutations (could be None or empty)
            mutations_list = result.mutations if result.mutations else []
            entity.mutations = [
                Mutation(
                    master_index=m.ma_position,
                    from_residue=m.wild_type,
                    to_residue=m.observed,
                    uniprot_id=(entity.uniprot_accessions[0] if entity.uniprot_accessions else "UNK"),
                    species=(entity.src_organism_names[0] if entity.src_organism_names else "UNK"),
                    tubulin_type=entity.family,
                    phenotype="Canonical",
                    database_source="PDB",
                    reference_link=f"https://rcsb.org/structure/{self.rcsb_id}",
                    keywords="canonical_mutation",
                )
                for m in mutations_list
            ]
            entity.alignment_stats = result.stats or {}

            ifile = Path(self.assets.paths.sequence_ingestion)
            ifile.parent.mkdir(parents=True, exist_ok=True)
            data = json.load(open(ifile)) if ifile.exists() else {}
            data[entity.entity_id] = {
                "processed_at": datetime.now().isoformat(),
                "family": entity.family.value,
                "data": asdict(result),
            }
            with open(ifile, "w") as f:
                json.dump(data, f, indent=2)
            
            coverage = result.stats.get('ma_coverage', 0) if result.stats else 0
            print(f"      {entity.entity_id}: mapped ({coverage} residues aligned)")
        except Exception as e:
            print(f"      Sequence Ingestion failed for {entity.entity_id}: {e}")

    async def generate_profile(
        self,
        overwrite: bool = False,
        run_sequence_ingestion: bool = True,
        run_ligand_extraction: bool = True,
        run_classification: bool = True,
        ligand_extraction_workers: int = 8,
    ) -> TubulinStructure:
        if os.path.exists(self.assets.paths.profile) and not overwrite:
            return self.assets.profile()
        print(f"\n{'─' * 60}\nCollecting {self.rcsb_id}\n{'─' * 60}")
        self.assets._verify_dir_exists()
        self.download_cif()

        print(f"  ◆ Fetching metadata...")
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
        all_entities = {**poly_entities, **ligand_entities}

        if run_classification:
            print(f"  ◆ HMM Classification:")
            self.run_classification(all_entities)

        if run_sequence_ingestion:
            print(f"  ◆ Sequence Ingestion:")
            for eid, ent in all_entities.items():
                if isinstance(ent, PolypeptideEntity):
                    self._sequence_ingestion_quiet(ent)

        org_info = self._infer_organisms(list(all_entities.values()))
        citation = (entry_data.get("citation") or [{}])[0]
        structure = TubulinStructure(
            rcsb_id=self.rcsb_id,
            expMethod=entry_data["exptl"][0]["method"],
            resolution=(entry_data["rcsb_entry_info"]["resolution_combined"] or [0.0])[0],
            deposition_date=entry_data["rcsb_accession_info"].get("deposit_date"),
            pdbx_keywords=(entry_data.get("struct_keywords") or {}).get("pdbx_keywords"),
            pdbx_keywords_text=(entry_data.get("struct_keywords") or {}).get("text"),
            rcsb_external_ref_id=[
                r.get("id") for r in (entry_data.get("rcsb_external_references") or [])
            ],
            rcsb_external_ref_type=[
                r.get("type") for r in (entry_data.get("rcsb_external_references") or [])
            ],
            rcsb_external_ref_link=[
                r.get("link") for r in (entry_data.get("rcsb_external_references") or [])
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

        # In collector.py, replace the ligand extraction block:

        if run_ligand_extraction:
            print(f"  ◆ Ligand Extraction & Augmentation:")
            asym_to_entity = {
                p.auth_asym_id: p.entity_id
                for p in polypeptides
                if isinstance(all_entities.get(p.entity_id), PolypeptideEntity)
            }
            
            # Filter ligand instances safely
            ligand_instances = []
            for n in (nonpolymers or []):
                entity = all_entities.get(n.entity_id)
                if entity and hasattr(entity, 'chemical_id') and entity.chemical_id not in SKIP_LIGANDS:
                    ligand_instances.append((entity.chemical_id, n.auth_asym_id))

            if ligand_instances:
                ingestion_path = Path(self.assets.paths.sequence_ingestion)
                results = extract_ligands_parallel(
                    self.rcsb_id,
                    Path(self.assets.paths.cif),
                    ligand_instances,
                    Path(self.assets.paths.base_dir),
                    EXTRACT_SCRIPT,
                    PROJECT_ROOT,
                    ligand_extraction_workers,
                    overwrite,
                    None,
                    asym_to_entity,
                    ingestion_path,  # Pass Path, not dict
                )
                print(f"    Summary: {sum(1 for r in results if r.success)}/{len(ligand_instances)} ligands processed.")


async def main(rcsb_id: str):
    await TubulinETLCollector(rcsb_id).generate_profile(overwrite=True)


if __name__ == "__main__":
    import sys

    asyncio.run(main(sys.argv[1] if len(sys.argv) > 1 else "6WVR"))
