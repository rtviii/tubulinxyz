import asyncio
import json
import os
import subprocess
import requests
from typing import Dict, List, Optional, Tuple, Union
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

from lib.etl.ligand_extraction import extract_ligands_parallel
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
SKIP_LIGANDS = {
    "HOH",
    "DOD",
    "NA",
    "CL",
    "K",
    "MG",
    "CA",
    "ZN",
    "SO4",
    "PO4",
    "GOL",
    "EDO",
    "ACT",
    "MES",
    "TRS",
    "CIT",
}

EXTRACT_SCRIPT = PROJECT_ROOT / "scripts_and_artifacts" / "extract_ixs.tsx"

class TubulinETLCollector:

    rcsb_id : str
    assets  : TubulinStructureAssets
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
        """Convert string value back to the appropriate enum."""
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
        entities: Dict[str, Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]],
    ) -> None:
        """
        Classify all polypeptide entities via HMM.
        Updates entity.family in place and saves full report to disk.
        """
        MIN_SCORE_THRESHOLD = 100.0
        MIN_DELTA_THRESHOLD = 50.0
        
        report = {
            "rcsb_id": self.rcsb_id,
            "summary": {
                "total_polypeptides": 0,
                "classified": 0,
                "unclassified": 0,
                "warnings": [],
            },
            "entities": {},
        }

        results_log = []  # Collect results for clean output

        for entity_id, entity in entities.items():
            if not isinstance(entity, PolypeptideEntity):
                continue

            report["summary"]["total_polypeptides"] += 1

            result = self.classifier.classify(
                sequence=entity.one_letter_code_can,
                rcsb_id=self.rcsb_id,
                auth_asym_id=entity_id,
            )

            sorted_hits = sorted(result.hits, key=lambda h: h.score, reverse=True)
            
            best_hit = sorted_hits[0] if sorted_hits else None
            second_hit = sorted_hits[1] if len(sorted_hits) > 1 else None
            
            best_score = best_hit.score if best_hit else None
            delta_to_second = (best_hit.score - second_hit.score) if (best_hit and second_hit) else None
            
            confidence = None
            entity_warnings = []
            
            if not best_hit:
                confidence = "none"
                entity_warnings.append("no HMM hits")
            elif best_score < MIN_SCORE_THRESHOLD:
                confidence = "low"
                entity_warnings.append(f"weak hit (score={best_score:.1f})")
            elif delta_to_second is not None and delta_to_second < MIN_DELTA_THRESHOLD:
                confidence = "ambiguous"
                entity_warnings.append(
                    f"ambiguous: {best_hit.family.value} vs {second_hit.family.value} "
                    f"(delta={delta_to_second:.1f})"
                )
            else:
                confidence = "high"

            entity.family = self._parse_family_enum(
                result.assigned_family.value if result.assigned_family else None
            )

            if entity.family:
                report["summary"]["classified"] += 1
            else:
                report["summary"]["unclassified"] += 1

            for w in entity_warnings:
                report["summary"]["warnings"].append(f"Entity {entity_id}: {w}")

            report["entities"][entity_id] = {
                "sequence_length": result.sequence_length,
                "assigned_family": result.assigned_family.value if result.assigned_family else None,
                "confidence": confidence,
                "best_score": best_score,
                "delta_to_second": delta_to_second,
                "warnings": entity_warnings if entity_warnings else None,
                "hits": [
                    {
                        "family": h.family.value,
                        "score": h.score,
                        "evalue": h.evalue,
                        "bias": h.bias,
                        "domains": [
                            {
                                "score": d.score,
                                "c_evalue": d.c_evalue,
                                "i_evalue": d.i_evalue,
                                "env_from": d.env_from,
                                "env_to": d.env_to,
                                "bias": d.bias,
                            }
                            for d in h.domains
                        ],
                    }
                    for h in sorted_hits
                ],
            }

            # Collect for batch output
            assigned = result.assigned_family
            status = assigned.value if assigned else "unclassified"
            warning_str = f" ({', '.join(entity_warnings)})" if entity_warnings else ""
            marker = "✓" if confidence == "high" else "⚠"
            results_log.append(f"     {entity_id}: {status}{warning_str} {marker}")

        # Print results in one block
        for line in results_log:
            print(line)
        
        summary = report["summary"]
        print(f"     [{summary['classified']}/{summary['total_polypeptides']} classified]")

        # Save report
        report_path = self.assets.paths.classification_report
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
    # -------------------------------------------------------------------------
    # CIF Download
    # -------------------------------------------------------------------------

    def download_cif(self) -> Path:
        """Download CIF file from RCSB PDB if not present."""
        cif_path = Path(self.assets.paths.cif)

        if cif_path.exists():
            print(f"CIF exists: {cif_path}")
            return cif_path

        cif_path.parent.mkdir(parents=True, exist_ok=True)

        url = f"https://files.rcsb.org/download/{self.rcsb_id}.cif"
        print(f"Downloading CIF from {url}...")

        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        cif_path.write_text(resp.text)

        print(f"Saved: {cif_path}")
        return cif_path

    # -------------------------------------------------------------------------
    # Ligand Extraction
    # -------------------------------------------------------------------------

    def _run_extraction_script(
        self,
        cif_path: Path,
        comp_id: str,
        auth_asym_id: str,
        output_path: Path,
    ) -> bool:
        """Run the TypeScript extraction script for a single ligand instance."""
        cmd = [
            "npx",
            "tsx",
            str(EXTRACT_SCRIPT),
            str(cif_path),
            comp_id,
            auth_asym_id,
            str(output_path),
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=PROJECT_ROOT,
            )

            if result.returncode != 0:
                print(f"  Extraction failed for {comp_id}/{auth_asym_id}:")
                print(f"  {result.stderr[:500]}")
                return False
            return True

        except subprocess.TimeoutExpired:
            print(f"  Extraction timed out for {comp_id}/{auth_asym_id}")
            return False
        except Exception as e:
            print(f"  Extraction error: {e}")
            return False

    def extract_ligand_interactions(
        self,
        structure: TubulinStructure,
        overwrite: bool = False,
    ) -> List[Path]:
        """
        Extract interactions for all ligand instances in the structure.

        Returns list of paths to generated report files.
        """
        cif_path = Path(self.assets.paths.cif)
        if not cif_path.exists():
            print(f"CIF not found, skipping ligand extraction")
            return []

        generated = []

        for nonpoly in structure.nonpolymers:
            entity = structure.entities.get(nonpoly.entity_id)
            if entity is None:
                continue

            comp_id = entity.chemical_id
            if comp_id in SKIP_LIGANDS:
                continue

            auth_asym_id = nonpoly.auth_asym_id
            output_name = f"{self.rcsb_id}_{comp_id}_{auth_asym_id}.json"
            output_path = Path(self.assets.paths.base_dir) / output_name

            if output_path.exists() and not overwrite:
                print(f"  Skipping (exists): {output_name}")
                generated.append(output_path)
                continue

            print(f"  Extracting: {comp_id} / {auth_asym_id}")

            if self._run_extraction_script(
                cif_path, comp_id, auth_asym_id, output_path
            ):
                generated.append(output_path)

        return generated

    # -------------------------------------------------------------------------
    # Existing methods (unchanged)
    # -------------------------------------------------------------------------

    def _infer_organisms(
        self,
        entities: List[
            Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]
        ],
    ) -> dict[str, List]:
        all_src_ids = []
        all_host_ids = []
        all_src_names = []
        all_host_names = []

        for ent in entities:
            if hasattr(ent, "src_organism_ids"):
                all_src_ids.extend(ent.src_organism_ids)
            if hasattr(ent, "host_organism_ids"):
                all_host_ids.extend(ent.host_organism_ids)
            if hasattr(ent, "src_organism_names"):
                all_src_names.extend(ent.src_organism_names)
            if hasattr(ent, "host_organism_names"):
                all_host_names.extend(ent.host_organism_names)

        return {
            "src_organism_ids": sorted(list(set(all_src_ids))),
            "src_organism_names": sorted(list(set(all_src_names))),
            "host_organism_ids": sorted(list(set(all_host_ids))),
            "host_organism_names": sorted(list(set(all_host_names))),
        }

    # def _classify_tubulin_family(self, desc: str) -> Optional[TubulinFamily]:
    #     desc = (desc or "").lower()
    #     if "alpha" in desc or "tuba" in desc:
    #         return TubulinFamily.ALPHA
    #     if "beta" in desc or "tubb" in desc:
    #         return TubulinFamily.BETA
    #     if "gamma" in desc or "tubg" in desc:
    #         return TubulinFamily.GAMMA
    #     return None

    def _process_polymers(
        self, polymers_data: dict
    ) -> Tuple[
        Dict[str, Union[PolypeptideEntity, PolynucleotideEntity]],
        List[Polypeptide],
        List[Polynucleotide],
    ]:
        entity_map = {}
        polypeptides = []
        polynucleotides = []

        if not polymers_data:
            return {}, [], []

        raw_entities = polymers_data.get("polymer_entities") or []

        all_auth_ids = []
        for raw in raw_entities:
            ids = raw.get("rcsb_polymer_entity_container_identifiers") or {}
            all_auth_ids.extend(ids.get("auth_asym_ids") or [])
        assembly_lookup = self._get_assembly_mappings(all_auth_ids)

        for raw in raw_entities:
            entity_id = raw["rcsb_polymer_entity_container_identifiers"]["entity_id"]
            desc = raw.get("rcsb_polymer_entity", {}).get("pdbx_description")

            src_orgs = raw.get("rcsb_entity_source_organism", []) or []
            host_orgs = raw.get("rcsb_entity_host_organism", []) or []

            poly_type = raw["entity_poly"]["rcsb_entity_polymer_type"]
            strand_ids_str = raw["entity_poly"]["pdbx_strand_id"]
            strand_ids = (
                [s.strip() for s in strand_ids_str.split(",")] if strand_ids_str else []
            )

            entity_obj = None

            if poly_type == "Protein":
                uniprots = [
                    u.get("rcsb_id")
                    for u in (raw.get("uniprots") or [])
                    if u.get("rcsb_id")
                ]

                entity_obj = PolypeptideEntity(
                    entity_id=entity_id,
                    pdbx_description=desc,
                    pdbx_strand_ids=strand_ids,
                    formula_weight=raw["rcsb_polymer_entity"].get("formula_weight"),
                    one_letter_code=raw["entity_poly"]["pdbx_seq_one_letter_code"],
                    one_letter_code_can=raw["entity_poly"][
                        "pdbx_seq_one_letter_code_can"
                    ],
                    sequence_length=raw["entity_poly"]["rcsb_sample_sequence_length"],
                    src_organism_names=[o.get("scientific_name") for o in src_orgs],
                    src_organism_ids=[o.get("ncbi_taxonomy_id") for o in src_orgs],
                    host_organism_names=[o.get("scientific_name") for o in host_orgs],
                    host_organism_ids=[o.get("ncbi_taxonomy_id") for o in host_orgs],
                    family=None,
                    uniprot_accessions=uniprots,
                )
            else:
                entity_obj = PolynucleotideEntity(
                    entity_id=entity_id,
                    pdbx_description=desc,
                    pdbx_strand_ids=strand_ids,
                    polymer_type=poly_type,
                    one_letter_code=raw["entity_poly"]["pdbx_seq_one_letter_code"],
                    one_letter_code_can=raw["entity_poly"][
                        "pdbx_seq_one_letter_code_can"
                    ],
                    sequence_length=raw["entity_poly"]["rcsb_sample_sequence_length"],
                    src_organism_names=[o.get("scientific_name") for o in src_orgs],
                    src_organism_ids=[o.get("ncbi_taxonomy_id") for o in src_orgs],
                )

            entity_map[entity_id] = entity_obj

            identifiers = raw.get("rcsb_polymer_entity_container_identifiers", {})
            auth_asym_ids = identifiers.get("auth_asym_ids") or []
            asym_ids = identifiers.get("asym_ids") or []

            for auth_id, asym_id in zip(auth_asym_ids, asym_ids):
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

    def _process_nonpolymers(
        self, nonpoly_data: dict
    ) -> Tuple[Dict[str, NonpolymerEntity], List[Nonpolymer]]:
        entity_map = {}
        instances = []

        if not nonpoly_data:
            return {}, []
        raw_list = nonpoly_data.get("nonpolymer_entities") or []
        if not raw_list:
            return {}, []

        comp_ids = list(set(r["pdbx_entity_nonpoly"]["comp_id"] for r in raw_list))
        chem_info_map = self._fetch_chem_info(comp_ids)

        all_auth_ids = []
        for raw in raw_list:
            ids = raw.get("rcsb_nonpolymer_entity_container_identifiers") or {}
            all_auth_ids.extend(ids.get("auth_asym_ids") or [])
        assembly_lookup = self._get_assembly_mappings(all_auth_ids)

        for raw in raw_list:
            entity_id = raw["rcsb_nonpolymer_entity_container_identifiers"]["entity_id"]
            comp_id = raw["pdbx_entity_nonpoly"]["comp_id"]
            chem = chem_info_map.get(comp_id, {})

            nonpolymer_comp_data = raw.get("nonpolymer_comp")

            identifiers = raw.get("rcsb_nonpolymer_entity_container_identifiers", {})
            asym_ids = identifiers.get("asym_ids") or []

            entity_obj = NonpolymerEntity(
                entity_id=entity_id,
                chemical_id=comp_id,
                chemical_name=raw["pdbx_entity_nonpoly"]["name"],
                pdbx_description=raw["rcsb_nonpolymer_entity"]["pdbx_description"],
                formula_weight=raw["rcsb_nonpolymer_entity"].get("formula_weight"),
                pdbx_strand_ids=asym_ids,
                SMILES=chem.get("SMILES"),
                SMILES_stereo=chem.get("SMILES_stereo"),
                InChI=chem.get("InChI"),
                InChIKey=chem.get("InChIKey"),
                nonpolymer_comp=nonpolymer_comp_data,
            )
            entity_map[entity_id] = entity_obj

            auth_asym_ids = identifiers.get("auth_asym_ids") or []

            for auth_id, asym_id in zip(auth_asym_ids, asym_ids):
                asm_id = assembly_lookup.get(auth_id, 0)
                instances.append(
                    Nonpolymer(
                        parent_rcsb_id=self.rcsb_id,
                        auth_asym_id=auth_id,
                        asym_id=asym_id,
                        entity_id=entity_id,
                        assembly_id=asm_id,
                    )
                )

        return entity_map, instances

    def _fetch_chem_info(self, comp_ids: List[str]) -> dict:
        if not comp_ids:
            return {}
        try:
            q = LigandsChemInfo.replace("$COMP_IDS", json.dumps(comp_ids))
            data = query_rcsb_api(q)
            res = {}
            if data and "chem_comps" in data:
                for c in data["chem_comps"]:
                    res[c["chem_comp"]["id"]] = c.get("rcsb_chem_comp_descriptor", {})
            return res
        except Exception:
            return {}

    def _get_assembly_mappings(self, auth_asym_ids: List[str]) -> Dict[str, int]:
        if not self.asm_maps:
            return {}
        mapping = {}
        target_ids = set(auth_asym_ids)

        for asm in self.asm_maps:
            try:
                asm_idx = int(asm.rcsb_id.split("-")[1])
            except:
                asm_idx = 0

            if asm.polymer_entity_instances:
                for wrapper in asm.polymer_entity_instances:
                    identifier = wrapper.get(
                        "rcsb_polymer_entity_instance_container_identifiers"
                    )
                    if identifier and identifier.auth_asym_id in target_ids:
                        mapping[identifier.auth_asym_id] = asm_idx

            if asm.nonpolymer_entity_instances:
                for wrapper in asm.nonpolymer_entity_instances:
                    identifier = wrapper.get(
                        "rcsb_nonpolymer_entity_instance_container_identifiers"
                    )
                    if identifier:
                        aid = identifier.auth_asym_id
                        if aid in target_ids and aid not in mapping:
                            mapping[aid] = asm_idx
        return mapping

    def sequence_ingestion(self, entity: PolypeptideEntity) -> None:
        if not entity.family or entity.family not in [
            TubulinFamily.ALPHA,
            TubulinFamily.BETA,
        ]:
            return

        try:
            from api.services.seq_aligner import TubulinIngestor
            from dataclasses import asdict

            MUSCLE_BINARY = str(PROJECT_ROOT / "muscle3.8.1")

            if entity.family == TubulinFamily.ALPHA:
                MASTER_PROFILE = str(
                    PROJECT_ROOT / "data" / "alpha_tubulin" / "alpha_tubulin.afasta"
                )
                profile_type = "Alpha"
            else:
                MASTER_PROFILE = str(
                    PROJECT_ROOT / "data" / "beta_tubulin" / "beta_tubulin.afasta"
                )
                profile_type = "Beta"

            print(
                f"SEQUENCE INGESTION: {self.rcsb_id} ENTITY {entity.entity_id} ({profile_type})"
            )

            ingestor = TubulinIngestor(MASTER_PROFILE, MUSCLE_BINARY)

            result = ingestor.process_sequence(
                sequence=entity.one_letter_code_can,
                identifier=f"{self.rcsb_id}_{entity.entity_id}",
                family=entity.family.value,
            )

            pydantic_mutations = []
            for m in result.mutations:
                pydantic_mutations.append(
                    Mutation(
                        master_index=m.ma_position,
                        from_residue=m.wild_type,
                        to_residue=m.observed,
                        uniprot_id=entity.uniprot_accessions[0]
                        if entity.uniprot_accessions
                        else "UNKNOWN",
                        species=entity.src_organism_names[0]
                        if entity.src_organism_names
                        else "UNKNOWN",
                        tubulin_type=entity.family.value,
                        phenotype="Canonical Sequence",
                        database_source="PDB/Canonical",
                        reference_link=f"https://rcsb.org/structure/{self.rcsb_id}",
                        keywords="canonical_mutation",
                        notes=f"Detected at entity index {m.pdb_auth_id}",
                    )
                )

            entity.mutations = pydantic_mutations
            entity.alignment_stats = result.stats

            ingestion_file = Path(self.assets.paths.sequence_ingestion)
            ingestion_file.parent.mkdir(parents=True, exist_ok=True)

            if ingestion_file.exists():
                with open(ingestion_file, "r") as f:
                    data = json.load(f)
            else:
                data = {}

            data[entity.entity_id] = {
                "processed_at": datetime.now().isoformat(),
                "family": profile_type,
                "data": asdict(result),
            }

            with open(ingestion_file, "w") as f:
                json.dump(data, f, indent=2)

            if pydantic_mutations:
                print(
                    f"  > Entity {entity.entity_id}: Found {len(pydantic_mutations)} mutations."
                )

        except Exception as e:
            print(f"INGESTION ERROR Entity {entity.entity_id}: {e}")
            import traceback

            traceback.print_exc()

    def _sequence_ingestion_quiet(self, entity: PolypeptideEntity) -> None:
        """Same as sequence_ingestion but without verbose output."""
        if not entity.family or entity.family not in [TubulinFamily.ALPHA, TubulinFamily.BETA]:
            return

        try:
            from api.services.seq_aligner import TubulinIngestor
            from dataclasses import asdict

            MUSCLE_BINARY = str(PROJECT_ROOT / "muscle3.8.1")

            if entity.family == TubulinFamily.ALPHA:
                MASTER_PROFILE = str(PROJECT_ROOT / "data" / "alpha_tubulin" / "alpha_tubulin.afasta")
                profile_type = "Alpha"
            else:
                MASTER_PROFILE = str(PROJECT_ROOT / "data" / "beta_tubulin" / "beta_tubulin.afasta")
                profile_type = "Beta"

            ingestor = TubulinIngestor(MASTER_PROFILE, MUSCLE_BINARY)
            result = ingestor.process_sequence(
                sequence=entity.one_letter_code_can,
                identifier=f"{self.rcsb_id}_{entity.entity_id}",
                family=entity.family.value,
            )

            pydantic_mutations = []
            for m in result.mutations:
                pydantic_mutations.append(
                    Mutation(
                        master_index=m.ma_position,
                        from_residue=m.wild_type,
                        to_residue=m.observed,
                        uniprot_id=entity.uniprot_accessions[0] if entity.uniprot_accessions else "UNKNOWN",
                        species=entity.src_organism_names[0] if entity.src_organism_names else "UNKNOWN",
                        tubulin_type=entity.family.value,
                        phenotype="Canonical Sequence",
                        database_source="PDB/Canonical",
                        reference_link=f"https://rcsb.org/structure/{self.rcsb_id}",
                        keywords="canonical_mutation",
                        notes=f"Detected at entity index {m.pdb_auth_id}",
                    )
                )

            entity.mutations = pydantic_mutations
            entity.alignment_stats = result.stats

            ingestion_file = Path(self.assets.paths.sequence_ingestion)
            ingestion_file.parent.mkdir(parents=True, exist_ok=True)

            if ingestion_file.exists():
                with open(ingestion_file, "r") as f:
                    data = json.load(f)
            else:
                data = {}

            data[entity.entity_id] = {
                "processed_at": datetime.now().isoformat(),
                "family": profile_type,
                "data": asdict(result),
            }

            with open(ingestion_file, "w") as f:
                json.dump(data, f, indent=2)

        except Exception as e:
            pass  # Silently fail for quiet mode
    # -------------------------------------------------------------------------
    # Main entry point
    # -------------------------------------------------------------------------

    async def generate_profile(
        self,
        overwrite: bool = False,
        run_sequence_ingestion: bool = True,
        run_ligand_extraction: bool = True,
        run_classification: bool = True,
        ligand_extraction_workers: int = 4,
    ) -> TubulinStructure:
        profile_path = self.assets.paths.profile
        if os.path.exists(profile_path) and not overwrite:
            return self.assets.profile()

        print(f"\n{'─'*60}")
        print(f"Collecting {self.rcsb_id}")
        print(f"{'─'*60}")

        self.assets._verify_dir_exists()
        
        # CIF
        cif_path = Path(self.assets.paths.cif)
        if cif_path.exists():
            print(f"  ◆ CIF: cached")
        else:
            print(f"  ◆ CIF: downloading...")
            self.download_cif()

        # API queries
        print(f"  ◆ Fetching metadata from RCSB...")
        entry_data = query_rcsb_api(EntryInfoString.replace("$RCSB_ID", self.rcsb_id))["entry"]
        asm_data = query_rcsb_api(AssemblyIdentificationString.replace("$RCSB_ID", self.rcsb_id))
        self.asm_maps = [
            AssemblyInstancesMap.model_validate(a)
            for a in asm_data.get("entry", {}).get("assemblies", [])
        ]

        poly_data = query_rcsb_api(PolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id))["entry"]
        poly_entities, polypeptides, polynucleotides = self._process_polymers(poly_data)

        nonpoly_data = query_rcsb_api(NonpolymerEntitiesString.replace("$RCSB_ID", self.rcsb_id))["entry"]
        ligand_entities, nonpolymers = self._process_nonpolymers(nonpoly_data)

        all_entities = {**poly_entities, **ligand_entities}

        # Classification
        if run_classification:
            print(f"  ◆ HMM Classification:")
            self.run_classification(all_entities)

        # Sequence ingestion
        if run_sequence_ingestion:
            ingested = []
            for eid, ent in all_entities.items():
                if isinstance(ent, PolypeptideEntity):
                    if ent.family in [TubulinFamily.ALPHA, TubulinFamily.BETA]:
                        self._sequence_ingestion_quiet(ent)
                        if ent.mutations:
                            ingested.append(f"{eid}:{len(ent.mutations)} mut")
            if ingested:
                print(f"  ◆ Sequence Ingestion: {', '.join(ingested)}")

        organism_info = self._infer_organisms(list(all_entities.values()))

        citation = (entry_data.get("citation") or [{}])[0]
        external_refs = entry_data.get("rcsb_external_references", []) or []

        structure = TubulinStructure(
            rcsb_id=entry_data["rcsb_id"],
            expMethod=entry_data["exptl"][0]["method"],
            resolution=(entry_data["rcsb_entry_info"]["resolution_combined"] or [None])[0] or 0.0,
            deposition_date=entry_data["rcsb_accession_info"].get("deposit_date"),
            pdbx_keywords=entry_data.get("struct_keywords", {}).get("pdbx_keywords"),
            pdbx_keywords_text=entry_data.get("struct_keywords", {}).get("text"),
            rcsb_external_ref_id=[ref.get("id") for ref in external_refs],
            rcsb_external_ref_type=[ref.get("type") for ref in external_refs],
            rcsb_external_ref_link=[ref.get("link") for ref in external_refs],
            citation_title=citation.get("title"),
            citation_year=citation.get("year"),
            citation_rcsb_authors=citation.get("rcsb_authors"),
            citation_pdbx_doi=citation.get("pdbx_database_id_DOI"),
            **organism_info,
            entities=all_entities,
            polypeptides=polypeptides,
            polynucleotides=polynucleotides,
            nonpolymers=nonpolymers,
            assembly_map=self.asm_maps,
        )

        with open(profile_path, "w") as f:
            f.write(structure.model_dump_json(indent=4))

        # Ligand extraction
        if run_ligand_extraction:
            ligand_instances = []
            for nonpoly in structure.nonpolymers:
                entity = structure.entities.get(nonpoly.entity_id)
                if entity:
                    ligand_instances.append((entity.chemical_id, nonpoly.auth_asym_id))
            
            if ligand_instances:
                results = extract_ligands_parallel(
                    rcsb_id=self.rcsb_id,
                    cif_path=Path(self.assets.paths.cif),
                    ligand_instances=ligand_instances,
                    output_dir=Path(self.assets.paths.base_dir),
                    script_path=PROJECT_ROOT / "scripts_and_artifacts" / "extract_ixs.tsx",
                    project_root=PROJECT_ROOT,
                    max_workers=ligand_extraction_workers,
                    overwrite=overwrite,
                )
                succeeded = sum(1 for r in results if r.success)
                print(f"  ◆ Ligand Extraction: {succeeded}/{len(results)}")

        print(f"  ◆ Saved: {profile_path}")
        print(f"{'─'*60}\n")

        return structure


async def main(rcsb_id: str):
    c = TubulinETLCollector(rcsb_id)
    await c.generate_profile(overwrite=True)


if __name__ == "__main__":
    asyncio.run(main("6WVR"))
