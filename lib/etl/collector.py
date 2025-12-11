import asyncio
import json
import os
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path
from datetime import datetime
from pydantic import ValidationError

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
from lib.models.types_tubulin import (
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


class TubulinETLCollector:
    rcsb_id: str
    assets: TubulinStructureAssets
    asm_maps: List[AssemblyInstancesMap] = []

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id.upper()
        self.assets = TubulinStructureAssets(self.rcsb_id)

    def _classify_tubulin_family(self, desc: str) -> Optional[TubulinFamily]:
        desc = (desc or "").lower()
        if "alpha" in desc or "tuba" in desc:
            return TubulinFamily.ALPHA
        if "beta" in desc or "tubb" in desc:
            return TubulinFamily.BETA
        if "gamma" in desc or "tubg" in desc:
            return TubulinFamily.GAMMA
        return None

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
            # --- A. Build Entity ---
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
                fam = self._classify_tubulin_family(desc)
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
                    family=fam,
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

            # --- B. Build Instances ---
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

            # 1. CAPTURE NONPOLYMER_COMP
            nonpolymer_comp_data = raw.get("nonpolymer_comp")  # Dict from GQL

            identifiers = raw.get("rcsb_nonpolymer_entity_container_identifiers", {})
            asym_ids = identifiers.get("asym_ids") or []

            entity_obj = NonpolymerEntity(
                entity_id=entity_id,
                chemicalId=comp_id,
                chemicalName=raw["pdbx_entity_nonpoly"]["name"],
                pdbx_description=raw["rcsb_nonpolymer_entity"]["pdbx_description"],
                formula_weight=raw["rcsb_nonpolymer_entity"].get("formula_weight"),
                pdbx_strand_ids=asym_ids,
                SMILES=chem.get("SMILES"),
                SMILES_stereo=chem.get("SMILES_stereo"),
                InChI=chem.get("InChI"),
                InChIKey=chem.get("InChIKey"),
                # 2. ASSIGN IT (Pydantic will parse the dict)
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
        """
        Runs alignment logic ONCE per Entity.
        Saves results to 'sequence_ingestion.json' and updates Entity object.
        """
        if not entity.family or entity.family not in [
            TubulinFamily.ALPHA,
            TubulinFamily.BETA,
        ]:
            return

        try:
            from api.services.alignment import TubulinIngestor
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

            # Map mutations for in-memory object
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

            # Store on Entity (for DB builder if it uses this object directly)
            entity.mutations = pydantic_mutations
            entity.alignment_stats = result.stats

            # 3. SAVE TO sequence_ingestion.json (The "Clean" way)
            ingestion_file = Path(self.assets.paths.sequence_ingestion)

            # Read existing
            if ingestion_file.exists():
                with open(ingestion_file, "r") as f:
                    data = json.load(f)
            else:
                data = {}

            # Update
            data[entity.entity_id] = {
                "processed_at": datetime.now().isoformat(),
                "family": profile_type,
                "data": asdict(result),  # Save the full dataclass result
            }

            # Write
            with open(ingestion_file, "w") as f:
                json.dump(data, f, indent=2)

            if pydantic_mutations:
                print(
                    f"  > Entity {entity.entity_id}: Found {len(pydantic_mutations)} mutations."
                )

        except Exception as e:
            print(f"⚠️ INGESTION ERROR Entity {entity.entity_id}: {e}")
            import traceback

            traceback.print_exc()

    async def generate_profile(
        self, overwrite: bool = False, run_sequence_ingestion: bool = True
    ) -> TubulinStructure:
        profile_path = self.assets.paths.profile
        if os.path.exists(profile_path) and not overwrite:
            return self.assets.profile()

        print(f"Generating profile for {self.rcsb_id}...")

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

        # 4. RUN SEQUENCE INGESTION (On Entities)
        if run_sequence_ingestion:
            for eid, ent in all_entities.items():
                if isinstance(ent, PolypeptideEntity):
                    self.sequence_ingestion(ent)

        citation = (entry_data.get("citation") or [{}])[0]
        external_refs = entry_data.get("rcsb_external_references", []) or []

        structure = TubulinStructure(
            rcsb_id=entry_data["rcsb_id"],
            expMethod=entry_data["exptl"][0]["method"],
            resolution=(entry_data["rcsb_entry_info"]["resolution_combined"] or [None])[
                0
            ]
            or 0.0,
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
            entities=all_entities,
            polypeptides=polypeptides,
            polynucleotides=polynucleotides,
            nonpolymers=nonpolymers,
            assembly_map=self.asm_maps,
        )

        self.assets._verify_dir_exists()

        # 5. SAVE PROFILE (EXCLUDING HEAVY MUTATIONS)
        # We exclude 'mutations' and 'alignment_stats' from all entities in the dictionary
        with open(profile_path, "w") as f:
            f.write(
                structure.model_dump_json(
                    indent=4,
                    exclude={"entities": {"__all__": {"mutations", "alignment_stats"}}},
                )
            )

        print(f"Saved: {profile_path}")
        return structure


async def main(rcsb_id: str):
    c = TubulinETLCollector(rcsb_id)
    await c.generate_profile(overwrite=True)


if __name__ == "__main__":
    asyncio.run(main("6WVR"))
