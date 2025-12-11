# tubexyz/lib/schema/types_tubulin.py
from typing import Dict, Optional, List, Literal, Union, Any
from enum import Enum
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, Field

# --- Enums ---


class TubulinFamily(str, Enum):
    ALPHA = "alpha"
    BETA = "beta"
    GAMMA = "gamma"
    DELTA = "delta"
    EPSILON = "epsilon"


class MutationType(str, Enum):
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


class ModificationType(str, Enum):
    ACETYLATION = "acetylation"
    PHOSPHORYLATION = "phosphorylation"
    METHYLATION = "methylation"
    UBIQUITINATION = "ubiquitination"
    SUMOYLATION = "sumoylation"
    PALMITOYLATION = "palmitoylation"
    NITROSYLATION = "nitrosylation"
    GLUTAMYLATION = "glutamylation"
    GLYCYLATION = "glycylation"
    TYROSINATION = "tyrosination"
    DETYROSINATION = "detyrosination"


# --- Support Models ---


class Mutation(BaseModel):
    master_index: int
    utn_position: Optional[int] = None
    from_residue: str
    to_residue: str

    uniprot_id: str
    species: str
    tubulin_type: str
    phenotype: str
    database_source: str
    reference_link: str
    keywords: str
    notes: Optional[str] = None


# --- Complex Property Models (Restored) ---


class NonpolymerComp(BaseModel):
    class Drugbank(BaseModel):
        class DrugbankInfo(BaseModel):
            cas_number: Optional[str] = None
            description: Optional[str] = None
            indication: Optional[str] = None
            mechanism_of_action: Optional[str] = None
            name: Optional[str] = None
            pharmacology: Optional[str] = None

        class DrugbankContainerIdentifiers(BaseModel):
            drugbank_id: str

        drugbank_container_identifiers: Optional[DrugbankContainerIdentifiers] = None
        drugbank_info: Optional[DrugbankInfo] = None

    class RcsbChemCompTarget(BaseModel):
        interaction_type: Optional[str] = None
        name: Optional[str] = None
        provenance_source: Optional[str] = None
        reference_database_accession_code: Optional[str] = None
        reference_database_name: Optional[str] = None

    drugbank: Optional[Drugbank] = None
    rcsb_chem_comp_target: Optional[List[RcsbChemCompTarget]] = None


# --- 1. ENTITIES (The Blueprints) ---


class BaseEntity(BaseModel):
    entity_id: str
    type: Literal["polymer", "non-polymer", "water", "branched"]
    pdbx_description: Optional[str] = None
    formula_weight: Optional[float] = None
    pdbx_strand_ids: List[str] = []


class PolypeptideEntity(BaseEntity):
    type: Literal["polymer"] = "polymer"
    polymer_type: Literal["Protein"] = "Protein"

    one_letter_code: str
    one_letter_code_can: str
    sequence_length: int

    src_organism_names: List[str] = []
    host_organism_names: List[str] = []
    src_organism_ids: List[int] = []
    host_organism_ids: List[int] = []

    family: Optional[TubulinFamily] = None
    uniprot_accessions: List[str] = []

    # Ingestion Results
    # (These will be populated in memory but excluded from profile.json save)
    mutations: List[Mutation] = []
    alignment_stats: Dict[str, Any] = {}

    def to_SeqRecord(self, rcsb_id: str) -> SeqRecord:
        return SeqRecord(
            seq=Seq(self.one_letter_code_can),
            id=f"{rcsb_id}_{self.entity_id}",
            description=self.pdbx_description or "",
            name=f"{rcsb_id}_entity_{self.entity_id}",
        )


class PolynucleotideEntity(BaseEntity):
    type: Literal["polymer"] = "polymer"
    polymer_type: str  # "DNA", "RNA", "NA-hybrid"

    one_letter_code: str
    one_letter_code_can: str
    sequence_length: int

    src_organism_names: List[str] = []
    src_organism_ids: List[int] = []


class NonpolymerEntity(BaseEntity):
    type: Literal["non-polymer"] = "non-polymer"

    chemicalId: str
    chemicalName: str

    SMILES: Optional[str] = None
    SMILES_stereo: Optional[str] = None
    InChI: Optional[str] = None
    InChIKey: Optional[str] = None

    # RESTORED: The rich annotation object
    nonpolymer_comp: Optional[NonpolymerComp] = None


# --- 2. INSTANCES (The Physical Objects) ---


class BaseInstance(BaseModel):
    parent_rcsb_id: str
    auth_asym_id: str
    asym_id: str
    entity_id: str
    assembly_id: int

    def __hash__(self):
        return hash(self.asym_id + self.parent_rcsb_id)


class Polypeptide(BaseInstance):
    pass


class Polynucleotide(BaseInstance):
    pass


class Nonpolymer(BaseInstance):
    pass


# --- 3. STRUCTURE ROOT ---


class AssemblyInstancesMap(BaseModel):
    class InstanceIdentifier(BaseModel):
        entity_id: str
        auth_asym_id: str
        asym_id: Optional[str] = None

    rcsb_id: str
    nonpolymer_entity_instances: Optional[List[Dict[str, InstanceIdentifier]]] = None
    polymer_entity_instances: List[Dict[str, InstanceIdentifier]]


class RCSBStructureMetadata(BaseModel):
    model_config = {"json_encoders": {Enum: lambda v: v.value}}
    rcsb_id: str
    expMethod: str
    resolution: float
    deposition_date: Optional[str] = None

    pdbx_keywords: Optional[str] = None
    pdbx_keywords_text: Optional[str] = None

    rcsb_external_ref_id: List[str]
    rcsb_external_ref_type: List[str]
    rcsb_external_ref_link: List[str]

    citation_year: Optional[int] = None
    citation_rcsb_authors: Optional[List[str]] = None
    citation_title: Optional[str] = None
    citation_pdbx_doi: Optional[str] = None


class TubulinStructure(RCSBStructureMetadata):
    entities: Dict[
        str, Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]
    ]

    polypeptides: List[Polypeptide]
    polynucleotides: List[Polynucleotide]
    nonpolymers: List[Nonpolymer]

    assembly_map: Optional[List[AssemblyInstancesMap]] = None
    polymerization_state: Optional[
        Literal["monomer", "dimer", "oligomer", "filament", "unknown"]
    ] = None
