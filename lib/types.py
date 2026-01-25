# lib/types.py
"""
Core type definitions for the tubulin ETL pipeline and database.
"""

import json
from typing import Dict, Optional, List, Literal, Union, Any
from enum import Enum
from pydantic import BaseModel, Field


# ============================================================
# Enums
# ============================================================


class VariantType(str, Enum):
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


class TubulinFamily(str, Enum):
    ALPHA = "tubulin_alpha"
    BETA = "tubulin_beta"
    GAMMA = "tubulin_gamma"
    DELTA = "tubulin_delta"
    EPSILON = "tubulin_epsilon"


class MapFamily(str, Enum):
    """Microtubule Associated Proteins (MAPs) and related enzymes."""

    ATAT1 = "map_atat1"
    CAMSAP1 = "map_camsap1"
    CAMSAP2 = "map_camsap2"
    CAMSAP3 = "map_camsap3"
    CCP = "map_ccp_deglutamylase"
    CFAP53 = "map_cfap53"
    CKAP5 = "map_ckap5_chtog"
    CLASP = "map_clasp"
    CLIP115 = "map_clip115"
    CLIP170 = "map_clip170"
    DOUBLECORTIN = "map_doublecortin"
    EB_FAMILY = "map_eb_family"
    FAP20 = "map_fap20_cfap20"
    GCP2_3 = "map_gcp2_3"
    GCP4 = "map_gcp4"
    GCP5_6 = "map_gcp5_6"
    KATANIN = "map_katanin_p60"
    KINESIN13 = "map_kinesin13"
    MAP1_HEAVY = "map_map1_heavy"
    MAP1S = "map_map1s"
    MAP2 = "map_map2"
    MAP4 = "map_map4"
    MAP7 = "map_map7"
    NME7 = "map_nme7"
    NME8 = "map_nme8"
    NUMA = "map_numa"
    PACRG = "map_pacrg"
    PRC1 = "map_prc1"
    RIB72 = "map_rib72_efhc"
    SPAG6 = "map_spag6"
    SPASTIN = "map_spastin"
    STATHMIN = "map_stathmin"
    TACC = "map_tacc"
    TAU = "map_tau"
    TPX2 = "map_tpx2"
    TTLL_LONG = "map_ttll_glutamylase_long"
    TTLL_SHORT = "map_ttll_glutamylase_short"
    VASH = "map_vash_detyrosinase"


PolymerClass = Union[TubulinFamily, MapFamily]


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


# ============================================================
# Index Mapping
# ============================================================


class IndexMappingData(BaseModel):
    """Bidirectional mapping between observed (auth_seq_id) and master alignment indices."""

    observed_to_master: Dict[
        int, Optional[int]
    ]  # auth_seq_id -> MA index (1-based) or None
    master_to_observed: Dict[
        int, Optional[int]
    ]  # MA index (1-based) -> auth_seq_id or None

    def get_master_index(self, auth_seq_id: int) -> Optional[int]:
        return self.observed_to_master.get(auth_seq_id)

    def get_observed_index(self, master_index: int) -> Optional[int]:
        return self.master_to_observed.get(master_index)

    def to_json_dict(self) -> Dict[str, str]:
        """Serialize for Neo4j storage as JSON strings."""
        return {
            "observed_to_master_json": json.dumps(
                {str(k): v for k, v in self.observed_to_master.items()}
            ),
            "master_to_observed_json": json.dumps(
                {str(k): v for k, v in self.master_to_observed.items()}
            ),
        }

    @classmethod
    def from_json_dict(cls, data: Dict[str, str]) -> "IndexMappingData":
        """Deserialize from Neo4j JSON strings."""
        return cls(
            observed_to_master={
                int(k): v
                for k, v in json.loads(data["observed_to_master_json"]).items()
            },
            master_to_observed={
                int(k): v
                for k, v in json.loads(data["master_to_observed_json"]).items()
            },
        )


# ============================================================
# Sequence Variants (Substitutions, Insertions, Deletions)
# ============================================================


class SequenceVariant(BaseModel):
    """
    A sequence variant detected by alignment or from literature.
    Covers substitutions, insertions, and deletions.
    """

    type: VariantType

    # Source of this annotation: "structural" (from alignment) or "literature"
    source: str = "structural"

    # For substitutions and deletions: the MA position (1-based)
    # For insertions: None
    master_index: Optional[int] = None

    # For substitutions and insertions: the auth_seq_id in the structure
    # For deletions: None
    observed_index: Optional[int] = None

    # For substitutions: wild_type -> observed
    # For insertions: observed only
    # For deletions: wild_type only
    wild_type: Optional[str] = None
    observed: Optional[str] = None

    # Optional metadata (primarily for literature-sourced variants)
    uniprot_id: Optional[str] = None
    phenotype: Optional[str] = None
    reference: Optional[str] = None  # DOI, PMID, or database link

    @classmethod
    def substitution(
        cls,
        master_index: int,
        observed_index: int,
        wild_type: str,
        observed: str,
        source: str = "structural",
        **kwargs,
    ):
        return cls(
            type=VariantType.SUBSTITUTION,
            master_index=master_index,
            observed_index=observed_index,
            wild_type=wild_type,
            observed=observed,
            source=source,
            **kwargs,
        )

    @classmethod
    def insertion(
        cls, observed_index: int, residue: str, source: str = "structural", **kwargs
    ):
        return cls(
            type=VariantType.INSERTION,
            observed_index=observed_index,
            observed=residue,
            source=source,
            **kwargs,
        )

    @classmethod
    def deletion(
        cls, master_index: int, expected: str, source: str = "structural", **kwargs
    ):
        return cls(
            type=VariantType.DELETION,
            master_index=master_index,
            wild_type=expected,
            source=source,
            **kwargs,
        )


# ============================================================
# Ligand Binding Sites
# ============================================================


class BindingSiteResidue(BaseModel):
    """A residue in a ligand's binding site."""

    auth_asym_id: str
    observed_index: int  # auth_seq_id
    comp_id: str  # 3-letter residue code
    master_index: Optional[int] = None  # 1-based MA index, if mapped

    def to_dict(self) -> Dict[str, Any]:
        """For Neo4j storage on relationship properties."""
        return {
            "auth_asym_id": self.auth_asym_id,
            "observed_index": self.observed_index,
            "comp_id": self.comp_id,
            "master_index": self.master_index,
        }


class LigandBindingSite(BaseModel):
    """Binding site for a single ligand instance."""

    ligand_comp_id: str
    ligand_auth_asym_id: str
    ligand_auth_seq_id: int
    residues: List[BindingSiteResidue]

    @property
    def residue_count(self) -> int:
        return len(self.residues)

    def residues_for_chain(self, auth_asym_id: str) -> List[BindingSiteResidue]:
        """Get residues belonging to a specific chain."""
        return [r for r in self.residues if r.auth_asym_id == auth_asym_id]


# ============================================================
# Modifications (PTMs) - for future use
# ============================================================


class Modification(BaseModel):
    """Post-translational modification from literature/databases."""

    master_index: int
    amino_acid: str
    modification_type: str

    uniprot_id: str
    species: str
    tubulin_type: str

    phenotype: str
    database_source: str
    database_link: str
    keywords: str
    notes: Optional[str] = None


# ============================================================
# Chemical / Nonpolymer Types
# ============================================================


class NonpolymerComp(BaseModel):
    """Drugbank and target info for a chemical."""

    class Drugbank(BaseModel):
        class DrugbankInfo(BaseModel):
            cas_number: Optional[str] = None
            description: Optional[str] = None

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


# ============================================================
# Instance Types (physical copies in a structure)
# ============================================================


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
    auth_seq_id: int  # Add this
    pass


# ============================================================
# Entity Types (blueprints - one per unique sequence/chemical per structure)
# ============================================================


class BaseEntity(BaseModel):
    entity_id: str
    type: Literal["polymer", "non-polymer", "water", "branched"]
    pdbx_description: Optional[str] = None
    formula_weight: Optional[float] = None
    pdbx_strand_ids: List[str] = []


class PolypeptideEntity(BaseModel):
    """Polypeptide entity with index mappings and variants."""

    type: Literal["polymer"] = "polymer"
    polymer_type: Literal["Protein"] = "Protein"
    entity_id: str
    pdbx_description: Optional[str] = None
    pdbx_strand_ids: List[str] = []

    one_letter_code: str
    one_letter_code_can: str
    sequence_length: int

    src_organism_names: List[str] = []
    host_organism_names: List[str] = []
    src_organism_ids: List[int] = []
    host_organism_ids: List[int] = []

    family: Optional[PolymerClass] = None
    uniprot_accessions: List[str] = []

    # Index mappings (for classified tubulin chains)
    index_mapping: Optional[IndexMappingData] = None

    # Sequence variants
    variants: List[SequenceVariant] = []

    # Alignment stats
    alignment_stats: Dict[str, Any] = {}

    @property
    def substitutions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.SUBSTITUTION]

    @property
    def insertions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.INSERTION]

    @property
    def deletions(self) -> List[SequenceVariant]:
        return [v for v in self.variants if v.type == VariantType.DELETION]


class PolynucleotideEntity(BaseEntity):
    type: Literal["polymer"] = "polymer"
    polymer_type: str
    one_letter_code: str
    one_letter_code_can: str
    sequence_length: int
    src_organism_names: List[str] = []
    src_organism_ids: List[int] = []


class NonpolymerEntity(BaseEntity):
    type: Literal["non-polymer"] = "non-polymer"

    chemical_id: str
    chemical_name: str

    pdbx_description: Optional[str] = None
    formula_weight: Optional[float] = None

    nonpolymer_comp: Optional[NonpolymerComp] = None

    SMILES: Optional[str] = None
    SMILES_stereo: Optional[str] = None
    InChI: Optional[str] = None
    InChIKey: Optional[str] = None

    num_instances: int = 0


# ============================================================
# Structure Root
# ============================================================


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

    src_organism_ids: List[int] = []
    src_organism_names: List[str] = []
    host_organism_ids: List[int] = []
    host_organism_names: List[str] = []


class TubulinStructure(RCSBStructureMetadata):
    entities: Dict[
        str, Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]
    ]

    polypeptides: List[Polypeptide]
    polynucleotides: List[Polynucleotide]
    nonpolymers: List[Nonpolymer]

    assembly_map: Optional[List[AssemblyInstancesMap]] = None
    ligand_binding_sites: List[LigandBindingSite] = []

    polymerization_state: Optional[
        Literal["monomer", "dimer", "oligomer", "filament", "unknown"]
    ] = None


# ============================================================
# Molstar Extraction Types
# ============================================================


class ObservedResidue(BaseModel):
    auth_seq_id: int
    label_seq_id: int
    comp_id: str
    one_letter: str


class ObservedSequenceData(BaseModel):
    auth_asym_id: str
    entity_id: str
    residues: List[ObservedResidue]

    @property
    def sequence(self) -> str:
        return "".join(r.one_letter for r in self.residues)

    @property
    def auth_seq_ids(self) -> List[int]:
        return [r.auth_seq_id for r in self.residues]


class MolstarExtractionResult(BaseModel):
    """Complete extraction result from Molstar for a structure."""

    rcsb_id: str
    sequences: List[ObservedSequenceData]
    ligand_neighborhoods: List[LigandBindingSite]

    def get_sequence_for_chain(
        self, auth_asym_id: str
    ) -> Optional[ObservedSequenceData]:
        for seq in self.sequences:
            if seq.auth_asym_id == auth_asym_id:
                return seq
        return None


# ============================================================
# Classification Report
# ============================================================


class EntityClassificationResult(BaseModel):
    """Classification result for a single entity."""

    entity_id: str
    auth_asym_ids: List[str]
    sequence_length: int
    assigned_family: Optional[str]
    best_hit_score: Optional[float] = None
    best_hit_evalue: Optional[float] = None
    all_hits: List[Dict[str, Any]] = []


class ClassificationReport(BaseModel):
    """File: {RCSB_ID}_classification_report.json"""

    rcsb_id: str
    generated_at: str
    summary: Dict[str, int]
    entities: Dict[str, EntityClassificationResult]


# ============================================================
# Output File Schemas
# ============================================================


class VariantsFile(BaseModel):
    """File: {RCSB_ID}_variants.json"""

    rcsb_id: str
    generated_at: str
    entities: Dict[str, List[SequenceVariant]]


class LigandBindingSitesFile(BaseModel):
    """File: {RCSB_ID}_ligand_binding_sites.json"""

    rcsb_id: str
    generated_at: str
    binding_sites: List[LigandBindingSite]
