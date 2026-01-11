# types_tubulin.py
from typing import Dict, Optional, List, Literal, Tuple, Union, Any
from enum import Enum
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, Field

class TubulinFamily(str, Enum):
    ALPHA   = "tubulin_alpha"
    BETA    = "tubulin_beta"
    GAMMA   = "tubulin_gamma"
    DELTA   = "tubulin_delta"
    EPSILON = "tubulin_epsilon"

class MapFamily(str, Enum):
    """
    Microtubule Associated Proteins (MAPs) and related enzymes.
    Values correspond to the filename base.
    """
    ATAT1           = "map_atat1"
    CAMSAP1         = "map_camsap1"
    CAMSAP2         = "map_camsap2"
    CAMSAP3         = "map_camsap3"
    CCP             = "map_ccp_deglutamylase"
    CFAP53          = "map_cfap53"
    CKAP5           = "map_ckap5_chtog"
    CLASP           = "map_clasp"
    CLIP115         = "map_clip115"
    CLIP170         = "map_clip170"
    DOUBLECORTIN    = "map_doublecortin"
    EB_FAMILY       = "map_eb_family"
    FAP20           = "map_fap20_cfap20"
    GCP2_3          = "map_gcp2_3"
    GCP4            = "map_gcp4"
    GCP5_6          = "map_gcp5_6"
    KATANIN         = "map_katanin_p60"
    KINESIN13       = "map_kinesin13"
    MAP1_HEAVY      = "map_map1_heavy"
    MAP1S           = "map_map1s"
    MAP2            = "map_map2"
    MAP4            = "map_map4"
    MAP7            = "map_map7"
    NME7            = "map_nme7"
    NME8            = "map_nme8"
    NUMA            = "map_numa"
    PACRG           = "map_pacrg"
    PRC1            = "map_prc1"
    RIB72           = "map_rib72_efhc"
    SPAG6           = "map_spag6"
    SPASTIN         = "map_spastin"
    STATHMIN        = "map_stathmin"
    TACC            = "map_tacc"
    TAU             = "map_tau"
    TPX2            = "map_tpx2"
    TTLL_LONG       = "map_ttll_glutamylase_long"
    TTLL_SHORT      = "map_ttll_glutamylase_short"
    VASH            = "map_vash_detyrosinase"

PolymerClass = Union[TubulinFamily, MapFamily]

# --- Ligand Interaction Models ---

class InteractionType(str, Enum):

    UNKNOWN            = "Unknown"
    IONIC              = "Ionic"
    CATION_PI          = "Cation-Pi Interaction"
    PI_STACKING        = "Pi Stacking"
    HYDROGEN_BOND      = "Hydrogen Bond"
    HALOGEN_BOND       = "Halogen Bond"
    HYDROPHOBIC        = "Hydrophobic Contact"
    METAL_COORDINATION = "Metal Coordination"
    WEAK_HYDROGEN_BOND = "Weak Hydrogen Bond"

class InteractionParticipant(BaseModel):
    """An atom participating in an interaction."""
    auth_asym_id: str
    auth_seq_id : int
    auth_comp_id: str
    atom_id     : str
    is_ligand   : bool
    master_index: Optional[int] = None  # Added field, made optional

    @classmethod
    def from_tuple(cls, t: list) -> "InteractionParticipant":
        return cls(
            auth_asym_id=t[0],
            auth_seq_id=t[1],
            auth_comp_id=t[2],
            atom_id=t[3],
            is_ligand=t[4],
            # Safely handle both 5-tuples and 6-tuples
            master_index=t[5] if len(t) > 5 else None 
        )
    
    def to_tuple(self) -> list:
        """Helper to convert back to the list format for JSON storage."""
        base = [
            self.auth_asym_id,
            self.auth_seq_id,
            self.auth_comp_id,
            self.atom_id,
            self.is_ligand
        ]
        if self.master_index is not None:
            base.append(self.master_index)
        return base

class LigandInteraction(BaseModel):
    """A single interaction between ligand and polymer."""
    type: str
    participants: Tuple[InteractionParticipant, InteractionParticipant]

    @classmethod
    def from_raw(cls, raw: dict) -> "LigandInteraction":
        return cls(
            type=raw["type"],
            participants=(
                InteractionParticipant.from_tuple(raw["participants"][0]),
                InteractionParticipant.from_tuple(raw["participants"][1]),
            ),
        )

class NeighborResidue(BaseModel):
    """A residue in the ligand's neighborhood."""
    auth_asym_id: str
    auth_seq_id: int
    auth_comp_id: str
    master_index: Optional[int] = None # Added

    @classmethod
    def from_tuple(cls, t: list) -> "NeighborResidue":
        return cls(
            auth_asym_id=t[0], 
            auth_seq_id=t[1], 
            auth_comp_id=t[2],
            master_index=t[3] if len(t) > 3 else None
        )
    
    def to_tuple(self) -> list:
        base = [self.auth_asym_id, self.auth_seq_id, self.auth_comp_id]
        if self.master_index is not None:
            base.append(self.master_index)
        return base

class LigandNeighborhood(BaseModel):
    ligand_auth_asym_id: str
    ligand_auth_seq_id: int
    ligand_comp_id: str
    interactions: List[LigandInteraction]
    neighborhood: List[NeighborResidue]

    @classmethod
    def from_raw(cls, raw: Any) -> "LigandNeighborhood":
        # If it's the raw list from Molstar (TSX), take the first entry
        if isinstance(raw, list):
            if not raw:
                raise ValueError("Empty ligand data list.")
            raw = raw[0]
            
        # The schema uses the "ligand" key which is [auth_asym_id, auth_seq_id, comp_id]
        ligand_info = raw["ligand"]
        return cls(
            ligand_auth_asym_id = ligand_info[0],
            ligand_auth_seq_id  = ligand_info[1],
            ligand_comp_id      = ligand_info[2],
            interactions        = [LigandInteraction.from_raw(i) for i in raw["interactions"]],
            neighborhood        = [NeighborResidue.from_tuple(n) for n in raw["neighborhood"]],
        )


# --- Enums ---
class MasterAlignment(BaseModel):
    """Versioned canonical reference MSA for a tubulin family"""
    version      : str
    family       : TubulinFamily
    fasta_content: str
    created_date : str
    description  : Optional[str] = None

class AlignmentMapping(BaseModel):
    seqres_to_master: str  # JSON array: list[int] (seqres_idx -> master_idx | -1)
    master_to_seqres: str  # JSON array: list[int] (master_idx -> seqres_idx | -1)



class MutationType(str, Enum):
    SUBSTITUTION = "substitution"
    INSERTION    = "insertion"
    DELETION     = "deletion"


class ModificationType(str, Enum):
    ACETYLATION     = "acetylation"
    PHOSPHORYLATION = "phosphorylation"
    METHYLATION     = "methylation"
    UBIQUITINATION  = "ubiquitination"
    SUMOYLATION     = "sumoylation"
    PALMITOYLATION  = "palmitoylation"
    NITROSYLATION   = "nitrosylation"
    GLUTAMYLATION   = "glutamylation"
    GLYCYLATION     = "glycylation"
    TYROSINATION    = "tyrosination"
    DETYROSINATION  = "detyrosination"


# --- Support Models ---

class Modification(BaseModel):
    """Post-translational modification from literature/databases"""
    
    master_index: int
    utn_position: Optional[int] = None
    
    amino_acid       : str
    modification_type: str
    
    uniprot_id  : str
    species     : str
    tubulin_type: str
    
    phenotype      : str
    database_source: str
    database_link  : str
    keywords       : str
    notes          : Optional[str] = None

class NonpolymericLigand(BaseModel):
    """Ligand model - unchanged from riboxyz"""
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    
    class NonpolymerComp(BaseModel):
        class Drugbank(BaseModel):
            class DrugbankInfo(BaseModel):
                cas_number: Optional[str] = None
                description: Optional[str] = None

class Mutation(BaseModel):
    master_index: int
    utn_position: Optional[int] = None
    from_residue: str
    to_residue: str

    uniprot_id     : str
    species        : str
    tubulin_type   : TubulinFamily
    phenotype      : str
    database_source: str
    reference_link : str
    keywords       : str
    notes          : Optional[str] = None


class NonpolymerComp(BaseModel):
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

# This one is for db only.
class ChemicalCompound(BaseModel):
    """
    Global chemical identity - shared across all structures.
    One instance per unique chemical_id in the entire PDB.
    """
    chemical_id  : str  # Primary key: "TAX", "PMM", "GTP"
    chemical_name: str
    
    SMILES        : Optional[str]   = None
    SMILES_stereo : Optional[str]   = None
    InChI         : Optional[str]   = None
    InChIKey      : Optional[str]   = None
    formula_weight: Optional[float] = None
    
    nonpolymer_comp: Optional[NonpolymerComp] = None



class BaseInstance(BaseModel):
    parent_rcsb_id: str
    auth_asym_id  : str
    asym_id       : str
    entity_id     : str
    assembly_id   : int

    def __hash__(self):
        return hash(self.asym_id + self.parent_rcsb_id)


class BaseEntity(BaseModel):
    entity_id       : str
    type            : Literal["polymer", "non-polymer", "water", "branched"]
    pdbx_description: Optional[str] = None
    formula_weight  : Optional[float] = None
    pdbx_strand_ids : List[str] = []




class NonpolymerEntity(BaseEntity):
    type: Literal["non-polymer"] = "non-polymer"
    
    # Reference to the shared chemical
    chemical_id: str
    chemical_name: str
    
    pdbx_description: Optional[str] = None
    formula_weight: Optional[float] = None
    
    # --- MISSING FIELDS ADDED BELOW ---
    # These are needed so node_ligand.py can create the Global Chemical Node
    nonpolymer_comp: Optional[NonpolymerComp] = None
    
    SMILES        : Optional[str] = None
    SMILES_stereo : Optional[str] = None
    InChI         : Optional[str] = None
    InChIKey      : Optional[str] = None

    num_instances : int = 0

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

    family: Optional[PolymerClass] = None  # Changed from Optional[TubulinFamily]
    uniprot_accessions: List[str] = []

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

    type               : Literal["polymer"] = "polymer"
    polymer_type       : str                             
    one_letter_code    : str
    one_letter_code_can: str
    sequence_length    : int
    src_organism_names: List[str] = []
    src_organism_ids  : List[int] = []




class Polypeptide(BaseInstance):
    pass


class Polynucleotide(BaseInstance):
    pass


class Nonpolymer(BaseInstance):
    """
    Represents a single nonpolymer instance (one molecule copy) in the structure.
    
    Inherits from BaseInstance:

        - parent_rcsb_id: str - The PDB ID (e.g., "6WVR")
        - auth_asym_id  : str - Author chain ID (e.g., "E")
        - asym_id       : str - Internal chain ID (e.g., "E")
        - entity_id     : str - References the NonpolymerEntity (e.g., "4")
        - assembly_id   : int - Which biological assembly (e.g., 1)
    
    Example:
        Three Taxol molecules in 6WVR would create three Nonpolymer instances,
        all pointing to the same NonpolymerEntity (which points to the same
        ChemicalCompound).
    """
    
    # Currently no additional fields beyond BaseInstance
    # Could add instance-specific data later if needed:
    # occupancy: Optional[float] = None
    # b_factor: Optional[float] = None




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
    
    # --- RESTORED FIELDS ---
    src_organism_ids   : List[int] = []
    src_organism_names : List[str] = []
    host_organism_ids  : List[int] = []
    host_organism_names: List[str] = []


class TubulinStructure(RCSBStructureMetadata):
    entities: Dict[
        str, Union[PolypeptideEntity, PolynucleotideEntity, NonpolymerEntity]
    ]

    polypeptides   : List[Polypeptide]
    polynucleotides: List[Polynucleotide]
    nonpolymers    : List[Nonpolymer]

    assembly_map: Optional[List[AssemblyInstancesMap]] = None
    polymerization_state: Optional[
        Literal["monomer", "dimer", "oligomer", "filament", "unknown"]
    ] = None




class MutationEntryData(BaseModel):
    """Mutation entry as stored in sequence ingestion results."""
    ma_position: int
    wild_type: str
    observed: str
    pdb_auth_id: int


class ProcessedChainData(BaseModel):
    """Result of sequence alignment/ingestion for a single entity."""
    pdb_id: str
    chain_id: str
    tubulin_class: str
    sequence: str
    
    # ma_to_auth_map[ma_idx] = auth_seq_id (or -1 if missing)
    ma_to_auth_map: List[int]
    
    # observed_to_ma_map[obs_idx] = MA position (1-based) or -2 for insertions
    observed_to_ma_map: List[int]
    
    mutations: List[MutationEntryData]
    stats: Dict[str, Any]


class SequenceIngestionEntry(BaseModel):
    """A single entity's ingestion record."""
    processed_at: str
    family: str
    data: ProcessedChainData
    
    def build_auth_to_ma_map(self) -> Dict[int, int]:
        """
        Build reverse lookup: auth_seq_id -> MA index (1-based).
        
        The ma_to_auth_map stores: ma_to_auth_map[ma_idx] = auth_seq_id
        We invert this to get: auth_seq_id -> ma_idx + 1
        """
        auth_to_ma: Dict[int, int] = {}
        for ma_idx, auth_id in enumerate(self.data.ma_to_auth_map):
            if auth_id != -1:
                auth_to_ma[auth_id] = ma_idx + 1  # MA positions are 1-based
        return auth_to_ma