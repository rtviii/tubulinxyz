# tubexyz/lib/schema/types_tubulin.py
from typing import Dict, Optional, List, Literal, Union, Any
from enum import Enum
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, Field

class TubulinFamily(str, Enum):
    ALPHA   = "alpha"
    BETA    = "beta"
    GAMMA   = "gamma"
    DELTA   = "delta"
    EPSILON = "epsilon"
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

    family: Optional[TubulinFamily] = None
    uniprot_accessions: List[str] = []

    # Mutations detected vis a vis the Master Alignmnet
    mutations      : List[Mutation] = []
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




# --- 2. INSTANCES (The Physical Objects) ---



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




# ------
