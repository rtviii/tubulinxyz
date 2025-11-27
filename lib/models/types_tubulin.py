# tubexyz/lib/schema/types_tubulin.py
from typing import Dict, Optional, List, Literal
from enum import Enum
import typing
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, Field

# --- Enums ---

class TubulinFamily(str, Enum):
    ALPHA   = "alpha"
    BETA    = "beta"
    GAMMA   = "gamma"
    DELTA   = "delta"
    EPSILON = "epsilon"
    # ZETA    = "zeta"
    # ETA     = "eta"
    # THETA   = "theta"
    # IOTA    = "iota"
    # KAPPA   = "kappa"

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


class Polymer(BaseModel):
    """Base polymer model (Instance-centric)"""
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    def to_SeqRecord(self) -> SeqRecord:
        return SeqRecord(
            seq         = Seq(self.entity_poly_seq_one_letter_code_can),
            id          = f"{self.src_organism_ids[0]}",
            description = f'{self.parent_rcsb_id}.{self.auth_asym_id}',
            name        = f'{self.parent_rcsb_id}.{self.auth_asym_id}'
        )

    assembly_id     : int       
    auth_asym_id    : str      
    parent_rcsb_id  : str
    entity_id       : str     
    asym_ids        : List[str] 

    src_organism_names : List[str]
    host_organism_names: List[str]
    src_organism_ids   : List[int]
    host_organism_ids  : List[int]

    rcsb_pdbx_description: Optional[str] = None

    entity_poly_strand_id           : str
    entity_poly_seq_one_letter_code : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length          : int
    entity_poly_polymer_type        : str
    entity_poly_entity_type         : str

class TubulinProtein(Polymer):
    """Tubulin-specific protein with family classification"""
    
    # Tubulin-specific field
    family: Optional[TubulinFamily] = None
    
    pfam_accessions  : List[str] = []
    pfam_comments    : List[str] = []
    pfam_descriptions: List[str] = []
    uniprot_accession: List[str] = []

class MasterAlignment(BaseModel):
    """Versioned canonical reference MSA for a tubulin family"""
    version: str  
    family: TubulinFamily
    fasta_content: str
    created_date: str 
    description: Optional[str] = None


class Mutation(BaseModel):
    """
    Mutation from literature/clinical databases"""
    
    master_index: int
    utn_position: Optional[int] = None
    
    from_residue: str  # Wild-type
    to_residue  : str  # Mutant
    
    uniprot_id  : str
    species     : str
    tubulin_type: str
    
    phenotype      : str
    database_source: str
    reference_link : str
    keywords       : str
    notes          : Optional[str] = None

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

            class DrugbankContainerIdentifiers(BaseModel):
                drugbank_id: str

            drugbank_container_identifiers: Optional[DrugbankContainerIdentifiers] = None
            drugbank_info: Optional[DrugbankInfo] = None

        class RcsbChemCompTarget(BaseModel):

            interaction_type                 : Optional[str] = None
            name                             : Optional[str] = None
            provenance_source                : Optional[str] = None
            reference_database_accession_code: Optional[str] = None
            reference_database_name          : Optional[str] = None

        drugbank: Optional[Drugbank] = None
        rcsb_chem_comp_target: Optional[List[RcsbChemCompTarget]] = None

    chemicalId         : str
    chemicalName       : str
    formula_weight     : Optional[float] = None
    pdbx_description   : str
    number_of_instances: int
    nonpolymer_comp    : Optional[NonpolymerComp] = None
    
    SMILES       : Optional[str] = None
    SMILES_stereo: Optional[str] = None
    InChI        : Optional[str] = None
    InChIKey     : Optional[str] = None

class RCSBStructureMetadata(BaseModel):
    """Standard RCSB structure metadata"""
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    
    rcsb_id: str
    expMethod: str
    resolution: float
    deposition_date: Optional[str] = None
    
    pdbx_keywords: Optional[str] = None
    pdbx_keywords_text: Optional[str] = None
    
    rcsb_external_ref_id  : List[str]
    rcsb_external_ref_type: List[str]
    rcsb_external_ref_link: List[str]
    
    citation_year        : Optional[int]       = None
    citation_rcsb_authors: Optional[List[str]] = None
    citation_title       : Optional[str]       = None
    citation_pdbx_doi    : Optional[str]       = None
    
    src_organism_ids   : List[int]
    src_organism_names : List[str]
    host_organism_ids  : List[int]
    host_organism_names: List[str]
    
    assembly_map: Optional[List["AssemblyInstancesMap"]] = None

class TubulinStructure(RCSBStructureMetadata):
    """Complete tubulin structure with all components"""
    proteins: List[TubulinProtein]
    other_polymers: List[Polymer]  # For any non-tubulin chains
    nonpolymeric_ligands: List[NonpolymericLigand]
    
    polymerization_state: Optional[Literal["monomer", "dimer", "oligomer", "filament", "unknown"]] = None

class AssemblyInstancesMap(BaseModel):
    """Maps polymer/ligand instances to assemblies"""
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    
    class NonpolymerEntityInstance(BaseModel):
        class NonpolymerEntityInstanceContainerIdentifiers(BaseModel):
            entity_id: str
            auth_asym_id: str
            auth_seq_id: str
        rcsb_nonpolymer_entity_instance_container_identifiers: NonpolymerEntityInstanceContainerIdentifiers

    class PolymerEntityInstance(BaseModel):
        class PolymerEntityInstanceContainerIdentifiers(BaseModel):
            entity_id: str
            auth_asym_id: str
        rcsb_polymer_entity_instance_container_identifiers: PolymerEntityInstanceContainerIdentifiers

    rcsb_id: str
    nonpolymer_entity_instances: Optional[List[NonpolymerEntityInstance]] = None
    polymer_entity_instances: List[PolymerEntityInstance]


# --- Relationship Models ---

class AlignmentMapping(BaseModel):
    seqres_to_master: str  # JSON array: list[int] (seqres_idx -> master_idx | -1)
    master_to_seqres: str  # JSON array: list[int] (master_idx -> seqres_idx | -1)