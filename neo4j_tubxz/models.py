# neo4j_tubxz/models.py
"""
Pydantic models for API request/response contracts.
These define the typed interface between frontend and backend.
"""
from typing import Optional, List, Generic, TypeVar
from pydantic import BaseModel, Field
from enum import Enum


# -----------------------------------------------------------------------------
# Enums (for validated filter options)
# -----------------------------------------------------------------------------

class ExpMethod(str, Enum):
    XRAY = "X-RAY DIFFRACTION"
    EM = "ELECTRON MICROSCOPY"  
    NMR = "SOLUTION NMR"
    NEUTRON = "NEUTRON DIFFRACTION"
    FIBER = "FIBER DIFFRACTION"
    OTHER = "OTHER"


class PolymerizationState(str, Enum):
    MONOMER = "monomer"
    DIMER = "dimer"
    OLIGOMER = "oligomer"
    FILAMENT = "filament"
    UNKNOWN = "unknown"


# -----------------------------------------------------------------------------
# Filter Request Models
# -----------------------------------------------------------------------------

class StructureFilters(BaseModel):
    """
    Filter parameters for Structure queries.
    All filters are optional and cumulative (AND logic).
    """
    # Pagination
    cursor: Optional[str] = Field(
        default=None,
        description="rcsb_id cursor for keyset pagination. Returns results where rcsb_id < cursor."
    )
    limit: int = Field(
        default=25, 
        ge=1, 
        le=500,
        description="Max results per page"
    )
    
    # Text search (searches across multiple fields)
    search: Optional[str] = Field(
        default=None,
        description="Full-text search across title, keywords, rcsb_id, organism names"
    )
    
    # Exact match filters
    rcsb_ids: Optional[List[str]] = Field(
        default=None,
        description="Filter to specific PDB IDs"
    )
    
    # Range filters
    resolution_min: Optional[float] = Field(default=None, description="Min resolution in Å")
    resolution_max: Optional[float] = Field(default=None, description="Max resolution in Å")
    year_min: Optional[int] = Field(default=None, description="Min publication year")
    year_max: Optional[int] = Field(default=None, description="Max publication year")
    
    # Categorical filters  
    exp_method: Optional[List[ExpMethod]] = Field(
        default=None,
        description="Filter by experimental method(s)"
    )
    polymerization_state: Optional[List[PolymerizationState]] = Field(
        default=None,
        description="Filter by polymerization state"
    )
    
    # Taxonomy filters (via PhylogenyNode relationships)
    source_organism_ids: Optional[List[int]] = Field(
        default=None,
        description="NCBI taxonomy IDs for source organisms (includes descendants)"
    )
    host_organism_ids: Optional[List[int]] = Field(
        default=None,
        description="NCBI taxonomy IDs for host organisms (includes descendants)"
    )
    
    # Related entity filters
    has_ligand_ids: Optional[List[str]] = Field(
        default=None,
        description="Filter to structures containing specific chemical IDs (e.g., ['TAX', 'GTP'])"
    )
    has_polymer_family: Optional[List[str]] = Field(
        default=None,
        description="Filter to structures with specific tubulin families"
    )
    has_uniprot: Optional[List[str]] = Field(
        default=None,
        description="Filter to structures with specific UniProt accessions"
    )

    # Polymer/mutation filters
    has_mutations: Optional[bool] = Field(default=None, description="Has any mutations")
    mutation_family: Optional[str] = Field(default=None, description="Family to scope mutation filters (e.g., tubulin_beta)")
    mutation_position_min: Optional[int] = Field(default=None, description="Min master alignment position")
    mutation_position_max: Optional[int] = Field(default=None, description="Max master alignment position")
    mutation_from: Optional[str] = Field(default=None, description="Wild-type residue")
    mutation_to: Optional[str] = Field(default=None, description="Mutant residue")
    mutation_phenotype: Optional[str] = Field(default=None, description="Phenotype keyword search")


class PolypeptideEntityFilters(BaseModel):
    """
    Filter parameters for PolypeptideEntity queries.
    Can include structure-level filters to scope the search.
    """
    # Pagination (compound cursor: rcsb_id + entity_id)
    cursor: Optional[str] = Field(
        default=None,
        description="Compound cursor 'RCSB_ID:ENTITY_ID' for pagination"
    )
    limit: int = Field(default=25, ge=1, le=500)
    
    # Structure-level filters (inherited)
    structure_filters: Optional[StructureFilters] = Field(
        default=None,
        description="Apply structure-level filters first"
    )
    
    # Entity-specific filters
    family: Optional[List[str]] = Field(
        default=None,
        description="Tubulin family filter"
    )
    uniprot_accession: Optional[str] = Field(
        default=None,
        description="Filter by UniProt accession"
    )
    sequence_contains: Optional[str] = Field(
        default=None,
        description="Filter by sequence motif"
    )
    sequence_length_min: Optional[int] = None
    sequence_length_max: Optional[int] = None
    has_mutations: Optional[bool] = Field(
        default=None,
        description="Filter to entities with/without mutations"
    )


class LigandFilters(BaseModel):
    """Filter parameters for Chemical/Ligand queries"""
    cursor: Optional[str] = None
    limit: int = Field(default=25, ge=1, le=500)
    
    search: Optional[str] = Field(
        default=None,
        description="Search chemical name or ID"
    )
    chemical_ids: Optional[List[str]] = None
    has_drugbank: Optional[bool] = None
    
    # Filter to ligands present in specific structures
    in_structures: Optional[List[str]] = Field(
        default=None,
        description="Filter to ligands present in these PDB IDs"
    )


# -----------------------------------------------------------------------------
# Response Models
# -----------------------------------------------------------------------------

T = TypeVar('T')


class PaginatedResponse(BaseModel, Generic[T]):
    """
    Standard paginated response wrapper.
    Frontend RTK-Query can use this shape directly.
    """
    data: List[T]
    total_count: int = Field(description="Total matching results (before pagination)")
    next_cursor: Optional[str] = Field(
        default=None,
        description="Cursor for next page. Null if no more results."
    )
    has_more: bool = Field(description="Whether more results exist")


class StructureSummary(BaseModel):
    """Lightweight structure representation for list views"""
    rcsb_id: str
    resolution: Optional[float] = None
    exp_method: Optional[str] = Field(None, alias="expMethod")
    citation_title: Optional[str] = None
    citation_year: Optional[int] = None
    deposition_date: Optional[str] = None
    src_organism_names: List[str] = []
    pdbx_keywords: Optional[str] = None
    
    # Counts for UI badges
    entity_count: Optional[int] = None
    ligand_count: Optional[int] = None
    
    class Config:
        populate_by_name = True


class PolypeptideEntitySummary(BaseModel):
    """Lightweight entity representation"""
    parent_rcsb_id: str
    entity_id: str
    pdbx_description: Optional[str] = None
    family: Optional[str] = None
    sequence_length: Optional[int] = None
    src_organism_names: List[str] = []
    uniprot_accessions: List[str] = []
    mutation_count: Optional[int] = None


class LigandSummary(BaseModel):
    """Lightweight ligand/chemical representation"""
    chemical_id: str
    chemical_name: Optional[str] = None
    drugbank_id: Optional[str] = None
    formula_weight: Optional[float] = None
    structure_count: Optional[int] = Field(
        None,
        description="Number of structures containing this ligand"
    )


# Concrete response types for the API
class StructureListResponse(PaginatedResponse[StructureSummary]):
    pass


class PolypeptideListResponse(PaginatedResponse[PolypeptideEntitySummary]):
    pass


class LigandListResponse(PaginatedResponse[LigandSummary]):
    pass


# Add to neo4j_tubxz/models.py

class FacetValue(BaseModel):
    """Single facet option with count"""
    value: str
    count: int


class LigandFacet(BaseModel):
    """Ligand facet for top ligands dropdown"""
    chemical_id: str
    chemical_name: Optional[str] = None
    count: int


class RangeValue(BaseModel):
    """Min/max range for numeric filters"""
    min: Optional[float] = None
    max: Optional[float] = None


class MutationByFamily(BaseModel):
    family: str
    mutation_count: int
    structure_count: int


class CommonMutation(BaseModel):
    family: Optional[str] = None
    position: int
    from_residue: str
    to_residue: str
    count: int


class MutationPositionRange(BaseModel):
    family: str
    min_position: int
    max_position: int


class FilterFacets(BaseModel):
    total_structures: int
    exp_methods: List[FacetValue]
    tubulin_families: List[FacetValue]
    year_range: RangeValue
    resolution_range: RangeValue
    top_ligands: List[LigandFacet]
    mutations_by_family: List[MutationByFamily]
    common_mutations: List[CommonMutation]
    mutation_position_ranges: List[MutationPositionRange]