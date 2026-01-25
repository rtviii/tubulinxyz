# neo4j_tubxz/models.py
"""
Pydantic models for API request/response contracts.
"""

from typing import Any, Optional, List, Generic, TypeVar
from pydantic import BaseModel, Field, field_validator
from enum import Enum


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


class VariantTypeFilter(str, Enum):
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


class StructureFilters(BaseModel):
    """Filter parameters for Structure queries."""

    class Config:
        populate_by_name = True

    cursor: Optional[str] = Field(
        default=None, description="rcsb_id cursor for keyset pagination"
    )
    limit: int = Field(default=25, ge=1, le=500)

    search: Optional[str] = Field(default=None, description="Full-text search")
    rcsb_ids: Optional[List[str]] = Field(
        default=None, description="Filter to specific PDB IDs"
    )

    resolution_min: Optional[float] = None
    resolution_max: Optional[float] = None
    year_min: Optional[int] = None
    year_max: Optional[int] = None

    exp_method: Optional[List[ExpMethod]] = None
    polymerization_state: Optional[List[PolymerizationState]] = None

    source_organism_ids: Optional[List[int]] = Field(default=None, alias="sourceTaxa")
    host_organism_ids: Optional[List[int]] = None

    has_ligand_ids: Optional[List[str]] = None
    has_polymer_family: Optional[List[str]] = None
    has_uniprot: Optional[List[str]] = None

    # Variant filters (renamed from mutation)
    has_variants: Optional[bool] = Field(default=None, description="Has any variants")
    variant_family: Optional[str] = Field(
        default=None, description="Family to scope variant filters"
    )
    variant_type: Optional[VariantTypeFilter] = Field(
        default=None, description="Filter by variant type"
    )
    variant_position_min: Optional[int] = Field(
        default=None, description="Min master alignment position"
    )
    variant_position_max: Optional[int] = Field(
        default=None, description="Max master alignment position"
    )
    variant_wild_type: Optional[str] = Field(
        default=None, description="Wild-type residue"
    )
    variant_observed: Optional[str] = Field(
        default=None, description="Observed/mutant residue"
    )
    variant_source: Optional[str] = Field(
        default=None, description="Source: structural or literature"
    )
    variant_phenotype: Optional[str] = Field(
        default=None, description="Phenotype keyword search"
    )

    @field_validator("source_organism_ids", mode="before")
    @classmethod
    def decode_taxa_list(cls, v: Any):
        if isinstance(v, str):
            return [int(x.strip()) for x in v.split(",") if x.strip()]
        return v

    @field_validator(
        "rcsb_ids",
        "exp_method",
        "polymerization_state",
        "has_ligand_ids",
        "has_polymer_family",
        mode="before",
    )
    @classmethod
    def decode_strings(cls, v: Any):
        if isinstance(v, str):
            return [x.strip() for x in v.split(",") if x.strip()]
        return v


class PolypeptideEntityFilters(BaseModel):
    """Filter parameters for PolypeptideEntity queries."""

    cursor: Optional[str] = Field(
        default=None, description="Compound cursor 'RCSB_ID:ENTITY_ID'"
    )
    limit: int = Field(default=25, ge=1, le=500)

    structure_filters: Optional[StructureFilters] = None

    family: Optional[List[str]] = None
    uniprot_accession: Optional[str] = None
    sequence_contains: Optional[str] = None
    sequence_length_min: Optional[int] = None
    sequence_length_max: Optional[int] = None
    has_variants: Optional[bool] = Field(
        default=None, description="Filter to entities with/without variants"
    )


class LigandFilters(BaseModel):
    """Filter parameters for Chemical/Ligand queries."""

    cursor: Optional[str] = None
    limit: int = Field(default=25, ge=1, le=500)

    search: Optional[str] = Field(
        default=None, description="Search chemical name or ID"
    )
    chemical_ids: Optional[List[str]] = None
    has_drugbank: Optional[bool] = None
    in_structures: Optional[List[str]] = Field(
        default=None, description="Filter to ligands in these PDB IDs"
    )


T = TypeVar("T")


class PaginatedResponse(BaseModel, Generic[T]):
    """Standard paginated response wrapper."""

    data: List[T]
    total_count: int = Field(description="Total matching results before pagination")
    next_cursor: Optional[str] = Field(default=None, description="Cursor for next page")
    has_more: bool = Field(description="Whether more results exist")


class StructureSummary(BaseModel):
    """Lightweight structure representation for list views."""

    rcsb_id: str
    resolution: Optional[float] = None
    exp_method: Optional[str] = Field(None, alias="expMethod")
    citation_title: Optional[str] = None
    citation_year: Optional[int] = None
    deposition_date: Optional[str] = None
    src_organism_names: List[str] = []
    pdbx_keywords: Optional[str] = None
    entity_count: Optional[int] = None
    ligand_count: Optional[int] = None

    class Config:
        populate_by_name = True


class PolypeptideEntitySummary(BaseModel):
    """Lightweight entity representation."""

    parent_rcsb_id: str
    entity_id: str
    pdbx_description: Optional[str] = None
    family: Optional[str] = None
    sequence_length: Optional[int] = None
    src_organism_names: List[str] = []
    uniprot_accessions: List[str] = []
    variant_count: Optional[int] = None


class LigandSummary(BaseModel):
    """Lightweight ligand/chemical representation."""

    chemical_id: str
    chemical_name: Optional[str] = None
    drugbank_id: Optional[str] = None
    formula_weight: Optional[float] = None
    structure_count: Optional[int] = Field(
        None, description="Number of structures containing this ligand"
    )


class StructureListResponse(PaginatedResponse[StructureSummary]):
    pass


class PolypeptideListResponse(PaginatedResponse[PolypeptideEntitySummary]):
    pass


class LigandListResponse(PaginatedResponse[LigandSummary]):
    pass


# =============================================================================
# Facets Response Models
# =============================================================================


class FacetValue(BaseModel):
    """Single facet option with count."""

    value: str
    count: int


class LigandFacet(BaseModel):
    """Ligand facet for top ligands dropdown."""

    chemical_id: str
    chemical_name: Optional[str] = None
    count: int


class RangeValue(BaseModel):
    """Min/max range for numeric filters."""

    min: Optional[float] = None
    max: Optional[float] = None


class VariantsByFamily(BaseModel):
    """Variant stats by family."""

    family: str
    variant_count: int
    structure_count: int


class CommonVariant(BaseModel):
    """Common variant for quick filters."""

    family: Optional[str] = None
    position: Optional[int] = None
    wild_type: Optional[str] = None
    observed: Optional[str] = None
    variant_type: str
    count: int


class VariantPositionRange(BaseModel):
    """Variant position range by family."""

    family: str
    min_position: int
    max_position: int


class FilterFacets(BaseModel):
    """Available filter options for the UI."""

    total_structures: int
    exp_methods: List[FacetValue] = []
    tubulin_families: List[FacetValue] = []
    year_range: RangeValue
    resolution_range: RangeValue
    top_ligands: List[LigandFacet] = []
    variants_by_family: List[VariantsByFamily] = []
    common_variants: List[CommonVariant] = []
    variant_position_ranges: List[VariantPositionRange] = []


# =============================================================================
# Ligand Neighborhood Response Models
# =============================================================================


class BindingSiteResidue(BaseModel):
    """A residue in a ligand binding site."""

    auth_asym_id: str
    observed_index: int
    comp_id: str
    master_index: Optional[int] = None


class LigandNeighborhood(BaseModel):
    """Ligand neighborhood for a specific polymer chain."""

    ligand_id: str
    ligand_name: Optional[str] = None
    ligand_auth_asym_id: str
    residues: List[BindingSiteResidue] = []
    residue_count: int
    drugbank_id: Optional[str] = None


class PolymerNeighborhoodsResponse(BaseModel):
    """All ligand neighborhoods for a polymer chain."""

    rcsb_id: str
    auth_asym_id: str
    neighborhoods: List[LigandNeighborhood]
    total_ligands: int
    total_residues: int
