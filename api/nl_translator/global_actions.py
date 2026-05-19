"""Pydantic models for the /nl_query/global endpoint.

The "global" endpoint sits on the landing page and accepts open-ended user
questions ("where does taxol bind", "what kinds of PTMs are out there"). It
produces a single response containing:

- An optional short blurb (LLM-generated, hard-capped).
- Zero or more `queries` (filter specs the frontend will run via existing
  list endpoints; cards may reference them by id).
- A small ranked list of `ActionCard` objects — one click each, navigating
  the user to a page with state preloaded via URL params.

The schema is intentionally flat. Both `EntityRef` and `ActionCard` use a
discriminator field (`kind` / `action`) plus a wide set of optional payload
fields, instead of an anyOf-of-nested-types. Haiku handles flat schemas with
nullable fields much more reliably than deeply nested discriminated unions.
Backend code interprets which fields are meaningful for each kind/action.
"""
from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import BaseModel, Field

from neo4j_tubxz.models import (
    StructureFilters,
    PolypeptideEntityFilters,
    LigandFilters,
)


# ---------------------------------------------------------------------------
# Entity references — what the LLM can point at in the database.
# ---------------------------------------------------------------------------

EntityKind = Literal[
    "structure",
    "chain",
    "polymer_entity",
    "family",
    "ligand",
    "variant",
    "residue_range",
]


class EntityRef(BaseModel):
    """A reference to one domain entity. Fields are optional and the backend
    interprets only those relevant for `kind`. Validation runs in hydration.py.
    """
    kind: EntityKind

    # Structure / chain / polymer_entity
    rcsb_id: Optional[str] = Field(default=None, description="PDB id, uppercase")
    auth_asym_id: Optional[str] = Field(default=None, description="Chain id, e.g. 'A'")
    entity_id: Optional[str] = Field(default=None, description="Polymer entity id, e.g. '1'")

    # Ligand
    chemical_id: Optional[str] = Field(default=None, description="PDB chemical component id, e.g. 'TA1'")
    auth_seq_id: Optional[int] = Field(default=None, description="Ligand instance position (when bound)")

    # Family
    family: Optional[str] = Field(
        default=None,
        description="One of the tubulin_families literal values, e.g. 'tubulin_beta'.",
    )

    # Variant
    master_index: Optional[int] = Field(default=None, description="Position in master alignment")
    wild_type: Optional[str] = Field(default=None, description="Single-letter WT residue")
    observed: Optional[str] = Field(default=None, description="Single-letter observed residue")

    # Residue range
    start: Optional[int] = Field(default=None, description="Range start (inclusive, auth_seq_id)")
    end: Optional[int] = Field(default=None, description="Range end (inclusive, auth_seq_id)")


# ---------------------------------------------------------------------------
# Action cards — what the LLM tells the UI to render.
# ---------------------------------------------------------------------------

ActionKind = Literal[
    "open_catalogue",
    "open_structure",
    "open_expert",
    "inspect_ligand",
    "view_variants",
    "clarify",
]


class AlignedRef(BaseModel):
    """One chain to load into the MSA alongside the primary chain."""
    rcsb_id: str
    auth_asym_id: str


class RangeRef(BaseModel):
    """A residue range (inclusive)."""
    start: int
    end: int


class ActionCard(BaseModel):
    """One clickable card. The frontend builds a URL from `action` + payload
    and routes the user when clicked.
    """
    action: ActionKind
    label: str = Field(..., max_length=80, description="Short chip text shown to user")
    description: Optional[str] = Field(
        default=None,
        max_length=160,
        description="Optional 1-line subtitle under the label",
    )

    # open_catalogue
    query_ref: Optional[str] = Field(
        default=None,
        description="Id of an entry in queries[] whose filters drive this catalogue view",
    )

    # open_structure / open_expert / inspect_ligand
    rcsb_id: Optional[str] = None
    focus_chains: Optional[List[str]] = None      # open_structure
    focus_ligands: Optional[List[str]] = None     # open_structure

    # open_expert
    primary_chain: Optional[str] = Field(
        default=None, description="auth_asym_id of the monomer-view primary chain"
    )
    aligned: Optional[List[AlignedRef]] = Field(
        default=None,
        description="Additional chains to load into the MSA alongside the primary",
    )
    focus_range: Optional[RangeRef] = None

    # inspect_ligand
    chemical_id: Optional[str] = None
    suggested_chain: Optional[str] = Field(
        default=None,
        description="Chain that contacts this ligand in the chosen structure",
    )

    # view_variants  (family also reused by open_catalogue as a filter)
    family: Optional[str] = None
    position_min: Optional[int] = None
    position_max: Optional[int] = None
    variant_type: Optional[str] = Field(
        default=None,
        description="substitution | insertion | deletion",
    )

    # open_catalogue direct-filter shortcuts (when query_ref isn't worth it).
    # Mirror the filter fields on UiFilters / StructureFilters.
    source_organism_ids: Optional[List[int]] = Field(
        default=None,
        description="NCBI tax ids to filter the catalogue by (e.g. [5811] for Toxoplasma gondii). Use when you want a catalogue card scoped to an organism without spelling out a full query in queries[].",
    )

    # clarify
    question: Optional[str] = Field(
        default=None, description="One-sentence clarification when intent is ambiguous"
    )


# ---------------------------------------------------------------------------
# Query specs the frontend will run against existing list endpoints.
# ---------------------------------------------------------------------------

QueryTarget = Literal["structures", "polymers", "ligands"]


class QuerySpec(BaseModel):
    """One filter query the frontend will execute via the existing list
    endpoints. `id` is referenced by cards' `query_ref`.
    """
    id: str = Field(..., pattern=r"^q\d+$", description="Stable id, e.g. 'q1'")
    target: QueryTarget
    filters_structures: Optional[StructureFilters] = None
    filters_polymers: Optional[PolypeptideEntityFilters] = None
    filters_ligands: Optional[LigandFilters] = None


# ---------------------------------------------------------------------------
# Top-level response envelope (the single tool the LLM calls).
# ---------------------------------------------------------------------------

MAX_BLURB_CHARS = 180
MAX_CARDS = 6
MAX_QUERIES = 3


class GlobalNLResponse(BaseModel):
    """The single tool call the LLM emits for a global NL query.

    `validation` is populated server-side after hydration; the LLM never sets
    it. Kept on the same envelope so the frontend sees one shape.
    """
    blurb: str = Field(
        default="",
        max_length=MAX_BLURB_CHARS,
        description=(
            "Optional 1-sentence interpretation/summary; <=180 chars. "
            "Leave empty if cards speak for themselves."
        ),
    )
    queries: List[QuerySpec] = Field(default_factory=list, max_length=MAX_QUERIES)
    cards: List[ActionCard] = Field(default_factory=list, max_length=MAX_CARDS)
    # Server-populated; LLM should not set this.
    validation: Dict[str, Dict[str, Any]] = Field(default_factory=dict)


# Tool registration metadata used by the translator.
EMIT_GLOBAL_TOOL_NAME = "emit_global_response"
CLARIFY_GLOBAL_TOOL_NAME = "request_clarification_global"


__all__ = [
    "EntityKind",
    "EntityRef",
    "ActionKind",
    "AlignedRef",
    "RangeRef",
    "ActionCard",
    "QueryTarget",
    "QuerySpec",
    "GlobalNLResponse",
    "MAX_BLURB_CHARS",
    "MAX_CARDS",
    "MAX_QUERIES",
    "EMIT_GLOBAL_TOOL_NAME",
    "CLARIFY_GLOBAL_TOOL_NAME",
]
