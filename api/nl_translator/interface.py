"""Provider-agnostic interface for NL -> filter translation.

The translator's only responsibility is: take a natural-language request and
produce a validated Pydantic filter model (one of StructureFilters,
PolypeptideEntityFilters, LigandFilters) OR a clarification string.

Implementations MUST NOT touch the database.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, Protocol, Union, Dict, Any, List

from pydantic import BaseModel

from neo4j_tubxz.models import (
    StructureFilters,
    PolypeptideEntityFilters,
    LigandFilters,
)

Target = Literal["structures", "polymers", "ligands"]

FilterUnion = Union[StructureFilters, PolypeptideEntityFilters, LigandFilters]


@dataclass
class ViewContext:
    """Snapshot of the molstar viewer's current state, sent by the frontend.

    The LLM uses it to validate arguments (refuse unknown chains/ligands) and
    pick context-appropriate actions. Keep it small — goes into every prompt.
    """
    rcsb_id: Optional[str] = None
    chain_ids: List[str] = field(default_factory=list)
    ligand_keys: List[str] = field(default_factory=list)
    view_mode: Optional[str] = None  # "structure" | "monomer"
    active_monomer_chain: Optional[str] = None


@dataclass
class ViewerTranslationResult:
    """Result of one viewer-action translation call.

    - `actions`: validated Pydantic action instances, in emission order.
    - `clarification`: set instead of actions when the model can't commit.
    - `summary`: short human readback (empty when clarification is set).
    """
    actions: List[BaseModel] = field(default_factory=list)
    summary: str = ""
    clarification: Optional[str] = None


@dataclass
class FacetContext:
    """Snapshot of valid values the LLM can choose from.

    Rendered into the system prompt. Keep top-N values; full exhaustive lists
    blow up token budgets and aren't needed for the PoC canonical query set.
    """
    exp_methods: List[str] = field(default_factory=list)
    tubulin_families: List[str] = field(default_factory=list)
    isotypes: List[str] = field(default_factory=list)
    top_ligands: List[Dict[str, str]] = field(default_factory=list)  # [{chemical_id, chemical_name}]
    common_source_organisms: List[Dict[str, Any]] = field(default_factory=list)  # [{tax_id, name}]
    year_range: Dict[str, Optional[int]] = field(default_factory=dict)
    resolution_range: Dict[str, Optional[float]] = field(default_factory=dict)


@dataclass
class TranslationResult:
    """Result of one translation call.

    Exactly one of `filters` or `clarification` should be non-None.
    `summary` is a short plain-English readback of the chosen filters
    (empty when clarification is set).
    """
    target: Optional[Target] = None
    filters: Optional[FilterUnion] = None
    summary: str = ""
    clarification: Optional[str] = None


class NLTranslator(Protocol):
    """Implementations: api.nl_translator.anthropic_impl.AnthropicNLTranslator
    and (later) api.nl_translator.openai_compat_impl.OpenAICompatNLTranslator.
    """

    def translate(
        self,
        text: str,
        target: Target,
        facets: FacetContext,
        current_filters: Optional[Dict[str, Any]] = None,
    ) -> TranslationResult:
        ...

    def translate_viewer(
        self,
        text: str,
        view_context: ViewContext,
    ) -> ViewerTranslationResult:
        ...
