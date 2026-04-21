"""Anthropic-backed NL -> filter translator.

Uses tool-use: each filter model becomes a tool, Claude fills its arguments.
No database access. No conversation state.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional, Type

from pydantic import BaseModel, ValidationError

from neo4j_tubxz.models import (
    StructureFilters,
    PolypeptideEntityFilters,
    LigandFilters,
)
from api.nl_translator.interface import (
    NLTranslator,
    TranslationResult,
    Target,
    FacetContext,
    FilterUnion,
    ViewContext,
    ViewerTranslationResult,
)
from api.nl_translator.viewer_actions import (
    VIEWER_ACTION_MODELS,
    VIEWER_ACTION_DESCRIPTIONS,
    CLARIFY_ACTION_NAME,
)


# Map of tool name -> (Pydantic model, target label). The tool name is what
# Claude emits in `tool_use.name`; we dispatch on it.
_TOOL_MODELS: Dict[str, tuple[Type[BaseModel], Target]] = {
    "filter_structures": (StructureFilters, "structures"),
    "filter_polymers": (PolypeptideEntityFilters, "polymers"),
    "filter_ligands": (LigandFilters, "ligands"),
}

_CLARIFY_TOOL = "request_clarification"

_DEFAULT_MODEL = "claude-haiku-4-5"
_DEFAULT_MAX_TOKENS = 1024


def _model_to_tool(name: str, description: str, model: Type[BaseModel]) -> dict:
    """Convert a Pydantic model into an Anthropic tool spec."""
    schema = model.model_json_schema(by_alias=False)
    # Anthropic rejects some schema fields; strip the metadata Pydantic adds.
    schema.pop("title", None)
    return {
        "name": name,
        "description": description,
        "input_schema": schema,
    }


def _build_tools() -> List[dict]:
    return [
        _model_to_tool(
            "filter_structures",
            "Build a StructureFilters object for querying the structures list. "
            "Use when the user asks about whole PDB structures/entries.",
            StructureFilters,
        ),
        _model_to_tool(
            "filter_polymers",
            "Build a PolypeptideEntityFilters object for querying polymer/protein entities. "
            "Use when the user asks about chains/subunits/sequences, not whole structures.",
            PolypeptideEntityFilters,
        ),
        _model_to_tool(
            "filter_ligands",
            "Build a LigandFilters object for querying chemical ligands. "
            "Use when the user asks about small molecules/drugs themselves, not structures that contain them.",
            LigandFilters,
        ),
        {
            "name": _CLARIFY_TOOL,
            "description": (
                "Call this tool when the user's request is ambiguous and you cannot "
                "confidently fill a filter object. Provide a short clarification question."
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "question": {
                        "type": "string",
                        "description": "One-sentence clarification question to ask the user.",
                    }
                },
                "required": ["question"],
            },
        },
    ]


def _build_system_prompt(target: Target, facets: FacetContext) -> str:
    """Render the FacetContext into a system prompt with domain rules."""
    return f"""You translate natural-language requests about a tubulin structural biology database into structured filter objects.

YOUR OUTPUT: Exactly one tool call. Pick the filter tool whose model best fits the user's intent.
- The user is currently on the {target!r} view; prefer the matching filter tool unless the request clearly belongs elsewhere.
- If the request is ambiguous, call `request_clarification` instead of guessing.
- Do NOT emit text alongside the tool call unless asking for clarification.

CRITICAL DOMAIN CONVENTIONS:
- Resolution in Angstroms: LOWER values mean BETTER resolution. "Higher resolution than 4A" or "better than 4A" => resolution_max=4.0. "Worse than 4A" => resolution_min=4.0.
- "Human" maps to NCBI tax_id 9606. Other organism mappings are in the facet table below.
- Ligand names and PDB chemical IDs differ. "Taxol" and "paclitaxel" share PDB chem ID "TA1". Use facet table below to resolve.
- If the user specifies a range with the word "between", treat it inclusively.

VALID VALUES (use these verbatim where applicable):
- exp_methods: {facets.exp_methods}
- tubulin_families: {facets.tubulin_families}
- isotypes: {facets.isotypes}
- resolution_range (Angstroms, global): {facets.resolution_range}
- year_range (global): {facets.year_range}
- common source organisms (tax_id -> name): {facets.common_source_organisms}
- top ligands (chemical_id -> chemical_name): {facets.top_ligands}

FILTER OBJECT NOTES:
- Omit any field you don't need; do not set it to null explicitly.
- For integer lists (organism IDs): pass integers, not strings.
- Do not set cursor or limit; the backend controls pagination.

If a concept in the user's query does not map to any field in the chosen filter model, call `request_clarification` and describe what the user asked for that you could not express.
"""


class AnthropicNLTranslator:
    """Concrete NLTranslator using the Anthropic Python SDK."""

    def __init__(
        self,
        model: str = _DEFAULT_MODEL,
        max_tokens: int = _DEFAULT_MAX_TOKENS,
        api_key: Optional[str] = None,
    ) -> None:
        # Import here so the backend boots even if anthropic isn't installed
        # (e.g. when running with NL_TRANSLATOR_PROVIDER=openai_compat).
        from anthropic import Anthropic

        self._client = Anthropic(api_key=api_key or os.environ.get("ANTHROPIC_API_KEY"))
        self._model = model
        self._max_tokens = max_tokens
        self._tools = _build_tools()
        self._viewer_tools = _build_viewer_tools()
        self._viewer_models_by_name: Dict[str, Type[BaseModel]] = {
            m.__name__: m for m in VIEWER_ACTION_MODELS
        }

    def translate(
        self,
        text: str,
        target: Target,
        facets: FacetContext,
        current_filters: Optional[Dict[str, Any]] = None,
    ) -> TranslationResult:
        system_prompt = _build_system_prompt(target, facets)

        user_content = text.strip()
        if current_filters:
            user_content = (
                f"Current active filters (treat as context; merge or overwrite as the user intends): "
                f"{json.dumps(current_filters)}\n\nRequest: {text.strip()}"
            )

        response = self._client.messages.create(
            model=self._model,
            max_tokens=self._max_tokens,
            system=system_prompt,
            tools=self._tools,
            tool_choice={"type": "auto"},
            messages=[{"role": "user", "content": user_content}],
        )

        return self._interpret(response)

    def translate_viewer(
        self,
        text: str,
        view_context: ViewContext,
    ) -> ViewerTranslationResult:
        system_prompt = _build_viewer_system_prompt(view_context)

        response = self._client.messages.create(
            model=self._model,
            max_tokens=self._max_tokens,
            system=system_prompt,
            tools=self._viewer_tools,
            tool_choice={"type": "auto"},
            messages=[{"role": "user", "content": text.strip()}],
        )

        return _interpret_viewer_anthropic(response, self._viewer_models_by_name)

    def _interpret(self, response) -> TranslationResult:
        tool_blocks = [b for b in response.content if getattr(b, "type", None) == "tool_use"]
        text_blocks = [b for b in response.content if getattr(b, "type", None) == "text"]

        if not tool_blocks:
            # Claude responded with text only - treat it as a clarification.
            msg = " ".join(getattr(b, "text", "") for b in text_blocks).strip()
            return TranslationResult(
                clarification=msg or "I could not interpret that request. Please rephrase."
            )

        block = tool_blocks[0]
        name = block.name
        inp = block.input if isinstance(block.input, dict) else {}

        if name == _CLARIFY_TOOL:
            return TranslationResult(clarification=inp.get("question", "Please clarify."))

        if name not in _TOOL_MODELS:
            return TranslationResult(
                clarification=f"Internal: model invoked unknown tool {name!r}."
            )

        model_cls, chosen_target = _TOOL_MODELS[name]
        try:
            filters = model_cls.model_validate(inp)
        except ValidationError as e:
            return TranslationResult(
                clarification=(
                    "The parsed filters were not valid; please rephrase. "
                    f"Details: {e.errors(include_url=False)[:3]}"
                )
            )

        return TranslationResult(
            target=chosen_target,
            filters=filters,
            summary=_summarize(filters),
        )


def _build_viewer_tools() -> List[dict]:
    tools: List[dict] = []
    for model in VIEWER_ACTION_MODELS:
        name = model.__name__
        description = VIEWER_ACTION_DESCRIPTIONS.get(name, (model.__doc__ or "").strip())
        tools.append(_model_to_tool(name, description, model))
    return tools


def _build_viewer_system_prompt(view_context: ViewContext) -> str:
    chains = ", ".join(view_context.chain_ids) or "(none loaded)"
    ligands = ", ".join(view_context.ligand_keys) or "(none)"
    return f"""You are a controller for a 3D molecular viewer (molstar) displaying a single PDB structure.

Each user turn, emit ONE OR MORE tool calls that together fulfill the request. The frontend will execute them in order against the current viewer.

ACTIONS AVAILABLE (tool names): FocusChain, FocusResidue, FocusResidueRange, ClearFocus, SetChainVisibility, IsolateChain, HighlightChain, HighlightResidueRange, ClearHighlight, RequestClarification.

CURRENT VIEWER STATE:
- rcsb_id: {view_context.rcsb_id}
- loaded chains (auth_asym_id): {chains}
- loaded ligand keys: {ligands}
- view_mode: {view_context.view_mode}
- active_monomer_chain: {view_context.active_monomer_chain}

RULES:
- auth_asym_id arguments MUST be one of the loaded chain ids above. If the user names a chain that is not loaded, emit exactly one `RequestClarification` tool call instead, naming the available chains.
- Treat "residue 240" / "position 240" as auth_seq_id=240.
- For range expressions ("100-150", "between 100 and 150") include both endpoints.
- Prefer the minimum number of actions that accomplish the intent. Don't add ClearFocus / ClearHighlight unless the user implies reset.
- If the user asks to "hide everything except X" / "show only X" / "isolate X", use IsolateChain (not a sequence of SetChainVisibility).
- Do NOT emit any free-form text alongside tool calls.
"""


def _interpret_viewer_anthropic(response, viewer_models_by_name: Dict[str, Type[BaseModel]]) -> ViewerTranslationResult:
    tool_blocks = [b for b in response.content if getattr(b, "type", None) == "tool_use"]
    text_blocks = [b for b in response.content if getattr(b, "type", None) == "text"]

    if not tool_blocks:
        msg = " ".join(getattr(b, "text", "") for b in text_blocks).strip()
        return ViewerTranslationResult(
            clarification=msg or "I could not interpret that request. Please rephrase."
        )

    actions: List[BaseModel] = []
    clarification: Optional[str] = None

    for block in tool_blocks:
        name = block.name
        inp = block.input if isinstance(block.input, dict) else {}

        if name == CLARIFY_ACTION_NAME:
            clarification = inp.get("question", "Please clarify.")
            break

        model_cls = viewer_models_by_name.get(name)
        if model_cls is None:
            return ViewerTranslationResult(
                clarification=f"Internal: model invoked unknown viewer tool {name!r}."
            )
        try:
            actions.append(model_cls.model_validate(inp))
        except ValidationError as e:
            return ViewerTranslationResult(
                clarification=(
                    f"Arguments for {name} failed validation: "
                    f"{e.errors(include_url=False)[:2]}"
                )
            )

    if clarification is not None:
        return ViewerTranslationResult(clarification=clarification)

    return ViewerTranslationResult(
        actions=actions,
        summary=_summarize_actions(actions),
    )


def _summarize_actions(actions: List[BaseModel]) -> str:
    if not actions:
        return "(no actions)"
    parts = []
    for a in actions:
        kv = ", ".join(f"{k}={v}" for k, v in a.model_dump().items())
        parts.append(f"{type(a).__name__}({kv})" if kv else type(a).__name__)
    return " ; ".join(parts)


def _summarize(filters: FilterUnion) -> str:
    """Short plain-English readback of non-null filter fields."""
    d = filters.model_dump(exclude_none=True, exclude_defaults=True)
    # Drop pagination and nested structure_filters noise from the summary.
    for k in ("cursor", "limit"):
        d.pop(k, None)
    if not d:
        return "(no filters set)"
    pairs = []
    for k, v in d.items():
        if isinstance(v, list):
            pairs.append(f"{k}={', '.join(map(str, v))}")
        elif isinstance(v, dict):
            pairs.append(f"{k}={v}")
        else:
            pairs.append(f"{k}={v}")
    return "; ".join(pairs)
