"""OpenAI-compatible NLTranslator.

Works with any endpoint that implements the OpenAI /v1/chat/completions API
with tool-calling. Verified targets:
- OpenRouter (openrouter.ai)  -- hosted Claude, GPT, Llama, etc.
- vLLM / SGLang / Ollama / TGI self-hosted open-weights models (Qwen, Llama).

Config (env vars):
- OPENAI_API_KEY       -- required. For OpenRouter, your sk-or-v1-... key.
- OPENAI_BASE_URL      -- default: https://openrouter.ai/api/v1
                          For self-hosted vLLM, e.g. http://gpu-box:8000/v1
- OPENAI_MODEL         -- default: anthropic/claude-haiku-4-5   (OpenRouter id)
                          For vLLM/Qwen, e.g. Qwen/Qwen2.5-32B-Instruct

Important shape differences from anthropic_impl:
- System prompt goes in messages[0] role="system", not a top-level field.
- Tool spec shape:  [{type: "function", function: {name, description, parameters}}]
- Tool call response: message.tool_calls[] with function.arguments as a
  JSON *string* that must be parsed.
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
    TranslationResult,
    Target,
    FacetContext,
    FilterUnion,
)


_TOOL_MODELS: Dict[str, tuple[Type[BaseModel], Target]] = {
    "filter_structures": (StructureFilters, "structures"),
    "filter_polymers": (PolypeptideEntityFilters, "polymers"),
    "filter_ligands": (LigandFilters, "ligands"),
}

_CLARIFY_TOOL = "request_clarification"

_DEFAULT_BASE_URL = "https://openrouter.ai/api/v1"
_DEFAULT_MODEL = "anthropic/claude-haiku-4-5"
_DEFAULT_MAX_TOKENS = 1024


def _model_to_tool(name: str, description: str, model: Type[BaseModel]) -> dict:
    schema = model.model_json_schema(by_alias=False)
    schema.pop("title", None)
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": description,
            "parameters": schema,
        },
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
            "type": "function",
            "function": {
                "name": _CLARIFY_TOOL,
                "description": (
                    "Call this tool when the user's request is ambiguous and you cannot "
                    "confidently fill a filter object. Provide a short clarification question."
                ),
                "parameters": {
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
        },
    ]


def _build_system_prompt(target: Target, facets: FacetContext) -> str:
    """Same rules as the Anthropic impl -- keeps translations consistent across providers."""
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


class OpenAICompatNLTranslator:
    """Translator for any OpenAI-compatible chat completions endpoint with tools."""

    def __init__(
        self,
        model: Optional[str] = None,
        base_url: Optional[str] = None,
        api_key: Optional[str] = None,
        max_tokens: int = _DEFAULT_MAX_TOKENS,
    ) -> None:
        from openai import OpenAI  # lazy import so the backend boots without openai if unused

        resolved_key = api_key or os.environ.get("OPENAI_API_KEY") or os.environ.get("OPENROUTER_API_KEY")
        if not resolved_key:
            raise RuntimeError(
                "OpenAICompatNLTranslator requires OPENAI_API_KEY (or OPENROUTER_API_KEY) to be set."
            )

        resolved_base = base_url or os.environ.get("OPENAI_BASE_URL", _DEFAULT_BASE_URL)
        resolved_model = model or os.environ.get("OPENAI_MODEL", _DEFAULT_MODEL)

        # OpenRouter recommends (but does not require) these headers for usage
        # attribution. They're harmless against a self-hosted vLLM endpoint too.
        default_headers = {
            "HTTP-Referer": os.environ.get("OPENROUTER_REFERER", "http://localhost:8000"),
            "X-Title": os.environ.get("OPENROUTER_APP_TITLE", "tubxz-nl-query"),
        }

        self._client = OpenAI(
            api_key=resolved_key,
            base_url=resolved_base,
            default_headers=default_headers,
        )
        self._model = resolved_model
        self._max_tokens = max_tokens
        self._tools = _build_tools()

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

        response = self._client.chat.completions.create(
            model=self._model,
            max_tokens=self._max_tokens,
            tools=self._tools,
            tool_choice="auto",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
        )

        return self._interpret(response)

    def _interpret(self, response) -> TranslationResult:
        if not response.choices:
            return TranslationResult(clarification="Empty response from model.")

        message = response.choices[0].message
        tool_calls = getattr(message, "tool_calls", None) or []

        if not tool_calls:
            text = (message.content or "").strip()
            return TranslationResult(
                clarification=text or "I could not interpret that request. Please rephrase."
            )

        call = tool_calls[0]
        name = call.function.name
        raw_args = call.function.arguments or "{}"
        try:
            inp = json.loads(raw_args)
        except json.JSONDecodeError as e:
            return TranslationResult(
                clarification=f"Model returned unparseable tool arguments: {e}"
            )

        if not isinstance(inp, dict):
            return TranslationResult(
                clarification="Tool arguments were not a JSON object."
            )

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


def _summarize(filters: FilterUnion) -> str:
    d = filters.model_dump(exclude_none=True, exclude_defaults=True)
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
