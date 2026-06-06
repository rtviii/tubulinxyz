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
    ViewContext,
    ViewerTranslationResult,
    GlobalTranslationResult,
)
from api.nl_translator.viewer_actions import (
    VIEWER_ACTION_MODELS,
    VIEWER_ACTION_DESCRIPTIONS,
    CLARIFY_ACTION_NAME,
    MENTION_ENTITIES_NAME,
    EMIT_NAV_CARD_NAME,
)
from api.nl_translator.global_actions import (
    GlobalNLResponse,
    EMIT_GLOBAL_TOOL_NAME,
    CLARIFY_GLOBAL_TOOL_NAME,
    MAX_CARDS,
    MAX_QUERIES,
    MAX_BLURB_CHARS,
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

        # Treat empty strings the same as unset. docker-compose's `${VAR:-}`
        # passthrough exports the variable as "" when it isn't in .env, which
        # `os.environ.get(KEY, DEFAULT)` would return as "" instead of falling
        # back to DEFAULT. Use a chained `or` so empty strings hit the default.
        resolved_base = base_url or os.environ.get("OPENAI_BASE_URL") or _DEFAULT_BASE_URL
        resolved_model = model or os.environ.get("OPENAI_MODEL") or _DEFAULT_MODEL

        # OpenRouter recommends (but does not require) these headers for usage
        # attribution. They're harmless against a self-hosted vLLM endpoint too.
        default_headers = {
            "HTTP-Referer": os.environ.get("OPENROUTER_REFERER") or "http://localhost:8000",
            "X-Title": os.environ.get("OPENROUTER_APP_TITLE") or "tubxz-nl-query",
        }

        self._client = OpenAI(
            api_key=resolved_key,
            base_url=resolved_base,
            default_headers=default_headers,
        )
        self._model = resolved_model
        self._max_tokens = max_tokens
        self._tools = _build_tools()
        self._viewer_tools = _build_viewer_tools()
        self._viewer_models_by_name: Dict[str, Type[BaseModel]] = {
            m.__name__: m for m in VIEWER_ACTION_MODELS
        }
        self._global_tools = _build_global_tools()

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

    def translate_global(
        self,
        text: str,
        facets: FacetContext,
    ) -> GlobalTranslationResult:
        """Front-page entry: emit `GlobalNLResponse` (blurb + queries + cards)."""
        system_prompt = _build_global_system_prompt(facets)

        response = self._client.chat.completions.create(
            model=self._model,
            max_tokens=self._max_tokens,
            tools=self._global_tools,
            tool_choice="auto",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": text.strip()},
            ],
        )

        return _interpret_global_openai(response)

    def translate_viewer(
        self,
        text: str,
        view_context: ViewContext,
    ) -> ViewerTranslationResult:
        # Facets are cached; cheap to load on every viewer call. Needed so the
        # universal preamble can render the KNOWN ORGANISMS grounding table.
        from api.nl_translator.facets_loader import load_facet_context
        facets = load_facet_context()
        system_prompt = _build_viewer_system_prompt(view_context, facets)

        response = self._client.chat.completions.create(
            model=self._model,
            max_tokens=self._max_tokens,
            tools=self._viewer_tools,
            tool_choice="auto",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": text.strip()},
            ],
        )

        return _interpret_viewer_openai(response, self._viewer_models_by_name)

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


def _build_viewer_tools() -> List[dict]:
    """Render VIEWER_ACTION_MODELS as OpenAI tool specs."""
    tools: List[dict] = []
    for model in VIEWER_ACTION_MODELS:
        name = model.__name__
        description = VIEWER_ACTION_DESCRIPTIONS.get(name, (model.__doc__ or "").strip())
        tools.append(_model_to_tool(name, description, model))
    return tools


def _build_viewer_system_prompt(view_context: ViewContext, facets: FacetContext) -> str:
    chains = ", ".join(view_context.chain_ids) or "(none loaded)"
    ligands = ", ".join(view_context.ligand_keys) or "(none)"
    return _build_universal_preamble(facets) + f"""You are a controller for a 3D molecular viewer (molstar) displaying a single PDB structure.

Each user turn, emit ONE OR MORE tool calls that together fulfill the request. The frontend executes them in order against the current viewer.

ACTIONS AVAILABLE (tool names): FocusChain, FocusResidue, FocusResidueRange, ClearFocus, SetChainVisibility, IsolateChain, HighlightChain, HighlightResidueRange, ClearHighlight, AlignChain, AddAnnotationTrack, RemoveAnnotationTrack, FocusBindingSite, MentionEntities, EmitNavigationCard, RequestClarification.

CURRENT VIEWER STATE:
- rcsb_id: {view_context.rcsb_id}
- loaded chains (auth_asym_id): {chains}
- loaded ligand keys: {ligands}
- view_mode: {view_context.view_mode}
- active_monomer_chain: {view_context.active_monomer_chain}
- active_family: {view_context.active_family}

VIEWER-SPECIFIC RULES:
- auth_asym_id arguments MUST be one of the loaded chain ids above. If the user names a chain that is not loaded, emit exactly one `RequestClarification` tool call instead, naming the available chains.
- Treat "residue 240" / "position 240" as auth_seq_id=240.
- Prefer the minimum number of actions that accomplish the intent. Don't add ClearFocus / ClearHighlight unless the user implies reset.
- If the user asks to "hide everything except X" / "show only X" / "isolate X", use IsolateChain (not a sequence of SetChainVisibility).
- Do NOT use EmitNavigationCard for in-page operations (focus, highlight, isolate) — those are viewer actions.

SURFACING ENTITIES (MentionEntities tool):
- AFTER your action tool calls, call MentionEntities ONCE with the chains / residues / ranges / ligands the user should be able to hover and click in the side panel.
- Use this whenever you focused or highlighted something specific (residue, range, chain, ligand). Skip it for pure visibility toggles or camera resets.
- Cap at 6 entities. Keep them genuinely useful — do not duplicate every action argument verbatim.
- Entity kinds: 'chain' (auth_asym_id), 'residue_range' (auth_asym_id + start + end), 'ligand' (chemical_id, optional auth_asym_id + auth_seq_id).
- Example: user asks "highlight the M-loop on chain B" — you emit HighlightResidueRange(B, 217, 230) then MentionEntities(entities=[{{kind:"residue_range", auth_asym_id:"B", start:217, end:230}}]).

NAVIGATION INTENT (EmitNavigationCard tool):
- If the user's question is about OTHER structures, browsing the catalogue, or switching to a different PDB entry, call EmitNavigationCard with ONE ActionCard.
- Signs of navigation intent: "what other structures...", "find me X", "show structures with Y", "switch to Z", references to a different PDB id than the one currently loaded.
- The card uses the same vocabulary as the front-page: open_catalogue, open_structure, open_expert, inspect_ligand, view_variants. Always fill the `description` field (~60-100 chars).
- Catalogue cards from this endpoint CANNOT reference queries (set `query_ref:null`). Encode filter intent directly on the card payload:
  - For a ligand-filtered catalogue, set `focus_ligands: ["TA1"]`.
  - For a family-filtered catalogue, set `family: "tubulin_beta"`.
  - For a specific PDB id, set `rcsb_id`.
- Example: on 1JFF, user asks "other structures with taxol" — emit EmitNavigationCard(card={{action:"open_catalogue", label:"Browse structures with taxol", description:"Catalogue filtered to all PDB entries containing taxol.", query_ref:null, focus_ligands:["TA1"]}}).

ADDING A SEQUENCE TO THE CURRENT ALIGNMENT (AlignChain tool) — EXPERT MODE ONLY (view_mode == "monomer"):
- When the user asks to ADD / LOAD / INCLUDE / "also show" another organism's sequence or structure in the CURRENT alignment (e.g. "add a bovine sequence", "also load human alpha", "compare with yeast tubulin", "add cow"), call AlignChain. This adds a new aligned row to the MSA + a ghost overlay in 3D and KEEPS everything already loaded.
- Set `organism_id` to the NCBI tax id (see ORGANISM TAX IDS). Leave `family`, `rcsb_id`, `auth_asym_id` null — the backend resolves a real structure+chain of the active family (active_family={view_context.active_family}).
- Do NOT use EmitNavigationCard for this — that navigates away and loses the already-loaded sequences. Do NOT emit RequestClarification — the intent is clear: add to the current view.
- Emit ONLY the AlignChain tool call for this intent (optionally one MentionEntities after). Do not emit clarification or a nav card alongside it, or the action will be dropped.
- If view_mode is NOT "monomer", treat "add a sequence" as navigation instead: use EmitNavigationCard(open_expert).
- Example (expert mode, active_family=tubulin_alpha): user "add a bovine sequence" -> AlignChain(organism_id=9913). user "also show human" -> AlignChain(organism_id=9606).

ANNOTATION TRACKS (AddAnnotationTrack / RemoveAnnotationTrack) — EXPERT MODE ONLY (view_mode == "monomer"):
- A track is a chain-independent MSA aux row that paints master columns matching a typed FilterSpec, and colors the corresponding residues on the primary 3D chain. It is the right tool for "show me X residues" overlay questions — variants, PTMs, ligand binding sites — scoped by organism, family, or position set.
- The `spec.kind` discriminator picks the FilterSpec type:
  - 'variants' (VariantSpec): wild_type_aas, observed_aas, sources (['structural','literature']), species_tax_ids (NCBI tax id list), positions, position_range, co_occurs_with_mod_type, phenotype_contains.
  - 'modifications' (ModificationSpec): modification_types (e.g. ['phosphorylation','acetylation','polyglutamylation']), species_tax_ids, positions, position_range, co_occurs_with_variant, phenotype_contains.
  - 'binding_contacts' (BindingContactSpec): chemical_ids (e.g. ['GTP','TA1','LOC']), positions.
- ALWAYS set `spec.family` (e.g. 'tubulin_alpha' / 'tubulin_beta'). For the active chain use active_family from view_context. The track will paint only master columns of that family.
- ALWAYS pick a hex `color`. For COMPARATIVE queries emit MULTIPLE AddAnnotationTrack calls (one per organism / one per ligand) with visually distinct colors.
- `label` is what the user sees in the MSA aux panel — keep short (~30 chars). Encode any disambiguating dimension in the label (e.g. "GTP site · human", "GTP site · Toxoplasma"). Comparative tracks SHOULD share a common label prefix and differ only in the disambiguating tail.
- If view_mode is NOT "monomer", do NOT emit AddAnnotationTrack — emit EmitNavigationCard(open_expert) so the user enters expert mode first.
- Do NOT emit RequestClarification when the user's intent is clearly an annotation overlay; emit the track(s). Do NOT emit AddAnnotationTrack alongside EmitNavigationCard — they are mutually exclusive (the nav card preempts).

EXAMPLES (expert mode):
- "Show PTMs in human alpha tubulin associated with fibrosis":
  AddAnnotationTrack(label="PTMs · human · fibrosis", color="#E74C3C", spec={{kind:"modifications", family:"tubulin_alpha", species_tax_ids:[9606], phenotype_contains:["fibrosis"]}})
- "Where does GTP bind on alpha tubulin?" (when GTP is in ligand_keys — use the 3D focus verb):
  FocusBindingSite(chemical_id="GTP")
- "Where does GTP bind on alpha tubulin?" (when GTP is NOT in ligand_keys — use the family-wide track instead):
  AddAnnotationTrack(label="GTP site", color="#1ABC9C", spec={{kind:"binding_contacts", family:"tubulin_alpha", chemical_ids:["GTP"]}})
- "Compare the GTP site in human vs Toxoplasma" — AlignChain to add the other organism, then FocusBindingSite on the primary so the binding pocket is visible across both:
  AlignChain(organism_id=5811)
  FocusBindingSite(chemical_id="GTP")
- "Show structural variants in human beta tubulin between positions 200 and 280":
  AddAnnotationTrack(label="Variants · human · 200-280", color="#8E44AD", spec={{kind:"variants", family:"tubulin_beta", species_tax_ids:[9606], sources:["structural"], position_range:[200,280]}})
- "Hide the GTP site track" — RemoveAnnotationTrack(label_match="GTP site")

FOCUS BINDING SITE (FocusBindingSite tool):
- Use when the ligand is ACTUALLY BOUND on the loaded structure (it appears in ligand_keys above). This anchors the 3D camera + draws the contact representation + jumps the MSA to the contact span.
- Pair with AlignChain in comparative queries: align the other organism FIRST, then FocusBindingSite on the active chain — the binding site is then visible on the primary while the alignment row carries the comparison.
- Defaults auth_asym_id to active_monomer_chain; only set it explicitly if the user names a different chain that's also loaded.
- If the ligand is NOT in ligand_keys, do NOT call FocusBindingSite — the per-chain contact data doesn't exist. Fall back to AddAnnotationTrack(binding_contacts) for a family-wide MSA overlay, or use the HighlightResidueRange + nav-card pattern.

HANDLING "WHERE WOULD X BIND" / "SHOW ME X'S SITE" WHEN X IS NOT LOADED:
- The loaded structure does NOT have ligand X bound. The user wants you to project the binding site from biological knowledge onto the loaded chain.

- HARD RULE: For any "where" / "show" / "locate" question, you MUST visually answer with a highlight on the 3D viewer. Words alone are insufficient. The viewer must reflect the user's question.

- HARD RULE: The FIRST tool call you emit MUST be `HighlightResidueRange`. If your first tool call would be RequestClarification or EmitNavigationCard, STOP and emit HighlightResidueRange first. If you cannot identify residues to highlight even from biological knowledge, do not respond at all — re-read the question.

- REQUIRED tool call ORDER for "where would X bind" queries (emit ALL FOUR, in this exact order):
  1. **HighlightResidueRange** — mandatory. Pick the relevant binding residues from biological knowledge and call this tool with concrete numbers. For heterodimer structures, beta-tubulin is typically chain B (verify against view_context's loaded chains).
     Conserved binding residue references (approximate, well-conserved across species):
     - Taxol / paclitaxel (TA1): beta-tubulin M-loop, residues 217-230.
     - Colchicine (LOC): beta-tubulin T7 loop / H7-H8 region, residues 240-260.
     - GTP / GDP: beta-tubulin (E-site) or alpha-tubulin (N-site) nucleotide pocket, residues 1-30 plus loop T4 around 140-180.
     - Vinblastine / vincristine: beta1-alpha2 dimer interface region.
     - Maytansine (MYT): alpha-tubulin S6-H7 / T7 loop region.
  2. **MentionEntities** — list the residue_range you just highlighted so the user can hover/click it in the side panel: `entities=[{{kind:"residue_range", auth_asym_id:<chain>, start:<lo>, end:<hi>}}]`.
  3. **EmitNavigationCard** — one open_catalogue card with `focus_ligands:["<pdb_chem_id>"]` so the user can jump to structures where the ligand is actually bound. Always fill the `description` field.
  4. **RequestClarification** — brief honest note (≤120 chars): "X (CODE) is not bound here; binding region projected onto chain Y from biology. Card below jumps to bound structures."

- DO NOT emit only RequestClarification or only EmitNavigationCard for these queries. The user expects to SEE SOMETHING CHANGE in the 3D viewer; that's HighlightResidueRange.

- WORKED EXAMPLE — user asks "where would colchicine bind" on a 1JFF-style structure (chains A,B; ligand_keys=[TA1,GTP,GDP]). Your response is FOUR tool calls in this order:
  call 1: HighlightResidueRange(auth_asym_id="B", start=240, end=260)
  call 2: MentionEntities(entities=[{{kind:"residue_range", auth_asym_id:"B", start:240, end:260}}])
  call 3: EmitNavigationCard(card={{action:"open_catalogue", label:"Structures with colchicine bound", description:"Catalogue filtered to PDB entries containing colchicine (LOC).", query_ref:null, focus_ligands:["LOC"]}})
  call 4: RequestClarification(question="Colchicine (LOC) is not bound in 1JFF; the T7 loop binding region on chain B is highlighted from biology. Card jumps to structures where colchicine is actually bound.")
  Four tool calls. All of them. In that order. Do not skip call 1.
"""


class _OpenAICompatViewerMixin:
    """Kept separate only for readability; actual methods attached to the class below."""


def _interpret_viewer_openai(response, viewer_models_by_name: Dict[str, Type[BaseModel]]) -> ViewerTranslationResult:
    if not response.choices:
        return ViewerTranslationResult(clarification="Empty response from model.")

    message = response.choices[0].message
    tool_calls = getattr(message, "tool_calls", None) or []

    if not tool_calls:
        text = (message.content or "").strip()
        return ViewerTranslationResult(
            clarification=text or "I could not interpret that request. Please rephrase."
        )

    actions: List[BaseModel] = []
    entities: List[Any] = []
    nav_card: Optional[Any] = None
    clarification: Optional[str] = None

    for call in tool_calls:
        name = call.function.name
        raw_args = call.function.arguments or "{}"
        try:
            inp = json.loads(raw_args)
        except json.JSONDecodeError as e:
            return ViewerTranslationResult(
                clarification=f"Model returned unparseable tool arguments: {e}"
            )
        if not isinstance(inp, dict):
            return ViewerTranslationResult(
                clarification="Tool arguments were not a JSON object."
            )

        if name == CLARIFY_ACTION_NAME:
            # Accumulate; don't break. The LLM may emit RequestClarification
            # PLUS EmitNavigationCard so the user sees an explanation AND a
            # concrete action to take (e.g. "this structure lacks colchicine,
            # but here's a card to browse ones that have it").
            clarification = inp.get("question", "Please clarify.")
            continue

        if name == MENTION_ENTITIES_NAME:
            # Companion tool — doesn't enter actions[], surfaces side-panel pills.
            raw_entities = inp.get("entities", [])
            if isinstance(raw_entities, list):
                entities = raw_entities
            continue

        if name == EMIT_NAV_CARD_NAME:
            # Navigation intent — overrides any viewer actions. We accept the
            # last emitted nav_card if multiple slip through.
            raw_card = inp.get("card")
            if isinstance(raw_card, dict):
                nav_card = raw_card
            continue

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
        # Clarification + (optionally) a nav card. The frontend renders the
        # clarification message in the chat input popover and the card in the
        # side panel.
        return ViewerTranslationResult(
            clarification=clarification,
            nav_card=nav_card,
        )

    if nav_card is not None:
        # Nav card overrides actions — they may not make sense in context.
        return ViewerTranslationResult(
            nav_card=nav_card,
            summary="",
        )

    return ViewerTranslationResult(
        actions=actions,
        entities=entities,
        summary=_summarize_actions(actions),
    )


def _summarize_actions(actions: List[BaseModel]) -> str:
    if not actions:
        return "(no actions)"
    parts = []
    for a in actions:
        if type(a).__name__ == "AlignChain":
            parts.append("Adding a sequence to the alignment")
            continue
        kv = ", ".join(f"{k}={v}" for k, v in a.model_dump().items())
        parts.append(f"{type(a).__name__}({kv})" if kv else type(a).__name__)
    return " ; ".join(parts)


def _build_universal_preamble(facets: FacetContext) -> str:
    """Shared rules prepended to BOTH global and viewer system prompts.

    Single source of truth for anti-hallucination, ligand-naming, domain
    conventions, and the organism→PDB-ID grounding table.
    """
    organism_lines: List[str] = []
    for o in facets.common_source_organisms:
        rep = o.get("rep_pdbs") or []
        rep_str = (
            f" [PDB ids: {', '.join(rep)}]"
            if rep
            else " [no entries indexed — DO NOT GUESS, use catalogue + tax_id filter]"
        )
        organism_lines.append(f"  - {o.get('tax_id')} -> {o.get('name')}{rep_str}")
    organism_table = "\n".join(organism_lines) if organism_lines else "  (none loaded)"

    return f"""=== UNIVERSAL RULES (apply to all responses) ===

GROUNDING — DO NOT AUTHOR PDB IDS FROM MEMORY (highest-priority rule):
- You will get literal PDB ids wrong — either inventing a non-existent id, or citing a real id with the wrong organism/family. So for entity cards you generally must NOT write `rcsb_id` yourself.
- Instead EXPRESS INTENT and let the backend resolve a real structure from the database:
  - open_structure / open_expert / inspect_ligand: set `family` (e.g. "tubulin_alpha") and `primary_organism_id` (NCBI tax id, e.g. 9606=human). The backend picks the best real structure + chain of that family from that organism. Leave `rcsb_id` null.
  - "compare A vs B" (open_expert): additionally set `aligned_organism_ids` (e.g. [5811]=Toxoplasma). The backend resolves each to a real aligned structure. Leave `aligned` null.
  - inspect_ligand: set `chemical_id` (+ optional `family` / `primary_organism_id`). The backend finds a real structure where that ligand binds.
- The ONLY time you write a literal `rcsb_id` is when the id came from (a) the user's query text, or (b) view_context. Then set `rcsb_id` directly and leave the organism selectors null.
- The backend DROPS any card whose intent resolves to nothing, so expressing intent is always safe — never substitute a guess.

ORGANISM TAX IDS:
- Map species to NCBI tax id for the selectors: human=9606, mouse=10090, rat=10116, pig/Sus scrofa=9823, cow/Bos taurus=9913, Toxoplasma gondii=5811, yeast/S. cerevisiae=4932. Use your taxonomy knowledge for others.
- You may also emit an open_catalogue card with `source_organism_ids:[tax_id]` so the user can browse an organism in the catalogue.

BLURB CONSISTENCY:
- The blurb MUST ONLY mention organisms/structures that correspond to cards in your response. Don't name a species or PDB id the cards don't back.

KNOWN ORGANISMS (organisms we have indexed, with tax ids + a few example PDB ids — use the TAX ID in the *_organism_id(s) selectors):
{organism_table}

LIGAND NAMING:
- Trivial names and PDB chemical IDs differ. Always use the chemical ID (3-4 char code).
- Common: taxol/paclitaxel=TA1, colchicine=LOC, vinblastine=VLB, vincristine=VCR, GTP=GTP, GDP=GDP, maytansine=MYT.

DOMAIN CONVENTIONS:
- Resolution in Angstroms: LOWER = BETTER. "Better than 3 A" -> resolution_max=3.0.
- Range expressions ("between X and Y", "X-Y") are INCLUSIVE on both ends.
- Tubulin heterodimer chains: alpha typically chain A, beta typically chain B. Verify against view_context's loaded chains when present.

OUTPUT FORMAT:
- Tool calls only. Do not emit free-form prose alongside tool calls unless the tool itself accepts text as an argument.

=== END UNIVERSAL RULES ===

"""


def _build_global_tools() -> List[dict]:
    """Two tools: emit the global response, or ask for clarification."""
    return [
        _model_to_tool(
            EMIT_GLOBAL_TOOL_NAME,
            (
                "Emit the structured global response (blurb + queries + cards). "
                "Call this for any query you can act on. The response routes the user "
                "to the right page with state preloaded."
            ),
            GlobalNLResponse,
        ),
        {
            "type": "function",
            "function": {
                "name": CLARIFY_GLOBAL_TOOL_NAME,
                "description": (
                    "Call instead of emit_global_response only when the request is so "
                    "ambiguous no action card would be useful. Provide one short question."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "question": {
                            "type": "string",
                            "description": "One-sentence clarification question.",
                        }
                    },
                    "required": ["question"],
                },
            },
        },
    ]


_GLOBAL_FEW_SHOTS = """\
FEW-SHOT EXAMPLES (shape only — express organism/family INTENT via selectors; the backend fills real ids. Note descriptions never name a specific PDB id, since you don't know which structure will be resolved):

Example 1 — "where does taxol bind"
  blurb: "Taxol (PDB chem id TA1) binds the beta-tubulin taxane pocket; structures with TA1 below."
  queries: [{ id: "q1", target: "structures", filters_structures: { has_ligand_ids: ["TA1"] } }]
  cards:
    - { action: "open_catalogue", query_ref: "q1",
        label: "Browse structures with taxol",
        description: "Catalogue filtered to all PDB entries containing chem id TA1." }
    - { action: "inspect_ligand", chemical_id: "TA1", family: "tubulin_beta",
        label: "View taxol in its binding pocket",
        description: "Opens a real structure where TA1 binds beta-tubulin, ligand focused." }
    - { action: "open_expert", family: "tubulin_beta", focus_range: { start: 217, end: 230 },
        label: "Inspect the M-loop taxane pocket",
        description: "Expert mode on a beta-tubulin chain; residues 217-230 (M-loop) highlighted." }

Example 2 — "what kinds of PTMs are in tubulin"
  blurb: "Tubulin PTMs (acetylation, glutamylation, detyrosination) are annotated per-residue in expert mode."
  queries: []
  cards:
    - { action: "open_expert", family: "tubulin_alpha",
        label: "Open an alpha-tubulin chain with PTM annotations",
        description: "Expert mode on an alpha-tubulin chain; PTM annotations appear in the side panel." }
    - { action: "view_variants", family: "tubulin_alpha",
        label: "Browse alpha-tubulin sequence variants",
        description: "Catalogue filtered to structures with alpha-tubulin variant annotations." }

Example 3 — "differences between the GTP binding site in human and toxoplasma tubulin"
  (Express both organisms as tax-id selectors on ONE open_expert card. The backend resolves
   each to a real structure of the family; if an organism has no structures the backend drops
   that aligned side automatically. Never write rcsb_id / aligned here.)
  blurb: "GTP binds alpha-tubulin's N-terminal pocket (~res 140-180); compare human and Toxoplasma chains aligned in expert mode."
  queries: []
  cards:
    - { action: "open_expert", family: "tubulin_alpha",
        primary_organism_id: 9606, aligned_organism_ids: [5811],
        focus_range: { start: 140, end: 180 },
        label: "Compare GTP site: human vs Toxoplasma alpha-tubulin",
        description: "Loads a human and a Toxoplasma alpha-tubulin chain aligned in expert mode; residues 140-180 highlighted." }
    - { action: "view_variants", family: "tubulin_alpha", position_min: 140, position_max: 180,
        label: "Browse alpha-tubulin variants near GTP pocket",
        description: "Catalogue filtered to alpha-tubulin variants in residues 140-180." }
"""


def _build_global_system_prompt(facets: FacetContext) -> str:
    return _build_universal_preamble(facets) + f"""You are the front-page query router for tube.xyz, a tubulin structural biology database. The user types an open-ended question; you produce a single structured response that routes them to the right page with state preloaded.

YOUR OUTPUT: Exactly one tool call to `{EMIT_GLOBAL_TOOL_NAME}`. Only call `{CLARIFY_GLOBAL_TOOL_NAME}` if the request is so ambiguous that no card would be useful.

THE RESPONSE HAS THREE PARTS:
1. `blurb` (optional, <={MAX_BLURB_CHARS} chars). One sentence stating what the cards collectively show. No preamble ("Here are..."), state the substance. Leave empty if cards are self-evident.
2. `queries` (0-{MAX_QUERIES} entries). Filter specs the frontend will run via the existing list endpoints. Each has an `id` (q1, q2, ...) and exactly one of `filters_structures` / `filters_polymers` / `filters_ligands`. Used only by `open_catalogue` cards via `query_ref`.
3. `cards` (1-{MAX_CARDS} entries, ranked most-specific first). One click each, routing to a page.

CARD TYPES (for entity cards set `family` + organism selectors, NOT `rcsb_id` — see GROUNDING):
- open_catalogue: aggregate browse. Set `query_ref` -> queries[].id, or `source_organism_ids` directly. Use for "list", "browse", "find all".
- open_structure: one structure in easy mode. Set `family` + `primary_organism_id`. Optional `focus_ligands`.
- open_expert: chain-level / alignment / residue ranges. Set `family` + `primary_organism_id`; for comparisons add `aligned_organism_ids`. Optional `focus_range`. Use for "compare", "align", "residue X", "binding site Y".
- inspect_ligand: one chemical. Set `chemical_id` (+ optional `family` / `primary_organism_id`).
- view_variants: variant browse. Set `family`. Optional `position_min`/`max`, `variant_type`.
- clarify: ask a question. Set `question`. (Single-card responses only.)

INTENT HEURISTICS:
- "list" / "browse" / "find structures with X" -> open_catalogue
- Specific PDB id named by the user -> open_structure (default) or open_expert (if alignment/range), with `rcsb_id` set to that id (the one case you set rcsb_id directly).
- "compare" / "align" / "differences between" -> open_expert with `primary_organism_id` + `aligned_organism_ids` (one tax id per organism). The backend resolves each side to a real structure; if an organism has none, that side is dropped automatically.
- Specific residue numbers / "binding site" / "pocket" -> open_expert with `focus_range`
- Ligand by trivial name (e.g. "taxol", "paclitaxel" both = TA1) -> inspect_ligand and/or open_expert
- Mutation / variant / substitution language -> view_variants (+ open_expert if specific structure)

RANKING / BUDGET:
- 2-5 cards typical. Never more than 6. Never duplicate intent.
- Never put `open_expert` as card #1 unless the query contains "align", "compare", "residue", "position", "binding site", or "pocket".
- `clarify` is always a single-card response.
- For nucleotides/ions (GTP, GDP, ATP, ADP, MG, ZN, etc.), do NOT emit an `inspect_ligand` card — they are too generic. Fall through to other cards.

DESCRIPTION FIELD (REQUIRED on every non-clarify card):
- Every card MUST have a `description` of ~60-100 chars stating concretely "what you'll see when you click" — the organism, family, chain, range, ligand, or filter that the click reveals.
- Do NOT name a specific PDB id in the description: for entity cards you set selectors, not ids, so you don't know which structure the backend will pick. Describe it by organism + family + what's highlighted.
- Bad: "Detailed structural analysis." (vague)
- Good: "Expert mode on a human alpha-tubulin chain; residues 140-180 (GTP pocket) highlighted."
- The `label` is the headline; `description` is the precise consequence. Don't repeat the label.

FACET VALUES (use verbatim where applicable):
- tubulin_families: {facets.tubulin_families}
- isotypes: {facets.isotypes}
- exp_methods: {facets.exp_methods}
- top ligands (chemical_id -> chemical_name): {facets.top_ligands}
- year_range: {facets.year_range}, resolution_range: {facets.resolution_range}

{_GLOBAL_FEW_SHOTS}"""


def _interpret_global_openai(response) -> GlobalTranslationResult:
    if not response.choices:
        return GlobalTranslationResult(clarification="Empty response from model.")

    message = response.choices[0].message
    tool_calls = getattr(message, "tool_calls", None) or []

    if not tool_calls:
        text = (message.content or "").strip()
        return GlobalTranslationResult(
            clarification=text or "I could not interpret that request. Please rephrase."
        )

    call = tool_calls[0]
    name = call.function.name
    raw_args = call.function.arguments or "{}"
    try:
        inp = json.loads(raw_args)
    except json.JSONDecodeError as e:
        return GlobalTranslationResult(
            clarification=f"Model returned unparseable tool arguments: {e}"
        )
    if not isinstance(inp, dict):
        return GlobalTranslationResult(
            clarification="Tool arguments were not a JSON object."
        )

    if name == CLARIFY_GLOBAL_TOOL_NAME:
        return GlobalTranslationResult(
            clarification=inp.get("question", "Please clarify.")
        )

    if name != EMIT_GLOBAL_TOOL_NAME:
        return GlobalTranslationResult(
            clarification=f"Internal: model invoked unknown tool {name!r}."
        )

    # The LLM must not set `validation`; strip it if present.
    inp.pop("validation", None)

    try:
        resp = GlobalNLResponse.model_validate(inp)
    except ValidationError as e:
        return GlobalTranslationResult(
            clarification=(
                "The model's response did not validate. "
                f"Details: {e.errors(include_url=False)[:3]}"
            )
        )

    return GlobalTranslationResult(response=resp)


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
