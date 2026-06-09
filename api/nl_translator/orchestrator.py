"""Grounded assistant orchestrator.

A single bounded tool-use loop that replaces the three ad-hoc NL endpoints. The
model is given:
  - READ tools (api.nl_translator.retrieval) — look up real data; loop continues.
  - COMMIT tools — emit_answer / emit_actions / request_clarification /
    cannot_answer — exactly one ends the loop.

The hard rule baked into the system prompt: the model must NEVER state a residue
number, a count, a chain id, a structure id, or a binding site from memory. It
looks everything up with a read tool first. Empty results lead to an honest
`cannot_answer` rather than a hallucinated action.

The loop is provider-agnostic over any OpenAI-compatible chat-completions
endpoint with tool calling (OpenRouter by default).
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional, Tuple, Type

from pydantic import BaseModel, Field, ValidationError

from api.nl_translator.facets_loader import load_facet_context
from api.nl_translator.global_actions import ActionCard, EntityRef, QuerySpec
from api.nl_translator.interface import FacetContext
from api.nl_translator.retrieval import RETRIEVAL_TOOLS, run_retrieval_tool
from api.nl_translator.viewer_actions import VIEWER_ACTION_MODELS

_DEFAULT_BASE_URL = "https://openrouter.ai/api/v1"
_DEFAULT_MODEL = "anthropic/claude-haiku-4-5"
_DEFAULT_MAX_TOKENS = 1500
_MAX_STEPS = 5  # model<->DB round trips before we force a terminal

_VIEWER_MODELS_BY_NAME: Dict[str, Type[BaseModel]] = {m.__name__: m for m in VIEWER_ACTION_MODELS}


# ---------------------------------------------------------------------------
# Request / result envelopes
# ---------------------------------------------------------------------------

class PageContext(BaseModel):
    """What the frontend tells us about where the user is and what's loaded."""
    page: str = Field("landing", description="'landing' | 'catalogue' | 'structure'")
    # structure-page viewer state (mirrors the old ViewContext)
    rcsb_id: Optional[str] = None
    chain_ids: List[str] = Field(default_factory=list)
    ligand_keys: List[str] = Field(default_factory=list)
    view_mode: Optional[str] = None  # 'structure' | 'monomer'
    active_monomer_chain: Optional[str] = None
    active_family: Optional[str] = None
    # Labels of annotation tracks already loaded in the MSA, so the model can
    # avoid creating duplicates (or RemoveAnnotationTrack first).
    loaded_tracks: List[str] = Field(default_factory=list)
    # catalogue-page current filters (freeform; merged/overwritten as user intends)
    current_filters: Optional[Dict[str, Any]] = None


class ViewerActionCall(BaseModel):
    """One molstar action as the frontend dispatcher consumes it: {type, args}.
    Validated server-side against VIEWER_ACTION_MODELS before being returned."""
    type: str = Field(..., description="Action type, e.g. 'FocusBindingSite', 'AddAnnotationTrack'.")
    args: Dict[str, Any] = Field(default_factory=dict, description="Arguments for the action; shape depends on type.")


class SuggestedAction(BaseModel):
    """A viewer action OFFERED to the user as a clickable chip (not auto-applied).
    Use for informational answers where visualizing the data would help but the
    user didn't explicitly ask to change the view."""
    label: str = Field(..., description="Short chip text, e.g. 'Show these PTMs on the MSA'.")
    action: ViewerActionCall = Field(..., description="The molstar action to run when clicked.")


class TraceEntry(BaseModel):
    tool: str
    args: Dict[str, Any] = Field(default_factory=dict)
    ok: bool = True
    summary: str = ""


class AssistantResult(BaseModel):
    """Result the endpoint serializes for the frontend. kind ∈
    {'respond','clarify','cannot'}. A 'respond' may carry ANY combination of
    answer text, auto-applied viewer actions, offered suggested actions, routing
    cards, and entity pills — the frontend renders whatever is present."""
    kind: str  # 'respond' | 'clarify' | 'cannot'
    answer_markdown: Optional[str] = None
    data: Optional[Dict[str, Any]] = None
    summary: Optional[str] = None
    clarification: Optional[str] = None
    reason: Optional[str] = None
    cards: List[ActionCard] = Field(default_factory=list)
    queries: List[QuerySpec] = Field(default_factory=list)
    entities: List[EntityRef] = Field(default_factory=list)
    viewer_actions: List[ViewerActionCall] = Field(default_factory=list)
    suggested_actions: List[SuggestedAction] = Field(default_factory=list)
    dropped_actions: List[Dict[str, Any]] = Field(default_factory=list)
    trace: List[TraceEntry] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Commit (terminal) tool argument models
# ---------------------------------------------------------------------------

class RespondArgs(BaseModel):
    """Finish the turn. Fill WHATEVER is relevant — these combine freely:
      - answer_markdown: a grounded textual answer (data/overview/comparison).
      - viewer_actions: molstar/MSA actions to AUTO-APPLY now (use when the user
        explicitly asked to change the view, e.g. 'focus/show/highlight X').
      - suggested_actions: molstar/MSA actions OFFERED as clickable chips (use to
        let the user visualize data you just reported, without mutating their
        view unasked).
      - cards: routing cards to another page (with query_ref into queries[]).
      - queries: filter specs the cards reference by id.
      - entities: chains/residues/ligands to surface as interactive pills.
    Every fact in answer_markdown MUST come from a read-tool result you received."""
    answer_markdown: Optional[str] = Field(None, description="Grounded answer in light markdown.")
    summary: Optional[str] = Field(None, description="One-line human readback (used when there's no answer text).")
    data: Optional[Dict[str, Any]] = Field(None, description="Optional structured data backing the answer.")
    viewer_actions: Optional[List[ViewerActionCall]] = Field(None, description="Molstar actions to auto-apply (structure page).")
    suggested_actions: Optional[List[SuggestedAction]] = Field(None, description="Molstar actions offered as clickable chips (structure page).")
    cards: Optional[List[ActionCard]] = Field(None, description="Routing cards.")
    queries: Optional[List[QuerySpec]] = Field(None, description="Filter specs referenced by open_catalogue cards via query_ref.")
    entities: Optional[List[EntityRef]] = Field(None, description="Entities to surface as pills.")


class ClarifyArgs(BaseModel):
    question: str = Field(..., description="One-sentence clarification question.")


class CannotAnswerArgs(BaseModel):
    reason: str = Field(..., description="Brief honest reason we can't answer (e.g. no data for that organism).")


_COMMIT_TOOLS: Dict[str, Type[BaseModel]] = {
    "respond": RespondArgs,
    "request_clarification": ClarifyArgs,
    "cannot_answer": CannotAnswerArgs,
}

_COMMIT_DESCRIPTIONS: Dict[str, str] = {
    "respond": "Finish the turn. Combine answer_markdown + viewer_actions + suggested_actions + cards as appropriate. THE primary terminal.",
    "request_clarification": "Ask one clarifying question instead of acting, when intent is genuinely ambiguous.",
    "cannot_answer": "Honestly decline when the data isn't there or the request is out of scope. Prefer this over guessing.",
}


# ---------------------------------------------------------------------------
# Tool-spec construction
# ---------------------------------------------------------------------------

def _model_to_tool(name: str, description: str, model: Type[BaseModel]) -> dict:
    schema = model.model_json_schema(by_alias=False)
    schema.pop("title", None)
    return {"type": "function", "function": {"name": name, "description": description, "parameters": schema}}


def _build_tools() -> List[dict]:
    tools: List[dict] = []
    for t in RETRIEVAL_TOOLS:
        tools.append(_model_to_tool(t.name, t.description, t.args_model))
    for name, model in _COMMIT_TOOLS.items():
        tools.append(_model_to_tool(name, _COMMIT_DESCRIPTIONS[name], model))
    return tools


# ---------------------------------------------------------------------------
# System prompt
# ---------------------------------------------------------------------------

def _organism_table(facets: FacetContext) -> str:
    lines = []
    for o in facets.common_source_organisms:
        rep = o.get("rep_pdbs") or []
        tail = f" [e.g. {', '.join(rep)}]" if rep else " [none indexed]"
        lines.append(f"  - {o.get('tax_id')} -> {o.get('name')}{tail}")
    return "\n".join(lines) if lines else "  (none)"


_VIEWER_ACTION_VOCAB = """\
VIEWER ACTION VOCABULARY (fill viewer_actions as a list of {type, args}; structure page only):
- FocusChain {auth_asym_id}
- FocusResidue {auth_asym_id, auth_seq_id}
- FocusResidueRange {auth_asym_id, start, end}
- HighlightResidueRange {auth_asym_id, start, end}
- HighlightChain {auth_asym_id}
- ClearFocus {} / ClearHighlight {}
- SetChainVisibility {auth_asym_id, visible} / IsolateChain {auth_asym_id, keep_ligands}
- FocusBindingSite {chemical_id, auth_asym_id?}  -- only for a ligand actually bound on that chain (check get_structure_chains first)
- AddAnnotationTrack {label, color (#RRGGBB), spec}  -- expert mode only. spec.family is REQUIRED. spec is one of:
    {kind:"binding_contacts", family, chemical_ids:[...]}            <- note PLURAL chemical_ids, always a list
    {kind:"modifications", family, modification_types?:[...], species_tax_ids?:[...], position_range?:[min,max], positions?:[...]}
    {kind:"variants", family, sources?:[...], species_tax_ids?:[...], position_range?:[min,max], wild_type_aas?:[...], observed_aas?:[...]}
- RemoveAnnotationTrack {label_match}
- AlignChain {organism_id}  -- expert mode only; adds another organism's chain to the alignment
"""


_CATALOGUE_CARDS = """\
CATALOGUE CARDS (open_catalogue) — encode the FULL filter set; never drop a field:
- An open_catalogue card MUST reference a query. Put the filters in queries[] and set the card's query_ref to that query's id.
- The query's filters_structures MUST contain EVERY filter you used to compute the count. If you called find_structures(has_ligand_ids=["TA1"], source_organism_ids=[9606]), the query MUST set BOTH has_ligand_ids AND source_organism_ids — dropping the ligand is a bug.
- StructureFilters fields you can set: search, rcsb_ids, has_ligand_ids, has_polymer_family, has_isotype, source_organism_ids, exp_method, resolution_min, resolution_max, year_min, year_max, has_variants.
- Worked example — user: "show me human structures with taxol":
    queries = [{"id":"q1","target":"structures","filters_structures":{"has_ligand_ids":["TA1"],"source_organism_ids":[9606]}}]
    cards   = [{"action":"open_catalogue","query_ref":"q1","label":"Browse human taxol-bound structures",
                "description":"Catalogue filtered to human (9606) structures containing TA1."}]
"""


def _build_system_prompt(ctx: PageContext, facets: FacetContext) -> str:
    page = ctx.page or "landing"
    grounding = f"""You are the assistant for tube.xyz, a tubulin structural-biology database. You answer questions and drive the UI by calling tools.

ABSOLUTE GROUNDING RULE — this is the most important instruction:
- NEVER state a residue number, a count, a chain id, a PDB id, a binding-site position, or "which organism/family has X" from your own memory. You WILL be wrong.
- To know any such fact you MUST call a READ tool and use its result. If you haven't looked it up, you don't know it.
- If a read tool returns nothing (found:false / count 0 / empty), say so honestly via cannot_answer or in your answer. Do NOT invent a fallback fact.

HOW THE LOOP WORKS:
- Call READ tools to gather the real data you need (you may call several, across multiple turns; you'll see each result before deciding the next step).
- When you have what you need, finish with EXACTLY ONE COMMIT tool: `respond`, `request_clarification`, or `cannot_answer`. ALWAYS finish by CALLING one — never reply in plain prose.
- `respond` combines freely: answer_markdown (text), viewer_actions (auto-apply to the viewer/MSA now), suggested_actions (offer as clickable chips), cards (route to another page), queries (filters the cards reference). Fill ALL that apply in ONE respond call.
- Be economical: don't call read tools whose result you won't use.

VISUALIZE, DON'T JUST DESCRIBE (this is the point of the app):
- This site exists to SHOW data on a 3D structure and an MSA. A wall of text is a weak answer here.
- Whenever your answer concerns residues / positions / sites / modifications / variants / a ligand's binding site, and the viewer can display them, you MUST pair the text with a way to SEE it:
  - If the user explicitly asked to change the view ("focus / show / highlight / color X") → put the action in viewer_actions (auto-apply).
  - If you're answering a data/overview question ("how many PTMs...", "what variants...") → offer the visualization in suggested_actions (chips the user clicks), so you don't mutate their view unasked.
- The grounded data you fetched (e.g. count_modifications positions) maps directly to an AddAnnotationTrack with the SAME filter spec — that paints the MSA and the 3D chain. Prefer that over listing numbers alone.

READ TOOLS: {", ".join(t.name for t in RETRIEVAL_TOOLS)}.
- get_structure_chains: which chain is which family and which ligands actually contact it — call this before focusing/aligning a ligand.
- get_binding_site / get_binding_contacts: real binding residues (master positions). Never recall these.
- count_modifications / count_variants: real counts for 'how many ...' questions.
- resolve_structure: turn organism+family[+ligand] into a real (rcsb_id, chain) instead of guessing an id.
- find_structures: catalogue counts + samples. get_facets: valid family/ligand/organism vocabulary.

LIGAND NAMING: trivial names != PDB chem ids. taxol/paclitaxel=TA1, colchicine=LOC, vinblastine=VLB, GTP=GTP, GDP=GDP, maytansine=MYT. Use chem ids in tool args.
DOMAIN: resolution in Angstroms, LOWER=better. Ranges inclusive. Tubulin heterodimer: alpha often chain A, beta often chain B — but VERIFY with get_structure_chains, don't assume.

KNOWN ORGANISMS (tax id -> name):
{_organism_table(facets)}

VALID FAMILIES: {facets.tubulin_families}
"""

    if page == "structure":
        view = f"""
CURRENT VIEWER STATE:
- rcsb_id: {ctx.rcsb_id}
- loaded chains: {", ".join(ctx.chain_ids) or "(none)"}
- loaded ligand keys: {", ".join(ctx.ligand_keys) or "(none)"}
- view_mode: {ctx.view_mode}   active_monomer_chain: {ctx.active_monomer_chain}   active_family: {ctx.active_family}
- annotation tracks already loaded: {", ".join(ctx.loaded_tracks) or "(none)"}

PAGE GUIDANCE (structure page):
- Pair text with visualization (see VISUALIZE above). "how many PTMs in human alpha tubulin" → answer_markdown with the counts AND a suggested_action AddAnnotationTrack(modifications, family=tubulin_alpha, species_tax_ids=[9606]) labelled e.g. "Show these PTMs on the MSA".
- "focus / show the <ligand> binding site": call get_structure_chains. If the ligand IS bound on a loaded chain, viewer_actions=[FocusBindingSite(chemical_id, auth_asym_id=<that chain>)]. If it is NOT bound on any loaded chain, DO NOT decline — call get_binding_site(chemical_id, family) and emit viewer_actions=[AddAnnotationTrack(binding_contacts, family=<the family that binds it>, chemical_ids=[chemical_id])] so the canonical site is shown, and say in the text that the ligand isn't bound here so the conserved site is projected from {{n}} structures.
- AddAnnotationTrack / AlignChain require view_mode == 'monomer'. If view_mode is 'structure' (easy mode), do NOT emit them — instead offer a card to open expert mode, or answer in text + a suggested card.
- Don't recreate a track that's already in "annotation tracks already loaded"; if the user wants to replace it, RemoveAnnotationTrack(label_match) first.
- For questions about OTHER structures, respond with cards (+ queries), no viewer_actions.

{_VIEWER_ACTION_VOCAB}
(suggested_actions use the SAME {{type, args}} vocabulary, wrapped as {{label, action}}.)

{_CATALOGUE_CARDS}"""
    elif page == "catalogue":
        cf = json.dumps(ctx.current_filters) if ctx.current_filters else "(none)"
        view = f"""
PAGE GUIDANCE (catalogue page). Current filters: {cf}
- For "show/filter/find structures with X", respond with ONE open_catalogue card that references a query (see CATALOGUE CARDS below).
- For "how many ..." answer in answer_markdown (use find_structures / count_* first), optionally + a card.

{_CATALOGUE_CARDS}"""
    else:  # landing
        view = f"""
PAGE GUIDANCE (landing page).
- Route the user with cards: open_catalogue (browse — see CATALOGUE CARDS), open_structure / open_expert (a real structure — use resolve_structure or organism+family selectors, NEVER guess rcsb_id), inspect_ligand, view_variants.
- For data/overview/comparison questions, answer in answer_markdown (look up the real numbers first), optionally attaching cards.
- Rank cards most-specific first; never emit two cards with the same intent.

{_CATALOGUE_CARDS}"""

    return grounding + view


# ---------------------------------------------------------------------------
# Client
# ---------------------------------------------------------------------------

def _make_client() -> Tuple[Any, str, int]:
    from openai import OpenAI

    key = os.environ.get("OPENAI_API_KEY") or os.environ.get("OPENROUTER_API_KEY")
    if not key:
        raise RuntimeError("Assistant orchestrator requires OPENAI_API_KEY (or OPENROUTER_API_KEY).")
    base = os.environ.get("OPENAI_BASE_URL") or _DEFAULT_BASE_URL
    model = os.environ.get("OPENAI_ASSISTANT_MODEL") or os.environ.get("OPENAI_MODEL") or _DEFAULT_MODEL
    headers = {
        "HTTP-Referer": os.environ.get("OPENROUTER_REFERER") or "http://localhost:8000",
        "X-Title": os.environ.get("OPENROUTER_APP_TITLE") or "tubxz-assistant",
    }
    client = OpenAI(api_key=key, base_url=base, default_headers=headers)
    return client, model, _DEFAULT_MAX_TOKENS


# ---------------------------------------------------------------------------
# Action validation (Phase 0: validate type+args; full chain-map grounding is Phase 1)
# ---------------------------------------------------------------------------

def _normalize_action_args(call: ViewerActionCall) -> None:
    """Repair common, unambiguous shape mistakes in-place before validation.
    Kept deliberately narrow — only well-known singular/plural slips, never
    semantic guesses."""
    if call.type == "AddAnnotationTrack":
        spec = call.args.get("spec")
        if isinstance(spec, dict) and spec.get("kind") == "binding_contacts":
            if "chemical_ids" not in spec and "chemical_id" in spec:
                spec["chemical_ids"] = [spec.pop("chemical_id")]


def _validate_viewer_actions(
    raw: List[ViewerActionCall],
) -> Tuple[List[ViewerActionCall], List[Dict[str, Any]]]:
    kept: List[ViewerActionCall] = []
    dropped: List[Dict[str, Any]] = []
    for call in raw:
        model_cls = _VIEWER_MODELS_BY_NAME.get(call.type)
        if model_cls is None:
            dropped.append({"type": call.type, "reason": "unknown action type"})
            continue
        _normalize_action_args(call)
        try:
            model_cls.model_validate(call.args)
        except ValidationError as e:
            dropped.append({"type": call.type, "reason": f"invalid args: {e.errors(include_url=False)[:1]}"})
            continue
        kept.append(call)
    return kept, dropped


# ---------------------------------------------------------------------------
# The loop
# ---------------------------------------------------------------------------

def _result_summary(obj: Any) -> str:
    try:
        s = json.dumps(obj)
    except (TypeError, ValueError):
        s = str(obj)
    return s if len(s) <= 200 else s[:197] + "..."


def run_assistant(text: str, page_context: Optional[Dict[str, Any]] = None) -> AssistantResult:
    ctx = PageContext.model_validate(page_context or {})
    facets = load_facet_context()
    system_prompt = _build_system_prompt(ctx, facets)
    client, model, max_tokens = _make_client()
    tools = _build_tools()

    messages: List[Dict[str, Any]] = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": text.strip()},
    ]
    trace: List[TraceEntry] = []

    for step in range(_MAX_STEPS):
        force_terminal = step == _MAX_STEPS - 1
        kwargs: Dict[str, Any] = dict(model=model, max_tokens=max_tokens, tools=tools, messages=messages)
        # On the last allowed step, require a tool call so we don't end empty.
        kwargs["tool_choice"] = "required" if force_terminal else "auto"
        if force_terminal:
            messages.append({
                "role": "system",
                "content": "You have used all lookups. Finish now with `respond` (include any visualization in viewer_actions/suggested_actions) or `cannot_answer`, based only on results you already have.",
            })

        response = client.chat.completions.create(**kwargs)
        message = response.choices[0].message
        tool_calls = getattr(message, "tool_calls", None) or []

        if not tool_calls:
            # The model replied in prose instead of calling a commit tool. If it
            # had already gathered data, treat that prose as the answer (it IS
            # answering); otherwise it's a genuine clarification/refusal.
            content = (message.content or "").strip()
            if trace and content:
                return AssistantResult(kind="respond", answer_markdown=content, trace=trace)
            return AssistantResult(kind="clarify", clarification=content or "Please rephrase that.", trace=trace)

        # Record the assistant turn (with its tool calls) for protocol correctness.
        messages.append({
            "role": "assistant",
            "content": message.content or "",
            "tool_calls": [
                {"id": c.id, "type": "function",
                 "function": {"name": c.function.name, "arguments": c.function.arguments or "{}"}}
                for c in tool_calls
            ],
        })

        # If any commit tool is present, that ends the loop — build the result.
        commit = next((c for c in tool_calls if c.function.name in _COMMIT_TOOLS), None)
        if commit is not None:
            return _build_terminal_result(commit, trace)

        # Otherwise: execute all retrieval calls, append results, continue.
        for c in tool_calls:
            name = c.function.name
            try:
                raw_args = json.loads(c.function.arguments or "{}")
            except json.JSONDecodeError:
                raw_args = {}
            try:
                result = run_retrieval_tool(name, raw_args)
                trace.append(TraceEntry(tool=name, args=raw_args, ok=True, summary=_result_summary(result)))
            except Exception as e:  # KeyError (unknown tool) or ValidationError or DB error
                result = {"error": f"{type(e).__name__}: {str(e)[:200]}"}
                trace.append(TraceEntry(tool=name, args=raw_args, ok=False, summary=_result_summary(result)))
            messages.append({
                "role": "tool",
                "tool_call_id": c.id,
                "content": json.dumps(result),
            })

    # Cap reached without a terminal (shouldn't happen with tool_choice=required).
    return AssistantResult(kind="cannot", reason="Ran out of steps without a conclusion.", trace=trace)


def _build_terminal_result(commit, trace: List[TraceEntry]) -> AssistantResult:
    name = commit.function.name
    try:
        raw = json.loads(commit.function.arguments or "{}")
    except json.JSONDecodeError:
        return AssistantResult(kind="cannot", reason="Model produced unparseable commit arguments.", trace=trace)

    if name == "request_clarification":
        return AssistantResult(kind="clarify", clarification=raw.get("question", "Please clarify."), trace=trace)
    if name == "cannot_answer":
        return AssistantResult(kind="cannot", reason=raw.get("reason", "I can't answer that."), trace=trace)

    # respond — the unified terminal
    try:
        args = RespondArgs.model_validate(raw)
    except ValidationError as e:
        return AssistantResult(kind="cannot", reason=f"response failed validation: {e.errors(include_url=False)[:1]}", trace=trace)

    kept, dropped = _validate_viewer_actions(args.viewer_actions or [])
    # Suggested actions: validate the embedded action; drop the chip if invalid.
    suggested: List[SuggestedAction] = []
    for s in (args.suggested_actions or []):
        sk, sd = _validate_viewer_actions([s.action])
        if sk:
            suggested.append(SuggestedAction(label=s.label, action=sk[0]))
        else:
            dropped.extend([{**d, "suggested": True, "label": s.label} for d in sd])

    proposed_auto = len(args.viewer_actions or [])
    has_text = bool(args.answer_markdown)
    has_cards = bool(args.cards)
    # Honesty: if the model ONLY tried to auto-apply viewer actions and they all
    # failed (no text, no cards, no surviving chips), surface the failure rather
    # than a silent no-op.
    if proposed_auto > 0 and not kept and not has_text and not has_cards and not suggested:
        reasons = "; ".join(f"{d['type']}: {d['reason']}" for d in dropped)
        return AssistantResult(
            kind="cannot",
            reason=f"Could not apply the requested viewer action(s): {reasons}",
            dropped_actions=dropped,
            trace=trace,
        )

    return AssistantResult(
        kind="respond",
        answer_markdown=args.answer_markdown,
        summary=args.summary,
        data=args.data,
        viewer_actions=kept,
        suggested_actions=suggested,
        cards=args.cards or [],
        queries=args.queries or [],
        entities=args.entities or [],
        dropped_actions=dropped,
        trace=trace,
    )


__all__ = ["run_assistant", "AssistantResult", "PageContext", "ViewerActionCall", "SuggestedAction"]
