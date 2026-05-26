"""Natural-language query translation route.

Takes a free-form user string and returns a structured filter object
(one of StructureFilters / PolypeptideEntityFilters / LigandFilters) OR a
clarification question. Does NOT execute any DB queries -- the frontend
applies the returned filters to the existing filter state and triggers
the normal list endpoint.
"""
from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from api.nl_translator import get_translator
from api.nl_translator.facets_loader import load_facet_context
from api.nl_translator.interface import ViewContext
from api.nl_translator.global_actions import GlobalNLResponse, ActionCard
from api.nl_translator.hydration import hydrate_response
from api.nl_translator.resolve import resolve_response, resolve_card, resolve_representative


router_nl_query = APIRouter()


class NLQueryRequest(BaseModel):
    text: str = Field(..., min_length=1, max_length=2000)
    target: Literal["structures", "polymers", "ligands"] = "structures"
    current_filters: Optional[Dict[str, Any]] = None


class NLQueryResponse(BaseModel):
    target: Optional[Literal["structures", "polymers", "ligands"]] = None
    filters: Optional[Dict[str, Any]] = None
    summary: str = ""
    clarification: Optional[str] = None


@router_nl_query.post(
    "/filters",
    response_model=NLQueryResponse,
    operation_id="nl_query_to_filters",
)
def nl_query_to_filters(req: NLQueryRequest) -> NLQueryResponse:
    """Translate a natural-language request into a filter object."""
    try:
        translator = get_translator()
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Translator unavailable: {e}")

    facets = load_facet_context()

    try:
        result = translator.translate(
            text=req.text,
            target=req.target,
            facets=facets,
            current_filters=req.current_filters,
        )
    except Exception as e:
        raise HTTPException(status_code=502, detail=f"Translator error: {e}")

    filters_dict = (
        result.filters.model_dump(exclude_none=True, exclude_defaults=True)
        if result.filters is not None
        else None
    )

    return NLQueryResponse(
        target=result.target,
        filters=filters_dict,
        summary=result.summary,
        clarification=result.clarification,
    )


# ---------------------------------------------------------------------------
# Viewer-action translation
# ---------------------------------------------------------------------------


class ViewContextBody(BaseModel):
    rcsb_id: Optional[str] = None
    chain_ids: List[str] = Field(default_factory=list)
    ligand_keys: List[str] = Field(default_factory=list)
    view_mode: Optional[str] = None
    active_monomer_chain: Optional[str] = None
    active_family: Optional[str] = None


class NLViewerRequest(BaseModel):
    text: str = Field(..., min_length=1, max_length=2000)
    view_context: ViewContextBody = Field(default_factory=ViewContextBody)


class ViewerActionEnvelope(BaseModel):
    """One action to execute, typed by `type` (the Pydantic class name)."""
    type: str
    args: Dict[str, Any] = Field(default_factory=dict)


class NLViewerResponse(BaseModel):
    kind: Literal["viewer_actions", "clarify", "nav_card"]
    actions: List[ViewerActionEnvelope] = Field(default_factory=list)
    entities: List[Dict[str, Any]] = Field(
        default_factory=list,
        description=(
            "EntityRef list surfaced via the MentionEntities tool. Frontend "
            "renders them as interactive pills in the side panel."
        ),
    )
    card: Optional[ActionCard] = Field(
        default=None,
        description=(
            "Set when kind='nav_card': a single ActionCard from the global "
            "vocabulary, rendered as a navigation chip in the side panel."
        ),
    )
    summary: str = ""
    clarification: Optional[str] = None


@router_nl_query.post(
    "/viewer",
    response_model=NLViewerResponse,
    operation_id="nl_query_to_viewer_actions",
)
def nl_query_to_viewer_actions(req: NLViewerRequest) -> NLViewerResponse:
    """Translate a natural-language request into an ordered list of viewer actions.

    The frontend is responsible for executing the actions against its
    MolstarInstance. This route only validates arguments and never touches
    the viewer or the database.
    """
    try:
        translator = get_translator()
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Translator unavailable: {e}")

    view_context = ViewContext(
        rcsb_id=req.view_context.rcsb_id,
        chain_ids=req.view_context.chain_ids,
        ligand_keys=req.view_context.ligand_keys,
        view_mode=req.view_context.view_mode,
        active_monomer_chain=req.view_context.active_monomer_chain,
        active_family=req.view_context.active_family,
    )

    try:
        result = translator.translate_viewer(text=req.text, view_context=view_context)
    except Exception as e:
        raise HTTPException(status_code=502, detail=f"Translator error: {e}")

    if result.clarification is not None:
        # Pass through nav_card if the LLM also emitted one alongside the
        # clarification. Best-effort validation; bad cards just get dropped.
        card: Optional[ActionCard] = None
        if result.nav_card is not None:
            try:
                card = ActionCard.model_validate(result.nav_card)
                # Resolve organism/family selectors to real ids (same as the
                # global flow); drop the nav card if it can't be resolved.
                if not resolve_card(card):
                    card = None
            except Exception:
                card = None
        return NLViewerResponse(
            kind="clarify",
            clarification=result.clarification,
            summary=result.summary,
            card=card,
        )

    if result.nav_card is not None:
        # Navigation intent — try to validate the card; if it doesn't fit the
        # ActionCard schema, downgrade to a clarification rather than 500.
        try:
            card = ActionCard.model_validate(result.nav_card)
        except Exception as e:
            return NLViewerResponse(
                kind="clarify",
                clarification=f"Could not parse navigation card: {e}",
            )
        # Resolve organism/family selectors to real ids; if the card's entity
        # can't be resolved, downgrade to a clarification rather than ship a
        # broken nav card.
        if not resolve_card(card):
            return NLViewerResponse(
                kind="clarify",
                clarification=result.summary or "I couldn't find a matching structure for that.",
            )
        return NLViewerResponse(kind="nav_card", card=card, summary=result.summary)

    # Resolve AlignChain selector actions (organism + active family) to real
    # (rcsb_id, auth_asym_id), same anti-hallucination path as the cards. Drop
    # any that can't be resolved rather than ship a guessed/empty alignment.
    resolved_actions = []
    for a in result.actions:
        if type(a).__name__ == "AlignChain":
            if getattr(a, "organism_id", None) is not None:
                fam = getattr(a, "family", None) or view_context.active_family
                rep = resolve_representative(organism_id=a.organism_id, family=fam)
                if rep is None:
                    continue
                a.rcsb_id, a.auth_asym_id = rep
            elif not (getattr(a, "rcsb_id", None) and getattr(a, "auth_asym_id", None)):
                # No organism selector and no concrete chain — nothing usable.
                continue
        resolved_actions.append(a)

    envelopes = [
        ViewerActionEnvelope(type=type(a).__name__, args=a.model_dump())
        for a in resolved_actions
    ]
    # Entities arrive as dicts (validated loosely in the interpreter); pass
    # through. Per-entity validation against EntityRef is best-effort —
    # malformed entries get filtered.
    entities: List[Dict[str, Any]] = []
    for e in result.entities:
        if isinstance(e, dict) and isinstance(e.get("kind"), str):
            entities.append(e)
    return NLViewerResponse(
        kind="viewer_actions",
        actions=envelopes,
        entities=entities,
        summary=result.summary,
    )


# ---------------------------------------------------------------------------
# Global front-page query
# ---------------------------------------------------------------------------


class NLGlobalRequest(BaseModel):
    text: str = Field(..., min_length=1, max_length=2000)


class NLGlobalResponseBody(BaseModel):
    kind: Literal["global", "clarify"]
    response: Optional[GlobalNLResponse] = None
    clarification: Optional[str] = None


@router_nl_query.post(
    "/global",
    response_model=NLGlobalResponseBody,
    operation_id="nl_query_global",
)
def nl_query_global(req: NLGlobalRequest) -> NLGlobalResponseBody:
    """Front-page entry point. Translates an open-ended question into a
    `GlobalNLResponse` (blurb + queries + ranked action cards), then validates
    entity references against the DB before returning. The frontend executes
    `queries[]` via the existing list endpoints and renders `cards[]` as
    clickable chips that route to the right page with state preloaded.
    """
    try:
        translator = get_translator()
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Translator unavailable: {e}")

    facets = load_facet_context()

    try:
        result = translator.translate_global(text=req.text, facets=facets)
    except Exception as e:
        raise HTTPException(status_code=502, detail=f"Translator error: {e}")

    if result.clarification is not None:
        return NLGlobalResponseBody(kind="clarify", clarification=result.clarification)

    if result.response is None:
        # Translator returned neither — defensive fallback.
        return NLGlobalResponseBody(
            kind="clarify",
            clarification="The model returned no response. Please rephrase.",
        )

    # Resolve LLM-expressed intent (organism + family + ligand) to REAL
    # structure ids BEFORE hydration. This is what stops the model from
    # citing guessed/wrong-organism PDB ids: entity cards get their ids from
    # the DB, not from the model. Hydration then existence-checks whatever
    # survives (user-named ids, resolved ids) as a second line of defense.
    resolved = resolve_response(result.response)
    hydrated = hydrate_response(resolved, known_families=facets.tubulin_families)
    return NLGlobalResponseBody(kind="global", response=hydrated)
