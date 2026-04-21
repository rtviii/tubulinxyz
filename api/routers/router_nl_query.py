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


class NLViewerRequest(BaseModel):
    text: str = Field(..., min_length=1, max_length=2000)
    view_context: ViewContextBody = Field(default_factory=ViewContextBody)


class ViewerActionEnvelope(BaseModel):
    """One action to execute, typed by `type` (the Pydantic class name)."""
    type: str
    args: Dict[str, Any] = Field(default_factory=dict)


class NLViewerResponse(BaseModel):
    kind: Literal["viewer_actions", "clarify"]
    actions: List[ViewerActionEnvelope] = Field(default_factory=list)
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
    )

    try:
        result = translator.translate_viewer(text=req.text, view_context=view_context)
    except Exception as e:
        raise HTTPException(status_code=502, detail=f"Translator error: {e}")

    if result.clarification is not None:
        return NLViewerResponse(
            kind="clarify",
            clarification=result.clarification,
            summary=result.summary,
        )

    envelopes = [
        ViewerActionEnvelope(type=type(a).__name__, args=a.model_dump())
        for a in result.actions
    ]
    return NLViewerResponse(
        kind="viewer_actions",
        actions=envelopes,
        summary=result.summary,
    )
