"""Assistant endpoint — the single grounded orchestrator entrypoint.

POST /assistant/query  { text, page_context? }  ->  AssistantResult

This is the consolidation target for the three legacy /nl_query/* endpoints.
They remain live during migration; the frontend is cut over page-by-page.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, Optional

from fastapi import APIRouter
from pydantic import BaseModel

from api.nl_translator.orchestrator import AssistantResult, run_assistant

logger = logging.getLogger(__name__)

router_assistant = APIRouter()


class AssistantQueryRequest(BaseModel):
    text: str
    page_context: Optional[Dict[str, Any]] = None


@router_assistant.post("/query", response_model=AssistantResult, operation_id="assistant_query")
def assistant_query(req: AssistantQueryRequest) -> AssistantResult:
    """Run the grounded assistant loop. Always returns an AssistantResult — on
    any internal failure we return kind='cannot' with a reason so the UI has a
    usable shape rather than a 500."""
    try:
        return run_assistant(req.text, req.page_context)
    except Exception as e:  # config error, provider error, etc.
        logger.exception("assistant_query failed")
        return AssistantResult(kind="cannot", reason=f"Assistant error: {type(e).__name__}: {str(e)[:200]}")
