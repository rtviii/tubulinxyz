"""Pydantic models for molstar viewer actions.

Each class is one tool the LLM can call. The translator walks the list in
`VIEWER_ACTION_MODELS`, renders each as a tool spec, and lets the model emit
ONE OR MORE tool calls per user turn (so "hide B and zoom to A" becomes two
actions). The frontend executes them in order against its `MolstarInstance`.

Field names mirror the frontend method arguments verbatim, so the frontend
dispatcher is a dumb switch on `type`.
"""
from __future__ import annotations

from typing import List, Type

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Camera / focus
# ---------------------------------------------------------------------------

class FocusChain(BaseModel):
    """Zoom the camera to a polymer chain."""
    auth_asym_id: str = Field(..., description="Author asym id, e.g. 'A'.")


class FocusResidue(BaseModel):
    """Zoom the camera to a single residue on a chain."""
    auth_asym_id: str
    auth_seq_id: int = Field(..., description="Author sequence id of the residue.")


class FocusResidueRange(BaseModel):
    """Zoom the camera to a contiguous residue range (inclusive)."""
    auth_asym_id: str
    start: int
    end: int


class ClearFocus(BaseModel):
    """Reset the camera to the default full-structure view."""
    pass


# ---------------------------------------------------------------------------
# Visibility
# ---------------------------------------------------------------------------

class SetChainVisibility(BaseModel):
    """Show or hide a single polymer chain."""
    auth_asym_id: str
    visible: bool


class IsolateChain(BaseModel):
    """Hide every chain except the one named. Ligands kept by default."""
    auth_asym_id: str
    keep_ligands: bool = True


# ---------------------------------------------------------------------------
# Highlight (transient glow, non-mutating)
# ---------------------------------------------------------------------------

class HighlightChain(BaseModel):
    """Glow-highlight an entire chain. Pair with ClearHighlight to remove."""
    auth_asym_id: str


class HighlightResidueRange(BaseModel):
    """Glow-highlight a contiguous residue range."""
    auth_asym_id: str
    start: int
    end: int


class ClearHighlight(BaseModel):
    """Remove all transient highlights."""
    pass


# ---------------------------------------------------------------------------
# Clarification — no-op action, frontend renders the question.
# ---------------------------------------------------------------------------

class RequestClarification(BaseModel):
    """Ask the user a clarifying question instead of taking any action."""
    question: str = Field(..., description="One-sentence question for the user.")


# ---------------------------------------------------------------------------
# Registry — drives tool generation in both translator impls.
# ---------------------------------------------------------------------------

VIEWER_ACTION_MODELS: List[Type[BaseModel]] = [
    FocusChain,
    FocusResidue,
    FocusResidueRange,
    ClearFocus,
    SetChainVisibility,
    IsolateChain,
    HighlightChain,
    HighlightResidueRange,
    ClearHighlight,
    RequestClarification,
]

# Description shown to the LLM for each tool. Kept terse — the field docstrings
# already carry most of the semantics.
VIEWER_ACTION_DESCRIPTIONS: dict[str, str] = {
    "FocusChain": "Zoom the camera to a polymer chain.",
    "FocusResidue": "Zoom the camera to one residue on a chain.",
    "FocusResidueRange": "Zoom the camera to a residue range on a chain (inclusive).",
    "ClearFocus": "Reset the camera to the default full-structure view.",
    "SetChainVisibility": "Show or hide a single polymer chain.",
    "IsolateChain": "Hide every chain except the named one; ligands kept by default.",
    "HighlightChain": "Glow-highlight an entire chain.",
    "HighlightResidueRange": "Glow-highlight a residue range.",
    "ClearHighlight": "Remove all transient highlights.",
    "RequestClarification": (
        "Ask the user a clarifying question. Use when the required chain/ligand "
        "is not present in view_context, or when the request is ambiguous."
    ),
}

CLARIFY_ACTION_NAME = "RequestClarification"
