"""Pydantic models for molstar viewer actions.

Each class is one tool the LLM can call. The translator walks the list in
`VIEWER_ACTION_MODELS`, renders each as a tool spec, and lets the model emit
ONE OR MORE tool calls per user turn (so "hide B and zoom to A" becomes two
actions). The frontend executes them in order against its `MolstarInstance`.

Field names mirror the frontend method arguments verbatim, so the frontend
dispatcher is a dumb switch on `type`.
"""
from __future__ import annotations

from typing import List, Optional, Type

from pydantic import BaseModel, Field

from api.nl_translator.global_actions import EntityRef, ActionCard


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
# Alignment — add another organism's chain to the current expert-mode view.
# ---------------------------------------------------------------------------

class AlignChain(BaseModel):
    """Add another organism's chain to the CURRENT expert-mode alignment — a new
    aligned row in the MSA plus a ghost overlay in 3D, keeping everything already
    loaded. Use when the user asks to add / load / include another organism's
    sequence or structure while in expert (monomer) mode (e.g. "add a bovine
    sequence", "also show human alpha", "compare with yeast tubulin"). Do NOT
    navigate away and do NOT ask for clarification — this adds in place.

    Express the organism as an NCBI tax id; the backend resolves it to a real
    structure + chain of the active family. Leave rcsb_id / auth_asym_id null.
    """
    organism_id: Optional[int] = Field(
        default=None,
        description="NCBI tax id of the organism to add (e.g. 9606 human, 9913 cattle/bovine, 4932 yeast). Backend resolves it to a real structure+chain of the active family.",
    )
    family: Optional[str] = Field(
        default=None,
        description="Tubulin family; leave null to use the active chain's family from view_context.",
    )
    # Server-filled after resolution; the model leaves these null.
    rcsb_id: Optional[str] = None
    auth_asym_id: Optional[str] = None


# ---------------------------------------------------------------------------
# Clarification — no-op action, frontend renders the question.
# ---------------------------------------------------------------------------

class RequestClarification(BaseModel):
    """Ask the user a clarifying question instead of taking any action."""
    question: str = Field(..., description="One-sentence question for the user.")


# ---------------------------------------------------------------------------
# Entity surfacing — companion tool that names the things the LLM acted on.
# Not a viewer action (it changes nothing in molstar). The frontend renders
# these as interactive pills with bidirectional sync to the viewer.
# ---------------------------------------------------------------------------

class MentionEntities(BaseModel):
    """Surface the chains / residues / ranges / ligands the user should see as
    interactive pills. Call ONCE per response, alongside your viewer-action
    tool calls. Do not duplicate entities already implied by actions — just
    list the things you want surfaced for hover/click in the side panel.
    """
    entities: List[EntityRef] = Field(
        default_factory=list,
        description="Up to 6 entities to render as pills.",
    )


class EmitNavigationCard(BaseModel):
    """Use INSTEAD of viewer actions when the user's intent is navigation —
    they're asking about other structures, want to browse the catalogue, or
    want to switch to a different PDB entry. The card uses the same vocabulary
    as the global front-page endpoint and routes the user away from this
    structure page. Pick ONE card and skip all action tools.
    """
    card: ActionCard = Field(
        ...,
        description=(
            "A single ActionCard. Typical choices on a structure page: "
            "open_catalogue (browse), open_structure (different PDB), or "
            "open_expert (different chain or alignment). Pick the most "
            "specific that fits the user's intent."
        ),
    )


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
    AlignChain,
    MentionEntities,
    EmitNavigationCard,
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
    "AlignChain": (
        "Expert mode only. Add another organism's chain to the CURRENT alignment "
        "(new MSA row + 3D ghost overlay) without navigating away. Use for 'add a "
        "<organism> sequence/structure'. Set organism_id (NCBI tax id); backend "
        "resolves it to a real structure+chain of the active family."
    ),
    "MentionEntities": (
        "Surface the entities (chains, residues, ranges, ligands) the user "
        "should see as interactive pills in the side panel. Call ONCE per "
        "response if any entities are worth pinning."
    ),
    "EmitNavigationCard": (
        "Use INSTEAD of viewer actions when the user's intent is navigation "
        "(browse catalogue, switch to another structure, open a different "
        "chain in expert mode). Emits a single ActionCard the frontend turns "
        "into a clickable chip. Skip all action tools when using this."
    ),
    "RequestClarification": (
        "Ask the user a clarifying question. Use when the required chain/ligand "
        "is not present in view_context, or when the request is ambiguous."
    ),
}

CLARIFY_ACTION_NAME = "RequestClarification"
MENTION_ENTITIES_NAME = "MentionEntities"
EMIT_NAV_CARD_NAME = "EmitNavigationCard"
