"""Pydantic models for molstar viewer actions.

Each class is one tool the LLM can call. The translator walks the list in
`VIEWER_ACTION_MODELS`, renders each as a tool spec, and lets the model emit
ONE OR MORE tool calls per user turn (so "hide B and zoom to A" becomes two
actions). The frontend executes them in order against its `MolstarInstance`.

Field names mirror the frontend method arguments verbatim, so the frontend
dispatcher is a dumb switch on `type`.
"""
from __future__ import annotations

import json as _json
from typing import Annotated, Any, List, Literal, Optional, Tuple, Type, Union

from pydantic import BaseModel, Field, field_validator

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
# Annotation tracks — chain-independent aux rows painted on MSA master columns.
#
# TrackSpec mirrors api/routers/router_annotations.py:VariantFilterSpec /
# ModificationFilterSpec / BindingContactFilterSpec. We re-declare here (rather
# than import) to keep viewer_actions.py independent of routers, matching the
# pattern used for EntityRef/ActionCard. Kept narrow on purpose — the LLM
# should not need to think about every filter; the truly load-bearing fields
# are family + (modification_types | chemical_ids | species_tax_ids | etc.).
# ---------------------------------------------------------------------------

class VariantSpec(BaseModel):
    kind: Literal['variants'] = 'variants'
    family: str = Field(..., description="e.g. 'tubulin_alpha', 'tubulin_beta'")
    wild_type_aas: Optional[List[str]] = None
    observed_aas: Optional[List[str]] = None
    substitution_pairs: Optional[List[Tuple[str, str]]] = None
    indel_present: Optional[bool] = None
    sources: Optional[List[Literal['structural', 'literature']]] = None
    uniprot_ids: Optional[List[str]] = None
    species_names: Optional[List[str]] = Field(
        None, description="Text species match (e.g. 'H. sapiens'). Prefer species_tax_ids for structural variants."
    )
    species_tax_ids: Optional[List[int]] = Field(
        None, description="NCBI tax ids; restricts structural variants by entity src_organism_ids (e.g. [9606])."
    )
    position_range: Optional[Tuple[int, int]] = None
    positions: Optional[List[int]] = Field(
        None, description="Restrict to specific master indices (e.g. residues of a binding site)."
    )
    co_occurs_with_mod_type: Optional[List[str]] = None
    phenotype_contains: Optional[List[str]] = Field(
        None, description="Case-insensitive substring search; OR within array. Useful with literature variants."
    )


class ModificationSpec(BaseModel):
    kind: Literal['modifications'] = 'modifications'
    family: str
    modification_types: Optional[List[str]] = Field(
        None, description="e.g. ['phosphorylation', 'acetylation', 'polyglutamylation']"
    )
    uniprot_ids: Optional[List[str]] = None
    species_tax_ids: Optional[List[int]] = Field(
        None, description="NCBI tax ids; modifications carry tax_id directly."
    )
    species_names: Optional[List[str]] = None
    position_range: Optional[Tuple[int, int]] = None
    positions: Optional[List[int]] = None
    co_occurs_with_variant: Optional[bool] = None
    evidence_source: Optional[List[str]] = None
    phenotype_contains: Optional[List[str]] = None


class BindingContactSpec(BaseModel):
    kind: Literal['binding_contacts'] = 'binding_contacts'
    family: str
    chemical_ids: List[str] = Field(..., description="Ligand chem comp ids (e.g. ['TA1', 'GTP', 'TXL']).")
    structure_ids: Optional[List[str]] = None
    positions: Optional[List[int]] = None


TrackSpec = Annotated[
    Union[VariantSpec, ModificationSpec, BindingContactSpec],
    Field(discriminator='kind'),
]


class AddAnnotationTrack(BaseModel):
    """Add a chain-independent annotation row to the MSA aux panel. The track
    paints master columns matching `spec` (variants / PTMs / ligand contacts)
    on every displayed chain of `spec.family` and colors the corresponding 3D
    residues on the primary chain. Expert (monomer) view only.

    Use for organism-scoped, family-scoped, or site-scoped annotation overlays:
    "Show PTMs in human alpha tubulin associated with fibrosis", "Where does GTP
    bind on alpha tubulin?", "Compare the GTP site in human vs Toxoplasma" (emit
    two tracks with distinct colors).

    Do NOT navigate away. Do NOT ask for clarification when the request is a
    direct annotation overlay request — emit the track. Multiple AddAnnotationTrack
    calls in a single response are encouraged for comparative queries (one per
    organism / one per ligand / etc.) — pick distinct colors.

    Pick a short human-readable label (e.g. "PTMs · human · fibrosis", "GTP site"
    "GTP site · Toxoplasma"). The user sees this in the MSA labels.
    """
    label: str = Field(..., description="Short label shown in the MSA aux panel (<= 40 chars).")
    spec: TrackSpec = Field(..., description="The typed filter spec; discriminated by 'kind'.")
    color: str = Field(
        ..., description="Hex color (e.g. '#E74C3C'). For comparative tracks, pick visually distinct colors."
    )
    description: Optional[str] = Field(
        None, description="One-sentence rationale for grounding (not shown in UI)."
    )

    @field_validator("spec", mode="before")
    @classmethod
    def _parse_spec_if_string(cls, v: Any) -> Any:
        # Some LLM/tool runtimes hand nested-object tool args back as escaped
        # JSON strings instead of dicts. Decode here so the discriminated
        # union still picks the right variant.
        if isinstance(v, str):
            try:
                return _json.loads(v)
            except (ValueError, TypeError):
                pass
        return v


class RemoveAnnotationTrack(BaseModel):
    """Remove an existing annotation track from the MSA. Substring-matches
    against existing track labels (case-insensitive). Use when the user asks to
    'hide', 'remove', 'drop', or 'clear' a specific track."""
    label_match: str = Field(
        ...,
        description=(
            "Substring of the track label to remove (case-insensitive). E.g. 'GTP site' "
            "matches both 'GTP site' and 'GTP site · Toxoplasma'."
        ),
    )


class FocusBindingSite(BaseModel):
    """Focus the camera + draw the binding-site representation for a named
    ligand on a loaded chain, using contact data already fetched per chain
    (no extra round trip). Also jumps the MSA viewport to the contact range.

    Use whenever the user asks to:
      - "focus / show / highlight the <ligand> binding site"
      - "compare X site in human vs Y" (after AlignChain, focus the site on the
        primary/active chain so the alignment context is visually anchored)
      - "where does <ligand> bind" — when the ligand IS bound on the loaded
        structure (if it isn't, use the highlight + nav-card pattern instead).

    Required: chemical_id (e.g. 'GTP', 'TA1', 'LOC'). Optional auth_asym_id —
    defaults to the active monomer chain from view_context.
    """
    chemical_id: str = Field(..., description="Ligand chem comp id, e.g. 'GTP', 'TA1', 'LOC'.")
    auth_asym_id: Optional[str] = Field(
        default=None,
        description="Chain to focus the binding site on. Defaults to active_monomer_chain.",
    )


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
    AddAnnotationTrack,
    RemoveAnnotationTrack,
    FocusBindingSite,
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
    "AddAnnotationTrack": (
        "Expert mode only. Add a chain-independent MSA aux row painted at master "
        "columns matching a typed FilterSpec (variants / modifications / "
        "binding_contacts). Use for annotation overlays scoped by organism, family, "
        "or site (e.g. 'PTMs in human alpha tubulin', 'GTP binding site', 'variants "
        "in the taxol pocket'). Emit MULTIPLE calls for comparative queries (one "
        "per organism / per ligand) with distinct colors."
    ),
    "RemoveAnnotationTrack": (
        "Remove an existing annotation track by case-insensitive label substring "
        "match. Use when the user asks to hide/clear/drop a track."
    ),
    "FocusBindingSite": (
        "Focus the camera on the binding site of a named ligand on a loaded chain "
        "(uses already-fetched per-chain contact data; no extra round trip). Also "
        "jumps the MSA to the contact span. Pair with AlignChain for comparative "
        "queries so the binding-site context is anchored on the primary chain "
        "after the alignment is loaded."
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
