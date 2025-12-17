"""
HMM module configuration.
Handles path resolution for Tubulin and MAP families.
"""
from pathlib import Path
from typing import Union
from lib.models.types_tubulin import TubulinFamily, MapFamily, HmmFamily

# --- Base Paths ---
_PROJECT_ROOT = Path(__file__).parent.parent.parent  # tubulinxyz/

_DATA_DIR = _PROJECT_ROOT / "data"

# Sequence Directories
_SEQ_TUBULIN_DIR = _DATA_DIR / "sequences" / "tubulin"
_SEQ_MAPS_DIR    = _DATA_DIR / "sequences" / "maps"

# HMM Directories
_HMM_TUBULIN_DIR = _DATA_DIR / "hmms" / "tubulin"
_HMM_MAPS_DIR    = _DATA_DIR / "hmms" / "maps"

# Binaries
_MUSCLE_BIN      = _PROJECT_ROOT / "muscle3.8.1"


# --- Getters ---

def get_muscle_bin() -> Path:
    """Path to MUSCLE binary."""
    return _MUSCLE_BIN

def _resolve_family_dirs(family: HmmFamily) -> tuple[Path, Path]:
    """
    Internal helper to route to the correct sequence/hmm folders 
    based on the family type.
    Returns: (sequence_dir, hmm_dir)
    """
    if isinstance(family, TubulinFamily):
        return _SEQ_TUBULIN_DIR, _HMM_TUBULIN_DIR
    elif isinstance(family, MapFamily):
        return _SEQ_MAPS_DIR, _HMM_MAPS_DIR
    else:
        raise ValueError(f"Unknown family type: {type(family)}")

def get_fasta_path(family: HmmFamily) -> Path:
    """
    Get path to input FASTA (seed sequences).
    Expected format: data/sequences/{type}/{family_value}.fasta
    """
    seq_dir, _ = _resolve_family_dirs(family)
    # Ensure directory exists when getting path
    seq_dir.mkdir(parents=True, exist_ok=True)
    return seq_dir / f"{family.value}.fasta"

def get_aligned_fasta_path(family: HmmFamily) -> Path:
    """
    Get path to output aligned FASTA.
    Format: data/sequences/{type}/{family_value}.afasta
    """
    seq_dir, _ = _resolve_family_dirs(family)
    seq_dir.mkdir(parents=True, exist_ok=True)
    return seq_dir / f"{family.value}.afasta"

def get_hmm_path(family: HmmFamily) -> Path:
    """
    Get path to HMM file.
    Format: data/hmms/{type}/{family_value}.hmm
    """
    _, hmm_dir = _resolve_family_dirs(family)
    hmm_dir.mkdir(parents=True, exist_ok=True)
    return hmm_dir / f"{family.value}.hmm"

def get_all_hmm_dirs() -> list[Path]:
    """Returns a list of all directories where HMMs might reside."""
    return [_HMM_TUBULIN_DIR, _HMM_MAPS_DIR]