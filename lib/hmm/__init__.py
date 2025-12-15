"""
HMM module configuration.
Modify these getters when relocating seed data / binaries.
"""
from pathlib import Path

# Base paths - modify these when relocating
_PROJECT_ROOT = Path(__file__).parent.parent.parent  # tubulinxyz/
_PROCESSED_DIR = _PROJECT_ROOT / "processed"
_RAW_DIR = _PROJECT_ROOT / "raw_fastas"
_HMM_DIR = _PROJECT_ROOT / "hmms"
_MUSCLE_BIN = _PROJECT_ROOT / "muscle3.8.1"


def get_processed_dir() -> Path:
    """Directory containing balanced seed FASTAs."""
    return _PROCESSED_DIR


def get_raw_dir() -> Path:
    """Directory containing raw unfiltered FASTAs from UniProt."""
    return _RAW_DIR


def get_hmm_dir() -> Path:
    """Directory for storing/loading HMM files."""
    return _HMM_DIR


def get_muscle_bin() -> Path:
    """Path to MUSCLE binary."""
    return _MUSCLE_BIN


def get_seed_fasta_path(family: str) -> Path:
    """Get path to seed FASTA for a tubulin family."""
    return get_processed_dir() / f"tubulin_{family}_final.fasta"


def get_raw_fasta_path(family: str) -> Path:
    """Get path to raw FASTA for a tubulin family."""
    return get_raw_dir() / f"tubulin_{family}_raw.fasta"


def get_hmm_path(family: str) -> Path:
    """Get path to HMM file for a tubulin family."""
    return get_hmm_dir() / f"tubulin_{family}.hmm"