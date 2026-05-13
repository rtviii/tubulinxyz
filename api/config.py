import platform
import sys
from pathlib import Path

# Calculate project root (assuming tubulin_api/ is one level deep)
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent

sys.path.insert(0, str(PROJECT_ROOT))


def resolve_muscle_binary(project_root: Path) -> Path:
    """Pick the MUSCLE 3.8.1 binary that matches the running platform.

    - Linux (Docker, PSI's VM): bin/muscle3.8.1 (ELF x86_64)
    - macOS (local dev):        bin/muscle3.8.1-darwin (Mach-O x86_64; runs via
                                Rosetta on Apple Silicon)

    The image baked at GHCR only ever runs Linux, so the darwin binary is dead
    weight inside the container (~430 KB) -- acceptable cost for keeping local
    dev workflows symmetric with prod.
    """
    if platform.system() == "Darwin":
        return project_root / "bin" / "muscle3.8.1-darwin"
    return project_root / "bin" / "muscle3.8.1"


class Settings:
    PROJECT_ROOT  : Path = PROJECT_ROOT
    MUSCLE_BINARY : str  = str(resolve_muscle_binary(PROJECT_ROOT))
    MASTER_PROFILE: str  = str(PROJECT_ROOT / "data" / "alpha_tubulin" / "alpha_tubulin.afasta")
    CACHE_DIR     : Path = PROJECT_ROOT / "cache"
    CACHE_DIR.mkdir(exist_ok=True)

settings = Settings()