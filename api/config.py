import sys
from pathlib import Path

# Calculate project root (assuming tubulin_api/ is one level deep)
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent

# Add root to sys.path so we can import sibling modules like tubulin_analyzer
sys.path.insert(0, str(PROJECT_ROOT))

class Settings:

    MUSCLE_BINARY : str  = str(PROJECT_ROOT / "muscle3.8.1")
    MASTER_PROFILE: str  = str(PROJECT_ROOT / "data" / "alpha_tubulin" / "alpha_tubulin.afasta")
    CACHE_DIR     : Path = PROJECT_ROOT / "cache"

    CACHE_DIR.mkdir(exist_ok=True)

settings = Settings()