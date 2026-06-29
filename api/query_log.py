"""Best-effort capture of every user query that hits an assistant/NL endpoint.

Appends one JSON object per line to a JSONL file. Default path
$TUBETL_DATA/query_log.jsonl, overridable via TUBXZ_QUERY_LOG. Never raises -- a
logging failure must not break a user's query.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


def _log_path() -> Path:
    override = os.environ.get("TUBXZ_QUERY_LOG")
    if override:
        return Path(override)
    return Path(os.environ.get("TUBETL_DATA") or ".") / "query_log.jsonl"


def log_query(source: str, text: str, context: Optional[dict[str, Any]] = None) -> None:
    """Append one query record to the JSONL log. Swallows all errors."""
    try:
        rec = {
            "ts": datetime.now(timezone.utc).isoformat(),
            "source": source,
            "text": text,
            "context": context or {},
        }
        p = _log_path()
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    except Exception:
        pass
