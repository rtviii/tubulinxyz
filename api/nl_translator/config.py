"""Factory that picks a translator based on the NL_TRANSLATOR_PROVIDER env var."""
from __future__ import annotations

import os
from functools import lru_cache

from api.nl_translator.interface import NLTranslator


@lru_cache(maxsize=1)
def get_translator() -> NLTranslator:
    provider = os.environ.get("NL_TRANSLATOR_PROVIDER", "anthropic").lower().strip()

    if provider == "anthropic":
        from api.nl_translator.anthropic_impl import AnthropicNLTranslator
        return AnthropicNLTranslator()

    if provider in ("openai_compat", "openai-compat", "vllm"):
        from api.nl_translator.openai_compat_impl import OpenAICompatNLTranslator
        return OpenAICompatNLTranslator()

    raise ValueError(
        f"Unknown NL_TRANSLATOR_PROVIDER={provider!r}. "
        "Expected one of: anthropic, openai_compat."
    )
