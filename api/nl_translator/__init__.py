"""Natural-language-to-filter translator package.

Exposes a provider-agnostic NLTranslator interface and factory.
"""
from api.nl_translator.interface import (
    NLTranslator,
    TranslationResult,
    Target,
    FacetContext,
)
from api.nl_translator.config import get_translator

__all__ = [
    "NLTranslator",
    "TranslationResult",
    "Target",
    "FacetContext",
    "get_translator",
]
