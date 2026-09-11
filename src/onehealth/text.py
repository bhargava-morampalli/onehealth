"""Normalization helpers for isolation-source matching."""

from __future__ import annotations

import re

from onehealth.constants import UNKNOWN_EXACT, UNKNOWN_GENERIC_TAILS

_SPLIT_RE = re.compile(r"[^a-z0-9]+")


def normalize_text(value: str | None) -> str:
    text = (value or "").strip().lower()
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text)
    return text


def tokenize(value: str | None) -> list[str]:
    normalized = normalize_text(value)
    if not normalized:
        return []
    return [part for part in _SPLIT_RE.split(normalized) if part]


def contains_phrase(haystack: list[str], phrase: list[str] | tuple[str, ...]) -> bool:
    needle = list(phrase)
    if not needle or len(needle) > len(haystack):
        return False
    n = len(needle)
    for i in range(len(haystack) - n + 1):
        if haystack[i : i + n] == needle:
            return True
    return False


def is_blank(value: str | None) -> bool:
    return normalize_text(value) == ""


def is_unknown_value(value: str | None) -> bool:
    """True when a field is empty, an unknown marker, or 'unknown <generic>'."""
    normalized = normalize_text(value)
    if normalized in UNKNOWN_EXACT:
        return True
    tokens = tokenize(normalized)
    if not tokens:
        return True
    if tokens[0] in {"unknown", "unk"}:
        return True
    if tokens[-1] in {"unknown", "unk"} and set(tokens[:-1]).issubset(UNKNOWN_GENERIC_TAILS):
        return True
    return False
