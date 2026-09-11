"""Load and apply the versioned seed lexicon."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml

from onehealth.constants import HEALTHCARE_NO, HEALTHCARE_YES, L1_VALUES, SEED_LEXICON_VERSION
from onehealth.text import contains_phrase, normalize_text, tokenize


@dataclass(frozen=True)
class LexiconRule:
    rule_id: str
    fields: tuple[str, ...]
    phrases: tuple[tuple[str, ...], ...]
    exact: tuple[str, ...]
    l1: str
    healthcare_environment: str
    priority: int
    note: str = ""


@dataclass(frozen=True)
class Lexicon:
    version: str
    rules: tuple[LexiconRule, ...]


FIELD_ALIASES = {
    "isolation_source": ("isolation_source",),
    "host": ("host",),
    "host_disease": ("host_disease",),
    "sample_type": ("sample_type",),
    "any": ("isolation_source", "host", "host_disease", "sample_type"),
}


def _healthcare(value: Any) -> str:
    if value is True:
        return HEALTHCARE_YES
    if value is False or value is None:
        return HEALTHCARE_NO
    return str(value).strip().lower()


def default_lexicon_path() -> Path:
    return Path(str(files("onehealth").joinpath("data/lexicon/v1.yaml")))


def load_lexicon(path: Path | None = None) -> Lexicon:
    lexicon_path = path or default_lexicon_path()
    payload = yaml.safe_load(lexicon_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Lexicon is not a mapping: {lexicon_path}")
    version = str(payload.get("version") or SEED_LEXICON_VERSION)
    rules = tuple(_parse_rule(item) for item in payload.get("rules") or [])
    return Lexicon(version=version, rules=rules)


def _parse_rule(item: dict[str, Any]) -> LexiconRule:
    rule_id = str(item["id"])
    l1 = str(item["l1"])
    if l1 not in L1_VALUES:
        raise ValueError(f"Lexicon rule {rule_id} has invalid l1 {l1!r}")
    healthcare = _healthcare(item.get("healthcare_environment"))
    if healthcare not in {HEALTHCARE_YES, HEALTHCARE_NO}:
        raise ValueError(f"Lexicon rule {rule_id} has invalid healthcare_environment {healthcare!r}")
    fields = tuple(item.get("fields") or ["isolation_source"])
    phrases = tuple(tuple(tokenize(phrase)) for phrase in item.get("phrases") or [])
    exact = tuple(normalize_text(value) for value in item.get("exact") or [])
    return LexiconRule(
        rule_id=rule_id,
        fields=fields,
        phrases=phrases,
        exact=exact,
        l1=l1,
        healthcare_environment=healthcare,
        priority=int(item.get("priority") or 0),
        note=str(item.get("note") or ""),
    )


@dataclass(frozen=True)
class RuleHit:
    rule: LexiconRule
    field: str
    matched: str


def match_rules(record: dict[str, str], lexicon: Lexicon) -> list[RuleHit]:
    hits: list[RuleHit] = []
    for rule in lexicon.rules:
        for field in rule.fields:
            for actual in FIELD_ALIASES.get(field, (field,)):
                value = record.get(actual, "")
                tokens = tokenize(value)
                normalized = normalize_text(value)
                matched = ""
                if rule.exact and normalized in rule.exact:
                    matched = normalized
                else:
                    for phrase in rule.phrases:
                        if contains_phrase(tokens, phrase):
                            matched = " ".join(phrase)
                            break
                if matched:
                    hits.append(RuleHit(rule=rule, field=actual, matched=matched))
                    break
    return hits


def best_hits(hits: list[RuleHit]) -> list[RuleHit]:
    if not hits:
        return []
    top = max(hit.rule.priority for hit in hits)
    winners = [hit for hit in hits if hit.rule.priority == top]
    # Deduplicate identical outcomes.
    seen: set[tuple[str, str, str]] = set()
    unique: list[RuleHit] = []
    for hit in winners:
        key = (hit.rule.l1, hit.rule.healthcare_environment, hit.rule.rule_id)
        if key in seen:
            continue
        seen.add(key)
        unique.append(hit)
    return unique
