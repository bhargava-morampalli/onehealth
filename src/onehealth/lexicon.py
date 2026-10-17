"""Load and apply the versioned seed lexicon."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml

from onehealth.constants import HEALTHCARE_NO, HEALTHCARE_YES, L1_VALUES, SEED_LEXICON_VERSION
from onehealth.text import contains_phrase, tokenize

# Metadata fields a rule may read (keys of io_tsv.metadata_view).
RULE_FIELDS = ("isolation_source", "host", "sample_type", "host_disease")


@dataclass(frozen=True)
class LexiconRule:
    rule_id: str
    fields: tuple[str, ...]
    phrases: tuple[tuple[str, ...], ...]
    exact: tuple[tuple[str, ...], ...]
    l1: str
    healthcare_environment: str
    priority: int
    # A body site or specimen (blood, urine, faeces, "clinical") that occurs in any
    # host. L1 then comes from host evidence; `l1` is used only when there is none.
    host_dependent: bool = False
    note: str = ""


@dataclass(frozen=True)
class Lexicon:
    version: str
    rules: tuple[LexiconRule, ...]


def _healthcare(value: Any) -> str:
    if value is True:
        return HEALTHCARE_YES
    if value is False or value is None:
        return HEALTHCARE_NO
    return str(value).strip().lower()


def default_lexicon_path() -> Path:
    return Path(str(files("onehealth").joinpath(f"data/lexicon/v{SEED_LEXICON_VERSION}.yaml")))


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
    unknown = [name for name in fields if name not in RULE_FIELDS]
    if unknown:
        raise ValueError(f"Lexicon rule {rule_id} reads unknown field(s) {unknown}; use {list(RULE_FIELDS)}")
    return LexiconRule(
        rule_id=rule_id,
        fields=fields,
        phrases=tuple(tuple(tokenize(phrase)) for phrase in item.get("phrases") or []),
        exact=tuple(tuple(tokenize(value)) for value in item.get("exact") or []),
        l1=l1,
        healthcare_environment=healthcare,
        priority=int(item.get("priority") or 0),
        host_dependent=bool(item.get("host_dependent", False)),
        note=str(item.get("note") or ""),
    )


@dataclass(frozen=True)
class RuleHit:
    rule: LexiconRule
    field: str
    matched: str


def match_rules(record: dict[str, str], lexicon: Lexicon) -> list[RuleHit]:
    """One hit per rule: the first of the rule's fields that matches."""
    hits: list[RuleHit] = []
    for rule in lexicon.rules:
        for field in rule.fields:
            tokens = tuple(tokenize(record.get(field, "")))
            if not tokens:
                continue
            matched = tokens if tokens in rule.exact else None
            if matched is None:
                matched = next((phrase for phrase in rule.phrases if contains_phrase(list(tokens), phrase)), None)
            if matched is not None:
                hits.append(RuleHit(rule=rule, field=field, matched=" ".join(matched)))
                break
    return hits


def best_hits(hits: list[RuleHit]) -> list[RuleHit]:
    """Hits that share the highest priority."""
    if not hits:
        return []
    top = max(hit.rule.priority for hit in hits)
    return [hit for hit in hits if hit.rule.priority == top]
