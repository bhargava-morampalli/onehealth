"""Load and apply the lexicon. onehealth ships one, and its `version` is recorded with every result."""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml

from onehealth.constants import HEALTHCARE_NO, HEALTHCARE_YES, L1_VALUES, SEED_LEXICON_VERSION
from onehealth.errors import InputError
from onehealth.io_tsv import read_text
from onehealth.taxonomy import KINDS, taxon_kind
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
    # "animal" or "human": the rule also matches a field whose whole value is a name
    # NCBI Taxonomy gives to an animal, or to people (onehealth.taxonomy).
    taxonomy: str = ""
    # The rule says whose specimen it is ("human feces") and nothing else: it decides
    # Human or Animal for a specimen rule that matched the same row, and is ignored
    # everywhere else ("human contact surface" is not a person).
    specimen_only: bool = False


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
    """The lexicon that ships with onehealth. There is one; its history is in the repository."""
    return Path(str(files("onehealth").joinpath("data/lexicon.yaml")))


def _lexicon_file(path: Path | None) -> Path:
    """The shipped lexicon, or a YAML file of the user's.

    A version number in place of a file ("v3") used to select an older lexicon that
    shipped alongside the current one. Only the current number still works.
    """
    if path is None:
        return default_lexicon_path()
    version = re.fullmatch(r"v?(\d+)", str(path))
    if version and not path.exists():
        if version.group(1) == SEED_LEXICON_VERSION:
            return default_lexicon_path()
        raise InputError(
            f"no lexicon {path}: onehealth ships one lexicon, version {SEED_LEXICON_VERSION}, which is the default. "
            "To reproduce results made with an earlier version, run the onehealth release or commit that made "
            "them, or pass that version's YAML file from the repository history"
        )
    return path


def load_lexicon(path: Path | None = None) -> Lexicon:
    lexicon_path = _lexicon_file(path)
    text = read_text(lexicon_path)
    try:
        payload = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise InputError(f"{lexicon_path}: not valid YAML ({' '.join(str(exc).split())})") from None
    if not isinstance(payload, dict):
        raise InputError(f"Lexicon is not a mapping: {lexicon_path}")
    version = _recorded_version(payload.get("version"), text)
    items = payload.get("rules") or []
    if not isinstance(items, list):
        raise InputError(f"{lexicon_path}: 'rules' must be a list of rules")
    rules = tuple(_parse_rule(item, number, lexicon_path) for number, item in enumerate(items, start=1))
    return Lexicon(version=version, rules=rules)


def _recorded_version(declared: Any, text: str) -> str:
    """What goes in the lexicon_version column, so a table always says which rules made it.

    A file's own `version` is used as written, with two exceptions. A file with no
    version is `custom-<digest>`; it used to be recorded as the shipped version. A file
    that carries the shipped lexicon's version but has different contents (an edited
    copy) is `<version>-modified-<digest>`. The digest is the first 8 hex digits of the
    file's SHA-256, with line endings normalised.
    """
    content = text.replace("\r\n", "\n").replace("\r", "\n")
    digest = hashlib.sha256(content.encode("utf-8")).hexdigest()[:8]
    version = "" if declared is None else str(declared).strip()
    if not version:
        return f"custom-{digest}"
    if version == SEED_LEXICON_VERSION:
        shipped = read_text(default_lexicon_path()).replace("\r\n", "\n").replace("\r", "\n")
        if shipped != content:
            return f"{version}-modified-{digest}"
    return version


def _parse_rule(item: Any, number: int, path: Path) -> LexiconRule:
    if not isinstance(item, dict):
        raise InputError(f"{path}: rule {number} is not a mapping (id, l1, phrases, ...)")
    for key in ("id", "l1"):
        if item.get(key) in (None, ""):
            raise InputError(f"{path}: rule {number} has no {key}")
    rule_id = str(item["id"])
    where = f"{path}: rule {rule_id}"
    l1 = str(item["l1"])
    if l1 not in L1_VALUES:
        raise InputError(f"{where} has invalid l1 {l1!r}; use one of {', '.join(L1_VALUES)}")
    healthcare = _healthcare(item.get("healthcare_environment"))
    if healthcare not in {HEALTHCARE_YES, HEALTHCARE_NO}:
        raise InputError(f"{where} has invalid healthcare_environment {healthcare!r}")
    fields = tuple(_words(item.get("fields"), where, "fields") or ["isolation_source"])
    unknown = [name for name in fields if name not in RULE_FIELDS]
    if unknown:
        raise InputError(f"{where} reads unknown field(s) {unknown}; use {list(RULE_FIELDS)}")
    try:
        priority = int(item.get("priority") or 0)
    except (TypeError, ValueError):
        raise InputError(f"{where}: priority must be a whole number, not {item.get('priority')!r}") from None
    return LexiconRule(
        rule_id=rule_id,
        fields=fields,
        phrases=tuple(tuple(tokenize(phrase)) for phrase in _words(item.get("phrases"), where, "phrases")),
        exact=tuple(tuple(tokenize(value)) for value in _words(item.get("exact"), where, "exact")),
        l1=l1,
        healthcare_environment=healthcare,
        priority=priority,
        host_dependent=bool(item.get("host_dependent", False)),
        note=str(item.get("note") or ""),
        taxonomy=_taxonomy(item.get("taxonomy"), where),
        specimen_only=bool(item.get("specimen_only", False)),
    )


def _taxonomy(value: Any, where: str) -> str:
    if value in (None, ""):
        return ""
    kind = str(value).strip().lower()
    if isinstance(value, bool) or kind not in KINDS:
        raise InputError(f"{where}: taxonomy must be one of {', '.join(KINDS)}, not {value!r}")
    return kind


def _words(value: Any, where: str, key: str) -> list[str]:
    """A list of words or phrases. A single one may be written on its own ("phrases: blood")."""
    if value is None:
        return []
    items = [value] if isinstance(value, (str, int, float)) else value
    if not isinstance(items, list):
        raise InputError(f"{where}: {key} must be a list of words or phrases")
    for entry in items:
        # YAML reads unquoted yes/no/on/off as true/false, which is never what a lexicon means.
        if isinstance(entry, bool) or not isinstance(entry, (str, int, float)):
            raise InputError(f"{where}: {key} entry {entry!r} is not text; put it in quotes")
    return [str(entry) for entry in items]


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
            if matched is None and rule.taxonomy and taxon_kind(record.get(field, "")) == rule.taxonomy:
                matched = tokens
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
