"""User override TSV. A complete override always wins; gaps go to review."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from onehealth.constants import HEALTHCARE_NO, HEALTHCARE_YES, L1_VALUES, PAPER_ROLLUP_VALUES
from onehealth.ids import parse_id
from onehealth.errors import InputError
from onehealth.io_tsv import first_present, get_field, read_tsv, row_match_ids
from onehealth.text import normalize_text

_OVERRIDE_ID_COLUMNS = ("id", "query_id", "accession", "biosample", "assembly", "sra", "pdt")
_L1_COLUMNS = ("l1", "source_l1", "classification")


@dataclass(frozen=True)
class OverrideRow:
    key: str
    l1: str
    healthcare_environment: str
    paper_rollup: str
    note: str
    complete: bool


@dataclass(frozen=True)
class OverrideHit:
    status: str  # applied | empty | conflict | invalid
    reason: str
    row: OverrideRow | None = None


def _canonical(value: str, allowed: tuple[str, ...]) -> str:
    """Case-insensitive match to an allowed label; unknown text is returned as typed."""
    text = value.strip()
    return {item.lower(): item for item in allowed}.get(text.lower(), text)


def load_overrides(path: Path | None) -> dict[str, list[OverrideRow]]:
    if path is None:
        return {}
    grouped: dict[str, list[OverrideRow]] = defaultdict(list)
    rows = read_tsv(path)
    # Without these columns every override would be ignored, or sent to review, with no sign why.
    if rows and first_present(rows[0], _OVERRIDE_ID_COLUMNS) is None:
        raise InputError(f"{path}: no ID column ({', '.join(_OVERRIDE_ID_COLUMNS)})")
    if rows and first_present(rows[0], _L1_COLUMNS) is None:
        raise InputError(f"{path}: no l1 column ({', '.join(_L1_COLUMNS)})")
    for raw in rows:
        key = _override_key(raw)
        if not key:
            continue
        l1 = _canonical(get_field(raw, _L1_COLUMNS), L1_VALUES)
        grouped[key].append(
            OverrideRow(
                key=key,
                l1=l1,
                healthcare_environment=normalize_text(get_field(raw, ("healthcare_environment", "healthcare_env"))),
                paper_rollup=_canonical(get_field(raw, ("paper_rollup", "paper")), PAPER_ROLLUP_VALUES),
                note=get_field(raw, ("note", "reason", "comment")),
                complete=bool(l1),
            )
        )
    return grouped


def lookup_override(record: dict[str, str], overrides: dict[str, list[OverrideRow]]) -> OverrideHit | None:
    if not overrides:
        return None
    # row_match_ids includes GCA<->GCF twins and the versionless assembly stem.
    matches = [row for key in sorted(row_match_ids(record)) for row in overrides.get(key, [])]
    if not matches:
        return None

    if any(not item.complete for item in matches):
        return OverrideHit(status="empty", reason="override_empty", row=matches[0])

    for item in matches:
        if item.l1 not in L1_VALUES:
            return OverrideHit(status="invalid", reason="override_invalid_l1", row=item)
        if item.healthcare_environment not in {"", HEALTHCARE_YES, HEALTHCARE_NO}:
            return OverrideHit(status="invalid", reason="override_invalid_healthcare", row=item)
        if item.paper_rollup and item.paper_rollup not in PAPER_ROLLUP_VALUES:
            return OverrideHit(status="invalid", reason="override_invalid_paper_rollup", row=item)

    l1_values = {item.l1 for item in matches}
    healthcare_values = {item.healthcare_environment or HEALTHCARE_NO for item in matches}
    paper_values = {item.paper_rollup for item in matches if item.paper_rollup}
    if len(l1_values) > 1 or len(healthcare_values) > 1 or len(paper_values) > 1:
        return OverrideHit(status="conflict", reason="override_conflict", row=matches[0])
    return OverrideHit(status="applied", reason="override", row=matches[0])


def _override_key(row: dict[str, str]) -> str:
    raw = get_field(row, _OVERRIDE_ID_COLUMNS)
    if not raw:
        return ""
    parsed = parse_id(raw)
    return (parsed.normalized or raw).strip().upper()
