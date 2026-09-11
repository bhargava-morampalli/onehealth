"""User override TSV. A complete override always wins; gaps go to review."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from onehealth.constants import HEALTHCARE_NO, HEALTHCARE_YES, L1_VALUES
from onehealth.ids import parse_id
from onehealth.io_tsv import first_present, get_field, read_tsv, row_match_ids
from onehealth.text import is_blank, normalize_text


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
    status: str  # applied | empty | conflict | invalid | missing
    reason: str
    row: OverrideRow | None = None


def load_overrides(path: Path | None) -> dict[str, list[OverrideRow]]:
    if path is None:
        return {}
    rows = read_tsv(path)
    grouped: dict[str, list[OverrideRow]] = defaultdict(list)
    for raw in rows:
        key = _override_key(raw)
        if not key:
            continue
        l1 = get_field(raw, ("l1", "L1", "source_l1", "classification"))
        healthcare = normalize_text(get_field(raw, ("healthcare_environment", "healthcare_env")))
        paper = get_field(raw, ("paper_rollup", "paper"))
        note = get_field(raw, ("note", "reason", "comment"))
        complete = not is_blank(l1)
        grouped[key].append(
            OverrideRow(
                key=key,
                l1=l1.strip(),
                healthcare_environment=healthcare,
                paper_rollup=paper.strip(),
                note=note,
                complete=complete,
            )
        )
    return grouped


def lookup_override(record: dict[str, str], overrides: dict[str, list[OverrideRow]]) -> OverrideHit | None:
    if not overrides:
        return None
    keys = {item.upper() for item in row_match_ids(record)}
    matches: list[OverrideRow] = []
    for key in keys:
        matches.extend(overrides.get(key, []))
    if not matches:
        return None

    if any(not item.complete for item in matches):
        return OverrideHit(status="empty", reason="override_empty", row=matches[0])

    invalid = [item for item in matches if item.l1 not in L1_VALUES]
    if invalid:
        return OverrideHit(status="invalid", reason="override_invalid_l1", row=invalid[0])

    healthcare_values = set()
    for item in matches:
        if item.healthcare_environment and item.healthcare_environment not in {
            HEALTHCARE_YES,
            HEALTHCARE_NO,
        }:
            return OverrideHit(status="invalid", reason="override_invalid_healthcare", row=item)
        healthcare_values.add(item.healthcare_environment or HEALTHCARE_NO)

    l1_values = {item.l1 for item in matches}
    paper_values = {item.paper_rollup for item in matches if item.paper_rollup}
    if len(l1_values) > 1 or len(healthcare_values) > 1 or len(paper_values) > 1:
        return OverrideHit(status="conflict", reason="override_conflict", row=matches[0])

    return OverrideHit(status="applied", reason="override", row=matches[0])


def _override_key(row: dict[str, str]) -> str:
    raw = get_field(row, ("id", "query_id", "accession", "biosample", "assembly", "sra", "pdt"))
    if not raw:
        # last resort: first column
        key = first_present(row, row.keys())
        raw = row.get(key, "") if key else ""
    parsed = parse_id(raw)
    return (parsed.normalized or raw).strip().upper()
