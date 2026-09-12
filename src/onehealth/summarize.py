"""Count classifications. Cluster/plasmid overlap only if those columns exist."""

from __future__ import annotations

from collections import Counter

from onehealth.constants import CLUSTER_ALIASES, ISOLATION_SOURCE_ALIASES, PLASMID_ALIASES
from onehealth.io_tsv import find_column, get_field
from onehealth.text import is_blank

# Always printed so classify_status=ok is not read as "source known".
SOURCE_EVIDENCE_CORE = (
    "ok_with_source",
    "ok_empty_source",
    "unmatched",
    "conflict",
)


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    out.extend(_counts(rows, "l1", "l1"))
    out.extend(_counts(rows, "paper_rollup", "paper_rollup"))
    out.extend(_counts(rows, "healthcare_environment", "healthcare_environment"))
    out.extend(_counts(rows, "classify_status", "classify_status"))
    out.extend(_counts(rows, "classify_reason", "classify_reason"))
    out.extend(_source_evidence(rows))

    cluster_col = find_column(rows, CLUSTER_ALIASES)
    if cluster_col:
        out.extend(_overlap(rows, cluster_col, "cluster_overlap"))

    plasmid_col = find_column(rows, PLASMID_ALIASES)
    if plasmid_col:
        out.extend(_overlap(rows, plasmid_col, "plasmid_overlap"))

    return out


def summary_fieldnames() -> list[str]:
    return ["section", "key", "count", "l1_set", "mixed"]


def _counts(rows: list[dict[str, str]], field: str, section: str) -> list[dict[str, str]]:
    counter: Counter[str] = Counter((row.get(field) or "") for row in rows)
    result = []
    for key, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
        result.append(
            {
                "section": section,
                "key": key,
                "count": str(count),
                "l1_set": "",
                "mixed": "",
            }
        )
    return result


def _source_evidence(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Split ok/review so empty NCBI source is not confused with a lexicon miss."""
    counter: Counter[str] = Counter(_evidence_key(row) for row in rows)
    keys = list(SOURCE_EVIDENCE_CORE)
    extras = sorted(key for key in counter if key not in SOURCE_EVIDENCE_CORE)
    keys.extend(extras)
    return [
        {
            "section": "source_evidence",
            "key": key,
            "count": str(counter.get(key, 0)),
            "l1_set": "",
            "mixed": "",
        }
        for key in keys
    ]


def _evidence_key(row: dict[str, str]) -> str:
    source = get_field(row, ISOLATION_SOURCE_ALIASES)
    status = (row.get("classify_status") or "").strip()
    reason = (row.get("classify_reason") or "").strip()
    resolve_status = (row.get("resolve_status") or "").strip()

    if "conflict" in reason:
        return "conflict"
    if reason == "unmatched":
        return "unmatched"
    if is_blank(source):
        if status == "ok" or resolve_status == "ok":
            return "ok_empty_source"
        if reason == "missing_source":
            return "missing_source"
        return "empty_source"
    if status == "ok":
        return "ok_with_source"
    if reason:
        return reason
    return "review"


def _overlap(rows: list[dict[str, str]], column: str, section: str) -> list[dict[str, str]]:
    groups: dict[str, set[str]] = {}
    for row in rows:
        key = (row.get(column) or "").strip()
        if not key:
            continue
        groups.setdefault(key, set()).add(row.get("l1") or "")
    result = []
    for key, l1s in sorted(groups.items()):
        cleaned = sorted(item for item in l1s if item)
        result.append(
            {
                "section": section,
                "key": key,
                "count": str(len(cleaned)),
                "l1_set": ",".join(cleaned),
                "mixed": "yes" if len(cleaned) > 1 else "no",
            }
        )
    return result
