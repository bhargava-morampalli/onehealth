"""Count classifications. Cluster/plasmid overlap only if those columns exist."""

from __future__ import annotations

from collections import Counter

from onehealth.constants import CLUSTER_ALIASES, ISOLATION_SOURCE_ALIASES, PLASMID_ALIASES
from onehealth.io_tsv import find_column, get_field, isolate_key
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
    out.extend(_dataset(rows))
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


def _row(section: str, key: str, count: int, l1_set: str = "", mixed: str = "") -> dict[str, str]:
    return {"section": section, "key": key, "count": str(count), "l1_set": l1_set, "mixed": mixed}


def _dataset(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Counts below are per input row. Rows that share a BioSample (or other ID) are the same isolate.

    A row with no identifier at all cannot be matched to another, so it counts as its own isolate.
    """
    keys = [isolate_key(row) for row in rows]
    isolates = len({key for key in keys if key}) + keys.count("")
    return [
        _row("dataset", "rows", len(rows)),
        _row("dataset", "unique_isolates", isolates),
        _row("dataset", "rows_sharing_a_biosample", len(rows) - isolates),
    ]


def _counts(rows: list[dict[str, str]], field: str, section: str) -> list[dict[str, str]]:
    counter: Counter[str] = Counter((row.get(field) or "") for row in rows)
    return [_row(section, key, count) for key, count in sorted(counter.items(), key=lambda item: (-item[1], item[0]))]


def _source_evidence(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """ok rows split by whether isolation_source was filled; review rows by reason."""
    counter: Counter[str] = Counter(_evidence_key(row) for row in rows)
    keys = list(SOURCE_EVIDENCE_CORE) + sorted(key for key in counter if key not in SOURCE_EVIDENCE_CORE)
    return [_row("source_evidence", key, counter.get(key, 0)) for key in keys]


def _evidence_key(row: dict[str, str]) -> str:
    status = (row.get("classify_status") or "").strip()
    reason = (row.get("classify_reason") or "").strip()
    if status == "ok":
        return "ok_empty_source" if is_blank(get_field(row, ISOLATION_SOURCE_ALIASES)) else "ok_with_source"
    if "conflict" in reason:
        return "conflict"
    return reason or "review"


def _overlap(rows: list[dict[str, str]], column: str, section: str) -> list[dict[str, str]]:
    """Per cluster: isolates (count) and the classified L1s present.

    Unclassified is not a source, so it never makes a cluster mixed.
    Shared cluster or plasmid is not evidence of transmission on its own.
    """
    members: dict[str, list[str]] = {}
    for row in rows:
        key = (row.get(column) or "").strip()
        if key:
            members.setdefault(key, []).append(row.get("l1") or "")
    result = []
    for key, l1s in sorted(members.items()):
        classified = sorted({l1 for l1 in l1s if l1 and l1 != "Unclassified"})
        result.append(_row(section, key, len(l1s), ",".join(classified), "yes" if len(classified) > 1 else "no"))
    return result
