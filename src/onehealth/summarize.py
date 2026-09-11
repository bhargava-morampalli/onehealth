"""Count classifications. Cluster/plasmid overlap only if those columns exist."""

from __future__ import annotations

from collections import Counter

from onehealth.constants import CLUSTER_ALIASES, PLASMID_ALIASES
from onehealth.io_tsv import find_column


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    out.extend(_counts(rows, "l1", "l1"))
    out.extend(_counts(rows, "paper_rollup", "paper_rollup"))
    out.extend(_counts(rows, "healthcare_environment", "healthcare_environment"))
    out.extend(_counts(rows, "classify_status", "classify_status"))
    out.extend(_counts(rows, "classify_reason", "classify_reason"))

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
