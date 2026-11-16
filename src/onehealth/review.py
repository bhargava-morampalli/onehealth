"""The review queue by distinct source text, and a curator's decisions on it.

`classify` lists every row it would not guess in the review queue. Many of those rows
say the same thing: in 159,045 NCBI Pathogen Detection records, the 2,626 review rows
with any text hold 694 distinct texts, and the 100 most common cover 69% of them. Here
the queue is grouped, one line per distinct text with the number of rows that carry
it, and with empty columns for a decision. Filled in, the same table is what
`classify --decisions` reads: each decision applies to every row with that text.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from onehealth.constants import HEALTHCARE_NO, HEALTHCARE_YES, L1_VALUES, PAPER_ROLLUP_VALUES
from onehealth.errors import InputError
from onehealth.io_tsv import first_present, get_field, isolate_key, metadata_view, read_tsv
from onehealth.overrides import OverrideHit, OverrideRow
from onehealth.text import environmental_package, is_unknown_value, known_part, normalize_text, tokenize

# The fields the lexicon reads. Two rows with the same words in each are classified alike.
EVIDENCE_FIELDS = ("isolation_source", "host", "sample_type", "host_disease")
# A specimen with no host is classified differently in a record filed as environmental,
# food or other, so that is part of what a decision is about.
ENVIRONMENTAL_PACKAGE = "environmental_package"

GROUP_FIELDS = (
    "rows", "classify_reason", "lexicon_suggestion", *EVIDENCE_FIELDS, ENVIRONMENTAL_PACKAGE, "example_id",
    "l1", "healthcare_environment", "note",
)
_L1_COLUMNS = ("l1", "source_l1", "classification")

# (isolation_source, host, sample_type, host_disease, "yes" or "no" for an environmental package)
DecisionKey = tuple[str, str, str, str, str]


def _words(value: str | None) -> str:
    """A field as the lexicon sees it: its words, or nothing for a missing-value word.

    Rules match on words, so "Chicken Breast", "chicken-breast" and "chicken_breast" are
    one text. "unknown animal" still says something ("animal"), so it is kept.
    """
    if is_unknown_value(value) and not known_part(value):
        return ""
    return " ".join(tokenize(value))


def evidence_key(view: dict[str, str]) -> DecisionKey:
    """What a row says, for grouping and for matching decisions. `view` is io_tsv.metadata_view of the row."""
    package = HEALTHCARE_YES if environmental_package(view.get("package")) else HEALTHCARE_NO
    words = tuple(_words(view.get(name)) for name in EVIDENCE_FIELDS)
    return (words[0], words[1], words[2], words[3], package)


def review_groups(classified: list[dict[str, str]]) -> list[dict[str, str]]:
    """One row per distinct source text and reason among the review rows, largest first.

    The text shown is that of the first row in the group, as it was typed. The decision
    columns (l1, healthcare_environment, note) are left empty for the curator.
    """
    groups: dict[tuple[DecisionKey, str], dict[str, Any]] = {}  # "rows" is a count until the end
    for row in classified:
        if row.get("classify_status") != "review":
            continue
        view = metadata_view(row)
        key = (evidence_key(view), row.get("classify_reason") or "")
        group = groups.get(key)
        if group is None:
            group = {
                "rows": 0,
                "classify_reason": key[1],
                "lexicon_suggestion": row.get("lexicon_suggestion") or "",
                # As typed, except that a missing-value word ("NULL", "not collected") is shown as empty.
                **{name: view.get(name, "") if words else "" for name, words in zip(EVIDENCE_FIELDS, key[0])},
                ENVIRONMENTAL_PACKAGE: HEALTHCARE_YES if key[0][4] == HEALTHCARE_YES else "",
                "example_id": isolate_key(row),
                "l1": "",
                "healthcare_environment": "",
                "note": "",
            }
            groups[key] = group
        group["rows"] += 1
    ordered = sorted(groups.items(), key=lambda item: (-item[1]["rows"], item[0]))
    return [{**group, "rows": str(group["rows"])} for _, group in ordered]


def load_decisions(path: Path | None) -> dict[DecisionKey, OverrideRow]:
    """Decisions from a review-groups table whose l1 column has been filled in.

    A row with an empty l1 has not been decided and is skipped. Anything wrong with a
    decision stops the run: applied to hundreds of rows at once, a mistyped class or two
    answers for one text must not pass.
    """
    if path is None:
        return {}
    rows = read_tsv(path)
    if not rows:
        return {}
    if first_present(rows[0], _L1_COLUMNS) is None:
        raise InputError(f"{path}: no l1 column ({', '.join(_L1_COLUMNS)}) to hold the decisions")
    if all(not metadata_view(row)[name] for row in rows for name in EVIDENCE_FIELDS):
        raise InputError(f"{path}: no source text to match decisions on (columns {', '.join(EVIDENCE_FIELDS)})")
    has_package = first_present(rows[0], (ENVIRONMENTAL_PACKAGE,)) is not None
    classes = {name.lower(): name for name in L1_VALUES}
    rollups = {name.lower(): name for name in PAPER_ROLLUP_VALUES}
    decisions: dict[DecisionKey, OverrideRow] = {}
    for number, row in enumerate(rows, start=1):
        label = get_field(row, _L1_COLUMNS)
        if not label:
            continue
        where = f"{path}: decision {number}"
        view = metadata_view(row)
        key = evidence_key(view)
        shown = " | ".join(view[name] for name in EVIDENCE_FIELDS if view[name])
        if not any(key[:4]):
            raise InputError(f"{where} has a class ({label}) but no isolation_source, host, sample_type or "
                             "host_disease, so it would apply to every row without them")
        if label.lower() not in classes:
            raise InputError(f"{where} ({shown}): {label!r} is not a class; use one of {', '.join(L1_VALUES)}")
        healthcare = normalize_text(get_field(row, ("healthcare_environment", "healthcare_env")))
        if healthcare not in {"", HEALTHCARE_YES, HEALTHCARE_NO}:
            raise InputError(f"{where} ({shown}): healthcare_environment must be yes or no, not {healthcare!r}")
        rollup = get_field(row, ("paper_rollup", "paper"))
        if rollup and rollup.lower() not in rollups:
            raise InputError(f"{where} ({shown}): {rollup!r} is not a paper rollup; use one of {', '.join(PAPER_ROLLUP_VALUES)}")
        decision = OverrideRow(
            key=shown,
            l1=classes[label.lower()],
            healthcare_environment=healthcare,
            paper_rollup=rollups[rollup.lower()] if rollup else "",
            note=get_field(row, ("note", "reason", "comment")),
            complete=True,
        )
        if has_package:
            package = normalize_text(get_field(row, (ENVIRONMENTAL_PACKAGE,)))
            if package not in {"", HEALTHCARE_YES, HEALTHCARE_NO}:
                raise InputError(f"{where} ({shown}): {ENVIRONMENTAL_PACKAGE} must be yes or empty, not {package!r}")
            keys = [(*key[:4], HEALTHCARE_YES if package == HEALTHCARE_YES else HEALTHCARE_NO)]
        else:
            # A table written by hand, without the column: the decision is about the text alone.
            keys = [(*key[:4], HEALTHCARE_YES), (*key[:4], HEALTHCARE_NO)]
        for full_key in keys:
            earlier = decisions.get(full_key)
            if earlier is not None and _answer(earlier) != _answer(decision):
                raise InputError(f"{where} ({shown}) is {decision.l1}, but an earlier decision for the same text is "
                                 f"{earlier.l1}; keep one decision per text")
            decisions[full_key] = decision
    return decisions


def _answer(decision: OverrideRow) -> tuple[str, str, str]:
    return decision.l1, decision.healthcare_environment or HEALTHCARE_NO, decision.paper_rollup


def lookup_decision(view: dict[str, str], decisions: dict[DecisionKey, OverrideRow]) -> OverrideHit | None:
    decision = decisions.get(evidence_key(view))
    if decision is None:
        return None
    return OverrideHit(status="applied", reason="decision", row=decision)


def unused_decisions(rows: list[dict[str, str]], decisions: dict[DecisionKey, OverrideRow]) -> list[str]:
    """Texts with a decision that no row carries, so the decision was never applied."""
    if not decisions:
        return []
    used = {evidence_key(metadata_view(row)) for row in rows} & decisions.keys()
    applied = {id(decisions[key]) for key in used}
    return sorted({decision.key for decision in decisions.values() if id(decision) not in applied})
