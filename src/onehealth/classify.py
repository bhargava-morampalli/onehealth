"""Classify isolation-source metadata. Gaps and conflicts go to review."""

from __future__ import annotations

from dataclasses import dataclass, asdict

from onehealth.constants import (
    CLASSIFY_EXTRA_FIELDS,
    HEALTHCARE_NO,
    HEALTHCARE_YES,
    derive_paper_rollup,
)
from onehealth.io_tsv import metadata_view
from onehealth.lexicon import Lexicon, RuleHit, best_hits, match_rules
from onehealth.overrides import OverrideHit, OverrideRow, lookup_override
from onehealth.text import is_blank, is_unknown_value


@dataclass(frozen=True)
class Classification:
    l1: str
    healthcare_environment: str
    paper_rollup: str
    classify_status: str
    classify_reason: str
    lexicon_version: str
    override_applied: str
    lexicon_suggestion: str

    def as_row(self) -> dict[str, str]:
        return asdict(self)


def classify_record(
    record: dict[str, str],
    lexicon: Lexicon,
    overrides: dict[str, list[OverrideRow]] | None = None,
) -> Classification:
    view = metadata_view(record)
    hits = match_rules(view, lexicon)
    suggestion = _from_hits(hits, view, lexicon.version)

    override = lookup_override(record, overrides or {})
    if override is not None:
        return _apply_override(override, suggestion, lexicon.version)

    if _missing_source(view):
        return _review(
            "Unclassified",
            HEALTHCARE_NO,
            "missing_source",
            lexicon.version,
            suggestion,
            override_applied="no",
        )

    if _unknown_blocks_environment(view):
        return _review(
            "Unclassified",
            HEALTHCARE_NO,
            "unknown_source",
            lexicon.version,
            suggestion,
            override_applied="no",
        )

    winners = best_hits(hits)
    if not winners:
        return _review(
            "Unclassified",
            HEALTHCARE_NO,
            "unmatched",
            lexicon.version,
            suggestion,
            override_applied="no",
        )

    l1_set = {hit.rule.l1 for hit in winners}
    hc_set = {hit.rule.healthcare_environment for hit in winners}
    if len(l1_set) > 1 or len(hc_set) > 1:
        if _allowed_host_env_split(winners):
            env_hit = next(hit for hit in winners if hit.rule.healthcare_environment == HEALTHCARE_YES)
            return _ok(env_hit.rule.l1, env_hit.rule.healthcare_environment, env_hit.rule.rule_id, lexicon.version)
        return _review(
            "Unclassified",
            HEALTHCARE_NO,
            "rule_conflict:" + ",".join(sorted(hit.rule.rule_id for hit in winners)),
            lexicon.version,
            suggestion,
            override_applied="no",
        )

    winner = winners[0]
    # Host vs source disagreement at different priorities is resolved by priority.
    # Same-priority leftovers already handled. Check complementary host/source conflict
    # when both sides produced different L1 at any priority and source is not healthcare-env.
    if _host_source_conflict(hits):
        return _review(
            "Unclassified",
            HEALTHCARE_NO,
            "host_source_conflict",
            lexicon.version,
            suggestion,
            override_applied="no",
        )

    return _ok(winner.rule.l1, winner.rule.healthcare_environment, winner.rule.rule_id, lexicon.version)


def classify_rows(
    rows: list[dict[str, str]],
    lexicon: Lexicon,
    overrides: dict[str, list[OverrideRow]] | None = None,
) -> list[dict[str, str]]:
    classified: list[dict[str, str]] = []
    for row in rows:
        result = classify_record(row, lexicon, overrides)
        merged = dict(row)
        merged.update(result.as_row())
        classified.append(merged)
    return classified


def classify_fieldnames(rows: list[dict[str, str]]) -> list[str]:
    extra = list(CLASSIFY_EXTRA_FIELDS)
    if not rows:
        return extra
    names: list[str] = []
    for key in rows[0]:
        if key not in extra:
            names.append(key)
    names.extend(extra)
    return names


def _apply_override(override: OverrideHit, suggestion: str, version: str) -> Classification:
    if override.status == "applied" and override.row is not None:
        row = override.row
        healthcare = row.healthcare_environment or HEALTHCARE_NO
        paper = row.paper_rollup or derive_paper_rollup(row.l1, healthcare)
        return Classification(
            l1=row.l1,
            healthcare_environment=healthcare,
            paper_rollup=paper,
            classify_status="ok",
            classify_reason="override",
            lexicon_version=version,
            override_applied="yes",
            lexicon_suggestion=suggestion,
        )
    reason = override.reason
    return _review(
        "Unclassified",
        HEALTHCARE_NO,
        reason,
        version,
        suggestion,
        override_applied="no",
    )


def _from_hits(hits: list[RuleHit], view: dict[str, str], version: str) -> str:
    if _missing_source(view) or _unknown_blocks_environment(view):
        return "Unclassified"
    winners = best_hits(hits)
    if len(winners) == 1:
        return winners[0].rule.l1
    if winners and len({hit.rule.l1 for hit in winners}) == 1:
        return winners[0].rule.l1
    return ""


def _missing_source(view: dict[str, str]) -> bool:
    return is_blank(view.get("isolation_source")) and is_blank(view.get("host")) and is_blank(
        view.get("sample_type")
    )


def _unknown_blocks_environment(view: dict[str, str]) -> bool:
    """Unknown isolation_source stays Unclassified; never map to Environmental."""
    source = view.get("isolation_source", "")
    host = view.get("host", "")
    if not is_unknown_value(source):
        return False
    if not is_blank(host) and not is_unknown_value(host):
        return False
    return True


def _allowed_host_env_split(winners: list[RuleHit]) -> bool:
    l1s = {hit.rule.l1 for hit in winners}
    healthcare = {hit.rule.healthcare_environment for hit in winners}
    if HEALTHCARE_YES in healthcare and l1s <= {"Human", "Environmental"}:
        return True
    return False


def _host_source_conflict(hits: list[RuleHit]) -> bool:
    source_hits = best_hits([hit for hit in hits if hit.field == "isolation_source"])
    host_hits = best_hits([hit for hit in hits if hit.field == "host"])
    if not source_hits or not host_hits:
        return False
    source_l1 = {hit.rule.l1 for hit in source_hits}
    host_l1 = {hit.rule.l1 for hit in host_hits}
    if source_l1 == host_l1:
        return False
    if any(hit.rule.healthcare_environment == HEALTHCARE_YES for hit in source_hits) and host_l1 == {"Human"}:
        return False
    # Specimen from an animal host is Animal, not a conflict with a human-anatomy word
    # only when the source rule is a generic human specimen AND host is Animal.
    # Those still need review if both are confident and disagree (e.g. human stool vs chicken).
    return True


def _ok(l1: str, healthcare: str, reason: str, version: str) -> Classification:
    return Classification(
        l1=l1,
        healthcare_environment=healthcare,
        paper_rollup=derive_paper_rollup(l1, healthcare),
        classify_status="ok",
        classify_reason=reason,
        lexicon_version=version,
        override_applied="no",
        lexicon_suggestion=l1,
    )


def _review(
    l1: str,
    healthcare: str,
    reason: str,
    version: str,
    suggestion: str,
    override_applied: str,
) -> Classification:
    return Classification(
        l1=l1,
        healthcare_environment=healthcare,
        paper_rollup=derive_paper_rollup(l1, healthcare),
        classify_status="review",
        classify_reason=reason,
        lexicon_version=version,
        override_applied=override_applied,
        lexicon_suggestion=suggestion,
    )
