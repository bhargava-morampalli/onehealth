"""Classify isolation-source metadata. Gaps and conflicts go to review."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from onehealth.constants import (
    CLASSIFY_EXTRA_FIELDS,
    HEALTHCARE_NO,
    HEALTHCARE_YES,
    derive_paper_rollup,
)
from onehealth.io_tsv import metadata_view
from onehealth.lexicon import Lexicon, LexiconRule, RuleHit, best_hits, match_rules
from onehealth.overrides import OverrideHit, OverrideRow, lookup_override
from onehealth.text import is_blank, is_taxon_only_source, is_unknown_value, known_part, tokenize

# Unclassified rows carry no healthcare-environment call.
HEALTHCARE_UNKNOWN = ""


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
    auto = _classify_with_lexicon(metadata_view(record), lexicon)
    override = lookup_override(record, overrides or {})
    if override is None:
        return auto
    return _apply_override(override, auto.lexicon_suggestion, lexicon.version)


def classify_rows(
    rows: list[dict[str, str]],
    lexicon: Lexicon,
    overrides: dict[str, list[OverrideRow]] | None = None,
) -> list[dict[str, str]]:
    classified: list[dict[str, str]] = []
    for row in rows:
        merged = dict(row)
        merged.update(classify_record(row, lexicon, overrides).as_row())
        classified.append(merged)
    return classified


def classify_fieldnames(rows: list[dict[str, str]]) -> list[str]:
    extra = list(CLASSIFY_EXTRA_FIELDS)
    if not rows:
        return extra
    return [key for key in rows[0] if key not in extra] + extra


def _classify_with_lexicon(view: dict[str, str], lexicon: Lexicon) -> Classification:
    version = lexicon.version
    # "unknown animal" leaves the species open, not the class: match on "animal".
    # Other fields keep all their words anyway, so a host of "unknown animal" matches too.
    source = view.get("isolation_source", "")
    matched = {**view, "isolation_source": known_part(source) or source}
    source = matched["isolation_source"]
    hits = match_rules(matched, lexicon)

    # Any other unknown marker ("unknown", "missing: lab stock") or a bare taxon name
    # ("Escherichia coli", "Bos taurus") in isolation_source is not source evidence.
    # Match again as if the field were blank: a rule stops at the first of its fields
    # that matches, so its other fields (source_type, host) must get their turn.
    if any(hit.field == "isolation_source" for hit in hits) and (is_unknown_value(source) or _taxon_only(source, hits)):
        hits = match_rules({**matched, "isolation_source": ""}, lexicon)

    winners = best_hits(hits)
    if not winners:
        return _review(_no_evidence_reason(view), version)

    l1s = {hit.rule.l1 for hit in winners}
    if len(l1s) > 1 or len({hit.rule.healthcare_environment for hit in winners}) > 1:
        reason = "rule_conflict:" + ",".join(sorted(hit.rule.rule_id for hit in winners))
        return _review(reason, version, l1s)

    rule = winners[0].rule
    if rule.l1 == "Unclassified":
        return _review(rule.rule_id, version)
    if all(hit.rule.host_dependent for hit in winners):
        return _classify_specimen(rule, hits, view, version)

    conflict = _host_source_conflict(hits)
    if conflict:
        return _review("host_source_conflict", version, conflict)
    return _ok(rule.l1, rule.healthcare_environment, rule.rule_id, version)


def _taxon_only(source: str, hits: list[RuleHit]) -> bool:
    """isolation_source is a bare scientific name, not a description of the sample.

    The word-shape test in is_taxon_only_source cannot tell "Bos taurus" from
    "bedpan washer", so a lexicon rule written for source text overrules it.
    Rules that also read `host` list host species ("bos taurus", "homo sapiens"),
    so a match from one of those is still treated as a taxon name.
    """
    if not is_taxon_only_source(source):
        return False
    return not any(hit.field == "isolation_source" and "host" not in hit.rule.fields for hit in hits)


def _classify_specimen(rule: LexiconRule, hits: list[RuleHit], view: dict[str, str], version: str) -> Classification:
    """Blood, urine, faeces, wound, "clinical": the host decides Human vs Animal.

    Anatomical entities are species-neutral (Uberon), so "pig faeces" or a
    Bos taurus host with "faeces" is Animal, not Human.
    """
    deciders = [hit for hit in hits if not hit.rule.host_dependent and hit.rule.l1 in {"Human", "Animal"}]
    candidates = {hit.rule.l1 for hit in deciders}
    # Another field that names a class no host can have (source_type "food" or
    # "environmental") contradicts the specimen. Words in the specimen's own text do
    # not: "stool from a food handler" describes the person, not the sample.
    specimen_fields = {hit.field for hit in hits if hit.rule.host_dependent}
    non_host = {
        hit.rule.l1
        for hit in hits
        if not hit.rule.host_dependent
        and hit.field not in specimen_fields
        and hit.rule.l1 not in {"Human", "Animal", "Unclassified"}
    }
    if non_host:
        return _review("host_source_conflict", version, (candidates or {rule.l1}) | non_host)
    if len(candidates) > 1:
        return _review("host_source_conflict", version, candidates)
    if deciders:
        decider = best_hits(deciders)[0]
        return _ok(decider.rule.l1, HEALTHCARE_NO, f"{rule.rule_id}+{decider.rule.rule_id}", version)
    if not is_unknown_value(view.get("host")):
        # A host the lexicon cannot place (plant, unlisted species): do not assume Human.
        return _review("host_unrecognised", version, {rule.l1})
    if _environmental_package(view.get("package")):
        # No host, and the submitter filed the sample as environmental, food or other, not
        # as clinical or host-associated: "feces" or "wound" here is often an animal's.
        return _review("no_host_env_package", version, {rule.l1})
    return _ok(rule.l1, rule.healthcare_environment, f"{rule.rule_id}+no_host", version)


def _environmental_package(value: str | None) -> bool:
    """The record was submitted under NCBI's package for non-clinical pathogen samples.

    True for the BioSample package "Pathogen.env.1.0", its display name "Pathogen:
    environmental/food/other", and NCBI Pathogen Detection's epi_type "environmental/other".
    """
    tokens = tokenize(value)
    return tokens[:2] in (["pathogen", "env"], ["pathogen", "environmental"]) or tokens == ["environmental", "other"]


def _host_source_conflict(hits: list[RuleHit]) -> set[str]:
    """L1s in conflict between the host field and the source fields, or empty.

    A rule that sends its matches to review (l1 Unclassified) names no class, so it
    cannot conflict with the host; it decides only when it is the best hit.
    """
    classes = [hit for hit in hits if hit.rule.l1 != "Unclassified"]
    source_hits = best_hits([hit for hit in classes if hit.field != "host"])
    host_hits = best_hits([hit for hit in classes if hit.field == "host"])
    source_l1 = {hit.rule.l1 for hit in source_hits}
    host_l1 = {hit.rule.l1 for hit in host_hits}
    if not source_l1 or not host_l1 or source_l1 == host_l1:
        return set()
    # Submitters often record host = Homo sapiens for hospital sink or surface swabs.
    if host_l1 == {"Human"} and any(hit.rule.healthcare_environment == HEALTHCARE_YES for hit in source_hits):
        return set()
    return source_l1 | host_l1


def _no_evidence_reason(view: dict[str, str]) -> str:
    source = view.get("isolation_source", "")
    # A missing-value term first: "not collected" also has the word shape of a binomial.
    if not is_blank(source) and is_unknown_value(source):
        return "unknown_source"
    if is_taxon_only_source(source):
        return "taxon_only_source"
    if is_blank(source) and is_unknown_value(view.get("host")):
        # NCBI sample_type ("whole organism", "pure culture") is material, not source.
        return "missing_source"
    return "unmatched"


def _apply_override(override: OverrideHit, suggestion: str, version: str) -> Classification:
    if override.status != "applied" or override.row is None:
        return _review(override.reason, version, {suggestion} if suggestion else set())
    row = override.row
    default_hc = HEALTHCARE_UNKNOWN if row.l1 == "Unclassified" else HEALTHCARE_NO
    healthcare = row.healthcare_environment or default_hc
    return Classification(
        l1=row.l1,
        healthcare_environment=healthcare,
        paper_rollup=row.paper_rollup or derive_paper_rollup(row.l1, healthcare),
        classify_status="ok",
        classify_reason="override",
        lexicon_version=version,
        override_applied="yes",
        lexicon_suggestion=suggestion,
    )


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


def _review(reason: str, version: str, candidates: set[str] | None = None) -> Classification:
    """Unclassified. lexicon_suggestion lists the competing L1s, if any, for the curator."""
    return Classification(
        l1="Unclassified",
        healthcare_environment=HEALTHCARE_UNKNOWN,
        paper_rollup="Unclassified",
        classify_status="review",
        classify_reason=reason,
        lexicon_version=version,
        override_applied="no",
        lexicon_suggestion="|".join(sorted(candidates or ())),
    )
