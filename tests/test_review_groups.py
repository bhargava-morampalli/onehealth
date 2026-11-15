"""The review queue by distinct source text, and decisions on it (review.py)."""

from pathlib import Path

import pytest

from onehealth import cli
from onehealth.classify import classify_rows
from onehealth.io_tsv import read_tsv, write_tsv
from onehealth.overrides import OverrideRow
from onehealth.review import GROUP_FIELDS, load_decisions, review_groups, unused_decisions


def _rows(*sources, **common):
    return [{"id": f"ISO{i}", "isolation_source": source, **common} for i, source in enumerate(sources)]


def _decisions(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "decisions.tsv"
    path.write_text(text, encoding="utf-8")
    return path


# --- grouping --------------------------------------------------------------------------------------


def test_review_rows_are_grouped_by_what_they_say(lexicon):
    rows = _rows("Bed Rails", "bed rails", "bed_rails", "bed-rails!", "ECG monitors", "blood", "", "not collected")
    groups = review_groups(classify_rows(rows, lexicon))
    assert [(g["rows"], g["isolation_source"], g["classify_reason"]) for g in groups] == [
        ("4", "Bed Rails", "taxon_only_source"),  # one text however it is typed; shown as first typed
        ("1", "", "missing_source"),  # an empty cell
        ("1", "", "unknown_source"),  # a missing-value word is shown as nothing too; its reason differs
        ("1", "ECG monitors", "taxon_only_source"),
    ]
    assert groups[0]["example_id"] == "ISO0"
    assert all(g["l1"] == g["healthcare_environment"] == g["note"] == "" for g in groups)  # left for the curator
    assert all(list(g) == list(GROUP_FIELDS) for g in groups)


def test_groups_are_split_by_host_reason_and_environmental_package(lexicon):
    rows = [
        {"id": "A", "isolation_source": "feces", "host": "Zea mays"},
        {"id": "B", "isolation_source": "feces", "host": "Zea mays"},
        {"id": "C", "isolation_source": "feces", "host": "", "package": "Pathogen.env.1.0"},
        {"id": "D", "isolation_source": "feces", "host": "", "epi_type": "environmental/other"},
        {"id": "E", "isolation_source": "feces", "host": ""},  # Human by the no-host default: not in review
        {"id": "F", "isolation_source": "river water", "host": "Homo sapiens"},
    ]
    groups = review_groups(classify_rows(rows, lexicon))
    assert [(g["rows"], g["isolation_source"], g["host"], g["environmental_package"], g["classify_reason"], g["lexicon_suggestion"])
            for g in groups] == [
        ("2", "feces", "", "yes", "no_host_env_package", "Human"),
        ("2", "feces", "Zea mays", "", "host_unrecognised", "Human"),
        ("1", "river water", "Homo sapiens", "", "host_source_conflict", "Environmental|Human"),
    ]


# --- decisions -------------------------------------------------------------------------------------


def test_a_decision_applies_to_every_row_with_that_text(lexicon, tmp_path: Path):
    rows = _rows("Bed Rails", "bed rails", "ECG monitors", "blood")
    decisions = load_decisions(_decisions(tmp_path, (
        "rows\tclassify_reason\tisolation_source\thost\tsample_type\thost_disease\tenvironmental_package\tl1\thealthcare_environment\tnote\n"
        "2\ttaxon_only_source\tBed Rails\t\t\t\t\tenvironmental\tyes\thospital ward\n"
        "1\ttaxon_only_source\tECG monitors\t\t\t\t\t\t\t\n"  # not decided yet: skipped
    )))
    out = classify_rows(rows, lexicon, decisions=decisions)
    assert [(r["l1"], r["healthcare_environment"], r["paper_rollup"], r["classify_status"], r["classify_reason"], r["override_applied"])
            for r in out] == [
        ("Environmental", "yes", "Environmental", "ok", "decision", "yes"),
        ("Environmental", "yes", "Environmental", "ok", "decision", "yes"),
        ("Unclassified", "", "Unclassified", "review", "taxon_only_source", "no"),
        ("Human", "no", "Human", "ok", "specimens+no_host", "no"),
    ]
    assert unused_decisions(rows, decisions) == []


def test_a_decision_also_corrects_a_text_the_lexicon_classifies(lexicon, tmp_path: Path):
    """A food-safety study can say that a bare "chicken" is meat. The lexicon's answer is kept as the suggestion."""
    decisions = load_decisions(_decisions(tmp_path, "isolation_source\tl1\nchicken\tFood\n"))
    out = classify_rows(_rows("Chicken", "chicken breast", "chicken", host=""), lexicon, decisions=decisions)
    assert [(r["l1"], r["classify_reason"], r["lexicon_suggestion"]) for r in out] == [
        ("Food", "decision", "Animal"), ("Food", "food-items", "Food"), ("Food", "decision", "Animal")]
    # The host is part of the text: "chicken" from a chicken is another text, and the decision does not reach it.
    with_host = classify_rows(_rows("chicken", host="Gallus gallus"), lexicon, decisions=decisions)
    assert with_host[0]["classify_reason"] == "animal-hosts"


def test_an_override_for_the_isolate_wins_over_a_decision_for_its_text(lexicon, tmp_path: Path):
    decisions = load_decisions(_decisions(tmp_path, "isolation_source\tl1\nbed rails\tEnvironmental\n"))
    overrides = {"ISO1": [OverrideRow(key="ISO1", l1="Laboratory", healthcare_environment="", paper_rollup="", note="", complete=True)]}
    out = classify_rows(_rows("bed rails", "bed rails"), lexicon, overrides, decisions)
    assert [(r["l1"], r["classify_reason"]) for r in out] == [("Environmental", "decision"), ("Laboratory", "override")]


def test_environmental_package_is_part_of_a_decision_when_the_table_has_the_column(lexicon, tmp_path: Path):
    rows = [
        {"id": "A", "isolation_source": "feces", "package": "Pathogen.env.1.0"},
        {"id": "B", "isolation_source": "feces", "package": "Pathogen.cl.1.0"},
    ]
    from_groups = load_decisions(_decisions(tmp_path, "isolation_source\tenvironmental_package\tl1\nfeces\tyes\tAnimal\n"))
    assert [(r["l1"], r["classify_reason"]) for r in classify_rows(rows, lexicon, decisions=from_groups)] == [
        ("Animal", "decision"), ("Human", "specimens+no_host")]
    # A table written by hand, without the column, is about the text alone.
    by_hand = load_decisions(_decisions(tmp_path, "isolation_source\tl1\nfeces\tAnimal\n"))
    assert [r["l1"] for r in classify_rows(rows, lexicon, decisions=by_hand)] == ["Animal", "Animal"]


@pytest.mark.parametrize(
    "table, message",
    [
        ("isolation_source\tl1\nbed rails\tHospital\n", "'Hospital' is not a class"),
        ("isolation_source\tl1\thealthcare_environment\nbed rails\tEnvironmental\tmaybe\n", "healthcare_environment must be yes or no"),
        ("isolation_source\tl1\tpaper_rollup\nbed rails\tEnvironmental\tFarm\n", "'Farm' is not a paper rollup"),
        ("isolation_source\tl1\nbed rails\tEnvironmental\nBed-Rails\tHuman\n", "an earlier decision for the same text is Environmental"),
        ("isolation_source\tl1\nbed rails\tEnvironmental\nnot collected\tHuman\n", "would apply to every row without them"),
        ("isolation_source\thost\nbed rails\t\n", "no l1 column"),
        ("id\tl1\nISO1\tHuman\n", "no source text to match decisions on"),
        ("isolation_source\tenvironmental_package\tl1\nfeces\tperhaps\tAnimal\n", "environmental_package must be yes or empty"),
    ],
)
def test_a_wrong_decision_stops_the_run(tmp_path: Path, table, message):
    with pytest.raises(ValueError, match=message):
        load_decisions(_decisions(tmp_path, table))


def test_the_same_decision_twice_is_fine_and_class_names_are_case_insensitive(tmp_path: Path):
    decisions = load_decisions(_decisions(tmp_path, "isolation_source\tl1\nbed rails\tenvironmental\nBed Rails\tEnvironmental\n"))
    assert {d.l1 for d in decisions.values()} == {"Environmental"}
    assert load_decisions(None) == {} and load_decisions(_decisions(tmp_path, "isolation_source\tl1\n")) == {}


# --- the command -----------------------------------------------------------------------------------


def test_classify_writes_the_groups_and_reads_them_back_as_decisions(tmp_path: Path, capsys):
    table = tmp_path / "metadata.tsv"
    rows = _rows("Bed Rails", "bed rails", '=HYPERLINK("http://x","swab")', '=hyperlink("http://x","swab")', "blood", host="")
    write_tsv(table, rows, ["id", "isolation_source", "host"])
    out = tmp_path / "classified.tsv"
    assert cli.main(["classify", str(table), "-o", str(out)]) == 0
    assert "the review rows hold 2 distinct source texts" in capsys.readouterr().err
    groups_path = tmp_path / "classified.review_groups.tsv"
    groups = read_tsv(groups_path)
    assert [(g["rows"], g["classify_reason"]) for g in groups] == [("2", "taxon_only_source"), ("2", "unmatched")]

    # The curator fills in l1 and passes the same table back. The formula text was written
    # with a leading apostrophe; it still matches, because texts are compared by their words.
    for group in groups:
        group["l1"] = "Environmental"
    decided = tmp_path / "decided.tsv"
    write_tsv(decided, groups, list(groups[0]))
    again = tmp_path / "again.tsv"
    assert cli.main(["classify", str(table), "-o", str(again), "--decisions", str(decided),
                     "--review-groups", str(tmp_path / "left.tsv")]) == 0
    err = capsys.readouterr().err
    assert "classified 4 rows" in err and "distinct source texts" not in err
    assert [r["classify_reason"] for r in read_tsv(again)] == ["decision"] * 4 + ["specimens+no_host"]
    assert read_tsv(tmp_path / "left.tsv") == []  # nothing left to review: a header and no rows
    assert (tmp_path / "left.tsv").read_text(encoding="utf-8").split("\t")[0] == "rows"


def test_classify_warns_about_decisions_that_matched_no_row(tmp_path: Path, capsys):
    table = tmp_path / "metadata.tsv"
    write_tsv(table, _rows("blood", host=""), ["id", "isolation_source", "host"])
    decided = _decisions(tmp_path, "isolation_source\tl1\nbed rails\tEnvironmental\nblood\tHuman\n")
    assert cli.main(["classify", str(table), "-o", str(tmp_path / "o.tsv"), "--decisions", str(decided)]) == 0
    assert "1 decision in" in (err := capsys.readouterr().err) and "matched no row (e.g. bed rails)" in err
    # A broken decisions table is one line on stderr and status 1, like any other input problem.
    assert cli.main(["classify", str(table), "-o", str(tmp_path / "o.tsv"), "--decisions",
                     str(_decisions(tmp_path, "isolation_source\tl1\nblood\tPerson\n"))]) == 1
    assert "'Person' is not a class" in capsys.readouterr().err
