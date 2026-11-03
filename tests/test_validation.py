"""onehealth validate: exact metrics on a hand-built table, then baselines on real curated data."""

import csv
from pathlib import Path

import pytest

from onehealth.classify import classify_rows
from onehealth.cli import main
from onehealth.io_tsv import read_tsv
from onehealth.validate import read_truth, validate

DATA = Path(__file__).parent / "data"
IMP = DATA / "imp_carbapenemase_sources" / "samples.tsv"
KLEBSIELLA = DATA / "ncbi_pd_klebsiella" / "metadata.tsv"


# --- exact metrics (hand-built: the expected numbers must be known) -----------------------


def _classified(*predictions):
    return [{"query_id": f"ISO{i}", "paper_rollup": p} for i, p in enumerate(predictions)]


def test_coverage_accuracy_and_agreement_are_exact():
    classified = _classified("Human", "Human", "Animal", "Unclassified", "Environmental", "Unclassified")
    truth = {"ISO0": "Human", "ISO1": "Animal", "ISO2": "Animal", "ISO3": "Human", "ISO4": "Environmental",
             "ISO5": "Unclassified", "MISSING": "Human"}
    result = validate(classified, truth)
    assert result.isolates == 6 and result.unmatched == ["MISSING"]
    assert result.coverage == pytest.approx(4 / 5)  # ISO3 went to review
    assert result.accuracy == pytest.approx(3 / 4)  # ISO1 called Human
    assert result.agreement == pytest.approx(4 / 6)  # ISO0, ISO2, ISO4, ISO5
    human = next(r for r in result.per_class() if r["class"] == "Human")
    assert (human["truth"], human["predicted"], human["correct"]) == (2, 2, 1)
    assert [d["id"] for d in result.disagreements] == ["ISO1", "ISO3"]


def test_truth_labels_must_be_classes(tmp_path: Path):
    path = tmp_path / "truth.tsv"
    path.write_text("id\ttruth\nISO0\thospital\n", encoding="utf-8")
    with pytest.raises(ValueError, match="map them before validating"):
        read_truth(path, "truth", "paper_rollup")
    path.write_text("id\ttruth\nISO0\tenvironment\nISO1\tHUMAN\n", encoding="utf-8")
    assert read_truth(path, "truth", "paper_rollup") == {"ISO0": "Environmental", "ISO1": "Human"}


# --- real data -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def imp_validation(lexicon):
    classified = classify_rows(read_tsv(IMP), lexicon)
    return validate(classified, read_truth(IMP, "truth", "paper_rollup"))


def test_imp_benchmark_baseline(imp_validation):
    """Floors just under the lexicon v2 baseline: lexicon changes must not make these worse."""
    result = imp_validation
    assert result.isolates == 4351 and not result.unmatched
    assert result.coverage >= 0.84
    assert result.accuracy >= 0.99
    precision = {r["class"]: r["precision"] for r in result.per_class()}
    assert precision["Human"] >= 0.99 and precision["Environmental"] >= 0.99


def test_imp_benchmark_unclassifiable_drains_are_documented(imp_validation):
    # The paper could not classify 9 "drain drainage" samples; with host Homo sapiens,
    # onehealth calls them Human. Recorded so a lexicon change that alters this is noticed.
    assert imp_validation.count(truth="Unclassified") == 9
    assert imp_validation.count("Unclassified", "Human") == 9


def test_ncbi_clinical_isolates_are_never_environmental_or_food(lexicon):
    """NCBI epi_type "clinical" (Pathogen.cl package) means host-associated."""
    classified = classify_rows(read_tsv(KLEBSIELLA), lexicon)
    clinical = [r for r in classified if r["epi_type"] == "clinical"]
    assert len(clinical) == 234
    assert sum(r["l1"] in ("Human", "Animal") for r in clinical) >= 233
    assert not [r for r in clinical if r["l1"] in ("Environmental", "Food")]


def test_validate_command_writes_all_tables(tmp_path: Path):
    classified = tmp_path / "classified.tsv"
    assert main(["classify", str(IMP), "-o", str(classified)]) == 0
    out = tmp_path / "validation.tsv"
    assert main(["validate", str(classified), "--truth", str(IMP), "-o", str(out)]) == 0
    summary = {r["metric"]: r["value"] for r in csv.DictReader(open(out), delimiter="\t")}
    assert summary["isolates"] == "4351"
    for table in ("by_class", "confusion", "disagreements"):
        assert (tmp_path / f"validation.{table}.tsv").exists()
