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
    """Floors just under the lexicon v4 baseline (84.5% coverage, 99.75% accuracy): changes must not make these worse."""
    result = imp_validation
    assert result.isolates == 4351 and not result.unmatched
    assert result.coverage >= 0.845
    assert result.accuracy >= 0.997
    precision = {r["class"]: r["precision"] for r in result.per_class()}
    assert precision["Human"] >= 0.99 and precision["Environmental"] >= 0.99


def test_imp_benchmark_drains_without_a_setting_go_to_review(imp_validation):
    # The paper could not classify 9 "drain drainage" samples, and v3 called them Human from
    # their host. v4 sends a fixture without a setting to review, so they now agree with the
    # curation, and only 3 environmental samples (stool, host Homo sapiens) are called Human.
    assert imp_validation.count(truth="Unclassified") == 9
    assert imp_validation.count("Unclassified", "Unclassified") == 9
    assert imp_validation.count("Environmental", "Human") == 3


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


def test_conflicting_truth_labels_are_an_error_not_a_coin_toss(tmp_path: Path, capsys):
    classified = tmp_path / "c.tsv"
    classified.write_text("query_id\tbiosample\tassembly\tpaper_rollup\nSAMN1\tSAMN1\tGCA_000000001.1\tHuman\n", encoding="utf-8")
    same_id = tmp_path / "t1.tsv"
    same_id.write_text("id\ttruth\nSAMN1\tHuman\nSAMN1\tAnimal\n", encoding="utf-8")
    two_ids = tmp_path / "t2.tsv"
    two_ids.write_text("id\ttruth\nSAMN1\tHuman\nGCA_000000001.1\tAnimal\n", encoding="utf-8")
    agree = tmp_path / "t3.tsv"
    agree.write_text("id\ttruth\nSAMN1\tHuman\nSAMN1\thuman\nGCA_000000001.1\tHuman\n", encoding="utf-8")
    for truth, message in ((same_id, "two different truth labels"), (two_ids, "listed under two IDs")):
        assert main(["validate", str(classified), "--truth", str(truth), "-o", str(tmp_path / "v.tsv")]) == 1
        assert message in capsys.readouterr().err
    assert main(["validate", str(classified), "--truth", str(agree), "-o", str(tmp_path / "v.tsv")]) == 0


def test_validate_with_no_overlap_says_n_a(tmp_path: Path, capsys):
    classified = tmp_path / "c.tsv"
    classified.write_text("query_id\tpaper_rollup\nSAMN1\tHuman\n", encoding="utf-8")
    truth = tmp_path / "t.tsv"
    truth.write_text("id\ttruth\nSAMN9\tHuman\n", encoding="utf-8")
    assert main(["validate", str(classified), "--truth", str(truth), "-o", str(tmp_path / "v.tsv")]) == 0
    assert "coverage n/a, accuracy n/a, agreement n/a" in capsys.readouterr().err


def test_lexicon_v2_still_ships_and_reproduces_its_benchmark():
    from onehealth.lexicon import load_lexicon

    v2 = load_lexicon(Path("v2"))
    result = validate(classify_rows(read_tsv(IMP), v2), read_truth(IMP, "truth", "paper_rollup"))
    assert v2.version == "2"
    assert (round(result.coverage, 4), round(result.accuracy, 4)) == (0.8448, 0.9918)  # 84.5% and 99.2%


def test_lexicon_v3_still_ships_and_reproduces_its_benchmark():
    from onehealth.lexicon import load_lexicon

    v3 = load_lexicon(Path("v3"))
    result = validate(classify_rows(read_tsv(IMP), v3), read_truth(IMP, "truth", "paper_rollup"))
    assert v3.version == "3"
    assert (round(result.coverage, 4), round(result.accuracy, 4)) == (0.8478, 0.9951)  # 84.8% and 99.5%
    assert result.count("Environmental", "Human") == 12


def test_lexicon_v3_classifies_more_of_the_ncbi_data_without_new_errors():
    """v3 classifies 10 more NCBI Pathogen Detection rows (waste water spellings) and gets none of them wrong."""
    from onehealth.lexicon import load_lexicon

    rows = read_tsv(KLEBSIELLA)
    expected = {"clinical": {"Human", "Animal"}, "environmental/other": {"Environmental"}}
    for version, classified_floor in (("v2", 253), ("v3", 263), ("v4", 263)):
        classified = [r for r in classify_rows(rows, load_lexicon(Path(version))) if r["classify_status"] == "ok"]
        labelled = [r for r in classified if r["epi_type"] in expected]
        wrong = [r for r in labelled if r["paper_rollup"] not in expected[r["epi_type"]]]
        assert len(classified) >= classified_floor and len(wrong) <= 2, version


def test_validate_matches_strain_tables_and_deduplicates_shared_biosamples(tmp_path: Path):
    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text("strain\ttruth\nKP1\tHuman\n", encoding="utf-8")
    truth = read_truth(truth_path, "truth", "paper_rollup")
    res = validate([{"strain": "KP1", "paper_rollup": "Human"}], truth)
    assert res.isolates == 1 and res.agreement == 1.0

    rows = [
        {"biosample": "SAMN1", "assembly": "GCA_000000001.1", "paper_rollup": "Human"},
        {"biosample": "SAMN1", "assembly": "GCA_000000002.1", "paper_rollup": "Human"},
    ]
    assert validate(rows, {"GCA_000000001.1": "Human", "GCA_000000002.1": "Human"}).isolates == 1
    with pytest.raises(ValueError, match="listed under two IDs"):
        validate(rows, {"GCA_000000001.1": "Human", "GCA_000000002.1": "Animal"})
