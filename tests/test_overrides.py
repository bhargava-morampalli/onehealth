from pathlib import Path

from onehealth.classify import classify_record
from onehealth.overrides import load_overrides


def test_override_always_wins(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\thealthcare_environment\nSAMN1\tFood\tno\n", encoding="utf-8")
    overrides = load_overrides(path)
    result = classify_record(
        {"query_id": "SAMN1", "biosample": "SAMN1", "isolation_source": "blood", "host": "Homo sapiens"},
        lexicon,
        overrides,
    )
    assert result.override_applied == "yes"
    assert result.l1 == "Food"
    assert result.classify_status == "ok"
    assert result.lexicon_suggestion == "Human"


def test_empty_override_goes_to_review_not_lexicon(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\nSAMN1\t\n", encoding="utf-8")
    overrides = load_overrides(path)
    result = classify_record(
        {"query_id": "SAMN1", "isolation_source": "blood", "host": "Homo sapiens"},
        lexicon,
        overrides,
    )
    assert result.classify_status == "review"
    assert result.classify_reason == "override_empty"
    assert result.l1 == "Unclassified"
    assert result.lexicon_suggestion == "Human"


def test_conflict_override_goes_to_review(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\nSAMN1\tHuman\nSAMN1\tAnimal\n", encoding="utf-8")
    overrides = load_overrides(path)
    result = classify_record(
        {"query_id": "SAMN1", "isolation_source": "blood"},
        lexicon,
        overrides,
    )
    assert result.classify_status == "review"
    assert result.classify_reason == "override_conflict"
    assert result.l1 == "Unclassified"


def test_invalid_override_l1_goes_to_review(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\nSAMN1\tHospital\n", encoding="utf-8")
    result = classify_record({"query_id": "SAMN1", "isolation_source": "blood"}, lexicon, load_overrides(path))
    assert result.classify_status == "review"
    assert result.classify_reason == "override_invalid_l1"
