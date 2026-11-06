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


def test_override_on_gcf_applies_to_gca_twin(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\nGCF_000005845.2\tLaboratory\n", encoding="utf-8")
    record = {"query_id": "GCA_000005845.2", "assembly": "GCA_000005845.2", "isolation_source": "urine"}
    result = classify_record(record, lexicon, load_overrides(path))
    assert (result.l1, result.override_applied) == ("Laboratory", "yes")


def test_override_on_other_assembly_version_does_not_apply(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\nGCA_000005845.2\tLaboratory\n", encoding="utf-8")
    record = {"query_id": "GCA_000005845.1", "isolation_source": "urine"}
    assert classify_record(record, lexicon, load_overrides(path)).override_applied == "no"


def test_override_l1_is_case_insensitive(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\nSAMN1\tanimal\n", encoding="utf-8")
    result = classify_record({"query_id": "SAMN1", "isolation_source": "blood"}, lexicon, load_overrides(path))
    assert (result.l1, result.classify_status) == ("Animal", "ok")


def test_invalid_override_paper_rollup_goes_to_review(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\tpaper_rollup\nSAMN1\tHuman\tHospital\n", encoding="utf-8")
    result = classify_record({"query_id": "SAMN1", "isolation_source": "blood"}, lexicon, load_overrides(path))
    assert result.classify_reason == "override_invalid_paper_rollup"
    assert result.paper_rollup == "Unclassified"


def test_override_blank_healthcare_equals_no(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("id\tl1\thealthcare_environment\nSAMN1\tHuman\t\nSAMN1\tHuman\tno\n", encoding="utf-8")
    result = classify_record({"query_id": "SAMN1"}, lexicon, load_overrides(path))
    assert result.classify_reason == "override"


def test_overrides_saved_by_excel_with_a_byte_order_mark_are_applied(tmp_path, lexicon):
    from onehealth.classify import classify_rows
    from onehealth.overrides import load_overrides

    path = tmp_path / "overrides.tsv"
    path.write_text("﻿id\tl1\nSAMN1\tAnimal\n", encoding="utf-8")
    row = classify_rows([{"id": "SAMN1", "isolation_source": "blood", "host": "Homo sapiens"}], lexicon, load_overrides(path))[0]
    assert (row["l1"], row["override_applied"]) == ("Animal", "yes")


def test_an_overrides_file_without_id_or_l1_columns_is_an_error(tmp_path):
    import pytest

    from onehealth.overrides import load_overrides

    no_id = tmp_path / "a.tsv"
    no_id.write_text("sample\tl1\nSAMN1\tAnimal\n", encoding="utf-8")
    no_l1 = tmp_path / "b.tsv"
    no_l1.write_text("id\tclass\nSAMN1\tAnimal\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no ID column"):
        load_overrides(no_id)
    with pytest.raises(ValueError, match="no l1 column"):
        load_overrides(no_l1)


def test_overrides_apply_to_tables_keyed_by_strain_or_sample(tmp_path: Path, lexicon):
    path = tmp_path / "overrides.tsv"
    path.write_text("strain\tl1\nKP1\tAnimal\n", encoding="utf-8")
    overrides = load_overrides(path)
    assert classify_record({"strain": "KP1", "isolation_source": "blood"}, lexicon, overrides).l1 == "Animal"
    assert classify_record({"sample": "KP1", "isolation_source": "blood"}, lexicon, overrides).l1 == "Animal"
