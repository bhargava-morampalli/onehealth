"""Specimen words are host-neutral; host and source_type decide."""

import pytest

from onehealth.classify import classify_record


def _cls(lexicon, isolation_source="", host="", **extra):
    return classify_record({"isolation_source": isolation_source, "host": host, **extra}, lexicon)


@pytest.mark.parametrize(
    "source, host",
    [
        ("pig feces", ""),
        ("bovine faeces", ""),
        ("chicken tissue", ""),
        ("dog urine", ""),
        ("swine nasal swab", ""),
        ("feces", "Bos taurus"),
        ("urine", "Canis lupus familiaris"),
        ("blood", "Gallus gallus"),
        ("clinical", "Sus scrofa"),
    ],
)
def test_specimen_from_an_animal_is_animal(lexicon, source, host):
    result = _cls(lexicon, source, host=host)
    assert (result.l1, result.classify_status, result.paper_rollup) == ("Animal", "ok", "Animal")


def test_specimen_with_human_host_is_human(lexicon):
    result = _cls(lexicon, "urine", host="Homo sapiens")
    assert result.l1 == "Human"
    assert result.classify_reason == "specimens+human-host"


def test_specimen_without_host_uses_the_documented_human_default(lexicon):
    result = _cls(lexicon, "blood")
    assert result.l1 == "Human"
    assert result.classify_reason == "specimens+no_host"


def test_specimen_with_human_host_and_animal_word_is_conflict(lexicon):
    result = _cls(lexicon, "pig bite wound", host="Homo sapiens")
    assert result.classify_status == "review"
    assert result.classify_reason == "host_source_conflict"
    assert result.lexicon_suggestion == "Animal|Human"


def test_specimen_with_unrecognised_host_is_review_not_human(lexicon):
    # Klebsiella variicola is a maize endophyte; a plant host must not default to Human.
    result = _cls(lexicon, "plant tissue", host="Zea mays")
    assert result.classify_status == "review"
    assert result.classify_reason == "host_unrecognised"


def test_taxon_name_in_source_falls_back_to_host(lexicon):
    result = _cls(lexicon, "Klebsiella pneumoniae", host="Homo sapiens")
    assert (result.l1, result.classify_status, result.classify_reason) == ("Human", "ok", "human-host")


def test_unknown_source_falls_back_to_source_type(lexicon):
    result = _cls(lexicon, "unknown", source_type="Food")
    assert (result.l1, result.classify_status) == ("Food", "ok")


def test_insdc_missing_value_is_unknown(lexicon):
    result = _cls(lexicon, "missing: lab stock")
    assert result.classify_reason == "unknown_source"


@pytest.mark.parametrize(
    "source_type, l1",
    [("human", "Human"), ("animal", "Animal"), ("environmental", "Environmental"), ("food", "Food")],
)
def test_ncbi_source_type_vocabulary_classifies(lexicon, source_type, l1):
    result = _cls(lexicon, source_type=source_type)
    assert (result.l1, result.classify_status) == (l1, "ok")


def test_ncbi_sample_type_material_is_not_source_evidence(lexicon):
    result = _cls(lexicon, sample_type="whole organism")
    assert result.classify_reason == "missing_source"


def test_source_type_wins_over_material_sample_type(lexicon):
    result = _cls(lexicon, sample_type="pure culture", source_type="animal")
    assert result.l1 == "Animal"


def test_source_type_disagreeing_with_host_is_conflict(lexicon):
    result = _cls(lexicon, host="Homo sapiens", source_type="food")
    assert result.classify_reason == "host_source_conflict"


@pytest.mark.parametrize(
    "source, host_disease, l1",
    [
        ("mastitis milk", "", "Animal"),
        ("milk", "bovine mastitis", "Animal"),
        ("breast milk", "", "Human"),
        ("retail chicken", "", "Food"),
        ("chicken air sac", "", "Animal"),
        ("marine fish", "", "Animal"),
        ("dairy farm", "", "Animal"),
        ("dairy products", "", "Food"),
    ],
)
def test_v2_lexicon_corrections(lexicon, source, host_disease, l1):
    result = _cls(lexicon, source, host_disease=host_disease)
    assert (result.l1, result.classify_status) == (l1, "ok")


def test_ventilator_associated_pneumonia_is_not_a_fixture(lexicon):
    result = _cls(lexicon, "ventilator-associated pneumonia")
    assert result.healthcare_environment != "yes"
    assert _cls(lexicon, "ventilator").healthcare_environment == "yes"


def test_wound_drainage_and_sink_drain_together_is_conflict_not_environmental(lexicon):
    result = _cls(lexicon, "wound drainage, sink drain")
    assert result.classify_status == "review"
    assert result.classify_reason.startswith("rule_conflict:")


def test_review_rows_make_no_healthcare_claim(lexicon):
    result = _cls(lexicon, "")
    assert result.healthcare_environment == ""
    assert result.paper_rollup == "Unclassified"


def test_hospital_sink_with_human_host_stays_healthcare_environment(lexicon):
    result = _cls(lexicon, "hospital sink", host="Homo sapiens")
    assert (result.l1, result.healthcare_environment) == ("Environmental", "yes")
