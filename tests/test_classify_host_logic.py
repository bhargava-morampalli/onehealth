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


def test_source_type_disagreeing_with_host_is_settled_by_source_type(lexicon):
    """Until v7 this went to review as host_source_conflict. NCBI's own category now decides."""
    result = _cls(lexicon, host="Homo sapiens", source_type="food")
    assert (result.l1, result.classify_status, result.classify_reason) == ("Food", "ok", "source_type")
    assert result.lexicon_suggestion == "Food|Human"  # what the free text left open


@pytest.mark.parametrize(
    "source, host_disease, l1",
    [
        ("mastitis milk", "", "Animal"),
        ("milk", "bovine mastitis", "Animal"),
        ("breast milk", "", "Human"),
        ("retail chicken", "", "Food"),
        ("chicken air sac", "", "Animal"),
        ("marine fish", "", "Animal"),
        ("dairy farm", "", "Environmental"),  # v7: a farm is environmental, as NCBI records it; never Food
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


def test_taxon_or_unknown_in_source_does_not_shadow_sample_type(lexicon):
    assert _cls(lexicon, "Bos taurus", source_type="animal").l1 == "Animal"
    result = _cls(lexicon, "unknown tissue", sample_type="blood")
    assert (result.l1, result.classify_reason) == ("Human", "specimens+no_host")


def test_an_unknown_animal_is_still_an_animal(lexicon):
    # "unknown animal" leaves the species open, not the class, in every field.
    assert _cls(lexicon, host="unknown animal").l1 == "Animal"
    assert _cls(lexicon, "unknown animal").l1 == "Animal"
    assert _cls(lexicon, "Unknown animal faeces", host="Homo sapiens").classify_reason == "host_source_conflict"
    assert (_cls(lexicon, "unknown food").l1, _cls(lexicon, "unk human").l1) == ("Food", "Human")
    # "unknown" with only general words, or an INSDC "missing: ..." term, still says nothing.
    for value in ("unknown", "unknown source", "Unknown environment", "unknown environmental sample", "missing: lab stock"):
        result = _cls(lexicon, value)
        assert (result.classify_status, result.classify_reason) == ("review", "unknown_source"), value
    assert _cls(lexicon, "unknown environment", host="Homo sapiens").l1 == "Human"  # the host decides
    faeces = _cls(lexicon, "faeces", host="unknown animal")
    assert (faeces.l1, faeces.classify_reason) == ("Animal", "specimens+animal-hosts")
    assert _cls(lexicon, sample_type="unknown food").l1 == "Food"
    no_host = _cls(lexicon, "faeces", host="missing")
    assert (no_host.l1, no_host.classify_reason) == ("Human", "specimens+no_host")


def test_specimen_with_conflicting_non_host_source_type_goes_to_review(lexicon):
    # The lexicon alone cannot settle these; until v7 they went to review as
    # host_source_conflict. NCBI's own category now decides, and the open question stays
    # in lexicon_suggestion.
    food = _cls(lexicon, "urine", host="Homo sapiens", source_type="food")
    assert (food.l1, food.classify_status, food.classify_reason, food.lexicon_suggestion) == (
        "Food",
        "ok",
        "source_type",
        "Food|Human",
    )
    env = _cls(lexicon, "blood", source_type="environmental")
    assert (env.l1, env.classify_status, env.classify_reason, env.lexicon_suggestion) == (
        "Environmental",
        "ok",
        "source_type",
        "Environmental|Human",
    )


def test_words_in_the_specimen_text_do_not_conflict_with_it(lexicon):
    # Only another field (source_type, host) naming food, environment or laboratory is a conflict.
    cases = [("Stool sample from food handler", "Homo sapiens", "Human"), ("urine, laboratory", "Homo sapiens", "Human"),
             ("faeces of food animal", "Sus scrofa", "Animal")]
    for source, host, expected in cases:
        result = _cls(lexicon, source, host=host)
        assert (result.classify_status, result.l1) == ("ok", expected), source
    assert _cls(lexicon, sample_type="blood", isolation_source="hospital food").classify_reason == "host_source_conflict"

