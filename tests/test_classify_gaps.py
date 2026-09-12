"""Gaps 2–4: food-producing / feed, lexicon holes, taxon-only source."""

from onehealth.classify import classify_record
from onehealth.text import is_taxon_only_source


def _cls(lexicon, isolation_source="", host="", **extra):
    row = {"isolation_source": isolation_source, "host": host, **extra}
    return classify_record(row, lexicon)


def test_animal_food_producing_is_animal_not_food(lexicon):
    result = _cls(lexicon, "Animal - food producing")
    assert result.l1 == "Animal"
    assert result.classify_status == "ok"
    assert result.l1 != "Food"
    assert result.paper_rollup == "Animal"


def test_food_producing_animal_is_animal(lexicon):
    result = _cls(lexicon, "food-producing animal")
    assert result.l1 == "Animal"
    assert result.classify_status == "ok"


def test_retail_meat_still_food(lexicon):
    result = _cls(lexicon, "retail meat")
    assert result.l1 == "Food"
    assert result.classify_status == "ok"


def test_animal_feed_is_review_not_food(lexicon):
    result = _cls(lexicon, "animal feed")
    assert result.l1 == "Unclassified"
    assert result.classify_status == "review"
    assert result.classify_reason == "animal_feed"
    assert result.l1 != "Food"


def test_livestock_feed_is_review_not_food(lexicon):
    result = _cls(lexicon, "livestock feed")
    assert result.classify_status == "review"
    assert result.classify_reason == "animal_feed"


def test_livestock_is_animal(lexicon):
    result = _cls(lexicon, "livestock")
    assert result.l1 == "Animal"
    assert result.classify_status == "ok"


def test_companion_is_animal(lexicon):
    result = _cls(lexicon, "Animal - companion")
    assert result.l1 == "Animal"
    assert result.classify_status == "ok"


def test_fresh_water_is_environmental_not_healthcare(lexicon):
    result = _cls(lexicon, "fresh water")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "no"
    assert result.classify_status == "ok"


def test_freshwater_token_is_environmental(lexicon):
    result = _cls(lexicon, "freshwater")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "no"


def test_experimentally_derived_is_laboratory(lexicon):
    result = _cls(lexicon, "experimentally derived")
    assert result.l1 == "Laboratory"
    assert result.classify_status == "ok"
    assert result.paper_rollup == "Environmental"


def test_bronchial_lavage_is_human(lexicon):
    result = _cls(lexicon, "bronchial lavage")
    assert result.l1 == "Human"
    assert result.healthcare_environment == "no"
    assert result.classify_status == "ok"


def test_hospital_carriage_is_human_not_fixture(lexicon):
    result = _cls(lexicon, "hospital carriage")
    assert result.l1 == "Human"
    assert result.healthcare_environment == "no"
    assert result.paper_rollup == "Human"
    assert result.classify_status == "ok"


def test_hospital_sink_is_healthcare_environment(lexicon):
    result = _cls(lexicon, "hospital sink")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"
    assert result.classify_status == "ok"


def test_bare_sink_is_not_healthcare_environment(lexicon):
    result = _cls(lexicon, "sink")
    assert result.classify_status == "review"
    assert result.classify_reason == "unmatched"
    assert result.healthcare_environment == "no"
    assert result.l1 == "Unclassified"


def test_hospital_faucet_is_healthcare_environment(lexicon):
    result = _cls(lexicon, "hospital faucet")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"
    assert result.classify_status == "ok"


def test_bare_faucet_is_not_healthcare_environment(lexicon):
    result = _cls(lexicon, "faucet")
    assert result.classify_status == "review"
    assert result.classify_reason == "unmatched"
    assert result.healthcare_environment == "no"
    assert result.l1 != "Environmental"


def test_escherichia_coli_source_is_review_not_specimen(lexicon):
    result = _cls(lexicon, "Escherichia coli")
    assert result.classify_status == "review"
    assert result.classify_reason == "taxon_only_source"
    assert result.l1 == "Unclassified"
    assert result.lexicon_suggestion == ""


def test_e_coli_source_is_review(lexicon):
    result = _cls(lexicon, "E. coli")
    assert result.classify_reason == "taxon_only_source"
    assert result.classify_status == "review"


def test_bos_taurus_isolation_source_is_review_not_animal(lexicon):
    result = _cls(lexicon, "Bos taurus")
    assert result.classify_status == "review"
    assert result.classify_reason == "taxon_only_source"
    assert result.l1 != "Animal"


def test_bos_taurus_host_without_source_is_still_animal(lexicon):
    result = _cls(lexicon, "", host="Bos taurus")
    assert result.l1 == "Animal"
    assert result.classify_status == "ok"


def test_taxon_plus_specimen_still_classifies(lexicon):
    result = _cls(lexicon, "Escherichia coli from urine")
    assert result.l1 == "Human"
    assert result.classify_status == "ok"
    assert result.classify_reason != "taxon_only_source"


def test_chicken_isolation_source_is_still_animal(lexicon):
    result = _cls(lexicon, "chicken")
    assert result.l1 == "Animal"
    assert result.classify_status == "ok"


def test_is_taxon_only_source_accepts_binomials_and_rejects_specimens():
    assert is_taxon_only_source("Escherichia coli") is True
    assert is_taxon_only_source("Bos taurus") is True
    assert is_taxon_only_source("Homo sapiens") is True
    assert is_taxon_only_source("E. coli") is True
    assert is_taxon_only_source("canis lupus familiaris") is True
    assert is_taxon_only_source("urine") is False
    assert is_taxon_only_source("chicken") is False
    assert is_taxon_only_source("surface water") is False
    assert is_taxon_only_source("retail meat") is False
    assert is_taxon_only_source("wound drainage") is False
    assert is_taxon_only_source("Escherichia coli from urine") is False
    assert is_taxon_only_source("experimentally derived") is False
    assert is_taxon_only_source("") is False
