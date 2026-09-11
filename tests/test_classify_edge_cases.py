"""Paper edge cases: drainage vs drain, clinical environment, unknown."""

from onehealth.classify import classify_record
from onehealth.constants import derive_paper_rollup
from onehealth.overrides import load_overrides
from onehealth.text import contains_phrase, tokenize


def _cls(lexicon, isolation_source="", host="", **extra):
    row = {"isolation_source": isolation_source, "host": host, **extra}
    return classify_record(row, lexicon)


def test_wound_drainage_is_human_not_hospital_drain(lexicon):
    result = _cls(lexicon, "wound drainage")
    assert result.l1 == "Human"
    assert result.healthcare_environment == "no"
    assert result.paper_rollup == "Human"
    assert result.classify_status == "ok"


def test_hospital_drain_is_healthcare_environment(lexicon):
    result = _cls(lexicon, "hospital drain")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"
    assert result.paper_rollup == "Environmental"


def test_sink_drain_is_healthcare_environment(lexicon):
    result = _cls(lexicon, "ICU sink drain")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"


def test_drain_token_does_not_match_drainage():
    assert contains_phrase(tokenize("wound drainage"), tokenize("drain")) is False
    assert contains_phrase(tokenize("wound drainage"), tokenize("drainage")) is True
    assert contains_phrase(tokenize("hospital drain"), tokenize("drain")) is True
    assert contains_phrase(tokenize("wound drainage"), tuple(tokenize("wound drainage"))) is True


def test_clinical_is_healthcare_env_not_human(lexicon):
    result = _cls(lexicon, "clinical")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"
    assert result.paper_rollup == "Environmental"
    assert result.l1 != "Human"


def test_clinical_environment_is_healthcare_env_not_human(lexicon):
    result = _cls(lexicon, "clinical environment", host="Homo sapiens")
    assert result.healthcare_environment == "yes"
    assert result.l1 == "Environmental"
    assert result.paper_rollup == "Environmental"
    assert result.l1 != "Human"


def test_unknown_stays_unclassified_never_environmental(lexicon):
    result = _cls(lexicon, "unknown")
    assert result.l1 == "Unclassified"
    assert result.paper_rollup == "Unclassified"
    assert result.classify_status == "review"
    assert result.l1 != "Environmental"


def test_unknown_environment_never_maps_to_environmental(lexicon):
    result = _cls(lexicon, "unknown environment")
    assert result.l1 == "Unclassified"
    assert result.paper_rollup != "Environmental"
    assert result.classify_status == "review"


def test_empty_source_goes_to_review(lexicon):
    result = _cls(lexicon, "", host="")
    assert result.l1 == "Unclassified"
    assert result.classify_status == "review"
    assert result.classify_reason == "missing_source"


def test_human_blood_is_human(lexicon):
    result = _cls(lexicon, "blood", host="Homo sapiens")
    assert result.l1 == "Human"
    assert result.healthcare_environment == "no"
    assert result.paper_rollup == "Human"


def test_human_host_hospital_sink_is_healthcare_env(lexicon):
    result = _cls(lexicon, "hospital sink", host="Homo sapiens")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"
    assert result.paper_rollup == "Environmental"


def test_food_paper_rollup_is_environmental(lexicon):
    result = _cls(lexicon, "retail meat")
    assert result.l1 == "Food"
    assert result.paper_rollup == "Environmental"


def test_laboratory_paper_rollup_is_environmental(lexicon):
    result = _cls(lexicon, "type strain")
    assert result.l1 == "Laboratory"
    assert result.paper_rollup == "Environmental"


def test_animal_host(lexicon):
    result = _cls(lexicon, "cloacal swab", host="Gallus gallus")
    assert result.l1 == "Animal"
    assert result.paper_rollup == "Animal"


def test_derive_paper_rollup_never_promotes_unknown():
    assert derive_paper_rollup("Unclassified", "no") == "Unclassified"
    assert derive_paper_rollup("Unclassified", "yes") == "Unclassified"


def test_derive_paper_rollup_maps_food_lab_and_healthcare_to_environmental():
    assert derive_paper_rollup("Food", "no") == "Environmental"
    assert derive_paper_rollup("Laboratory", "no") == "Environmental"
    assert derive_paper_rollup("Human", "yes") == "Environmental"
