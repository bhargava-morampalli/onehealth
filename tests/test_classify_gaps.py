"""Gaps 2–4: food-producing / feed, lexicon holes, taxon-only source."""


from onehealth.classify import classify_record
from onehealth.lexicon import load_lexicon
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
    assert result.classify_reason == "animal-feed"
    assert result.l1 != "Food"


def test_livestock_feed_is_review_not_food(lexicon):
    result = _cls(lexicon, "livestock feed")
    assert result.classify_status == "review"
    assert result.classify_reason == "animal-feed"


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
    assert result.classify_reason == "fixture-no-setting"  # v4: review, with the reason named
    assert result.healthcare_environment == ""
    assert result.l1 == "Unclassified"


def test_hospital_faucet_is_healthcare_environment(lexicon):
    result = _cls(lexicon, "hospital faucet")
    assert result.l1 == "Environmental"
    assert result.healthcare_environment == "yes"
    assert result.classify_status == "ok"


def test_bare_faucet_is_not_healthcare_environment(lexicon):
    result = _cls(lexicon, "faucet")
    assert result.classify_status == "review"
    assert result.classify_reason == "fixture-no-setting"  # v4: review, with the reason named
    assert result.healthcare_environment == ""
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


def test_chicken_on_its_own_is_meat_and_chicken_from_a_chicken_is_the_animal(lexicon):
    """Until v7 a bare "chicken" was Animal. With no host it is Food, as in NCBI's records; see test_ncbi_conventions.py."""
    result = _cls(lexicon, "chicken")
    assert (result.l1, result.classify_status, result.classify_reason) == ("Food", "ok", "meat-species-no-host")
    assert _cls(lexicon, "chicken", host="Gallus gallus").l1 == "Animal"


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


def _custom_lexicon(tmp_path, fields):
    path = tmp_path / "lexicon.yaml"
    path.write_text(
        'version: "test"\nrules:\n'
        f"  - id: hospital-fixtures-extra\n    fields: {fields}\n"
        '    phrases: ["bedpan washer", "door handle", "bed rail"]\n'
        "    l1: Environmental\n    healthcare_environment: true\n    priority: 90\n"
        "  - id: animal-hosts\n    fields: [host, isolation_source]\n"
        '    phrases: ["bos taurus"]\n    l1: Animal\n    priority: 40\n',
        encoding="utf-8",
    )
    return load_lexicon(path)


def test_lexicon_rule_for_source_text_beats_the_taxon_guess(tmp_path):
    # "bedpan washer" has the shape of a binomial; a rule written for source text must still match it.
    lexicon = _custom_lexicon(tmp_path, "[isolation_source]")
    for source in ("bedpan washer", "Door handle", "bed rail"):
        result = _cls(lexicon, source)
        assert (result.l1, result.healthcare_environment, result.classify_reason) == (
            "Environmental", "yes", "hospital-fixtures-extra"), source


def test_binomial_matched_by_a_host_rule_is_still_ignored_in_isolation_source(tmp_path):
    lexicon = _custom_lexicon(tmp_path, "[isolation_source]")
    result = _cls(lexicon, "Bos taurus")
    assert (result.l1, result.classify_reason) == ("Unclassified", "taxon_only_source")
    # A rule that also reads host lists species names, so it does not overrule the guess.
    host_rule = _custom_lexicon(tmp_path, "[isolation_source, host]")
    assert _cls(host_rule, "door handle").classify_reason == "taxon_only_source"


def test_real_sources_are_not_reported_as_taxon_names(lexicon):
    # Every string the guess flagged in the IMP and NCBI Pathogen Detection data; none is a species.
    for source in ("bedpan washer", "respiratory tract infection", "urinary tract infection", "Cleaning car",
                   "WWTP influent"):
        assert is_taxon_only_source(source) is False, source
        assert _cls(lexicon, source).classify_reason != "taxon_only_source", source
    assert is_taxon_only_source("Klebsiella pneumoniae") is True
    assert is_taxon_only_source("Pseudomonas aeruginosa") is True


# --- lexicon v3 additions ------------------------------------------------------------------------


import pytest  # noqa: E402


@pytest.mark.parametrize(
    "source, host, expected",
    [
        ("bedpan washer", "", ("Environmental", "yes")),
        ("Bedpan", "", ("Environmental", "yes")),
        ("WC sink drains", "", ("Environmental", "yes")),
        ("Hospital waste water", "", ("Environmental", "yes")),
        ("waste water", "", ("Environmental", "no")),
        ("WWTP influent", "", ("Environmental", "no")),
        ("water", "", ("Environmental", "no")),
        ("environmental swab", "", ("Environmental", "no")),
        ("urinary tract infection", "", ("Human", "no")),
        ("respiratory tract infection", "Bos taurus", ("Animal", "no")),
        ("catheter", "", ("Human", "no")),
        ("indwelling urinary catheter", "Sus scrofa", ("Animal", "no")),
        ("human", "", ("Human", "no")),
        ("patient", "Homo sapiens", ("Human", "no")),
        ("passaging", "", ("Laboratory", "no")),
    ],
)
def test_lexicon_v3_additions(lexicon, source, host, expected):
    result = _cls(lexicon, source, host=host)
    assert (result.l1, result.healthcare_environment) == expected, result.classify_reason


@pytest.mark.parametrize(
    "source, reason",
    [
        ("Shower 3", "fixture-no-setting"),  # bare fixtures stay in review; v4 names the reason
        ("Sink swab", "fixture-no-setting"),
        ("infection control screening", "unmatched"),  # only "tract infection" is a specimen word
    ],
)
def test_lexicon_v3_leaves_these_for_review(lexicon, source, reason):
    result = _cls(lexicon, source)
    assert (result.classify_status, result.classify_reason) == ("review", reason)


def test_human_as_a_source_with_an_animal_host_is_a_conflict(lexicon):
    assert _cls(lexicon, "human", host="Gallus gallus").classify_reason == "host_source_conflict"
    assert _cls(lexicon, "patient room surface").l1 != "Human"


# --- lexicon v4: a fixture without a setting ---------------------------------------------------------


@pytest.mark.parametrize(
    "source, host, expected",
    [
        # A fixture alone says neither whose sample it is nor where: review, whatever a human host says.
        ("Rm 14 Drain", "Homo sapiens", ("review", "fixture-no-setting", "Unclassified")),
        ("Sink swab", "Homo sapiens", ("review", "fixture-no-setting", "Unclassified")),
        # v7: a floor is environmental, so with this host it is a conflict; still review.
        ("Floor Trap Bathroom Room 1 Swab", "Homo sapiens", ("review", "host_source_conflict", "Unclassified")),
        ("kitchen sink", "Homo sapiens", ("review", "fixture-no-setting", "Unclassified")),
        ("drain drainage", "Homo sapiens", ("review", "fixture-no-setting", "Unclassified")),
        ("Shower 3", "", ("review", "fixture-no-setting", "Unclassified")),
        # A setting word, a specimen or an animal host still decides, as in v3.
        ("hospital drain", "Homo sapiens", ("ok", "hospital-drain", "Environmental")),
        ("Sink drain, ward 4", "Homo sapiens", ("ok", "hospital-drain", "Environmental")),
        ("ICU sink", "Homo sapiens", ("ok", "healthcare-fixtures", "Environmental")),
        ("drain tip", "Homo sapiens", ("ok", "specimens+human-host", "Human")),
        ("drain fluid", "Homo sapiens", ("ok", "wound-drainage+human-host", "Human")),
        ("wound drain", "Homo sapiens", ("ok", "specimens+human-host", "Human")),
        ("drainage", "Homo sapiens", ("ok", "human-host", "Human")),  # whole words: drainage is not drain
        ("drain", "Sus scrofa", ("ok", "animal-hosts", "Animal")),
    ],
)
def test_lexicon_v4_fixture_without_a_setting(lexicon, source, host, expected):
    result = _cls(lexicon, source, host=host)
    assert (result.classify_status, result.classify_reason, result.l1) == expected
