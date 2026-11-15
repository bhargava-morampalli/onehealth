"""Food, food premises and food handlers (lexicon v5): retail meat is Food, a swab of food premises is Environmental,
and a person who handles food is Human. v4 called the first Animal and the other two Food."""

from pathlib import Path

import pytest

from onehealth.classify import classify_record
from onehealth.constants import SEED_LEXICON_VERSION as SEED
from onehealth.lexicon import load_lexicon


def _classify(lexicon, source="", host="", **more):
    return classify_record({"isolation_source": source, "host": host, **more}, lexicon)


@pytest.mark.parametrize(
    "source",
    [
        "raw intact chicken",
        "comminuted chicken",
        "ground turkey",
        "chicken breast",
        "skinless boneless chicken breast",
        "chicken wings",
        "smoked salmon",
        "fish fillet",
        "pork chop",
        "Product-Raw-Intact-Beef",
        "RTE product",
        "peanut butter",  # two plain words: v4 took it for a species name
        "ice cream",
        "enoki mushroom",
        "avocado",
    ],
)
def test_meat_cuts_and_foods_are_food_not_animal(lexicon, source):
    result = _classify(lexicon, source)
    assert (result.l1, result.classify_status, result.classify_reason) == ("Food", "ok", "food-items")
    assert result.paper_rollup == "Environmental"


@pytest.mark.parametrize(
    "source, host, l1, reason",
    [
        # The species with a host, or in the plural, still means the animal. On its own
        # with no host, "chicken" is meat since v7 (test_ncbi_conventions.py).
        ("chicken", "Gallus gallus", "Animal", "animal-hosts"),
        ("chickens", "", "Animal", "animal-hosts"),
        ("cecum", "Gallus gallus", "Animal", "animal-hosts"),
        # An organ is not a cut: retail giblets or a necropsy sample.
        ("chicken liver", "", "Animal", "animal-hosts"),
        ("chicken liver", "Gallus gallus", "Animal", "animal-hosts"),
        ("chicken liver (giblets)", "", "Food", "food-items"),
        # A bare body part can be a patient's.
        ("breast", "Homo sapiens", "Human", "human-host"),
        ("breast abscess", "Homo sapiens", "Human", "specimens+human-host"),
        ("breast milk", "", "Human", "specimens+no_host"),
        # Livestock named with a food word.
        ("food animal", "", "Animal", "animal-food-producing"),
        ("beef cattle", "", "Animal", "animal-hosts"),
    ],
)
def test_what_v5_leaves_alone(lexicon, source, host, l1, reason):
    result = _classify(lexicon, source, host)
    assert (result.l1, result.classify_reason) == (l1, reason)


def test_a_cut_of_meat_with_an_animal_host_goes_to_review(lexicon):
    """As "chicken meat" with host chicken always did: the text says food, the host says animal."""
    result = _classify(lexicon, "ground chicken", "Gallus gallus")
    assert (result.classify_status, result.classify_reason) == ("review", "host_source_conflict")
    assert result.lexicon_suggestion == "Animal|Food"


@pytest.mark.parametrize(
    "source",
    [
        "environmental swab sponge non food contact surface",
        "Environmental: food-contact surface",
        "food processing environment",
        "production environment",
        "poultry processing plant",
        "environmental swab from a pet food plant",
        "environmental sponge",
    ],
)
def test_swabs_of_food_premises_are_environmental_not_food(lexicon, source):
    result = _classify(lexicon, source)
    assert (result.l1, result.healthcare_environment, result.classify_reason) == ("Environmental", "no", "environmental-sampling")


def test_hospital_and_specimen_rules_still_outrank_environmental_sampling(lexicon):
    assert _classify(lexicon, "hospital surface environmental swab").healthcare_environment == "yes"
    assert _classify(lexicon, "environmental swab", "Homo sapiens").classify_reason == "host_source_conflict"
    assert _classify(lexicon, "stool from processing plant worker").classify_reason == "specimens+no_host"


def test_a_food_handler_is_a_person(lexicon):
    assert (_classify(lexicon, "food handler").l1, _classify(lexicon, "food handler").classify_reason) == ("Human", "food-handler")
    # A specimen from one: still Human, and the reason now names the person in place of "no_host".
    assert (_classify(lexicon, "food handler stool").l1, _classify(lexicon, "food handler stool").classify_reason) == (
        "Human", "specimens+food-handler")
    assert _classify(lexicon, "food handler", "Bos taurus").classify_reason == "host_source_conflict"


@pytest.mark.parametrize("source", ["fecal", "faecal sample", "fecal swab", "stools"])
def test_other_spellings_of_faeces_follow_the_host(lexicon, source):
    assert _classify(lexicon, source, "Homo sapiens").classify_reason == "specimens+human-host"
    assert _classify(lexicon, source, "Sus scrofa").classify_reason == "specimens+animal-hosts"


def test_one_lexicon_ships_and_an_older_version_number_says_where_to_find_it():
    """`--lexicon v3` used to select an older lexicon shipped next to the current one."""
    assert load_lexicon(Path(f"v{SEED}")).version == load_lexicon(Path(SEED)).version == load_lexicon().version == SEED
    for old in ("v4", "v3", "2"):
        with pytest.raises(ValueError, match=f"onehealth ships one lexicon, version {SEED}, which is the default"):
            load_lexicon(Path(old))


# --- a specimen with no host, filed by its submitter as not clinical ---------------------------------


@pytest.mark.parametrize(
    "column, value",
    [
        ("package", "Pathogen.env.1.0"),
        ("attribute_package", "Pathogen: environmental/food/other"),
        ("epi_type", "environmental/other"),
    ],
)
def test_specimen_without_host_in_an_environmental_package_goes_to_review(lexicon, column, value):
    result = _classify(lexicon, "feces", **{column: value})
    assert (result.l1, result.classify_status, result.classify_reason) == ("Unclassified", "review", "no_host_env_package")
    assert result.lexicon_suggestion == "Human"  # what the default would have been


@pytest.mark.parametrize(
    "more",
    [
        {},  # no package column at all: as before
        {"package": "Pathogen.cl.1.0"},
        {"attribute_package": "Pathogen: clinical or host-associated"},
        {"epi_type": "clinical"},
        {"package": "Microbe.1.0"},
        {"package": "MIGS.ba.human-associated.6.0"},
    ],
)
def test_specimen_without_host_keeps_its_default_in_any_other_package(lexicon, more):
    result = _classify(lexicon, "feces", **more)
    assert (result.l1, result.classify_reason) == ("Human", "specimens+no_host")


def test_the_package_never_overrides_a_host_or_a_source(lexicon):
    env = {"package": "Pathogen.env.1.0"}
    assert _classify(lexicon, "feces", "Homo sapiens", **env).classify_reason == "specimens+human-host"
    assert _classify(lexicon, "pig feces", **env).classify_reason == "specimens+animal-hosts"
    assert _classify(lexicon, "river water", **env).l1 == "Environmental"
    assert _classify(lexicon, "", "Homo sapiens", **env).l1 == "Human"
