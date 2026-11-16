"""Where onehealth follows NCBI (lexicon v7 and classify): documentation/guides/definitions.md.

Four conventions, each measured in NCBI Pathogen Detection records whose submitters filled in
NCBI's own source_type (tests/data/ncbi_pd_source_types): a record's source_type is its class;
what animals leave behind and where they are kept is environmental; a species that is also
the name of its meat is the meat when there is no host; and a specimen with no host is a person's.
"""

from collections import Counter
from pathlib import Path

import pytest

from onehealth import cli
from onehealth.classify import classify_record, classify_rows
from onehealth.io_tsv import read_tsv, write_tsv
from onehealth.lexicon import load_lexicon
from onehealth.ncbi import attributes_from_biosample
from onehealth.overrides import OverrideRow
from onehealth.resolve import resolve_ids
from onehealth.review import load_decisions
from onehealth.taxonomy import taxon_kind
from onehealth.text import contains_phrase, is_unknown_value, tokenize

RECORDS = Path(__file__).parent / "data" / "ncbi_pd_source_types" / "records.tsv"


def _classify(lexicon, source="", host="", **more):
    return classify_record({"isolation_source": source, "host": host, **more}, lexicon)


def _labels(phrase: str, exact: bool = False, no_host: bool = False) -> Counter:
    """NCBI source_type of the benchmark records whose isolation_source has this phrase."""
    words, labels = phrase.split(), Counter()
    for row in read_tsv(RECORDS):
        tokens = tokenize(row["isolation_source"])
        if (tokens == words if exact else contains_phrase(tokens, words)) and (not no_host or is_unknown_value(row["host"])):
            labels[row["ncbi_source_type"]] += int(row["records"])
    return labels


# --- 1. a record's source_type is its class ---------------------------------------------------------


@pytest.mark.parametrize(
    "source, host, source_type, l1, reason, suggestion",
    [
        # The free text is silent or agrees: the lexicon's own reason is kept.
        ("", "", "food", "Food", "food-products", "Food"),
        ("river water", "", "environmental", "Environmental", "environment-outdoor", "Environmental"),
        # The free text cannot settle it: these went to review until v7.
        ("animal-cattle-steer", "bovine", "food", "Food", "source_type", "Animal|Food"),
        ("bovine carcass", "bovine", "animal", "Animal", "source_type", "Animal|Food"),
        ("stool", "Gallus gallus", "human", "Human", "source_type", "Animal|Human"),
        ("animal feed", "", "food", "Food", "source_type", ""),
        ("feces", "Zea mays", "environmental", "Environmental", "source_type", "Environmental|Human"),
        # The free text says something else: NCBI's category is used, and the reason says so.
        ("raw milk", "", "environmental", "Environmental", "source_type_over_text", "Food"),
        ("laboratory", "", "human", "Human", "source_type_over_text", "Laboratory"),
        ("ground turkey", "", "animal", "Animal", "source_type_over_text", "Food"),
    ],
)
def test_source_type_is_the_class(lexicon, source, host, source_type, l1, reason, suggestion):
    result = _classify(lexicon, source, host, source_type=source_type)
    assert (result.l1, result.classify_status, result.classify_reason) == (l1, "ok", reason)
    assert result.lexicon_suggestion == suggestion
    assert result.paper_rollup == {"Food": "Environmental"}.get(l1, l1)


@pytest.mark.parametrize("value", ["other", "Other", "not provided", "missing", "clinical", ""])
def test_a_source_type_that_names_no_class_declares_nothing(lexicon, value):
    result = _classify(lexicon, "blood", "Homo sapiens", source_type=value)
    assert (result.l1, result.classify_reason) == ("Human", "specimens+human-host")


def test_source_type_keeps_the_healthcare_flag_and_loses_to_a_curator(lexicon, tmp_path: Path):
    hospital = _classify(lexicon, "hospital sink", "Homo sapiens", source_type="environmental")
    assert (hospital.l1, hospital.healthcare_environment, hospital.classify_reason) == ("Environmental", "yes", "healthcare-fixtures")
    # Only the column named source_type declares a class. "tissue" and NCBI's sample_type are free text.
    assert _classify(lexicon, "chicken breast", tissue="animal").classify_reason == "food-items"
    row = {"id": "ISO1", "isolation_source": "stool", "host": "Gallus gallus", "source_type": "human"}
    override = {"ISO1": [OverrideRow(key="ISO1", l1="Animal", healthcare_environment="", paper_rollup="", note="", complete=True)]}
    assert classify_record(row, lexicon, override).classify_reason == "override"
    decided = tmp_path / "d.tsv"
    decided.write_text("isolation_source\thost\tsample_type\tl1\nstool\tGallus gallus\thuman\tAnimal\n", encoding="utf-8")
    assert classify_record(row, lexicon, None, load_decisions(decided)).classify_reason == "decision"


def test_resolve_writes_source_type_in_a_column_of_its_own(tmp_path: Path):
    """sample_type holds source_type, or tissue, or NCBI's material type: classify could not tell which."""
    report = {"accession": "SAMN1", "attributes": [
        {"name": "isolation_source", "value": "animal-cattle-steer"}, {"name": "host", "value": "bovine"},
        {"name": "source_type", "value": "Food"}, {"name": "sample_type", "value": "whole organism"}]}
    fields = attributes_from_biosample(report)
    assert (fields["source_type"], fields["sample_type"]) == ("Food", "Food")
    assert attributes_from_biosample({"accession": "SAMN2", "attributes": [{"name": "tissue", "value": "human"}]})["source_type"] == ""
    local = resolve_ids(["SAMN3"], local_rows=[{"biosample": "SAMN3", "isolation_source": "bovine carcass", "host": "bovine",
                                                "source_type": "animal"}])
    assert local[0]["source_type"] == "animal"
    table = tmp_path / "resolved.tsv"
    write_tsv(table, local, list(local[0]))
    out = tmp_path / "classified.tsv"
    assert cli.main(["classify", str(table), "-o", str(out)]) == 0
    row = read_tsv(out)[0]
    assert (row["l1"], row["classify_reason"], row["lexicon_suggestion"]) == ("Animal", "source_type", "Animal|Food")


# --- 2. farm settings ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "source, host, l1, reason",
    [
        # What animals leave behind is environmental, whatever the host.
        ("manure", "", "Environmental", "animal-waste"),
        ("dairy manure", "", "Environmental", "animal-waste"),
        ("poultry litter", "", "Environmental", "animal-waste"),
        ("dairy cattle manure", "Bos taurus", "Environmental", "animal-waste"),
        # A place: with no host the sample is of the place, with a host it is of the host.
        ("farm", "", "Environmental", "farm-environment"),
        ("dairy farm calf barn", "", "Environmental", "farm-environment"),
        ("poultry shed", "", "Environmental", "farm-environment"),
        ("pig farm", "Sus scrofa domesticus", "Animal", "farm-environment+animal-hosts"),
        ("organic meat chicken farm", "Gallus gallus", "Animal", "farm-environment+animal-hosts"),
        ("farm", "Homo sapiens", "Human", "farm-environment+human-host"),
        ("farm", "Larus dominicanus", "Animal", "farm-environment+host-taxonomy"),
        # A farm is a farm even when its text names a food, and a farmed food is a food.
        ("water: shrimp farm", "", "Environmental", "farm-environment"),
        ("frozen farm raised large mouth bass", "", "Food", "farmed-food"),
        # Swabs of a poultry house.
        ("drag swab poultry", "", "Environmental", "environmental-sampling"),
        ("boot swab", "", "Environmental", "environmental-sampling"),
        # A specimen still follows its host, and a bird's body site is the bird.
        ("feces, calf barn", "Bos taurus", "Animal", "specimens+animal-hosts"),
        ("feces from self-breeding pig farm", "pig", "Animal", "specimens+animal-hosts"),
        ("farm worker stool", "", "Human", "specimens+no_host"),  # until v7 "farm" made this Animal
        ("cloaca swab farm", "", "Animal", "animal-sites"),
    ],
)
def test_farm_settings_follow_ncbi(lexicon, source, host, l1, reason):
    result = _classify(lexicon, source, host)
    assert (result.l1, result.classify_status, result.classify_reason) == (l1, "ok", reason)


def test_ncbi_records_farm_settings_as_environmental():
    """The figures quoted in documentation/guides/definitions.md, from all 139,554 benchmark records."""
    assert _labels("manure") == {"environmental": 646, "animal": 101}
    assert _labels("manure", no_host=True) == {"environmental": 567, "animal": 42}
    for phrase, n in (("litter", 77), ("bedding", 68), ("barn", 19), ("poultry shed", 69), ("drag swab", 152)):
        assert _labels(phrase) == {"environmental": n}, phrase
    farm_no_host, farm = _labels("farm", no_host=True), _labels("farm")
    assert farm_no_host == {"environmental": 141, "food": 30, "animal": 16}
    assert farm - farm_no_host == {"animal": 507, "environmental": 125, "human": 9}
    assert _labels("floor drain") == {"environmental": 88}


# --- 3. a species that is also the name of its meat -------------------------------------------------


def test_chicken_and_turkey_with_no_host_are_food(lexicon):
    for word in ("chicken", "Chicken", "turkey"):
        result = _classify(lexicon, word)
        assert (result.l1, result.classify_reason, result.paper_rollup) == ("Food", "meat-species-no-host", "Environmental")
    # With a host the rule does not apply: the host decides, as before.
    assert _classify(lexicon, "chicken", "Gallus gallus").classify_reason == "animal-hosts"
    assert _classify(lexicon, "chicken", "Homo sapiens").classify_reason == "host_source_conflict"
    # A missing-value word in host is no host.
    assert _classify(lexicon, "turkey", "not collected").classify_reason == "meat-species-no-host"
    # Only the bare word, and only these two: no other has the evidence.
    for source, reason in (("live chicken", "animal-hosts"), ("chicken cecum", "animal-hosts"), ("duck", "animal-hosts"),
                           ("fish", "animal-hosts"), ("pig", "animal-hosts")):
        assert (_classify(lexicon, source).l1, _classify(lexicon, source).classify_reason) == ("Animal", reason), source


def test_ncbi_records_a_bare_chicken_or_turkey_as_food():
    chicken, turkey = _labels("chicken", exact=True, no_host=True), _labels("turkey", exact=True, no_host=True)
    assert chicken == {"food": 810, "animal": 56, "human": 1}
    assert turkey == {"food": 255, "animal": 1}


# --- other v7 changes ------------------------------------------------------------------------------


# --- 4. a specimen with no host is a person's -----------------------------------------------------


def test_ncbi_records_a_specimen_with_no_host_as_human(lexicon):
    """The default is older than v7. It was measured against NCBI's categories in November 2026 and kept."""
    labels, named_animal, shellfish = Counter(), 0, 0
    for row in read_tsv(RECORDS):
        result = classify_record({k: row[k] for k in ("isolation_source", "host", "host_disease")}, lexicon)
        if result.classify_reason.endswith("+no_host"):
            assert result.l1 == "Human"
            n = int(row["records"])
            labels[row["ncbi_source_type"]] += n
            words = tokenize(row["isolation_source"])
            shellfish += n if {"oyster", "catfish"} & set(words) else 0
            named_animal += n if any(
                taxon_kind(" ".join(words[i:j])) == "animal" for i in range(len(words)) for j in range(i + 1, min(i + 4, len(words)) + 1)
            ) and row["ncbi_source_type"] == "animal" else 0
    assert sum(labels.values()) == 15235
    assert labels == {"human": 14735, "animal": 345, "environmental": 117, "food": 38}
    assert round(100 * labels["human"] / 15235, 1) == 96.7
    # The two kinds of exception a text could settle, listed in the documentation as not yet handled.
    assert named_animal == 78 and shellfish == 120


@pytest.mark.parametrize(
    "source, host, l1, healthcare, reason",
    [
        # A floor drain is a food plant's as often as a hospital's: environmental, no healthcare flag.
        ("floor drain in produce prep room", "", "Environmental", "no", "environmental-sampling"),
        ("floor drain", "", "Environmental", "no", "environmental-sampling"),
        ("hospital floor drain", "", "Environmental", "yes", "healthcare-fixtures"),
        ("sink drain", "", "Environmental", "yes", "hospital-drain"),
        # Fixtures and equipment of food premises, even when the text names the food.
        ("drain-raw meat", "", "Environmental", "no", "environmental-sampling"),
        ("deli case-interior deli case (by potato and egg salad)", "", "Environmental", "no", "environmental-sampling"),
        ("utensils-raw chicken/egg white cutting boards", "", "Environmental", "no", "environmental-sampling"),
        # Water in any phrase, and water bodies.
        ("pond water", "", "Environmental", "no", "environment-outdoor"),
        ("creek", "", "Environmental", "no", "environment-outdoor"),
        ("water bath", "", "Environmental", "no", "environment-outdoor"),
        ("bottled water", "", "Food", "no", "food-items"),
        ("water buffalo", "", "Animal", "no", "animal-food-producing"),
        # Hospital equipment.
        ("Bed rails", "", "Environmental", "yes", "healthcare-fixtures"),
        ("ECG monitors", "Homo sapiens", "Environmental", "yes", "healthcare-fixtures"),
        # A scope is how a patient sample was taken: the IMP curators call these human.
        ("Bronchoscope", "Homo sapiens", "Human", "no", "human-host"),
        # More body sites and fluids, and "human ..." as a whole value.
        ("rectal infection", "", "Human", "no", "specimens+no_host"),
        ("respiratory secretion", "Canis lupus familiaris", "Animal", "no", "specimens+animal-hosts"),
        ("Human Listeriosis", "", "Human", "no", "human-source"),
    ],
)
def test_other_v7_vocabulary(lexicon, source, host, l1, healthcare, reason):
    result = _classify(lexicon, source, host)
    assert (result.l1, result.healthcare_environment, result.classify_reason) == (l1, healthcare, reason)


def test_a_specimen_word_in_a_food_is_not_a_patient(lexicon):
    """"blood" is a specimen, and with no host a specimen is Human. In "pork blood sausage" it is an ingredient."""
    sausage = _classify(lexicon, "pork blood sausage")
    assert (sausage.classify_status, sausage.classify_reason, sausage.lexicon_suggestion) == ("review", "specimen_or_food", "Food|Human")
    # A host, or "human" in the text, still says whose specimen it is.
    assert _classify(lexicon, "pork blood sausage", "Homo sapiens").classify_reason == "specimens+human-host"
    assert _classify(lexicon, "stool from a food handler").classify_reason == "specimens+food-handler"
    # A food word that is part of the specimen's own name is the specimen.
    assert _classify(lexicon, "breast milk").classify_reason == "specimens+no_host"
    assert _classify(lexicon, "blood").classify_reason == "specimens+no_host"


def test_accents_do_not_split_words(lexicon):
    assert tokenize("Pâté de foie, Hôpital Saint-André") == ["pate", "de", "foie", "hopital", "saint", "andre"]
    assert tokenize("jalapeño") == ["jalapeno"]
    assert _classify(lexicon, "jalapeño").classify_reason == "food-items"


def test_the_new_rule_keys(tmp_path: Path):
    rules = (
        'version: "mine"\nrules:\n'
        "  - id: place\n    phrases: [farm]\n    l1: Environmental\n    priority: 70\n    yields_to_host: true\n"
        "  - id: waste\n    phrases: [manure]\n    l1: Environmental\n    priority: 60\n    setting_of: animal\n"
        "  - id: meat\n    exact: [chicken]\n    l1: Food\n    priority: 45\n    no_host: true\n"
        "  - id: animals\n    fields: [host, isolation_source]\n    phrases: [chicken, pig]\n    l1: Animal\n    priority: 40\n"
    )
    path = tmp_path / "lex.yaml"
    path.write_text(rules, encoding="utf-8")
    lexicon = load_lexicon(path)
    rows = classify_rows([
        {"isolation_source": "farm", "host": ""}, {"isolation_source": "farm", "host": "pig"},
        {"isolation_source": "manure", "host": "pig"}, {"isolation_source": "chicken", "host": ""},
        {"isolation_source": "chicken", "host": "chicken"},
    ], lexicon)
    assert [(r["l1"], r["classify_reason"]) for r in rows] == [
        ("Environmental", "place"), ("Animal", "place+animals"), ("Environmental", "waste"), ("Food", "meat"), ("Animal", "animals")]
    path.write_text(rules.replace("setting_of: animal", "setting_of: human"), encoding="utf-8")
    with pytest.raises(ValueError, match="setting_of must be animal"):
        load_lexicon(path)
