"""Hosts the word lists do not name are recognised by the names NCBI Taxonomy gives to animals (lexicon v6)."""

import gzip
from importlib.resources import files
from pathlib import Path

import pytest

from onehealth.classify import classify_record, classify_rows
from onehealth.io_tsv import read_tsv
from onehealth.lexicon import load_lexicon
from onehealth.taxonomy import taxon_kind
from onehealth.validate import read_truth, validate

IMP = Path(__file__).parent / "data" / "imp_carbapenemase_sources" / "samples.tsv"


def _classify(lexicon, source="", host="", **more):
    return classify_record({"isolation_source": source, "host": host, **more}, lexicon)


@pytest.mark.parametrize(
    "value",
    [
        "Larus dominicanus",  # a species
        "Macaca",  # a genus on its own
        "Cervidae",  # a family
        "Aves",  # a class
        "Bubalus bubalis",
        "Bos primigenius taurus",  # a species with a subspecies or breed
        "Sus scrofa domesticus",
        "Meleagris sp.",  # a genus with sp.
        "Corvus spp",
        "gull",  # a common name; NCBI lists "gulls"
        "Gulls",
        "sika deer",
        "Rhesus macaque",
        "American crow",
        "mud crab",  # common names cover animals that are not chordates
        "  SIKA_DEER ",
    ],
)
def test_names_of_animals_are_recognised(value):
    assert taxon_kind(value) == "animal"


@pytest.mark.parametrize("value", ["Homo", "Homo sapiens", "Homo sapiens sapiens", "human", "Humans"])
def test_names_of_people_are_recognised(value):
    assert taxon_kind(value) == "human"


@pytest.mark.parametrize(
    "value",
    [
        "",
        None,
        "sea water",  # Sea is a genus of fish
        "colon cancer",  # Colon and Cancer are genera of animals that are not chordates
        "sepsis",  # a genus of flies
        "cancer patient",
        "Environment",
        "Water",
        "patient",
        "old lady",
        "soldier",
        "pure culture",
        "white stork nestling",  # an animal, but not a name: left for review, never guessed from a word inside
        "Zea mays",  # a plant
        "Allium cepa",
        "Escherichia coli",  # a bacterium
        "Salmonella enterica",
        "environmental samples",  # NCBI's name for a group that is not an organism
        "unclassified Aves",
        "not collected",
    ],
)
def test_other_values_are_not_names_of_animals_or_people(value):
    assert taxon_kind(value) == ""


@pytest.mark.parametrize(
    "source, host, l1, status, reason",
    [
        ("", "Gull", "Animal", "ok", "host-taxonomy"),
        ("", "Bubalus bubalis", "Animal", "ok", "host-taxonomy"),
        ("feces", "Macaca mulatta", "Animal", "ok", "specimens+host-taxonomy"),  # was host_unrecognised
        ("wound", "sika deer", "Animal", "ok", "specimens+host-taxonomy"),
        ("blood", "Homo", "Human", "ok", "specimens+host-taxonomy-human"),
        # A host the word lists already know keeps its reason.
        ("feces", "Bos taurus", "Animal", "ok", "specimens+animal-hosts"),
        ("", "chicken", "Animal", "ok", "animal-hosts"),
        ("", "Homo sapiens", "Human", "ok", "human-host"),
        # An animal host with a source that says otherwise is a conflict, as for a host on the word list.
        ("river water", "gull", "Unclassified", "review", "host_source_conflict"),
        ("chicken breast", "gull", "Unclassified", "review", "host_source_conflict"),
        # Not an animal: as before.
        ("urine", "Zea mays", "Unclassified", "review", "host_unrecognised"),
        ("", "Environment", "Unclassified", "review", "unmatched"),
        ("", "sea water", "Unclassified", "review", "unmatched"),
        ("feces", "white stork nestling", "Unclassified", "review", "host_unrecognised"),
    ],
)
def test_classify_places_hosts_by_taxonomy(lexicon, source, host, l1, status, reason):
    result = _classify(lexicon, source, host)
    assert (result.l1, result.classify_status, result.classify_reason) == (l1, status, reason)


def test_a_species_name_in_isolation_source_is_still_not_a_description_of_the_sample(lexicon):
    """The taxonomy rules read host only. "Larus dominicanus" typed into isolation_source stays in review."""
    result = _classify(lexicon, "Larus dominicanus")
    assert (result.classify_status, result.classify_reason) == ("review", "taxon_only_source")
    assert _classify(lexicon, "gull").classify_reason == "unmatched"


def test_taxonomy_is_a_rule_so_a_lexicon_without_it_behaves_as_before(tmp_path: Path):
    rules = (
        "rules:\n"
        "  - id: specimens\n    phrases: [feces]\n    l1: Human\n    host_dependent: true\n    priority: 80\n"
        "  - id: animal-hosts\n    fields: [host]\n    phrases: [chicken]\n    l1: Animal\n    priority: 40\n"
    )
    without = tmp_path / "without.yaml"
    without.write_text('version: "mine"\n' + rules, encoding="utf-8")
    assert _classify(load_lexicon(without), "feces", "Gull").classify_reason == "host_unrecognised"
    with_rule = tmp_path / "with.yaml"
    with_rule.write_text('version: "mine"\n' + rules + "  - id: any-animal\n    fields: [host]\n    taxonomy: animal\n"
                         "    l1: Animal\n    priority: 40\n", encoding="utf-8")
    result = _classify(load_lexicon(with_rule), "feces", "Gull")
    assert (result.l1, result.classify_reason) == ("Animal", "specimens+any-animal")


@pytest.mark.parametrize("value", ["plant", "yes", "true", "animals"])
def test_an_unknown_taxonomy_value_is_rejected(tmp_path: Path, value):
    path = tmp_path / "lex.yaml"
    path.write_text(f'version: "x"\nrules:\n  - id: r\n    fields: [host]\n    taxonomy: {value}\n    l1: Animal\n', encoding="utf-8")
    with pytest.raises(ValueError, match="taxonomy must be one of animal, human"):
        load_lexicon(path)


def test_the_name_table_holds_animals_only():
    scientific, common = set(), set()
    path = files("onehealth").joinpath("data/animal_names.tsv.gz")
    with path.open("rb") as raw, gzip.open(raw, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert lines[0].startswith("# Names NCBI Taxonomy gives to animals") and "md5" in lines[1]
    for line in lines:
        if not line.startswith("#"):
            kind, name = line.split("\t")
            (scientific if kind == "s" else common).add(name)
    assert len(scientific) > 80_000 and len(common) > 30_000
    names = scientific | common
    # People are left to the lexicon's own rules, and groups that are not organisms are left out.
    assert not names & {"homo", "homo sapiens", "human", "humans", "environmental samples"}
    assert not [n for n in names if {"unclassified", "environmental", "uncultured"} & set(n.split())]
    # Words a host field can hold that are also genera of animals that are not chordates.
    assert not names & {"sepsis", "cancer", "colon", "trachea", "tibia", "water", "environment", "food", "patient"}
    assert all(name == name.lower().strip() for name in names)


def test_imp_benchmark_animal_hosts_are_all_found(lexicon):
    """13 animal isolates (American crow, a silver gull) had hosts no word list named; every curated animal is now Animal."""
    result = validate(classify_rows(read_tsv(IMP), lexicon), read_truth(IMP, "truth", "paper_rollup"))
    assert result.count(truth="Animal") == result.count("Animal", "Animal") == 51
    assert (round(result.coverage, 4), round(result.accuracy, 4)) == (0.8485, 0.9976)


# --- "human" in a specimen's own text ---------------------------------------------------------------


@pytest.mark.parametrize(
    "source, host, more, l1, reason",
    [
        ("human feces", "", {}, "Human", "specimens+human-word"),  # was specimens+no_host
        ("Human, blood culture", "", {}, "Human", "specimens+human-word"),
        # Filed as environmental, food or other: the text still says whose specimen it is.
        ("human feces", "", {"package": "Pathogen.env.1.0"}, "Human", "specimens+human-word"),
        ("feces", "", {"package": "Pathogen.env.1.0"}, "Unclassified", "no_host_env_package"),
        ("human feces", "Homo sapiens", {}, "Human", "specimens+human-host"),
        ("human stool", "Sus scrofa", {}, "Unclassified", "host_source_conflict"),
        ("pig farm worker human stool", "", {}, "Unclassified", "host_source_conflict"),
        # On its own the word decides nothing.
        ("human contact surface", "", {}, "Environmental", "environmental-sampling"),
        ("human listeriosis", "", {}, "Unclassified", "unmatched"),
        ("human sewage", "", {}, "Environmental", "environment-outdoor"),
        ("human food", "", {}, "Food", "food-products"),
        # A veterinary record has patients too, so "patient" is not read this way.
        ("patient feces", "", {"package": "Pathogen.env.1.0"}, "Unclassified", "no_host_env_package"),
    ],
)
def test_human_in_a_specimens_text_says_whose_it_is(lexicon, source, host, more, l1, reason):
    result = _classify(lexicon, source, host, **more)
    assert (result.l1, result.classify_reason) == (l1, reason)
