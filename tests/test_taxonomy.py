"""Hosts the word lists do not name are recognised by the names NCBI Taxonomy gives to animals (lexicon v6)."""

import gzip
import importlib.util
import io
import re
import tarfile
import time
from importlib.resources import files
from pathlib import Path

import pytest

from onehealth import cli, taxonomy
from onehealth.classify import classify_record, classify_rows
from onehealth.io_tsv import read_tsv
from onehealth.lexicon import load_lexicon
from onehealth.taxonomy import animal_name_within, taxon_kind
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
        "nestling of a white stork in a zoo",  # the name is in the middle, not at the end
        "bologna, hot dogs, sausage-like products",  # a list of things
        "kitchen sponge",  # one word at the end must name a vertebrate
        "sponge cake",
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
        # v8: a name at the end of the value, or before a word for the animal's age.
        ("feces", "white stork nestling", "Animal", "ok", "specimens+host-taxonomy"),
        ("feces", "European rook", "Animal", "ok", "specimens+host-taxonomy"),
        ("", "free-range white-tailed deer", "Animal", "ok", "host-taxonomy"),
        # v8: an oyster, a shrimp or a fly by its species name.
        ("tissue", "Magallana gigas", "Animal", "ok", "specimens+host-taxonomy"),
        ("", "Litopenaeus vannamei", "Animal", "ok", "host-taxonomy"),
        ("", "Musca domestica", "Animal", "ok", "host-taxonomy"),
        # A list of several things is not one host.
        ("feces", "bologna, hot dogs, sausage-like products", "Unclassified", "review", "host_unrecognised"),
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
    assert len(scientific) > 200_000 and len(common) > 30_000
    assert {line.split("\t")[0] for line in lines if not line.startswith("#")} == {"s", "v", "c"}
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
    assert (round(result.coverage, 4), round(result.accuracy, 4)) == (0.8489, 0.9976)


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
        ("human contact", "", {}, "Unclassified", "unmatched"),
        ("human sewage", "", {}, "Environmental", "environment-outdoor"),
        ("human food", "", {}, "Food", "food-products"),
        # A veterinary record has patients too, so "patient" is not read this way.
        ("patient feces", "", {"package": "Pathogen.env.1.0"}, "Unclassified", "no_host_env_package"),
    ],
)
def test_human_in_a_specimens_text_says_whose_it_is(lexicon, source, host, more, l1, reason):
    result = _classify(lexicon, source, host, **more)
    assert (result.l1, result.classify_reason) == (l1, reason)


# --- v8: a name at the end of a host value, and inside a specimen's text ------------------------------


@pytest.mark.parametrize(
    "value",
    ["white stork nestling", "European rook", "ground squirrel", "forest musk deer", "wild seal",
     "free-range white-tailed deer", "others - Daboia russelii", "adult female macaque",
     "Magallana gigas", "Crassostrea virginica", "Procambarus clarkii", "Litopenaeus vannamei", "Musca domestica"],
)
def test_a_name_at_the_end_of_a_host_value_and_invertebrate_species(value):
    assert taxon_kind(value) == "animal"


@pytest.mark.parametrize(
    "text, name",
    [
        ("stool (Odocoileus hemionus)", "odocoileus hemionus"),
        ("feces (Ursus americanus)", "ursus americanus"),
        ("giant panda feces", "giant panda"),
        ("tibetan antelope feces", "tibetan antelope"),
        ("pigeon stool", "pigeon"),
        ("roe deer faeces", "deer"),
        # Ordinary words that are also animals' names.
        ("catheter tip", ""),  # a moth
        ("white blood cells", ""),  # a butterfly
        ("wound, sole of foot", ""),  # a fish
        ("ECG monitor swab", ""),  # a lizard
        ("tick bite wound", ""),
        ("sponge", ""),
        ("maternal placental", ""),  # "placentals", a clade of mammals
        ("blood", ""),
        ("rectal swab", ""),
        ("", ""),
    ],
)
def test_an_animal_name_inside_a_text(text, name):
    assert animal_name_within(text) == name


def test_an_animal_named_in_a_specimen_text_decides_whose_it_is(lexicon):
    """Until v8 these were Human: the word lists do not know the animal, and with no host a specimen is a person's."""
    for source in ("stool (Odocoileus hemionus)", "pigeon stool", "giant panda feces", "rectal swab (Ursus americanus)"):
        result = _classify(lexicon, source, "")
        assert (result.l1, result.classify_reason) == ("Animal", "specimens+animal-name-in-specimen"), source
    # Only for a specimen, and only with no host.
    assert _classify(lexicon, "stool", "").classify_reason == "specimens+no_host"
    assert _classify(lexicon, "catheter tip", "").classify_reason == "specimens+no_host"
    assert _classify(lexicon, "giant panda feces", "Homo sapiens").classify_reason == "specimens+human-host"
    assert _classify(lexicon, "giant panda enclosure", "").classify_status == "review"
    # A bare species name in isolation_source is still not a description of a sample.
    assert _classify(lexicon, "Bos taurus", "").classify_reason == "taxon_only_source"


def test_within_needs_taxonomy_animal(tmp_path: Path):
    path = tmp_path / "lex.yaml"
    path.write_text('version: "x"\nrules:\n  - id: r\n    within: true\n    phrases: [a]\n    l1: Animal\n', encoding="utf-8")
    with pytest.raises(ValueError, match="within needs taxonomy: animal"):
        load_lexicon(path)


# --- the table, its date, and refreshing it -----------------------------------------------------------

TOOLS = Path(__file__).resolve().parent.parent / "tools"


def _tool(name: str):
    spec = importlib.util.spec_from_file_location(f"tool_{name}", TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _dmp(rows: list[tuple[str, ...]]) -> bytes:
    """NCBI's .dmp format: fields separated by <TAB>|<TAB>, each line ending with <TAB>|."""
    return "".join("\t|\t".join(row) + "\t|\n" for row in rows).encode("utf-8")


def _dump(names: list[tuple[str, str, str]], when: float = 1_763_078_400.0) -> bytes:
    """A new_taxdump.tar.gz with names.dmp and rankedlineage.dmp only, written at `when` (2025-11-14 UTC)."""
    # tax_id, name, species, genus, family, order, class, phylum, kingdom, domain
    lineage = [
        ("1", "Larus", "", "Larus", "Laridae", "Charadriiformes", "Aves", "Chordata", "Metazoa", "Eukaryota"),
        ("2", "Homo sapiens", "", "Homo", "Hominidae", "Primates", "Mammalia", "Chordata", "Metazoa", "Eukaryota"),
        ("3", "Magallana gigas", "", "Magallana", "Ostreidae", "Ostreida", "Bivalvia", "Mollusca", "Metazoa", "Eukaryota"),
        ("4", "Magallana", "", "", "Ostreidae", "Ostreida", "Bivalvia", "Mollusca", "Metazoa", "Eukaryota"),
        ("5", "Musca domestica", "", "Musca", "Muscidae", "Diptera", "Insecta", "Arthropoda", "Metazoa", "Eukaryota"),
        ("6", "Sepsis fulgens", "", "Sepsis", "Sepsidae", "Diptera", "Insecta", "Arthropoda", "Metazoa", "Eukaryota"),
        ("7", "Penaeus vannamei", "", "Penaeus", "Penaeidae", "Decapoda", "Malacostraca", "Arthropoda", "Metazoa", "Eukaryota"),
        ("8", "Zea mays", "", "Zea", "Poaceae", "Poales", "", "Streptophyta", "Viridiplantae", "Eukaryota"),
        ("9", "Rattus sp. ABC-2011", "", "Rattus", "Muridae", "Rodentia", "Mammalia", "Chordata", "Metazoa", "Eukaryota"),
    ]
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w:gz") as tar:
        for member, data in (("names.dmp", _dmp([(t, n, "", c) for t, n, c in names])), ("rankedlineage.dmp", _dmp(lineage))):
            info = tarfile.TarInfo(member)
            info.size, info.mtime = len(data), when
            tar.addfile(info, io.BytesIO(data))
    return archive.getvalue()


NAMES = [
    ("1", "Larus", "scientific name"), ("1", "gull", "genbank common name"),
    ("2", "Homo sapiens", "scientific name"), ("2", "human", "genbank common name"),
    ("3", "Magallana gigas", "scientific name"), ("3", "Crassostrea gigas", "synonym"),
    ("4", "Magallana", "scientific name"),
    ("5", "Musca domestica", "scientific name"), ("5", "house fly", "genbank common name"),
    ("6", "Sepsis fulgens", "scientific name"), ("6", "Sepsis", "synonym"),
    ("7", "Penaeus vannamei", "scientific name"), ("7", "Litopenaeus vannamei", "synonym"),
    ("8", "Zea mays", "scientific name"), ("8", "maize", "common name"),
    ("9", "Rattus sp. ABC-2011", "scientific name"),
]


def test_the_builder_keeps_the_names_it_should():
    builder = _tool("build_animal_names")
    scientific, vertebrate, common = builder.build(_dump(NAMES))
    # A chordate at any rank; another animal by a species name, if it is a mollusc, a crustacean or
    # has a common name; synonyms of those. Not people, plants, unnamed entries or a lone genus.
    assert scientific == {"larus", "magallana gigas", "crassostrea gigas", "musca domestica", "penaeus vannamei",
                          "litopenaeus vannamei"}
    assert vertebrate == {"gull"} and common == {"house fly"}
    assert builder.dump_date(_dump(NAMES)) == "2025-11-14"


def test_a_table_is_rebuilt_byte_for_byte_and_records_its_date(tmp_path):
    builder, refresh = _tool("build_animal_names"), _tool("refresh_animal_names")
    dump = tmp_path / "new_taxdump.tar.gz"
    dump.write_bytes(_dump(NAMES))
    first, second = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    assert builder.main([str(dump), "--out", str(first)]) == 0
    time.sleep(1.1)  # a gzip time stamp would differ now
    assert builder.main([str(dump), "--out", str(second)]) == 0
    assert first.read_bytes() == second.read_bytes()
    names, date = refresh.read_table(first)
    assert date == "2025-11-14" and names["v"] == {"gull"} and len(names["s"]) == 6


def _shipped_date() -> str:
    """The date line of the table in the package, read here without onehealth's own code.

    Not written into the tests: the monthly refresh replaces the table, and with it the date.
    """
    path = files("onehealth").joinpath("data/animal_names.tsv.gz")
    with path.open("rb") as raw, gzip.open(raw, "rt", encoding="utf-8") as handle:
        dated = [line for line in handle.read(4000).splitlines() if line.startswith("# NCBI Taxonomy of ")]
    assert len(dated) == 1, "the table must record the date of its NCBI dump"
    return dated[0].split()[4]


def test_the_shipped_table_says_which_dump_it_is_from_and_classify_prints_it(tmp_path, capsys):
    shipped = _shipped_date()
    assert re.fullmatch(r"20\d\d-\d\d-\d\d", shipped) and shipped >= "2026-11-14"
    assert taxonomy.table_date() == shipped
    table = tmp_path / "t.tsv"
    table.write_text("query_id\thost\tisolation_source\nA\tgull\tfeces\n", encoding="utf-8")
    assert cli.main(["classify", str(table), "-o", str(tmp_path / "out.tsv")]) == 0
    assert f"animal names from NCBI Taxonomy of {shipped}" in capsys.readouterr().err


def test_a_refresh_compares_the_names_and_the_classes_and_puts_the_shipped_table_back(tmp_path, lexicon):
    builder, refresh = _tool("build_animal_names"), _tool("refresh_animal_names")
    without = [row for row in NAMES if row[1] != "gull"]
    old, new = tmp_path / "old.tsv.gz", tmp_path / "new.tsv.gz"
    old.write_bytes(builder.render(builder.build(_dump(without)), "old", "0", "2025-11-01"))
    new.write_bytes(builder.render(builder.build(_dump(NAMES)), "new", "1", "2025-11-14"))

    names = refresh.compare(refresh.read_table(old)[0], refresh.read_table(new)[0])
    assert names == {"s": ([], []), "v": (["gull"], []), "c": ([], [])}

    rows = [{"isolation_source": "feces", "host": "gull", "label": "Animal"},
            {"isolation_source": "blood", "host": "Homo sapiens", "label": "Human"}]
    data = [("tiny", rows, lambda row: 3, lambda row: row["label"])]
    (result,) = refresh.class_changes(old, new, lexicon, data)
    assert result["rows"] == 6 and result["toward"] == 3 and result["away"] == 0
    assert result["changes"] == {("feces", "gull", "Unclassified", "Animal"): 3}

    text = refresh.report("2025-11-01", "2025-11-14", names, [result])
    assert "Animal names from NCBI Taxonomy of 2025-11-14 (was 2025-11-01)" in text
    assert "| common names of chordates | 1 | 0 |" in text and "`feces`, host `gull`: Unclassified to Animal" in text
    assert "If any class changes above, also raise the package version" in text
    entry = refresh.changelog_entry("2025-11-01", "2025-11-14", names, [result])
    assert entry.startswith("## 14 November 2025: animal names from NCBI Taxonomy of 2025-11-14\n")
    assert "1 names added and 0 removed" in entry and "3 records change class" in entry
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text("# Changelog\n\nIntro.\n\n## Older\n\nText.\n", encoding="utf-8")
    refresh.add_to_changelog(changelog, entry)
    assert changelog.read_text(encoding="utf-8") == "# Changelog\n\nIntro.\n\n" + entry + "## Older\n\nText.\n"
    # The shipped table is in use again.
    assert taxonomy.table_date() == _shipped_date() and taxon_kind("European rook") == "animal"


def test_a_refresh_with_the_same_names_changes_nothing(tmp_path, monkeypatch):
    builder, refresh = _tool("build_animal_names"), _tool("refresh_animal_names")
    dump = tmp_path / "new_taxdump.tar.gz"
    dump.write_bytes(_dump(NAMES))
    table = tmp_path / "animal_names.tsv.gz"
    table.write_bytes(builder.render(builder.build(_dump(NAMES)), "earlier", "0", "2025-11-01"))
    before = table.read_bytes()
    monkeypatch.setattr(refresh, "TABLE", table)
    outputs = tmp_path / "github_output"
    monkeypatch.setenv("GITHUB_OUTPUT", str(outputs))
    assert refresh.main([str(dump), "--summary", str(tmp_path / "report.md")]) == 0
    assert table.read_bytes() == before  # a newer dump with the same names is not a change
    assert outputs.read_text(encoding="utf-8") == "changed=false\ndate=2025-11-14\n"
    assert "Nothing to do" in (tmp_path / "report.md").read_text(encoding="utf-8")
