"""The classifier against NCBI's own source_type (see tests/data/ncbi_pd_source_types/README.md).

`source_type` is hidden from the classifier and used as the answer. Each lexicon version was
written with the earlier splits in view, so the figure to quote for a version is the one from
the split that was drawn after it was finished: the fresh split, for v8.
"""

from collections import Counter
from pathlib import Path

import pytest

from onehealth.classify import classify_record
from onehealth.io_tsv import read_tsv

RECORDS = Path(__file__).parent / "data" / "ncbi_pd_source_types" / "records.tsv"
LABEL = {"human": "Human", "animal": "Animal", "food": "Food", "environmental": "Environmental"}


@pytest.fixture(scope="module")
def records():
    return read_tsv(RECORDS)


def _score(rows, lexicon):
    """Coverage, accuracy of what was classified, and precision and recall per class, weighted by records."""
    confusion: Counter = Counter()
    for row in rows:
        evidence = {name: row[name] for name in ("isolation_source", "host", "host_disease")}
        confusion[(LABEL[row["ncbi_source_type"]], classify_record(evidence, lexicon).l1)] += int(row["records"])
    total = sum(confusion.values())
    classified = sum(n for (_, predicted), n in confusion.items() if predicted != "Unclassified")
    correct = sum(n for (truth, predicted), n in confusion.items() if truth == predicted)

    def share(label, of):
        return confusion[(label, label)] / sum(n for pair, n in confusion.items() if pair[of] == label)

    return {
        "records": total,
        "coverage": classified / total,
        "accuracy": correct / classified,
        "precision": {label: share(label, 1) for label in LABEL.values()},
        "recall": {label: share(label, 0) for label in LABEL.values()},
    }


def test_the_fixture_is_what_the_readme_says(records):
    by_split = Counter()
    for row in records:
        by_split[row["split"]] += int(row["records"])
    assert by_split == {"development": 33372, "heldout": 39722, "validation": 32208, "final": 34252, "fresh": 26604}
    assert sum(by_split.values()) == 166158 and sum(n for split, n in by_split.items() if split != "fresh") == 139554
    assert {row["ncbi_source_type"] for row in records} == set(LABEL)
    # Every record has something to classify from, and the answer is not in a column classify reads.
    assert all(row["isolation_source"] or row["host"] for row in records)
    assert "source_type" not in records[0] and "sample_type" not in records[0]
    groups = {split: {row["organism_group"] for row in records if row["split"] == split} for split in by_split}
    assert groups["heldout"] - groups["development"] == {
        "Staphylococcus_aureus", "Vibrio_parahaemolyticus", "Bacillus_cereus_group", "Cronobacter"}
    assert groups["validation"] - groups["development"] - groups["heldout"] == {
        "Clostridium_perfringens", "Vibrio_cholerae", "Yersinia_enterocolitica"}
    assert groups["final"] - groups["development"] - groups["heldout"] - groups["validation"] == {
        "Vibrio_vulnificus", "Listeria_innocua", "Staphylococcus_pseudintermedius", "Clostridium_botulinum", "Aeromonas"}
    # Listeria and Acinetobacter had no unread part left for the fresh split.
    assert groups["fresh"] == {
        "Salmonella", "Escherichia_coli_Shigella", "Campylobacter", "Klebsiella",
        "Vibrio_alginolyticus", "Enterococcus_faecalis", "Enterococcus_faecium", "Citrobacter_freundii",
        "Aeromonas_veronii", "Pseudomonas_aeruginosa", "Enterobacter_hormaechei", "Serratia"}


def test_fresh_baseline(records, lexicon):
    """Floors just under the shipped lexicon on the fresh split: 97.8% classified, 97.9% of those right.

    The figure to quote for v8 is a little lower: 97.8% and 96.7%. That is what v8 scored when the
    split was read, once, after v8 was finished (v7: 97.1% and 95.4%; Animal precision 90.7% to
    97.4%, Food recall 91.7% to 97.4%). The reading showed one thing: 296 records described only
    as "healthcare facility", human to their submitter, which v7 had taken for a hospital fixture
    with no labelled record to go by. That one value was corrected afterwards, so the split is no
    longer unseen and the 97.9% is not a fair test. Nothing else was changed.

    Animal recall is 92.0%: the records missed name the animal in isolation_source and give no
    host ("Larus cachinnans", "gull", "fox"), or are faeces with no host at all.
    """
    result = _score([r for r in records if r["split"] == "fresh"], lexicon)
    assert result["records"] == 26604
    assert result["coverage"] >= 0.975
    assert result["accuracy"] >= 0.975
    assert result["precision"]["Human"] >= 0.985 and result["precision"]["Food"] >= 0.98
    assert result["precision"]["Animal"] >= 0.97
    assert result["precision"]["Environmental"] >= 0.935
    assert result["recall"]["Human"] >= 0.98 and result["recall"]["Food"] >= 0.97
    assert result["recall"]["Animal"] >= 0.915
    assert result["recall"]["Environmental"] >= 0.845


def test_final_baseline(records, lexicon):
    """Floors just under lexicon v8 on the final split: 98.3% classified, 98.3% of those right.

    The final split was read once, after v7 was finished: v7 scored 95.5% and 96.5% (v6: 94.6% and
    95.3%), and that was the figure to quote for v7. v8 was written with it in view, so for v8
    it is a seen split: Food recall 83.3% to 94.4%, Animal precision 91.7% to 97.2%.
    """
    result = _score([r for r in records if r["split"] == "final"], lexicon)
    assert result["coverage"] >= 0.98
    assert result["accuracy"] >= 0.98
    assert result["precision"]["Human"] >= 0.99 and result["precision"]["Food"] >= 0.98
    assert result["precision"]["Environmental"] >= 0.97
    assert result["precision"]["Animal"] >= 0.965
    assert result["recall"]["Animal"] >= 0.97
    assert result["recall"]["Environmental"] >= 0.945
    assert result["recall"]["Food"] >= 0.94


def test_validation_baseline(records, lexicon):
    """Floors just under lexicon v8 on the validation split: 98.8% classified, 98.1% of those right.

    This split was meant to be read once, for v7. It showed that one v7 rule went against NCBI's
    records (a farm with an animal host is the animal, not a conflict), and the rule was corrected.
    v7 scored 96.7% and 96.1%.
    """
    result = _score([r for r in records if r["split"] == "validation"], lexicon)
    assert result["coverage"] >= 0.985
    assert result["accuracy"] >= 0.98
    assert result["precision"]["Human"] >= 0.99 and result["precision"]["Environmental"] >= 0.99
    assert result["precision"]["Food"] >= 0.99
    assert result["precision"]["Animal"] >= 0.94
    assert result["recall"]["Environmental"] >= 0.89
    assert result["recall"]["Animal"] >= 0.97
    assert result["recall"]["Food"] >= 0.955


def test_heldout_baseline(records, lexicon):
    """Floors just under lexicon v8 on the held-out records: 97.1% classified, 97.9% of those right.

    Lexicon v4, the last one before this benchmark existed, classified 76.9% and got 80.9% of
    those right: fewer than half of its Animal calls agreed with NCBI (45.4%), because it called
    retail meat Animal. v5 fixed that (88.7% and 93.9%), v6 found more animal hosts, and v7
    followed NCBI on farm settings and on "chicken" (90.3% and 96.1%). v8 does the same for every
    meat species and names the foods v7 left in review: Food recall 73.9% to 92.9%, Animal
    precision 84.2% to 92.4%.
    """
    result = _score([r for r in records if r["split"] == "heldout"], lexicon)
    assert result["coverage"] >= 0.97
    assert result["accuracy"] >= 0.975
    assert result["precision"]["Human"] >= 0.99 and result["precision"]["Environmental"] >= 0.98
    assert result["precision"]["Food"] >= 0.98
    assert result["precision"]["Animal"] >= 0.92
    assert result["recall"]["Animal"] >= 0.95
    assert result["recall"]["Food"] >= 0.925
    assert result["recall"]["Environmental"] >= 0.915


def test_development_baseline(records, lexicon):
    result = _score([r for r in records if r["split"] == "development"], lexicon)
    assert result["coverage"] >= 0.975 and result["accuracy"] >= 0.98


def test_an_animal_named_alone_in_isolation_source_is_the_documented_known_limitation(records, lexicon):
    """documentation/commands/classify.md: with no host, a species name in isolation_source goes to review."""
    from onehealth.taxonomy import taxon_kind
    from onehealth.text import is_unknown_value

    animal = alone = 0
    for row in records:
        if row["split"] != "fresh" or row["ncbi_source_type"] != "animal":
            continue
        animal += int(row["records"])
        result = classify_record({name: row[name] for name in ("isolation_source", "host", "host_disease")}, lexicon)
        if result.classify_status == "review" and is_unknown_value(row["host"]) and taxon_kind(row["isolation_source"]) == "animal":
            alone += int(row["records"])
    assert (alone, animal) == (151, 5697)
    page = (Path(__file__).parents[1] / "documentation" / "commands" / "classify.md").read_text(encoding="utf-8")
    assert f"{alone} of {animal:,} animal records to review" in page
