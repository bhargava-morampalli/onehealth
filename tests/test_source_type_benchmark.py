"""The classifier against NCBI's own source_type (see tests/data/ncbi_pd_source_types/README.md).

`source_type` is hidden from the classifier and used as the answer. The held-out split was
not looked at while lexicon v5 was written, so its numbers are the ones to quote.
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
    assert by_split == {"development": 33372, "heldout": 39722, "validation": 32208}
    assert {row["ncbi_source_type"] for row in records} == set(LABEL)
    # Every record has something to classify from, and the answer is not in a column classify reads.
    assert all(row["isolation_source"] or row["host"] for row in records)
    assert "source_type" not in records[0] and "sample_type" not in records[0]
    groups = {split: {row["organism_group"] for row in records if row["split"] == split} for split in by_split}
    assert groups["heldout"] - groups["development"] == {
        "Staphylococcus_aureus", "Vibrio_parahaemolyticus", "Bacillus_cereus_group", "Cronobacter"}
    assert groups["validation"] - groups["development"] - groups["heldout"] == {
        "Clostridium_perfringens", "Vibrio_cholerae", "Yersinia_enterocolitica"}


def test_validation_baseline(records, lexicon):
    """Floors just under lexicon v7 on records nobody looked at while it was written: 95.9% classified,
    96.0% of those right. This is the figure to quote for v7.

    The validation split was read once, after v7 was finished, and v7 was not changed afterwards.
    Before v7 the same records scored 95.5% and 94.8% (Animal precision 82.0%, Environmental recall 73.6%).
    """
    result = _score([r for r in records if r["split"] == "validation"], lexicon)
    assert result["coverage"] >= 0.955
    assert result["accuracy"] >= 0.955
    assert result["precision"]["Human"] >= 0.985 and result["precision"]["Environmental"] >= 0.985
    assert result["precision"]["Food"] >= 0.98
    assert result["precision"]["Animal"] >= 0.855
    assert result["recall"]["Environmental"] >= 0.84
    assert result["recall"]["Animal"] >= 0.895


def test_heldout_baseline(records, lexicon):
    """Floors just under lexicon v7 on the held-out records: 90.3% classified, 96.1% of those right.

    Lexicon v4, the last one before this benchmark existed, classified 76.9% and got 80.9% of
    those right: fewer than half of its Animal calls agreed with NCBI (45.4%), because it called
    retail meat Animal. v5 fixed that (88.7% and 93.9%), v6 found more animal hosts, and v7
    follows NCBI on farm settings and on "chicken" (Animal precision 74.3% to 84.2%, Environmental
    recall 77.6% to 86.9%). v7 was checked against this split, so it is no longer unseen; the
    validation split is. Animal recall is a little lower than in v6 (91.7% to 90.8%): a minority
    of submitters label manure and a bare "chicken" as animal.
    """
    result = _score([r for r in records if r["split"] == "heldout"], lexicon)
    assert result["coverage"] >= 0.90
    assert result["accuracy"] >= 0.955
    assert result["precision"]["Human"] >= 0.975 and result["precision"]["Environmental"] >= 0.975
    assert result["precision"]["Food"] >= 0.975
    assert result["precision"]["Animal"] >= 0.835
    assert result["recall"]["Animal"] >= 0.90
    assert result["recall"]["Food"] >= 0.73
    assert result["recall"]["Environmental"] >= 0.86


def test_development_baseline(records, lexicon):
    result = _score([r for r in records if r["split"] == "development"], lexicon)
    assert result["coverage"] >= 0.955 and result["accuracy"] >= 0.97
