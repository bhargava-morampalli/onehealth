"""The classifier against NCBI's own source_type (see tests/data/ncbi_pd_source_types/README.md).

`source_type` is hidden from the classifier and used as the answer. The held-out split was
not looked at while lexicon v5 was written, so its numbers are the ones to quote.
"""

from collections import Counter
from pathlib import Path

import pytest

from onehealth.classify import classify_record
from onehealth.io_tsv import read_tsv
from onehealth.lexicon import load_lexicon

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
    assert by_split == {"development": 33372, "heldout": 39722}
    assert {row["ncbi_source_type"] for row in records} == set(LABEL)
    # Every record has something to classify from, and the answer is not in a column classify reads.
    assert all(row["isolation_source"] or row["host"] for row in records)
    assert "source_type" not in records[0] and "sample_type" not in records[0]
    heldout_only = {row["organism_group"] for row in records if row["split"] == "heldout"} - {
        row["organism_group"] for row in records if row["split"] == "development"
    }
    assert heldout_only == {"Staphylococcus_aureus", "Vibrio_parahaemolyticus", "Bacillus_cereus_group", "Cronobacter"}


def test_heldout_baseline(records, lexicon):
    """Floors just under lexicon v5 on records nobody looked at: 88.7% classified, 93.9% of those right."""
    result = _score([r for r in records if r["split"] == "heldout"], lexicon)
    assert result["coverage"] >= 0.885
    assert result["accuracy"] >= 0.935
    assert result["precision"]["Human"] >= 0.975 and result["precision"]["Environmental"] >= 0.985
    assert result["precision"]["Food"] >= 0.97
    # Still the weak class: a bare "chicken" or "turkey" is called Animal, and NCBI's label is often food.
    assert result["precision"]["Animal"] >= 0.73
    assert result["recall"]["Food"] >= 0.70


def test_development_baseline(records, lexicon):
    result = _score([r for r in records if r["split"] == "development"], lexicon)
    assert result["coverage"] >= 0.935 and result["accuracy"] >= 0.945


def test_v5_fixes_what_v4_got_wrong_on_records_it_was_not_tuned_on(records, lexicon):
    """v4 called retail meat Animal: fewer than half of its Animal calls agreed with NCBI."""
    heldout = [r for r in records if r["split"] == "heldout"]
    v4, v5 = _score(heldout, load_lexicon(Path("v4"))), _score(heldout, lexicon)
    assert v4["precision"]["Animal"] < 0.5 and v4["recall"]["Food"] < 0.35
    assert (round(v4["coverage"], 3), round(v4["accuracy"], 3)) == (0.769, 0.809)
    assert v5["coverage"] > v4["coverage"] + 0.1 and v5["accuracy"] > v4["accuracy"] + 0.1
    for label in LABEL.values():
        assert v5["recall"][label] >= v4["recall"][label] - 0.01, label
        assert v5["precision"][label] >= v4["precision"][label] - 0.01, label
