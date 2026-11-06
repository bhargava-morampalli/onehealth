"""Compare onehealth's source classes with independent, curated labels.

Use this to measure the lexicon against a hand-curated set before trusting it on a
new dataset. Three numbers summarise the result:

* coverage  - of isolates with a real class in the truth set, the share onehealth
              classified (the rest went to the review queue);
* accuracy  - of those it classified, the share whose class matches the truth;
* agreement - the share of all isolates where onehealth and the truth agree,
              counting "Unclassified" on both sides as agreement.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from onehealth.constants import ID_COLUMN_ALIASES, L1_VALUES, PAPER_ROLLUP_VALUES
from onehealth.errors import InputError
from onehealth.io_tsv import (
    ISOLATE_NAME_ALIASES,
    first_present,
    isolate_key,
    read_tsv,
    row_match_ids,
)

LEVELS = {"l1": L1_VALUES, "paper_rollup": PAPER_ROLLUP_VALUES}
UNCLASSIFIED = "Unclassified"
_TRUTH_ID_COLUMNS = (*ID_COLUMN_ALIASES, "biosample", "biosample_acc", "assembly", "target_acc", *ISOLATE_NAME_ALIASES)


@dataclass
class Validation:
    level: str
    isolates: int  # isolates that have a truth label and appear in the classified table
    unmatched: list[str]  # truth IDs not found in the classified table
    confusion: Counter = field(default_factory=Counter)  # (truth, predicted) -> isolates
    disagreements: list[dict] = field(default_factory=list)

    @property
    def classes(self) -> list[str]:
        order = list(LEVELS[self.level])
        seen = {label for pair in self.confusion for label in pair}
        return [c for c in order if c in seen]

    def count(self, truth: str | None = None, predicted: str | None = None) -> int:
        """Isolates with this truth class and/or predicted class (None matches any)."""
        return sum(
            n
            for (t, p), n in self.confusion.items()
            if (truth is None or t == truth) and (predicted is None or p == predicted)
        )

    def _pairs(self) -> list[tuple[str, str]]:
        return list(self.confusion.elements())  # one (truth, predicted) pair per isolate

    @property
    def coverage(self) -> float:
        labelled = [(t, p) for t, p in self._pairs() if t != UNCLASSIFIED]
        classified = [(t, p) for t, p in labelled if p != UNCLASSIFIED]
        return len(classified) / len(labelled) if labelled else float("nan")

    @property
    def accuracy(self) -> float:
        classified = [(t, p) for t, p in self._pairs() if t != UNCLASSIFIED and p != UNCLASSIFIED]
        return sum(t == p for t, p in classified) / len(classified) if classified else float("nan")

    @property
    def agreement(self) -> float:
        agree = sum(n for (t, p), n in self.confusion.items() if t == p)
        return agree / self.isolates if self.isolates else float("nan")

    def per_class(self) -> list[dict]:
        """Recall: of isolates truly in a class, the share onehealth put there.
        Precision: of isolates onehealth put in a class, the share truly in it."""
        rows = []
        for label in self.classes:
            truth_n, predicted_n = self.count(truth=label), self.count(predicted=label)
            correct = self.count(label, label)
            rows.append(
                {
                    "class": label,
                    "truth": truth_n,
                    "predicted": predicted_n,
                    "correct": correct,
                    "recall": correct / truth_n if truth_n else float("nan"),
                    "precision": correct / predicted_n if predicted_n else float("nan"),
                }
            )
        return rows


def read_truth(path: Path, column: str, level: str) -> dict[str, str]:
    """Isolate ID -> truth class. Labels must name a class of `level` (case-insensitive)."""
    rows = read_tsv(path)
    if not rows:
        return {}
    id_column = first_present(rows[0], _TRUTH_ID_COLUMNS)
    if id_column is None or column not in rows[0]:
        raise InputError(f"{path}: needs an ID column ({', '.join(_TRUTH_ID_COLUMNS[:5])}, ...) and a {column!r} column")
    classes = {name.lower(): name for name in LEVELS[level]}
    classes["environment"] = "Environmental"
    truth: dict[str, str] = {}
    unknown, conflicts = set(), set()
    for row in rows:
        label = row[column].strip()
        if not row[id_column] or not label:
            continue
        if label.lower() not in classes:
            unknown.add(label)
            continue
        key, value = row[id_column].strip().upper(), classes[label.lower()]
        if truth.get(key, value) != value:
            conflicts.add(key)
        truth[key] = value
    if unknown:
        raise InputError(
            f"{path}: truth labels {sorted(unknown)} are not {level} classes ({', '.join(LEVELS[level])}); "
            "map them before validating"
        )
    if conflicts:
        # Keeping either label would bias accuracy without anyone noticing.
        raise InputError(f"{path}: {len(conflicts)} isolates have two different truth labels "
                         f"(e.g. {', '.join(sorted(conflicts)[:3])}); keep one label per isolate")
    return truth


def validate(classified: list[dict[str, str]], truth: dict[str, str], level: str = "paper_rollup") -> Validation:
    """Cross-tabulate truth against onehealth's class for every isolate in both tables."""
    if level not in LEVELS:
        raise InputError(f"level must be one of {list(LEVELS)}")
    by_id: dict[str, dict[str, str]] = {}
    for row in classified:
        for identifier in row_match_ids(row):
            by_id.setdefault(identifier, row)

    result = Validation(level=level, isolates=0, unmatched=[])
    seen: dict[str, tuple[str, str]] = {}  # isolate key -> first truth ID and label for it
    conflicts: list[str] = []
    for identifier, expected in truth.items():
        row = by_id.get(identifier)
        if row is None:
            result.unmatched.append(identifier)
            continue
        key = isolate_key(row) or str(id(row))
        if key in seen:
            # Two truth IDs for one isolate (say BioSample and assembly): count it once, if they agree.
            first_id, first_label = seen[key]
            if first_label != expected:
                conflicts.append(f"{first_id} and {identifier}")
            continue
        seen[key] = (identifier, expected)
        predicted = (row.get(level) or UNCLASSIFIED).strip() or UNCLASSIFIED
        result.isolates += 1
        result.confusion[(expected, predicted)] += 1
        if predicted != expected:
            result.disagreements.append(
                {
                    "id": identifier,
                    "truth": expected,
                    "predicted": predicted,
                    "host": row.get("host", ""),
                    "isolation_source": row.get("isolation_source", ""),
                    "sample_type": row.get("sample_type", ""),
                    "classify_reason": row.get("classify_reason", ""),
                }
            )
    if conflicts:
        raise InputError(f"the truth labels disagree for {len(conflicts)} isolates listed under two IDs "
                         f"(e.g. {'; '.join(conflicts[:3])}); keep one label per isolate")
    return result
