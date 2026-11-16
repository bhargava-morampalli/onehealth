"""Rebuild the table of animal names from a newer NCBI Taxonomy, and say what that changes.

    python tools/refresh_animal_names.py                          # downloads today's dump (160 MB)
    python tools/refresh_animal_names.py new_taxdump.tar.gz       # or reads a copy you already have
    python tools/refresh_animal_names.py --summary refresh.md     # writes the report to a file

NCBI adds and renames species every day, so the table that ships with onehealth
(src/onehealth/data/animal_names.tsv.gz) slowly falls behind. A newer table can change
classes: a host or a specimen text that names a newly listed animal is placed, where it
went to review before. So a refresh is a change to onehealth's results, and is made the
way any such change is: looked at, then released with a new version.

This tool does the looking. It builds a table from the dump and compares it with the one
in the repository:

1. The names: how many of each kind were added and removed, with examples.
2. The classes: every row of the three benchmark datasets in tests/data is classified with
   the old table and with the new one, and every row whose class changes is listed. For
   the NCBI source_type benchmark it also says whether each change agrees with NCBI's label.

If no name changed, nothing is written: a new dump with the same names is not a change.
Otherwise the new table replaces the old one, and the report says what to do next.

On GitHub the workflow .github/workflows/taxonomy.yml runs this once a month, runs the
tests with the new table, and opens a pull request with the report.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tempfile
from collections import Counter
from collections.abc import Callable, Iterable
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

import build_animal_names as builder  # noqa: E402

from onehealth import taxonomy  # noqa: E402
from onehealth.classify import classify_record  # noqa: E402
from onehealth.io_tsv import read_tsv  # noqa: E402
from onehealth.lexicon import Lexicon, load_lexicon  # noqa: E402

TABLE = builder.OUT
KINDS = {"s": "scientific names", "v": "common names of chordates", "c": "common names of other animals"}
LABEL = {"human": "Human", "animal": "Animal", "food": "Food", "environmental": "Environmental"}
EXAMPLES = 12
LISTED = 25


def read_table(path: Path) -> tuple[dict[str, set[str]], str]:
    """The names of each kind in a table, and the date of the dump it came from ("" if it does not say)."""
    names: dict[str, set[str]] = {kind: set() for kind in KINDS}
    date = ""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# NCBI Taxonomy of "):
                date = line.split()[4]
            elif not line.startswith("#"):
                kind, _, name = line.rstrip("\n").partition("\t")
                names[kind].add(name)
    return names, date


def compare(old: dict[str, set[str]], new: dict[str, set[str]]) -> dict[str, tuple[list[str], list[str]]]:
    """For each kind of name: (added, removed), sorted."""
    return {kind: (sorted(new[kind] - old[kind]), sorted(old[kind] - new[kind])) for kind in KINDS}


# --- the classes ---------------------------------------------------------------------------------------


# (name, rows, how much each row counts, NCBI's label or None)
Dataset = tuple[str, list[dict[str, str]], Callable[[dict[str, str]], int], Callable[[dict[str, str]], str | None]]


def datasets(root: Path = ROOT) -> list[Dataset]:
    data = root / "tests" / "data"
    benchmark = [{**row} for row in read_tsv(data / "ncbi_pd_source_types" / "records.tsv")]
    return [
        # source_type is not in this file (it is ncbi_source_type), so the text alone decides, as in the tests.
        ("NCBI source_type benchmark", benchmark, lambda row: int(row["records"]),
         lambda row: LABEL.get(row["ncbi_source_type"])),
        ("IMP carbapenemase isolates", read_tsv(data / "imp_carbapenemase_sources" / "samples.tsv"), lambda row: 1,
         lambda row: None),
        ("NCBI Klebsiella", read_tsv(data / "ncbi_pd_klebsiella" / "metadata.tsv"), lambda row: 1, lambda row: None),
    ]


def use_table(path: Path) -> None:
    """Classify with this table from now on."""
    taxonomy._names_path = lambda: path  # type: ignore[assignment,return-value]
    taxonomy.forget_tables()


def classes(rows: Iterable[dict[str, str]], lexicon: Lexicon) -> list[str]:
    return [classify_record(row, lexicon).l1 for row in rows]


def class_changes(old: Path, new: Path, lexicon: Lexicon, data: list[Dataset]) -> list[dict]:
    """For each dataset: the rows whose class differs between the two tables."""
    original = taxonomy._names_path
    results = []
    try:
        for name, rows, weight, label in data:
            use_table(old)
            before = classes(rows, lexicon)
            use_table(new)
            after = classes(rows, lexicon)
            changed: Counter = Counter()
            toward = away = 0
            for row, was, now in zip(rows, before, after):
                if was == now:
                    continue
                n = weight(row)
                changed[(row.get("isolation_source", ""), row.get("host", ""), was, now)] += n
                truth = label(row)
                toward += n if truth == now else 0
                away += n if truth == was else 0
            results.append({"dataset": name, "rows": sum(weight(row) for row in rows), "changes": changed,
                            "toward": toward, "away": away, "labelled": label(rows[0]) is not None if rows else False})
    finally:
        taxonomy._names_path = original  # type: ignore[assignment]
        taxonomy.forget_tables()
    return results


# --- the report --------------------------------------------------------------------------------------


def report(old_date: str, new_date: str, names: dict[str, tuple[list[str], list[str]]], changes: list[dict]) -> str:
    lines = [f"## Animal names from NCBI Taxonomy of {new_date} (was {old_date or 'not recorded'})", ""]
    lines += ["### Names", "", "| Kind | Added | Removed |", "|---|---|---|"]
    lines += [f"| {KINDS[kind]} | {len(added):,} | {len(removed):,} |" for kind, (added, removed) in names.items()]
    lines.append("")
    for kind, (added, removed) in names.items():
        for verb, listed in (("added", added), ("removed", removed)):
            if listed:
                more = f", and {len(listed) - EXAMPLES:,} more" if len(listed) > EXAMPLES else ""
                lines.append(f"- {KINDS[kind].capitalize()} {verb}: " + ", ".join(f"`{n}`" for n in listed[:EXAMPLES]) + more)
    lines += ["", "### Classes", "",
              "Every row of the benchmark datasets in `tests/data`, classified with the old table and with the new one.", "",
              "| Dataset | Records | Records whose class changes | Of those, now agreeing with NCBI | No longer agreeing |",
              "|---|---|---|---|---|"]
    for result in changes:
        total = sum(result["changes"].values())
        agree = (f"{result['toward']:,}", f"{result['away']:,}") if result["labelled"] else ("", "")
        lines.append(f"| {result['dataset']} | {result['rows']:,} | {total:,} | {agree[0]} | {agree[1]} |")
    for result in changes:
        if result["changes"]:
            lines += ["", f"**{result['dataset']}**, the texts that change (records):", ""]
            for (source, host, was, now), n in result["changes"].most_common(LISTED):
                lines.append(f"- {n:,}: `{source}`, host `{host}`: {was} to {now}")
            if len(result["changes"]) > LISTED:
                lines.append(f"- and {len(result['changes']) - LISTED:,} more texts")
    lines += ["", "### Before merging", "",
              "1. Read the class changes above. A change that disagrees with NCBI, or a name that is an everyday word "
              "(a moth called `tip`), may need an entry in `_WORDS_NOT_NAMES` in `src/onehealth/taxonomy.py`.",
              "2. Run the tests. Tests that quote exact counts from the benchmark fail when the table changes a class; "
              "update the counts, and the documentation pages that quote them, if the change is right.",
              "3. Add a line to `CHANGELOG.md` saying that the animal names now come from NCBI Taxonomy of this "
              "date. If any class changes above, also raise the package version in `pyproject.toml` and "
              "`src/onehealth/__init__.py`: results made before and after are then not comparable row for row.",
              "4. `classify` prints the date of the table it used, so a methods section can name it either way."]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Rebuild the table of animal names and report what it changes.")
    parser.add_argument("dump", nargs="?", help="new_taxdump.tar.gz (default: download it from NCBI)")
    parser.add_argument("--summary", type=Path, help="Write the report to this file (default: print it)")
    args = parser.parse_args(argv)

    data, source = builder.read_dump(args.dump)
    built = builder.build(data)
    new_date = builder.dump_date(data)
    old_names, old_date = read_table(TABLE)
    names = compare(old_names, dict(zip(KINDS, built)))
    changed = any(added or removed for added, removed in names.values())
    outputs = {"changed": str(changed).lower(), "date": new_date}
    github = os.environ.get("GITHUB_OUTPUT")
    if github:
        with open(github, "a", encoding="utf-8") as handle:
            handle.writelines(f"{key}={value}\n" for key, value in outputs.items())
    if not changed:
        text = f"The names in NCBI Taxonomy of {new_date} are the ones already in the table ({old_date}). Nothing to do.\n"
        print(text, end="", file=sys.stderr)
        if args.summary:
            args.summary.write_text(text, encoding="utf-8")
        return 0

    with tempfile.TemporaryDirectory() as folder:
        new_table = Path(folder) / "animal_names.tsv.gz"
        new_table.write_bytes(builder.render(built, source, hashlib.md5(data).hexdigest(), new_date))
        old_table = Path(folder) / "old_animal_names.tsv.gz"
        old_table.write_bytes(TABLE.read_bytes())
        changes = class_changes(old_table, new_table, load_lexicon(), datasets())
        TABLE.write_bytes(new_table.read_bytes())
    text = report(old_date, new_date, names, changes)
    if args.summary:
        args.summary.write_text(text, encoding="utf-8")
    else:
        print(text)
    print(f"wrote {TABLE}: animal names from NCBI Taxonomy of {new_date}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
