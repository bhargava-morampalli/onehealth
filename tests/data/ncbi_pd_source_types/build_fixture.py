"""Rebuild the NCBI source_type benchmark (standard library only).

    python tests/data/ncbi_pd_source_types/build_fixture.py

NCBI Pathogen Detection publishes one metadata table per organism group. Some
submitters fill in `source_type`, NCBI's own One Health category (human, animal, food,
environmental). Those records give a label that onehealth did not make and that the
classifier can be tested against: hide `source_type`, classify from `isolation_source`
and `host` alone, and compare.

The tables are large (Salmonella is 900 MB), so this reads evenly spaced parts of each
with HTTP range requests instead of whole files. Nothing is simulated.

* development   Eight 3 MB parts of each of six organism groups. Lexicon v5 was written
                while looking at these records.
* heldout       Eight other parts of the same six tables, half-way between the first
                eight, and four organism groups not in the development split. Nobody
                looked at these before v5 was finished.

Records are counted, not listed: one output row per distinct (organism group,
isolation_source, host, host_disease, source_type), with the number of records. A
record with neither an isolation_source nor a host is left out, because no text
classifier can place it.

Each release is fixed below, so the output is reproducible for as long as NCBI keeps
those releases online.
"""

from __future__ import annotations

import sys
import urllib.request
from collections import Counter
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results"
OUT = Path(__file__).resolve().parent / "records.tsv"
PARTS = 8
PART_BYTES = 3_000_000
LABELS = ("human", "animal", "food", "environmental")
COLUMNS = ("isolation_source", "host", "host_disease")

# Organism group -> release (all read on 14 November 2026).
DEVELOPMENT = {
    "Salmonella": "PDG000000002.4298",
    "Escherichia_coli_Shigella": "PDG000000004.6440",
    "Listeria": "PDG000000001.4986",
    "Campylobacter": "PDG000000003.2951",
    "Klebsiella": "PDG000000012.2579",
    "Acinetobacter": "PDG000000010.2002",
}
HELDOUT_ONLY = {
    "Staphylococcus_aureus": "PDG000000073.1331",
    "Vibrio_parahaemolyticus": "PDG000000023.1373",
    "Bacillus_cereus_group": "PDG000000088.787",
    "Cronobacter": "PDG000000043.488",
}


def url(group: str, release: str) -> str:
    return f"{BASE}/{group}/{release}/Metadata/{release}.metadata.tsv"


def size(address: str) -> int:
    request = urllib.request.Request(address, method="HEAD")
    with urllib.request.urlopen(request, timeout=120) as response:
        return int(response.headers["Content-Length"])


def read_range(address: str, start: int, end: int) -> str:
    """Bytes start..end inclusive, as text. A multi-byte character cut at either end is replaced."""
    request = urllib.request.Request(address, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(request, timeout=600) as response:
        return response.read().decode("utf-8", errors="replace")


def header(address: str) -> list[str]:
    return read_range(address, 0, 9_999).split("\n", 1)[0].lstrip("#").split("\t")


def offsets(total: int, split: str) -> list[int]:
    """Where each part starts. The held-out parts start half-way between the development parts."""
    step = total // PARTS
    shift = 5_000 if split == "development" else step // 2
    return [step * k + shift for k in range(PARTS)]


def records(group: str, release: str, split: str) -> Counter:
    """Distinct labelled records in the parts of one table, counted."""
    address = url(group, release)
    names = header(address)
    column = {name: i for i, name in enumerate(names)}
    total = size(address)
    if total <= PARTS * PART_BYTES:
        parts = [read_range(address, 0, total - 1).split("\n")[1:]]  # a small table: all of it, minus the header
    else:
        # The first and last line of a part are cut short, so they are dropped.
        parts = [read_range(address, start, start + PART_BYTES).split("\n")[1:-1] for start in offsets(total, split)]
    seen: set[str] = set()
    counts: Counter = Counter()
    for lines in parts:
        for line in lines:
            cells = line.rstrip("\r").split("\t")
            if len(cells) != len(names):
                continue

            def cell(name: str) -> str:
                value = cells[column[name]].strip()
                return "" if value == "NULL" else value

            target, label = cell("target_acc"), cell("source_type").lower()
            if label not in LABELS or target in seen:
                continue
            seen.add(target)
            if cell("isolation_source") or cell("host"):
                counts[(group, *(cell(name) for name in COLUMNS), label)] += 1
    return counts


def main() -> int:
    plan = [("development", group, release) for group, release in DEVELOPMENT.items()]
    plan += [("heldout", group, release) for group, release in {**DEVELOPMENT, **HELDOUT_ONLY}.items()]
    lines = ["# NCBI Pathogen Detection records that carry NCBI's own source_type. See README.md."]
    lines += [f"# {group} {release}" for group, release in {**DEVELOPMENT, **HELDOUT_ONLY}.items()]
    lines.append("\t".join(("split", "organism_group", *COLUMNS, "ncbi_source_type", "records")))
    for split, group, release in plan:
        counts = records(group, release, split)
        print(f"{split:12} {group:28} {release}: {sum(counts.values()):6} records, {len(counts):5} distinct", file=sys.stderr)
        for key, n in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append("\t".join((split, *key, str(n))))
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(lines) - len(DEVELOPMENT) - len(HELDOUT_ONLY) - 2} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
