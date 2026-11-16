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
                looked at these before v5 was finished. v7 was checked against them.
* validation    Eight more parts of the same six tables, right after the held-out
                parts, and three organism groups in neither other split. Meant to be
                read once, after v7; it showed that one v7 rule (a farm with a host)
                went against NCBI's records, and the rule was corrected.
* final         Eight more parts of the same six tables, right after the
                development parts, and five more organism groups. Read once, after
                v7 was corrected, and not used to change v7. Lexicon v8 was written
                with all four splits above in view.
* fresh         Drawn after v8 was finished and read once. Eight parts of the four
                tables that still had unread parts (right after the validation parts;
                Listeria and Acinetobacter were already read in full), and eight more
                organism groups, chosen by name before any of their records were seen.

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
# Read on 15 November 2026.
VALIDATION_ONLY = {
    "Clostridium_perfringens": "PDG000000061.381",
    "Vibrio_cholerae": "PDG000000055.1027",
    "Yersinia_enterocolitica": "PDG000000087.760",
}
FINAL_ONLY = {
    "Vibrio_vulnificus": "PDG000000058.490",
    "Listeria_innocua": "PDG000000091.235",
    "Staphylococcus_pseudintermedius": "PDG000000042.545",
    "Clostridium_botulinum": "PDG000000059.243",
    "Aeromonas": "PDG000000068.313",
}
# Read on 16 November 2026.
FRESH_ONLY = {
    "Vibrio_alginolyticus": "PDG000000136.164",
    "Enterococcus_faecalis": "PDG000000072.881",
    "Enterococcus_faecium": "PDG000000071.979",
    "Citrobacter_freundii": "PDG000000039.1355",
    "Aeromonas_veronii": "PDG000000076.214",
    "Pseudomonas_aeruginosa": "PDG000000036.1837",
    "Enterobacter_hormaechei": "PDG000000105.895",
    "Serratia": "PDG000000016.1070",
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


def offsets(total: int, split: str, shared: bool) -> list[tuple[int, int]]:
    """(start, length) of each part.

    The development parts start each eighth of the table, and the held-out parts half-way
    between them. In the six tables every split reads, the final parts come right after
    the development parts and the validation parts right after the held-out parts, each
    shortened where the table is too small for a full part. The fresh parts come right
    after the validation parts, where the table has room left; a table that has none
    gives no parts. The byte ranges do not overlap, so no record is in two splits. A
    table that only one of the later splits reads is read at the development positions.
    """
    step = total // PARTS
    if split == "development" or (split in ("validation", "final", "fresh") and not shared):
        shift, length = 5_000, PART_BYTES
    elif split == "heldout":
        shift, length = step // 2, PART_BYTES
    elif split == "validation":
        shift = step // 2 + PART_BYTES + 5_000
        length = min(PART_BYTES, step - shift - 5_000)
    elif split == "final":
        shift = PART_BYTES + 10_000
        length = min(PART_BYTES, step // 2 - shift - 5_000)
    else:
        shift = step // 2 + 2 * PART_BYTES + 10_000
        length = min(PART_BYTES, step - shift - 5_000)
    return [(step * k + shift, length) for k in range(PARTS)] if length > 0 else []


def _cell(cells: list[str], column: dict[str, int], name: str) -> str:
    """One value of a metadata line; NCBI writes NULL for empty."""
    value = cells[column[name]].strip()
    return "" if value == "NULL" else value


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
        parts = [read_range(address, start, start + length).split("\n")[1:-1]
                 for start, length in offsets(total, split, group in DEVELOPMENT)]
    seen: set[str] = set()
    counts: Counter = Counter()
    for lines in parts:
        for line in lines:
            cells = line.rstrip("\r").split("\t")
            if len(cells) != len(names):
                continue

            target, label = _cell(cells, column, "target_acc"), _cell(cells, column, "source_type").lower()
            if label not in LABELS or target in seen:
                continue
            seen.add(target)
            if _cell(cells, column, "isolation_source") or _cell(cells, column, "host"):
                counts[(group, *(_cell(cells, column, name) for name in COLUMNS), label)] += 1
    return counts


def main() -> int:
    plan = [("development", group, release) for group, release in DEVELOPMENT.items()]
    plan += [("heldout", group, release) for group, release in {**DEVELOPMENT, **HELDOUT_ONLY}.items()]
    plan += [("validation", group, release) for group, release in {**DEVELOPMENT, **VALIDATION_ONLY}.items()]
    plan += [("final", group, release) for group, release in {**DEVELOPMENT, **FINAL_ONLY}.items()]
    plan += [("fresh", group, release) for group, release in {**DEVELOPMENT, **FRESH_ONLY}.items()]
    every = {**DEVELOPMENT, **HELDOUT_ONLY, **VALIDATION_ONLY, **FINAL_ONLY, **FRESH_ONLY}
    lines = ["# NCBI Pathogen Detection records that carry NCBI's own source_type. See README.md."]
    lines += [f"# {group} {release}" for group, release in every.items()]
    lines.append("\t".join(("split", "organism_group", *COLUMNS, "ncbi_source_type", "records")))
    for split, group, release in plan:
        counts = records(group, release, split)
        print(f"{split:12} {group:28} {release}: {sum(counts.values()):6} records, {len(counts):5} distinct", file=sys.stderr)
        for key, n in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append("\t".join((split, *key, str(n))))
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(lines) - len(every) - 2} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
