"""TSV helpers and column aliasing."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from onehealth.constants import (
    ASSEMBLY_ALIASES,
    BIOSAMPLE_ALIASES,
    CLUSTER_ALIASES,
    HOST_ALIASES,
    HOST_DISEASE_ALIASES,
    ID_COLUMN_ALIASES,
    ISOLATION_SOURCE_ALIASES,
    PDT_ALIASES,
    PLASMID_ALIASES,
    SAMPLE_TYPE_ALIASES,
    SRA_ALIASES,
)
from onehealth.ids import parse_id


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if not sample.strip():
            return []
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        rows = []
        for raw in reader:
            cleaned = {(key or "").strip(): (value or "").strip() for key, value in raw.items()}
            first = next(iter(cleaned.values()), "")
            if first.startswith("#"):
                continue
            rows.append(cleaned)
        return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in names})


def read_id_list(path: Path) -> list[str]:
    """Accept one ID per line, or a TSV with an id-like column."""
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        return []
    first = next((line for line in text.splitlines() if line.strip() and not line.startswith("#")), "")
    if "\t" in first:
        rows = read_tsv(path)
        if rows:
            key = first_present(rows[0], ID_COLUMN_ALIASES)
            if key:
                return [row.get(key, "").strip() for row in rows if row.get(key, "").strip()]
    ids: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        ids.append(stripped.split("\t")[0].strip())
    return ids


def first_present(row: dict[str, str], aliases: Iterable[str]) -> str | None:
    lower = {key.lower(): key for key in row}
    for alias in aliases:
        if alias.lower() in lower:
            return lower[alias.lower()]
    return None


def get_field(row: dict[str, str], aliases: Iterable[str]) -> str:
    """Return the first non-empty alias. A blank isolation_source must not hide source."""
    lower = {key.lower(): key for key in row}
    for alias in aliases:
        key = lower.get(alias.lower())
        if key is None:
            continue
        value = (row.get(key) or "").strip()
        if value:
            return value
    return ""


def row_match_ids(row: dict[str, str]) -> set[str]:
    """All identifiers on a row that an override file may key on."""
    values = [
        get_field(row, ID_COLUMN_ALIASES),
        get_field(row, BIOSAMPLE_ALIASES),
        get_field(row, ASSEMBLY_ALIASES),
        get_field(row, SRA_ALIASES),
        get_field(row, PDT_ALIASES),
    ]
    ids = set()
    for value in values:
        parsed = parse_id(value)
        if parsed.normalized:
            ids.add(parsed.normalized.upper())
            ids.add(value.strip().upper())
    return {item for item in ids if item}


def find_column(rows: list[dict[str, str]], aliases: Iterable[str]) -> str | None:
    if not rows:
        return None
    return first_present(rows[0], aliases)


def optional_group_column(rows: list[dict[str, str]], aliases: Iterable[str]) -> str | None:
    return find_column(rows, aliases)


def metadata_view(row: dict[str, str]) -> dict[str, str]:
    return {
        "host": get_field(row, HOST_ALIASES),
        "isolation_source": get_field(row, ISOLATION_SOURCE_ALIASES),
        "host_disease": get_field(row, HOST_DISEASE_ALIASES),
        "sample_type": get_field(row, SAMPLE_TYPE_ALIASES),
        "biosample": get_field(row, BIOSAMPLE_ALIASES),
        "assembly": get_field(row, ASSEMBLY_ALIASES),
        "sra": get_field(row, SRA_ALIASES),
        "pdt": get_field(row, PDT_ALIASES),
        "cluster_id": get_field(row, CLUSTER_ALIASES),
        "plasmid_id": get_field(row, PLASMID_ALIASES),
    }
