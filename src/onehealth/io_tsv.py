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
from onehealth.ids import assembly_join_keys, parse_id

# Columns that can hold the input identifier of an ID list.
_ID_LIST_ALIASES = ID_COLUMN_ALIASES + BIOSAMPLE_ALIASES + ASSEMBLY_ALIASES + SRA_ALIASES + PDT_ALIASES


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV with a header row. Lines starting with '#' are comments, even above the header."""
    with path.open(newline="", encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    rows = []
    for raw in reader:
        extra = raw.pop(None, None)  # cells beyond the header
        if extra and any((cell or "").strip() for cell in extra):
            raise ValueError(f"{path}: data line {reader.line_num} has more cells than the header")
        rows.append({(key or "").strip(): (value or "").strip() for key, value in raw.items()})
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
    """One ID per line, or a TSV whose header has an id-like column (else its first column)."""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line and not line.startswith("#")]
    if not lines:
        return []
    if "\t" in lines[0]:
        rows = read_tsv(path)
        if not rows:
            return []
        key = first_present(rows[0], _ID_LIST_ALIASES) or next(iter(rows[0]))
        return [row[key] for row in rows if row.get(key)]
    return [line.split("\t")[0].strip() for line in lines]


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
    """All identifiers on a row that a local-metadata or override row may key on.

    Assemblies also yield the versionless stem and the GCA/GCF twin, which NCBI
    gives the same numeric part for the same assembly.
    """
    values = [
        get_field(row, ID_COLUMN_ALIASES),
        get_field(row, BIOSAMPLE_ALIASES),
        get_field(row, ASSEMBLY_ALIASES),
        get_field(row, SRA_ALIASES),
        get_field(row, PDT_ALIASES),
    ]
    ids: set[str] = set()
    for value in values:
        parsed = parse_id(value)
        if parsed.normalized:
            ids.add(parsed.normalized.upper())
            ids.add(value.strip().upper())
            ids.update(assembly_join_keys(parsed.normalized))
    return {item for item in ids if item}


def isolate_key(row: dict[str, str]) -> str:
    """One key per biological isolate: the BioSample when it is a single accession, else the query ID."""
    biosample = get_field(row, BIOSAMPLE_ALIASES)
    if biosample and "," not in biosample:
        return biosample.upper()
    return get_field(row, ID_COLUMN_ALIASES).upper()


def find_column(rows: list[dict[str, str]], aliases: Iterable[str]) -> str | None:
    if not rows:
        return None
    return first_present(rows[0], aliases)


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


# --- genomic inputs for `onehealth report` -------------------------------------

_FEATURE_ID_ALIASES = ("id", "name", "sample", "sample_id", "isolate", "strain", "genome", "assembly", "biosample", "query_id")
_FEATURE_NAME_ALIASES = ("feature", "gene", "gene_symbol", "gene symbol", "element_symbol", "element symbol")
_FEATURE_CLASS_ALIASES = ("class", "feature_class", "group", "drug_class", "category")


def read_features(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Long table: one row per (isolate, feature), optional class column.

    AMRFinderPlus output run with --name fits directly (Name, Element symbol, Class).
    Every isolate named in the file counts as genotyped; a row with an empty
    feature records an isolate with no hits.
    """
    rows = read_tsv(path)
    if not rows:
        return {}, {}
    id_col = first_present(rows[0], _FEATURE_ID_ALIASES)
    feature_col = first_present(rows[0], _FEATURE_NAME_ALIASES)
    if id_col is None or feature_col is None:
        raise ValueError(
            f"{path}: need an isolate column ({', '.join(_FEATURE_ID_ALIASES[:4])}, ...) and a feature column "
            f"({', '.join(_FEATURE_NAME_ALIASES)})"
        )
    class_col = first_present(rows[0], _FEATURE_CLASS_ALIASES)
    presence: dict[str, set[str]] = {}
    classes: dict[str, str] = {}
    for row in rows:
        name = row.get(id_col, "").strip()
        if not name:
            continue
        found = presence.setdefault(name, set())
        feature = row.get(feature_col, "").strip()
        if feature and feature != "-":
            found.add(feature)
            if class_col and row.get(class_col, "").strip():
                classes.setdefault(feature, row[class_col].strip())
    return presence, classes


def read_distance_matrix(path: Path) -> tuple[list[str], "np.ndarray"]:
    """Pairwise distances from snp-dists: square (default) or long (-m) format, TSV or CSV."""
    import numpy as np

    lines = [line.rstrip("\r\n") for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line.strip() and not line.startswith("#")]
    if not lines:
        return [], np.zeros((0, 0))
    delimiter = "\t" if "\t" in lines[0] else ","
    cells = [[cell.strip() for cell in line.split(delimiter)] for line in lines]
    header_names = cells[0][1:]
    if header_names and header_names == [row[0] for row in cells[1:]]:
        names = header_names
        try:
            matrix = np.array([[float(value) for value in row[1:]] for row in cells[1:]])
        except ValueError as exc:
            raise ValueError(f"{path}: non-numeric distance ({exc})") from exc
        if matrix.shape != (len(names), len(names)):
            raise ValueError(f"{path}: distance matrix is not square")
        return names, matrix
    if any(len(row) != 3 for row in cells):
        raise ValueError(f"{path}: expected a square snp-dists matrix or three-column (name, name, distance) rows")
    order: dict[str, int] = {}
    triples = []
    for a, b, value in cells:
        try:
            distance = float(value)
        except ValueError:
            continue  # a header line in long format
        triples.append((order.setdefault(a, len(order)), order.setdefault(b, len(order)), distance))
    matrix = np.full((len(order), len(order)), np.nan)
    np.fill_diagonal(matrix, 0.0)
    for i, j, distance in triples:
        matrix[i, j] = matrix[j, i] = distance
    return list(order), matrix


_LINEAGE_ID_ALIASES = ("id", "strain", "name", "taxon", "sample", "sample_id", "isolate", "genome", "assembly", "biosample")
_LINEAGE_VALUE_ALIASES = ("st", "mlst_st", "sequence_type", "lineage", "sublineage", "sl", "poppunk_cluster", "cluster")


def read_lineages(path: Path) -> dict[str, str]:
    """Isolate -> lineage from a typing table: Kleborate (strain, ST or *__mlst__ST),
    mlst, or PopPUNK clusters.csv (Taxon, Cluster). Tab- or comma-separated."""
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return {}
    delimiter = "\t" if "\t" in lines[0] else ","
    rows = list(csv.DictReader(lines, delimiter=delimiter))
    if not rows:
        return {}
    header = [key for key in rows[0] if key]
    id_col = first_present(rows[0], _LINEAGE_ID_ALIASES) or header[0]
    value_col = first_present(rows[0], _LINEAGE_VALUE_ALIASES) or next(
        (key for key in header if key.lower().endswith("__st")), None
    )
    if value_col is None:
        raise ValueError(f"{path}: no lineage column ({', '.join(_LINEAGE_VALUE_ALIASES)}, or a Kleborate *__ST column)")
    out = {}
    for row in rows:
        name, value = (row.get(id_col) or "").strip(), (row.get(value_col) or "").strip()
        if name and value:
            out[name] = value
    return out
