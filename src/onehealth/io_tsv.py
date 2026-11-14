"""TSV helpers and column aliasing."""

from __future__ import annotations

import csv
import io
import re
from array import array
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Iterator

from onehealth.constants import (
    ASSEMBLY_ALIASES,
    BIOSAMPLE_ALIASES,
    CLUSTER_ALIASES,
    HOST_ALIASES,
    HOST_DISEASE_ALIASES,
    ID_COLUMN_ALIASES,
    ISOLATION_SOURCE_ALIASES,
    PACKAGE_ALIASES,
    PDT_ALIASES,
    PLASMID_ALIASES,
    SAMPLE_TYPE_ALIASES,
    SRA_ALIASES,
)
from onehealth.errors import InputError
from onehealth.ids import assembly_join_keys, parse_id

if TYPE_CHECKING:  # numpy is imported lazily, only by the readers that need it
    import numpy as np

# Columns that can hold an isolate's identifier: the ID column of an ID list, and what
# resolve --metadata rows are matched on.
ID_LIST_ALIASES = ID_COLUMN_ALIASES + BIOSAMPLE_ALIASES + ASSEMBLY_ALIASES + SRA_ALIASES + PDT_ALIASES


def read_text(path: Path) -> str:
    """A text file as UTF-8, without the byte-order mark Excel adds. Any other encoding is reported, not guessed."""
    try:
        return path.read_bytes().decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputError(f"{path}: not UTF-8 text (byte {exc.start}); save it as UTF-8") from None


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV with a header row. Lines starting with '#' are comments, even above the header.

    The one exception is a header that itself starts with '#', as in the tables NCBI
    publishes ("#label" in Pathogen Detection metadata, "#assembly_accession" in
    assembly_summary.txt): see _hash_header.

    One line is one row. A cell wrapped in double quotes (how Excel and pandas write a
    cell that contains a quote) is unwrapped; any other quote is kept as typed. csv's
    quoting rules let a stray quote, as in `"wound swab`, swallow the lines after it.
    """
    # Split lines exactly as a file opened with newline="" does (\n, \r\n or \r only).
    numbered = [
        (number, line)
        for number, line in enumerate(io.StringIO(read_text(path), newline=""), start=1)
        if line.strip()
    ]
    lines = [(number, line) for number, line in numbered if not line.startswith("#")]
    if not lines:
        return []
    hash_header = _hash_header(numbered, lines[0][0])
    if hash_header is not None:
        lines.insert(0, hash_header)
    _reject_other_delimiters(path, lines[0][1])
    header = [cell.strip() for cell in _cells(lines[0][1])]
    _reject_repeated_columns(path, header)
    rows = []
    for number, line in lines[1:]:
        cells = _cells(line)
        if any(cell.strip() for cell in cells[len(header):]):
            raise InputError(f"{path}: line {number} has more cells than the header")
        cells += [""] * (len(header) - len(cells))
        rows.append({key: value.strip() for key, value in zip(header, cells)})
    return rows


# Column names that mark a '#' line as a table's header and not a comment.
_KNOWN_COLUMNS = frozenset(
    name.lower()
    for name in (
        *ID_COLUMN_ALIASES, *BIOSAMPLE_ALIASES, *ASSEMBLY_ALIASES, *SRA_ALIASES, *PDT_ALIASES,
        *HOST_ALIASES, *ISOLATION_SOURCE_ALIASES, *HOST_DISEASE_ALIASES, *SAMPLE_TYPE_ALIASES,
        "collection_date", "geo_loc_name", "organism", "organism_name",
    )
)


def _hash_header(numbered: list[tuple[int, str]], first_row: int) -> tuple[int, str] | None:
    """The header of a table whose header line starts with '#', or None.

    Only the last '#' line before the first row is considered, and only when it is
    tab-separated and names more columns onehealth reads ("host", "biosample_acc", ...)
    than the line below it does. A commented-out data row has values in its cells, not
    column names, and a commented-out copy of the header names no more columns than
    the real one below it, so both stay comments, as before.
    """
    above = [(number, line) for number, line in numbered if number < first_row]
    if not above or "\t" not in above[-1][1]:
        return None
    number, line = above[-1]
    stripped = line.lstrip("#").lstrip(" ")
    below = next(text for row, text in numbered if row == first_row)
    if _known_columns(stripped) <= _known_columns(below):
        return None
    return number, stripped


def _known_columns(line: str) -> int:
    return sum(1 for cell in _cells(line) if cell.strip().lower() in _KNOWN_COLUMNS)


def _reject_repeated_columns(path: Path, header: list[str]) -> None:
    """Two columns with one name cannot both be kept: a row is read as name -> value, so the
    second would silently replace the first (a second `host` column changes the class)."""
    seen: dict[str, int] = {}
    for position, name in enumerate(header, start=1):
        if not name:
            continue
        if name in seen:
            raise InputError(f"{path}: column {name!r} appears twice (columns {seen[name]} and {position}); "
                             "rename or remove one")
        seen[name] = position


_OTHER_DELIMITERS = ((",", "comma"), (";", "semicolon"))


def _reject_other_delimiters(path: Path, header: str) -> None:
    """A comma- or semicolon-separated file would be read as one column, and every row of it
    would then go to review for lack of a host or source. Say so instead.

    Only a header without a tab counts: in a TSV, a comma in a column name is just text.
    """
    if "\t" in header:
        return
    for delimiter, name in _OTHER_DELIMITERS:
        try:
            names = next(csv.reader([header.rstrip("\r\n")], delimiter=delimiter))
        except csv.Error:  # e.g. a NUL byte: not a table this check can say anything about
            return
        if sum(1 for cell in names if cell.strip()) >= 2:
            raise InputError(f"{path}: this looks {name}-separated (the header has {name}s but no tabs); "
                             "save it as tab-separated text")


def _cells(line: str) -> list[str]:
    return [_unquote(cell) for cell in line.rstrip("\r\n").split("\t")]


def _unquote(cell: str) -> str:
    if len(cell) >= 2 and cell[0] == '"' == cell[-1]:
        return cell[1:-1].replace('""', '"')
    return cell


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: _one_line(row.get(name, "")) for name in names})


def _one_line(value: object) -> object:
    """A TSV cell cannot hold a tab or a line break, so they become spaces.

    A cell that a spreadsheet would run as a formula gets a leading apostrophe, which
    Excel and LibreOffice show as text. The values come from public BioSample records,
    which anyone can submit.
    """
    if not isinstance(value, str):
        return value
    if "\t" in value or "\n" in value or "\r" in value:
        value = value.replace("\r\n", " ").replace("\t", " ").replace("\n", " ").replace("\r", " ")
    return f"'{value}" if _is_formula(value) else value


def _is_formula(value: str) -> bool:
    """'=' starts a formula in every spreadsheet. '+', '-' and '@' do too, but they also start
    ordinary values ("-", "-5", "+ve", "@home"), so those count only with a call or a DDE pipe."""
    if value.startswith("="):
        return True
    return value[:1] in ("+", "-", "@") and ("(" in value or "|" in value)


def require_columns(rows: list[dict[str, str]], path: Path, columns: Iterable[str], hint: str = "") -> None:
    """A table with rows must have these columns; a header-only table passes (it has nothing to misread)."""
    if not rows:
        return
    missing = [name for name in columns if name not in rows[0]]
    if missing:
        names = ", ".join(repr(name) for name in missing)
        found = ", ".join(list(rows[0])[:8]) + (", ..." if len(rows[0]) > 8 else "")
        raise InputError(f"{path}: no {names} column{'s' if len(missing) > 1 else ''} (found: {found}). {hint}".strip())


def read_id_list(path: Path) -> list[str]:
    """One ID per line, or a TSV whose header has an id-like column (else its first column).

    No accession contains a comma or semicolon, so a line with one is a CSV table or several
    IDs on one line, and is reported rather than looked up as one unknown ID.
    """
    numbered = [(number, line.strip()) for number, line in enumerate(read_text(path).splitlines(), start=1)]
    numbered = [(number, line) for number, line in numbered if line and not line.startswith("#")]
    if not numbered:
        return []
    if "\t" in numbered[0][1]:
        rows = read_tsv(path)
        if not rows:
            return []
        key = first_present(rows[0], ID_LIST_ALIASES) or next(iter(rows[0]))
        return [row[key] for row in rows if row.get(key)]
    ids = []
    for number, line in numbered:
        value = line.split("\t")[0].strip()
        for delimiter, name in _OTHER_DELIMITERS:
            if delimiter in value:
                raise InputError(f"{path}: line {number} has a {name}; put one ID per line, "
                                 "or save a table of IDs as tab-separated text")
        ids.append(value)
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
    """All identifiers on a row that a local-metadata or override row may key on.

    Assemblies also yield the versionless stem and the GCA/GCF twin, which NCBI
    gives the same numeric part for the same assembly. Rows without an accession or
    ID column fall back to their isolate key (strain, isolate, sample or name).
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
    matched = {item for item in ids if item}
    if not matched and (key := isolate_key(row)):
        return {key}
    return matched


# Columns that name an isolate in tables without accessions: Kleborate `strain`,
# AMRFinderPlus `Name`, lab sheets `sample`.
ISOLATE_NAME_ALIASES = ("strain", "isolate", "sample", "sample_name", "name")


def isolate_key(row: dict[str, str]) -> str:
    """One key per biological isolate, or "" when the row names none.

    The BioSample when it is a single accession; otherwise the first identifier
    found: the query/ID column, assembly, SRA, PDT, then a strain, isolate or
    sample name.
    """
    biosample = get_field(row, BIOSAMPLE_ALIASES)
    if biosample and "," not in biosample:
        return biosample.upper()
    for aliases in (ID_COLUMN_ALIASES, ASSEMBLY_ALIASES, SRA_ALIASES, PDT_ALIASES, ISOLATE_NAME_ALIASES):
        value = get_field(row, aliases)
        if value:
            return value.upper()
    return ""


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
        "package": get_field(row, PACKAGE_ALIASES),
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
        raise InputError(
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


def _is_float(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


def read_distance_matrix(path: Path) -> tuple[list[str], "np.ndarray"]:
    """Pairwise distances from snp-dists: square (default) or long (-m) format, TSV or CSV.

    The file is streamed and the numbers are parsed by numpy, so memory stays close to
    the matrix itself (8 bytes per cell) rather than a Python object per cell.
    """
    import numpy as np

    header = next(_data_lines(path), None)
    if header is None:
        return [], np.zeros((0, 0))
    delimiter = "\t" if "\t" in header else ","
    names = [cell.strip() for cell in header.split(delimiter)[1:]]
    first = _first_cells(path, delimiter) if names else []
    if names and first == names:
        return names, _read_square(path, delimiter, len(names))
    # Every row is named in the header, but the rows are not the columns in order:
    # a square matrix with rows missing, or rows sorted differently. Reading it
    # would put distances on the wrong pairs. (A long-format file's first line ends
    # in a distance, so its "header" is never a list of names; see _is_float.)
    looks_long = len(names) == 2 and _is_float(names[1])
    if first and set(first) <= set(names) and not looks_long:
        if len(first) != len(names):
            raise InputError(f"{path}: distance matrix is not square ({len(first)} rows for {len(names)} columns)")
        raise InputError(f"{path}: the rows of the distance matrix are not in the same order as its columns")
    return _read_long(path, delimiter)


def _data_lines(path: Path) -> Iterator[str]:
    """Lines without their newline, skipping blank lines and lines starting with '#'."""
    with path.open(encoding="utf-8-sig") as handle:
        try:
            for line in handle:
                line = line.rstrip("\r\n")
                if line.strip() and not line.startswith("#"):
                    yield line
        except UnicodeDecodeError as exc:
            raise InputError(f"{path}: not UTF-8 text ({exc.reason}); save it as UTF-8") from None


def _first_cells(path: Path, delimiter: str) -> list[str]:
    """First cell of every line after the header: a square matrix's row names."""
    lines = _data_lines(path)
    next(lines)
    return [line.split(delimiter, 1)[0].strip() for line in lines]


def _read_square(path: Path, delimiter: str, n: int) -> "np.ndarray":
    import numpy as np

    def values() -> Iterator[str]:
        lines = _data_lines(path)
        next(lines)
        for line in lines:
            if line.count(delimiter) != n:
                raise InputError(f"{path}: distance matrix is not square")
            yield line.split(delimiter, 1)[1]

    try:
        matrix = np.loadtxt(values(), delimiter=delimiter, comments=None, dtype=float, ndmin=2)
    except InputError:
        raise
    except ValueError as exc:
        raise InputError(f"{path}: non-numeric distance ({exc})") from exc
    if matrix.shape != (n, n):
        raise InputError(f"{path}: distance matrix is not square")
    _check_square(path, matrix)
    return matrix


def _check_square(path: Path, matrix: "np.ndarray", block: int = 512) -> None:
    """A distance cannot be negative, and A to B is the same as B to A.

    Only the upper triangle is used later, so a file holding just the lower triangle (zeros
    above the diagonal) would read as every pair being 0 apart. Checked a block of rows at a
    time, so nothing the size of the matrix is copied.
    """
    import numpy as np

    n = len(matrix)
    for start in range(0, n, block):
        rows = matrix[start : start + block]
        if (rows < 0).any():
            raise InputError(f"{path}: the distance matrix has a negative distance (row {start + int(np.argmax((rows < 0).any(axis=1))) + 1})")
        mirror = matrix[:, start : start + block].T
        differ = ~np.isclose(rows, mirror, rtol=1e-9, atol=1e-9, equal_nan=True)
        if differ.any():
            i, j = (int(k) for k in np.argwhere(differ)[0])
            raise InputError(
                f"{path}: the distance matrix is not symmetric (row {start + i + 1}, column {j + 1}: "
                f"{rows[i, j]:g}, but {mirror[i, j]:g} the other way round); onehealth needs the full square "
                "matrix, as snp-dists writes it, not one triangle"
            )


def _read_long(path: Path, delimiter: str) -> tuple[list[str], "np.ndarray"]:
    """Three columns: name, name, distance. Lines with a non-numeric distance (a header) are skipped."""
    import numpy as np

    order: dict[str, int] = {}
    first, second, distances = array("q"), array("q"), array("d")
    for line in _data_lines(path):
        cells = line.split(delimiter)
        if len(cells) != 3:
            raise InputError(f"{path}: expected a square snp-dists matrix or three-column (name, name, distance) rows")
        a, b, value = (cell.strip() for cell in cells)
        try:
            distance = float(value)
        except ValueError:
            continue  # a header line in long format
        if distance < 0:
            raise InputError(f"{path}: negative distance between {a} and {b}")
        first.append(order.setdefault(a, len(order)))
        second.append(order.setdefault(b, len(order)))
        distances.append(distance)
    n = len(order)
    matrix = np.full((n, n), np.nan)
    np.fill_diagonal(matrix, 0.0)
    i, j, d = np.array(first, dtype=np.int64), np.array(second, dtype=np.int64), np.array(distances)
    # A pair listed twice keeps its last distance, in both directions.
    pair = np.minimum(i, j) * n + np.maximum(i, j)
    _, last_reversed = np.unique(pair[::-1], return_index=True)
    keep = len(pair) - 1 - last_reversed
    matrix[i[keep], j[keep]] = d[keep]
    matrix[j[keep], i[keep]] = d[keep]
    return list(order), matrix


_LINEAGE_ID_ALIASES = ("id", "strain", "name", "taxon", "sample", "sample_id", "isolate", "genome", "assembly", "biosample")
_LINEAGE_VALUE_ALIASES = ("st", "mlst_st", "sequence_type", "lineage", "sublineage", "sl", "poppunk_cluster", "cluster")


def read_lineages(path: Path) -> dict[str, str]:
    """Isolate -> lineage from a typing table: Kleborate (strain, ST or *__mlst__ST),
    mlst (with or without a header row), or PopPUNK clusters.csv (Taxon, Cluster).
    Tab- or comma-separated."""
    lines = [line for line in read_text(path).splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return {}
    delimiter = "\t" if "\t" in lines[0] else ","
    if _is_mlst_without_header(lines[0], delimiter):
        pairs = _mlst_pairs(path, lines, delimiter)
    else:
        pairs = _table_pairs(path, lines, delimiter)
    out: dict[str, str] = {}
    for name, value in pairs:
        if name and value:
            if out.get(name, value) != value:
                raise InputError(f"{path}: {name} is listed twice, as {out[name]} and {value}; keep one lineage per isolate")
            out[name] = value
    return out


def _table_pairs(path: Path, lines: list[str], delimiter: str) -> list[tuple[str, str]]:
    if delimiter == "\t":
        rows = read_tsv(path)
    else:
        rows = [
            {(k or "").strip(): (v or "").strip() for k, v in row.items() if k}
            for row in csv.DictReader(lines, delimiter=delimiter)
        ]
    if not rows:
        return []
    header = [key for key in rows[0] if key]
    id_col = first_present(rows[0], _LINEAGE_ID_ALIASES) or header[0]
    value_col = first_present(rows[0], _LINEAGE_VALUE_ALIASES) or next(
        (key for key in header if key.lower().endswith("__st")), None
    )
    if value_col is None:
        raise InputError(f"{path}: no lineage column ({', '.join(_LINEAGE_VALUE_ALIASES)}, or a Kleborate *__ST column)")
    # Sequence type 0 means "no ST" (Kleborate writes it when the alleles match none), like
    # mlst's "-". Kept, every untyped isolate would share a lineage called 0.
    sequence_types = value_col.lower() in ("st", "mlst_st", "sequence_type") or value_col.lower().endswith("__st")
    pairs = [((row.get(id_col) or "").strip(), (row.get(value_col) or "").strip()) for row in rows]
    return [(name, "" if sequence_types and value == "0" else value) for name, value in pairs]


# mlst's default output has no header: FILE, SCHEME, ST, then one allele per locus as
# name(number), such as gapA(3), infB(~4) or pgi(-). ST is "-" when none is assigned, and a
# file with no matching scheme has only "FILE - -".
_MLST_ST = re.compile(r"^(\d+|-)$")
_MLST_ALLELE = re.compile(r"^[^\s()]+\(")


def _split(line: str, delimiter: str) -> list[str]:
    cells = line.split("\t") if delimiter == "\t" else next(csv.reader([line], delimiter=delimiter))
    return [cell.strip() for cell in cells]


def _is_mlst_without_header(line: str, delimiter: str) -> bool:
    """A header has names in it, never an ST number or "-" in the third column."""
    cells = _split(line, delimiter)
    return len(cells) >= 3 and bool(_MLST_ST.match(cells[2])) and (len(cells) == 3 or bool(_MLST_ALLELE.match(cells[3])))


def _mlst_pairs(path: Path, lines: list[str], delimiter: str) -> list[tuple[str, str]]:
    pairs = []
    for line in lines:
        cells = _split(line, delimiter)
        if len(cells) < 3 or not _MLST_ST.match(cells[2]):
            raise InputError(f"{path}: this looks like mlst output, but {line[:60]!r} is not an mlst line "
                             "(file, scheme, ST, alleles)")
        pairs.append((cells[0], "" if cells[2] == "-" else cells[2]))
    return pairs
