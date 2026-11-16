"""Standard One Health population analyses on a classified table.

The analysis set follows recent One Health genomic studies (cited in the report's methods):
source composition, lineage diversity and sharing, lineage-source enrichment,
resistance/virulence prevalence by source, and genomic linkage between sources
at a SNP threshold with a permutation null. Every function is pure and returns
plain rows (dicts) so each number in the report is also written to a table.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass, field, replace
from datetime import date

import numpy as np

from onehealth.constants import ASSEMBLY_ALIASES, HEALTHCARE_YES, PDT_ALIASES
from onehealth.errors import InputError
from onehealth.io_tsv import first_present, get_field, isolate_key, row_match_ids
from onehealth.text import is_unknown_value

# Fixed display order. Colours follow the source, never its rank (report.py).
SOURCE_GROUPS = ("Human", "Animal", "Environmental", "Food", "Healthcare environment", "Laboratory")
ROLLUP_GROUPS = ("Human", "Animal", "Environmental")
LINEAGE_ALIASES = (
    "st",
    "mlst_st",
    "sequence_type",
    "mlst",
    "lineage",
    "sublineage",
    "sl",
    "poppunk_cluster",
    "cluster_id",
    "cluster",
)
HOST_ID_ALIASES = ("host_id", "host_subject_id", "patient_id", "subject_id", "individual_id", "animal_id")
COUNTRY_ALIASES = ("geo_loc_name", "country", "geo_loc")
# Permutation tests can shuffle source labels only within these strata.
STRATA_FIELDS = ("year", "country")
# One-to-one renames from the INSDC list of historical geo_loc_names. Old records keep
# the old name ("Czech Republic"), new ones use the current one ("Czechia"); without
# this the same country would be two strata. Splits and mergers (USSR, Yugoslavia,
# Czechoslovakia, Serbia and Montenegro) have no single successor and are left as is.
HISTORICAL_COUNTRIES = {
    "belgian congo": "Democratic Republic of the Congo",
    "british guiana": "Guyana",
    "burma": "Myanmar",
    "czech republic": "Czechia",
    "east timor": "Timor-Leste",
    "macedonia": "North Macedonia",
    "siam": "Thailand",
    "swaziland": "Eswatini",
    "the former yugoslav republic of macedonia": "North Macedonia",
    "zaire": "Democratic Republic of the Congo",
}
# Thorpe et al. 2022 swept 0-20 SNPs; Hetland et al. 2025 used 22; David et al. 2019 21.
SWEEP_THRESHOLDS = (0, 1, 2, 5, 10, 20, 50, 100)
# Linkage keeps pairs up to this distance (or the threshold, if larger) for the
# histogram and the sweep; pairs further apart are only counted.
CLOSE_LIMIT = 100.0


@dataclass(frozen=True)
class Isolate:
    key: str
    source: str  # "" when Unclassified
    lineage: str
    year: float | None
    host_id: str
    ids: frozenset[str]
    country: str = ""


@dataclass
class Dataset:
    isolates: list[Isolate]
    groups: tuple[str, ...]
    lineage_column: str | None
    rows: int
    duplicates: int
    unclassified: int
    unnamed: int = 0  # rows with no identifier; counted, but never matched to other files
    untyped_keys: set[str] = field(default_factory=set)  # isolates whose only lineage is a "not typed" word
    unread_dates: int = 0  # isolates whose collection_date is in a format decimal_year cannot read

    @property
    def untyped_placeholders(self) -> int:
        """Isolates left untyped because their lineage was NF, Novel, ST 0 or the like."""
        return len(self.untyped_keys)

    @property
    def classified(self) -> list[Isolate]:
        return [iso for iso in self.isolates if iso.source]


# --- building the isolate table -------------------------------------------------


def build_dataset(rows: list[dict[str, str]], by: str = "source", lineage_column: str | None = None) -> Dataset:
    """One record per isolate (rows sharing a BioSample are collapsed)."""
    if by not in {"source", "paper_rollup"}:
        raise InputError("group-by must be 'source' or 'paper_rollup'")
    groups = SOURCE_GROUPS if by == "source" else ROLLUP_GROUPS
    lineage_col = lineage_column or (first_present(rows[0], LINEAGE_ALIASES) if rows else None)
    if lineage_column and rows and lineage_column not in rows[0]:
        raise InputError(f"lineage column {lineage_column!r} is not in the table")
    isolates: dict[str, Isolate] = {}
    unnamed = unread_dates = 0
    untyped: set[str] = set()
    sequence_types = is_sequence_type_column(lineage_col)
    for number, row in enumerate(rows, start=1):
        key = isolate_key(row)
        if not key:
            # No identifier to match on: still an isolate, just never a duplicate.
            unnamed += 1
            key = f"(row {number})"
        if key in isolates:
            # Same BioSample (e.g. two assemblies): keep the first row's values, but every
            # identifier, so distances or genes keyed by either assembly still match.
            isolates[key] = replace(isolates[key], ids=isolates[key].ids | frozenset(_row_names(row)))
            continue
        lineage = _clean(row.get(lineage_col, "")) if lineage_col else ""
        if lineage and is_untyped(lineage, sequence_types):
            untyped.add(key)
            lineage = ""
        year = decimal_year(row.get("collection_date", ""))
        if year is None and _clean(row.get("collection_date", "")):
            unread_dates += 1
        isolates[key] = Isolate(
            key=key,
            source=_source_group(row, by),
            lineage=lineage,
            year=year,
            host_id=_clean(get_field(row, HOST_ID_ALIASES)),
            ids=frozenset(_row_names(row)),
            country=_country(get_field(row, COUNTRY_ALIASES)),
        )
    records = list(isolates.values())
    return Dataset(
        isolates=records,
        groups=groups,
        lineage_column=lineage_col,
        rows=len(rows),
        duplicates=len(rows) - len(records),
        unclassified=sum(1 for iso in records if not iso.source),
        unnamed=unnamed,
        untyped_keys=untyped,
        unread_dates=unread_dates,
    )


# What typing tools write when they could not assign a type: SRST2 "NF" (not found) and
# "ND" (not determined), ARIBA "Novel" and "ND", and words people type by hand. Left as
# they are, every untyped isolate would join one lineage called "NF", and that lineage
# would be reported as shared between sources.
UNTYPED_WORDS = frozenset(
    {"nf", "nd", "nt", "novel", "novel st", "untyped", "untypeable", "untypable", "nontypeable", "nontypable",
     "non typeable", "non typable", "not typed", "not typeable", "not typable", "failed"}
)
# Sequence type 0 is not a sequence type: Kleborate writes 0 when the alleles match no ST.
# Only in an ST column, though; a cluster numbered from zero has a real cluster 0.
UNTYPED_SEQUENCE_TYPES = frozenset({"0"})
SEQUENCE_TYPE_COLUMNS = frozenset({"st", "mlst_st", "sequence_type", "mlst"})


def is_sequence_type_column(name: str | None) -> bool:
    """A column of MLST sequence types, by its name ("ST", "mlst", Kleborate's "..__mlst__ST")."""
    lower = (name or "").strip().lower()
    return lower in SEQUENCE_TYPE_COLUMNS or lower.endswith("__st")


def is_untyped(value: str, sequence_type: bool = False) -> bool:
    """True for a lineage value that says "no type was assigned" rather than naming one.

    SRST2 marks uncertain calls with a trailing * or ?; "NF*" is still not found.
    """
    word = " ".join((value or "").strip().lower().replace("_", " ").replace("-", " ").split()).rstrip("*?").strip()
    return word in UNTYPED_WORDS or (sequence_type and word in UNTYPED_SEQUENCE_TYPES)


def apply_lineages(dataset: Dataset, lineages: dict[str, str]) -> list[str]:
    """Set lineages from an external typing table (Kleborate, mlst, PopPUNK). Returns unmatched names."""
    mapping, unmatched = match_names(list(lineages), dataset)
    by_key: dict[str, str] = {}
    untyped: set[str] = set()
    for name, key in mapping.items():
        value = _clean(lineages[name])
        if not value:
            continue
        if is_untyped(value):
            untyped.add(key)
            continue
        if by_key.get(key, value) != value:
            raise InputError(
                f"isolate {key} has two different lineages ({by_key[key]} and {value}); keep one lineage per isolate"
            )
        by_key[key] = value
    dataset.isolates = [replace(iso, lineage=by_key.get(iso.key, iso.lineage)) for iso in dataset.isolates]
    dataset.lineage_column = dataset.lineage_column or "external table"
    # Counted once, and only while the isolate has no lineage from either the table or this file.
    typed = {iso.key for iso in dataset.isolates if iso.lineage}
    dataset.untyped_keys = (dataset.untyped_keys | untyped) - typed
    return unmatched


def _source_group(row: dict[str, str], by: str) -> str:
    l1 = (row.get("l1") or "").strip()
    if not l1 or l1 == "Unclassified":
        return ""
    if by == "paper_rollup":
        rollup = (row.get("paper_rollup") or "").strip()
        return rollup if rollup in ROLLUP_GROUPS else ""
    if (row.get("healthcare_environment") or "").strip() == HEALTHCARE_YES:
        return "Healthcare environment"
    return l1 if l1 in SOURCE_GROUPS else ""


def _country(geo_loc_name: str) -> str:
    """INSDC geo_loc_name is "Country: region, locality"; keep the current country name."""
    country = _clean(geo_loc_name.split(":", 1)[0])
    return HISTORICAL_COUNTRIES.get(country.lower(), country)


def _clean(value: str) -> str:
    value = (value or "").strip()
    return "" if is_unknown_value(value) else value


def _row_names(row: dict[str, str]) -> set[str]:
    names = set(row_match_ids(row))
    for column in ("strain", "isolate", "sample_name", "name"):
        value = (row.get(column) or "").strip()
        if value:
            names.add(value.upper())
    names.add(isolate_key(row))
    for value in (get_field(row, ASSEMBLY_ALIASES), get_field(row, PDT_ALIASES)):
        if value:
            names.add(value.upper())
    return {name for name in names if name}


_DATE = re.compile(r"^\s*(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?")
# The other INSDC collection_date formats, "15-Mar-2019" and "Mar-2019", and how people
# write the same by hand ("15 March 2019", "March 2019").
_DAY_MONTH_YEAR = re.compile(r"^\s*(?:(\d{1,2})[-\s/.]+)?([A-Za-z]{3,9})\.?[-\s/.,]+(\d{4})\b")
# 15/03/2019 or 03/15/2019: the year is certain, the order of day and month is not.
_NUMERIC_DATE = re.compile(r"^\s*\d{1,2}[/.-]\d{1,2}[/.-](\d{4})\s*$")
_MONTHS = {name: number for number, name in enumerate(
    ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), start=1)}


def decimal_year(text: str) -> float | None:
    """A collection date as a year with a fraction, or None when it has no readable year.

    '2019-03-15' -> 2019.20; '2019-03' -> 2019.21; '2019' -> 2019.5. The INSDC formats
    with a month name read the same way: '15-Mar-2019' -> 2019.20, 'Mar-2019' -> 2019.21.
    A numeric date that ends in the year ('15/03/2019', '03/15/2019') gives 2019.5, the
    middle of the year, because its day and month cannot be told apart.
    """
    text = text or ""
    match = _DATE.match(text)
    if match:
        year, month, day = int(match.group(1)), match.group(2), match.group(3)
    else:
        named = _DAY_MONTH_YEAR.match(text)
        numeric = None if named else _NUMERIC_DATE.match(text)
        if named and named.group(2)[:3].lower() in _MONTHS:
            year, month, day = int(named.group(3)), str(_MONTHS[named.group(2)[:3].lower()]), named.group(1)
        elif numeric:
            year, month, day = int(numeric.group(1)), None, None
        else:
            return None
    if not 1900 <= year <= 2100:
        return None
    try:
        if month and day:
            start = date(year, 1, 1).toordinal()
            length = date(year + 1, 1, 1).toordinal() - start
            return year + (date(year, int(month), int(day)).toordinal() - start + 0.5) / length
        if month and 1 <= int(month) <= 12:
            return year + (int(month) - 0.5) / 12
    except ValueError:
        pass
    return year + 0.5


# --- name matching for genotype and distance files -------------------------------

_ASSEMBLY_FILE = re.compile(r"^(GC[AF]_\d+\.\d+)_", re.IGNORECASE)
# Sequence files as typing tools name them (mlst's FILE column, for one), optionally compressed.
_SEQUENCE_FILE = re.compile(r"\.(fa|fas|fasta|fna|ffn|fsa|gb|gbk|gbff|genbank|embl)(\.(gz|bz2|xz|zip))?$", re.IGNORECASE)


def _name_candidates(name: str) -> list[str]:
    """What an external name may stand for, most specific first: the name itself, an NCBI
    assembly file's accession, a sequence file's name without its directory and extension,
    and the file's directory, as in shovill's ISO1/contigs.fa."""
    name = name.strip()
    parts = [part for part in re.split(r"[/\\]", name) if part]
    base = parts[-1] if parts else name
    candidates = [name]
    stem = _ASSEMBLY_FILE.match(base)
    if stem:
        candidates.append(stem.group(1))
    if _SEQUENCE_FILE.search(base):
        candidates.append(_SEQUENCE_FILE.sub("", base))
        if len(parts) > 1:
            candidates.append(parts[-2])
    return list(dict.fromkeys(c.upper() for c in candidates if c))


def match_names(names: list[str], dataset: Dataset) -> tuple[dict[str, str], list[str]]:
    """Map external sample names to isolate keys.

    The most specific candidate that names any isolate decides (see _name_candidates). If it
    names more than one, the name is ambiguous and left unmatched rather than guessed.
    """
    owners: dict[str, set[str]] = {}
    for iso in dataset.isolates:
        for ident in iso.ids:
            owners.setdefault(ident, set()).add(iso.key)
    mapping: dict[str, str] = {}
    unmatched: list[str] = []
    for name in names:
        found = next((owners[c] for c in _name_candidates(name) if c in owners), set())
        if len(found) == 1:
            mapping[name] = next(iter(found))
        else:
            unmatched.append(name)
    return mapping, unmatched


# --- statistics helpers ------------------------------------------------------------


# The 97.5% point of the standard normal distribution, for two-sided 95% intervals.
Z_95 = 1.959963984540054


def wilson(successes: int, n: int, z: float = Z_95) -> tuple[float, float]:
    """Wilson score interval for a proportion (95% by default), as R's prop.test(correct=FALSE)."""
    if n == 0:
        return math.nan, math.nan
    p = successes / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """BH-adjusted q-values, same order as the input. NaN stays NaN."""
    indexed = [(p, i) for i, p in enumerate(pvalues) if not math.isnan(p)]
    q = [math.nan] * len(pvalues)
    m = len(indexed)
    running = 1.0
    for rank, (p, i) in reversed(list(enumerate(sorted(indexed), start=1))):
        running = min(running, p * m / rank)
        q[i] = running
    return q


def perm_pvalue(exceed: int, permutations: int) -> float:
    """Phipson & Smyth (2010): never zero. NaN when no permutations ran, since then nothing was tested."""
    if permutations < 1:
        return math.nan
    return (exceed + 1) / (permutations + 1)


def strata_codes(isolates: list[Isolate], by: tuple[str, ...] = ()) -> np.ndarray | None:
    """Integer stratum per isolate from collection year and/or country; None when not stratifying.

    Isolates missing a value share an "unknown" stratum for that field.
    """
    if not by:
        return None
    unknown = [f for f in by if f not in STRATA_FIELDS]
    if unknown:
        raise InputError(f"cannot stratify by {unknown}; choose from {list(STRATA_FIELDS)}")

    def part(iso: Isolate, name: str) -> str:
        if name == "year":
            return str(int(iso.year)) if iso.year is not None else "?"
        return iso.country.lower() or "?"

    keys = ["|".join(part(iso, name) for name in by) for iso in isolates]
    index = {key: i for i, key in enumerate(sorted(set(keys)))}
    return np.array([index[key] for key in keys], dtype=np.int64)


def permute_within(labels: np.ndarray, strata: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    """Shuffle labels, but only among isolates of the same stratum."""
    if strata is None:
        return rng.permutation(labels)
    slots = np.argsort(strata, kind="stable")  # positions grouped by stratum
    draws = np.lexsort((rng.random(len(labels)), strata))  # same grouping, random order inside
    out = np.empty_like(labels)
    out[slots] = labels[draws]
    return out


def strata_summary(isolates: list[Isolate], by: tuple[str, ...]) -> dict:
    """How informative the strata are: isolates in strata with one source cannot change label."""
    codes = strata_codes([iso for iso in isolates if iso.source], by)
    if codes is None:
        return {"strata": 0, "isolates": 0, "isolates_in_single_source_strata": 0}
    classified = [iso for iso in isolates if iso.source]
    sources: dict[int, set[str]] = {}
    for iso, code in zip(classified, codes.tolist()):
        sources.setdefault(code, set()).add(iso.source)
    fixed = sum(1 for code in codes.tolist() if len(sources[code]) == 1)
    return {"strata": len(sources), "isolates": len(classified), "isolates_in_single_source_strata": fixed}


def simpson_variance(counts: list[int]) -> float:
    """Sampling variance of Simpson's index of diversity (Simpson 1949), exact for any sample size.

    The variance of the index over samples of n isolates drawn from a population in which
    the lineages have the frequencies seen here.
    """
    n = sum(counts)
    if n < 2:
        return math.nan
    p = [c / n for c in counts]
    s2, s3 = sum(x**2 for x in p), sum(x**3 for x in p)
    return max((4 * (n - 2) * s3 + 2 * s2 - 2 * (2 * n - 3) * s2**2) / (n * (n - 1)), 0.0)


def simpson(counts: list[int]) -> tuple[float, float, float]:
    """Simpson's index of diversity (Hunter & Gaston 1988) with an approximate 95% CI.

    The interval is D +/- 2 standard errors, as in Grundmann et al. (2001), but with the
    exact sampling variance of the index (Simpson 1949) and not its large-sample form
    4/n (sum p^3 - (sum p^2)^2). The two agree for large samples. The large-sample form
    is exactly zero whenever the lineages are equally common, so ten isolates of ten
    lineages had D = 1.00 with the interval 1.00-1.00.
    """
    n = sum(counts)
    if n < 2:
        return math.nan, math.nan, math.nan
    d = 1 - sum(c * (c - 1) for c in counts) / (n * (n - 1))
    half = 2 * math.sqrt(simpson_variance(counts))
    return d, max(0.0, d - half), min(1.0, d + half)


def rarefaction(counts: list[int], sizes: list[int]) -> list[float]:
    """Expected number of lineages in a random subsample of n isolates (Hurlbert 1971)."""
    total = sum(counts)

    def log_comb(a: int, b: int) -> float:
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

    expected = []
    for n in sizes:
        e = 0.0
        for c in counts:
            e += 1.0 if total - c < n else 1 - math.exp(log_comb(total - c, n) - log_comb(total, n))
        expected.append(e)
    return expected


# --- 1. sampling -------------------------------------------------------------------


def source_counts(dataset: Dataset) -> list[dict]:
    counts = Counter(iso.source for iso in dataset.classified)
    total = sum(counts.values())
    return [
        {"source": g, "isolates": counts[g], "percent": 100 * counts[g] / total if total else 0.0}
        for g in dataset.groups
        if counts[g]
    ]


def counts_by_year(dataset: Dataset) -> list[dict]:
    rows: dict[int, Counter] = {}
    for iso in dataset.classified:
        if iso.year is not None:
            rows.setdefault(int(iso.year), Counter())[iso.source] += 1
    return [{"year": year, **{g: rows[year][g] for g in dataset.groups}} for year in sorted(rows)]


# --- 2. lineage diversity, sharing and enrichment ----------------------------------


def _typed(dataset: Dataset) -> list[Isolate]:
    return [iso for iso in dataset.classified if iso.lineage]


def diversity(dataset: Dataset) -> list[dict]:
    typed = _typed(dataset)
    out = []
    for group in dataset.groups:
        counts = list(Counter(iso.lineage for iso in typed if iso.source == group).values())
        if not counts:
            continue
        d, lo, hi = simpson(counts)
        out.append(
            {"source": group, "isolates": sum(counts), "lineages": len(counts), "simpson": d, "ci_low": lo, "ci_high": hi}
        )
    return out


def rarefaction_curves(dataset: Dataset, points: int = 25) -> list[dict]:
    typed = _typed(dataset)
    out = []
    for group in dataset.groups:
        counts = list(Counter(iso.lineage for iso in typed if iso.source == group).values())
        total = sum(counts)
        if total < 2:
            continue
        sizes = sorted({max(1, round(total * k / points)) for k in range(1, points + 1)})
        for n, e in zip(sizes, rarefaction(counts, sizes)):
            out.append({"source": group, "isolates_sampled": n, "expected_lineages": e})
    return out


def lineage_sources(dataset: Dataset) -> list[dict]:
    """Per lineage: isolates per source, ordered by size."""
    table: dict[str, Counter] = {}
    for iso in _typed(dataset):
        table.setdefault(iso.lineage, Counter())[iso.source] += 1
    rows = []
    for lineage, counts in table.items():
        present = [g for g in dataset.groups if counts[g]]
        rows.append(
            {"lineage": lineage, "isolates": sum(counts.values()), "sources": len(present), **{g: counts[g] for g in dataset.groups}}
        )
    return sorted(rows, key=lambda r: (-r["isolates"], r["lineage"]))


def lineage_sharing_summary(dataset: Dataset) -> dict:
    """Hetland et al. 2025: share of lineages in >= 2 sources and the isolates they hold."""
    rows = lineage_sources(dataset)
    typed = sum(r["isolates"] for r in rows)
    shared = [r for r in rows if r["sources"] >= 2]
    return {
        "lineages": len(rows),
        "shared_lineages": len(shared),
        "typed_isolates": typed,
        "isolates_in_shared": sum(r["isolates"] for r in shared),
    }


def lineage_intersections(dataset: Dataset, top: int = 15) -> list[dict]:
    """UpSet-style: how many lineages occur in exactly each combination of sources."""
    combos: Counter = Counter()
    for row in lineage_sources(dataset):
        combos[tuple(g for g in dataset.groups if row[g])] += 1
    ordered = sorted(combos.items(), key=lambda item: (-item[1], len(item[0]), item[0]))
    return [{"sources": " + ".join(combo), "members": combo, "lineages": n} for combo, n in ordered[:top]]


def lineage_enrichment(
    dataset: Dataset, min_isolates: int = 5, permutations: int = 9999, seed: int = 1, strata: tuple[str, ...] = ()
) -> list[dict]:
    """Thorpe et al. 2022: is a lineage over- or under-represented in a source?

    Source labels are permuted across typed isolates (within strata, if given);
    two-sided p-values use |count - expected|, are never zero (Phipson & Smyth 2010)
    and are BH-adjusted. The expected count is exact: the sum over strata of
    lineage size x source size / stratum size.
    """
    typed = _typed(dataset)
    groups = [g for g in dataset.groups if any(iso.source == g for iso in typed)]
    sizes = Counter(iso.lineage for iso in typed)
    tested = sorted((lin for lin, n in sizes.items() if n >= min_isolates), key=lambda lin: (-sizes[lin], lin))
    if not tested or len(groups) < 2:
        return []
    lin_index = {lin: i for i, lin in enumerate(tested)}
    src_index = {g: i for i, g in enumerate(groups)}
    lin = np.array([lin_index.get(iso.lineage, -1) for iso in typed])
    src = np.array([src_index[iso.source] for iso in typed])
    keep = lin >= 0
    n_l, n_s = len(tested), len(groups)

    def counts(labels: np.ndarray) -> np.ndarray:
        return np.bincount(lin[keep] * n_s + labels[keep], minlength=n_l * n_s).reshape(n_l, n_s)

    observed = counts(src)
    codes = strata_codes(typed, strata)
    stratum = codes if codes is not None else np.zeros(len(typed), dtype=np.int64)
    expected = np.zeros((n_l, n_s))
    for k in np.unique(stratum):
        inside = stratum == k
        lineage_k = np.bincount(lin[inside & keep], minlength=n_l)
        source_k = np.bincount(src[inside], minlength=n_s)
        expected += np.outer(lineage_k, source_k) / inside.sum()
    deviation = np.abs(observed - expected) - 1e-9
    exceed = np.zeros_like(observed)
    rng = np.random.default_rng(seed)
    for _ in range(permutations):
        exceed += np.abs(counts(permute_within(src, codes, rng)) - expected) >= deviation
    rows = []
    for i, lineage in enumerate(tested):
        for j, group in enumerate(groups):
            rows.append(
                {
                    "lineage": lineage,
                    "source": group,
                    "observed": int(observed[i, j]),
                    "expected": float(expected[i, j]),
                    "ratio": float(observed[i, j] / expected[i, j]) if expected[i, j] else math.nan,
                    "p": perm_pvalue(int(exceed[i, j]), permutations),
                }
            )
    for row, q in zip(rows, benjamini_hochberg([r["p"] for r in rows])):
        row["q"] = q
        row["direction"] = "enriched" if row["observed"] > row["expected"] else "depleted"
    return rows


# --- 3. resistance / virulence / plasmid features by source --------------------------


@dataclass
class FeatureResult:
    groups: list[str]
    genotyped: dict[str, int]
    rows: list[dict] = field(default_factory=list)  # one per feature x source
    tests: list[dict] = field(default_factory=list)  # one per feature
    unmatched: list[str] = field(default_factory=list)


def feature_prevalence(
    dataset: Dataset,
    presence: dict[str, set[str]],
    classes: dict[str, str],
    permutations: int = 9999,
    seed: int = 1,
    strata: tuple[str, ...] = (),
) -> FeatureResult:
    """Prevalence (Wilson 95% CI) of each feature per source, and a permutation
    chi-squared test of independence from source with BH correction."""
    mapping, unmatched = match_names(list(presence), dataset)
    by_key: dict[str, set[str]] = {}
    for name, key in mapping.items():
        # Share the caller's set; only two names for one isolate need a new (merged) one.
        by_key[key] = by_key[key] | presence[name] if key in by_key else presence[name]
    isolates = [iso for iso in dataset.classified if iso.key in by_key]
    groups = [g for g in dataset.groups if any(iso.source == g for iso in isolates)]
    result = FeatureResult(groups=groups, genotyped={g: 0 for g in groups}, unmatched=unmatched)
    if not isolates:
        return result
    features = sorted({f for iso in isolates for f in by_key[iso.key]})
    column = {f: k for k, f in enumerate(features)}
    src_index = {g: i for i, g in enumerate(groups)}
    src = np.array([src_index[iso.source] for iso in isolates])
    # Float rather than int: numpy multiplies float matrices with BLAS, about 35x faster
    # in the permutation loop below, and whole-number counts stay exact.
    matrix = np.zeros((len(isolates), len(features)))
    rows = np.repeat(np.arange(len(isolates)), [len(by_key[iso.key]) for iso in isolates])
    matrix[rows, [column[f] for iso in isolates for f in by_key[iso.key]]] = 1.0
    n_s = np.bincount(src, minlength=len(groups))
    result.genotyped = {g: int(n_s[i]) for i, g in enumerate(groups)}
    onehot = np.eye(len(groups))
    observed = onehot[src].T @ matrix  # sources x features

    totals = matrix.sum(axis=0)
    testable = (totals > 0) & (totals < len(isolates)) & (len(groups) > 1)
    stat_obs = _chi2(observed, n_s, totals)
    exceed = np.zeros(len(features), dtype=np.int64)
    if testable.any():
        rng = np.random.default_rng(seed)
        codes = strata_codes(isolates, strata)
        sub = matrix[:, testable]
        # Many permutations per matrix product: one pass over the isolate x feature
        # matrix per batch instead of one per permutation. Permutations are drawn in
        # the same order, so results do not depend on the batch size.
        batch = max(1, min(256, 2_000_000 // (len(groups) * len(isolates))))
        done = 0
        while done < permutations:
            size = min(batch, permutations - done)
            labels = np.stack([permute_within(src, codes, rng) for _ in range(size)])  # size x isolates
            rows = onehot[labels].transpose(0, 2, 1).reshape(size * len(groups), len(isolates))
            counts = (rows @ sub).reshape(size, len(groups), -1)  # size x sources x features
            exceed[testable] += (_chi2(counts, n_s, totals[testable]) >= stat_obs[testable] - 1e-9).sum(axis=0)
            done += size
    pvalues = [perm_pvalue(int(exceed[k]), permutations) if testable[k] else math.nan for k in range(len(features))]
    qvalues = benjamini_hochberg(pvalues)
    for k, feature in enumerate(features):
        result.tests.append(
            {
                "feature": feature,
                "class": classes.get(feature, ""),
                "isolates_with_feature": int(totals[k]),
                "percent_overall": float(100 * totals[k] / len(isolates)),
                "chi2": float(stat_obs[k]),
                "p": pvalues[k],
                "q": qvalues[k],
            }
        )
        for i, group in enumerate(groups):
            lo, hi = wilson(int(observed[i, k]), int(n_s[i]))
            result.rows.append(
                {
                    "feature": feature,
                    "class": classes.get(feature, ""),
                    "source": group,
                    "present": int(observed[i, k]),
                    "genotyped": int(n_s[i]),
                    "percent": float(100 * observed[i, k] / n_s[i]) if n_s[i] else math.nan,
                    "ci_low": 100 * lo,
                    "ci_high": 100 * hi,
                }
            )
    return result


def _chi2(present: np.ndarray, n_s: np.ndarray, totals: np.ndarray) -> np.ndarray:
    """Pearson chi-squared for each feature's 2 x sources table.

    `present` is sources x features, or a stack of those (one per permutation).
    """
    n = n_s.sum()
    expected = np.outer(n_s, totals) / n
    absent = n_s[:, None] - present
    expected_absent = np.outer(n_s, n - totals) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = (present - expected) ** 2 / expected + (absent - expected_absent) ** 2 / expected_absent
    return np.nansum(np.where(np.isfinite(terms), terms, 0.0), axis=-2)


# --- 4. genomic linkage between sources ------------------------------------------------


@dataclass
class LinkageResult:
    threshold: float
    groups: list[str]
    matched: int
    unmatched: list[str]
    per_source: dict[str, int]
    events: list[dict]  # per source pair
    linked_fraction: list[dict]  # per (from source, to source)
    sweep: list[dict]
    clusters: list[dict]
    cross_edges: list[dict]
    close_within: np.ndarray  # same-source distances up to close_limit
    close_between: np.ndarray  # cross-source distances up to close_limit
    pairs_within: int  # every same-source pair with a distance
    pairs_between: int  # every cross-source pair with a distance
    close_limit: float


def linkage(
    dataset: Dataset,
    names: list[str],
    matrix: np.ndarray,
    threshold: float = 21,
    permutations: int = 9999,
    seed: int = 1,
    strata: tuple[str, ...] = (),
) -> LinkageResult:
    """Strain sharing within and between sources at a SNP threshold.

    Following Thorpe et al. 2022 and Hetland et al. 2025: isolates within the
    threshold are joined by single linkage; each unique source pair connected by
    at least one qualifying pair inside a cluster is one sharing event. Pairs from
    the same host individual are excluded. Events are normalised per 100 isolates
    of the two sources, and compared with a null in which source labels are
    permuted across isolates (sampling effort held fixed; cf. Muloi et al. 2022).
    With `strata`, labels are shuffled only among isolates from the same collection
    year and/or country, so sources sampled at different times or places are not
    expected to share strains.
    """
    isolates, index, unmatched = _matched_isolates(dataset, names, matrix)
    groups = [g for g in dataset.groups if any(iso.source == g for iso in isolates)]
    n_s = len(groups)
    src = np.array([groups.index(iso.source) if iso.source in groups else -1 for iso in isolates], dtype=np.int64)
    per_source = {g: int((src == k).sum()) for k, g in enumerate(groups)}

    # Only pairs within reach of a threshold are kept; holding all n^2/2 pairs at once
    # was what made large matrices run out of memory.
    close_limit = max(threshold, CLOSE_LIMIT, *SWEEP_THRESHOLDS)
    pair_a, pair_b, pair_d = _close_pairs(matrix, index, close_limit)
    pairs_within, pairs_between = _pairs_by_source(matrix, index, src)
    classified = (src[pair_a] >= 0) & (src[pair_b] >= 0)
    same = classified & (src[pair_a] == src[pair_b])
    host_codes: dict[str, int] = {}
    host = np.array([host_codes.setdefault(iso.host_id, len(host_codes) + 1) if iso.host_id else 0 for iso in isolates],
                    dtype=np.int64)
    same_host = (host[pair_a] > 0) & (host[pair_a] == host[pair_b])

    def edges_at(t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        hit = pair_d <= t
        # Clusters include pairs from one host individual, so a patient whose two
        # isolates link to different sources stays one cluster (one sharing event per
        # source pair); only the pair itself is not counted as sharing.
        clusters = _components(len(isolates), pair_a[hit], pair_b[hit])
        keep = hit & ~same_host
        return pair_a[keep], pair_b[keep], pair_d[keep], clusters

    def event_counts(labels: np.ndarray, a: np.ndarray, b: np.ndarray, clusters: np.ndarray) -> np.ndarray:
        la, lb = labels[a], labels[b]
        ok = (la >= 0) & (lb >= 0)
        lo, hi = np.minimum(la[ok], lb[ok]), np.maximum(la[ok], lb[ok])
        codes = np.unique(clusters[a[ok]] * n_s * n_s + lo * n_s + hi)
        return np.bincount(codes % (n_s * n_s), minlength=n_s * n_s).reshape(n_s, n_s)

    a, b, d, clusters = edges_at(threshold)
    observed = event_counts(src, a, b, clusters) if n_s else np.zeros((0, 0), dtype=np.int64)

    # Permutation null: shuffle labels among classified isolates, graph fixed.
    exceed = np.zeros_like(observed)
    perm_values = []
    labelled = np.flatnonzero(src >= 0)
    codes = strata_codes([isolates[i] for i in labelled], strata)
    if n_s and len(a) and permutations > 0:
        rng = np.random.default_rng(seed)
        draws = []
        for _ in range(permutations):
            labels = src.copy()
            labels[labelled] = permute_within(src[labelled], codes, rng)
            draws.append(event_counts(labels, a, b, clusters))
        perm_values = np.stack(draws)
        mean = perm_values.mean(axis=0)
        exceed = (np.abs(perm_values - mean) >= np.abs(observed - mean) - 1e-9).sum(axis=0)
    events = []
    for x in range(n_s):
        for y in range(x, n_s):
            denominator = per_source[groups[x]] + (per_source[groups[y]] if x != y else 0)
            row = {
                "source_a": groups[x],
                "source_b": groups[y],
                "kind": "within" if x == y else "between",
                "events": int(observed[x, y]),
                "isolates_in_pair": denominator,
                "events_per_100": float(100 * observed[x, y] / denominator) if denominator else math.nan,
            }
            if len(perm_values):
                column = perm_values[:, x, y]
                row.update(
                    expected=float(column.mean()),
                    expected_low=float(np.percentile(column, 2.5)),
                    expected_high=float(np.percentile(column, 97.5)),
                    p=perm_pvalue(int(exceed[x, y]), permutations),
                )
            else:
                row.update(expected=math.nan, expected_low=math.nan, expected_high=math.nan, p=math.nan)
            events.append(row)
    for row, q in zip(events, benjamini_hochberg([r["p"] for r in events])):
        row["q"] = q

    # Hetland et al. 2025: share of isolates of one source with a close relative in another.
    neighbours = [set() for _ in isolates]
    for x, y in zip(a, b):
        if src[y] >= 0:
            neighbours[x].add(int(src[y]))
        if src[x] >= 0:
            neighbours[y].add(int(src[x]))
    linked = []
    for k, group in enumerate(groups):
        members = np.flatnonzero(src == k)
        for m, other in enumerate(groups):
            count = sum(1 for i in members if m in neighbours[i])
            linked.append(
                {
                    "from_source": group,
                    "to_source": other,
                    "isolates": len(members),
                    "with_close_relative": count,
                    "percent": 100 * count / len(members) if len(members) else math.nan,
                }
            )

    sweep = []
    for t in sorted(set(SWEEP_THRESHOLDS) | {threshold}):
        ta, tb, _, tc = edges_at(t)
        counts = event_counts(src, ta, tb, tc) if n_s else np.zeros((0, 0))
        diagonal = int(np.trace(counts)) if n_s else 0
        sweep.append({"threshold": t, "within_events": diagonal, "between_events": int(counts.sum()) - diagonal})

    cross_edges = []
    for x, y, distance in zip(a, b, d):
        if src[x] >= 0 and src[y] >= 0 and src[x] != src[y]:
            gap = (
                abs(isolates[x].year - isolates[y].year)
                if isolates[x].year is not None and isolates[y].year is not None
                else math.nan
            )
            cross_edges.append(
                {
                    "cluster": int(clusters[x]),
                    "isolate_a": isolates[x].key,
                    "source_a": isolates[x].source,
                    "isolate_b": isolates[y].key,
                    "source_b": isolates[y].source,
                    "lineage": isolates[x].lineage or isolates[y].lineage,
                    "snps": float(distance),
                    "years_apart": gap,
                }
            )

    # Members and edges grouped by cluster in one pass (not one scan of all edges per cluster).
    cluster_members: dict[int, list[int]] = {}
    for i, cid in enumerate(clusters.tolist()):
        cluster_members.setdefault(cid, []).append(i)
    cluster_edges: dict[int, list[tuple[int, int, float]]] = {}
    for x, y, distance in zip(a.tolist(), b.tolist(), d.tolist()):
        cluster_edges.setdefault(int(clusters[x]), []).append((x, y, distance))
    cluster_rows = []
    for cid in sorted(cluster_members):
        members = np.array(cluster_members[cid])
        if len(members) < 2:
            continue
        sources = Counter(isolates[i].source or "Unclassified" for i in members)
        inside = cluster_edges.get(cid, [])
        cluster_rows.append(
            {
                "cluster": int(cid),
                "isolates": len(members),
                "sources": ", ".join(f"{g} ({sources[g]})" for g in [*groups, "Unclassified"] if sources[g]),
                "n_sources": len([g for g in groups if sources[g]]),
                "lineages": ", ".join(sorted({isolates[i].lineage for i in members if isolates[i].lineage})),
                "members": [isolates[i].key for i in members],
                "member_sources": [isolates[i].source for i in members],
                "edges": [(int(np.searchsorted(members, x)), int(np.searchsorted(members, y)), distance)
                          for x, y, distance in inside],
            }
        )
    cluster_rows.sort(key=lambda r: (-r["n_sources"], -r["isolates"], r["cluster"]))
    return LinkageResult(
        threshold=threshold,
        groups=groups,
        matched=len(isolates),
        unmatched=unmatched,
        per_source=per_source,
        events=events,
        linked_fraction=linked,
        sweep=sweep,
        clusters=cluster_rows,
        cross_edges=cross_edges,
        close_within=pair_d[same],
        close_between=pair_d[classified & ~same],
        pairs_within=pairs_within,
        pairs_between=pairs_between,
        close_limit=close_limit,
    )


def fractional_distances(matrix: np.ndarray, block: int = 512) -> bool:
    """True when the distances are fractions no larger than 1, as Mash or ANI distances are.

    SNP counts are whole numbers. A SNP threshold applied to such a matrix links every
    pair of isolates. Read a block of rows at a time, so the matrix is not copied.
    """
    fractional = False
    for start in range(0, len(matrix), block):
        rows = matrix[start : start + block]
        finite = rows[np.isfinite(rows)]
        if finite.size and finite.max() > 1:
            return False
        fractional = fractional or bool((finite != np.round(finite)).any())
    return fractional


def _matched_isolates(
    dataset: Dataset, names: list[str], matrix: np.ndarray
) -> tuple[list[Isolate], np.ndarray, list[str]]:
    """Isolates found in a distance matrix (first occurrence of each), their rows in it, and unmatched names."""
    mapping, unmatched = match_names(names, dataset)
    by_key = {iso.key: iso for iso in dataset.isolates}
    index: list[int] = []
    isolates: list[Isolate] = []
    seen: set[str] = set()
    for i, name in enumerate(names):
        key = mapping.get(name)
        if key is None or key in seen:
            continue
        seen.add(key)
        index.append(i)
        isolates.append(by_key[key])
    return isolates, np.array(index, dtype=np.int64), unmatched


def _upper_rows(matrix: np.ndarray, index: np.ndarray) -> Iterator[tuple[int, np.ndarray]]:
    """For each isolate i (a position in `index`), its distances to isolates i+1, i+2, ...

    Rows are read from the full matrix one at a time, so nothing of size n^2 is copied.
    """
    contiguous = np.array_equal(index, np.arange(len(index)))  # the matrix's first rows, in order
    for i in range(len(index) - 1):
        yield i, (matrix[i, i + 1 : len(index)] if contiguous else matrix[index[i], index[i + 1 :]])


def _close_pairs(matrix: np.ndarray, index: np.ndarray, cutoff: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pairs i < j with a finite distance <= cutoff, in row order, and their distances.

    Built row by row, so memory grows with the number of close pairs, not with n^2.
    """
    a_parts, b_parts, d_parts = [], [], []
    for i, row in _upper_rows(matrix, index):
        j = np.flatnonzero(np.isfinite(row) & (row <= cutoff))
        if len(j):
            a_parts.append(np.full(len(j), i, dtype=np.int64))
            b_parts.append(j + i + 1)
            d_parts.append(row[j])
    if not a_parts:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0)
    return np.concatenate(a_parts), np.concatenate(b_parts), np.concatenate(d_parts)


def _pairs_by_source(matrix: np.ndarray, index: np.ndarray, src: np.ndarray) -> tuple[int, int]:
    """Same-source and cross-source pairs that have a distance (both isolates classified)."""
    within = between = 0
    for i, row in _upper_rows(matrix, index):
        if src[i] < 0:
            continue
        others = src[i + 1 :]
        usable = np.isfinite(row) & (others >= 0)
        same = int((usable & (others == src[i])).sum())
        within += same
        between += int(usable.sum()) - same
    return within, between


@dataclass
class Dereplication:
    threshold: float
    isolates_before: int
    isolates_after: int
    clusters_collapsed: int  # clusters with more than one isolate
    not_in_matrix: int  # classified isolates kept as-is because they have no distances


def dereplicate(dataset: Dataset, names: list[str], matrix: np.ndarray, threshold: float) -> tuple[Dataset, Dereplication]:
    """Keep one isolate per genomic cluster per source.

    Isolates from the same source within `threshold` SNPs are joined by single
    linkage and represented by one isolate (the lowest key, so the choice is
    reproducible). This removes the pseudo-replication that outbreak clusters add
    to per-isolate tests, as in the IMP-cluster dereplication of Vezina et al. 2025.
    Isolates from different sources are never merged, so cross-source sharing is
    not hidden. Isolates without distances are kept unchanged.
    """
    isolates, index, _ = _matched_isolates(dataset, names, matrix)
    a, b, _ = _close_pairs(matrix, index, threshold)
    sources = np.array([iso.source for iso in isolates], dtype=object)
    link = (sources[a] != "") & (sources[a] == sources[b])
    clusters = _components(len(isolates), a[link], b[link])
    members: dict[int, list[str]] = {}
    for iso, cid in zip(isolates, clusters.tolist()):
        members.setdefault(cid, []).append(iso.key)
    dropped = {key for keys in members.values() for key in sorted(keys)[1:]}
    kept = [iso for iso in dataset.isolates if iso.key not in dropped]
    in_matrix = {iso.key for iso in isolates}
    result = Dataset(
        isolates=kept,
        groups=dataset.groups,
        lineage_column=dataset.lineage_column,
        rows=dataset.rows,
        duplicates=dataset.duplicates,
        unclassified=sum(1 for iso in kept if not iso.source),
        unnamed=dataset.unnamed,
        untyped_keys=set(dataset.untyped_keys),
        unread_dates=dataset.unread_dates,
    )
    summary = Dereplication(
        threshold=threshold,
        isolates_before=len(dataset.classified),
        isolates_after=len(result.classified),
        clusters_collapsed=sum(1 for keys in members.values() if len(keys) > 1),
        not_in_matrix=sum(1 for iso in dataset.classified if iso.key not in in_matrix),
    )
    return result, summary


def _components(n: int, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Connected components (single linkage) with union-find."""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for x, y in zip(a.tolist(), b.tolist()):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)
    return np.array([find(x) for x in range(n)], dtype=np.int64)
