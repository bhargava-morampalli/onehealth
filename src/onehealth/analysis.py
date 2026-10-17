"""Standard One Health population analyses on a classified table.

The analysis set follows recent One Health genomic studies (see ROADMAP.md):
source composition, lineage diversity and sharing, lineage-source enrichment,
resistance/virulence prevalence by source, and genomic linkage between sources
at a SNP threshold with a permutation null. Every function is pure and returns
plain rows (dicts) so each number in the report is also written to a table.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass, field
from datetime import date

import numpy as np

from onehealth.constants import ASSEMBLY_ALIASES, HEALTHCARE_YES, PDT_ALIASES
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
# Thorpe et al. 2022 swept 0-20 SNPs; Hetland et al. 2025 used 22; David et al. 2019 21.
SWEEP_THRESHOLDS = (0, 1, 2, 5, 10, 20, 50, 100)


@dataclass(frozen=True)
class Isolate:
    key: str
    source: str  # "" when Unclassified
    lineage: str
    year: float | None
    host_id: str
    ids: frozenset[str]


@dataclass
class Dataset:
    isolates: list[Isolate]
    groups: tuple[str, ...]
    lineage_column: str | None
    rows: int
    duplicates: int
    unclassified: int

    @property
    def classified(self) -> list[Isolate]:
        return [iso for iso in self.isolates if iso.source]


# --- building the isolate table -------------------------------------------------


def build_dataset(rows: list[dict[str, str]], by: str = "source", lineage_column: str | None = None) -> Dataset:
    """One record per isolate (rows sharing a BioSample are collapsed)."""
    if by not in {"source", "paper_rollup"}:
        raise ValueError("group-by must be 'source' or 'paper_rollup'")
    groups = SOURCE_GROUPS if by == "source" else ROLLUP_GROUPS
    lineage_col = lineage_column or (first_present(rows[0], LINEAGE_ALIASES) if rows else None)
    if lineage_column and rows and lineage_column not in rows[0]:
        raise ValueError(f"lineage column {lineage_column!r} is not in the table")
    isolates: dict[str, Isolate] = {}
    for row in rows:
        key = isolate_key(row)
        if not key or key in isolates:
            continue
        isolates[key] = Isolate(
            key=key,
            source=_source_group(row, by),
            lineage=_clean(row.get(lineage_col, "")) if lineage_col else "",
            year=decimal_year(row.get("collection_date", "")),
            host_id=_clean(get_field(row, HOST_ID_ALIASES)),
            ids=frozenset(_row_names(row)),
        )
    records = list(isolates.values())
    return Dataset(
        isolates=records,
        groups=groups,
        lineage_column=lineage_col,
        rows=len(rows),
        duplicates=len(rows) - len(records),
        unclassified=sum(1 for iso in records if not iso.source),
    )


def apply_lineages(dataset: Dataset, lineages: dict[str, str]) -> list[str]:
    """Set lineages from an external typing table (Kleborate, mlst, PopPUNK). Returns unmatched names."""
    mapping, unmatched = match_names(list(lineages), dataset)
    by_key = {key: _clean(lineages[name]) for name, key in mapping.items()}
    dataset.isolates = [
        Isolate(iso.key, iso.source, by_key.get(iso.key, iso.lineage), iso.year, iso.host_id, iso.ids)
        for iso in dataset.isolates
    ]
    dataset.lineage_column = dataset.lineage_column or "external table"
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


def decimal_year(text: str) -> float | None:
    """'2019-03-15' -> 2019.20; '2019-03' -> 2019.21; '2019' -> 2019.5; else None."""
    match = _DATE.match(text or "")
    if not match:
        return None
    year = int(match.group(1))
    if not 1900 <= year <= 2100:
        return None
    month, day = match.group(2), match.group(3)
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


def match_names(names: list[str], dataset: Dataset) -> tuple[dict[str, str], list[str]]:
    """Map external sample names to isolate keys. Ambiguous identifiers are not used."""
    owners: dict[str, set[str]] = {}
    for iso in dataset.isolates:
        for ident in iso.ids:
            owners.setdefault(ident, set()).add(iso.key)
    mapping: dict[str, str] = {}
    unmatched: list[str] = []
    for name in names:
        candidates = [name.strip().upper()]
        stem = _ASSEMBLY_FILE.match(name.strip())
        if stem:
            candidates.append(stem.group(1).upper())
        key = next((next(iter(owners[c])) for c in candidates if len(owners.get(c, ())) == 1), None)
        if key is None:
            unmatched.append(name)
        else:
            mapping[name] = key
    return mapping, unmatched


# --- statistics helpers ------------------------------------------------------------


def wilson(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a proportion."""
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
    """Phipson & Smyth (2010): never zero."""
    return (exceed + 1) / (permutations + 1)


def simpson(counts: list[int]) -> tuple[float, float, float]:
    """Simpson's index of diversity (Hunter & Gaston 1988) with Grundmann et al. (2001) 95% CI."""
    n = sum(counts)
    if n < 2:
        return math.nan, math.nan, math.nan
    d = 1 - sum(c * (c - 1) for c in counts) / (n * (n - 1))
    p = [c / n for c in counts]
    variance = 4 / n * (sum(x**3 for x in p) - sum(x**2 for x in p) ** 2)
    half = 2 * math.sqrt(max(variance, 0.0))
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
    dataset: Dataset, min_isolates: int = 5, permutations: int = 9999, seed: int = 1
) -> list[dict]:
    """Thorpe et al. 2022: is a lineage over- or under-represented in a source?

    Source labels are permuted across typed isolates; two-sided p-values use
    |count - expected|, are never zero (Phipson & Smyth 2010) and are BH-adjusted.
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
    expected = np.outer([sizes[t] for t in tested], np.bincount(src, minlength=n_s)) / len(typed)
    deviation = np.abs(observed - expected) - 1e-9
    exceed = np.zeros_like(observed)
    rng = np.random.default_rng(seed)
    for _ in range(permutations):
        exceed += np.abs(counts(rng.permutation(src)) - expected) >= deviation
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
) -> FeatureResult:
    """Prevalence (Wilson 95% CI) of each feature per source, and a permutation
    chi-squared test of independence from source with BH correction."""
    mapping, unmatched = match_names(list(presence), dataset)
    by_key: dict[str, set[str]] = {}
    for name, key in mapping.items():
        by_key.setdefault(key, set()).update(presence[name])
    isolates = [iso for iso in dataset.classified if iso.key in by_key]
    groups = [g for g in dataset.groups if any(iso.source == g for iso in isolates)]
    result = FeatureResult(groups=groups, genotyped={g: 0 for g in groups}, unmatched=unmatched)
    if not isolates:
        return result
    features = sorted({f for iso in isolates for f in by_key[iso.key]})
    src_index = {g: i for i, g in enumerate(groups)}
    src = np.array([src_index[iso.source] for iso in isolates])
    matrix = np.array([[f in by_key[iso.key] for f in features] for iso in isolates], dtype=np.int64)
    n_s = np.bincount(src, minlength=len(groups))
    result.genotyped = {g: int(n_s[i]) for i, g in enumerate(groups)}
    onehot = np.eye(len(groups), dtype=np.int64)
    observed = onehot[src].T @ matrix  # sources x features

    totals = matrix.sum(axis=0)
    testable = (totals > 0) & (totals < len(isolates)) & (len(groups) > 1)
    stat_obs = _chi2(observed, n_s, totals)
    exceed = np.zeros(len(features), dtype=np.int64)
    if testable.any():
        rng = np.random.default_rng(seed)
        sub = matrix[:, testable]
        for _ in range(permutations):
            perm = onehot[rng.permutation(src)].T @ sub
            exceed[testable] += _chi2(perm, n_s, totals[testable]) >= stat_obs[testable] - 1e-9
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
    """Pearson chi-squared for each feature's 2 x sources table."""
    n = n_s.sum()
    expected = np.outer(n_s, totals) / n
    absent = n_s[:, None] - present
    expected_absent = np.outer(n_s, n - totals) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = (present - expected) ** 2 / expected + (absent - expected_absent) ** 2 / expected_absent
    return np.nansum(np.where(np.isfinite(terms), terms, 0.0), axis=0)


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
    distances_within: np.ndarray
    distances_between: np.ndarray


def linkage(
    dataset: Dataset,
    names: list[str],
    matrix: np.ndarray,
    threshold: float = 21,
    permutations: int = 9999,
    seed: int = 1,
) -> LinkageResult:
    """Strain sharing within and between sources at a SNP threshold.

    Following Thorpe et al. 2022 and Hetland et al. 2025: isolates within the
    threshold are joined by single linkage; each unique source pair connected by
    at least one qualifying pair inside a cluster is one sharing event. Pairs from
    the same host individual are excluded. Events are normalised per 100 isolates
    of the two sources, and compared with a null in which source labels are
    permuted across isolates (sampling effort held fixed; cf. Muloi et al. 2022).
    """
    isolates, sub, unmatched = _matched_isolates(dataset, names, matrix)
    groups = [g for g in dataset.groups if any(iso.source == g for iso in isolates)]
    n_s = len(groups)
    src = np.array([groups.index(iso.source) if iso.source in groups else -1 for iso in isolates], dtype=np.int64)
    hosts = np.array([iso.host_id for iso in isolates], dtype=object)
    per_source = {g: int((src == k).sum()) for k, g in enumerate(groups)}

    iu, ju = np.triu_indices(len(isolates), k=1)
    values = sub[iu, ju] if len(isolates) else np.zeros(0)
    finite = np.isfinite(values)
    classified = (src[iu] >= 0) & (src[ju] >= 0) & finite
    same = classified & (src[iu] == src[ju])
    within = values[same]
    between = values[classified & ~same]

    def edges_at(t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        hit = finite & (values <= t)
        a, b = iu[hit], ju[hit]
        clusters = _components(len(isolates), a, b)
        same_host = np.array([bool(hosts[x]) and hosts[x] == hosts[y] for x, y in zip(a, b)], dtype=bool)
        keep = ~same_host
        return a[keep], b[keep], clusters

    def event_counts(labels: np.ndarray, a: np.ndarray, b: np.ndarray, clusters: np.ndarray) -> np.ndarray:
        la, lb = labels[a], labels[b]
        ok = (la >= 0) & (lb >= 0)
        lo, hi = np.minimum(la[ok], lb[ok]), np.maximum(la[ok], lb[ok])
        codes = np.unique(clusters[a[ok]] * n_s * n_s + lo * n_s + hi)
        return np.bincount(codes % (n_s * n_s), minlength=n_s * n_s).reshape(n_s, n_s)

    a, b, clusters = edges_at(threshold)
    observed = event_counts(src, a, b, clusters) if n_s else np.zeros((0, 0), dtype=np.int64)

    # Permutation null: shuffle labels among classified isolates, graph fixed.
    exceed = np.zeros_like(observed)
    perm_values = []
    labelled = np.flatnonzero(src >= 0)
    if n_s and len(a):
        rng = np.random.default_rng(seed)
        draws = []
        for _ in range(permutations):
            labels = src.copy()
            labels[labelled] = rng.permutation(src[labelled])
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
        ta, tb, tc = edges_at(t)
        counts = event_counts(src, ta, tb, tc) if n_s else np.zeros((0, 0))
        diagonal = int(np.trace(counts)) if n_s else 0
        sweep.append({"threshold": t, "within_events": diagonal, "between_events": int(counts.sum()) - diagonal})

    cross_edges = []
    for x, y in zip(a, b):
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
                    "snps": float(sub[x, y]),
                    "years_apart": gap,
                }
            )

    cluster_rows = []
    for cid in sorted(set(clusters.tolist())):
        members = np.flatnonzero(clusters == cid)
        if len(members) < 2:
            continue
        sources = Counter(isolates[i].source or "Unclassified" for i in members)
        inside = [(int(x), int(y)) for x, y in zip(a, b) if clusters[x] == cid]
        cluster_rows.append(
            {
                "cluster": int(cid),
                "isolates": len(members),
                "sources": ", ".join(f"{g} ({sources[g]})" for g in [*groups, "Unclassified"] if sources[g]),
                "n_sources": len([g for g in groups if sources[g]]),
                "lineages": ", ".join(sorted({isolates[i].lineage for i in members if isolates[i].lineage})),
                "members": [isolates[i].key for i in members],
                "member_sources": [isolates[i].source for i in members],
                "edges": [(int(np.searchsorted(members, x)), int(np.searchsorted(members, y)), float(sub[x, y])) for x, y in inside],
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
        distances_within=within,
        distances_between=between,
    )


def _matched_isolates(
    dataset: Dataset, names: list[str], matrix: np.ndarray
) -> tuple[list[Isolate], np.ndarray, list[str]]:
    """Isolates found in a distance matrix (first occurrence of each), their sub-matrix, and unmatched names."""
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
    sub = matrix[np.ix_(index, index)] if index else np.zeros((0, 0))
    return isolates, sub, unmatched


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
    isolates, sub, _ = _matched_isolates(dataset, names, matrix)
    iu, ju = np.triu_indices(len(isolates), k=1)
    values = sub[iu, ju] if len(isolates) else np.zeros(0)
    sources = np.array([iso.source for iso in isolates], dtype=object)
    link = np.isfinite(values) & (values <= threshold) & (sources[iu] == sources[ju])
    clusters = _components(len(isolates), iu[link], ju[link])
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
