"""The statistics in `report`, checked against independent implementations.

Three kinds of check, so that no number in a report rests on onehealth's code alone:

1. R and vegan. tests/data/stats_reference/reference.tsv holds values computed by
   reference.R on the inputs in cases.tsv (counts from the real Klebsiella fixture, and
   edge cases). The same inputs go through onehealth here.
2. scipy, on the tables `report` makes from the real Klebsiella data.
3. Every possible relabelling. A permutation test samples the ways of shuffling source
   labels; on a dataset of ten isolates all of them can be listed, which gives the exact
   p-value. The statistic is recomputed here in a few plain lines, not with onehealth's
   own code.

documentation/guides/statistics.md describes the result.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse, stats
from scipy.sparse.csgraph import connected_components

from onehealth import analysis as A
from onehealth.classify import classify_rows
from onehealth.io_tsv import read_distance_matrix, read_features, read_tsv

REFERENCE = Path(__file__).parent / "data" / "stats_reference"
KLEBSIELLA = Path(__file__).parent / "data" / "ncbi_pd_klebsiella"
SOURCES = ("Human", "Animal", "Environmental", "Food", "Healthcare environment")
PERMUTATIONS = 9999
# Deterministic quantities must agree with R to this many decimal places.
EXACT = 1e-12


def _table(name: str) -> list[list[str]]:
    lines = [line.rstrip("\n") for line in (REFERENCE / name).read_text(encoding="utf-8").splitlines()]
    return [line.split("\t") for line in lines if line and not line.startswith("#")][1:]


CASES = {(kind, case): values for kind, case, values in _table("cases.tsv")}
R = {(kind, case, quantity): (math.nan if value == "NA" else float(value)) for kind, case, quantity, value in _table("reference.tsv")}


def _cases(kind: str) -> list[str]:
    return [case for k, case in CASES if k == kind]


def _numbers(text: str) -> list[float]:
    return [math.nan if part == "NA" else float(part) for part in text.split(",")]


def _named(text: str) -> dict[str, list[int]]:
    parts = (part.rpartition("=") for part in text.split(";"))
    return {name: [int(x) for x in values.split(",")] for name, _, values in parts}


def _row(i: int, source: str, lineage: str = "", year: int = 2020, host_id: str = "") -> dict[str, str]:
    healthcare = source == "Healthcare environment"
    return {
        "query_id": f"ISO{i}",
        "biosample": f"SAMN{i}",
        "l1": "Environmental" if healthcare else source,
        "healthcare_environment": "yes" if healthcare else "no",
        "st": lineage,
        "collection_date": str(year),
        "host_subject_id": host_id,
    }


def _monte_carlo_margin(p: float, permutations: int) -> float:
    """How far a p-value from this many random permutations may be from the exact one: 4.5 standard errors."""
    p = min(max(p, 0.0), 1.0)  # a sum of probabilities can come out as 1.0000000000000002
    return 4.5 * math.sqrt(p * (1 - p) / permutations) + 2 / permutations


# --- 1. against R and vegan ----------------------------------------------------------------


def test_the_reference_file_covers_every_case():
    assert len(CASES) == 88 and len(R) == 567
    assert {kind for kind, _ in CASES} == {"counts", "proportion", "pvalues", "features", "lineages"}


@pytest.mark.parametrize("case", _cases("counts"))
def test_simpson_index_matches_vegan(case):
    counts = [int(x) for x in _numbers(CASES["counts", case])]
    ours, reference = A.simpson(counts)[0], R["counts", case, "simpson"]
    if math.isnan(reference):  # one isolate: no pair to compare
        assert math.isnan(ours)
    else:
        assert ours == pytest.approx(reference, abs=EXACT)


@pytest.mark.parametrize("case", [c for c in _cases("counts") if ("counts", c, "simpson_variance") in R])
def test_simpson_variance_matches_two_million_samples_drawn_in_r(case):
    """The variance behind the interval is the variance of the index over repeated samples."""
    counts = [int(x) for x in _numbers(CASES["counts", case])]
    simulated, error = R["counts", case, "simpson_variance"], R["counts", case, "simpson_variance_se"]
    assert A.simpson_variance(counts) == pytest.approx(simulated, abs=5 * error + 1e-15)


@pytest.mark.parametrize("case", _cases("counts"))
def test_rarefaction_matches_vegan(case):
    counts = [int(x) for x in _numbers(CASES["counts", case])]
    sizes = sorted(int(q.split("_")[1]) for k, c, q in R if (k, c) == ("counts", case) and q.startswith("rarefy_"))
    assert sizes or sum(counts) < 2
    for size, ours in zip(sizes, A.rarefaction(counts, sizes)):
        assert ours == pytest.approx(R["counts", case, f"rarefy_{size}"], abs=1e-9)


@pytest.mark.parametrize("case", _cases("proportion"))
def test_wilson_interval_matches_prop_test(case):
    successes, n = (int(x) for x in _numbers(CASES["proportion", case]))
    low, high = A.wilson(successes, n)
    assert low == pytest.approx(R["proportion", case, "wilson_low"], abs=EXACT)
    assert high == pytest.approx(R["proportion", case, "wilson_high"], abs=EXACT)


@pytest.mark.parametrize("case", _cases("pvalues"))
def test_benjamini_hochberg_matches_p_adjust(case):
    ours = A.benjamini_hochberg(_numbers(CASES["pvalues", case]))
    for k, q in enumerate(ours, start=1):
        reference = R["pvalues", case, f"q_{k}"]
        assert math.isnan(q) if math.isnan(reference) else q == pytest.approx(reference, abs=EXACT)


def _feature_dataset(case: str):
    """A dataset with the per-source counts of one `features` case: each feature is in the first isolates of a source."""
    parts = _named(CASES["features", case])
    sizes = parts.pop("sizes")
    rows, presence, i = [], {}, 0
    for source, size in zip(SOURCES, sizes):
        for k in range(size):
            rows.append(_row(i, source))
            presence[f"ISO{i}"] = {feature for feature, counts in parts.items() if k < counts[SOURCES.index(source)]}
            i += 1
    return A.build_dataset(rows), presence, list(parts)


@pytest.mark.parametrize("case", _cases("features"))
def test_feature_test_matches_chisq_test(case):
    """The statistic matches R exactly; the permutation p-value matches R's from a million random tables."""
    dataset, presence, features = _feature_dataset(case)
    result = A.feature_prevalence(dataset, presence, {}, permutations=PERMUTATIONS, seed=1)
    tests = {t["feature"]: t for t in result.tests}
    assert set(tests) == set(features)
    for feature in features:
        assert tests[feature]["chi2"] == pytest.approx(R["features", f"{case}|{feature}", "chi2"], rel=1e-10)
        reference = R["features", f"{case}|{feature}", "p_simulated"]
        assert abs(tests[feature]["p"] - reference) <= _monte_carlo_margin(reference, PERMUTATIONS), feature


@pytest.mark.parametrize("case", _cases("lineages"))
def test_lineage_enrichment_matches_the_hypergeometric_distribution(case):
    """With no strata, the isolates of a lineage that fall in a source follow a hypergeometric
    distribution, so the p-value the permutations estimate is known exactly."""
    parts = _named(CASES["lineages", case])
    rows = []
    for lineage, counts in parts.items():
        for source, n in zip(SOURCES, counts):
            rows += [_row(len(rows) + k, source, lineage) for k in range(n)]
    result = A.lineage_enrichment(A.build_dataset(rows), min_isolates=5, permutations=PERMUTATIONS, seed=1)
    compared = 0
    for row in result:
        key = f"{case}|{row['lineage']}|{SOURCES.index(row['source']) + 1}"
        assert row["expected"] == pytest.approx(R["lineages", key, "expected"], abs=1e-9)
        exact = R["lineages", key, "p_exact"]
        assert abs(row["p"] - exact) <= _monte_carlo_margin(exact, PERMUTATIONS), key
        compared += 1
    assert compared == sum(1 for k, c, q in R if k == "lineages" and c.startswith(case + "|") and q == "p_exact")


# --- 2. against scipy, on the tables made from real data -----------------------------------


@pytest.fixture(scope="module")
def klebsiella(lexicon):
    classified = classify_rows(read_tsv(KLEBSIELLA / "metadata.tsv"), lexicon)
    return A.build_dataset(classified, lineage_column="PDS_acc")


def test_real_diversity_table_matches_scipy(klebsiella):
    """Rarefaction is one hypergeometric probability per lineage: the chance that a subsample misses it."""
    curves = A.rarefaction_curves(klebsiella)
    assert len(curves) > 20
    typed = [iso for iso in klebsiella.classified if iso.lineage]
    for point in curves:
        counts = list(Counter(iso.lineage for iso in typed if iso.source == point["source"]).values())
        missed = stats.hypergeom.pmf(0, sum(counts), counts, point["isolates_sampled"])
        assert point["expected_lineages"] == pytest.approx(float((1 - missed).sum()), abs=1e-9)
    for row in A.diversity(klebsiella):
        counts = np.array(list(Counter(iso.lineage for iso in typed if iso.source == row["source"]).values()))
        if row["isolates"] >= 2:
            # Two isolates drawn without replacement are of different lineages with this probability.
            different = 1 - float(stats.hypergeom.pmf(2, row["isolates"], counts, 2).sum())
            assert row["simpson"] == pytest.approx(different, abs=EXACT)


def test_real_feature_table_matches_scipy(klebsiella):
    presence, classes = read_features(KLEBSIELLA / "features.tsv")
    result = A.feature_prevalence(klebsiella, presence, classes, permutations=199, seed=1)
    assert len(result.rows) > 500
    for row in result.rows:
        interval = stats.binomtest(row["present"], row["genotyped"]).proportion_ci(confidence_level=0.95, method="wilson")
        assert row["ci_low"] == pytest.approx(100 * interval.low, abs=1e-9)
        assert row["ci_high"] == pytest.approx(100 * interval.high, abs=1e-9)
    sizes = np.array([result.genotyped[g] for g in result.groups])
    by_feature: dict[str, list[int]] = {}
    for row in result.rows:
        by_feature.setdefault(row["feature"], []).append(row["present"])
    tested = [t for t in result.tests if not math.isnan(t["p"])]
    assert len(tested) > 100
    for test in tested:
        present = np.array(by_feature[test["feature"]])
        reference = stats.chi2_contingency(np.vstack([present, sizes - present]), correction=False).statistic
        assert test["chi2"] == pytest.approx(reference, rel=1e-10)
    q = stats.false_discovery_control([t["p"] for t in tested], method="bh")
    assert [t["q"] for t in tested] == pytest.approx(list(q), abs=EXACT)


@pytest.mark.parametrize("threshold", [0, 5, 21, 50])
def test_real_clusters_match_scipy_connected_components(klebsiella, threshold):
    """Single linkage at a threshold is the connected components of the graph of pairs within it."""
    names, matrix = read_distance_matrix(KLEBSIELLA / "snp_dists.tsv")
    result = A.linkage(klebsiella, names, matrix, threshold=threshold, permutations=0)
    ours = {frozenset(cluster["members"]) for cluster in result.clusters}

    mapping, _ = A.match_names(names, klebsiella)
    first: dict[str, int] = {}
    for i, name in enumerate(names):
        first.setdefault(mapping[name], i)  # two assemblies of one BioSample: the first is the isolate
    keys, index = list(first), np.array(list(first.values()))
    sub = matrix[np.ix_(index, index)]
    graph = sparse.csr_matrix(np.isfinite(sub) & (sub <= threshold) & ~np.eye(len(index), dtype=bool))
    count, labels = connected_components(graph, directed=False)
    members: dict[int, set[str]] = {}
    for key, label in zip(keys, labels):
        members.setdefault(int(label), set()).add(key)
    reference = {frozenset(m) for m in members.values() if len(m) > 1}
    assert ours == reference and len(reference) > 5


# --- 3. against every possible relabelling of a small dataset ------------------------------

# Ten isolates. Year is the stratum. Lineage L1 is mostly human, L2 mostly animal.
SMALL = [
    # source, lineage, year, patient
    ("Human", "L1", 2019, "p1"),
    ("Human", "L1", 2019, "p1"),
    ("Animal", "L2", 2019, ""),
    ("Environmental", "L1", 2019, ""),
    ("Environmental", "L2", 2019, ""),
    ("Human", "L1", 2020, ""),
    ("Human", "L1", 2020, ""),
    ("Animal", "L2", 2020, ""),
    ("Animal", "L2", 2020, ""),
    ("Environmental", "L2", 2020, ""),
]
# SNP distances: isolates 0-3, 4-6 and 7-9 are three chains of close relatives.
SMALL_CLOSE = [(0, 1, 2), (1, 2, 5), (2, 3, 8), (4, 5, 3), (5, 6, 9), (7, 8, 1), (8, 9, 4), (0, 3, 30), (4, 9, 40)]
SMALL_GENES = {"g1": [0, 1, 3, 5, 6], "g2": [2, 7, 8, 9], "g3": [0, 2, 4, 6, 8], "g4": [1]}
THRESHOLD = 10


def _relabellings(sources: list[str], strata: list[int] | None):
    """Every distinct way of handing the source labels back to the isolates (within strata, if given)."""
    if strata is None:
        yield from sorted(set(itertools.permutations(sources)))
        return
    blocks = sorted(set(strata))
    options = [sorted(set(itertools.permutations([s for s, k in zip(sources, strata) if k == block]))) for block in blocks]
    slots = [[i for i, k in enumerate(strata) if k == block] for block in blocks]
    for choice in itertools.product(*options):
        labels = [""] * len(sources)
        for positions, block_labels in zip(slots, choice):
            for position, label in zip(positions, block_labels):
                labels[position] = label
        yield tuple(labels)


def _exact_p(observed: float, values: list[float], centre: float | None = None) -> float:
    """Share of relabellings at least as extreme as the data: one-sided, or two-sided about `centre`."""
    if centre is None:
        return sum(v >= observed - 1e-9 for v in values) / len(values)
    return sum(abs(v - centre) >= abs(observed - centre) - 1e-9 for v in values) / len(values)


def _small_dataset():
    rows = [_row(i, source, lineage, year, patient) for i, (source, lineage, year, patient) in enumerate(SMALL)]
    return A.build_dataset(rows)


@pytest.mark.parametrize("strata", [(), ("year",)])
def test_enrichment_p_values_match_every_relabelling(strata):
    sources, lineages = [s for s, *_ in SMALL], [lin for _, lin, *_ in SMALL]
    years = [year for *_, year, _ in SMALL] if strata else None
    relabelled = list(_relabellings(sources, years))
    assert len(relabelled) == (4200 if not strata else 900)  # 10!/(4! 3! 3!), or 30 x 30 within the two years
    result = A.lineage_enrichment(_small_dataset(), min_isolates=5, permutations=PERMUTATIONS, seed=1, strata=strata)
    assert len(result) == 6
    for row in result:
        def count(labels):
            return sum(1 for label, lineage in zip(labels, lineages) if label == row["source"] and lineage == row["lineage"])

        values = [count(labels) for labels in relabelled]
        expected = sum(values) / len(values)
        assert row["observed"] == count(sources)
        assert row["expected"] == pytest.approx(expected, abs=1e-9)
        exact = _exact_p(row["observed"], values, centre=expected)
        assert abs(row["p"] - exact) <= _monte_carlo_margin(exact, PERMUTATIONS), row


@pytest.mark.parametrize("strata", [(), ("year",)])
def test_feature_p_values_match_every_relabelling(strata):
    sources = [s for s, *_ in SMALL]
    years = [year for *_, year, _ in SMALL] if strata else None
    presence = {f"ISO{i}": {gene for gene, carriers in SMALL_GENES.items() if i in carriers} for i in range(len(SMALL))}
    result = A.feature_prevalence(_small_dataset(), presence, {}, permutations=PERMUTATIONS, seed=1, strata=strata)
    relabelled = list(_relabellings(sources, years))
    groups = sorted(set(sources))

    def chi2(labels, carriers):
        table = np.array([[sum(1 for i, label in enumerate(labels) if label == g and (i in carriers) == present) for g in groups]
                          for present in (True, False)])
        return stats.chi2_contingency(table, correction=False).statistic

    assert {t["feature"] for t in result.tests} == set(SMALL_GENES)
    for test in result.tests:
        carriers = set(SMALL_GENES[test["feature"]])
        assert test["chi2"] == pytest.approx(chi2(sources, carriers), rel=1e-10)
        exact = _exact_p(test["chi2"], [chi2(labels, carriers) for labels in relabelled])
        assert abs(test["p"] - exact) <= _monte_carlo_margin(exact, PERMUTATIONS), test


def _sharing_events(labels, patients) -> Counter:
    """Sharing events, from the definition: in each cluster of isolates joined by pairs within the
    threshold, every pair of sources that has such a pair (not two isolates of one patient) is one event."""
    n = len(labels)
    close = [(a, b) for a, b, snps in SMALL_CLOSE if snps <= THRESHOLD]
    graph = sparse.csr_matrix((np.ones(len(close)), ([a for a, _ in close], [b for _, b in close])), shape=(n, n))
    _, cluster = connected_components(graph, directed=False)
    seen = {(int(cluster[a]), *sorted((labels[a], labels[b]))) for a, b in close if not (patients[a] and patients[a] == patients[b])}
    return Counter((x, y) for _, x, y in seen)


@pytest.mark.parametrize("strata", [(), ("year",)])
def test_linkage_events_and_p_values_match_every_relabelling(strata):
    sources, patients = [s for s, *_ in SMALL], [patient for *_, patient in SMALL]
    years = [year for *_, year, _ in SMALL] if strata else None
    n = len(SMALL)
    matrix = np.full((n, n), 500.0)
    np.fill_diagonal(matrix, 0.0)
    for a, b, snps in SMALL_CLOSE:
        matrix[a, b] = matrix[b, a] = snps
    names = [f"ISO{i}" for i in range(n)]
    result = A.linkage(_small_dataset(), names, matrix, threshold=THRESHOLD, permutations=PERMUTATIONS, seed=1, strata=strata)

    observed = _sharing_events(sources, patients)
    relabelled = [_sharing_events(labels, patients) for labels in _relabellings(sources, years)]
    assert sum(observed.values()) > 3
    for row in result.events:
        pair = tuple(sorted((row["source_a"], row["source_b"])))
        assert row["events"] == observed[pair], pair
        values = [events[pair] for events in relabelled]
        mean = sum(values) / len(values)
        assert row["expected"] == pytest.approx(mean, abs=0.05)
        exact = _exact_p(observed[pair], values, centre=mean)
        assert abs(row["p"] - exact) <= _monte_carlo_margin(exact, PERMUTATIONS), row


def test_shuffling_within_strata_gives_every_relabelling_equally_often():
    """The permutation itself: each of the 900 relabellings within year should come up equally often."""
    sources = np.array([s for s, *_ in SMALL], dtype=object)
    years = np.array([year for *_, year, _ in SMALL])
    rng = np.random.default_rng(1)
    draws = 90_000
    seen = Counter(tuple(A.permute_within(sources, years, rng)) for _ in range(draws))
    assert set(seen) == set(_relabellings(list(sources), list(years)))
    assert stats.chisquare(list(seen.values())).pvalue > 0.001
    plain = Counter(tuple(A.permute_within(sources, None, rng)) for _ in range(draws))
    assert len(plain) == 4200 and stats.chisquare(list(plain.values())).pvalue > 0.001


# --- the documentation ------------------------------------------------------------------------


def test_the_values_quoted_in_the_documentation_are_the_ones_computed():
    """documentation/guides/statistics.md shows a few values side by side; each is recomputed here."""
    page = (Path(__file__).parents[1] / "documentation" / "guides" / "statistics.md").read_text(encoding="utf-8")
    human = [int(x) for x in _numbers(CASES["counts", "klebsiella_human"])]
    assert (sum(human), len(human)) == (231, 22)
    low, high = A.wilson(151, 231)
    dataset, presence, _ = _feature_dataset("klebsiella")
    gene = next(t for t in A.feature_prevalence(dataset, presence, {}, permutations=PERMUTATIONS, seed=1).tests
                if t["feature"] == "aac(6')-Ib-cr5")
    rows = []
    for lineage, counts in _named(CASES["lineages", "klebsiella"]).items():
        for source, n in zip(SOURCES, counts):
            rows += [_row(len(rows) + k, source, lineage) for k in range(n)]
    cell = next(r for r in A.lineage_enrichment(A.build_dataset(rows), min_isolates=5, permutations=PERMUTATIONS, seed=1)
                if (r["lineage"], r["source"]) == ("PDS000045679.7", "Human"))
    assert cell["observed"] == 27 and f"{cell['expected']:.1f}" == "24.9"
    for line in (
        f"| Simpson's index, 231 human isolates in 22 lineages | {A.simpson(human)[0]:.6f} | {R['counts', 'klebsiella_human', 'simpson']:.6f} |",
        f"| Expected lineages in 58 of those isolates | {A.rarefaction(human, [58])[0]:.4f} | {R['counts', 'klebsiella_human', 'rarefy_58']:.4f} |",
        f"| Wilson 95% interval, `oqxB` in 151 of 231 human isolates | {100 * low:.4f}–{100 * high:.4f}% | "
        f"{100 * R['proportion', 'oqxB|Human', 'wilson_low']:.4f}–{100 * R['proportion', 'oqxB|Human', 'wilson_high']:.4f}% |",
        f"| Chi-squared, `aac(6')-Ib-cr5` by source | {gene['chi2']:.4f} | {R['features', 'klebsiella|aac(6' + chr(39) + ')-Ib-cr5', 'chi2']:.4f} |",
        f"| Its permutation p-value | {gene['p']:.4f} | {R['features', 'klebsiella|aac(6' + chr(39) + ')-Ib-cr5', 'p_simulated']:.4f} |",
        f"| Enrichment of lineage `PDS000045679.7` in Human (27 isolates, 24.9 expected), p | {cell['p']:.4f} | "
        f"{R['lineages', 'klebsiella|PDS000045679.7|1', 'p_exact']:.4f} (exact) |",
    ):
        assert line in page, line
    old, new = A.wilson(51, 150, 1.96), A.wilson(51, 150)
    assert f"{100 * old[0]:.4f}–{100 * old[1]:.4f}% became {100 * new[0]:.4f}–{100 * new[1]:.4f}%" in page
    assert f"{len(CASES)} inputs" in page and f"computes {len(R)} values" in page
