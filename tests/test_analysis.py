import math
from pathlib import Path

import numpy as np
import pytest

from onehealth import analysis as A
from onehealth.io_tsv import read_distance_matrix, read_features, read_lineages


def _row(i, l1="Human", st="ST1", hc="no", host_id="", date="2019-05-01", biosample=None):
    return {
        "query_id": f"ISO{i}",
        "biosample": biosample or f"SAMN{i}",
        "l1": l1,
        "healthcare_environment": hc,
        "paper_rollup": "Environmental" if hc == "yes" or l1 in {"Food", "Laboratory", "Environmental"} else l1,
        "st": st,
        "host_subject_id": host_id,
        "collection_date": date,
    }


# --- statistics checked against hand-computed values -------------------------------------


def test_simpson_matches_hunter_gaston():
    assert A.simpson([1, 1])[0] == 1.0
    assert A.simpson([3, 1])[0] == pytest.approx(1 - 6 / 12)
    assert math.isnan(A.simpson([1])[0])
    d, lo, hi = A.simpson([5, 5, 5, 5])
    assert lo <= d <= hi


def test_rarefaction_is_exact_for_simple_cases():
    assert A.rarefaction([1, 1, 1, 1], [1, 2, 4]) == pytest.approx([1, 2, 4])
    assert A.rarefaction([4], [1, 3]) == pytest.approx([1, 1])
    # two lineages of 2: E[S_2] = 2 - 2 * C(2,2)/C(4,2) = 2 - 2/6
    assert A.rarefaction([2, 2], [2]) == pytest.approx([2 - 2 / 6])


def test_wilson_interval():
    lo, hi = A.wilson(0, 10)
    assert lo == 0 and hi == pytest.approx(0.2775, abs=1e-3)
    lo, hi = A.wilson(5, 10)
    assert (lo, hi) == pytest.approx((0.2366, 0.7634), abs=1e-3)


def test_benjamini_hochberg():
    assert A.benjamini_hochberg([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.04, 0.04])
    q = A.benjamini_hochberg([0.5, math.nan])
    assert q[0] == 0.5 and math.isnan(q[1])


def test_permutation_p_is_never_zero():
    assert A.perm_pvalue(0, 999) == pytest.approx(0.001)


def test_decimal_year():
    assert A.decimal_year("2019") == 2019.5
    assert A.decimal_year("2019-07") == pytest.approx(2019.54, abs=0.01)
    assert 2019.0 < A.decimal_year("2019-01-02") < 2019.01
    assert A.decimal_year("missing") is None
    assert A.decimal_year("2019-13-40") == 2019.5


# --- dataset --------------------------------------------------------------------------------


def test_dataset_collapses_rows_of_the_same_biosample_and_splits_healthcare():
    rows = [_row(1), _row(2, biosample="SAMN1"), _row(3, "Environmental", hc="yes"), _row(4, "Unclassified")]
    ds = A.build_dataset(rows)
    assert (len(ds.isolates), ds.duplicates, ds.unclassified) == (3, 1, 1)
    assert {iso.source for iso in ds.classified} == {"Human", "Healthcare environment"}
    rollup = A.build_dataset(rows, by="paper_rollup")
    assert {iso.source for iso in rollup.classified} == {"Human", "Environmental"}


def test_lineage_sharing_summary_and_intersections():
    rows = [_row(1, "Human", "ST1"), _row(2, "Animal", "ST1"), _row(3, "Human", "ST2"), _row(4, "Food", "ST3")]
    ds = A.build_dataset(rows)
    summary = A.lineage_sharing_summary(ds)
    assert summary == {"lineages": 3, "shared_lineages": 1, "typed_isolates": 4, "isolates_in_shared": 2}
    combos = {r["sources"]: r["lineages"] for r in A.lineage_intersections(ds)}
    assert combos == {"Human": 1, "Food": 1, "Human + Animal": 1}


def test_lineage_enrichment_finds_a_source_restricted_lineage_and_is_reproducible():
    rows = [_row(i, "Animal", "ST661") for i in range(12)]
    rows += [_row(100 + i, "Human", f"ST{i % 4}") for i in range(30)]
    ds = A.build_dataset(rows)
    first = A.lineage_enrichment(ds, permutations=999, seed=3)
    again = A.lineage_enrichment(ds, permutations=999, seed=3)
    assert first == again
    hit = next(r for r in first if r["lineage"] == "ST661" and r["source"] == "Animal")
    assert hit["direction"] == "enriched" and hit["q"] < 0.01


# --- features ------------------------------------------------------------------------------


def test_feature_prevalence_counts_tests_and_unmatched(tmp_path: Path):
    rows = [_row(i, "Animal") for i in range(15)] + [_row(100 + i, "Human") for i in range(15)]
    path = tmp_path / "amr.tsv"
    lines = ["Name\tElement symbol\tClass"]
    lines += [f"ISO{i}\ttet(A)\tTETRACYCLINE" for i in range(15)]
    lines += [f"ISO{100 + i}\tblaSHV-1\tBETA-LACTAM" for i in range(15)]
    lines += ["GHOST\tmcr-1\tCOLISTIN"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    presence, classes = read_features(path)
    result = A.feature_prevalence(A.build_dataset(rows), presence, classes, permutations=999)
    assert result.genotyped == {"Human": 15, "Animal": 15}
    assert result.unmatched == ["GHOST"]
    tet = next(r for r in result.rows if r["feature"] == "tet(A)" and r["source"] == "Animal")
    assert (tet["present"], tet["genotyped"], tet["percent"]) == (15, 15, 100.0)
    test = next(t for t in result.tests if t["feature"] == "tet(A)")
    assert test["class"] == "TETRACYCLINE" and test["q"] < 0.01


# --- linkage -------------------------------------------------------------------------------


def _linkage_fixture():
    rows = [
        _row(0, "Human", host_id="P1"),
        _row(1, "Human", host_id="P1"),  # same patient as 0
        _row(2, "Animal"),
        _row(3, "Human", host_id="P2"),
        _row(4, "Food", st="ST9"),
        _row(5, "Animal", st="ST9"),
    ]
    names = [f"ISO{i}" for i in range(6)]
    d = np.full((6, 6), 5000.0)
    np.fill_diagonal(d, 0)
    for i, j, v in [(0, 1, 1), (0, 2, 3), (1, 2, 4), (2, 3, 2), (4, 5, 10)]:
        d[i, j] = d[j, i] = v
    return A.build_dataset(rows), names, d


def test_linkage_counts_each_source_pair_once_per_cluster_and_skips_same_host():
    ds, names, d = _linkage_fixture()
    link = A.linkage(ds, names, d, threshold=5, permutations=99)
    events = {(r["source_a"], r["source_b"]): r["events"] for r in link.events}
    # Cluster {0,1,2,3}: Human-Animal (several pairs) counts once; Human-Human only via 0-1 (same host) -> 0.
    assert events[("Human", "Animal")] == 1
    assert events[("Human", "Human")] == 0
    # Cluster {4,5} needs 10 SNPs, so Animal-Food is 0 at 5 SNPs.
    assert events[("Animal", "Food")] == 0
    sweep = {r["threshold"]: r for r in link.sweep}
    assert sweep[10]["between_events"] == 2
    linked = {(r["from_source"], r["to_source"]): r["with_close_relative"] for r in link.linked_fraction}
    assert linked[("Human", "Animal")] == 3 and linked[("Animal", "Human")] == 1


def test_linkage_reports_unmatched_names_and_cross_edges():
    ds, names, d = _linkage_fixture()
    link = A.linkage(ds, names + ["NOPE"], np.pad(d, ((0, 1), (0, 1)), constant_values=9999), threshold=5, permutations=9)
    assert link.unmatched == ["NOPE"] and link.matched == 6
    assert {(e["source_a"], e["source_b"]) for e in link.cross_edges} <= {("Human", "Animal"), ("Animal", "Human")}


def test_match_names_accepts_ncbi_assembly_file_names():
    ds = A.build_dataset([{**_row(1), "assembly": "GCA_000005845.2"}])
    mapping, unmatched = A.match_names(["GCA_000005845.2_ASM584v2_genomic", "other"], ds)
    assert mapping == {"GCA_000005845.2_ASM584v2_genomic": "SAMN1"} and unmatched == ["other"]


# --- readers -------------------------------------------------------------------------------


def test_read_distance_matrix_square_and_long(tmp_path: Path):
    square = tmp_path / "sq.tsv"
    square.write_text("snp-dists 0.8.2\tA\tB\nA\t0\t7\nB\t7\t0\n", encoding="utf-8")
    names, m = read_distance_matrix(square)
    assert names == ["A", "B"] and m[0, 1] == 7
    long = tmp_path / "long.tsv"
    long.write_text("A\tB\t7\nA\tC\t3\n", encoding="utf-8")
    names, m = read_distance_matrix(long)
    assert names == ["A", "B", "C"] and m[0, 2] == 3 and math.isnan(m[1, 2])


def test_read_distance_matrix_csv_comments_crlf_and_missing_values(tmp_path: Path):
    path = tmp_path / "d.csv"
    path.write_bytes(b"# snp-dists -c\r\n\r\nsnp-dists, A , B\r\n A ,0, nan\r\nB , nan,0\r\n")
    names, m = read_distance_matrix(path)
    assert names == ["A", "B"] and m[0, 0] == 0 and math.isnan(m[0, 1])


def test_read_distance_matrix_long_format_keeps_the_last_value_of_a_repeated_pair(tmp_path: Path):
    path = tmp_path / "long.tsv"
    path.write_text("seq1\tseq2\tdistance\nA\tB\t7\nB\tA\t9\nA\tC\t1\nA\tB\t4\n", encoding="utf-8")
    names, m = read_distance_matrix(path)
    assert names == ["A", "B", "C"]
    assert m[0, 1] == m[1, 0] == 4 and m[0, 2] == m[2, 0] == 1 and math.isnan(m[1, 2])


def test_read_distance_matrix_rejects_ragged_and_non_numeric_squares(tmp_path: Path):
    ragged = tmp_path / "ragged.tsv"
    ragged.write_text("x\tA\tB\nA\t0\nB\t7\t0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not square"):
        read_distance_matrix(ragged)
    text = tmp_path / "text.tsv"
    text.write_text("x\tA\tB\nA\t0\t-\nB\t7\t0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="non-numeric distance"):
        read_distance_matrix(text)


def test_read_lineages_from_kleborate_v3_and_poppunk(tmp_path: Path):
    kleborate = tmp_path / "k.tsv"
    kleborate.write_text("strain\tklebsiella_pneumo_complex__mlst__ST\nISO1\tST258\n", encoding="utf-8")
    assert read_lineages(kleborate) == {"ISO1": "ST258"}
    poppunk = tmp_path / "clusters.csv"
    poppunk.write_text("Taxon,Cluster\nISO1,4\n", encoding="utf-8")
    assert read_lineages(poppunk) == {"ISO1": "4"}


# --- dereplication -------------------------------------------------------------------------


def test_dereplicate_keeps_one_isolate_per_cluster_per_source():
    rows = [_row(i, "Human") for i in range(3)] + [_row(3, "Animal"), _row(4, "Human"), _row(5, "Food")]
    names = [f"ISO{i}" for i in range(5)]  # ISO5 has no distances
    d = np.full((5, 5), 5000.0)
    np.fill_diagonal(d, 0)
    # 0-1-2 form a human cluster by single linkage (0-2 is 30 SNPs apart); 3 is an animal 1 SNP from 0.
    for i, j, v in [(0, 1, 2), (1, 2, 4), (0, 2, 30), (0, 3, 1)]:
        d[i, j] = d[j, i] = v
    derep, summary = A.dereplicate(A.build_dataset(rows), names, d, threshold=5)
    kept = sorted(iso.key for iso in derep.classified)
    assert kept == ["SAMN0", "SAMN3", "SAMN4", "SAMN5"]  # lowest key of the human cluster; animal never merged
    assert (summary.isolates_before, summary.isolates_after) == (6, 4)
    assert (summary.clusters_collapsed, summary.not_in_matrix) == (1, 1)


def test_dereplication_removes_clonal_inflation_of_feature_tests():
    # Two sources of 10 clusters each; within a cluster all isolates are identical, so a feature
    # carried by 1 of 10 animal clusters (cloned 10x) looks source-associated until clusters collapse.
    rows, names, presence = [], [], {}
    for k in range(20):
        source = "Human" if k < 10 else "Animal"
        for copy in range(10 if k == 10 else 1):
            i = len(rows)
            rows.append(_row(i, source))
            names.append(f"ISO{i}")
            presence[f"ISO{i}"] = {"gene"} if k == 10 else set()
    n = len(rows)
    clone = [int(name[3:]) for name in names]
    cluster = [10 if 10 <= c < 20 else c for c in clone]
    d = np.array([[0.0 if a == b else 5000.0 for b in cluster] for a in cluster])
    ds = A.build_dataset(rows)
    raw = A.feature_prevalence(ds, presence, {}, permutations=999)
    derep, _ = A.dereplicate(ds, names, d, threshold=5)
    collapsed = A.feature_prevalence(derep, presence, {}, permutations=999)
    assert raw.tests[0]["p"] < 0.05
    assert collapsed.tests[0]["p"] > 0.2
    assert n == 29 and len(derep.classified) == 20


def test_collapsed_rows_keep_every_identifier():
    rows = [{**_row(1), "assembly": "GCA_000000001.1"}, {**_row(2, biosample="SAMN1"), "assembly": "GCA_000000002.1"}]
    ds = A.build_dataset(rows)
    assert len(ds.isolates) == 1
    mapping, unmatched = A.match_names(["GCA_000000002.1", "ISO2"], ds)
    assert set(mapping.values()) == {"SAMN1"} and not unmatched


# --- stratified permutation null (hand-built: exact expectations) -----------------------


def test_year_confounded_with_source_gives_no_evidence_under_year_strata():
    # Humans sampled in 2018, animals in 2019; close relatives only within each source.
    rows = [_row(i, "Human", date="2018-03-01") for i in range(4)] + [_row(4 + i, "Animal", date="2019-03-01") for i in range(4)]
    names = [f"ISO{i}" for i in range(8)]
    d = np.full((8, 8), 5000.0)
    np.fill_diagonal(d, 0)
    for i, j in [(0, 1), (2, 3), (4, 5), (6, 7)]:
        d[i, j] = d[j, i] = 2
    ds = A.build_dataset(rows)
    plain = A.linkage(ds, names, d, threshold=5, permutations=999, seed=1)
    strat = A.linkage(ds, names, d, threshold=5, permutations=999, seed=1, strata=("year",))
    cross_plain = next(r for r in plain.events if r["kind"] == "between")
    cross_strat = next(r for r in strat.events if r["kind"] == "between")
    assert cross_plain["events"] == cross_strat["events"] == 0
    assert cross_plain["expected"] > 0.5  # shuffling across years invents human-animal links
    assert (cross_strat["expected"], cross_strat["p"]) == (0.0, 1.0)  # within a year nothing can move


def test_stratified_lineage_expectation_is_the_sum_over_strata():
    # 2018: ST1 x2 (Human), ST2 x2 (Animal). 2019: ST1 x2 and ST2 x2, all Animal.
    rows = [
        _row(0, "Human", "ST1", date="2018"), _row(1, "Human", "ST1", date="2018"),
        _row(2, "Animal", "ST2", date="2018"), _row(3, "Animal", "ST2", date="2018"),
        _row(4, "Animal", "ST1", date="2019"), _row(5, "Animal", "ST1", date="2019"),
        _row(6, "Animal", "ST2", date="2019"), _row(7, "Animal", "ST2", date="2019"),
    ]
    rows = [{**r, "st": r["st"]} for r in rows]
    result = A.lineage_enrichment(A.build_dataset(rows), min_isolates=2, permutations=99, strata=("year",))
    st1_human = next(r for r in result if r["lineage"] == "ST1" and r["source"] == "Human")
    # 2018: 2 ST1 x 2 Human / 4 = 1; 2019: no humans -> 0.
    assert st1_human["expected"] == pytest.approx(1.0)
    unstratified = A.lineage_enrichment(A.build_dataset(rows), min_isolates=2, permutations=99)
    # Pooled: 4 ST1 x 2 Human / 8 = 1 as well, but the variance differs; check the other cell.
    st1_animal = next(r for r in unstratified if r["lineage"] == "ST1" and r["source"] == "Animal")
    assert st1_animal["expected"] == pytest.approx(3.0)


def test_unknown_stratum_field_is_rejected():
    with pytest.raises(ValueError, match="cannot stratify"):
        A.strata_codes(A.build_dataset([_row(1)]).isolates, ("month",))


def test_dataset_keeps_isolates_named_by_strain_and_rows_without_ids():
    rows = [
        {"strain": "KP1", "l1": "Human", "st": "ST1"},
        {"strain": "KP2", "l1": "Animal", "st": "ST1"},
        {"strain": "KP2", "l1": "Animal", "st": "ST1"},  # same isolate listed twice
        {"assembly": "GCA_000005845.2", "l1": "Human", "st": "ST2"},
        {"l1": "Environmental", "st": "ST3"},
        {"l1": "Environmental", "st": "ST3"},
    ]
    ds = A.build_dataset(rows)
    assert (len(ds.isolates), ds.duplicates, ds.unnamed) == (5, 1, 2)
    assert {iso.key for iso in ds.isolates} >= {"KP1", "KP2", "GCA_000005845.2"}
    # Strain names still match the genotype and distance files.
    mapping, unmatched = A.match_names(["KP1", "kp2"], ds)
    assert mapping == {"KP1": "KP1", "kp2": "KP2"} and unmatched == []


def test_feature_permutations_match_a_one_at_a_time_reference():
    """The feature test multiplies many permutations at once; p-values must match doing them one by one."""
    rng = np.random.default_rng(5)
    sources = ["Human", "Animal", "Environmental"]
    rows = [_row(i, sources[i % 3], date=f"{2015 + i % 4}-01-01") for i in range(90)]
    presence = {f"ISO{i}": {f"g{k}" for k in range(12) if rng.random() < 0.15 + 0.3 * (i % 3 == k % 3)} for i in range(90)}
    ds = A.build_dataset(rows)
    permutations = 300  # more than one batch
    result = A.feature_prevalence(ds, presence, {}, permutations=permutations, seed=11, strata=("year",))

    isolates = [iso for iso in ds.classified if iso.key in {f"SAMN{i}" for i in range(90)}]
    features = sorted({f for genes in presence.values() for f in genes})
    groups = [g for g in ds.groups if any(iso.source == g for iso in isolates)]
    src = np.array([groups.index(iso.source) for iso in isolates])
    matrix = np.array([[f in presence[f"ISO{iso.key[4:]}"] for f in features] for iso in isolates], dtype=float)
    n_s, totals, onehot = np.bincount(src), matrix.sum(axis=0), np.eye(len(groups))
    observed = A._chi2(onehot[src].T @ matrix, n_s, totals)
    codes, draw = A.strata_codes(isolates, ("year",)), np.random.default_rng(11)
    exceed = np.zeros(len(features))
    for _ in range(permutations):
        exceed += A._chi2(onehot[A.permute_within(src, codes, draw)].T @ matrix, n_s, totals) >= observed - 1e-9
    expected = {f: A.perm_pvalue(int(e), permutations) for f, e in zip(features, exceed) if 0 < totals[features.index(f)] < len(isolates)}
    assert {t["feature"]: t["p"] for t in result.tests if not math.isnan(t["p"])} == expected
