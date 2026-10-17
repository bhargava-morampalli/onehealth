"""Tests on real public data: NCBI Pathogen Detection Klebsiella (see tests/data/ncbi_pd_klebsiella/README.md)."""

import csv
from pathlib import Path

import numpy as np
import pytest

from onehealth import analysis as A
from onehealth.classify import classify_rows
from onehealth.cli import main
from onehealth.io_tsv import read_distance_matrix, read_features, read_tsv

DATA = Path(__file__).parent / "data" / "ncbi_pd_klebsiella"


@pytest.fixture(scope="module")
def classified(lexicon):
    return classify_rows(read_tsv(DATA / "metadata.tsv"), lexicon)


@pytest.fixture(scope="module")
def distances():
    return read_distance_matrix(DATA / "snp_dists.tsv")


def test_vcf_distances_reproduce_ncbi_snp_clusters(distances):
    """Single linkage at NCBI's 50-SNP cluster definition recovers 21 of 22 NCBI clusters."""
    names, matrix = distances
    cluster = {r["target_acc"]: r["PDS_acc"] for r in csv.DictReader(open(DATA / "metadata.tsv"), delimiter="\t")}
    iu, ju = np.triu_indices(len(names), k=1)
    linked = np.isfinite(matrix[iu, ju]) & (matrix[iu, ju] <= 50)
    components = A._components(len(names), iu[linked], ju[linked])
    parts: dict[str, set[int]] = {}
    for name, component in zip(names, components.tolist()):
        parts.setdefault(cluster[name], set()).add(component)
    assert len(parts) == 22
    assert sum(len(p) == 1 for p in parts.values()) >= 21


def test_historical_country_names_share_a_stratum(classified):
    dataset = A.build_dataset(classified)
    countries = {iso.country for iso in dataset.isolates}
    assert "Czechia" in countries and "Czech Republic" not in countries


def test_permutation_within_real_year_country_strata_keeps_each_stratum_composition(classified):
    isolates = A.build_dataset(classified).classified
    strata = A.strata_codes(isolates, ("year", "country"))
    labels = np.array([iso.source for iso in isolates], dtype=object)
    rng = np.random.default_rng(1)
    moved = 0
    for _ in range(50):
        shuffled = A.permute_within(labels, strata, rng)
        moved += int((shuffled != labels).sum())
        for k in np.unique(strata):
            inside = strata == k
            assert sorted(shuffled[inside]) == sorted(labels[inside])
    assert moved > 0


def test_stratified_linkage_changes_only_the_null(classified, distances):
    names, matrix = distances
    dataset = A.build_dataset(classified)
    plain = A.linkage(dataset, names, matrix, threshold=21, permutations=199, seed=7)
    strat = A.linkage(dataset, names, matrix, threshold=21, permutations=199, seed=7, strata=("year", "country"))
    again = A.linkage(dataset, names, matrix, threshold=21, permutations=199, seed=7, strata=("year", "country"))
    # 307 NCBI targets, but SAMN02996352 and SAMN02996353 each have two assemblies: 305 isolates.
    assert plain.matched == 305 and not plain.unmatched
    assert [r["events"] for r in plain.events] == [r["events"] for r in strat.events]
    assert [r["expected"] for r in plain.events] != [r["expected"] for r in strat.events]
    assert [r["p"] for r in strat.events] == [r["p"] for r in again.events]
    assert all(0 < r["p"] <= 1 for r in strat.events)


def test_ncbi_amrfinder_features_and_stratified_tests_run(classified):
    presence, classes = read_features(DATA / "features.tsv")
    dataset = A.build_dataset(classified)
    result = A.feature_prevalence(dataset, presence, classes, permutations=199, strata=("country",))
    assert not result.unmatched
    assert sum(result.genotyped.values()) >= 250
    kpc = [t for t in result.tests if t["feature"].startswith("blaKPC")]
    assert kpc and all(0 < t["p"] <= 1 for t in kpc)


def test_report_on_real_data_with_strata_and_dereplication(tmp_path: Path, classified):
    table = tmp_path / "classified.tsv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(classified[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(classified)
    out = tmp_path / "report"
    code = main([
        "report", str(table), "-o", str(out),
        "--features", str(DATA / "features.tsv"), "--distances", str(DATA / "snp_dists.tsv"),
        "--lineage-column", "PDS_acc", "--stratify-by", "year,country", "--dereplicate", "--permutations", "99",
    ])
    assert code == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "permutation strata: year, country" in page
    assert "shuffled source labels only among isolates with the same collection year and country" in page
