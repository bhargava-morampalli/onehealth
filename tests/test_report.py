import re
from pathlib import Path

import numpy as np
import pytest

from onehealth.cli import main


def _write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    classified = tmp_path / "classified.tsv"
    lines = ["query_id\tl1\thealthcare_environment\tpaper_rollup\tst\tcollection_date\tisolation_source"]
    sources = ["Human"] * 12 + ["Animal"] * 8 + ["Environmental"] * 6 + ["Unclassified"] * 2
    for i, l1 in enumerate(sources):
        hc = "yes" if l1 == "Environmental" and i % 2 else "no"
        lines.append(f"ISO{i}\t{l1}\t{hc}\t{l1}\tST{i % 5}\t{2018 + i % 3}-06-01\tx")
    classified.write_text("\n".join(lines) + "\n", encoding="utf-8")
    features = tmp_path / "amr.tsv"
    feature_lines = ["Name\tElement symbol\tClass", "ISO0\t<b>bad</b>\tBETA-LACTAM"]
    feature_lines += [f"ISO{i}\ttet(A)\tTETRACYCLINE" for i in range(12, 20)]
    feature_lines += [f"ISO{i}\tblaSHV-1\tBETA-LACTAM" for i in range(28)]
    features.write_text("\n".join(feature_lines) + "\n", encoding="utf-8")
    distances = tmp_path / "dists.tsv"
    n = len(sources)
    d = np.full((n, n), 4000)
    np.fill_diagonal(d, 0)
    for i, j in [(0, 1), (1, 12), (12, 13), (20, 21), (2, 20)]:
        d[i, j] = d[j, i] = 3
    names = [f"ISO{i}" for i in range(n)]
    rows = ["snp-dists\t" + "\t".join(names)] + [names[i] + "\t" + "\t".join(map(str, d[i])) for i in range(n)]
    distances.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return classified, features, distances


def test_report_writes_html_figures_and_tables(tmp_path: Path):
    classified, features, distances = _write_inputs(tmp_path)
    out = tmp_path / "report"
    code = main([
        "report", str(classified), "-o", str(out), "--features", str(features), "--distances", str(distances),
        "--permutations", "99", "--threshold", "5", "--note", "test <data>",
    ])
    assert code == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    for heading in ("Sampling by source", "Lineage diversity", "features by source", "Clusters that span sources", "Methods"):
        assert heading in page
    # User-supplied text is escaped, never injected as markup.
    assert "<b>bad</b>" not in page and "test <data>" not in page
    assert "test &lt;data&gt;" in page
    # Inline SVG is re-themed with CSS variables (dark mode), standalone files keep hex colours.
    assert "var(--s-human)" in page
    assert "#2a78d6" in (out / "figures" / "01_sampling.svg").read_text(encoding="utf-8")
    assert (out / "figures" / "06_features.png").stat().st_size > 0
    assert (out / "tables" / "09_linked_isolates.tsv").exists()


def test_report_without_optional_inputs_skips_those_sections(tmp_path: Path):
    classified, _, _ = _write_inputs(tmp_path)
    out = tmp_path / "report"
    assert main(["report", str(classified), "-o", str(out), "--permutations", "19"]) == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "Lineage diversity" in page
    assert "Clusters that span sources" not in page
    assert "features by source" not in page


def test_report_dereplicate_uses_representatives_and_needs_distances(tmp_path: Path, capsys):
    classified, features, distances = _write_inputs(tmp_path)
    out = tmp_path / "report"
    assert main(["report", str(classified), "-o", str(out), "--dereplicate", "--permutations", "19"]) == 2
    assert "--dereplicate needs --distances" in capsys.readouterr().err
    code = main([
        "report", str(classified), "-o", str(out), "--features", str(features), "--distances", str(distances),
        "--threshold", "5", "--permutations", "19", "--dereplicate",
    ])
    assert code == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "dereplicated: yes" in page and "After dereplication" in page
    assert "Vezina B, Morampalli BR" in page
    assert "one representative per genomic cluster per source" in page


def _kpi(page: str, label: str) -> str:
    match = re.search(rf'<div class="label">{re.escape(label)}</div><div class="value">([^<]*)</div>', page)
    assert match, f"no {label!r} tile"
    return match.group(1)


def test_report_counts_isolates_named_only_by_strain(tmp_path: Path):
    classified = tmp_path / "classified.tsv"
    lines = ["strain\tl1\thealthcare_environment\tpaper_rollup\tst\tcollection_date"]
    lines += [f"KP{i}\t{l1}\tno\t{l1}\tST{i % 3}\t2019-01-01" for i, l1 in enumerate(["Human"] * 5 + ["Animal"] * 4)]
    classified.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out = tmp_path / "report"
    assert main(["report", str(classified), "-o", str(out), "--permutations", "19"]) == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert _kpi(page, "Isolates") == "9"
    assert "Lineage diversity" in page
    assert "were collapsed" not in page


def test_diversity_headline_skips_a_source_with_one_typed_isolate(tmp_path: Path):
    # Human comes first in the fixed source order and has a single typed isolate, so no D.
    classified = tmp_path / "classified.tsv"
    lines = ["query_id\tl1\thealthcare_environment\tpaper_rollup\tst"]
    lines += ["H1\tHuman\tno\tHuman\tST1"] + [f"A{i}\tAnimal\tno\tAnimal\tST{i % 2}" for i in range(4)]
    classified.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out = tmp_path / "report"
    assert main(["report", str(classified), "-o", str(out), "--permutations", "19"]) == 0
    finding = re.search(r'id="02_diversity-title">.*?<p class="finding">([^<]*)</p>', (out / "report.html").read_text(
        encoding="utf-8"), re.S).group(1)
    assert "nan" not in finding
    assert finding.startswith("Only Animal has two or more typed isolates")


def test_zero_permutations_give_a_report_that_says_the_tests_were_skipped(tmp_path: Path):
    classified, features, distances = _write_inputs(tmp_path)
    out = tmp_path / "report"
    code = main(["report", str(classified), "-o", str(out), "--features", str(features), "--distances", str(distances),
                 "--threshold", "5", "--permutations", "0"])
    assert code == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "Permutation tests were skipped (--permutations 0)" in page
    assert "Enrichment was not tested" in page and "were not tested (--permutations 0)" in page
    assert "Is cross-source sharing more or less than chance?" not in page
    assert "No lineage was significantly enriched" not in page
    assert (out / "tables" / "08_sharing_events.tsv").exists()


@pytest.mark.parametrize("value, message", [("-1", "must be 0 or more"), ("x", "not a whole number")])
def test_report_rejects_a_bad_permutation_count(tmp_path: Path, capsys, value, message):
    classified, _, _ = _write_inputs(tmp_path)
    with pytest.raises(SystemExit) as exc:
        main(["report", str(classified), "-o", str(tmp_path / "r"), "--permutations", value])
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


# --- edge cases found by the sweep ---------------------------------------------------------------


def _quality(page: str) -> list[str]:
    import html as _html

    block = page.split('id="quality"')[1].split("</section>")[0]
    return [_html.unescape(re.sub("<[^>]+>", "", item)) for item in re.findall(r"<li>(.*?)</li>", block)]


def _table(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_report_on_a_single_isolate(tmp_path: Path):
    table = _table(tmp_path, "one.tsv", "id\tl1\tpaper_rollup\thealthcare_environment\tst\nI1\tHuman\tHuman\tno\tST1\n")
    assert main(["report", str(table), "-o", str(tmp_path / "r"), "--permutations", "9"]) == 0
    page = (tmp_path / "r" / "report.html").read_text(encoding="utf-8")
    assert "Lineage accumulation needs a source with at least two typed isolates" in " ".join(_quality(page))
    assert "Lineage accumulation</h2>" not in page


def test_a_threshold_beyond_the_histogram_axis_is_described_not_drawn(tmp_path: Path):
    classified, _, distances = _write_inputs(tmp_path)
    out = tmp_path / "r"
    assert main(["report", str(classified), "-o", str(out), "--distances", str(distances), "--threshold", "1e9",
                 "--permutations", "9"]) == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "threshold is beyond this axis" in page
    assert (out / "figures" / "07_distances.png").stat().st_size < 2_000_000


def test_skipped_report_sections_are_explained(tmp_path: Path):
    classified, _, _ = _write_inputs(tmp_path)
    empty = _table(tmp_path, "empty.tsv", "")
    same_gene = _table(tmp_path, "same.tsv", "Name\tElement symbol\n" + "".join(f"ISO{i}\tfosA\n" for i in range(28)))
    stranger = _table(tmp_path, "stranger.tsv", "x\tA\tB\nA\t0\t1\nB\t1\t0\n")
    cases = [
        (["--features", str(empty)], "lists no isolates, so the features section was skipped"),
        (["--features", str(same_gene)], "nothing to compare by source"),
        (["--distances", str(empty)], "genomic linkage sections were skipped"),
        (["--distances", str(stranger)], "genomic linkage sections were skipped"),
        (["--lineages", str(empty)], "lists no isolates"),
    ]
    for extra, note in cases:
        out = tmp_path / "r"
        assert main(["report", str(classified), "-o", str(out), "--permutations", "9", *extra]) == 0, extra
        assert note in " ".join(_quality((out / "report.html").read_text(encoding="utf-8"))), extra


def test_unclassified_isolates_with_lineages_get_the_right_reason(tmp_path: Path):
    table = _table(tmp_path, "u.tsv", "id\tl1\tpaper_rollup\thealthcare_environment\tst\n"
                   + "".join(f"I{i}\tUnclassified\tUnclassified\t\tST1\n" for i in range(4)))
    assert main(["report", str(table), "-o", str(tmp_path / "r"), "--permutations", "9"]) == 0
    notes = " ".join(_quality((tmp_path / "r" / "report.html").read_text(encoding="utf-8")))
    assert "No isolate with a source class has a lineage in 'st'" in notes
    assert "No lineage column was found" not in notes


def test_report_flags_mash_distances_unreadable_dates_and_untyped_lineages(tmp_path: Path):
    """Each of these used to pass without a word and change the results."""
    classified = tmp_path / "classified.tsv"
    lines = ["query_id\tl1\thealthcare_environment\tpaper_rollup\tst\tcollection_date"]
    for i in range(12):
        l1 = "Human" if i < 7 else "Animal"
        st = "NF" if i in (0, 8) else f"ST{i % 3}"
        date = "week 14" if i == 1 else "15-Mar-2019" if i == 2 else "2019-06-01"
        lines.append(f"ISO{i}\t{l1}\tno\t{l1}\t{st}\t{date}")
    classified.write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = [f"ISO{i}" for i in range(12)]
    mash = np.full((12, 12), 0.0123)
    np.fill_diagonal(mash, 0)
    distances = tmp_path / "mash.tsv"
    distances.write_text("\n".join(["\t" + "\t".join(names)] + [names[i] + "\t" + "\t".join(f"{v:g}" for v in mash[i])
                                                                 for i in range(12)]) + "\n", encoding="utf-8")
    out = tmp_path / "report"
    assert main(["report", str(classified), "-o", str(out), "--distances", str(distances), "--permutations", "19"]) == 0
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "The distances in mash.tsv are fractions between 0 and 1" in page
    assert "2 isolates have a lineage that means no type was assigned" in page
    assert "1 isolates have a collection date in a format that could not be read" in page
    assert "Simpson EH (1949) Measurement of diversity" in page
    # The two NF isolates are untyped, not a lineage shared by Human and Animal.
    lineages = (out / "tables" / "04_lineage_sources.tsv").read_text(encoding="utf-8")
    assert "NF" not in lineages

    # The same table with SNP counts: none of the three notes about the distances.
    snps = np.full((12, 12), 300)
    np.fill_diagonal(snps, 0)
    distances.write_text("\n".join(["\t" + "\t".join(names)] + [names[i] + "\t" + "\t".join(map(str, snps[i]))
                                                                 for i in range(12)]) + "\n", encoding="utf-8")
    assert main(["report", str(classified), "-o", str(out), "--distances", str(distances), "--permutations", "19"]) == 0
    assert "fractions between 0 and 1" not in (out / "report.html").read_text(encoding="utf-8")
