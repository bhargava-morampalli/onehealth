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
