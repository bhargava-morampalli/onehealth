from pathlib import Path

import numpy as np

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
