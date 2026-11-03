from pathlib import Path

import pytest

from onehealth.cli import main


def test_cli_help_exits_zero(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "resolve" in out
    assert "classify" in out
    assert "summarize" in out
    assert "plot" in out


def test_subcommand_help(capsys):
    for command in ("resolve", "classify", "summarize", "plot"):
        with pytest.raises(SystemExit) as exc:
            main([command, "--help"])
        assert exc.value.code == 0
        assert command in capsys.readouterr().out


def test_offline_pipeline(tmp_path: Path, fixtures: Path):
    ids = fixtures / "ids.txt"
    metadata = fixtures / "metadata.tsv"
    overrides = fixtures / "overrides.tsv"
    resolved = tmp_path / "resolved.tsv"
    classified = tmp_path / "classified.tsv"
    summary = tmp_path / "summary.tsv"
    plot = tmp_path / "counts.png"

    assert main(["resolve", str(ids), "-o", str(resolved), "--metadata", str(metadata), "--no-network"]) == 0
    assert main(["classify", str(resolved), "-o", str(classified), "--overrides", str(overrides)]) == 0
    assert main(["summarize", str(classified), "-o", str(summary)]) == 0
    assert main(["plot", str(classified), "-o", str(plot)]) == 0

    classified_text = classified.read_text(encoding="utf-8")
    assert "wound drainage" in classified_text
    assert plot.exists() and plot.stat().st_size > 0
    assert "cluster_overlap" in summary.read_text(encoding="utf-8")
    review = classified.with_name("classified.review.tsv")
    assert review.exists()
    assert "OVERRIDE_EMPTY" in review.read_text(encoding="utf-8")


def test_resolve_keeps_cluster_column_when_first_row_has_none(tmp_path: Path):
    ids = tmp_path / "ids.txt"
    ids.write_text("not-an-id\nSAMN1\n", encoding="utf-8")
    metadata = tmp_path / "meta.tsv"
    metadata.write_text("biosample\tisolation_source\tcluster_id\nSAMN1\turine\tST15\n", encoding="utf-8")
    out = tmp_path / "resolved.tsv"
    assert main(["resolve", str(ids), "-o", str(out), "--metadata", str(metadata), "--no-network"]) == 0
    header, _, second = out.read_text(encoding="utf-8").splitlines()
    assert "cluster_id" in header.split("\t")
    assert "ST15" in second.split("\t")


def test_help_lists_every_command_and_the_current_seed_lexicon(capsys):
    from onehealth.constants import SEED_LEXICON_VERSION

    with pytest.raises(SystemExit):
        main(["--help"])
    out = capsys.readouterr().out
    assert "validate" in out and "report" in out
    with pytest.raises(SystemExit):
        main(["classify", "--help"])
    assert f"seeded v{SEED_LEXICON_VERSION}" in " ".join(capsys.readouterr().out.split())
