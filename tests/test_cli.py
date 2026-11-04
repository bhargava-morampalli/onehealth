from pathlib import Path

import pytest

from onehealth.cli import main
from onehealth.io_tsv import read_tsv


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


# --- user errors get one line, not a traceback --------------------------------------------------


def _error(capsys, argv) -> tuple[int, str]:
    code = main(argv)
    err = capsys.readouterr().err
    assert "Traceback" not in err
    return code, err.strip().splitlines()[-1]


def test_missing_and_unreadable_inputs_give_one_line_errors(tmp_path: Path, capsys):
    code, line = _error(capsys, ["classify", str(tmp_path / "missing.tsv"), "-o", str(tmp_path / "o.tsv")])
    assert code == 1 and line.startswith("onehealth classify: error: No such file or directory:")
    code, line = _error(capsys, ["summarize", str(tmp_path), "-o", str(tmp_path / "s.tsv")])
    assert code == 1 and "Is a directory" in line
    latin1 = tmp_path / "latin1.tsv"
    latin1.write_bytes("id\tisolation_source\nI1\tsangre São Paulo\n".encode("latin-1"))
    code, line = _error(capsys, ["classify", str(latin1), "-o", str(tmp_path / "o.tsv")])
    assert code == 1 and "not UTF-8 text" in line and "save it as UTF-8" in line


def test_bad_output_paths_give_one_line_errors(tmp_path: Path, capsys, fixtures: Path):
    table = tmp_path / "t.tsv"
    table.write_text("id\tl1\tpaper_rollup\nI1\tHuman\tHuman\n", encoding="utf-8")
    code, line = _error(capsys, ["plot", str(table), "-o", str(tmp_path / "counts.txt")])
    assert code == 1 and "cannot write .txt images" in line
    code, line = _error(capsys, ["report", str(table), "-o", str(table)])
    assert code == 1 and "-o needs a directory" in line


@pytest.mark.parametrize(
    "rules, message",
    [
        ("rules:\n  - l1: Human\n    phrases: [blood]\n", "rule 1 has no id"),
        ("rules:\n  - id: r\n    phrases: [blood]\n", "rule 1 has no l1"),
        ("rules:\n  - id: r\n    l1: Human\n    priority: high\n    phrases: [blood]\n", "priority must be a whole number"),
        ("rules:\n  r: 1\n", "'rules' must be a list"),
        ("rules:\n  - just text\n", "rule 1 is not a mapping"),
        ("rules:\n  - id: r\n    l1: Human\n    phrases: [yes, blood]\n", "put it in quotes"),
        ("rules: [\n  - id: x\n  bad", "not valid YAML"),
    ],
)
def test_broken_lexicons_give_one_line_errors(tmp_path: Path, capsys, fixtures: Path, rules, message):
    lexicon = tmp_path / "lexicon.yaml"
    lexicon.write_text('version: "x"\n' + rules, encoding="utf-8")
    code, line = _error(capsys, ["classify", str(fixtures / "metadata.tsv"), "-o", str(tmp_path / "o.tsv"), "--lexicon", str(lexicon)])
    assert code == 1 and message in line


def test_a_single_phrase_may_be_written_without_a_list(tmp_path: Path):
    from onehealth.classify import classify_record
    from onehealth.lexicon import load_lexicon

    lexicon = tmp_path / "lexicon.yaml"
    lexicon.write_text('version: "x"\nrules:\n  - id: r\n    l1: Human\n    fields: isolation_source\n    phrases: blood\n',
                       encoding="utf-8")
    loaded = load_lexicon(lexicon)
    assert loaded.rules[0].phrases == (("blood",),) and loaded.rules[0].fields == ("isolation_source",)
    assert classify_record({"isolation_source": "blood"}, loaded).l1 == "Human"


@pytest.mark.parametrize("option, value", [("--threshold", "nan"), ("--threshold", "inf"), ("--threshold", "-5"),
                                           ("--threshold", "x"), ("--seed", "-1")])
def test_report_rejects_impossible_thresholds_and_seeds(tmp_path: Path, capsys, option, value):
    with pytest.raises(SystemExit) as exc:
        main(["report", str(tmp_path / "t.tsv"), "-o", str(tmp_path / "r"), option, value])
    assert exc.value.code == 2
    assert option in capsys.readouterr().err


def test_a_bug_still_shows_its_traceback(monkeypatch, tmp_path: Path):
    import onehealth.cli as cli

    def broken(rows):
        raise RuntimeError("a real bug")

    table = tmp_path / "t.tsv"
    table.write_text("id\tl1\n", encoding="utf-8")
    monkeypatch.setattr(cli, "summarize", broken)
    with pytest.raises(RuntimeError, match="a real bug"):
        cli.main(["summarize", str(table), "-o", str(tmp_path / "s.tsv")])


def test_commands_refuse_a_table_that_was_never_classified(tmp_path: Path, capsys, fixtures: Path):
    raw = fixtures / "metadata.tsv"  # resolve-style metadata, no l1 or paper_rollup
    for argv in (["summarize", str(raw), "-o", str(tmp_path / "s.tsv")],
                 ["plot", str(raw), "-o", str(tmp_path / "p.png")],
                 ["report", str(raw), "-o", str(tmp_path / "r")]):
        code, line = _error(capsys, argv)
        assert code == 1 and "no 'l1'" in line and "onehealth classify" in line, argv
    truth = tmp_path / "truth.tsv"
    truth.write_text("id\ttruth\nSAMN1\tHuman\n", encoding="utf-8")
    code, line = _error(capsys, ["validate", str(raw), "--truth", str(truth), "-o", str(tmp_path / "v.tsv")])
    assert code == 1 and "no 'paper_rollup'" in line


def test_a_header_only_classified_table_is_still_fine(tmp_path: Path):
    table = tmp_path / "t.tsv"
    table.write_text("id\tl1\tpaper_rollup\n", encoding="utf-8")
    assert main(["summarize", str(table), "-o", str(tmp_path / "s.tsv")]) == 0
    assert main(["report", str(table), "-o", str(tmp_path / "r")]) == 0


# --- tables that would otherwise send every row to review without saying why -----------------


def test_classify_names_the_problem_with_a_csv_or_an_id_list(tmp_path: Path, capsys):
    csv_file = tmp_path / "m.csv"
    csv_file.write_text("biosample,host,isolation_source\nSAMN1,Homo sapiens,blood\n", encoding="utf-8")
    code, line = _error(capsys, ["classify", str(csv_file), "-o", str(tmp_path / "o.tsv")])
    assert code == 1 and "looks comma-separated" in line
    ids = tmp_path / "ids.txt"
    ids.write_text("SAMN00000001\nSAMN00000002\n", encoding="utf-8")
    code, line = _error(capsys, ["classify", str(ids), "-o", str(tmp_path / "o.tsv")])
    assert code == 1 and "no column classify can read" in line and "found: SAMN00000001" in line
    assert "run onehealth resolve first" in line
    assert not (tmp_path / "o.tsv").exists()


def test_classify_without_evidence_columns_still_applies_overrides(tmp_path: Path, capsys):
    table = tmp_path / "ids.tsv"
    table.write_text("id\tnote\nSAMN1\tx\nSAMN2\ty\n", encoding="utf-8")
    overrides = tmp_path / "overrides.tsv"
    overrides.write_text("id\tl1\nSAMN1\tAnimal\n", encoding="utf-8")
    out = tmp_path / "o.tsv"
    assert main(["classify", str(table), "-o", str(out), "--overrides", str(overrides)]) == 0
    assert "warning:" in capsys.readouterr().err
    assert [(r["id"], r["l1"]) for r in read_tsv(out)] == [("SAMN1", "Animal"), ("SAMN2", "Unclassified")]
    # Without an ID column the overrides cannot match anything, so it is an error after all.
    ids = tmp_path / "ids.txt"
    ids.write_text("SAMN1\nSAMN2\n", encoding="utf-8")
    code, line = _error(capsys, ["classify", str(ids), "-o", str(out), "--overrides", str(overrides)])
    assert code == 1 and "no column classify can read" in line


def test_resolve_metadata_needs_an_id_column(tmp_path: Path, capsys):
    ids = tmp_path / "ids.txt"
    ids.write_text("SAMN1\n", encoding="utf-8")
    metadata = tmp_path / "m.tsv"
    metadata.write_text("strain\thost\tisolation_source\nKP1\tHomo sapiens\tblood\n", encoding="utf-8")
    code, line = _error(capsys, ["resolve", str(ids), "-o", str(tmp_path / "r.tsv"), "--metadata", str(metadata), "--no-network"])
    assert code == 1 and "no ID column to match IDs on" in line and "found: strain, host, isolation_source" in line


def test_truth_and_override_files_in_csv_say_so(tmp_path: Path, capsys, fixtures: Path):
    classified = tmp_path / "c.tsv"
    assert main(["classify", str(fixtures / "metadata.tsv"), "-o", str(classified)]) == 0
    truth = tmp_path / "truth.csv"
    truth.write_text("biosample,truth\nSAMN1,Human\n", encoding="utf-8")
    code, line = _error(capsys, ["validate", str(classified), "--truth", str(truth), "-o", str(tmp_path / "v.tsv")])
    assert code == 1 and "looks comma-separated" in line
    overrides = tmp_path / "o.csv"
    overrides.write_text("id,l1\nSAMN1,Animal\n", encoding="utf-8")
    code, line = _error(capsys, ["classify", str(fixtures / "metadata.tsv"), "-o", str(tmp_path / "o.tsv"), "--overrides", str(overrides)])
    assert code == 1 and "looks comma-separated" in line
