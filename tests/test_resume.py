"""A resolve run that was stopped, killed or refused by NCBI can be finished, and classify does not
take what was never looked up for a record with no source."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from onehealth import analysis as A
from onehealth import cli, resolve
from onehealth.classify import classify_record, classify_rows
from onehealth.io_tsv import read_tsv, write_tsv
from onehealth.ncbi import NCBIClient
from onehealth.overrides import OverrideRow
from onehealth.resolve import answered, earlier_answers, resolve_ids
from test_resolve import MIXED_IDS, DatasetsStub, PathogenDetectionStub, _pd_row, _stub

IDS = [f"SAMN{i}" for i in range(1, 8)] + ["SAMN99999", "GCA_000000002.1", "not-an-id"]


class Stopping(DatasetsStub):
    """NCBI, until the user presses Ctrl-C (or the machine dies) at the `stop_at`-th request."""

    def __init__(self, base: DatasetsStub, stop_at: int, error: type[BaseException] = KeyboardInterrupt):
        super().__init__(base.biosamples, base.genomes, base.max_page, base.fail)
        self.stop_at, self.error = stop_at, error

    def get_json(self, url, headers=None):
        if len(self.urls) + 1 >= self.stop_at:
            self.urls.append(url)
            raise self.error()
        return super().get_json(url, headers)


def _run(monkeypatch, tmp_path: Path, stub, ids=IDS, *extra: str) -> tuple[int, list[dict[str, str]]]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    listed = tmp_path / "ids.txt"
    listed.write_text("\n".join(ids) + "\n", encoding="utf-8")
    monkeypatch.setattr(cli, "NCBIClient", lambda **kw: NCBIClient(fetcher=stub))
    code = cli.main(["resolve", str(listed), "-o", str(tmp_path / "resolved.tsv"), "--email", "x@example.org", *extra])
    return code, read_tsv(tmp_path / "resolved.tsv")


def _asked(stub) -> set[str]:
    """Every accession any request named."""
    import re
    import urllib.parse

    return {a.upper() for url in stub.urls for a in re.findall(r"(?:SAMN|GC[AF]_)[\w.]+", urllib.parse.unquote(url), flags=re.I)}


# --- what an earlier run settled ---------------------------------------------------------------------


@pytest.mark.parametrize(
    "status, reason, groups, settled",
    [
        ("ok", "ncbi", [], True),
        ("ok", "local_metadata", [], True),
        ("ok", "ncbi_assembly", [], True),
        # NCBI answered: there is no such record, or more than one.
        ("review", "biosample_not_found", [], True),
        ("review", "assembly_not_found", [], True),
        ("review", "multiple_biosamples", [], True),
        ("review", "unknown_id", [], True),
        # The lookup did not happen, or NCBI did not reply.
        ("review", "interrupted", [], False),
        ("review", "ncbi_error:HTTP 500 for https://api.ncbi.nlm.nih.gov/x: Internal Server Error", [], False),
        ("review", "missing_metadata", [], False),
        ("review", "", [], False),
        # A PDT is worth another try only when there is now a table to look in.
        ("review", "pdt_unresolved", [], True),
        ("review", "pdt_unresolved", ["Klebsiella"], False),
    ],
)
def test_which_rows_of_an_earlier_run_are_answers(status, reason, groups, settled):
    assert answered({"query_id": "X", "resolve_status": status, "resolve_reason": reason}, groups) is settled


def test_earlier_answers_are_keyed_by_query_id_and_keep_the_first_row():
    rows = [
        {"query_id": "samn1 ", "resolve_status": "ok", "resolve_reason": "ncbi", "isolation_source": "blood"},
        {"query_id": "SAMN1", "resolve_status": "ok", "resolve_reason": "ncbi", "isolation_source": "second"},
        {"query_id": "SAMN2", "resolve_status": "review", "resolve_reason": "interrupted"},
        {"query_id": "", "resolve_status": "ok", "resolve_reason": "ncbi"},
    ]
    kept = earlier_answers(rows)
    assert list(kept) == ["SAMN1"] and kept["SAMN1"]["isolation_source"] == "blood"


# --- resume ------------------------------------------------------------------------------------------


def test_a_stopped_run_is_finished_by_resume_and_gives_the_table_of_one_run(monkeypatch, tmp_path: Path, capsys):
    whole_code, whole = _run(monkeypatch, tmp_path / "whole", _stub())
    assert whole_code == 0

    # Ctrl-C at the fifth request: some IDs are answered, the rest are marked.
    code, stopped = _run(monkeypatch, tmp_path, Stopping(_stub(), stop_at=5))
    assert code == 130
    left = [r["query_id"] for r in stopped if r["resolve_reason"] == "interrupted"]
    assert 0 < len(left) < len(IDS)
    assert "run the same command with --resume" in capsys.readouterr().err

    again = _stub()
    code, finished = _run(monkeypatch, tmp_path, again, IDS, "--resume")
    assert code == 0
    assert finished == whole
    err = capsys.readouterr().err
    assert f"{len(IDS) - len(left)} of {len(IDS)} IDs were answered by the earlier run" in err
    # Only what was left was asked for (and SAMN2, the BioSample of the assembly that was left).
    # "not-an-id" never goes to NCBI.
    assert _asked(again) <= {q.upper() for q in left} | {"SAMN2"}
    assert "SAMN1" not in _asked(again)
    assert not (tmp_path / ".resolved.tsv.partial").exists()


def test_resume_tries_an_ncbi_error_again_but_not_a_record_ncbi_does_not_have(monkeypatch, tmp_path: Path, capsys):
    code, first = _run(monkeypatch, tmp_path, _stub(fail=("SAMN5",)))
    assert code == 0  # the run finished; one ID got no answer
    by_id = {r["query_id"]: r for r in first}
    assert by_id["SAMN5"]["resolve_reason"].startswith("ncbi_error:HTTP 400")
    assert by_id["SAMN99999"]["resolve_reason"] == "biosample_not_found"
    err = capsys.readouterr().err
    assert "warning: NCBI gave no usable answer for" in err and "--resume" in err

    again = _stub()
    code, second = _run(monkeypatch, tmp_path, again, IDS, "--resume")
    assert code == 0
    assert {r["query_id"]: r["resolve_reason"] for r in second}["SAMN5"] == "ncbi"
    assert _asked(again) == {"SAMN5"}
    assert "warning: NCBI gave no usable answer" not in capsys.readouterr().err
    # Everything else is the row the first run wrote.
    assert [r for r in second if r["query_id"] != "SAMN5"] == [r for r in first if r["query_id"] != "SAMN5"]


def test_resume_gives_the_same_table_as_one_run_for_every_kind_of_id(monkeypatch, tmp_path: Path):
    _, whole = _run(monkeypatch, tmp_path / "whole", _stub(), MIXED_IDS)
    for stop_at in (1, 2, 4, 7, 11):
        folder = tmp_path / f"stop{stop_at}"
        code, _ = _run(monkeypatch, folder, Stopping(_stub(), stop_at=stop_at), MIXED_IDS)
        assert code == 130
        code, finished = _run(monkeypatch, folder, _stub(), MIXED_IDS, "--resume")
        assert code == 0 and finished == whole, stop_at


def test_resume_with_no_earlier_output_looks_everything_up(monkeypatch, tmp_path: Path, capsys):
    code, rows = _run(monkeypatch, tmp_path, _stub(), IDS, "--resume")
    assert code == 0 and len(rows) == len(IDS)
    assert "--resume found no earlier output" in capsys.readouterr().err


def test_without_resume_an_earlier_output_is_not_read(monkeypatch, tmp_path: Path):
    _run(monkeypatch, tmp_path, _stub())
    again = _stub()
    _run(monkeypatch, tmp_path, again)
    assert {"SAMN1", "SAMN7"} <= _asked(again)


def test_resume_looks_an_unresolved_pdt_up_again_once_an_organism_group_is_named(monkeypatch, tmp_path: Path):
    table = [_pd_row("PDT000000003.2", "SAMN3", "GCA_000000003.1")]
    code, first = _run(monkeypatch, tmp_path, PathogenDetectionStub({"Klebsiella": table}), ["PDT000000003.2"])
    assert (code, first[0]["resolve_reason"]) == (0, "pdt_unresolved")
    # The same command again: the PDT was answered ("not found by assembly name"), so nothing is asked.
    quiet = PathogenDetectionStub({"Klebsiella": table})
    _, same = _run(monkeypatch, tmp_path, quiet, ["PDT000000003.2"], "--resume")
    assert same == first and quiet.urls == []
    # With the organism group, the table is read and the PDT is found.
    _, found = _run(monkeypatch, tmp_path, PathogenDetectionStub({"Klebsiella": table}), ["PDT000000003.2"],
                    "--resume", "--pd-organism", "Klebsiella")
    assert (found[0]["resolve_status"], found[0]["resolve_reason"], found[0]["biosample"]) == ("ok", "ncbi_pdt", "SAMN3")


def test_resume_reads_the_pathogen_detection_table_only_for_pdts_still_open(monkeypatch, tmp_path: Path):
    table = [_pd_row(f"PDT00000000{i}.1", f"SAMN{i}", f"GCA_00000000{i}.1") for i in range(1, 6)]
    ids = [f"PDT00000000{i}.1" for i in range(1, 6)]
    _, first = _run(monkeypatch, tmp_path, PathogenDetectionStub({"Klebsiella": table}), ids, "--pd-organism", "Klebsiella")
    assert [r["resolve_status"] for r in first] == ["ok"] * 5
    again = PathogenDetectionStub({"Klebsiella": table})
    _, second = _run(monkeypatch, tmp_path, again, ids, "--resume", "--pd-organism", "Klebsiella")
    assert second == first
    assert again.urls == []  # every PDT was answered: the table of 200 MB or more is not downloaded again


# --- checkpoints: a killed run keeps what it had --------------------------------------------------------


def test_the_ids_are_looked_up_a_chunk_at_a_time_with_a_checkpoint_after_each(monkeypatch):
    monkeypatch.setattr(resolve, "CHUNK", 3)
    seen = []

    def checkpoint(rows, done, total):
        seen.append((done, total, sum(1 for r in rows if r["resolve_reason"]), len(rows)))

    ids = [f"SAMN{i}" for i in range(1, 9)] + ["not-an-id"]
    rows = resolve_ids(ids, client=NCBIClient(fetcher=_stub()), allow_network=True, checkpoint=checkpoint)
    # Eight IDs go to NCBI in chunks of three; "not-an-id" was settled before any request.
    assert seen == [(3, 8, 4, 9), (6, 8, 7, 9), (8, 8, 9, 9)]
    assert rows == resolve_ids(ids, client=NCBIClient(fetcher=_stub()), allow_network=True)


def test_a_killed_run_has_saved_all_but_its_last_chunk_and_resume_finishes_it(monkeypatch, tmp_path: Path, capsys):
    monkeypatch.setattr(resolve, "CHUNK", 3)
    monkeypatch.setattr(cli, "CHECKPOINT_SECONDS", 0.0)
    _, whole = _run(monkeypatch, tmp_path / "whole", _stub())

    class Killed(Exception):
        """Stands for a power cut: nothing in onehealth handles it."""

    # One request per chunk of three BioSamples; the third chunk has an assembly, which takes
    # a request of its own first. The fourth request is where the machine dies.
    with pytest.raises(Killed):
        _run(monkeypatch, tmp_path, Stopping(_stub(), stop_at=4, error=Killed))
    saved = read_tsv(tmp_path / "resolved.tsv")
    assert len(saved) == len(IDS)  # a whole table, not half of one
    assert [r["resolve_status"] for r in saved[:6]] == ["ok"] * 6
    open_ids = [r["query_id"] for r in saved if r["resolve_reason"] == "interrupted"]
    assert open_ids == ["SAMN7", "SAMN99999", "GCA_000000002.1"]
    assert [r["query_id"] for r in read_tsv(tmp_path / "resolved.review.tsv")] == [*open_ids, "not-an-id"]
    assert "IDs looked up at NCBI; saved to" in capsys.readouterr().err

    again = _stub()
    code, finished = _run(monkeypatch, tmp_path, again, IDS, "--resume")
    assert code == 0 and finished == whole
    assert not {"SAMN1", "SAMN3", "SAMN4", "SAMN5", "SAMN6"} & _asked(again)


def test_checkpoints_are_not_written_more_often_than_every_thirty_seconds(monkeypatch, tmp_path: Path, capsys):
    monkeypatch.setattr(resolve, "CHUNK", 2)
    assert cli.CHECKPOINT_SECONDS == 30.0
    code, _ = _run(monkeypatch, tmp_path, _stub())
    assert code == 0 and "saved to" not in capsys.readouterr().err  # a run of seconds writes once, at the end


def test_a_table_is_written_whole_or_not_at_all(monkeypatch, tmp_path: Path):
    path = tmp_path / "table.tsv"
    write_tsv(path, [{"a": "1"}, {"a": "2"}], ["a"])
    before = path.read_text(encoding="utf-8")
    real = csv.DictWriter.writerow

    def fails_on_the_second_row(self, row):
        if row["a"] == "second":
            raise OSError("No space left on device")
        return real(self, row)

    monkeypatch.setattr(csv.DictWriter, "writerow", fails_on_the_second_row)
    with pytest.raises(OSError):
        write_tsv(path, [{"a": "first"}, {"a": "second"}], ["a"])
    assert path.read_text(encoding="utf-8") == before
    assert [p.name for p in tmp_path.iterdir()] == ["table.tsv"]


# --- classify: not looked up is not "no source" ---------------------------------------------------------


@pytest.mark.parametrize(
    "resolve_reason, classify_reason",
    [
        ("ncbi_error:HTTP 500 for https://api.ncbi.nlm.nih.gov/x: Internal Server Error", "unresolved:ncbi_error"),
        ("interrupted", "unresolved:interrupted"),
        ("missing_metadata", "unresolved:missing_metadata"),
        ("biosample_not_found", "unresolved:biosample_not_found"),
        ("pdt_unresolved", "unresolved:pdt_unresolved"),
        ("multiple_biosamples", "unresolved:multiple_biosamples"),
        ("", "unresolved:review"),
    ],
)
def test_a_row_resolve_could_not_get_is_unresolved_not_missing_source(lexicon, resolve_reason, classify_reason):
    row = {"query_id": "SAMN1", "biosample": "SAMN1", "host": "", "isolation_source": "",
           "resolve_status": "review", "resolve_reason": resolve_reason}
    result = classify_record(row, lexicon)
    assert (result.l1, result.classify_status, result.classify_reason) == ("Unclassified", "review", classify_reason)


def test_a_record_that_was_fetched_and_has_no_source_is_still_missing_source(lexicon):
    fetched = {"biosample": "SAMN1", "host": "", "isolation_source": "", "resolve_status": "ok", "resolve_reason": "ncbi"}
    assert classify_record(fetched, lexicon).classify_reason == "missing_source"
    # A table that never went through resolve has no resolve columns.
    assert classify_record({"biosample": "SAMN1", "host": "", "isolation_source": ""}, lexicon).classify_reason == "missing_source"
    # A missing-value word is a statement by the submitter, not a failed lookup.
    unknown = {"isolation_source": "not collected", "host": "", "resolve_status": "ok", "resolve_reason": "ncbi"}
    assert classify_record(unknown, lexicon).classify_reason == "unknown_source"


def test_an_unresolved_row_with_evidence_or_an_override_is_classified(lexicon):
    # resolve sent it to review, but the table has a source (say, from --metadata): that is evidence.
    row = {"biosample": "SAMN1", "host": "Homo sapiens", "isolation_source": "blood",
           "resolve_status": "review", "resolve_reason": "ncbi_error:HTTP 500"}
    assert classify_record(row, lexicon).classify_reason == "specimens+human-host"
    bare = {"biosample": "SAMN2", "host": "", "isolation_source": "", "resolve_status": "review", "resolve_reason": "interrupted"}
    override = {"SAMN2": [OverrideRow(key="SAMN2", l1="Animal", healthcare_environment="", paper_rollup="", note="", complete=True)]}
    assert classify_record(bare, lexicon, override).classify_reason == "override"


def test_classify_warns_about_unresolved_rows_and_report_counts_them(tmp_path: Path, capsys):
    table = tmp_path / "resolved.tsv"
    rows = [
        {"query_id": "SAMN1", "biosample": "SAMN1", "host": "Homo sapiens", "isolation_source": "blood",
         "resolve_status": "ok", "resolve_reason": "ncbi"},
        {"query_id": "SAMN2", "biosample": "SAMN2", "host": "", "isolation_source": "",
         "resolve_status": "ok", "resolve_reason": "ncbi"},
        {"query_id": "SAMN3", "biosample": "SAMN3", "host": "", "isolation_source": "",
         "resolve_status": "review", "resolve_reason": "ncbi_error:HTTP 500"},
        {"query_id": "SAMN4", "biosample": "SAMN4", "host": "", "isolation_source": "",
         "resolve_status": "review", "resolve_reason": "interrupted"},
        {"query_id": "SAMN5", "biosample": "SAMN5", "host": "", "isolation_source": "",
         "resolve_status": "review", "resolve_reason": "biosample_not_found"},
    ]
    write_tsv(table, rows, list(rows[0]))
    out = tmp_path / "classified.tsv"
    assert cli.main(["classify", str(table), "-o", str(out)]) == 0
    err = capsys.readouterr().err
    assert "warning: 3 of 5 rows were not resolved at NCBI" in err
    assert "their source is unknown, not missing" in err
    assert "2 of them were never looked up or got an error; run resolve again with --resume" in err
    classified = read_tsv(out)
    assert [r["classify_reason"] for r in classified] == [
        "specimens+human-host", "missing_source", "unresolved:ncbi_error", "unresolved:interrupted",
        "unresolved:biosample_not_found"]
    dataset = A.build_dataset(classified)
    assert (dataset.unclassified, dataset.unresolved) == (4, 3)

    # Nothing unresolved: no warning.
    write_tsv(table, rows[:2], list(rows[0]))
    assert cli.main(["classify", str(table), "-o", str(out)]) == 0
    assert "were not resolved" not in capsys.readouterr().err


def test_the_whole_way_an_ncbi_error_then_resume_then_classify(monkeypatch, tmp_path: Path, capsys, lexicon):
    """What the fix is for: an NCBI outage in the middle of a study must not look like data."""
    _run(monkeypatch, tmp_path, _stub(fail=("SAMN5",)))
    first = classify_rows(read_tsv(tmp_path / "resolved.tsv"), lexicon)
    assert {r["query_id"]: r["classify_reason"] for r in first}["SAMN5"] == "unresolved:ncbi_error"
    _run(monkeypatch, tmp_path, _stub(), IDS, "--resume")
    second = classify_rows(read_tsv(tmp_path / "resolved.tsv"), lexicon)
    assert not [r for r in second if r["classify_reason"].startswith("unresolved:ncbi_error")]
    assert {r["query_id"]: r["classify_status"] for r in second}["SAMN5"] in {"ok", "review"}
    assert {r["query_id"]: r["isolation_source"] for r in second}["SAMN5"] == "source 5"
