from pathlib import Path

import pytest

from onehealth.io_tsv import read_id_list, read_tsv


def test_comment_above_header_is_skipped(tmp_path: Path):
    path = tmp_path / "m.tsv"
    path.write_text("# exported 2026-10-01\nid\thost\nSAMN1\tHomo sapiens\n", encoding="utf-8")
    assert read_tsv(path) == [{"id": "SAMN1", "host": "Homo sapiens"}]


def test_trailing_empty_cells_are_ignored_but_extra_values_raise(tmp_path: Path):
    path = tmp_path / "m.tsv"
    path.write_text("id\thost\nSAMN1\tHomo sapiens\t\n", encoding="utf-8")
    assert read_tsv(path)[0]["host"] == "Homo sapiens"
    path.write_text("id\thost\nSAMN1\tHomo sapiens\turine\n", encoding="utf-8")
    with pytest.raises(ValueError, match="more cells than the header"):
        read_tsv(path)


def test_id_list_tsv_with_biosample_column_does_not_return_the_header(tmp_path: Path):
    path = tmp_path / "ids.tsv"
    path.write_text("biosample\thost\nSAMN1\tHomo sapiens\n", encoding="utf-8")
    assert read_id_list(path) == ["SAMN1"]


def test_id_list_tsv_with_leading_comment(tmp_path: Path):
    path = tmp_path / "ids.tsv"
    path.write_text("# my ids\nid\tnote\nSAMN1\tx\n", encoding="utf-8")
    assert read_id_list(path) == ["SAMN1"]


# --- files written by Excel and other tools -------------------------------------------------------

BOM = "﻿"


def test_excel_byte_order_mark_is_not_part_of_the_first_column(tmp_path: Path):
    from onehealth.io_tsv import read_distance_matrix, read_lineages

    path = tmp_path / "m.tsv"
    path.write_text(BOM + "id\thost\nSAMN1\tHomo sapiens\n", encoding="utf-8")
    assert read_tsv(path) == [{"id": "SAMN1", "host": "Homo sapiens"}]
    ids = tmp_path / "ids.txt"
    ids.write_text(BOM + "SAMN1\nSAMN2\n", encoding="utf-8")
    assert read_id_list(ids) == ["SAMN1", "SAMN2"]
    lineages = tmp_path / "clusters.csv"
    lineages.write_text(BOM + "Taxon,Cluster\nI1,7\n", encoding="utf-8")
    assert read_lineages(lineages) == {"I1": "7"}
    matrix = tmp_path / "d.tsv"
    matrix.write_text(BOM + "snp-dists\tA\tB\nA\t0\t3\nB\t3\t0\n", encoding="utf-8")
    assert read_distance_matrix(matrix)[0] == ["A", "B"]


def test_a_stray_quote_stays_in_its_cell_and_does_not_swallow_later_rows(tmp_path: Path):
    path = tmp_path / "m.tsv"
    path.write_text('id\tisolation_source\thost\nI1\t"wound swab\tHomo sapiens\nI2\tblood\tHomo sapiens\n', encoding="utf-8")
    assert read_tsv(path) == [
        {"id": "I1", "isolation_source": '"wound swab', "host": "Homo sapiens"},
        {"id": "I2", "isolation_source": "blood", "host": "Homo sapiens"},
    ]


def test_quoted_cells_from_excel_and_pandas_are_unwrapped(tmp_path: Path):
    path = tmp_path / "m.tsv"
    path.write_text('id\tisolation_source\n"I1"\t"5"" catheter tip"\nI2\t"""clinical"""\n', encoding="utf-8")
    assert [row["isolation_source"] for row in read_tsv(path)] == ['5" catheter tip', '"clinical"']
    assert read_tsv(path)[0]["id"] == "I1"


def test_written_tables_read_back_unchanged_except_for_tabs_and_line_breaks(tmp_path: Path):
    from onehealth.io_tsv import write_tsv

    values = ['5" catheter', '"quoted"', 'ends with "', '"', "", "plain", "a\tb", "two\nlines", "crlf\r\nline"]
    path = tmp_path / "out.tsv"
    write_tsv(path, [{"id": str(i), "value": v} for i, v in enumerate(values)], ["id", "value"])
    expected = [v.replace("\r\n", " ").replace("\t", " ").replace("\n", " ") for v in values]
    assert [row["value"] for row in read_tsv(path)] == expected


def test_spreadsheet_missing_markers_count_as_missing(lexicon):
    from onehealth.classify import classify_record

    for marker in ("NaN", "nan", "#N/A", "not determined", "Undetermined"):
        result = classify_record({"isolation_source": marker}, lexicon)
        assert (result.classify_status, result.classify_reason) == ("review", "unknown_source"), marker
    assert classify_record({"isolation_source": "NaN", "host": "Homo sapiens"}, lexicon).l1 == "Human"


def test_every_missing_value_term_is_reported_as_unknown_not_as_a_species_name(lexicon):
    from onehealth.classify import classify_record
    from onehealth.constants import UNKNOWN_EXACT

    for term in sorted(t for t in UNKNOWN_EXACT if t):
        assert classify_record({"isolation_source": term}, lexicon).classify_reason == "unknown_source", term
