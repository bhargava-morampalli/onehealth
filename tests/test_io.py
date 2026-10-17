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
