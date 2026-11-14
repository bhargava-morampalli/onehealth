"""Input, NCBI and provenance fixes from the November 2026 review (see CHANGELOG.md)."""

import json
from pathlib import Path

import pytest

from onehealth import __version__, cli, ncbi
from onehealth.io_tsv import read_distance_matrix, read_id_list, read_lineages, read_tsv, write_tsv
from onehealth.lexicon import default_lexicon_path, load_lexicon
from onehealth.ncbi import NCBIClient, attributes_from_biosample, parse_biosample_xml
from onehealth.resolve import resolve_ids

# --- tables whose header starts with '#' ----------------------------------------------------------

PD_HEADER = "#label\tHHS_region\tasm_acc\tbiosample_acc\thost\tisolation_source\tsource_type\ttarget_acc"
PD_ROW = "PDT000000011.3|SAMN02265470\t5\tGCA_004104335.1\tSAMN02265470\tNULL\tblood\thuman\tPDT000000011.3"


def test_ncbi_pathogen_detection_metadata_is_read_as_downloaded(tmp_path: Path):
    """Its header starts with '#'. It was dropped as a comment and the first isolate became the header."""
    path = tmp_path / "PDG000000001.1.metadata.tsv"
    path.write_text(f"{PD_HEADER}\n{PD_ROW}\n", encoding="utf-8")
    rows = read_tsv(path)
    assert len(rows) == 1 and list(rows[0])[0] == "label"
    assert rows[0]["biosample_acc"] == "SAMN02265470" and rows[0]["isolation_source"] == "blood"
    assert read_id_list(path) == ["SAMN02265470"]
    out = tmp_path / "classified.tsv"
    assert cli.main(["classify", str(path), "-o", str(out)]) == 0
    assert read_tsv(out)[0]["l1"] == "Human"


def test_assembly_summary_header_after_a_comment_line_is_read(tmp_path: Path):
    path = tmp_path / "assembly_summary.txt"
    path.write_text(
        "## See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n"
        "#assembly_accession\tbioproject\tbiosample\torganism_name\n"
        "GCF_000005845.2\tPRJNA57779\tSAMN02604091\tEscherichia coli\n",
        encoding="utf-8",
    )
    assert read_tsv(path) == [
        {"assembly_accession": "GCF_000005845.2", "bioproject": "PRJNA57779", "biosample": "SAMN02604091",
         "organism_name": "Escherichia coli"}
    ]


@pytest.mark.parametrize(
    "above",
    [
        "# exported 2026-10-01",  # prose
        "#SAMN0\tBos taurus\ttissue",  # a commented-out row: values, not column names
        "#id\thost\tisolation_source",  # a commented-out copy of the header
        "# id\thost",  # an older, shorter header
    ],
)
def test_other_hash_lines_above_the_header_stay_comments(tmp_path: Path, above):
    path = tmp_path / "m.tsv"
    path.write_text(f"{above}\nid\thost\tisolation_source\nSAMN1\tHomo sapiens\tblood\n", encoding="utf-8")
    assert read_tsv(path) == [{"id": "SAMN1", "host": "Homo sapiens", "isolation_source": "blood"}]


# --- two columns with one name --------------------------------------------------------------------


def test_a_repeated_column_name_is_an_error_not_last_one_wins(tmp_path: Path):
    """With two `host` columns the second replaced the first, so blood from a patient became Animal."""
    path = tmp_path / "m.tsv"
    path.write_text("id\thost\tisolation_source\thost\nA\tHomo sapiens\tblood\tBos taurus\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"column 'host' appears twice \(columns 2 and 4\)"):
        read_tsv(path)
    # Unnamed columns (trailing tabs in the header) are not names, and names differ by case as before.
    path.write_text("id\thost\tHost\t\t\nA\tx\ty\t\t\n", encoding="utf-8")
    assert read_tsv(path)[0]["id"] == "A"


# --- cells a spreadsheet would run ----------------------------------------------------------------


@pytest.mark.parametrize(
    "value, written",
    [
        ('=HYPERLINK("http://x","blood")', "'=HYPERLINK(\"http://x\",\"blood\")"),
        ("=1+1", "'=1+1"),
        ("+cmd|' /C calc'!A0", "'+cmd|' /C calc'!A0"),
        ("@SUM(1+1)", "'@SUM(1+1)"),
        ("-2+3+cmd|' /C calc'!A0", "'-2+3+cmd|' /C calc'!A0"),
        # Ordinary values that start with the same characters are left exactly as they are.
        ("-", "-"),
        ("-5", "-5"),
        ("-0.25", "-0.25"),
        ("+ve", "+ve"),
        ("+/-", "+/-"),
        ("@home", "@home"),
        ("blood (=venous)", "blood (=venous)"),
    ],
)
def test_formula_cells_are_written_as_text(tmp_path: Path, value, written):
    path = tmp_path / "out.tsv"
    write_tsv(path, [{"id": "A", "isolation_source": value}], ["id", "isolation_source"])
    assert read_tsv(path)[0]["isolation_source"] == written


# --- distance matrices ----------------------------------------------------------------------------


def test_a_one_triangle_or_negative_distance_matrix_is_rejected(tmp_path: Path):
    path = tmp_path / "d.tsv"
    path.write_text("snp-dists\tA\tB\tC\nA\t0\t5\t9\nB\t5\t0\t7\nC\t9\t7\t0\n", encoding="utf-8")
    names, matrix = read_distance_matrix(path)
    assert names == ["A", "B", "C"] and matrix[0, 2] == 9
    # Lower triangle only: read as it was, every pair would be 0 SNPs apart.
    path.write_text("snp-dists\tA\tB\tC\nA\t0\t0\t0\nB\t5\t0\t0\nC\t9\t7\t0\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"not symmetric \(row 1, column 2: 0, but 5 the other way round\)"):
        read_distance_matrix(path)
    path.write_text("snp-dists\tA\tB\nA\t0\t-3\nB\t-3\t0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="negative distance"):
        read_distance_matrix(path)
    path.write_text("A\tB\t-3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="negative distance between A and B"):
        read_distance_matrix(path)


def test_sequence_type_zero_in_a_typing_table_is_untyped(tmp_path: Path):
    """Kleborate writes ST 0 for alleles that match no ST. A PopPUNK or other cluster 0 is a real cluster."""
    kleborate = tmp_path / "kleborate.tsv"
    kleborate.write_text("strain\tST\nISO1\tST258\nISO2\t0\nISO3\t0\n", encoding="utf-8")
    assert read_lineages(kleborate) == {"ISO1": "ST258"}
    v3 = tmp_path / "kleborate3.tsv"
    v3.write_text("strain\tklebsiella_pneumo_complex__mlst__ST\nISO1\tST11\nISO2\t0\n", encoding="utf-8")
    assert read_lineages(v3) == {"ISO1": "ST11"}
    clusters = tmp_path / "clusters.csv"
    clusters.write_text("Taxon,Cluster\nISO1,0\nISO2,1\n", encoding="utf-8")
    assert read_lineages(clusters) == {"ISO1": "0", "ISO2": "1"}


# --- NCBI lookups ---------------------------------------------------------------------------------

XML = """<BioSampleSet><BioSample accession="{accession}">
  <Ids><Id db="BioSample" is_primary="1">{accession}</Id><Id db="SRA">SRS1</Id></Ids>
  <Description><Organism taxonomy_name="Klebsiella pneumoniae"/></Description>
  <Package display_name="Pathogen: environmental/food/other; version 1.0">Pathogen.env.1.0</Package>
  <Attributes><Attribute harmonized_name="isolation_source">feces</Attribute></Attributes>
</BioSample></BioSampleSet>"""


class _Fetcher:
    def __init__(self, payloads):
        self.payloads, self.urls = payloads, []

    def _lookup(self, url):
        self.urls.append(url)
        for key, payload in self.payloads.items():
            if key in url:
                if isinstance(payload, BaseException):
                    raise payload
                return payload
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": []}}
        raise AssertionError(f"unexpected URL {url}")

    def get_json(self, url, headers=None):
        return self._lookup(url)

    def get_text(self, url, headers=None):
        payload = self._lookup(url)
        return payload if isinstance(payload, str) else json.dumps(payload)


def test_the_biosample_package_is_carried_into_the_resolved_table():
    report = {"accession": "SAMN1", "package": "Pathogen.env.1.0", "attributes": [{"name": "isolation_source", "value": "feces"}]}
    assert attributes_from_biosample(report)["package"] == "Pathogen.env.1.0"
    assert parse_biosample_xml(XML.format(accession="SAMN1"))["package"] == "Pathogen.env.1.0"
    rows = resolve_ids(["SAMN1"], client=NCBIClient(fetcher=_Fetcher({"biosample/accession/SAMN1": {"reports": [report]}})),
                       allow_network=True)
    assert (rows[0]["package"], rows[0]["isolation_source"], rows[0]["resolve_status"]) == ("Pathogen.env.1.0", "feces", "ok")
    # A local table's attribute_package or epi_type fills the same column.
    local = resolve_ids(["SAMN2"], local_rows=[{"biosample": "SAMN2", "isolation_source": "feces", "epi_type": "environmental/other"}])
    assert local[0]["package"] == "environmental/other"


def test_the_xml_route_searches_the_accession_field_and_checks_what_comes_back():
    payloads = {
        "biosample/accession/SAMN5": {"reports": [{}]},
        "term=SAMN5%5Baccn%5D": {"esearchresult": {"idlist": ["5"]}},
        "efetch.fcgi": XML.format(accession="SAMN5"),
    }
    fetcher = _Fetcher(payloads)
    rows = resolve_ids(["SAMN5"], client=NCBIClient(fetcher=fetcher), allow_network=True)
    assert (rows[0]["resolve_status"], rows[0]["isolation_source"], rows[0]["package"]) == ("ok", "feces", "Pathogen.env.1.0")
    assert sum("esearch.fcgi" in url for url in fetcher.urls) == 1  # found by accession: no second search

    # Nothing under the accession field, and the all-fields search lands on another sample
    # that only mentions SAMN5: its metadata must not be attached to SAMN5.
    other = {
        "biosample/accession/SAMN5": {"reports": [{}]},
        "term=SAMN5&": {"esearchresult": {"idlist": ["9"]}},
        "efetch.fcgi": XML.format(accession="SAMN9"),
    }
    fetcher = _Fetcher(other)
    rows = resolve_ids(["SAMN5"], client=NCBIClient(fetcher=fetcher), allow_network=True)
    assert (rows[0]["resolve_status"], rows[0]["resolve_reason"], rows[0]["isolation_source"]) == ("review", "biosample_not_found", "")
    assert rows[0]["biosample"] == "SAMN5"
    assert sum("esearch.fcgi" in url for url in fetcher.urls) == 2  # accession field, then all fields

    # A sample listed under one of its other IDs (here its SRA sample) is that sample.
    assert parse_biosample_xml(XML.format(accession="SAMN5"))["ids"] == ["SAMN5", "SRS1"]


def test_a_reply_of_the_wrong_shape_sends_one_id_to_review_and_the_run_continues():
    """`attributes` as text made the run stop with a traceback and write nothing."""
    good = {"accession": "SAMN2", "attributes": [{"name": "isolation_source", "value": "blood"}]}
    fetcher = _Fetcher({
        "accession/SAMN1,SAMN2": {"reports": []},  # the batch answers nothing, so each is looked up alone
        "accession/SAMN1/": {"reports": [{"accession": "SAMN1", "attributes": "not a list of attributes"}]},
        "accession/SAMN2/": {"reports": [good]},
    })
    rows = resolve_ids(["SAMN1", "SAMN2"], client=NCBIClient(fetcher=fetcher), allow_network=True)
    assert [(r["query_id"], r["resolve_status"]) for r in rows] == [("SAMN1", "review"), ("SAMN2", "ok")]
    assert rows[0]["resolve_reason"].startswith("ncbi_error:unexpected reply (AttributeError")
    assert rows[1]["isolation_source"] == "blood"


def test_ctrl_c_keeps_what_was_already_looked_up(tmp_path: Path, monkeypatch, capsys):
    good = {"accession": "SAMN1", "attributes": [{"name": "isolation_source", "value": "blood"}]}
    payloads = {
        "accession/SAMN1,SAMN2,SAMN3": {"reports": []},
        "accession/SAMN1/": {"reports": [good]},
        "accession/SAMN2/": KeyboardInterrupt(),
    }
    rows = resolve_ids(["SAMN1", "SAMN2", "SAMN3"], client=NCBIClient(fetcher=_Fetcher(payloads)), allow_network=True)
    assert [(r["resolve_status"], r["resolve_reason"]) for r in rows] == [
        ("ok", "ncbi"), ("review", "interrupted"), ("review", "interrupted")]
    assert rows[0]["isolation_source"] == "blood"

    ids = tmp_path / "ids.txt"
    ids.write_text("SAMN1\nSAMN2\nSAMN3\n", encoding="utf-8")
    out = tmp_path / "resolved.tsv"
    monkeypatch.setattr(cli, "NCBIClient", lambda **kw: NCBIClient(fetcher=_Fetcher(payloads)))
    assert cli.main(["resolve", str(ids), "-o", str(out), "--email", "x@example.org"]) == 130
    assert "stopped before 2 of 3 IDs were looked up" in capsys.readouterr().err
    written = read_tsv(out)
    assert [r["resolve_reason"] for r in written] == ["ncbi", "interrupted", "interrupted"]
    assert [r["query_id"] for r in read_tsv(tmp_path / "resolved.review.tsv")] == ["SAMN2", "SAMN3"]


def test_requests_name_onehealth_as_the_client(monkeypatch):
    import io

    seen = []

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def urlopen(request, timeout):
        seen.append(request)
        return _Response(b'{"reports": []}\n')

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", urlopen)
    fetcher = ncbi.UrllibFetcher(delay_s=0)
    fetcher.get_json("https://api.example/x", headers={"Accept": "application/json", "api-key": "k"})
    fetcher.get_text("https://api.example/y")
    list(fetcher.iter_lines("https://ftp.example/z"))
    assert [r.get_header("User-agent") for r in seen] == [f"onehealth/{__version__}"] * 3
    assert seen[0].get_header("Api-key") == "k" and seen[1].get_header("Accept") == "application/json"


def test_api_key_can_come_from_the_environment(tmp_path: Path, monkeypatch):
    ids = tmp_path / "ids.txt"
    ids.write_text("not_an_id\n", encoding="utf-8")
    made = []
    real = cli.NCBIClient
    monkeypatch.setattr(cli, "NCBIClient", lambda **kw: made.append(kw) or real(**kw))
    monkeypatch.setenv("NCBI_API_KEY", " from-env ")
    assert cli.main(["resolve", str(ids), "-o", str(tmp_path / "a.tsv"), "--email", "x@example.org"]) == 0
    assert cli.main(["resolve", str(ids), "-o", str(tmp_path / "b.tsv"), "--email", "x@example.org", "--api-key", "typed"]) == 0
    monkeypatch.delenv("NCBI_API_KEY")
    assert cli.main(["resolve", str(ids), "-o", str(tmp_path / "c.tsv"), "--email", "x@example.org"]) == 0
    assert [kw["api_key"] for kw in made] == ["from-env", "typed", ""]


# --- which lexicon made a table -------------------------------------------------------------------

RULES = "rules:\n  - id: r\n    l1: Human\n    phrases: [blood]\n"


def test_lexicon_version_tells_an_edited_or_unversioned_file_from_the_seed(tmp_path: Path):
    """A copy of the seed with a rule changed was recorded as the seed's version; so was a file with no version."""
    seed_text = default_lexicon_path().read_text(encoding="utf-8")
    seed = load_lexicon()

    copy = tmp_path / "copy.yaml"
    copy.write_text(seed_text, encoding="utf-8")
    assert load_lexicon(copy).version == seed.version == "5"
    copy.write_bytes(seed_text.replace("\n", "\r\n").encode("utf-8"))  # the same file saved on Windows
    assert load_lexicon(copy).version == "5"

    edited = tmp_path / "edited.yaml"
    edited.write_text(seed_text.replace("      - homo sapien\n", "      - martian\n"), encoding="utf-8")
    version = load_lexicon(edited).version
    assert version.startswith("5-modified-") and len(version) == len("5-modified-") + 8
    assert load_lexicon(edited).version == version  # the same file, the same name

    unversioned = tmp_path / "mine.yaml"
    unversioned.write_text(RULES, encoding="utf-8")
    assert load_lexicon(unversioned).version.startswith("custom-")
    unversioned.write_text('version: ""\n' + RULES, encoding="utf-8")
    assert load_lexicon(unversioned).version.startswith("custom-")

    named = tmp_path / "named.yaml"
    named.write_text('version: "mylab-2"\n' + RULES, encoding="utf-8")
    assert load_lexicon(named).version == "mylab-2"
    named.write_text('version: "4"\n' + RULES, encoding="utf-8")  # not the shipped version: taken at its word
    assert load_lexicon(named).version == "4"


def test_classify_reports_and_records_a_modified_lexicon(tmp_path: Path, capsys):
    table = tmp_path / "m.tsv"
    table.write_text("id\thost\nA\tmartian\n", encoding="utf-8")
    edited = tmp_path / "edited.yaml"
    edited.write_text(default_lexicon_path().read_text(encoding="utf-8").replace("      - homo sapien\n", "      - martian\n"),
                      encoding="utf-8")
    out = tmp_path / "out.tsv"
    assert cli.main(["classify", str(table), "-o", str(out), "--lexicon", str(edited)]) == 0
    row = read_tsv(out)[0]
    assert row["l1"] == "Human" and row["lexicon_version"].startswith("5-modified-")
    assert f"lexicon {row['lexicon_version']} (seed v5)" in capsys.readouterr().err
