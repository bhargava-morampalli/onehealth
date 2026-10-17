import json

from onehealth.ids import parse_id
from onehealth.io_tsv import get_field, metadata_view
from onehealth.ncbi import NCBIClient, attributes_from_biosample, biosample_from_genome, parse_biosample_xml
from onehealth.resolve import resolve_ids

SAMN63084674_XML = """<?xml version="1.0" ?>
<BioSampleSet><BioSample accession="SAMN63084674">
  <Description>
    <Organism taxonomy_name="Klebsiella pneumoniae">
      <OrganismName>Klebsiella pneumoniae</OrganismName>
    </Organism>
  </Description>
  <Attributes>
    <Attribute attribute_name="host" harmonized_name="host">Homo sapiens</Attribute>
    <Attribute attribute_name="isolation_source" harmonized_name="isolation_source">urine</Attribute>
    <Attribute attribute_name="sample_type" harmonized_name="sample_type">whole organism</Attribute>
    <Attribute attribute_name="geo_loc_name" harmonized_name="geo_loc_name">USA</Attribute>
  </Attributes>
</BioSample></BioSampleSet>
"""


class FakeFetcher:
    def __init__(self, payloads: dict[str, object]) -> None:
        self.payloads = payloads
        self.urls: list[str] = []

    def _lookup(self, url: str):
        self.urls.append(url)
        for key, payload in self.payloads.items():
            if key in url:
                return payload
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": []}}
        if "esummary.fcgi" in url:
            return {"result": {}}
        raise AssertionError(f"unexpected URL {url}")

    def get_json(self, url: str, headers=None):
        payload = self._lookup(url)
        if isinstance(payload, str):
            raise AssertionError(f"JSON expected, got text for {url}")
        return payload

    def get_text(self, url: str, headers=None):
        payload = self._lookup(url)
        if isinstance(payload, str):
            return payload
        return json.dumps(payload)


def _client(payloads: dict[str, object]) -> NCBIClient:
    return NCBIClient(fetcher=FakeFetcher(payloads))


def test_resolve_joins_local_metadata_without_network():
    rows = resolve_ids(
        ["SAMN1", "GCA_000000001.1", "mystery"],
        local_rows=[
            {
                "id": "SAMN1",
                "biosample": "SAMN1",
                "host": "Homo sapiens",
                "isolation_source": "wound drainage",
            },
            {
                "assembly": "GCA_000000001.1",
                "biosample": "SAMN9",
                "isolation_source": "hospital drain",
            },
        ],
        allow_network=False,
    )
    assert rows[0]["resolve_status"] == "ok"
    assert rows[0]["isolation_source"] == "wound drainage"
    assert rows[1]["resolve_status"] == "ok"
    assert rows[1]["biosample"] == "SAMN9"
    assert rows[2]["resolve_status"] == "review"
    assert rows[2]["resolve_reason"] == "unknown_id"


def test_missing_local_metadata_without_network_is_review():
    rows = resolve_ids(["SAMN99"], local_rows=[], allow_network=False)
    assert rows[0]["resolve_status"] == "review"
    assert rows[0]["resolve_reason"] == "missing_metadata"
    assert parse_id("SAMN99").id_type == "biosample"


def test_blank_isolation_source_does_not_hide_source():
    assert get_field({"isolation_source": "", "source": "hospital sink"}, ("isolation_source", "source")) == "hospital sink"
    view = metadata_view({"isolation_source": "", "source": "hospital sink"})
    assert view["isolation_source"] == "hospital sink"
    rows = resolve_ids(
        ["SAMN1"],
        local_rows=[{"biosample": "SAMN1", "isolation_source": "", "source": "hospital sink"}],
        allow_network=False,
    )
    assert rows[0]["resolve_status"] == "ok"
    assert rows[0]["isolation_source"] == "hospital sink"


def test_local_gca_joins_gcf_and_versionless():
    local = [{"assembly": "GCF_000005845.2", "biosample": "SAMN02604091", "isolation_source": "lab"}]
    rows = resolve_ids(
        ["GCA_000005845.2", "GCA_000005845", "gca_000005845.2"],
        local_rows=local,
        allow_network=False,
    )
    assert [row["resolve_status"] for row in rows] == ["ok", "ok", "ok"]
    assert all(row["biosample"] == "SAMN02604091" for row in rows)


def test_duplicate_ids_are_not_collapsed():
    rows = resolve_ids(
        ["SAMN1", "SAMN1"],
        local_rows=[{"id": "SAMN1", "host": "Homo sapiens", "isolation_source": "urine"}],
        allow_network=False,
    )
    assert len(rows) == 2
    assert [row["query_id"] for row in rows] == ["SAMN1", "SAMN1"]
    assert [row["resolve_status"] for row in rows] == ["ok", "ok"]


def test_srp_is_refused_without_network():
    rows = resolve_ids(["SRP020237"], local_rows=[], allow_network=True, client=_client({}))
    assert rows[0]["id_type"] == "sra"
    assert rows[0]["resolve_status"] == "review"
    assert rows[0]["resolve_reason"] == "sra_study_ambiguous"
    assert rows[0]["biosample"] == ""


def test_empty_datasets_report_falls_back_to_eutils_xml():
    client = _client(
        {
            "biosample/accession/SAMN63084674": {"reports": [{}]},
            "db=biosample&term=SAMN63084674": {"esearchresult": {"idlist": ["63084674"]}},
            "efetch.fcgi": SAMN63084674_XML,
        }
    )
    rows = resolve_ids(["SAMN63084674"], client=client, allow_network=True)
    assert rows[0]["biosample"] == "SAMN63084674"
    assert rows[0]["host"] == "Homo sapiens"
    assert rows[0]["isolation_source"] == "urine"
    assert rows[0]["sample_type"] == "whole organism"
    assert rows[0]["resolve_status"] == "ok"


def test_empty_datasets_report_does_not_wipe_biosample_when_xml_missing():
    client = _client(
        {
            "biosample/accession/SAMN63084674": {"reports": [{}]},
            "db=biosample&term=SAMN63084674": {"esearchresult": {"idlist": []}},
        }
    )
    rows = resolve_ids(["SAMN63084674"], client=client, allow_network=True)
    assert rows[0]["biosample"] == "SAMN63084674"
    assert rows[0]["resolve_status"] == "review"
    assert rows[0]["resolve_reason"] == "biosample_not_found"


def test_assembly_with_empty_source_is_ok():
    client = _client(
        {
            "genome/accession/GCA_000016305.1": {
                "reports": [
                    {
                        "accession": "GCA_000016305.1",
                        "assembly_info": {"biosample": {"accession": "SAMN02603941"}},
                    }
                ]
            },
            "biosample/accession/SAMN02603941": {
                "reports": [
                    {
                        "accession": "SAMN02603941",
                        "description": {
                            "organism": {"organism_name": "Klebsiella pneumoniae subsp. pneumoniae MGH 78578"}
                        },
                        "attributes": [{"name": "strain", "value": "ATCC 700721; MGH 78578"}],
                        "strain": "ATCC 700721; MGH 78578",
                    }
                ]
            },
        }
    )
    rows = resolve_ids(["GCA_000016305.1"], client=client, allow_network=True)
    assert rows[0]["resolve_status"] == "ok"
    assert rows[0]["resolve_reason"] == "ncbi_assembly"
    assert rows[0]["biosample"] == "SAMN02603941"
    assert rows[0]["host"] == ""
    assert rows[0]["isolation_source"] == ""
    assert rows[0]["strain"] == "ATCC 700721; MGH 78578"


def test_organism_is_not_copied_to_host():
    client = _client(
        {
            "biosample/accession/SAMN00630374": {
                "reports": [
                    {
                        "accession": "SAMN00630374",
                        "description": {"organism": {"organism_name": "Homo sapiens"}},
                        "attributes": [{"name": "cell_line", "value": "DB"}],
                    }
                ]
            }
        }
    )
    rows = resolve_ids(["SAMN00630374"], client=client, allow_network=True)
    assert rows[0]["resolve_status"] == "ok"
    assert rows[0]["organism"] == "Homo sapiens"
    assert rows[0]["host"] == ""


def test_missing_host_token_is_blanked():
    report = {
        "accession": "SAMN10527295",
        "description": {"organism": {"organism_name": "Klebsiella pneumoniae"}},
        "attributes": [
            {"name": "host", "value": "missing"},
            {"name": "isolation_source", "value": "urine"},
            {"name": "sample_type", "value": "pure culture"},
        ],
        "host": "missing",
        "isolation_source": "urine",
    }
    attrs = attributes_from_biosample(report)
    assert attrs["host"] == ""
    assert attrs["isolation_source"] == "urine"
    assert attrs["sample_type"] == "pure culture"


def test_pdt_resolves_via_assembly_name():
    client = _client(
        {
            "PDT000806148.1%5BAssembly+Name%5D": {"esearchresult": {"idlist": ["7853101"]}},
            "db=assembly&id=7853101": {
                "result": {"7853101": {"assemblyname": "PDT000806148.1", "biosampleaccn": "SAMN09405588"}}
            },
            "biosample/accession/SAMN09405588": {
                "reports": [
                    {
                        "accession": "SAMN09405588",
                        "host": "Homo sapiens",
                        "isolation_source": "blood",
                        "description": {"organism": {"organism_name": "Escherichia coli"}},
                    }
                ]
            },
        }
    )
    rows = resolve_ids(["PDT000806148.1"], client=client, allow_network=True)
    assert rows[0]["resolve_status"] == "ok"
    assert rows[0]["resolve_reason"] == "ncbi_pdt"
    assert rows[0]["biosample"] == "SAMN09405588"
    assert rows[0]["isolation_source"] == "blood"
    assert not any("db=biosample&term=PDT" in url for url in client.fetcher.urls)


def test_pdt_unresolved_when_assembly_name_missing():
    rows = resolve_ids(["PDT000000001.1"], client=_client({}), allow_network=True)
    assert rows[0]["resolve_status"] == "review"
    assert rows[0]["resolve_reason"] == "pdt_unresolved"
    assert rows[0]["biosample"] == ""


def test_parse_biosample_report_attributes():
    report = {
        "accession": "SAMN1",
        "host": "Homo sapiens",
        "isolation_source": "blood",
        "description": {"organism": {"organism_name": "Klebsiella pneumoniae"}},
        "attributes": [{"name": "geo_loc_name", "value": "USA"}],
    }
    attrs = attributes_from_biosample(report)
    assert attrs["biosample"] == "SAMN1"
    assert attrs["host"] == "Homo sapiens"
    assert attrs["organism"] == "Klebsiella pneumoniae"
    assert attrs["geo_loc_name"] == "USA"


def test_parse_biosample_xml_maps_harmonized_attributes():
    report = parse_biosample_xml(SAMN63084674_XML)
    assert report is not None
    attrs = attributes_from_biosample(report)
    assert attrs["biosample"] == "SAMN63084674"
    assert attrs["host"] == "Homo sapiens"
    assert attrs["isolation_source"] == "urine"
    assert attrs["sample_type"] == "whole organism"


def test_biosample_from_genome_report():
    report = {"accession": "GCA_1.1", "assembly_info": {"biosample": {"accession": "SAMN2"}}}
    assert biosample_from_genome(report) == "SAMN2"


def test_versioned_genome_empty_report_retries_stem():
    client = _client(
        {
            "genome/accession/GCA_000005845.1": {},
            "genome/accession/GCA_000005845/": {
                "reports": [
                    {
                        "accession": "GCA_000005845.2",
                        "assembly_info": {"biosample": {"accession": "SAMN02604091"}},
                    }
                ]
            },
            "biosample/accession/SAMN02604091": {
                "reports": [
                    {
                        "accession": "SAMN02604091",
                        "description": {"organism": {"organism_name": "Escherichia coli"}},
                        "attributes": [{"name": "strain", "value": "K-12"}],
                    }
                ]
            },
        }
    )
    rows = resolve_ids(["GCA_000005845.1"], client=client, allow_network=True)
    assert rows[0]["resolve_status"] == "ok"
    assert rows[0]["biosample"] == "SAMN02604091"
    assert rows[0]["assembly"] == "GCA_000005845.2"


def test_local_exact_accession_beats_versionless_match():
    local = [
        {"assembly": "GCA_000005845.1", "isolation_source": "urine"},
        {"assembly": "GCA_000005845.2", "isolation_source": "river water"},
    ]
    rows = resolve_ids(["GCA_000005845.1", "GCA_000005845.2"], local_rows=local)
    assert [row["isolation_source"] for row in rows] == ["urine", "river water"]


def test_conflicting_local_rows_for_same_key_go_to_review():
    local = [
        {"assembly": "GCA_000005845.2", "isolation_source": "urine"},
        {"assembly": "GCF_000005845.2", "isolation_source": "river water"},
    ]
    rows = resolve_ids(["GCA_000005845"], local_rows=local)
    assert rows[0]["resolve_status"] == "review"
    assert rows[0]["resolve_reason"] == "local_metadata_conflict"


def test_identical_local_rows_are_not_a_conflict():
    local = [
        {"id": "SAMN1", "biosample": "SAMN1", "isolation_source": "urine"},
        {"id": "PDT1.1", "biosample": "SAMN1", "isolation_source": "urine"},
    ]
    assert resolve_ids(["SAMN1"], local_rows=local)[0]["resolve_status"] == "ok"


def test_unversioned_pdt_tries_assembly_name_wildcard():
    client = _client(
        {
            "PDT000002407%2A": {"esearchresult": {"idlist": ["9"]}},
            "db=assembly&id=9": {"result": {"9": {"assemblyname": "PDT000002407.4", "biosampleaccn": "SAMN02344983"}}},
        }
    )
    assert client.pdt_to_biosample("PDT000002407") == ["SAMN02344983"]


def test_partial_unknown_text_is_kept_but_null_tokens_are_blanked():
    attrs = attributes_from_biosample(
        {"accession": "SAMN1", "host": "not collected", "isolation_source": "unknown river"}
    )
    assert attrs["host"] == ""
    assert attrs["isolation_source"] == "unknown river"


def test_ncbi_source_type_preferred_over_material_sample_type():
    attrs = attributes_from_biosample(
        {
            "accession": "SAMN1",
            "attributes": [
                {"name": "sample_type", "value": "whole organism"},
                {"name": "source_type", "value": "animal"},
            ],
        }
    )
    assert attrs["sample_type"] == "animal"


def test_fetcher_retries_rate_limit_then_succeeds(monkeypatch):
    import io
    import urllib.error

    from onehealth import ncbi

    calls = []

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def fake_urlopen(request, timeout):
        calls.append(request.full_url)
        if len(calls) == 1:
            raise urllib.error.HTTPError(request.full_url, 429, "Too Many Requests", {}, io.BytesIO(b""))
        return _Response(b'{"ok": 1}')

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    assert ncbi.UrllibFetcher(delay_s=0).get_json("https://example.org/x") == {"ok": 1}
    assert len(calls) == 2


def test_fetcher_does_not_retry_client_errors(monkeypatch):
    import io
    import urllib.error

    import pytest

    from onehealth import ncbi

    calls = []

    def fake_urlopen(request, timeout):
        calls.append(1)
        raise urllib.error.HTTPError(request.full_url, 400, "Bad Request", {}, io.BytesIO(b""))

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(ncbi.HttpError):
        ncbi.UrllibFetcher(delay_s=0).get_text("https://example.org/x")
    assert len(calls) == 1


def test_http_error_text_does_not_leak_api_key(monkeypatch):
    import io
    import urllib.error

    import pytest

    from onehealth import ncbi

    def fake_urlopen(request, timeout):
        raise urllib.error.HTTPError(request.full_url, 400, "Bad", {}, io.BytesIO(b"bad\trequest\n"))

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(ncbi.HttpError) as exc:
        ncbi.UrllibFetcher(delay_s=0).get_text("https://eutils.example/esearch.fcgi?term=x&api_key=SECRET")
    assert "SECRET" not in str(exc.value)
    assert "\t" not in str(exc.value) and "\n" not in str(exc.value)


def test_local_extra_columns_are_carried_through():
    local = [{"biosample": "SAMN1", "isolation_source": "urine", "ST": "ST258", "host_subject_id": "P1"}]
    row = resolve_ids(["SAMN1"], local_rows=local)[0]
    assert (row["ST"], row["host_subject_id"]) == ("ST258", "P1")
