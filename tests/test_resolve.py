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


def _fake_urlopen(responses: dict[str, list[bytes]], seen: list):
    """urlopen stand-in: each URL fragment answers with its next body (the last one repeats)."""
    import io

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def urlopen(request, timeout):
        seen.append(request)
        for fragment, bodies in responses.items():
            if fragment in request.full_url:
                return _Response(bodies.pop(0) if len(bodies) > 1 else bodies[0])
        raise AssertionError(f"unexpected URL {request.full_url}")

    return urlopen


def test_garbled_ncbi_answer_sends_one_id_to_review_and_the_run_continues(monkeypatch):
    from onehealth import ncbi

    good = b'{"reports": [{"accession": "SAMN1", "attributes": [{"name": "isolation_source", "value": "blood"}]}]}'
    maintenance = b"<html><body>Service temporarily unavailable</body></html>"
    seen: list = []
    monkeypatch.setattr(ncbi.urllib.request, "urlopen", _fake_urlopen({"SAMN1/": [good], "SAMN2/": [maintenance]}, seen))
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    client = NCBIClient(fetcher=ncbi.UrllibFetcher(delay_s=0, retries=2))
    rows = resolve_ids(["SAMN2", "SAMN1"], client=client, allow_network=True)
    assert [(r["query_id"], r["resolve_status"]) for r in rows] == [("SAMN2", "review"), ("SAMN1", "ok")]
    assert rows[0]["resolve_reason"].startswith("ncbi_error:Invalid JSON from https://api.ncbi.nlm.nih.gov/")
    assert rows[1]["isolation_source"] == "blood"
    assert sum("SAMN2/" in r.full_url for r in seen) == 3  # first try plus two retries


def test_fetcher_retries_a_garbled_body_then_succeeds(monkeypatch):
    from onehealth import ncbi

    seen: list = []
    monkeypatch.setattr(ncbi.urllib.request, "urlopen", _fake_urlopen({"x": [b"<html>", b'{"ok": 1}']}, seen))
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    assert ncbi.UrllibFetcher(delay_s=0).get_json("https://example.org/x") == {"ok": 1}
    assert len(seen) == 2


def test_fetcher_rejects_json_that_is_not_an_object(monkeypatch):
    import pytest

    from onehealth import ncbi

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", _fake_urlopen({"x": [b"[1, 2]"]}, []))
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    with pytest.raises(ncbi.HttpError, match="Expected a JSON object"):
        ncbi.UrllibFetcher(delay_s=0, retries=0).get_json("https://example.org/x")


def test_custom_fetcher_json_error_does_not_stop_the_run():
    class Broken(FakeFetcher):
        def get_json(self, url, headers=None):
            if "SAMN2" in url:
                return json.loads("<html>")
            return super().get_json(url, headers)

    good = {"reports": [{"accession": "SAMN1", "attributes": [{"name": "isolation_source", "value": "urine"}]}]}
    rows = resolve_ids(["SAMN1", "SAMN2"], client=NCBIClient(fetcher=Broken({"SAMN1/": good})), allow_network=True)
    assert [(r["resolve_status"], r["resolve_reason"][:11]) for r in rows] == [("ok", "ncbi"), ("review", "ncbi_error:")]


def test_rejected_api_key_is_not_written_to_the_table():
    from onehealth.ncbi import HttpError

    class Rejecting(FakeFetcher):
        def get_json(self, url, headers=None):
            raise HttpError('HTTP 400 for https://api.ncbi.nlm.nih.gov/datasets/v2/x: {"message":"API key invalid","api-key":"MYKEY123"}')

    client = NCBIClient(api_key="MYKEY123", fetcher=Rejecting({}))
    row = resolve_ids(["SAMN1"], client=client, allow_network=True)[0]
    assert row["resolve_status"] == "review"
    assert "MYKEY123" not in row["resolve_reason"] and "***" in row["resolve_reason"]


def test_api_key_raises_the_request_rate_and_reaches_datasets_as_a_header():
    from onehealth.ncbi import UrllibFetcher

    assert NCBIClient().fetcher.delay_s == 0.35  # 3 requests/s without a key
    assert NCBIClient(api_key="K").fetcher.delay_s == 0.11  # 10 requests/s with one
    assert NCBIClient(api_key="K", fetcher=UrllibFetcher(delay_s=1.0)).fetcher.delay_s == 1.0

    class Recording(FakeFetcher):
        headers: list = []

        def get_json(self, url, headers=None):
            Recording.headers.append(headers or {})
            return super().get_json(url, headers)

    fetcher = Recording({"biosample_report": {"reports": [{"accession": "SAMN1"}]}})
    NCBIClient(api_key="K", fetcher=fetcher).biosample_report("SAMN1")
    assert Recording.headers[-1].get("api-key") == "K"
    assert all("K" not in url.split("?")[-1] for url in fetcher.urls)


# --- batched lookups --------------------------------------------------------------------------


class DatasetsStub:
    """In-memory NCBI that answers like the live API did when probed (November 2026).

    BioSample batches return one report per accession asked, an empty one for an unknown
    record. Genome batches come back in no particular order, 20 per page unless page_size
    says otherwise, answer an unversioned accession with its current version, and return
    nothing for an outdated version. E-utilities finds nothing, so the XML fallback is a miss.
    """

    def __init__(self, biosamples: dict[str, dict], genomes: dict[str, dict], max_page: int = 3, fail: tuple = ()):
        self.biosamples, self.genomes, self.max_page, self.fail = biosamples, genomes, max_page, fail
        self.urls: list[str] = []

    def get_json(self, url, headers=None):
        import urllib.parse

        from onehealth.ncbi import HttpError

        self.urls.append(url)
        path, _, query = url.partition("?")
        params = dict(urllib.parse.parse_qsl(query))
        if any(bad in path for bad in self.fail):
            raise HttpError(f"HTTP 400 for {path}: bad request")
        if "/biosample/accession/" in path:
            asked = urllib.parse.unquote(path.split("/biosample/accession/")[1].split("/")[0]).split(",")
            return {"reports": [self.biosamples.get(a.upper(), {}) for a in asked]}
        if "/genome/accession/" in path:
            asked = urllib.parse.unquote(path.split("/genome/accession/")[1].split("/")[0]).split(",")
            hits: dict[str, dict] = {}
            for a in (x.upper() for x in asked):
                if "." in a:
                    report = self.genomes.get(a)
                else:
                    versions = sorted(k for k in self.genomes if k.split(".")[0] == a)
                    report = self.genomes[versions[-1]] if versions else None
                if report is not None and report.get("current", True):
                    hits[report["accession"]] = report
            found = list(reversed(hits.values()))  # not the order asked
            size = min(int(params.get("page_size", 20)), self.max_page)
            start = int(params.get("page_token", 0))
            page = {"reports": found[start : start + size], "total_count": len(found)}
            if start + size < len(found):
                page["next_page_token"] = str(start + size)
            return page
        if "esearch.fcgi" in path:
            return {"esearchresult": {"idlist": []}}
        raise AssertionError(f"unexpected URL {url}")

    def get_text(self, url, headers=None):
        raise AssertionError(f"unexpected URL {url}")


def _biosample(accession: str, source: str) -> dict:
    return {"accession": accession, "attributes": [{"name": "isolation_source", "value": source}]}


def _genome(accession: str, biosample: str, current: bool = True) -> dict:
    return {"accession": accession, "assembly_info": {"biosample": {"accession": biosample}}, "current": current}


def _stub(**kwargs) -> DatasetsStub:
    biosamples = {f"SAMN{i}": _biosample(f"SAMN{i}", f"source {i}") for i in range(1, 251)}
    biosamples["SAMN900"] = {}  # Datasets has an empty record for it
    genomes = {
        "GCA_000000001.2": _genome("GCA_000000001.2", "SAMN1"),
        "GCA_000000001.1": _genome("GCA_000000001.1", "SAMN1", current=False),  # outdated version
        "GCF_000000001.1": _genome("GCF_000000001.1", "SAMN1"),
        **{f"GCA_00000000{i}.1": _genome(f"GCA_00000000{i}.1", f"SAMN{i}") for i in range(2, 10)},
        "GCA_000000010.1": {"accession": "GCA_000000010.1", "assembly_info": {}},  # no BioSample
    }
    return DatasetsStub(biosamples, genomes, **kwargs)


MIXED_IDS = (
    [f"SAMN{i}" for i in (3, 1, 250, 120)]
    + ["GCA_000000001.2", "GCA_000000001", "gcf_000000001.1", "GCA_000000001.1", "GCA_000000010.1", "GCA_999999999.1"]
    + [f"GCA_00000000{i}.1" for i in range(2, 10)]
    + ["SAMN900", "SAMN99999", "SAMN3", "GCA_000000004.1", "not-an-id", "SRP000001"]
)


def _one_at_a_time(client: NCBIClient) -> NCBIClient:
    """The same client with batching switched off: every ID takes the single-ID path."""
    client.biosample_reports = lambda accessions: {}
    client.genome_reports = lambda accessions: {}
    return client


def test_batched_lookups_give_exactly_the_rows_of_one_at_a_time_lookups():
    batched_stub, single_stub = _stub(), _stub()
    batched = resolve_ids(MIXED_IDS, client=NCBIClient(fetcher=batched_stub), allow_network=True)
    single = resolve_ids(MIXED_IDS, client=_one_at_a_time(NCBIClient(fetcher=single_stub)), allow_network=True)
    assert batched == single
    assert [row["query_id"] for row in batched] == MIXED_IDS
    reasons = {row["query_id"]: (row["resolve_status"], row["resolve_reason"], row["assembly"], row["biosample"]) for row in batched}
    assert reasons["GCA_000000001"] == ("ok", "ncbi_assembly", "GCA_000000001.2", "SAMN1")  # unversioned: current
    assert reasons["GCA_000000001.1"] == ("ok", "ncbi_assembly", "GCA_000000001.2", "SAMN1")  # outdated: via the stem
    assert reasons["gcf_000000001.1"][:2] == ("ok", "ncbi_assembly")
    assert reasons["GCA_000000010.1"][:2] == ("review", "assembly_no_biosample")
    assert reasons["GCA_999999999.1"][:2] == ("review", "assembly_not_found")
    assert reasons["SAMN900"][:2] == reasons["SAMN99999"][:2] == ("review", "biosample_not_found")
    assert len(batched_stub.urls) < len(single_stub.urls) / 2


def test_genome_batches_follow_every_page():
    stub = _stub(max_page=3)  # 12 assemblies come back 3 to a page
    rows = resolve_ids([f"GCA_00000000{i}.1" for i in range(2, 10)], client=NCBIClient(fetcher=stub), allow_network=True)
    assert all(row["resolve_status"] == "ok" for row in rows)
    single_genome_urls = [u for u in stub.urls if "/genome/" in u and "," not in u]
    assert single_genome_urls == []  # every assembly came from the paged batch


def test_many_biosamples_take_one_request_per_hundred():
    stub = _stub()
    ids = [f"SAMN{i}" for i in range(1, 251)]
    rows = resolve_ids(ids, client=NCBIClient(fetcher=stub), allow_network=True)
    assert [row["isolation_source"] for row in rows] == [f"source {i}" for i in range(1, 251)]
    assert len(stub.urls) == 3


def test_a_rejected_batch_falls_back_to_one_at_a_time():
    # Only the comma-joined batch URLs fail; the same accessions asked one by one succeed.
    stub = _stub(fail=(",",))
    rows = resolve_ids(MIXED_IDS, client=NCBIClient(fetcher=stub), allow_network=True)
    assert rows == resolve_ids(MIXED_IDS, client=_one_at_a_time(NCBIClient(fetcher=_stub())), allow_network=True)


def test_a_value_that_is_not_an_accession_is_never_put_in_a_batch():
    from onehealth.ncbi import _chunks

    assert list(_chunks(["SAMN1", "12345", " samn2 ", "SAMN1", "GCA_000000001.1"], "biosample")) == [["SAMN1", "SAMN2"]]
    assert list(_chunks(["GCA_000000001.1", "SAMN1", "gca_000000002"], "assembly")) == [["GCA_000000001.1", "GCA_000000002"]]


def test_timeouts_and_dropped_connections_are_retried(monkeypatch):
    import io

    from onehealth import ncbi

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    failures = [TimeoutError("The read operation timed out"), ConnectionResetError(104, "Connection reset by peer")]

    def urlopen(request, timeout):
        if failures:
            raise failures.pop(0)
        return _Response(b'{"ok": 1}')

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    assert ncbi.UrllibFetcher(delay_s=0).get_json("https://example.org/x") == {"ok": 1}


def test_a_connection_that_keeps_timing_out_sends_only_that_id_to_review(monkeypatch):
    import io

    from onehealth import ncbi

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    good = b'{"reports": [{"accession": "SAMN1", "attributes": [{"name": "isolation_source", "value": "blood"}]}]}'

    def urlopen(request, timeout):
        if "SAMN2" in request.full_url:
            raise TimeoutError("The read operation timed out")
        return _Response(good)

    monkeypatch.setattr(ncbi.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    rows = resolve_ids(["SAMN2", "SAMN1"], client=NCBIClient(fetcher=ncbi.UrllibFetcher(delay_s=0, retries=1)), allow_network=True)
    assert [(r["resolve_status"], r["resolve_reason"][:11]) for r in rows] == [("review", "ncbi_error:"), ("ok", "ncbi")]
