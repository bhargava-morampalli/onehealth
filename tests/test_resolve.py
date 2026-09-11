from onehealth.ids import parse_id
from onehealth.ncbi import attributes_from_biosample, biosample_from_genome
from onehealth.resolve import resolve_ids


class FakeFetcher:
    def __init__(self, payloads: dict[str, dict]) -> None:
        self.payloads = payloads

    def get_json(self, url: str, headers=None):
        for key, payload in self.payloads.items():
            if key in url:
                return payload
        raise AssertionError(f"unexpected URL {url}")

    def get_text(self, url: str, headers=None):
        return ""


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


def test_biosample_from_genome_report():
    report = {"accession": "GCA_1.1", "assembly_info": {"biosample": {"accession": "SAMN2"}}}
    assert biosample_from_genome(report) == "SAMN2"
