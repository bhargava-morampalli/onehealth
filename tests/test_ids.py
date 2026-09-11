from onehealth.ids import parse_id, parse_ids


def test_parse_mixed_public_ids():
    parsed = parse_ids(
        [
            "GCA_000005845.2",
            "GCF_000005845.2",
            "SAMN12345678",
            "SAMEA1234567",
            "SRR1234567",
            "PDT000123456.1",
            "not-an-id",
        ]
    )
    assert [item.id_type for item in parsed] == [
        "assembly",
        "assembly",
        "biosample",
        "biosample",
        "sra",
        "pdt",
        "unknown",
    ]


def test_parse_assembly_without_version():
    assert parse_id("gca_000005845").normalized == "GCA_000005845"
    assert parse_id("gca_000005845").id_type == "assembly"


def test_blank_is_unknown():
    assert parse_id("").id_type == "unknown"
