from onehealth.summarize import SOURCE_EVIDENCE_CORE, summarize


def test_summarize_skips_overlap_without_cluster_or_plasmid_columns():
    rows = [
        {"l1": "Human", "paper_rollup": "Human", "healthcare_environment": "no", "classify_status": "ok", "classify_reason": "blood"},
        {"l1": "Environmental", "paper_rollup": "Environmental", "healthcare_environment": "yes", "classify_status": "ok", "classify_reason": "hospital-drain"},
    ]
    summary = summarize(rows)
    sections = {row["section"] for row in summary}
    assert "l1" in sections
    assert "paper_rollup" in sections
    assert "cluster_overlap" not in sections
    assert "plasmid_overlap" not in sections


def test_summarize_cluster_overlap_only_when_column_present():
    rows = [
        {"l1": "Human", "cluster_id": "A", "paper_rollup": "Human", "healthcare_environment": "no", "classify_status": "ok", "classify_reason": "x"},
        {"l1": "Environmental", "cluster_id": "A", "paper_rollup": "Environmental", "healthcare_environment": "yes", "classify_status": "ok", "classify_reason": "y"},
        {"l1": "Animal", "cluster_id": "B", "paper_rollup": "Animal", "healthcare_environment": "no", "classify_status": "ok", "classify_reason": "z"},
    ]
    summary = summarize(rows)
    overlap = [row for row in summary if row["section"] == "cluster_overlap"]
    mixed = {row["key"]: row["mixed"] for row in overlap}
    assert mixed["A"] == "yes"
    assert mixed["B"] == "no"


def test_summarize_plasmid_overlap_only_when_column_present():
    rows = [
        {"l1": "Human", "plasmid_id": "IncF", "paper_rollup": "Human", "healthcare_environment": "no", "classify_status": "ok", "classify_reason": "x"},
        {"l1": "Food", "plasmid_id": "IncF", "paper_rollup": "Environmental", "healthcare_environment": "no", "classify_status": "ok", "classify_reason": "y"},
    ]
    summary = summarize(rows)
    overlap = [row for row in summary if row["section"] == "plasmid_overlap"]
    assert overlap[0]["mixed"] == "yes"
    assert "Human" in overlap[0]["l1_set"]
    assert "Food" in overlap[0]["l1_set"]


def test_summarize_source_evidence_splits_ok_empty_from_unmatched_and_conflict():
    rows = [
        {
            "l1": "Human",
            "isolation_source": "blood",
            "classify_status": "ok",
            "classify_reason": "human-specimens",
            "paper_rollup": "Human",
            "healthcare_environment": "no",
        },
        {
            "l1": "Animal",
            "isolation_source": "",
            "host": "Bos taurus",
            "classify_status": "ok",
            "classify_reason": "animal-hosts",
            "paper_rollup": "Animal",
            "healthcare_environment": "no",
        },
        {
            "l1": "Unclassified",
            "isolation_source": "",
            "resolve_status": "ok",
            "classify_status": "review",
            "classify_reason": "missing_source",
            "paper_rollup": "Unclassified",
            "healthcare_environment": "no",
        },
        {
            "l1": "Unclassified",
            "isolation_source": "sink",
            "classify_status": "review",
            "classify_reason": "unmatched",
            "paper_rollup": "Unclassified",
            "healthcare_environment": "no",
        },
        {
            "l1": "Unclassified",
            "isolation_source": "urine",
            "host": "chicken",
            "classify_status": "review",
            "classify_reason": "host_source_conflict",
            "paper_rollup": "Unclassified",
            "healthcare_environment": "no",
        },
    ]
    summary = summarize(rows)
    evidence = {row["key"]: row["count"] for row in summary if row["section"] == "source_evidence"}
    for key in SOURCE_EVIDENCE_CORE:
        assert key in evidence
    assert evidence["ok_with_source"] == "1"
    # Only classify ok rows are "ok_*". The resolve-ok / classify-review row is missing_source.
    assert evidence["ok_empty_source"] == "1"
    assert evidence["missing_source"] == "1"
    assert evidence["unmatched"] == "1"
    assert evidence["conflict"] == "1"


def test_summarize_always_emits_source_evidence_core_keys():
    rows = [
        {
            "l1": "Human",
            "isolation_source": "blood",
            "classify_status": "ok",
            "classify_reason": "human-specimens",
            "paper_rollup": "Human",
            "healthcare_environment": "no",
        }
    ]
    summary = summarize(rows)
    evidence = {row["key"]: int(row["count"]) for row in summary if row["section"] == "source_evidence"}
    assert evidence["ok_with_source"] == 1
    assert evidence["ok_empty_source"] == 0
    assert evidence["unmatched"] == 0
    assert evidence["conflict"] == 0


def test_unclassified_does_not_make_a_cluster_mixed_and_count_is_isolates():
    rows = [
        {"l1": "Human", "cluster_id": "A"},
        {"l1": "Human", "cluster_id": "A"},
        {"l1": "Unclassified", "cluster_id": "A"},
    ]
    overlap = [row for row in summarize(rows) if row["section"] == "cluster_overlap"]
    assert overlap == [{"section": "cluster_overlap", "key": "A", "count": "3", "l1_set": "Human", "mixed": "no"}]


def test_dataset_section_reports_rows_that_share_a_biosample():
    rows = [
        {"query_id": "SAMN1", "biosample": "SAMN1"},
        {"query_id": "PDT1.1", "biosample": "SAMN1"},
        {"query_id": "GCA_1.1", "biosample": "SAMN2"},
        {"query_id": "not-an-id", "biosample": ""},
    ]
    dataset = {row["key"]: row["count"] for row in summarize(rows) if row["section"] == "dataset"}
    assert dataset == {"rows": "4", "unique_isolates": "3", "rows_sharing_a_biosample": "1"}
