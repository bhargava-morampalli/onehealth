from onehealth.summarize import summarize


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
