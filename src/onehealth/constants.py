"""Shared vocabulary for v1 classification."""

from __future__ import annotations

L1_VALUES = (
    "Human",
    "Animal",
    "Environmental",
    "Food",
    "Laboratory",
    "Unclassified",
)

PAPER_ROLLUP_VALUES = ("Human", "Environmental", "Animal", "Unclassified")

HEALTHCARE_YES = "yes"
HEALTHCARE_NO = "no"

SEED_LEXICON_VERSION = "3"

UNKNOWN_EXACT = frozenset(
    {
        "",
        "-",
        ".",
        "na",
        "n/a",
        "n.a",
        "n.a.",
        "none",
        "null",
        "nil",
        "missing",
        "unknown",
        "unk",
        "not applicable",
        "not collected",
        "not provided",
        "not available",
        "not recorded",
        "not known",
        "not specified",
        "restricted access",
        "unspecified",
        "unavailable",
        "missing data",
        "no data",
        "none provided",
        # Spreadsheet and pandas exports, and lab shorthand.
        "nan",
        "#n/a",
        "not determined",
        "undetermined",
    }
)

# If a value is "unknown" plus only these tails, it stays Unclassified.
UNKNOWN_GENERIC_TAILS = frozenset(
    {
        "environment",
        "environmental",
        "source",
        "origin",
        "host",
        "sample",
        "isolate",
        "type",
        "location",
        "site",
    }
)

ID_COLUMN_ALIASES = (
    "query_id",
    "id",
    "accession",
    "isolate_id",
    "sample_id",
)

# The spaced names are what NCBI Datasets (`dataformat tsv genome`) writes.
BIOSAMPLE_ALIASES = ("biosample", "biosample_acc", "biosample_accession", "biosample accession", "assembly biosample accession")
ASSEMBLY_ALIASES = ("assembly", "asm_acc", "assembly_accession", "genome", "assembly accession")
SRA_ALIASES = ("sra", "srr", "run", "run_accession")
PDT_ALIASES = ("pdt", "target_acc", "pdt_accession")

HOST_ALIASES = ("host",)
ISOLATION_SOURCE_ALIASES = (
    "isolation_source",
    "isolation source",
    "source",
    "sample_source",
)
HOST_DISEASE_ALIASES = ("host_disease", "host disease")
# NCBI `source_type` is a controlled One Health category (human, animal, food,
# environmental, other). NCBI `sample_type` is the material ("whole organism",
# "pure culture"), so it is read last.
SAMPLE_TYPE_ALIASES = ("source_type", "tissue", "sample_type")

CLUSTER_ALIASES = ("cluster_id", "cluster", "poppunk_cluster")
PLASMID_ALIASES = ("plasmid_id", "plasmid", "plasmid_replicon")

RESOLVE_FIELDS = (
    "query_id",
    "id_type",
    "biosample",
    "assembly",
    "sra",
    "pdt",
    "organism",
    "host",
    "isolation_source",
    "host_disease",
    "sample_type",
    "geo_loc_name",
    "collection_date",
    "strain",
    "isolate",
    "resolve_status",
    "resolve_reason",
)

CLASSIFY_EXTRA_FIELDS = (
    "l1",
    "healthcare_environment",
    "paper_rollup",
    "classify_status",
    "classify_reason",
    "lexicon_version",
    "override_applied",
    "lexicon_suggestion",
)


def derive_paper_rollup(l1: str, healthcare_environment: str) -> str:
    """Map L1 + healthcare-environment onto the paper's 3-class rollup.

    Food, Laboratory, and healthcare-environment roll to Environmental.
    Unknown / Unclassified stays Unclassified and never becomes Environmental.
    """
    if l1 == "Unclassified":
        return "Unclassified"
    if healthcare_environment == HEALTHCARE_YES:
        return "Environmental"
    if l1 == "Human":
        return "Human"
    if l1 == "Animal":
        return "Animal"
    if l1 in {"Environmental", "Food", "Laboratory"}:
        return "Environmental"
    return "Unclassified"
