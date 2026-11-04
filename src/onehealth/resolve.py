"""Resolve mixed NCBI IDs to BioSample isolation-source metadata."""

from __future__ import annotations

from pathlib import Path

from onehealth.constants import (
    ASSEMBLY_ALIASES,
    BIOSAMPLE_ALIASES,
    ID_COLUMN_ALIASES,
    PDT_ALIASES,
    RESOLVE_FIELDS,
    SRA_ALIASES,
)
from onehealth.ids import ParsedId, assembly_join_keys, is_sra_study, parse_id, pdt_stem
from onehealth.io_tsv import ID_LIST_ALIASES, first_present, get_field, metadata_view, read_tsv, row_match_ids
from onehealth.errors import InputError
from onehealth.ncbi import (
    HttpError,
    NCBIClient,
    attributes_from_biosample,
    assembly_accession,
    biosample_from_genome,
)


def load_local_metadata(path: Path | None) -> list[dict[str, str]]:
    """The --metadata table. Without an ID column no row could match, and every ID would
    quietly be looked up at NCBI (or sent to review with --no-network), so that is an error."""
    if path is None:
        return []
    rows = read_tsv(path)
    if rows and first_present(rows[0], ID_LIST_ALIASES) is None:
        found = ", ".join(list(rows[0])[:8]) + (", ..." if len(rows[0]) > 8 else "")
        raise InputError(f"{path}: no ID column to match IDs on ({', '.join(ID_LIST_ALIASES[:6])}, ...; found: {found})")
    return rows


def resolve_ids(
    ids: list[str],
    local_rows: list[dict[str, str]] | None = None,
    client: NCBIClient | None = None,
    allow_network: bool = False,
    pd_groups: list[str] | None = None,
) -> list[dict[str, str]]:
    """One row per ID, in order. `pd_groups` (e.g. ["Klebsiella"]) looks PDTs up in those
    Pathogen Detection organism groups' metadata before trying anything else for them."""
    index = _index_local(local_rows or [])
    resolved: list[dict[str, str]] = []
    online: list[tuple[dict[str, str], ParsedId]] = []  # rows to look up at NCBI, filled in place below
    for raw in ids:
        parsed = parse_id(raw)
        row = _blank_row(parsed)
        local, conflict = _lookup_local(parsed, index)
        if conflict:
            row["resolve_status"] = "review"
            row["resolve_reason"] = "local_metadata_conflict"
            resolved.append(row)
            continue
        if local is not None:
            row.update(local)
            row["resolve_status"] = "ok"
            row["resolve_reason"] = "local_metadata"
            resolved.append(row)
            continue
        if parsed.id_type == "unknown" or not parsed.normalized:
            row["resolve_status"] = "review"
            row["resolve_reason"] = "unknown_id"
            resolved.append(row)
            continue
        if parsed.id_type == "sra" and is_sra_study(parsed.normalized):
            row["resolve_status"] = "review"
            row["resolve_reason"] = "sra_study_ambiguous"
            resolved.append(row)
            continue
        if not allow_network or client is None:
            row["resolve_status"] = "review"
            row["resolve_reason"] = "missing_metadata"
            resolved.append(row)
            continue
        resolved.append(row)
        online.append((row, parsed))
    if online and client is not None:
        _resolve_online(online, client, pd_groups or [])
    return resolved


def _pathogen_detection(online: list[tuple[dict[str, str], ParsedId]], client: NCBIClient, groups: list[str]) -> dict:
    """PDT stem -> current Pathogen Detection record, or {} (then PDTs take the assembly-name route)."""
    pdts = [parsed.normalized for _, parsed in online if parsed.id_type == "pdt"]
    if not pdts or not groups:
        return {}
    try:
        names = client.pathogen_detection_groups(groups)  # an unknown group name stops the run here
        return client.pathogen_detection_targets(pdts, names)
    except InputError:
        raise
    except (HttpError, ValueError) as exc:
        client.note(f"Pathogen Detection lookup failed ({exc}); PDTs are looked up by assembly name instead")
        return {}


def _resolve_online(online: list[tuple[dict[str, str], ParsedId]], client: NCBIClient, pd_groups: list[str]) -> None:
    """Look rows up at NCBI, filling each row in place.

    Assembly and BioSample reports are fetched many per request. Whatever a batch does
    not answer cleanly is looked up on its own, exactly as a single ID always was, so
    batching changes the number of requests, not the result.
    """
    targets = _pathogen_detection(online, client, pd_groups)
    genomes = client.genome_reports(parsed.normalized for _, parsed in online if parsed.id_type == "assembly")
    found: list[tuple[dict[str, str], dict[str, str], str]] = []
    for row, parsed in online:
        try:
            fetched, biosample = _find_biosample(parsed, client, genomes, targets)
        except (HttpError, ValueError) as exc:
            _ncbi_error(row, exc, client)
            continue
        if biosample:
            found.append((row, fetched, biosample))
        else:
            _finish(row, fetched)
    reports = client.biosample_reports(biosample for _, _, biosample in found)
    for row, fetched, biosample in found:
        try:
            _finish(row, _add_biosample(fetched, biosample, reports.get(biosample.upper()), client))
        except (HttpError, ValueError) as exc:
            _ncbi_error(row, exc, client)


def _ncbi_error(row: dict[str, str], exc: Exception, client: NCBIClient) -> None:
    # One bad answer must not end a run of thousands of IDs (ValueError: JSON from a
    # custom fetcher). NCBI echoes a rejected API key in its error body, and this text
    # lands in the table, so the key is masked.
    message = str(exc).replace(client.api_key, "***") if client.api_key else str(exc)
    row["resolve_status"] = "review"
    row["resolve_reason"] = f"ncbi_error:{message}"


def _finish(row: dict[str, str], fetched: dict[str, str]) -> None:
    row.update(fetched)
    if fetched.get("resolve_status") == "review":
        row["resolve_status"] = "review"
        row["resolve_reason"] = fetched.get("resolve_reason") or "ncbi_incomplete"
    elif _single_biosample(row.get("biosample")):
        row["resolve_status"] = "ok"
        row["resolve_reason"] = fetched.get("resolve_reason") or row.get("resolve_reason") or "ncbi"
    else:
        row["resolve_status"] = "review"
        row["resolve_reason"] = row.get("resolve_reason") or "ncbi_incomplete"


def _blank_row(parsed: ParsedId) -> dict[str, str]:
    row = {name: "" for name in RESOLVE_FIELDS}
    row["query_id"] = parsed.raw or parsed.normalized
    row["id_type"] = parsed.id_type
    if parsed.id_type == "assembly":
        row["assembly"] = parsed.normalized
    elif parsed.id_type == "biosample":
        row["biosample"] = parsed.normalized
    elif parsed.id_type == "sra":
        row["sra"] = parsed.normalized
    elif parsed.id_type == "pdt":
        row["pdt"] = parsed.normalized
    row["resolve_status"] = "review"
    row["resolve_reason"] = ""
    return row


def _single_biosample(value: str | None) -> bool:
    text = (value or "").strip()
    if not text or "," in text:
        return False
    return True


_ID_COLUMNS = {
    alias.lower()
    for alias in (*ID_COLUMN_ALIASES, *BIOSAMPLE_ALIASES, *ASSEMBLY_ALIASES, *SRA_ALIASES, *PDT_ALIASES)
}


def _local_view(row: dict[str, str]) -> dict[str, str]:
    """Standard fields plus any other local column (ST, patient ID, ...), minus identifiers."""
    view = metadata_view(row)
    standard = {
        "biosample": view["biosample"],
        "assembly": view["assembly"],
        "sra": view["sra"],
        "pdt": view["pdt"],
        "organism": get_field(row, ("organism", "organism_name", "tax_name")),
        "host": view["host"],
        "isolation_source": view["isolation_source"],
        "host_disease": view["host_disease"],
        "sample_type": view["sample_type"],
        "geo_loc_name": get_field(row, ("geo_loc_name", "geo_loc", "country")),
        "collection_date": get_field(row, ("collection_date", "collection date")),
        "strain": get_field(row, ("strain",)),
        "isolate": get_field(row, ("isolate",)),
        "cluster_id": view["cluster_id"],
        "plasmid_id": view["plasmid_id"],
    }
    extras = {
        key: value
        for key, value in row.items()
        if key and value and key not in standard and key not in RESOLVE_FIELDS and key.lower() not in _ID_COLUMNS
    }
    return {**standard, **extras}


def _index_local(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    """Map every identifier on a local row to the distinct metadata rows carrying it."""
    index: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        view = _local_view(row)
        for key in row_match_ids(row):
            bucket = index.setdefault(key, [])
            if view not in bucket:
                bucket.append(view)
    return index


def _lookup_local(
    parsed: ParsedId, index: dict[str, list[dict[str, str]]]
) -> tuple[dict[str, str] | None, bool]:
    """Return (metadata, conflict). Exact accession first, then versionless / GCA-GCF twin.

    Two different local rows for the same key are a conflict, not "last row wins".
    """
    exact = [parsed.normalized.upper(), parsed.raw.strip().upper()]
    related = sorted(assembly_join_keys(parsed.normalized))
    for key in [*exact, *related]:
        views = index.get(key) if key else None
        if not views:
            continue
        if len(views) > 1:
            return None, True
        hit = dict(views[0])
        column = parsed.id_type if parsed.id_type in {"assembly", "biosample", "sra", "pdt"} else ""
        if column:
            hit[column] = hit.get(column) or parsed.normalized
        return {k: v for k, v in hit.items() if v or k in RESOLVE_FIELDS}, False
    return None, False


def _find_biosample(
    parsed: ParsedId, client: NCBIClient, genomes: dict[str, dict], targets: dict[str, dict[str, str]]
) -> tuple[dict[str, str], str]:
    """The one BioSample an ID points to, with what was learnt on the way.

    Returns (fields, "") when the ID ends here, with a review reason in fields.
    """
    out: dict[str, str] = {}
    biosamples: list[str] = []
    if parsed.id_type == "biosample":
        biosamples = [parsed.normalized]
    elif parsed.id_type == "assembly":
        # From the batch when it answered; otherwise on its own (which also tries the unversioned stem).
        report = genomes.get(parsed.normalized.upper()) or client.genome_report(parsed.normalized)
        if report is None:
            return {"resolve_status": "review", "resolve_reason": "assembly_not_found"}, ""
        out["assembly"] = assembly_accession(report) or parsed.normalized
        biosample = biosample_from_genome(report)
        if not biosample:
            return {**out, "resolve_status": "review", "resolve_reason": "assembly_no_biosample"}, ""
        biosamples = [biosample]
        out["resolve_reason"] = "ncbi_assembly"
    elif parsed.id_type == "sra":
        biosamples = client.sra_to_biosample(parsed.normalized)
        if not biosamples:
            return {"resolve_status": "review", "resolve_reason": "sra_no_biosample"}, ""
        out["resolve_reason"] = "ncbi_sra"
    elif parsed.id_type == "pdt":
        record = targets.get(pdt_stem(parsed.normalized))
        if record and record["biosample"]:
            # Pathogen Detection's own table: the current PDT version, its BioSample and assembly.
            biosamples = [record["biosample"]]
            out["pdt"] = record["target"]
            if record["assembly"]:
                out["assembly"] = record["assembly"]
        else:
            # Only works when NCBI assembled the genome itself and named it after the PDT.
            biosamples = client.pdt_to_biosample(parsed.normalized)
        if not biosamples:
            return {"resolve_status": "review", "resolve_reason": "pdt_unresolved"}, ""
        out["resolve_reason"] = "ncbi_pdt"
    else:
        return {"resolve_status": "review", "resolve_reason": "unknown_id"}, ""

    unique = list(dict.fromkeys(item.strip() for item in biosamples if item and item.strip()))
    if len(unique) > 1:
        return {
            **out,
            "biosample": ",".join(unique),
            "resolve_status": "review",
            "resolve_reason": "multiple_biosamples",
        }, ""
    if not unique:
        return {**out, "resolve_status": "review", "resolve_reason": "ncbi_incomplete"}, ""
    return out, unique[0]


def _add_biosample(
    out: dict[str, str], known: str, report: dict | None, client: NCBIClient
) -> dict[str, str]:
    """BioSample attributes from the batch report, or looked up on its own (Datasets, then E-utilities XML)."""
    if report is None:
        report = client.biosample_report(known)
    if report is None:
        report = client.biosample_xml_report(known)
    if report is None:
        return {**out, "biosample": known, "resolve_status": "review", "resolve_reason": "biosample_not_found"}
    attrs = attributes_from_biosample(report)
    if not attrs.get("biosample"):
        attrs.pop("biosample", None)
    out.update(attrs)
    out["biosample"] = out.get("biosample") or known
    return out
