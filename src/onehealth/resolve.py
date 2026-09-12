"""Resolve mixed NCBI IDs to BioSample isolation-source metadata."""

from __future__ import annotations

from pathlib import Path

from onehealth.constants import RESOLVE_FIELDS
from onehealth.ids import ParsedId, assembly_join_keys, is_sra_study, parse_id
from onehealth.io_tsv import get_field, metadata_view, read_tsv, row_match_ids
from onehealth.ncbi import (
    HttpError,
    NCBIClient,
    attributes_from_biosample,
    assembly_accession,
    biosample_from_genome,
)


def load_local_metadata(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    return read_tsv(path)


def resolve_ids(
    ids: list[str],
    local_rows: list[dict[str, str]] | None = None,
    client: NCBIClient | None = None,
    allow_network: bool = False,
) -> list[dict[str, str]]:
    index = _index_local(local_rows or [])
    resolved: list[dict[str, str]] = []
    for raw in ids:
        parsed = parse_id(raw)
        row = _blank_row(parsed)
        local = _lookup_local(parsed, index)
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
        try:
            fetched = _fetch(parsed, client)
        except HttpError as exc:
            row["resolve_status"] = "review"
            row["resolve_reason"] = f"ncbi_error:{exc}"
            resolved.append(row)
            continue
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
        resolved.append(row)
    return resolved


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


def _index_local(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        view = {
            "biosample": metadata_view(row)["biosample"],
            "assembly": metadata_view(row)["assembly"],
            "sra": metadata_view(row)["sra"],
            "pdt": metadata_view(row)["pdt"],
            "organism": get_field(row, ("organism", "organism_name", "tax_name")),
            "host": metadata_view(row)["host"],
            "isolation_source": metadata_view(row)["isolation_source"],
            "host_disease": metadata_view(row)["host_disease"],
            "sample_type": metadata_view(row)["sample_type"],
            "geo_loc_name": get_field(row, ("geo_loc_name", "geo_loc", "country")),
            "collection_date": get_field(row, ("collection_date", "collection date")),
            "strain": get_field(row, ("strain",)),
            "isolate": get_field(row, ("isolate",)),
            "cluster_id": metadata_view(row)["cluster_id"],
            "plasmid_id": metadata_view(row)["plasmid_id"],
        }
        keys = row_match_ids(row)
        query = get_field(row, ("query_id", "id", "accession"))
        if query:
            keys.add(parse_id(query).normalized.upper() or query.upper())
        expanded = set(keys)
        for key in list(keys):
            expanded.update(assembly_join_keys(key))
        for key in expanded:
            index[key.upper()] = view
    return index


def _lookup_local(parsed: ParsedId, index: dict[str, dict[str, str]]) -> dict[str, str] | None:
    candidates = {parsed.normalized.upper(), parsed.raw.strip().upper()}
    candidates.update(assembly_join_keys(parsed.normalized))
    for key in filter(None, candidates):
        if key in index:
            hit = dict(index[key])
            if parsed.id_type == "assembly":
                hit["assembly"] = hit.get("assembly") or parsed.normalized
            elif parsed.id_type == "biosample":
                hit["biosample"] = hit.get("biosample") or parsed.normalized
            elif parsed.id_type == "sra":
                hit["sra"] = hit.get("sra") or parsed.normalized
            elif parsed.id_type == "pdt":
                hit["pdt"] = hit.get("pdt") or parsed.normalized
            return {k: v for k, v in hit.items() if v or k in RESOLVE_FIELDS}
    return None


def _fetch(parsed: ParsedId, client: NCBIClient) -> dict[str, str]:
    out: dict[str, str] = {}
    biosamples: list[str] = []
    if parsed.id_type == "biosample":
        biosamples = [parsed.normalized]
    elif parsed.id_type == "assembly":
        report = client.genome_report(parsed.normalized)
        if report is None:
            return {"resolve_status": "review", "resolve_reason": "assembly_not_found"}
        out["assembly"] = assembly_accession(report) or parsed.normalized
        biosample = biosample_from_genome(report)
        if not biosample:
            return {**out, "resolve_status": "review", "resolve_reason": "assembly_no_biosample"}
        biosamples = [biosample]
        out["resolve_reason"] = "ncbi_assembly"
    elif parsed.id_type == "sra":
        biosamples = client.sra_to_biosample(parsed.normalized)
        if not biosamples:
            return {"resolve_status": "review", "resolve_reason": "sra_no_biosample"}
        out["resolve_reason"] = "ncbi_sra"
    elif parsed.id_type == "pdt":
        biosamples = client.pdt_to_biosample(parsed.normalized)
        if not biosamples:
            return {"resolve_status": "review", "resolve_reason": "pdt_unresolved"}
        out["resolve_reason"] = "ncbi_pdt"
    else:
        return {"resolve_status": "review", "resolve_reason": "unknown_id"}

    unique = list(dict.fromkeys(item.strip() for item in biosamples if item and item.strip()))
    if len(unique) > 1:
        return {
            **out,
            "biosample": ",".join(unique),
            "resolve_status": "review",
            "resolve_reason": "multiple_biosamples",
        }
    if not unique:
        return {**out, "resolve_status": "review", "resolve_reason": "ncbi_incomplete"}

    known = unique[0]
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
