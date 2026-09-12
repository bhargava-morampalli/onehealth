"""NCBI Datasets and E-utilities lookups. No genome download."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any, Protocol

from onehealth.ids import parse_id, pdt_stem
from onehealth.text import is_unknown_value

DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_TEXT_FIELDS = (
    "organism",
    "host",
    "isolation_source",
    "host_disease",
    "sample_type",
    "geo_loc_name",
    "collection_date",
    "strain",
    "isolate",
)


class HttpError(RuntimeError):
    pass


class Fetcher(Protocol):
    def get_json(self, url: str, headers: dict[str, str] | None = None) -> dict[str, Any]: ...

    def get_text(self, url: str, headers: dict[str, str] | None = None) -> str: ...


class UrllibFetcher:
    def __init__(self, delay_s: float = 0.35) -> None:
        self.delay_s = delay_s
        self._last = 0.0

    def get_json(self, url: str, headers: dict[str, str] | None = None) -> dict[str, Any]:
        return json.loads(self.get_text(url, headers=headers))

    def get_text(self, url: str, headers: dict[str, str] | None = None) -> str:
        self._pause()
        request = urllib.request.Request(url, headers=headers or {"Accept": "application/json"})
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            raise HttpError(f"HTTP {exc.code} for {url}: {body[:200]}") from exc
        except urllib.error.URLError as exc:
            raise HttpError(f"Network error for {url}: {exc.reason}") from exc

    def _pause(self) -> None:
        now = time.monotonic()
        wait = self.delay_s - (now - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()


@dataclass
class NCBIClient:
    email: str = ""
    api_key: str = ""
    fetcher: Fetcher = field(default_factory=UrllibFetcher)

    def biosample_report(self, accession: str) -> dict[str, Any] | None:
        url = f"{DATASETS_BASE}/biosample/accession/{urllib.parse.quote(accession)}/biosample_report"
        payload = self.fetcher.get_json(url)
        return _usable_datasets_report(payload, require_accession=True)

    def genome_report(self, accession: str) -> dict[str, Any] | None:
        report = self._genome_report_once(accession)
        if report is not None:
            return report
        parsed = parse_id(accession)
        if parsed.id_type == "assembly" and "." in parsed.normalized:
            stem = parsed.normalized.rsplit(".", 1)[0]
            if stem != parsed.normalized:
                return self._genome_report_once(stem)
        return None

    def _genome_report_once(self, accession: str) -> dict[str, Any] | None:
        url = f"{DATASETS_BASE}/genome/accession/{urllib.parse.quote(accession)}/dataset_report"
        payload = self.fetcher.get_json(url)
        return _usable_datasets_report(payload, require_accession=False)

    def biosample_xml_report(self, accession: str) -> dict[str, Any] | None:
        uid = self._esearch_id("biosample", accession)
        if not uid:
            return None
        params = {"db": "biosample", "id": uid, "retmode": "xml"}
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        text = self.fetcher.get_text(url, headers={"Accept": "application/xml"})
        return parse_biosample_xml(text)

    def sra_to_biosample(self, accession: str) -> list[str]:
        uid = self._esearch_id("sra", accession)
        if not uid:
            return []
        links = self._elink(dbfrom="sra", db="biosample", uid=uid)
        biosamples: list[str] = []
        for link_id in links:
            accession_name = self._esummary_accession("biosample", link_id)
            if accession_name:
                biosamples.append(accession_name)
        return biosamples

    def pdt_to_biosample(self, accession: str) -> list[str]:
        """Resolve PDT as an assembly *name*, not a BioSample accession."""
        found = self._assembly_name_biosamples(accession)
        if found:
            return found
        stem = pdt_stem(accession)
        if stem and stem.upper() != accession.upper():
            found = self._assembly_name_biosamples(stem)
            if found:
                return found
            return self._assembly_wildcard_biosamples(stem)
        return []

    def _assembly_name_biosamples(self, name: str) -> list[str]:
        return self._biosamples_from_assembly_uids(
            self._esearch_ids("assembly", f"{name}[Assembly Name]", retmax=5)
        )

    def _assembly_wildcard_biosamples(self, stem: str) -> list[str]:
        uids = self._esearch_ids("assembly", f"{stem}*", retmax=5)
        biosamples: list[str] = []
        stem_upper = stem.upper()
        for uid in uids:
            record = self._esummary_record("assembly", uid)
            name = str(record.get("assemblyname") or "")
            if not name.upper().startswith(stem_upper):
                continue
            acc = str(record.get("biosampleaccn") or record.get("biosampleacc") or "")
            if acc:
                biosamples.append(acc)
        return list(dict.fromkeys(biosamples))

    def _biosamples_from_assembly_uids(self, uids: list[str]) -> list[str]:
        biosamples: list[str] = []
        for uid in uids:
            record = self._esummary_record("assembly", uid)
            acc = str(record.get("biosampleaccn") or record.get("biosampleacc") or "")
            if acc:
                biosamples.append(acc)
        return list(dict.fromkeys(biosamples))

    def _esearch_id(self, db: str, term: str) -> str:
        ids = self._esearch_ids(db, term, retmax=1)
        return ids[0] if ids else ""

    def _esearch_ids(self, db: str, term: str, retmax: int = 1) -> list[str]:
        params = {
            "db": db,
            "term": term,
            "retmode": "json",
            "retmax": str(retmax),
        }
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
        payload = self.fetcher.get_json(url)
        ids = (((payload.get("esearchresult") or {}).get("idlist")) or [])
        return [str(item) for item in ids if item]

    def _elink(self, dbfrom: str, db: str, uid: str) -> list[str]:
        params = {"dbfrom": dbfrom, "db": db, "id": uid, "retmode": "json"}
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/elink.fcgi?{urllib.parse.urlencode(params)}"
        payload = self.fetcher.get_json(url)
        ids: list[str] = []
        for linkset in payload.get("linksets") or []:
            for dbset in linkset.get("linksetdbs") or []:
                ids.extend(str(item) for item in dbset.get("links") or [])
        return ids

    def _esummary_record(self, db: str, uid: str) -> dict[str, Any]:
        params = {"db": db, "id": uid, "retmode": "json"}
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/esummary.fcgi?{urllib.parse.urlencode(params)}"
        payload = self.fetcher.get_json(url)
        result = payload.get("result") or {}
        record = result.get(str(uid)) or {}
        return record if isinstance(record, dict) else {}

    def _esummary_accession(self, db: str, uid: str) -> str:
        record = self._esummary_record(db, uid)
        for key in ("accession", "biosampleacc", "caption", "uid"):
            value = record.get(key)
            if value:
                return str(value)
        return ""

    def _add_eutils(self, params: dict[str, str]) -> None:
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        params["tool"] = "onehealth"


def _usable_datasets_report(payload: dict[str, Any], *, require_accession: bool) -> dict[str, Any] | None:
    if not isinstance(payload, dict) or payload == {}:
        return None
    reports = payload.get("reports") or []
    if not reports:
        return None
    report = reports[0]
    if not isinstance(report, dict) or report == {}:
        return None
    if require_accession and not str(report.get("accession") or "").strip():
        return None
    return report


def attributes_from_biosample(report: dict[str, Any]) -> dict[str, str]:
    attrs = {str(item.get("name", "")).lower(): str(item.get("value") or "") for item in report.get("attributes") or []}
    organism = ((report.get("description") or {}).get("organism") or {}).get("organism_name", "")
    fields = {
        "biosample": str(report.get("accession") or ""),
        "organism": str(organism or ""),
        "host": str(report.get("host") or attrs.get("host") or ""),
        "isolation_source": str(report.get("isolation_source") or attrs.get("isolation_source") or ""),
        "host_disease": str(report.get("host_disease") or attrs.get("host_disease") or ""),
        "sample_type": str(
            report.get("sample_type")
            or attrs.get("sample_type")
            or report.get("source_type")
            or attrs.get("source_type")
            or report.get("tissue")
            or attrs.get("tissue")
            or ""
        ),
        "geo_loc_name": str(report.get("geo_loc_name") or attrs.get("geo_loc_name") or ""),
        "collection_date": str(report.get("collection_date") or attrs.get("collection_date") or ""),
        "strain": str(report.get("strain") or attrs.get("strain") or ""),
        "isolate": str(report.get("isolate") or attrs.get("isolate") or ""),
    }
    return _blank_unknown_fields(fields)


def parse_biosample_xml(text: str) -> dict[str, Any] | None:
    if not (text or "").strip():
        return None
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return None
    sample = root if root.tag == "BioSample" else root.find("BioSample")
    if sample is None:
        return None
    accession = (sample.get("accession") or "").strip()
    organism_el = sample.find("./Description/Organism/OrganismName")
    if organism_el is None:
        organism_el = sample.find("./Description/Organism")
    organism = ""
    if organism_el is not None:
        organism = (organism_el.text or organism_el.get("taxonomy_name") or "").strip()
    attributes: list[dict[str, str]] = []
    for attr in sample.findall("./Attributes/Attribute"):
        name = (attr.get("harmonized_name") or attr.get("attribute_name") or "").strip()
        value = (attr.text or "").strip()
        if name:
            attributes.append({"name": name, "value": value})
    if not accession and not attributes:
        return None
    top = {item["name"].lower(): item["value"] for item in attributes}
    return {
        "accession": accession,
        "description": {"organism": {"organism_name": organism}},
        "attributes": attributes,
        "host": top.get("host", ""),
        "isolation_source": top.get("isolation_source", ""),
        "host_disease": top.get("host_disease", ""),
        "sample_type": top.get("sample_type", ""),
        "geo_loc_name": top.get("geo_loc_name", ""),
        "collection_date": top.get("collection_date", ""),
        "strain": top.get("strain", ""),
        "isolate": top.get("isolate", ""),
    }


def _blank_unknown_fields(fields: dict[str, str]) -> dict[str, str]:
    cleaned = dict(fields)
    for name in _TEXT_FIELDS:
        value = cleaned.get(name, "")
        cleaned[name] = "" if is_unknown_value(value) else str(value or "").strip()
    cleaned["biosample"] = str(cleaned.get("biosample") or "").strip()
    return cleaned


def biosample_from_genome(report: dict[str, Any]) -> str:
    info = report.get("assembly_info") or {}
    biosample = info.get("biosample") or {}
    if isinstance(biosample, dict):
        return str(biosample.get("accession") or "")
    return str(biosample or "")


def assembly_accession(report: dict[str, Any]) -> str:
    return str(report.get("accession") or (report.get("assembly_info") or {}).get("assembly_accession") or "")
