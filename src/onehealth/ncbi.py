"""NCBI Datasets and E-utilities lookups. No genome download."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Protocol

DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


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
        reports = payload.get("reports") or []
        if not reports:
            return None
        return reports[0]

    def genome_report(self, accession: str) -> dict[str, Any] | None:
        url = f"{DATASETS_BASE}/genome/accession/{urllib.parse.quote(accession)}/dataset_report"
        payload = self.fetcher.get_json(url)
        reports = payload.get("reports") or []
        if not reports:
            return None
        return reports[0]

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
        # Pathogen Detection isolate IDs are often stored as BioSample identifiers.
        uid = self._esearch_id("biosample", accession)
        if uid:
            name = self._esummary_accession("biosample", uid)
            return [name] if name else []
        uid = self._esearch_id("biosample", f"{accession}[All Fields]")
        if not uid:
            return []
        name = self._esummary_accession("biosample", uid)
        return [name] if name else []

    def _esearch_id(self, db: str, term: str) -> str:
        params = {
            "db": db,
            "term": term,
            "retmode": "json",
            "retmax": "1",
        }
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
        payload = self.fetcher.get_json(url)
        ids = (((payload.get("esearchresult") or {}).get("idlist")) or [])
        return str(ids[0]) if ids else ""

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

    def _esummary_accession(self, db: str, uid: str) -> str:
        params = {"db": db, "id": uid, "retmode": "json"}
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/esummary.fcgi?{urllib.parse.urlencode(params)}"
        payload = self.fetcher.get_json(url)
        result = payload.get("result") or {}
        record = result.get(str(uid)) or {}
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


def attributes_from_biosample(report: dict[str, Any]) -> dict[str, str]:
    attrs = {str(item.get("name", "")).lower(): str(item.get("value") or "") for item in report.get("attributes") or []}
    organism = ((report.get("description") or {}).get("organism") or {}).get("organism_name", "")
    fields = {
        "biosample": str(report.get("accession") or ""),
        "organism": str(organism or ""),
        "host": str(report.get("host") or attrs.get("host") or ""),
        "isolation_source": str(report.get("isolation_source") or attrs.get("isolation_source") or ""),
        "host_disease": str(report.get("host_disease") or attrs.get("host_disease") or ""),
        "sample_type": str(report.get("source_type") or report.get("tissue") or attrs.get("source_type") or attrs.get("tissue") or ""),
        "geo_loc_name": str(report.get("geo_loc_name") or attrs.get("geo_loc_name") or ""),
        "collection_date": str(report.get("collection_date") or attrs.get("collection_date") or ""),
        "strain": str(report.get("strain") or attrs.get("strain") or ""),
        "isolate": str(report.get("isolate") or attrs.get("isolate") or ""),
    }
    return fields


def biosample_from_genome(report: dict[str, Any]) -> str:
    info = report.get("assembly_info") or {}
    biosample = info.get("biosample") or {}
    if isinstance(biosample, dict):
        return str(biosample.get("accession") or "")
    return str(biosample or "")


def assembly_accession(report: dict[str, Any]) -> str:
    return str(report.get("accession") or (report.get("assembly_info") or {}).get("assembly_accession") or "")
