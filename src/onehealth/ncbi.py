"""NCBI Datasets and E-utilities lookups. No genome download."""

from __future__ import annotations

import difflib
import http.client
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Iterator, Protocol

from onehealth import __version__
from onehealth.constants import UNKNOWN_EXACT
from onehealth.errors import InputError
from onehealth.ids import parse_id, pdt_stem
from onehealth.text import normalize_text

DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# NCBI Pathogen Detection results, one folder per organism group (Klebsiella, Salmonella, ...).
PD_BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results"
# Accessions per batched Datasets request. 300 worked when tested; 100 keeps URLs short.
DATASETS_BATCH = 100
# The genome endpoint returns 20 reports per page unless asked for more.
_DATASETS_PAGE_SIZE = 1000
_DATASETS_MAX_PAGES = 20
# NCBI asks automated clients to say who they are; Python's default name is sometimes blocked.
USER_AGENT = f"onehealth/{__version__}"

_TEXT_FIELDS = (
    "organism",
    "host",
    "isolation_source",
    "host_disease",
    "sample_type",
    "source_type",
    "package",
    "geo_loc_name",
    "collection_date",
    "strain",
    "isolate",
)


class HttpError(RuntimeError):
    pass


class _Retryable(Exception):
    pass


class Fetcher(Protocol):
    def get_json(self, url: str, headers: dict[str, str] | None = None) -> dict[str, Any]: ...

    def get_text(self, url: str, headers: dict[str, str] | None = None) -> str: ...

    def iter_lines(self, url: str) -> Iterator[str]: ...


class UrllibFetcher:
    """Rate-limited GET. NCBI allows 3 requests/s without an API key and 10/s with one.

    HTTP 429 (rate limit), 5xx, network errors and garbled JSON are retried a few
    times, so a transient NCBI hiccup does not turn into a random review row.
    """

    DEFAULT_DELAY_S = 0.35
    API_KEY_DELAY_S = 0.11

    def __init__(self, delay_s: float = DEFAULT_DELAY_S, retries: int = 3, backoff_s: float = 2.0) -> None:
        self.delay_s = delay_s
        self.retries = retries
        self.backoff_s = backoff_s
        self._last = 0.0

    def get_json(self, url: str, headers: dict[str, str] | None = None) -> dict[str, Any]:
        """A JSON object. A body that is not one (an HTML maintenance page sent with
        HTTP 200, a cut-off transfer) is retried, then raised as HttpError."""
        for attempt in range(self.retries + 1):
            text = self.get_text(url, headers=headers)
            try:
                payload = json.loads(text)
            except ValueError as exc:
                problem = f"Invalid JSON from {_endpoint(url)}: {exc}"
            else:
                if isinstance(payload, dict):
                    return payload
                problem = f"Expected a JSON object from {_endpoint(url)}, got {type(payload).__name__}"
            if attempt == self.retries:
                raise HttpError(problem)
            time.sleep(self.backoff_s * (attempt + 1))
        raise AssertionError("unreachable")

    def get_text(self, url: str, headers: dict[str, str] | None = None) -> str:
        for attempt in range(self.retries + 1):
            try:
                return self._get_once(url, headers)
            except _Retryable as exc:
                if attempt == self.retries:
                    raise HttpError(str(exc)) from exc
                time.sleep(self.backoff_s * (attempt + 1))
        raise AssertionError("unreachable")

    def iter_lines(self, url: str) -> Iterator[str]:
        """A large text file line by line, never held in memory whole.

        Opening is retried like any request. A failure part-way through raises
        HttpError, because the read cannot resume where it stopped.
        """
        response = None
        for attempt in range(self.retries + 1):
            try:
                with _network_errors(url):
                    self._pause()
                    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                    response = urllib.request.urlopen(request, timeout=60)
                break
            except _Retryable as exc:
                if attempt == self.retries:
                    raise HttpError(str(exc)) from exc
                time.sleep(self.backoff_s * (attempt + 1))
        assert response is not None
        with response:
            try:
                with _network_errors(url):
                    for raw in response:
                        yield raw.decode("utf-8", errors="replace").rstrip("\r\n")
            except _Retryable as exc:
                raise HttpError(f"{exc} (part-way through the file)") from exc

    def _get_once(self, url: str, headers: dict[str, str] | None) -> str:
        with _network_errors(url):
            self._pause()
            request = urllib.request.Request(
                url, headers={"User-Agent": USER_AGENT, **(headers or {"Accept": "application/json"})}
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.read().decode("utf-8")

    def _pause(self) -> None:
        now = time.monotonic()
        wait = self.delay_s - (now - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()


def _endpoint(url: str) -> str:
    return url.split("?", 1)[0]


@contextmanager
def _network_errors(url: str) -> Iterator[None]:
    """HTTP 429, 5xx and network failures become _Retryable; other HTTP errors HttpError."""
    try:
        yield
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        # No query string: it carries api_key and email, and this text lands in the TSV.
        message = " ".join(f"HTTP {exc.code} for {_endpoint(url)}: {body[:200]}".split())
        if exc.code == 429 or exc.code >= 500:
            raise _Retryable(message) from exc
        raise HttpError(message) from exc
    except urllib.error.URLError as exc:
        raise _Retryable(f"Network error for {_endpoint(url)}: {exc.reason}") from exc
    except (OSError, http.client.HTTPException) as exc:
        # A read that times out, a connection dropped mid-reply or a TLS error escapes urlopen's URLError.
        raise _Retryable(f"Network error for {_endpoint(url)}: {type(exc).__name__} {exc}".strip()) from exc


@dataclass
class NCBIClient:
    email: str = ""
    api_key: str = ""
    fetcher: Fetcher = field(default_factory=UrllibFetcher)
    log: Callable[[str], None] | None = field(default=None, repr=False)  # progress notes for the user

    def __post_init__(self) -> None:
        # With a key NCBI allows 10 requests/s, not 3. Keep any pace the caller chose.
        if (
            self.api_key
            and isinstance(self.fetcher, UrllibFetcher)
            and self.fetcher.delay_s == UrllibFetcher.DEFAULT_DELAY_S
        ):
            self.fetcher.delay_s = UrllibFetcher.API_KEY_DELAY_S

    def _datasets_headers(self) -> dict[str, str]:
        """Datasets takes the key as a header (E-utilities as a parameter); it keeps the key out of URLs."""
        headers = {"Accept": "application/json"}
        if self.api_key:
            headers["api-key"] = self.api_key
        return headers

    def biosample_report(self, accession: str) -> dict[str, Any] | None:
        url = f"{DATASETS_BASE}/biosample/accession/{urllib.parse.quote(accession)}/biosample_report"
        payload = self.fetcher.get_json(url, headers=self._datasets_headers())
        return _usable_datasets_report(payload, require_accession=True)

    def biosample_reports(self, accessions: Iterable[str]) -> dict[str, dict[str, Any]]:
        """Datasets BioSample reports for many accessions, DATASETS_BATCH per request, keyed by upper-case accession.

        A fast path only. An accession NCBI leaves out (Datasets answers an unknown or
        empty record with an empty report), or one in a request that fails, is missing
        from the result; the caller then looks it up on its own, as for a single ID.
        """
        found: dict[str, dict[str, Any]] = {}
        for chunk in _chunks(accessions, "biosample"):
            for report in self._datasets_reports(f"{DATASETS_BASE}/biosample/accession/{_joined(chunk)}/biosample_report"):
                accession = str(report.get("accession") or "").strip().upper()
                if accession:
                    found.setdefault(accession, report)
        return found

    def genome_reports(self, accessions: Iterable[str]) -> dict[str, dict[str, Any]]:
        """Datasets assembly reports for many accessions, keyed by the (upper-case) accession asked for.

        Replies come back in no particular order. A versioned accession must come back as
        itself; an unversioned one takes the single report that shares its stem (Datasets
        answers with the current version). Anything else, such as an outdated version, is
        missing from the result, and the caller looks it up with genome_report, which also
        tries the unversioned stem.
        """
        found: dict[str, dict[str, Any]] = {}
        for chunk in _chunks(accessions, "assembly"):
            reports = self._datasets_reports(f"{DATASETS_BASE}/genome/accession/{_joined(chunk)}/dataset_report")
            by_accession = {str(r.get("accession") or "").strip().upper(): r for r in reports}
            by_stem: dict[str, list[dict[str, Any]]] = {}
            for accession, report in by_accession.items():
                if accession:
                    by_stem.setdefault(accession.split(".", 1)[0], []).append(report)
            for asked in chunk:
                if "." in asked:
                    report = by_accession.get(asked)
                else:
                    matches = by_stem.get(asked, [])
                    report = matches[0] if len(matches) == 1 else None
                if report is not None:
                    found[asked] = report
        return found

    def _datasets_reports(self, url: str) -> list[dict[str, Any]]:
        """Every non-empty report from a Datasets report URL, across pages. [] if any request fails."""
        reports: list[dict[str, Any]] = []
        token = ""
        try:
            for _ in range(_DATASETS_MAX_PAGES):
                query = {"page_size": str(_DATASETS_PAGE_SIZE), **({"page_token": token} if token else {})}
                payload = self.fetcher.get_json(f"{url}?{urllib.parse.urlencode(query)}", headers=self._datasets_headers())
                if not isinstance(payload, dict):
                    return []
                reports.extend(r for r in payload.get("reports") or [] if isinstance(r, dict) and r)
                token = str(payload.get("next_page_token") or "")
                if not token:
                    return reports
        except (HttpError, ValueError):
            pass  # the whole batch falls back to one-at-a-time lookups
        return []

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
        payload = self.fetcher.get_json(url, headers=self._datasets_headers())
        return _usable_datasets_report(payload, require_accession=False)

    def biosample_xml_report(self, accession: str) -> dict[str, Any] | None:
        uid = self._accession_uid("biosample", accession)
        if not uid:
            return None
        params = {"db": "biosample", "id": uid, "retmode": "xml"}
        self._add_eutils(params)
        url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        text = self.fetcher.get_text(url, headers={"Accept": "application/xml"})
        report = parse_biosample_xml(text)
        # A search of all fields can land on another sample whose text mentions this
        # accession. Its metadata must not be attached to this one.
        if report is not None and accession.strip().upper() not in report["ids"]:
            return None
        return report

    def _accession_uid(self, db: str, accession: str) -> str:
        """The record with this accession: the accession field first, then all fields (a
        replaced accession is found only there)."""
        return self._esearch_id(db, f"{accession}[accn]") or self._esearch_id(db, accession)

    def sra_to_biosample(self, accession: str) -> list[str]:
        uid = self._accession_uid("sra", accession)
        if not uid:
            return []
        links = self._elink(dbfrom="sra", db="biosample", uid=uid)
        biosamples: list[str] = []
        for link_id in links:
            accession_name = self._esummary_accession("biosample", link_id)
            if accession_name:
                biosamples.append(accession_name)
        return biosamples

    def note(self, message: str) -> None:
        """Pass a progress note to `log`, if there is one."""
        if self.log is not None:
            self.log(message)

    def pathogen_detection_groups(self, requested: Iterable[str]) -> list[str]:
        """Organism group names as NCBI spells them (case-insensitive). An unknown name is an InputError."""
        listing = self.fetcher.get_text(f"{PD_BASE}/", headers={"Accept": "text/html"})
        known = {name.lower(): name for name in re.findall(r'href="([A-Za-z0-9_]+)/"', listing)}
        groups = []
        for name in requested:
            match = known.get(name.strip().lower())
            if match is None:
                close = difflib.get_close_matches(name.strip(), list(known.values()), n=3)
                hint = f" Did you mean {', '.join(close)}?" if close else ""
                raise InputError(f"{name!r} is not a Pathogen Detection organism group.{hint} The list is at {PD_BASE}/")
            groups.append(match)
        return list(dict.fromkeys(groups))

    def pathogen_detection_targets(self, accessions: Iterable[str], groups: Iterable[str]) -> dict[str, dict[str, str]]:
        """Each PDT's current Pathogen Detection record (target, biosample, assembly), keyed by PDT stem.

        Each group's metadata from its newest release is streamed and not kept: about
        200 MB for Klebsiella, 900 MB for Salmonella. Reading stops as soon as every PDT
        is found. A release lists only the current version of each PDT, so an outdated or
        unversioned PDT gets the current one.
        """
        wanted = {pdt_stem(a) for a in accessions if parse_id(a).id_type == "pdt"}
        found: dict[str, dict[str, str]] = {}
        for group in groups:
            if not wanted - found.keys():
                break
            release, url = self._pd_metadata_url(group)
            before = len(found)
            lines = self.fetcher.iter_lines(url)
            try:
                header = next(lines, "").lstrip("#").split("\t")
                col = {name: i for i, name in enumerate(header)}
                if "target_acc" not in col or "biosample_acc" not in col:
                    raise HttpError(f"{url}: no target_acc and biosample_acc columns")
                for line in lines:
                    cells = line.split("\t")
                    stem = _pd_cell(cells, col, "target_acc").split(".", 1)[0].upper()
                    if stem in wanted and stem not in found:
                        found[stem] = {
                            "target": _pd_cell(cells, col, "target_acc"),
                            "biosample": _pd_cell(cells, col, "biosample_acc"),
                            "assembly": _pd_cell(cells, col, "asm_acc"),
                        }
                        if not wanted - found.keys():
                            break
            finally:
                close = getattr(lines, "close", None)
                if close is not None:
                    close()
            self.note(f"Pathogen Detection {group} {release}: {len(found) - before} PDTs found "
                       f"({len(found)} of {len(wanted)} in all)")
        return found

    def _pd_metadata_url(self, group: str) -> tuple[str, str]:
        """The newest release with a metadata table: latest_kmer, or else latest_snps (which lags behind)."""
        for link in ("latest_kmer", "latest_snps"):
            folder = f"{PD_BASE}/{group}/{link}/Metadata/"
            try:
                listing = self.fetcher.get_text(folder, headers={"Accept": "text/html"})
            except HttpError:
                continue
            releases = re.findall(r'href="(PDG\d+\.\d+)\.metadata\.tsv"', listing)
            if releases:
                return releases[0], f"{folder}{releases[0]}.metadata.tsv"
        raise HttpError(f"no metadata table for {group} under {PD_BASE}/{group}/")

    def pdt_to_biosample(self, accession: str) -> list[str]:
        """Resolve PDT as an assembly *name*, not a BioSample accession.

        Exact name first. Otherwise any assembly named PDTstem.N (a retired or
        unversioned PDT); several distinct BioSamples become multiple_biosamples.
        """
        found = self._assembly_name_biosamples(accession)
        if found:
            return found
        stem = pdt_stem(accession)
        return self._assembly_wildcard_biosamples(stem) if stem else []

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


def _pd_cell(cells: list[str], col: dict[str, int], name: str) -> str:
    """A Pathogen Detection metadata cell; NCBI writes NULL for empty."""
    index = col.get(name)
    value = cells[index].strip() if index is not None and index < len(cells) else ""
    return "" if value == "NULL" else value


def _chunks(accessions: Iterable[str], kind: str, size: int = DATASETS_BATCH) -> Iterator[list[str]]:
    """Distinct well-formed accessions of one kind, upper-case and in order, `size` at a time.

    Anything else is left for one-at-a-time lookup, so one odd value cannot make NCBI
    reject a whole batch.
    """
    parsed = (parse_id(item) for item in accessions if item and item.strip())
    unique = list(dict.fromkeys(p.normalized.upper() for p in parsed if p.id_type == kind))
    for start in range(0, len(unique), size):
        yield unique[start : start + size]


def _joined(accessions: list[str]) -> str:
    return urllib.parse.quote(",".join(accessions), safe=",")


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
        # source_type is NCBI's One Health category (human/animal/food/environmental);
        # sample_type is the material (whole organism, pure culture), so it comes last.
        "sample_type": str(
            report.get("source_type")
            or attrs.get("source_type")
            or report.get("tissue")
            or attrs.get("tissue")
            or report.get("sample_type")
            or attrs.get("sample_type")
            or ""
        ),
        # On its own too: where the submitter named a One Health class, classify uses it as is.
        "source_type": str(report.get("source_type") or attrs.get("source_type") or ""),
        # The BioSample package ("Pathogen.cl.1.0", "Pathogen.env.1.0"): what kind of
        # sample the submitter said this was.
        "package": str(report.get("package") or ""),
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
    ids = {accession.upper()} | {(el.text or "").strip().upper() for el in sample.findall("./Ids/Id")}
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
        "ids": sorted(ids - {""}),
        "package": (sample.findtext("./Package") or "").strip(),
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
    """Blank whole-value null tokens ("missing", "not collected"). Keep any other text.

    Partial values such as "unknown river" are kept verbatim for the curator;
    classify treats them as unknown.
    """
    cleaned = dict(fields)
    for name in _TEXT_FIELDS:
        value = str(cleaned.get(name) or "").strip()
        cleaned[name] = "" if normalize_text(value) in UNKNOWN_EXACT else value
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
