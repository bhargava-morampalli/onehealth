"""Parse mixed public NCBI identifiers."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable

ASSEMBLY_RE = re.compile(r"^(GC[AF])_(\d+)(?:\.(\d+))?$", re.IGNORECASE)
BIOSAMPLE_RE = re.compile(r"^(SAMN|SAMD|SAMEA|SAME)\d+$", re.IGNORECASE)
SRA_RE = re.compile(r"^([SED]R[RXPX])(\d+)$", re.IGNORECASE)
# PDT000000000 or PDT000000000.1
PDT_RE = re.compile(r"^(PDT\d+)(?:\.(\d+))?$", re.IGNORECASE)

ID_TYPES = ("assembly", "biosample", "sra", "pdt", "unknown")


@dataclass(frozen=True)
class ParsedId:
    raw: str
    normalized: str
    id_type: str


def normalize_id(value: str) -> str:
    return (value or "").strip()


def parse_id(value: str) -> ParsedId:
    raw = normalize_id(value)
    if not raw:
        return ParsedId(raw=value, normalized="", id_type="unknown")

    upper = raw.upper()
    if ASSEMBLY_RE.match(upper):
        return ParsedId(raw=raw, normalized=_canonical_assembly(upper), id_type="assembly")
    if BIOSAMPLE_RE.match(upper):
        return ParsedId(raw=raw, normalized=upper, id_type="biosample")
    if SRA_RE.match(upper):
        return ParsedId(raw=raw, normalized=upper, id_type="sra")
    if PDT_RE.match(upper):
        return ParsedId(raw=raw, normalized=_canonical_pdt(upper), id_type="pdt")
    return ParsedId(raw=raw, normalized=raw, id_type="unknown")


def parse_ids(values: Iterable[str]) -> list[ParsedId]:
    return [parse_id(value) for value in values]


def _canonical_assembly(value: str) -> str:
    match = ASSEMBLY_RE.match(value)
    assert match is not None
    prefix, digits, version = match.group(1).upper(), match.group(2), match.group(3)
    if version:
        return f"{prefix}_{digits}.{version}"
    return f"{prefix}_{digits}"


def _canonical_pdt(value: str) -> str:
    match = PDT_RE.match(value)
    assert match is not None
    stem, version = match.group(1).upper(), match.group(2)
    if version:
        return f"{stem}.{version}"
    return stem
