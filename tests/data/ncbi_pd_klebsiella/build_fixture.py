"""Rebuild the NCBI Pathogen Detection Klebsiella test fixture (standard library only).

    python tests/data/ncbi_pd_klebsiella/build_fixture.py

Everything comes from one fixed NCBI Pathogen Detection release, so the output is
reproducible. Nothing is simulated:

* metadata.tsv       BioSample metadata as NCBI publishes it (host, isolation_source,
                     source_type, epi_type, collection date, country) plus the NCBI
                     SNP cluster (PDS_acc) of each isolate.
* features.tsv       NCBI's own AMRFinderPlus calls (AMR_genotypes, virulence_genotypes).
* snp_dists.tsv      Pairwise SNP differences within each selected cluster, counted from
                     NCBI's per-cluster SNP calls (variation.vcf.gz) the way snp-dists
                     counts an alignment: sites called in both isolates that differ.

Selection: SNP clusters of 4-30 isolates spanning at least two collection years.
The first MIXED clusters (by accession) that contain both NCBI epi_type "clinical"
and "environmental/other" isolates, plus the first CLINICAL_ONLY clusters for contrast.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import sys
import tarfile
import urllib.request
from pathlib import Path

RELEASE = "PDG000000012.2554"  # Klebsiella, released 2026-10-17
BASE = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/{RELEASE}"
OUT = Path(__file__).resolve().parent
MIN_SIZE, MAX_SIZE = 4, 30
MIXED, CLINICAL_ONLY = 16, 6
KEEP = (
    "target_acc",
    "biosample_acc",
    "asm_acc",
    "collection_date",
    "geo_loc_name",
    "host",
    "host_disease",
    "isolation_source",
    "source_type",
    "epi_type",
)


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=300) as response:
        return response.read()


def stream_lines(url: str):
    with urllib.request.urlopen(url, timeout=600) as response:
        for raw in io.TextIOWrapper(response, encoding="utf-8", newline=""):
            yield raw.rstrip("\r\n")


def clean(value: str) -> str:
    return "" if value in ("NULL", "") else value


def year(value: str) -> str:
    match = re.match(r"(\d{4})", value or "")
    return match.group(1) if match else ""


def main() -> None:
    print("cluster list ...", file=sys.stderr)
    cluster_of: dict[str, str] = {}
    for row in csv.DictReader(io.StringIO(fetch(f"{BASE}/Clusters/{RELEASE}.reference_target.cluster_list.tsv").decode()), delimiter="\t"):
        cluster_of[row["target_acc"]] = row["PDS_acc"]

    print("metadata (streamed) ...", file=sys.stderr)
    lines = stream_lines(f"{BASE}/Metadata/{RELEASE}.metadata.tsv")
    header = next(lines).lstrip("#").split("\t")
    col = {name: i for i, name in enumerate(header)}
    records: dict[str, dict[str, str]] = {}
    for line in lines:
        cells = line.split("\t")
        target = cells[col["target_acc"]]
        if target not in cluster_of:
            continue
        record = {name: clean(cells[col[name]]) for name in KEEP}
        record["AMR_genotypes"] = clean(cells[col["AMR_genotypes"]])
        record["virulence_genotypes"] = clean(cells[col["virulence_genotypes"]])
        record["amrfinder_version"] = clean(cells[col["amrfinder_version"]])
        records[target] = record

    members: dict[str, list[str]] = {}
    for target, pds in cluster_of.items():
        if target in records:
            members.setdefault(pds, []).append(target)

    def accession_number(pds: str) -> int:
        return int(re.match(r"PDS(\d+)", pds).group(1))

    mixed, clinical = [], []
    for pds in sorted(members, key=accession_number):
        targets = members[pds]
        if not MIN_SIZE <= len(targets) <= MAX_SIZE:
            continue
        if len({year(records[t]["collection_date"]) for t in targets} - {""}) < 2:
            continue
        epi = {records[t]["epi_type"] for t in targets}
        if {"clinical", "environmental/other"} <= epi and len(mixed) < MIXED:
            mixed.append(pds)
        elif epi == {"clinical"} and len(clinical) < CLINICAL_ONLY:
            clinical.append(pds)
        if len(mixed) == MIXED and len(clinical) == CLINICAL_ONLY:
            break
    selected = mixed + clinical
    print(f"selected {len(mixed)} mixed and {len(clinical)} clinical-only clusters", file=sys.stderr)

    pairs: list[tuple[str, str, int]] = []
    for pds in selected:
        archive = tarfile.open(fileobj=io.BytesIO(fetch(f"{BASE}/SNP_trees/{pds}.tar.gz")), mode="r:gz")
        vcf = next(m for m in archive.getmembers() if m.name.endswith(".variation.vcf.gz"))
        text = gzip.decompress(archive.extractfile(vcf).read()).decode()
        pairs.extend(vcf_distances(text, set(members[pds])))

    targets = sorted(t for pds in selected for t in members[pds])
    with (OUT / "metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([*KEEP, "PDS_acc"])
        for t in targets:
            writer.writerow([*(records[t][name] for name in KEEP), cluster_of[t]])
    with (OUT / "features.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["biosample", "feature", "class"])
        for t in targets:
            r = records[t]
            if not r["amrfinder_version"]:
                continue  # not genotyped by NCBI
            hits = [(g, "AMR") for g in genes(r["AMR_genotypes"])] + [(g, "VIRULENCE") for g in genes(r["virulence_genotypes"])]
            for gene, klass in hits or [("", "")]:
                writer.writerow([r["biosample_acc"], gene, klass])
    with (OUT / "snp_dists.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_a", "target_b", "snps"])
        writer.writerows(sorted(pairs))
    print(f"wrote {len(targets)} isolates, {len(pairs)} pairs to {OUT}", file=sys.stderr)


def genes(value: str) -> list[str]:
    """'blaKPC-3=COMPLETE,fosA=COMPLETE' -> ['blaKPC-3', 'fosA']."""
    return sorted({item.split("=", 1)[0] for item in value.split(",") if item})


def vcf_distances(text: str, wanted: set[str]) -> list[tuple[str, str, int]]:
    """Pairwise differences at sites called (GT present, FT PASS) in both isolates."""
    samples: list[str] = []
    calls: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith("##"):
            continue
        cells = line.split("\t")
        if line.startswith("#"):
            samples = cells[9:]
            continue
        keys = cells[8].split(":")
        gt, ft = keys.index("GT"), keys.index("FT") if "FT" in keys else None
        row = []
        for cell in cells[9:]:
            values = cell.split(":")
            ok = values[gt] not in (".", "") and (ft is None or values[ft] == "PASS")
            row.append(values[gt] if ok else "")
        calls.append(row)
    keep = [i for i, s in enumerate(samples) if s in wanted]
    out = []
    for x in range(len(keep)):
        for y in range(x + 1, len(keep)):
            i, j = keep[x], keep[y]
            d = sum(1 for row in calls if row[i] and row[j] and row[i] != row[j])
            a, b = sorted((samples[i], samples[j]))
            out.append((a, b, d))
    return out


if __name__ == "__main__":
    main()
