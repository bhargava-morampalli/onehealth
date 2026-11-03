"""Rebuild the IMP carbapenemase source-classification benchmark.

    python tests/data/imp_carbapenemase_sources/build_fixture.py

Truth labels are the curated source categories published with
Vezina, Morampalli et al. (2025) Nat Commun 17:183 (Supplementary Data, Table S1,
column isolation_type_nenad). The classifier input is what NCBI publishes for each
BioSample today, fetched from NCBI Datasets and mapped with onehealth's own
attributes_from_biosample, so the benchmark sees exactly what `onehealth resolve`
would produce. Nothing is simulated.

Output: samples.tsv, one row per public BioSample.
"""

from __future__ import annotations

import csv
import io
import json
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from onehealth.ncbi import attributes_from_biosample

SUPPLEMENT = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41467-025-66874-7/MediaObjects/41467_2025_66874_MOESM4_ESM.xlsx"
)
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/biosample/accession/{}/biosample_report"
BATCH = 200
OUT = Path(__file__).resolve().parent

# The paper's three source categories. "drain" ("drain drainage") was left unclassified
# in the paper (its nine unclassified samples), so the correct outcome is review.
TRUTH = {
    "human": "Human",
    "clinical": "Human",
    "animal": "Animal",
    "clinical environment": "Environmental",
    "clinical/clinical environment": "Environmental",
    "environmental/other": "Environmental",
    "drain": "Unclassified",
}
HEALTHCARE = {"clinical environment", "clinical/clinical environment"}


def main() -> None:
    print("supplementary data ...", file=sys.stderr)
    request = urllib.request.Request(SUPPLEMENT, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=300) as response:
        workbook = response.read()
    rows = list(read_sheet(workbook, "Table_S1"))
    header, body = rows[0], rows[1:]
    col = {name: i for i, name in enumerate(header)}

    curated: dict[str, set[str]] = {}
    paper_source: dict[str, str] = {}
    for row in body:
        biosample = row[col["BioSample"]]
        if not biosample.startswith("SAM"):
            continue  # no public BioSample: nothing for the classifier to read
        curated.setdefault(biosample, set()).add(row[col["isolation_type_nenad"]])
        paper_source.setdefault(biosample, row[col["isolation_source"]])

    accessions = sorted(curated)
    print(f"{len(accessions)} BioSamples; fetching NCBI metadata ...", file=sys.stderr)
    ncbi: dict[str, dict] = {}
    for start in range(0, len(accessions), BATCH):
        chunk = accessions[start : start + BATCH]
        url = DATASETS.format("%2C".join(chunk))
        with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=300) as r:
            for report in json.load(r).get("reports", []):
                if report.get("accession"):
                    ncbi[report["accession"]] = report

    columns = ["biosample", "host", "isolation_source", "host_disease", "sample_type", "package",
               "curated_label", "truth", "truth_healthcare", "paper_isolation_source"]
    written = 0
    with (OUT / "samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for accession in accessions:
            labels = curated[accession]
            truths = {TRUTH[label] for label in labels}
            if len(truths) != 1 or accession not in ncbi:
                continue  # conflicting curation, or NCBI no longer returns the record
            report = ncbi[accession]
            fields = attributes_from_biosample(report)
            writer.writerow([
                accession,
                fields["host"],
                fields["isolation_source"],
                fields["host_disease"],
                fields["sample_type"],
                str((report.get("package") or "")),
                ";".join(sorted(labels)),
                truths.pop(),
                "yes" if labels & HEALTHCARE else "no",
                paper_source[accession],
            ])
            written += 1
    print(f"wrote {written} BioSamples to {OUT / 'samples.tsv'} ({len(accessions) - written} skipped)", file=sys.stderr)


def read_sheet(workbook: bytes, name: str):
    """Yield the rows of one worksheet of an .xlsx file as lists of strings."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    archive = zipfile.ZipFile(io.BytesIO(workbook))
    sheet_rid = next(s.get(f"{{{ns['r']}}}id") for s in ET.fromstring(archive.read("xl/workbook.xml")).iter(f"{{{ns['m']}}}sheet")
                     if s.get("name") == name)
    rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    target = next(rel.get("Target") for rel in rels if rel.get("Id") == sheet_rid)
    shared = ["".join(t.text or "" for t in item.iter(f"{{{ns['m']}}}t"))
              for item in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("m:si", ns)]
    width = 0
    for _, element in ET.iterparse(archive.open("xl/" + target.lstrip("/").removeprefix("xl/"))):
        if element.tag != f"{{{ns['m']}}}row":
            continue
        cells: dict[int, str] = {}
        for cell in element.findall("m:c", ns):
            value = cell.find("m:v", ns)
            if value is None:
                continue
            text = shared[int(value.text)] if cell.get("t") == "s" else value.text or ""
            cells[column_number(cell.get("r"))] = text
        width = max(width, max(cells, default=-1) + 1)
        yield [cells.get(i, "") for i in range(width)]
        element.clear()


def column_number(reference: str) -> int:
    """'A1' -> 0, 'AB7' -> 27."""
    number = 0
    for letter in re.match(r"[A-Z]+", reference).group(0):
        number = number * 26 + ord(letter) - 64
    return number - 1


if __name__ == "__main__":
    main()
