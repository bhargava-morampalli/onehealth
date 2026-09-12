# onehealth CLI test report

Tested `/home/bmorampa/onehealth-pipeline` as it stands (v0.1.0). No product-code changes. Existing pytest plus a throwaway harness under `/tmp/onehealth-extensive-test` (not added to the package).

`/cursor/stores/bc-aa037deb-470e-4509-8057-709134b82ef9/docs/` is not writable on this machine, so this file is the user-facing copy.

## Results and usefulness verdict

**Unit suite: 33 passed, 0 failed** (`uv run pytest -q`). Covers ID parse, overrides, offline resolve join, paper edge cases (wound drainage vs drain; clinical / clinical environment; unknown; Food/Laboratory/healthcare rollup), cluster/plasmid overlap only when columns exist, and CLI help + offline pipeline.

**Throwaway / live pass:** classify of published KlebNET + Pavia specimen strings and NCBI-messy phrasings mostly behaved as the current rules say. Live NCBI metadata calls **succeeded** (HTTP 200). The CLI still marked every live row `resolve_status=review`, and several public records came back with empty host/isolation_source. A BioSample that E-utilities lists as `host=Homo sapiens`, `isolation_source=urine` came back from NCBI Datasets as an empty report; the CLI wrote a blank row.

**Verdict:** `classify` on a **curated host/source TSV** is useful as a first pass and a review queue — it would cut hand-merging for clear clinical specimens, hospital drains, and unknown markers. It is **not** yet a replacement for that merge. Silent wrong or incomplete calls remain (see below). `resolve` is only useful if you read the metadata columns and ignore `resolve_status` on live runs, or if you supply `--metadata`. `summarize` / `plot` are fine count dumps once classify is trustworthy; on live NCBI output they currently show almost everything Unclassified, which looks like “no signal” rather than “metadata missing / flagged review.”

---

## What passed

| Area | Result |
| --- | --- |
| Existing pytest | 33/33 |
| Mixed IDs parsed | GCA/GCF (with/without version, mixed case), SAMN, SRR, PDT → correct `id_type`. `PRJNA225` and `not-an-id` → `unknown`. |
| Duplicate IDs | Two `SAMN02604091` lines both resolved from local metadata; not collapsed. |
| Empty / unknown / not collected / n/a / missing / unknown environment / unknown river | Unclassified + `review` (`unknown_source` or `missing_source`). Never Environmental. |
| Wound drainage vs hospital/sink/shower drain | Human vs Environmental + healthcare yes. Token `drain` ≠ `drainage`. Mixed case and `Wound-Drainage` OK. |
| `clinical`, `clinical environment`, `clinical/clinical environment` | Environmental + healthcare yes, not Human. |
| `clinical blood` + human host | Human (specimen outranks bare clinical). |
| Hospital sewage / hospital wastewater / hospital sink | Environmental + healthcare yes (when the source field is actually read). |
| Food (`retail meat`, `Food`, `carcass`, `milk`) | L1 Food, paper rollup Environmental. |
| Laboratory (`type strain`, `laboratory strain`, `reference strain`) | L1 Laboratory, paper rollup Environmental. |
| Override complete | Always wins (retail meat override stayed Food). |
| Override empty / conflict / invalid L1 (`Hospital`) | Review; L1 Unclassified; lexicon not applied as L1. |
| Host/source conflict (chicken+urine, chicken+blood, cow+wound, dog+stool, human+river) | Review `host_source_conflict` when both fields are seen. |
| Human host + hospital sink | Environmental + healthcare yes (allowed split). |
| Paper rollup matrix | Unclassified stays Unclassified even if healthcare=yes. Food, Laboratory, healthcare-env → Environmental. |
| Cluster/plasmid columns absent | No overlap sections. Present: ST15 marked mixed Human+Environmental. |
| Alternate column `source` (alone) | Treated as isolation_source. |
| Offline `--metadata` join | Works when the TSV actually contains matching id/assembly/biosample columns. |
| Plot | PNG written (offline ~30 KB, live ~18 KB). |
| Review sidecars | Written next to `-o` outputs. |

Pavia 2022 published diagnostic specimen words (Nature Microbiology 13, 1520–1531; “urine (n = 417) … blood (n = 74) … rectal swab (n = 70) … bronchial (n = 45) … sputum (n = 20)”) classify as Human when host is Homo sapiens. Those n= values are the paper’s, not counts from this CLI run.

KlebNET scheme strings that classified cleanly: `human blood`; `river water grab sample`; source type `Food` / `Environmental`; `Clinical environment: sinks and drains`; `Clinical environment: surfaces`; `Hospital wastewater`; `Wastewater (not hospital)`; `Seawater`; `Soil`; tissue Blood/Urine/Sputum/BAL/Rectal swab/Wound/Faeces.

---

## What failed or looked wrong (current behavior, not a patch list)

**Live resolve status is misleading.** NCBI Datasets returned assembly→BioSample, organism, and strain for `GCF_000005845.2` and `GCA_000016305.1`, and SRA→`SAMN00630374` for `SRR390728`. Every row stayed `resolve_status=review` (`ncbi_assembly` / `ncbi_sra` / empty reason). A scientist sorting on status would throw away successful lookups.

**Datasets missed attributes that E-utilities has.** `SAMN63084674` (public 2026-09-11): E-utilities summary has `host=Homo sapiens`, `isolation_source=urine`. Datasets `biosample_report` was `{"reports":[{}]}`. CLI output: empty biosample/host/source, `review`. Classify: `missing_source` / Unclassified.

**Reference / lab genomes have no isolation_source.** `SAMN02604091` (E. coli K-12 MG1655) and `SAMN02603941` (K. pneumoniae MGH 78578): strain only. Classify cannot see Laboratory from `strain=K-12`. Organism is not used.

**`organism=Homo sapiens` is not `host`.** `SRR390728` → BioSample `SAMN00630374` (cell line DB). CLI put `Homo sapiens` in `organism`. Classify saw empty host/source → `missing_source`. NCBI attributes were DSMZ / cell_line only.

**Empty synonym column wins.** A TSV with both `source` (filled) and `isolation source` (empty) classified as if source were missing: shower drain, hospital sink, type strain, clinical/clinical environment, Food, animal feed → `missing_source`. Same rows with only `source` classified correctly. First matching alias is used even when blank.

**Local GCA vs GCF pairing.** Offline metadata for `GCF_000005845.2` did not satisfy query `gca_000005845.2` (normalized `GCA_000005845.2`). Missing BioSample / unmatched assembly → `missing_metadata`. PDT `PDT000002407.1` → `pdt_unresolved` (E-utils count 0). `PRJNA225` is not an accepted ID type.

**Silent or shaky classes**

| Input | CLI | Why it matters |
| --- | --- | --- |
| KlebNET `Animal - food producing` | Food → paper Environmental | Phrase `food` inside an **animal** source-type. |
| `faucet` | Environmental, healthcare yes | Community tap looks like a hospital fixture. |
| Bare `sink` | unmatched / review | Inconsistent with `hospital sink` / `sink drain`. |
| `animal feed` | unmatched, or Human if host is human | KlebNET example isolate source; not in seed lexicon. |
| `Fresh water`, `Livestock`, `Companion animals`, source type `Human`, `Experimentally derived` | unmatched or review | Published vocab the lexicon does not cover. |
| `hospital carriage` / `hospital disease` + human host | Human via `human-host` only | Right for Pavia clinical isolates; the phrases themselves are not recognized. |
| `bronchial` + human host | Human via host only | Paper specimen; word `bronchial` is not in the specimen list. |
| Tissue `Unknown` + Homo sapiens | Human ok (`human-host`) | Unknown tissue does not force review if host is known. |
| Invalid override `l1=Hospital` on `Food` | Unclassified review | Correctly not silent, but the Food suggestion is easy to miss. |

---

## Live NCBI vs classify (exactly)

Network: allowed. No genomes/SRA downloads. Tool: `onehealth resolve` → Datasets v2 + E-utils, then `classify`.

| Query | NCBI / CLI metadata actually stored | classify |
| --- | --- | --- |
| `GCF_000005845.2` | BioSample `SAMN02604091`, organism E. coli K-12 MG1655, strain K-12. No host, no isolation_source. Status **review** / `ncbi_assembly`. | Unclassified, `missing_source` |
| `GCA_000016305.1` | BioSample `SAMN02603941`, organism K. pneumoniae MGH 78578, strain ATCC 700721; MGH 78578. No host/source. Status **review** / `ncbi_assembly`. | Unclassified, `missing_source` |
| `SAMN02604091` | Same K-12 record. Status **review**, empty reason. | Unclassified, `missing_source` |
| `SRR390728` | BioSample `SAMN00630374`. Datasets: cell_line DB, no isolation_source. CLI `organism=Homo sapiens`. Status **review** / `ncbi_sra`. | Unclassified, `missing_source` |
| `PDT000002407.1` | No BioSample. `pdt_unresolved`. | Unclassified, `missing_source` |
| `not-an-id` | `unknown_id`. | Unclassified, `missing_source` |
| `SAMN63084674` | E-utils: host Homo sapiens, isolation_source urine. Datasets: empty report. CLI: blank row, **review**. | Unclassified, `missing_source` |

---

## Usefulness by command

### resolve
**Clear:** ID type, linked BioSample/assembly/SRA when Datasets fills them, organism/strain on older complete genomes, `unknown_id` / `pdt_unresolved` / `missing_metadata` reasons, review sidecar.

**Still manual:** Almost all real NCBI isolation_source recovery (Datasets empty or fields absent). GCA↔GCF. PDT. Dedup. Anything in E-utils but not Datasets.

**Misleading:** `resolve_status=review` on rows that did fetch BioSample/organism. Blank BioSample row for a live accession that E-utils describes fully. `organism` looks like host and is not classified.

### classify
**Clear:** L1 + healthcare_environment + paper_rollup + reason + `lexicon_suggestion` + override flag. Review queue for unknown, empty override, conflict, unmatched. Wound vs drain and clinical vs clinical-environment work when the source column is the one being read.

**Still manual:** Every `unmatched` / `host_source_conflict` / unknown row; KlebNET gaps (animal feed, livestock, companion animals, fresh water, experimentally derived, Animal-food-producing); animal hosts with human specimen words; lab strains identified only by strain name; hospital carriage/disease/bronchial without relying on host.

**Misleading:** Empty alias column hides a real source. `Animal - food producing` → Food. `faucet` → healthcare-env. Live NCBI table → wall of Unclassified even after a “successful” resolve. Tissue Unknown + human host → confident Human.

### summarize
**Clear:** Counts by L1, paper_rollup, healthcare_environment, status, reason. Cluster/plasmid overlap and `mixed=yes` only if those columns exist (ST15 Human+Environmental).

**Still manual:** No per-isolate drill-down; mixed clusters still need the classified TSV.

**Misleading:** Live summary is “Unclassified / review” only — reads like classify failed globally. Does not say how many rows have usable NCBI text that was unused (organism, strain).

### plot
**Clear:** Two bars (L1 and paper rollup). File is a real PNG.

**Still manual:** No review vs ok, no healthcare split, no labels for why.

**Misleading:** Live PNG is almost all Unclassified. `--by` only retitles a panel; both series are always drawn.

---

## Where the lexicon still forces hand work

1. NCBI records with empty isolation_source (typical of reference genomes and many SRA-linked samples).
2. Published vocab not in the seed: animal feed, livestock, companion animals, fresh water, source type Human, experimentally derived, Animal - food producing (and the last one is worse because it auto-maps to Food).
3. Animal isolate + human-anatomy words (blood/urine/wound/stool) → review every time.
4. Clinical-environment vs Human when the only text is `clinical`, `clinical isolate`, or `hospital carriage`/`hospital disease` without a specimen.
5. Overrides that are empty or invalid: correct to queue, but then nothing is filled in.
6. Conflicting or unused columns (`source` vs empty `isolation source`; `organism` vs `host`).
7. PDT and unpaired GCA/GCF.

No feature work was done. Harness artifacts: `/tmp/onehealth-extensive-test/`.
