# NCBI `resolve`: can it be made airtight?

**No product code was changed.** Live calls were metadata only (Datasets v2 reports, E-utilities search/summary/fetch, SRA runinfo). No FASTA, SRA reads, or genome packages were downloaded.

Assigned store path `/cursor/stores/bc-aa037deb-470e-4509-8057-709134b82ef9/docs/` is not writable on this machine. This file is the user-facing copy.

Related test log: `docs/cli-test-report.md`. Fresh probes: `/tmp/onehealth-resolve-analysis-probe/` (not in the package).

---

## Lead answer

**Yes, for identity and honesty. No, for filling every isolation_source.**

`resolve` can be made airtight as an **ID → BioSample identity** step: one output row per input ID, the linked accessions are the ones NCBI actually has, status is `ok` only when that identity fetch succeeded, and host/source cells are either copied from a named NCBI attribute or left empty. It cannot be made airtight as “every public ID comes back with a usable host and isolation_source.” Many real records have no those fields. That gap belongs to `classify` (Unclassified + review) and to a curator override TSV, not to silent invention in `resolve`.

**Airtight here means:**

1. Every input line becomes exactly one output row (duplicates are not collapsed).
2. The parser’s accepted ID shapes are documented and deterministic (case-folded; version kept if the user typed one).
3. If NCBI has a **single** BioSample for that ID, the row has that accession. If it has none, or more than one, or the fetch failed, the row says so and is `review`.
4. `resolve_status=ok` means **identity resolved** (BioSample known, or a local `--metadata` hit), not “host and isolation_source are filled.”
5. Empty isolation_source on a real BioSample is still `ok`. `classify` then emits `missing_source`.
6. Fields are copied from named NCBI attributes only. `organism` never becomes `host`. Strain, package, title, and Pathogen Detection labels are never turned into host or isolation_source.
7. NCBI tokens that mean “no value” (`missing`, `not collected`, `n/a`, …) are stored as empty, not as a host named “missing.”
8. A Datasets `{"reports":[{}]}` is **not** a record. It must not wipe a BioSample the caller already knew.

**Recommended cascade (smallest path that meets that contract):**

1. **Parse** the line (mixed-case → canonical).
2. **Local `--metadata` first.** Join on any ID column. Skip blank alias columns. Also index the versionless assembly stem and the GCA↔GCF twin.
3. **Network, by type:**
   - **Assembly (GCA/GCF):** Datasets `genome/.../dataset_report`. Read `assembly_info.biosample.accession`, `current_accession`, `paired_accession`. If the body is `{}` and the query was versioned, retry the unversioned stem (Datasets already maps `GCA_000005845` → `.2`). If there is still no BioSample, stop: `assembly_no_biosample`.
   - **SRA run/experiment (SRR/SRX and E/D equivalents):** E-utils `esearch db=sra` → `elink db=biosample`. If that is empty, SRA runinfo `BioSample` column. If several distinct BioSamples: `multiple_biosamples`. Study-level SRP/ERP/DRP: do not take `retmax=1` and pretend it is one sample.
   - **BioSample (SAMN/SAMD/SAMEA/SAME):** Datasets `biosample/.../biosample_report`. If the report is missing, empty `{}`, or has no `accession`, E-utils `esearch` + `efetch retmode=xml` and read `<Attribute harmonized_name=...>`.
   - **PDT:** E-utils `esearch db=assembly` on `[Assembly Name]`, not `db=biosample`. Datasets `genome/accession/PDT…` returns `{}`. If the exact version misses, optionally try the current `PDTstem.N` assembly name and record the version move. If still none: `pdt_unresolved`. Do **not** use nucleotide `All Fields` (it can hit the right WGS comment, and it can also look like a hit when you have not checked the assembly name).
4. **Map:** `host` ← host only; `isolation_source` ← isolation_source only; `organism` stays organism; `sample_type` ← `sample_type` / `source_type` / `tissue`. Blank NCBI unknown tokens.
5. **Status:** `ok` if BioSample is non-empty after a successful fetch, or if local metadata hit. Otherwise `review` with a specific reason.

Pathogen Detection BigQuery / FTP isolate tables can resolve every live `target_acc`, but that is a new service, not the smallest change.

---

## What the code does today

### Parse (`src/onehealth/ids.py`)

| Input shape | Regex | Normalized | `id_type` |
| --- | --- | --- | --- |
| `GCA_` / `GCF_` + digits, optional `.version` | `^(GC[AF])_(\d+)(?:\.(\d+))?$` | prefix upper, version kept if present | `assembly` |
| `SAMN` / `SAMD` / `SAMEA` / `SAME` + digits | `SAMEA` before `SAME` | upper | `biosample` |
| `S/E/D` + `R` + `R\|X\|P` + digits | `^([SED]R[RXPX])(\d+)$` | upper | `sra` |
| `PDT` + digits, optional `.version` | `^(PDT\d+)(?:\.(\d+))?$` | stem upper, version kept if present | `pdt` |
| anything else, including blank | — | raw or `""` | `unknown` |

Mixed case is folded (`gca_…` → `GCA_…`). Leading/trailing space is stripped.

**Parse gaps (live-checked):**

| Typed | Today | Why it matters |
| --- | --- | --- |
| `SRS212581`, `ERS1932678`, `DRS000001` | `unknown` | Comment says sample (`S`) is allowed; the class is `[RXPX]` (R, X, P only). SRA *sample* accessions never enter the SRA cascade. |
| `SRP020237` | `sra` | Study ID. E-utils count **230**. CLI `retmax=1` then links **one** BioSample (`SAMN05182466`). Silent wrong isolate. |
| `SAMN02604091.1` | `unknown` | Version suffix on BioSample is not a real accession form. |
| `GCA000005845.2` (no underscore) | `unknown` | |
| `PRJNA225`, `PRJNA288601` | `unknown` | BioProject is not an isolate ID. |
| `GCA_000005845` vs `GCA_000005845.2` | different normalized strings | Local join is exact. Network Datasets already resolves the unversioned stem to `.2`. |
| `GCA_000005845.2` vs `GCF_000005845.2` | different strings | Network finds the same BioSample. Local `--metadata` does not cross the twin. |

### Network calls (`src/onehealth/ncbi.py`)

Default client: `UrllibFetcher` at 0.35 s between requests, 30 s timeout. E-utils params add `tool=onehealth` and optional `--email` / `--api-key`. **No genome or SRA download.**

| `id_type` | Calls | What is kept from the hit |
| --- | --- | --- |
| `assembly` | `GET https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{id}/dataset_report` | `accession` / `assembly_info.assembly_accession`; `assembly_info.biosample.accession` |
| `sra` | E-utils `esearch db=sra` (`retmax=1`) → `elink db=biosample` → `esummary` keys `accession`, `biosampleacc`, `caption`, `uid` | list of BioSample accessions |
| `pdt` | E-utils `esearch db=biosample` on the PDT string, then the same with `{id}[All Fields]` | **does not work** for live PDTs (0 hits) |
| `biosample` (always, after a BioSample is known) | `GET …/datasets/v2/biosample/accession/{id}/biosample_report` | first `reports[]` element, even if `{}` |

`reports` empty or missing → `None`. `{"reports":[{}]}` is treated as a real report.

Unused by the CLI today, but present on Datasets genome reports: `current_accession`, `paired_accession`, `assembly_info.assembly_name` (often the PDT string), nested BioSample attributes (host, isolation_source).

### Field mapping (`attributes_from_biosample`)

Attribute names are lowercased. Values come from top-level report keys first, then `attributes[]`.

| Output column | Source | Not read |
| --- | --- | --- |
| `biosample` | `report.accession` | |
| `organism` | `description.organism.organism_name` | |
| `host` | `report.host` or attr `host` | `organism` (correctly not copied) |
| `isolation_source` | `report.isolation_source` or attr `isolation_source` | `source`, `env_medium`, `ifsac_category`, title |
| `host_disease` | `host_disease` | |
| `sample_type` | `source_type` / `tissue` only | attr **`sample_type`** (present on SAMN63084674 and SAMN13648599) |
| `geo_loc_name` | `geo_loc_name` | |
| `collection_date` | `collection_date` | |
| `strain` | `strain` | `substrain` |
| `isolate` | `isolate` | |

Host is not filled from organism. That part is already correct.

### Status (`src/onehealth/resolve.py`)

`_blank_row` always sets `resolve_status=review`, `resolve_reason=""`, and pre-fills the typed ID column (`assembly` / `biosample` / `sra` / `pdt`).

Successful `_fetch` usually **omits** `resolve_status`. After `row.update(fetched)` the blank-row `review` is still there. The post-check is:

```text
if status != "review" and (host or isolation_source or biosample):
    status = ok
elif status != "review":
    status = review / ncbi_incomplete
```

So a live success never becomes `ok`. That is why every live row in `docs/cli-test-report.md` was `review` even when BioSample, organism, and strain were filled.

`_fetch` *does* set `review` itself for: `assembly_not_found`, `assembly_no_biosample`, `sra_no_biosample`, `pdt_unresolved`, `unknown_id`, `multiple_biosamples`, `biosample_not_found`, plus CLI-level `unknown_id`, `missing_metadata`, `ncbi_error:…`.

Local `--metadata` sets `ok` / `local_metadata` and skips the network.

Empty Datasets report wipe: `attributes_from_biosample({})` sets `biosample=""`. `out.setdefault("biosample", known)` does not restore it. For a BioSample query, `_blank_row` had the accession; `update` then erases it. That is `SAMN63084674` → blank row.

### Local join and alias shadowing (`io_tsv.get_field` / `first_present`)

`first_present` returns the first alias **whose column exists**, even if every cell is blank.

`ISOLATION_SOURCE_ALIASES = isolation_source, "isolation source", source, sample_source`.

A TSV with empty `isolation_source` and filled `source` joins as isolation_source empty. Offline resolve still marks `ok` (local hit). Classify then sees no source.

Local index keys are exact `parse_id(...).normalized.upper()` plus raw upper. No versionless stem, no GCA↔GCF twin.

---

## Live NCBI (metadata only, 2026-09-12)

HTTP 200 unless noted. Datasets TSV `Accept: text/tab-separated-values` on BioSample reports returned **500** (not a usable fallback).

### Assemblies

| Query | Datasets genome | BioSample | Notes |
| --- | --- | --- | --- |
| `GCF_000005845.2` | `GCF_000005845.2`, paired `GCA_000005845.2` | `SAMN02604091` | K-12. Nested BioSample attrs: strain, substrain only. |
| `GCA_000005845.2` | twin of the above | same | |
| `GCA_000005845` (no version) | **current** `GCA_000005845.2` | same | Datasets version-joins. CLI local join does not. |
| `GCA_000005845.1` | `{}` (2 bytes) | — | Superseded version. CLI: `assembly_not_found`. |
| `GCA_000016305.1` / `GCF_000016305.1` | paired to each other | `SAMN02603941` | Kp MGH 78578, strain only. |
| `GCA_014136285.1` | current, name **`PDT000806148.1`** | `SAMN09405588` | Nested attrs include host + isolation_source. |
| `GCF_000001405.40` (GRCh38.p14) | current, paired `GCA_000001405.29` | **none** | Reference human. CLI: `assembly_no_biosample`. |
| `GCA_000001405.29` | same pair | **none** | |
| `GCA_000002135.3` | current | **none** | Another assembly without BioSample. |
| `GCA_000001215.4` | Drosophila | `SAMN02803731` | Has BioSample; tissue/strain, not isolation_source. |
| `GCF_000009045.1` | *B. subtilis* 168 | `SAMEA3138188` | ENA-style BioSample, package-ish attrs. |

### BioSamples

| Query | Datasets `biosample_report` | E-utils `efetch` XML | CLI today |
| --- | --- | --- | --- |
| `SAMN63084674` | `{"reports":[{}]}` (16 bytes) | host=`Homo sapiens`, isolation_source=`urine`, sample_type=`whole organism`, isolate, geo=USA, date=2026. Package Microbe.1.0. Submitted **2026-09-11** | Blank row, `review`, biosample wiped |
| `SAMN02604091` | Generic.1.0, strain K-12, substrain MG1655, organism E. coli K-12 MG1655 | same; no host/source | Fields kept; status stuck `review` |
| `SAMN02603941` | Generic.1.0, strain only | same | same |
| `SAMN00630374` | Generic.1.0, organism **Homo sapiens**, attrs DSMZ / cell_line | same; no host/isolation_source | `organism=Homo sapiens`, empty host |
| `SAMN10527295` | Pathogen.cl.1.0, host=`missing`, isolation_source=`urine` | (not re-fetched) | Would store host=`missing` as a real host string |
| `SAMEA104307700` | Generic.1.0, host=`missing`, isolation_source=`urine` | — | SAMEA parse works |
| `SAMN09405588` | Pathogen.cl.1.0, host=`Homo sapiens`, isolation_source=`blood` | same | Full clinical record; Datasets works |
| `SAMD00060925` | Generic.1.0, strain ORS 571, type-material; no host/source | same | SAMD parse works |
| `SAMN00001583` | HapMap / Coriell human; no host/isolation_source | same | SRA-linked 1000G sample |
| `SAMN13648599` | Microbe.1.0, host=`chicken`, isolation_source=`anal swab`, **sample_type=`pure culture`** | isolation in sampledata | CLI drops `sample_type` |
| `SAMN12706423` | Pathogen.env.1.0, isolation_source=`wastewater`, no host | yes | Env package, empty host is real |
| `SAMN13151938` | Pathogen.env.1.0, isolation_source=`Wastewater Effluent` | yes | |
| `SAMN02344983` | Pathogen.env.1.0, isolation_source=`papaya`, `source_type` present, **no host** | — | Current BioSample for PDT000002407.**4** |

### SRA

| Query | esearch count | elink BioSample | runinfo `BioSample` |
| --- | --- | --- | --- |
| `SRR390728` | 1 | `SAMN00630374` | `SAMN00630374` (runinfo also lists sister run `SRR292241`, same BioSample) |
| `SRX079566` | 1 (same uid as the run) | same | same |
| `SRR000001` | 1 | `SAMN00001583` | same (file also lists other runs of the experiment) |
| `ERR036001` | 1 | uid → `SAMEA958765` | `SAMEA958765` |
| `SRP020237` | **230** | CLI would link **one** uid → `SAMN05182466` | one experiment’s runs, not the study |

elink and runinfo agreed on the run/experiment IDs above. Study-level IDs are the hazard.

### PDT / Pathogen Detection

PDT is **`target_acc` / assembly name**, not an Entrez BioSample accession. Official bulk metadata is BigQuery `ncbi-pathogen-detect.pdbrowser.isolates` (`target_acc`, `biosample_acc`, `asm_acc`). Isolates Browser HTML is not an API. Datasets `genome/accession/PDT…` is `{}`. FTP `ftp.ncbi.nlm.nih.gov/pathogen/Results/` did not resolve from this host.

| Query | `db=biosample` | `db=assembly` | Datasets genome | What actually exists |
| --- | --- | --- | --- | --- |
| `PDT000002407.1` | 0 | 0 | `{}` | Retired version. CLI: `pdt_unresolved` (correct given today’s search). |
| `PDT000002407` (no version) | 0 | 0 | — | Nucleotide All Fields hits WGS `AAHGAG000000000` because the **comment** says assembly name `PDT000002407.4`. |
| `PDT000002407.4` | 0 | **1** → `GCA_007182035.1`, name `PDT000002407.4`, BioSample `SAMN02344983` | `{}` | Salmonella Mbandaka, isolation_source **papaya**. Not Klebsiella. |
| `PDT000806148.1` | 0 | **1** → `GCA_014136285.1`, BioSample `SAMN09405588` | `{}` | *E. coli*, host Homo sapiens, blood. |
| `PDT000000001.1` | 0 | 0 | — | No Entrez hit. |

Current `pdt_to_biosample` only searches `db=biosample`. Live PDT000806148.1 therefore stays unresolved even though assembly name search is one E-utils call.

Nucleotide `All Fields` on an unversioned PDT stem is **not** a safe fallback: it worked here only because the WGS annotation contains `PDT000002407.4`. It must not be used without checking that the assembly name equals the query (or the current version of that stem).

---

## Failure modes a user will hit

| # | What the user has | What happens now | What airtight does |
| --- | --- | --- | --- |
| 1 | Live NCBI success, empty host/source (K-12, MGH 78578) | BioSample filled, **status review** (`ncbi_assembly` or empty reason) | `ok` + empty host/source. Classify: `missing_source`. Do not invent Laboratory from `strain=K-12`. |
| 2 | Empty Datasets report, full E-utils XML (`SAMN63084674`) | `{"reports":[{}]}` treated as a record; biosample **wiped**; review | Treat as missing; efetch XML; keep `SAMN63084674`; map host/urine; `ok`. |
| 3 | `organism=Homo sapiens`, no host (`SRR390728` / `SAMN00630374`) | Organism filled; classify `missing_source` | Same fields. Never copy organism → host. Cell-line / HapMap records stay empty host. |
| 4 | NCBI host=`missing` (`SAMN10527295`, `SAMEA104307700`) | Literal `missing` stored in `host` | Blank it (same unknown token list classify already uses). Keep `isolation_source=urine`. |
| 5 | Package-only / Generic.1.0 (K-12, SAMD type strain) | Strain + organism only | `ok` if BioSample known. Package name is not isolation_source. |
| 6 | Pathogen.cl vs Pathogen.env vs Microbe.1.0 | Package ignored (good) | Keep ignoring package as a class. Env package + `wastewater` is still just the isolation_source text. |
| 7 | Pathogen Detection PDT vs BioSample | `db=biosample` 0 hits; `pdt_unresolved` | Assembly-name lookup. Exact version first. Optional current `PDTstem.N`. Else `pdt_unresolved` + `--metadata`. No BigQuery in the small change. |
| 8 | PDT version suffix (`.1` vs `.4`) | `.1` unresolved; user never sees papaya / `SAMN02344983` | Either `pdt_unresolved` (exact) or documented version-follow with reason `pdt_version_updated` and `pdt` column = current name. Do not silently treat nucleotide All Fields as identity. |
| 9 | Unversioned vs versioned GCA | Network: unversioned works. Local: `GCA_000005845` ≠ `.2`. Superseded `.1` → `{}` | Network retry stem on `{}`. Local index versionless + current. Write `current_accession` if you add one column; not required for airtight status. |
| 10 | GCA ↔ GCF twin | Network same BioSample. Local `GCA_…` vs metadata `GCF_…` → `missing_metadata` | Index `paired_accession` / swap GC[AF]. |
| 11 | Mixed case | Parse OK. Local join OK if indexed upper. | Keep uppercasing. |
| 12 | Duplicate IDs on two lines | Two rows (already) | Keep two rows. |
| 13 | SRA-only run (`SRR390728`) | BioSample via elink; status review | Same link; `ok` if one BioSample. |
| 14 | SRA sample `SRS` / `ERS` | `unknown_id` | Parse as `sra` (fix `[RXPX]` → include `S`) then same cascade. |
| 15 | SRA study `SRP` | Silent one-of-many BioSample | `review` / `sra_study_ambiguous` (smallest). Expanding every run is a feature, not the small fix. |
| 16 | Assembly without BioSample (GRCh38, `GCA_000002135.3`) | `assembly_no_biosample` if `_fetch` returns that; then still stuck review | Keep the reason; status stays `review`. Do not invent a sample. |
| 17 | Several BioSamples for one ID | `multiple_biosamples`, comma-joined accessions, review | Keep. Do not pick the first. |
| 18 | Column alias shadowing | Empty `isolation_source` hides filled `source` | `get_field`: first **non-blank** alias. |
| 19 | `sample_type` attribute | Dropped (`source_type`/`tissue` only) | Read `sample_type` too. Still not a substitute for isolation_source. |
| 20 | No `--metadata`, `--no-network` | `missing_metadata` | Keep. |
| 21 | HTTP / DNS failure | `ncbi_error:…` | Keep. |
| 22 | BioProject / junk / blank | `unknown_id` | Keep. Do not search those strings as All Fields. |
| 23 | New BioSample, Datasets lag (same day as `SAMN63084674`) | Empty report | E-utils fallback. Datasets TSV accept is 500; do not use it. |
| 24 | Food / env records Datasets already has (`chicken`+`anal swab`, `wastewater`, `papaya`) | Mapping works when the report is non-empty; status still review | Status `ok`; leave classify to read the text. |

---

## Airtight contract

### Inputs (one token per line, or a TSV id-like column)

Accepted:

- Assembly: `GCA_` / `GCF_` + digits, optional `.version`
- BioSample: `SAMN` / `SAMD` / `SAMEA` / `SAME` + digits
- SRA run/experiment/sample: `S/E/D` + `R` + `R|X|S` + digits (SRR/SRX/SRS and E/D twins)
- PDT: `PDT` + digits, optional `.version`

Rejected (`unknown_id`): BioProject, WGS master (`AAHGAG…`), versioned BioSample (`SAMN….1`), underscore-less GCA, PDT/SRA with junk, blank.

SRA **project** (`SRP`/`ERP`/`DRP`): accepted as a type but **not** resolved to a single BioSample (`sra_study_ambiguous`).

Duplicates: one row per input line.

### Fallback cascade

```text
parse
  → local --metadata (non-blank aliases; versionless GCA/GCF + twin)
    → network:
         assembly: Datasets genome
                   if {} and versioned → retry stem
                   biosample accession / else assembly_no_biosample
         sra run|exp|sample: elink biosample → else runinfo BioSample
                   0 → sra_no_biosample; >1 distinct → multiple_biosamples
         sra study: stop, sra_study_ambiguous
         biosample: Datasets report
                   if empty/no accession → E-utils efetch XML
         pdt: assembly esearch [Assembly Name] exact version
                   optional: current PDTstem.N
                   else pdt_unresolved
    → map named attributes only
    → status
```

Local always wins over network for that ID. Local does not call NCBI to “improve” empty cells.

### Which NCBI field wins

| Output | Winner | Never |
| --- | --- | --- |
| `biosample` | Linked accession from genome / elink / runinfo / assembly summary / the query itself | Empty Datasets `accession` overwriting a known ID |
| `assembly` | Datasets `accession` (current if Datasets substituted) | Inventing GCA from PDT without an assembly hit |
| `sra` / `pdt` | The normalized query (and current PDT name if version-followed) | |
| `organism` | BioSample `organism_name` | Using it as host |
| `host` | Harmonized `host` only, after unknown-token blanking | organism, tax name, cell_line, title |
| `isolation_source` | Harmonized `isolation_source` only | strain, package, isolate name, `ifsac_category`, food_origin, PDT label |
| `sample_type` | `sample_type` else `source_type` else `tissue` | Treating `pure culture` / `whole organism` as the isolation source |
| `strain` / `isolate` / `geo_loc_name` / `collection_date` / `host_disease` | Same-named harmonized attributes | |

If Datasets and E-utils disagree on a filled host/source, **E-utils XML wins** when Datasets was empty; if both are filled and differ, keep Datasets and set `review` / `ncbi_field_conflict` (rare; not seen in this probe set). Do not merge strings.

### `ok` vs `review`

**`ok`** when any of:

- Local metadata keyed to this ID (even if host/source empty).
- Network identity succeeded and `biosample` is a single non-empty accession.

**`review`** when:

| Reason | When |
| --- | --- |
| `unknown_id` | Parser reject |
| `missing_metadata` | No local row and network off |
| `ncbi_error:…` | HTTP / transport |
| `assembly_not_found` | Genome report missing/`{}` after stem retry |
| `assembly_no_biosample` | Assembly exists, no BioSample |
| `sra_no_biosample` | Run/experiment/sample has no BioSample link |
| `sra_study_ambiguous` | SRP/ERP/DRP |
| `multiple_biosamples` | Two or more distinct BioSamples |
| `biosample_not_found` | Both Datasets and E-utils failed to return a record |
| `pdt_unresolved` | No assembly-name hit (after optional current version) |
| `pdt_version_updated` | Optional: resolved via a newer PDT version — may be `ok` with this reason, or `review` if you want a human to see the version move. Pick one and keep it. Recommendation: `ok` + reason `pdt_version_updated` so status stays an identity bit. |

Empty isolation_source, empty host, host/source that classify will later conflict, package-only records: **not** resolve-review.

### Must never be invented

- `host` from `organism` (including `Homo sapiens` on human cell-line / HapMap BioSamples).
- `isolation_source` from strain, substrain, package (`Generic.1.0`, `Pathogen.cl`, `Pathogen.env`, `Microbe.1.0`), isolate name, title, or PDT stem.
- Laboratory / Food / Environmental as a resolve field.
- A BioSample taken from the first of 230 SRA study hits.
- A BioSample taken from a nucleotide All Fields hit that was not confirmed as assembly name.
- Lexicon or override logic inside `resolve`.
- Filling blank cells from the other GCA/GCF twin’s *attributes* unless that twin’s BioSample is the same accession (identity join only).

---

## Smallest code change (do not implement)

Four edits meet the contract for GCA/GCF/BioSample/SRA. PDT and parse/join extras are listed after; they are still small, but they are not required to stop the “always review / blank SAMN” failures.

### Must (identity + status + Datasets hole)

1. **Status.** Stop treating blank-row `review` as final. After `_fetch`, set `ok` when `biosample` is non-empty and `_fetch` did not set `review`. Do not require host or isolation_source. File: `resolve.py` (`_blank_row` / post-check around lines 61–68).

2. **Empty Datasets report.** In `biosample_report` / `genome_report`, treat `not reports or reports[0] == {}` or missing `accession` as `None`. In `_fetch`, never `update` from `attributes_from_biosample` when the report is empty. Do not write `biosample=""`. File: `ncbi.py`, `resolve.py`.

3. **E-utils XML fallback.** If Datasets BioSample is `None`, `esearch db=biosample` + `efetch retmode=xml` and map `harmonized_name` attributes (same output columns). This is the only way `SAMN63084674` (and other Datasets-lag records) get host/source. File: `ncbi.py`.

4. **Skip blank aliases.** `get_field` / `first_present`: first alias with a **non-empty** value. File: `io_tsv.py`. Fixes local resolve + classify shadowing.

### Should (still small; closes PDT / join / SRA study)

5. **PDT via assembly name.** Replace `pdt_to_biosample`’s `db=biosample` search with `esearch db=assembly term="{id}[Assembly Name]"`, then `esummary.biosampleaccn`. Optional: if exact version misses, search `PDTstem` and accept a single current `PDTstem.N`. Do not search nucleotide All Fields. File: `ncbi.py`.

6. **Local GCA/GCF + versionless keys.** When indexing/looking up assemblies, add the stem (`GCA_000005845`) and the GC[AF] twin. File: `resolve.py`.

7. **SRA parse + study.** Allow `S` in the third letter (`SRS`/`ERS`/`DRS`). If the third letter is `P`, do not `retmax=1` resolve; set `sra_study_ambiguous`. File: `ids.py`, `resolve.py`.

8. **Map `sample_type` attr; blank unknown tokens** using the existing `UNKNOWN_EXACT` list when writing resolve fields. File: `ncbi.py` (and maybe `text.py`).

### Do not add in this slice

- Pathogen Detection BigQuery or FTP dumps.
- Datasets TSV accept (500).
- Copying organism → host.
- Inferring Laboratory from `strain=K-12`.
- Deduplicating input IDs.
- New output columns unless you want `paired_accession` / `current_accession` (optional, not required for airtight status).

### Tests that would lock the contract (also not written here)

- Fake Datasets `{"reports":[{}]}` plus fake E-utils XML for `SAMN63084674` → biosample kept, host/source filled, `ok`.
- Assembly fetch with BioSample and empty host/source → `ok` / `ncbi_assembly`.
- `organism=Homo sapiens`, no host attr → host empty.
- host attr `missing` → host empty, isolation_source kept.
- Local TSV empty `isolation_source` + filled `source` → source used.
- `GCA_*` query vs local `GCF_*` same digits/version → local hit.
- PDT assembly-name fixture → BioSample; PDT with 0 assembly hits → `pdt_unresolved`.
- `SRP…` → review, not a single BioSample.
- Two input lines of the same SAMN → two rows.

---

## What “airtight” does not fix

`classify` will still send K-12, cell-line SRA, HapMap, and any BioSample with no isolation_source to Unclassified / `missing_source`. That is correct. Airtight resolve makes that queue honest: `ok` + empty source means “NCBI had nothing,” not “lookup failed.”

PDT IDs that never appear as an assembly name (retired stems with no current `.N`, or isolates only in Isolates Browser / BigQuery) stay `pdt_unresolved` until the user passes `--metadata`. That is the honest ceiling without a PD-specific client.
