# resolve

Turns a list of NCBI IDs into one row of BioSample metadata per ID: host, isolation source, collection date, location and the other fields `classify` and `report` use. IDs of different kinds can be mixed in one list.

```bash
onehealth resolve ids.txt -o work/resolved.tsv --email you@example.com
```

`resolve` reads metadata only. It never downloads sequences or reads.

## The ID list

One ID per line, or a TSV with a header whose ID column is `query_id`, `id`, `accession`, `biosample`, `assembly`, `sra`, `pdt` or a similar name (otherwise its first column is used). Blank lines and lines starting with `#` are ignored. Case does not matter. No accession contains a comma or semicolon, so a line with one, as in a CSV file or several IDs on one line, stops with an error naming the line.

| ID type | Examples | Notes |
|---|---|---|
| BioSample | `SAMN02604091`, `SAMEA3375632`, `SAMD00020338` | |
| Assembly | `GCA_000005845.2`, `GCF_000005845`, `gcf_000005845.2` | An unversioned accession gets the current version. An outdated version falls back to the current one. |
| SRA run, experiment or sample | `SRR390728`, `ERR1217372`, `DRR000001`, `SRX…`, `SRS…` | Study accessions (`SRP`, `ERP`, `DRP`) cover many samples and go to review as `sra_study_ambiguous`. |
| Pathogen Detection target | `PDT000046563.2`, `PDT000046563` | Pass `--pd-organism`; see below. |

Anything else, such as a BioProject (`PRJNA…`), goes to review as `unknown_id`.

## Options

| Option | What it does |
|---|---|
| `-o`, `--out` | Output TSV (required). |
| `--metadata` | A local metadata TSV to join first. IDs found there are not looked up at NCBI. It needs an ID column; see [Input files](../guides/input-files.md#identifier-columns). |
| `--no-network` | Use only `--metadata`; IDs not found there go to review as `missing_metadata`. |
| `--email` | Your contact address, sent to NCBI E-utilities as NCBI asks. |
| `--api-key` | An [NCBI API key](https://account.ncbi.nlm.nih.gov/settings/). Raises the request rate from 3 to 10 per second. Without the option, the key is taken from the `NCBI_API_KEY` environment variable, which keeps it out of your shell history. |
| `--pd-organism` | Pathogen Detection organism group(s) of your PDT IDs, comma-separated, for example `Klebsiella` or `Escherichia_coli_Shigella`. |
| `--review` | Where to write the review queue (default `<out>.review.tsv`). |

## How IDs are looked up

1. **Local metadata.** With `--metadata`, each ID is matched against every ID column of that table. An assembly also matches its GCA/GCF twin and its unversioned form. Two different rows for the same ID are a conflict, not "last row wins". Every other column of the matched row (sequence type, patient ID, dates, ...) is carried into the output.
2. **NCBI**, for IDs not found locally:
    - **Assemblies** are looked up in NCBI Datasets, 100 per request, to find their BioSample.
    - **SRA accessions** are looked up one at a time through E-utilities.
    - **PDTs** are looked up in NCBI Pathogen Detection's metadata table when you pass `--pd-organism`; see below.
    - **BioSample reports** for all IDs are fetched from NCBI Datasets, 100 per request. A record Datasets returns empty is read from E-utilities instead.

Anything a batch does not answer cleanly is looked up on its own. If NCBI sends an error or a garbled reply for one ID, after a few retries that ID goes to review as `ncbi_error:<message>` and the run continues. The same goes for a reply that is well-formed but not shaped as NCBI documents it (`ncbi_error:unexpected reply ...`). Your API key is never written into the table.

A BioSample read from E-utilities is searched for by its accession, and the record that comes back is used only if it carries that accession. A search of all fields can land on another sample that merely mentions it.

**Stopping a run.** Ctrl-C stops the lookups, not the run: what NCBI has already answered is written to the output, the IDs not yet looked up are marked `interrupted` and listed in the review queue, and `resolve` exits with status 130. Before, a run stopped part-way wrote nothing.

Standard missing-value words in NCBI's fields (`missing`, `not collected`, `not provided`, ...) are left blank in the output.

### PDT IDs

NCBI names an assembly after its PDT only when it assembled the genome itself, so few PDTs can be found by assembly name. Pass the PDT's organism group:

```bash
onehealth resolve pdts.txt -o work/resolved.tsv --email you@example.com --pd-organism Klebsiella
```

This streams the group's current metadata table from NCBI's FTP site, which lists every PDT with its BioSample and assembly. The table is not saved: it is about 200 MB for Klebsiella and 900 MB for Salmonella, and reading stops as soon as every PDT has been found. An outdated or unversioned PDT gets its current version, and the `pdt` and `assembly` columns are filled in from the table.

- **Group names** are as listed at <https://ftp.ncbi.nlm.nih.gov/pathogen/Results/>, in any case. A misspelt name stops the run with suggestions.
- **PDTs the table does not list** still get the assembly-name lookup.

## Output

`resolve` writes one row per input ID, in the input order:

`query_id`, `id_type`, `biosample`, `assembly`, `sra`, `pdt`, `organism`, `host`, `isolation_source`, `host_disease`, `sample_type`, `source_type`, `package`, `geo_loc_name`, `collection_date`, `strain`, `isolate`, `resolve_status`, `resolve_reason`, then any extra columns from `--metadata`.

`source_type` is NCBI's One Health category of the record (`human`, `animal`, `food`, `environmental`), where the submitter gave one. `classify` uses it as the class; see [Definitions](../guides/definitions.md#1-a-records-own-source_type-is-its-class). `sample_type` holds the same value when there is one, and otherwise the `tissue` or NCBI's material type, which are only text.

`package` is the BioSample package the record was submitted under, such as `Pathogen.cl.1.0` (clinical or host-associated) or `Pathogen.env.1.0` (environmental, food or other). [`classify`](classify.md#specimens-follow-the-host) uses it for specimens that have no host.

Rows with `resolve_status` = `review` are also written to the review queue. Duplicate IDs are kept as separate rows.

### resolve_status and resolve_reason

| Status | Reason | Meaning |
|---|---|---|
| ok | `local_metadata` | Found in `--metadata` |
| ok | `ncbi` | BioSample found at NCBI |
| ok | `ncbi_assembly` | Assembly found, then its BioSample |
| ok | `ncbi_sra` | SRA accession found, then its BioSample |
| ok | `ncbi_pdt` | PDT found, then its BioSample |
| review | `local_metadata_conflict` | Two different rows in `--metadata` for this ID |
| review | `unknown_id` | Not an ID type `resolve` knows |
| review | `sra_study_ambiguous` | An SRA study covers many samples |
| review | `missing_metadata` | Not in `--metadata`, and `--no-network` was given |
| review | `assembly_not_found` | NCBI has no such assembly |
| review | `assembly_no_biosample` | The assembly has no BioSample |
| review | `sra_no_biosample` | The SRA accession has no BioSample |
| review | `pdt_unresolved` | The PDT was found neither in the Pathogen Detection table nor by assembly name |
| review | `multiple_biosamples` | The ID points to several BioSamples, listed in `biosample` |
| review | `biosample_not_found` | NCBI has no record for the BioSample |
| review | `ncbi_incomplete` | NCBI's answer had no single BioSample |
| review | `ncbi_error:<message>` | NCBI failed for this ID; the message says how |
| review | `interrupted` | The run was stopped (Ctrl-C) before this ID was looked up |

## Rate limits

NCBI allows 3 requests a second without an API key and 10 with one; `resolve` keeps to that. Because BioSample and assembly reports come 100 per request, a list of 917 mixed IDs took 30 seconds in testing with `--pd-organism Klebsiella`. Without that option it took 83 seconds, because each PDT then needs searches of its own. SRA accessions are looked up one at a time and are the other slow part.
