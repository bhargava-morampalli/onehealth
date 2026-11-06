# Input files

`onehealth` reads tab-separated tables. Any of these work as they are:

- the output of [`resolve`](../commands/resolve.md);
- a metadata table of your own, with the columns below;
- NCBI Pathogen Detection metadata tables (`*.metadata.tsv`).

## Format

- **Tab-separated**, with a header row. A comma- or semicolon-separated file, such as Excel's CSV export, stops with an error that says so, rather than being read as one column.
- **UTF-8.** A byte-order mark at the start is fine. Any other encoding stops with an error naming the byte, rather than being guessed. Excel's tab-delimited export may not be UTF-8, so a table with accented or other non-ASCII letters can fail this check; LibreOffice and Google Sheets can save tab-separated UTF-8.
- **One line is one row.** Blank lines are skipped, and lines starting with `#` are comments, even above the header.
- **Quotes:** a cell wrapped in double quotes, as Excel and pandas write a cell that contains a quote, is unwrapped. Any other quote is kept as typed.
- **Cells** are trimmed of surrounding spaces. A row with fewer cells than the header is padded with empty cells; a row with more stops with an error naming the line, unless the extra cells are empty.
- **Column names** are matched in any case, so `Host` and `host` are the same column.

## Columns `classify` reads

| Field | Taken from the first non-empty of |
|---|---|
| isolation_source | `isolation_source`, `isolation source`, `source`, `sample_source` |
| host | `host` |
| sample_type | `source_type`, `tissue`, `sample_type` |
| host_disease | `host_disease`, `host disease` |

A table with none of these columns stops with an error, because every row would go to review; see [`classify`](../commands/classify.md#what-the-classes-are-based-on), which also explains how the fields are used.

## Identifier columns

Identifiers link rows to overrides, truth labels, lineages, features and distances, and tell which rows are the same isolate. A `--metadata` table for `resolve` must have at least one of them.

| Identifier | Columns |
|---|---|
| Your own ID | `query_id`, `id`, `accession`, `isolate_id`, `sample_id` |
| BioSample | `biosample`, `biosample_acc`, `biosample_accession`, `biosample accession`, `assembly biosample accession` |
| Assembly | `assembly`, `asm_acc`, `assembly_accession`, `genome`, `assembly accession` |
| SRA | `sra`, `srr`, `run`, `run_accession` |
| Pathogen Detection target | `pdt`, `target_acc`, `pdt_accession` |

The names with spaces are the ones NCBI Datasets (`dataformat tsv genome`) writes. An assembly also matches its GCA/GCF twin and its unversioned form, so `GCF_000005845` finds `GCA_000005845.2`.

**Rows that share a BioSample are the same isolate**, for example two assemblies of one sample, or a BioSample and a PDT for it. `classify` keeps every row, `summarize` reports how many isolates the rows are, and `report` counts each isolate once. A row without a BioSample is identified by its ID, assembly, SRA or PDT accession, or its `strain`, `isolate`, `sample` or `name` column. A row with no identifier at all counts as its own isolate, but cannot be matched to anything.

## Columns `summarize` and `report` use

| Column | Used by |
|---|---|
| `cluster_id`, `cluster`, `poppunk_cluster` | `summarize`: which classes share each cluster |
| `plasmid_id`, `plasmid`, `plasmid_replicon` | `summarize`: which classes share each plasmid |
| `st`, `mlst_st`, `sequence_type`, `mlst`, `lineage`, `sublineage`, `sl`, `poppunk_cluster`, `cluster_id`, `cluster` | `report`: the lineage |
| `collection_date` | `report`: year of collection |
| `geo_loc_name`, `country`, `geo_loc` | `report`: country, for `--stratify-by country` |
| `host_id`, `host_subject_id`, `patient_id`, `subject_id`, `individual_id`, `animal_id` | `report`: pairs from the same host individual are not counted as sharing |

Every other column is carried through to the output unchanged. [`report`](../commands/report.md#input-files) describes the features, distances and lineage files it reads alongside the table.

## Missing values

These words mean "no information", in any case and with underscores read as spaces. `onehealth` treats them as an empty cell and never as a class: `unknown` is never Environmental.

`-` · `.` · `na` · `n/a` · `n.a` · `n.a.` · `nan` · `#n/a` · `none` · `null` · `nil` · `missing` · `missing data` · `unknown` · `unk` · `unspecified` · `unavailable` · `undetermined` · `no data` · `none provided` · `not applicable` · `not available` · `not collected` · `not determined` · `not known` · `not provided` · `not recorded` · `not specified` · `restricted access`

Also missing:

- any value that starts with `unknown`, `unk` or `missing`, such as `missing: control sample` (the INSDC form) or `unknown environment`;
- `unknown` after only general words (`environment`, `environmental`, `source`, `origin`, `host`, `sample`, `isolate`, `type`, `location`, `site`), such as `source unknown` or `sample type unknown`.

In `isolation_source`, a missing value is ignored and the host or `source_type` decides; with nothing else, the row goes to review as `unknown_source`. Other fields keep their words, so a host of `unknown animal` still says the host is an animal. A bare scientific name (`Escherichia coli`, `Bos taurus`) is ignored the same way, as taxonomy belongs in `host`.
