# summarize and plot

Two quick looks at a classified table: `summarize` counts the classes and the evidence behind them, and `plot` draws the counts.

```bash
onehealth summarize work/classified.tsv -o work/summary.tsv
onehealth plot work/classified.tsv -o work/counts.png
```

Both read the output of [`classify`](classify.md). For the population analyses papers report, use [`report`](report.md).

## summarize

### Options

| Option | What it does |
|---|---|
| `-o`, `--out` | Output TSV (required). |

### Output

One row per count, with the columns `section`, `key`, `count`, `l1_set` and `mixed`. The sections, in order:

| Section | Counts |
|---|---|
| `dataset` | `rows`, `unique_isolates` and `rows_sharing_a_biosample` |
| `l1` | Rows per L1 class |
| `paper_rollup` | Rows per rollup class |
| `healthcare_environment` | Rows per `yes`, `no` and empty |
| `classify_status` | Rows per `ok` and `review` |
| `classify_reason` | Rows per reason |
| `source_evidence` | What the classes rest on (below) |
| `cluster_overlap` | Only if the table has a `cluster_id`, `cluster` or `poppunk_cluster` column (below) |
| `plasmid_overlap` | Only if the table has a `plasmid_id`, `plasmid` or `plasmid_replicon` column (below) |

Counts are per row. Rows that share a BioSample are the same isolate, so the `dataset` section says how many isolates there are. An isolate is its BioSample; a row without one is matched by its ID, assembly, SRA or PDT accession, or strain or sample name, and a row with no identifier at all counts as its own isolate.

### source_evidence

`classify_status = ok` does not mean the isolation source was known: a row classified from its host or `source_type` alone, or by an override, is `ok` with an empty isolation source. `source_evidence` splits the rows so you can tell:

| Key | Rows |
|---|---|
| `ok_with_source` | Classified, with an isolation source |
| `ok_empty_source` | Classified, with no isolation source |
| `unmatched` | In review because no rule matched |
| `conflict` | In review because of any conflict: host and source, rules, or overrides |
| any other review reason | In review for that reason, for example `missing_source` |

The first four are always listed, even at zero.

### Cluster and plasmid overlap

With a cluster or plasmid column, each cluster or plasmid gets one row: `count` is the number of rows in it, `l1_set` the L1 classes present, and `mixed` is `yes` when there is more than one. `Unclassified` is not a source, so it never makes a cluster mixed.

A shared cluster or plasmid is not evidence of transmission on its own. For linkage by SNP distance, use [`report`](report.md#genomic-linkage).

## plot

### Options

| Option | What it does |
|---|---|
| `-o`, `--out` | Output image (required). The extension sets the format: `.png`, `.svg` or `.pdf`. |

### Output

Two bar charts side by side: rows per L1 class and rows per paper rollup class. The classes are always in the same order, and classes with no rows are drawn at zero, so plots of different datasets line up.
