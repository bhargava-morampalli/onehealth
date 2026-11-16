# Tutorial

This walks through the whole pipeline on data that ships with `onehealth`, so nothing needs downloading. Run the commands from the root of the repository. They write into `work/`, which git ignores.

## Part 1: from IDs to classes

`examples/ids.txt` holds four IDs of different kinds: a BioSample, an assembly, an SRA run and a Pathogen Detection PDT. `examples/metadata.tsv` holds their metadata, so this part runs offline.

### Resolve the IDs

```bash
onehealth resolve examples/ids.txt -o work/resolved.tsv --metadata examples/metadata.tsv --no-network
```

```text
resolve: 4 rows (4 ok, 0 review)
```

Each ID is matched to its row in the metadata, whatever kind of ID it is. `work/resolved.tsv` has one row per ID:

| query_id | biosample | isolation_source | resolve_status | resolve_reason |
|---|---|---|---|---|
| SAMN00000001 | SAMN00000001 | wound drainage | ok | local_metadata |
| GCA_000000001.1 | SAMN00000002 | hospital drain | ok | local_metadata |
| SRR0000001 | SAMN00000003 | cloacal swab | ok | local_metadata |
| PDT000000001.1 | SAMN00000001 | wound drainage | ok | local_metadata |

Without `--metadata` and `--no-network`, the same command looks each ID up at NCBI instead; see [`resolve`](commands/resolve.md).

### Classify them

```bash
onehealth classify work/resolved.tsv -o work/classified.tsv --overrides examples/overrides.tsv
```

```text
classify: 4 rows (4 ok, 0 review)
lexicon v8 (seed v8)
```

`classify` adds the class columns:

| query_id | isolation_source | l1 | paper_rollup | classify_reason |
|---|---|---|---|---|
| SAMN00000001 | wound drainage | Human | Human | wound-drainage+human-host |
| GCA_000000001.1 | hospital drain | Environmental | Environmental | hospital-drain |
| SRR0000001 | cloacal swab | Animal | Animal | override |
| PDT000000001.1 | wound drainage | Human | Human | wound-drainage+human-host |

`classify_reason` says why: `wound-drainage+human-host` means the specimen rule matched and the human host decided Human. `override` means the row in `examples/overrides.tsv` decided. Rows that cannot be classified would be in `work/classified.review.tsv`, and each distinct source text among them once, with its row count, in `work/classified.review_groups.tsv`.

### Count and plot

```bash
onehealth summarize work/classified.tsv -o work/summary.tsv
onehealth plot work/classified.tsv -o work/counts.png
```

The summary starts with the dataset itself:

| section | key | count |
|---|---|---|
| dataset | rows | 4 |
| dataset | unique_isolates | 3 |
| dataset | rows_sharing_a_biosample | 1 |

There are 4 rows but only 3 isolates: the BioSample and the PDT are the same isolate, `SAMN00000001`. Counts are per row, so the summary says so.

## Part 2: the report

`examples/demo/make_demo.py` writes a synthetic dataset of 316 isolates: metadata, AMRFinderPlus-style gene calls and a SNP distance matrix.

```bash
python examples/demo/make_demo.py work/demo
onehealth classify work/demo/metadata.tsv -o work/demo/classified.tsv
onehealth report work/demo/classified.tsv -o work/demo/report \
    --features work/demo/amrfinder.tsv --distances work/demo/snp_dists.tsv \
    --note "Synthetic demonstration data"
```

```text
classify: 316 rows (307 ok, 9 review)
wrote work/demo/report/report.html
```

Open `work/demo/report/report.html` in a browser. It is one self-contained file. Next to it:

- `figures/` has every figure as SVG and PNG, for papers.
- `tables/` has every number behind every figure as TSV, from `01_sampling.tsv` to `11_cross_source_pairs.tsv`.

With 9,999 permutations per test (the default) the report takes a few seconds here. [`report`](commands/report.md) describes each section and its options.

!!! note
    The demonstration data is synthetic. It shows how the report looks and behaves, not real biology.

## Next steps

- Run on your own IDs: [`resolve`](commands/resolve.md).
- Run on your own metadata table: [Input files](guides/input-files.md), then [`classify`](commands/classify.md).
- Check the classes against labels you trust before relying on them: [`validate`](commands/validate.md).
