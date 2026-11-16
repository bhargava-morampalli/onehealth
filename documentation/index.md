# onehealth

`onehealth` is a command-line tool for One Health genomic studies of bacteria. It does the step every study does by hand: it turns the messy isolation-source metadata of public isolates into consistent One Health categories. It then runs the population analyses those studies publish and writes them as a report.

```text
NCBI IDs or local metadata ──► resolve ──► classify ──► report ──► report.html, figures/, tables/
                                              │
                                              ├──► summarize, plot
                                              └──► validate (against labels you trust)
```

!!! warning "The classes follow NCBI"
    Where NCBI's records show a convention, `onehealth` follows it: a record's own `source_type` is its class, farm settings (manure, a barn) are Environmental, `chicken` on its own with no host is Food, and a specimen with no host is Human. These definitions decide how many isolates are in each group. Read [Definitions](guides/definitions.md) before you compare sources.

## What it does

| Command | What it does |
|---|---|
| [`resolve`](commands/resolve.md) | Turns a mixed list of NCBI IDs (BioSample, GCA/GCF assembly, SRA run, Pathogen Detection PDT) into one row of BioSample metadata per ID |
| [`classify`](commands/classify.md) | Assigns each isolate a One Health class (Human, Animal, Environmental, Food, Laboratory) from its host and isolation source, and sends anything unclear to a review queue |
| [`summarize`, `plot`](commands/summarize-plot.md) | Count the classes and the evidence behind them |
| [`validate`](commands/validate.md) | Compares the classes with labels you trust: coverage, accuracy, confusion matrix, every disagreement |
| [`report`](commands/report.md) | Lineage diversity and sharing, resistance and virulence by source, and genomic linkage between sources, as one self-contained HTML report with figures and tables |

## What it does not do

`onehealth` does not download genomes, assemble, type or call SNPs. Established tools do that well. `report` reads their output: Kleborate, mlst or PopPUNK for lineages, AMRFinderPlus for genes, and snp-dists for SNP distances.

## How it decides

- **Nothing is guessed silently.** Missing, unknown or conflicting evidence goes to a review queue (`<output>.review.tsv`) with a reason, never to a default class. `unknown` is never treated as Environmental.
- **The rules are a versioned table.** Classification comes from a [lexicon](guides/lexicon.md) of phrases, shipped as a YAML file. Every classified row records the lexicon version, and older versions still ship, so results can be reproduced.
- **Your corrections win.** A curator [overrides file](commands/classify.md#overrides) always takes precedence over the lexicon.
- **Every number is in a table.** Each figure in the report has its underlying TSV, and each statistic has its method and reference in the report's Methods section.

## Where to start

- New here: [install](install.md), then follow the [tutorial](tutorial.md). It runs on example data shipped with the tool.
- Bringing your own table: read [Input files](guides/input-files.md) for what `onehealth` accepts.
- Checking it on your data first: see [Benchmarks](benchmarks.md) and [`validate`](commands/validate.md).
