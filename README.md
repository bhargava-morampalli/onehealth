# onehealth

Python CLI that resolves mixed public NCBI isolate IDs and classifies isolation-source metadata into One Health categories.

v1 does **not** download genomes, assemble, run SKA/PopPUNK/MOB-suite, or start a database. Cluster or plasmid overlap is reported only when those columns are already in the input table.

## Install

```bash
cd ~/onehealth-pipeline
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Or with uv:

```bash
cd ~/onehealth-pipeline
uv sync --extra dev
uv run onehealth --help
```

## Commands

| Command | What it does |
| --- | --- |
| `onehealth resolve` | Map GCA/GCF, BioSample, SRA, and optional PDT IDs to BioSample metadata |
| `onehealth classify` | Assign L1, `healthcare_environment`, and a paper rollup |
| `onehealth summarize` | Count classes; add cluster/plasmid overlap only if those columns exist |
| `onehealth plot` | Write L1 and paper-rollup bar charts |

```bash
onehealth resolve examples/ids.txt -o work/resolved.tsv --metadata examples/metadata.tsv --no-network
onehealth classify work/resolved.tsv -o work/classified.tsv --overrides examples/overrides.tsv
onehealth summarize work/classified.tsv -o work/summary.tsv
onehealth plot work/classified.tsv -o work/counts.png
```

Live NCBI lookup (metadata only; no FASTA/SRA download):

```bash
onehealth resolve examples/ids.txt -o work/resolved.tsv --email you@example.com
```

Rows that cannot be resolved or classified are written next to the output as `<stem>.review.tsv`.

## Classification

Default **L1**: `Human | Animal | Environmental | Food | Laboratory | Unclassified`.

`healthcare_environment` is `yes` or `no`. Hospital fixtures and bare `clinical` / `clinical environment` labels are healthcare-environment, not Human.

**Paper rollup** (Human / Environmental / Animal):

| L1 | healthcare_environment | paper_rollup |
| --- | --- | --- |
| Human | no | Human |
| Animal | * | Animal |
| Environmental, Food, or Laboratory | * | Environmental |
| any L1 | yes | Environmental |
| Unclassified | * | Unclassified |

Unknown, empty, or conflicting evidence stays `Unclassified` and goes to the review queue. `unknown` is never mapped to Environmental.

Matching is token-based, so `wound drainage` is Human and `hospital drain` is healthcare-environment. The word `drain` does not match `drainage`.

A complete user override TSV always wins. An override row with a missing/empty/invalid `l1`, or two disagreeing rows for the same ID, is queued for review. The lexicon is not applied silently in those cases.

Override columns: `id` (or BioSample/assembly/SRA/PDT), `l1`, optional `healthcare_environment`, `paper_rollup`, `note`.

## Lexicon

The seeded lexicon is versioned at `src/onehealth/data/lexicon/v1.yaml`. Pass `--lexicon` to use another YAML. This is a small maintained table, not a wrapper of an external R notebook.

## Tests

```bash
pip install -e ".[dev]"
pytest
```
