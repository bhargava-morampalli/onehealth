# onehealth

Python CLI that resolves mixed public NCBI isolate IDs and classifies isolation-source metadata into One Health categories.

v1 does **not** download genomes, assemble, run SKA/PopPUNK/MOB-suite, or start a database. Cluster or plasmid overlap is reported only when those columns are already in the input table.

## Install

```bash
git clone https://github.com/bhargava-morampalli/onehealth.git
cd onehealth
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Or with uv:

```bash
cd onehealth
uv sync --extra dev
uv run onehealth --help
```

## Commands

| Command | What it does |
| --- | --- |
| `onehealth resolve` | Map GCA/GCF, BioSample, SRA, and optional PDT IDs to BioSample metadata |
| `onehealth classify` | Assign L1, `healthcare_environment`, and a paper rollup |
| `onehealth summarize` | Count classes and source evidence (`ok` with a source vs `ok` with an empty source vs each review reason); rows vs unique isolates; cluster/plasmid overlap only if those columns exist |
| `onehealth plot` | Write L1 and paper-rollup bar charts |
| `onehealth validate` | Compare the classes with a curated label table: coverage, accuracy, confusion matrix and every disagreement |
| `onehealth report` | Standard One Health analyses as an HTML report, figures (SVG/PNG) and tables (TSV) |

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

NCBI allows 3 requests a second, or 10 with a free API key (create one in your [NCBI account settings](https://account.ncbi.nlm.nih.gov/settings/)); pass it with `--api-key` for large ID lists. If NCBI sends back an error or a garbled reply for an ID, that ID goes to review as `ncbi_error` and the run continues.

Rows that cannot be resolved or classified are written next to the output as `<stem>.review.tsv`.

With `--metadata`, every extra local column (ST, `host_subject_id`, dates, ...) is carried into the resolved table. You can also run `classify` directly on your own metadata TSV.

## Validate

Check the classifier against labels you trust before relying on it for a new dataset. The truth file needs an ID column and a label column naming a class (Human, Animal, Environmental, Unclassified; or the L1 classes with `--level l1`).

```bash
onehealth validate work/classified.tsv --truth curated.tsv --truth-column truth -o work/validation.tsv
```

It writes a summary (coverage, accuracy, agreement), per-class recall and precision, a confusion matrix, and a list of every disagreement with the host, source and rule involved.

Against the manually curated source labels of 4,351 public BioSamples from the IMP carbapenemase study (Vezina, Morampalli et al. 2025), lexicon v2 classifies 84.5% of isolates, and 99.2% of those agree with the curation. Most of the rest have no isolation source or host in their BioSample. The benchmark lives in `tests/data/imp_carbapenemase_sources/`.

## Report

`onehealth report` runs the analyses that recent One Health genomic studies publish, and writes them as one self-contained `report.html` plus `figures/` (SVG and PNG) and `tables/` (TSV). [ROADMAP.md](ROADMAP.md) lists which studies each analysis comes from and the statistics used.

```bash
onehealth report work/classified.tsv -o work/report \
    --features amrfinder.tsv --distances snp_dists.tsv --threshold 21
```

| Input | Enables | Format |
| --- | --- | --- |
| classified TSV (required) | Sampling, lineage diversity, rarefaction, lineage overlap and enrichment | Output of `classify`; lineage from a column named `st`, `mlst`, `lineage`, `cluster_id`, ... or `--lineages` |
| `--lineages` | Lineages from a typing tool | Kleborate v2/v3, mlst, PopPUNK `clusters.csv` |
| `--features` | Resistance / virulence / plasmid prevalence by source | Long table: isolate, feature, optional class. AMRFinderPlus run with `--name` fits as is. Only isolates listed count as genotyped. |
| `--distances` | Genomic linkage between sources | snp-dists square matrix or `-m` long format |

Add `--stratify-by year`, `country` or `year,country` to shuffle source labels only among isolates from the same collection year and country in every permutation test. Sources sampled at different times or places then no longer look like barriers to sharing. Historical INSDC country names (for example "Czech Republic") are mapped to their current names first.

Add `--dereplicate` (with `--distances`) to run the lineage and feature analyses on one isolate per genomic cluster per source. Isolates from the same source within `--threshold` SNPs collapse to one representative, so outbreak clusters do not inflate the tests; isolates from different sources are never merged. Sampling and linkage sections always use every isolate.

The SNP threshold defaults to 21, the published *K. pneumoniae* value. Choose one for your organism and pipeline; the report always shows a sweep from 0 to 100 SNPs. Tests use 9,999 permutations with a fixed `--seed`, so results reproduce exactly.

Try it on synthetic data:

```bash
python examples/demo/make_demo.py work/demo
onehealth classify work/demo/metadata.tsv -o work/demo/classified.tsv
onehealth report work/demo/classified.tsv -o work/demo/report \
    --features work/demo/amrfinder.tsv --distances work/demo/snp_dists.tsv --note "Synthetic demonstration data"
```

## Classification

Default **L1**: `Human | Animal | Environmental | Food | Laboratory | Unclassified`.

`healthcare_environment` is `yes` or `no` for classified rows and empty for `Unclassified` rows. Hospital fixtures and `clinical environment` labels are healthcare-environment, not Human.

**Paper rollup** (Human / Environmental / Animal):

| L1 | healthcare_environment | paper_rollup |
| --- | --- | --- |
| Human | no | Human |
| Animal | * | Animal |
| Environmental, Food, or Laboratory | * | Environmental |
| any L1 | yes | Environmental |
| Unclassified | * | Unclassified |

Unknown, empty, or conflicting evidence stays `Unclassified` and goes to the review queue. `unknown` is never mapped to Environmental.

**Specimens follow the host.** Body sites and specimens (`blood`, `urine`, `faeces`, `wound`, `tissue`) occur in any host. A recognised animal host or animal word makes them Animal (`pig faeces`, `urine` + `Canis lupus familiaris`). A human host makes them Human. Both at once is a `host_source_conflict`. A host the lexicon does not recognise, such as a plant, goes to review as `host_unrecognised`. With no host at all they default to Human. The reason then ends in `+no_host`, so you can count those calls.

Bare `clinical` / `clinical isolate` is a patient sample, as in the NCBI `Pathogen.cl` package ("clinical or host-associated"). It is not an environmental swab, and it follows the host like any specimen.

NCBI `source_type` (`human`, `animal`, `food`, `environmental`) is read into `sample_type` before `tissue` and the NCBI material type (`whole organism`, `pure culture`). The material type is never treated as source evidence.

KlebNET `Animal - food producing` is Animal, not Food. `animal feed` is not retail Food and has no L1 here, so it goes to review (`animal-feed`). A taxon-only `isolation_source` (`Escherichia coli`, `Bos taurus`) or an unknown marker (`unknown`, `missing: lab stock`) is ignored. Taxon-only is a guess from word shape, so a lexicon rule written for source text (one that does not read `host`) always wins over it: add `bedpan washer` to such a rule and it is classified. The host or `source_type` then decides, and with nothing else the row goes to review (`taxon_only_source` / `unknown_source`). Bare `sink` or `faucet` is not a hospital fixture; `hospital sink` / `hospital faucet` are.

Matching is token-based, so `wound drainage` is Human and `hospital drain` is healthcare-environment. The word `drain` does not match `drainage`.

A complete user override TSV always wins. An override row with a missing/empty/invalid `l1` or `paper_rollup`, or two disagreeing rows for the same ID, is queued for review. The lexicon is not applied silently in those cases. Class names are case-insensitive. An override on a `GCF_` accession also applies to its `GCA_` twin, and an unversioned accession covers every version.

Override columns: `id` (or BioSample/assembly/SRA/PDT), `l1`, optional `healthcare_environment`, `paper_rollup`, `note`.

**Counting.** Summary and plot counts are per input row. If the same isolate was listed twice (for example by BioSample and by PDT), the summary's `dataset` section reports `rows_sharing_a_biosample`. Remove duplicates before reporting isolate counts. A shared cluster or plasmid (`mixed=yes`) is a lead, not evidence of transmission. `Unclassified` never makes a cluster mixed.

## Lexicon

The seeded lexicon is versioned at `src/onehealth/data/lexicon/v2.yaml`, and every classified row records `lexicon_version`. Pass `--lexicon` to use another YAML. The file header explains `priority`, `exact` vs `phrases`, `host_dependent: true` (specimen rules), and `l1: Unclassified` (send to review). This is a small maintained table, not a wrapper of an external R notebook.

What changed from v1 and why, with literature, is in [REVIEW.md](REVIEW.md).

## Tests

```bash
pip install -e ".[dev]"
pytest
```

Besides small hand-built cases, the tests run on real public data: 307 *Klebsiella* isolates from NCBI Pathogen Detection with NCBI's own metadata, SNP calls and AMRFinderPlus results. Provenance and the rebuild script are in `tests/data/ncbi_pd_klebsiella/`.
