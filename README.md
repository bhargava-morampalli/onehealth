# onehealth

Python CLI that resolves mixed public NCBI isolate IDs and classifies isolation-source metadata into One Health categories.

v1 does **not** download genomes, assemble, run SKA/PopPUNK/MOB-suite, or start a database. Cluster or plasmid overlap is reported only when those columns are already in the input table.

The full documentation, with a tutorial and a page for each command, is published at <https://bhargava-morampalli.github.io/onehealth/>. Its source is in [`documentation/`](documentation/index.md).

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

NCBI allows 3 requests a second, or 10 with a free API key (create one in your [NCBI account settings](https://account.ncbi.nlm.nih.gov/settings/)); pass it with `--api-key` for large ID lists, or set the `NCBI_API_KEY` environment variable to keep it out of your shell history. BioSample and assembly reports are fetched 100 per request, so a list of thousands of BioSamples or assemblies takes a few dozen requests. If NCBI sends back an error or a garbled reply for an ID, that ID goes to review as `ncbi_error` and the run continues. Ctrl-C stops the lookups but still writes what NCBI had answered; the IDs not reached are marked `interrupted`.

PDT IDs need the Pathogen Detection organism group:

```bash
onehealth resolve pdts.txt -o work/resolved.tsv --email you@example.com --pd-organism Klebsiella
```

This reads the group's current metadata table from NCBI, which lists every PDT with its BioSample and assembly. The table is streamed, not saved: about 200 MB for Klebsiella, 900 MB for Salmonella, and reading stops once every PDT is found. Several groups can be given, comma-separated. Group names are as at <https://ftp.ncbi.nlm.nih.gov/pathogen/Results/>. An outdated PDT version gets the current one. Without `--pd-organism`, a PDT is found only if NCBI named its assembly after it, which few genomes are.

A run that was stopped, was killed, or met an NCBI error is finished with the same command and `--resume`:

```bash
onehealth resolve examples/ids.txt -o work/resolved.tsv --email you@example.com --resume
```

It keeps every ID the earlier run settled and looks up only the rest. While it works, `resolve` saves the table every 30 seconds, so a killed run loses at most that. An ID that was never looked up (`interrupted`, `ncbi_error`) is not a record with no source: `classify` sends it to review as `unresolved:<reason>` and warns, where a record that was fetched and has no source is `missing_source`.

Rows that cannot be resolved or classified are written next to the output as `<stem>.review.tsv`. `classify` also writes `<stem>.review_groups.tsv`: each distinct source text in the review queue once, with the number of rows that carry it and empty columns for a decision. Fill in `l1` for the texts you can place and pass the table back with `--decisions`; each decision then applies to every row with that text. In a sample of 159,045 NCBI records, the 2,626 review rows that say anything hold 694 distinct texts.

With `--metadata`, every extra local column (ST, `host_subject_id`, dates, ...) is carried into the resolved table. You can also run `classify` directly on your own metadata TSV.

Input tables are tab-separated UTF-8 text with a header row; a byte-order mark is fine, and cells may contain quote characters. NCBI Pathogen Detection metadata tables work as downloaded, although their header starts with `#`. A column name that appears twice stops with an error rather than letting the second column replace the first. A comma- or semicolon-separated file, such as Excel's CSV export, stops with an error rather than being misread. A problem with an input file or option is reported in one line, and the command exits with status 1. Changes between versions are listed in [CHANGELOG.md](CHANGELOG.md).

## Validate

Check the classifier against labels you trust before relying on it for a new dataset. The truth file needs an ID column and a label column naming a class (Human, Animal, Environmental, Unclassified; or the L1 classes with `--level l1`).

```bash
onehealth validate work/classified.tsv --truth curated.tsv --truth-column truth -o work/validation.tsv
```

It writes a summary (coverage, accuracy, agreement), per-class recall and precision, a confusion matrix, and a list of every disagreement with the host, source and rule involved.

Against the manually curated source labels of 4,351 public BioSamples from the IMP carbapenemase study (Vezina, Morampalli et al. 2025), the lexicon classifies 84.9% of isolates, and 99.8% of those agree with the curation. Most of the rest have no isolation source or host in their BioSample. The benchmark lives in `tests/data/imp_carbapenemase_sources/`.

That study is almost all patient isolates, and the lexicon was tuned on it. For food, animal and environmental sources there is a second benchmark, `tests/data/ncbi_pd_source_types/`: 166,158 NCBI Pathogen Detection records from 26 organism groups whose submitters filled in NCBI's own `source_type`, which is hidden from the classifier and used as the answer. On the 26,604 records drawn after the current lexicon (v8) was finished, and read once, it classifies 97.8%, and 96.7% of those match NCBI's label. A Human call matches 99% of the time, Food 98% and Animal 97%. Environmental matched 81%, because of one description (`healthcare facility`) that has since been corrected; without it, 94%. Both benchmarks are described in [documentation/benchmarks.md](documentation/benchmarks.md).

## Report

`onehealth report` runs the analyses that recent One Health genomic studies publish, and writes them as one self-contained `report.html` plus `figures/` (SVG and PNG) and `tables/` (TSV). The report's Methods section names the study each analysis follows and the statistics used, with references.

```bash
onehealth report work/classified.tsv -o work/report \
    --features amrfinder.tsv --distances snp_dists.tsv --threshold 21
```

| Input | Enables | Format |
| --- | --- | --- |
| classified TSV (required) | Sampling, lineage diversity, rarefaction, lineage overlap and enrichment | Output of `classify`; lineage from a column named `st`, `mlst`, `lineage`, `cluster_id`, ... or `--lineages` |
| `--lineages` | Lineages from a typing tool | Kleborate v2/v3, mlst (with or without a header row), PopPUNK `clusters.csv` |
| `--features` | Resistance / virulence / plasmid prevalence by source | Long table: isolate, feature, optional class. AMRFinderPlus run with `--name` fits as is. Only isolates listed count as genotyped. |
| `--distances` | Genomic linkage between sources | snp-dists square matrix or `-m` long format |

Add `--stratify-by year`, `country` or `year,country` to shuffle source labels only among isolates from the same collection year and country in every permutation test. Sources sampled at different times or places then no longer look like barriers to sharing. Historical INSDC country names (for example "Czech Republic") are mapped to their current names first.

Add `--dereplicate` (with `--distances`) to run the lineage and feature analyses on one isolate per genomic cluster per source. Isolates from the same source within `--threshold` SNPs collapse to one representative, so outbreak clusters do not inflate the tests; isolates from different sources are never merged. Sampling and linkage sections always use every isolate.

The SNP threshold defaults to 21, the published *K. pneumoniae* value. Choose one for your organism and pipeline; the report always shows a sweep from 0 to 100 SNPs. Tests use 9,999 permutations with a fixed `--seed`, so results reproduce exactly. `--permutations 0` skips them for a quick look: p-values are left blank and the observed-vs-expected section is omitted.

Every statistic in the report is compared with R, the vegan package or scipy, or with the exact answer, each time the tests run. [documentation/guides/statistics.md](documentation/guides/statistics.md) shows the comparison.

Try it on synthetic data:

```bash
python examples/demo/make_demo.py work/demo
onehealth classify work/demo/metadata.tsv -o work/demo/classified.tsv
onehealth report work/demo/classified.tsv -o work/demo/report \
    --features work/demo/amrfinder.tsv --distances work/demo/snp_dists.tsv --note "Synthetic demonstration data"
```

## Classification

> **The classes follow NCBI.** A One Health class is a definition, and definitions decide how many isolates are in each group. Where NCBI's records show a convention, `onehealth` follows it:
>
> 1. **A record's own `source_type` is its class.** Where the submitter gave NCBI the One Health category (`human`, `animal`, `food`, `environmental`), it is used as it is. The lexicon reads free text and is not allowed to disagree with it. A third of NCBI records have one.
> 2. **Farm settings are Environmental.** Manure, litter, a barn and a farm were Animal until lexicon v7. NCBI submitters record manure as environmental in 646 records of 747, and a farm with no host in 141 of 187. A farm *with* a host is the host: `pig farm`, host *Sus scrofa*, is Animal.
> 3. **A meat species with no host is Food.** `chicken`, `turkey`, `fish`, `poultry`, `salmon`, `beef` or `lamb` in the text, with no host, is the meat: food in 3,079 NCBI records of 3,190 (96.5%). Until v7 these were Animal; v7 changed the whole values `chicken` and `turkey`, v8 the rest (`fish`, `chicken liver`). With a host, the host decides.
> 4. **A specimen with no host is Human.** `blood` or `stool` with no host named is a person's in 14,982 NCBI records of 15,290 (98.0%). This default is older than v7; it was measured afterwards and kept.
>
> [documentation/guides/definitions.md](documentation/guides/definitions.md) has the figures behind each and what is still `onehealth`'s own judgement. If your study defines a source differently, one line in a `--decisions` file changes it for every row with that text.

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

**Specimens follow the host.** Body sites and specimens (`blood`, `urine`, `faeces`, `wound`, `tissue`) occur in any host. A recognised animal host or animal word makes them Animal (`pig faeces`, `urine` + `Canis lupus familiaris`). A human host makes them Human. Both at once is a `host_source_conflict`. The word `human` in the specimen's own text makes them Human too (`human feces`), and with no host the name of any animal makes them Animal (`stool (Odocoileus hemionus)`, `pigeon stool`). A host the lexicon does not recognise, such as a plant, goes to review as `host_unrecognised`. With no host at all they default to Human. The reason then ends in `+no_host`, so you can count those calls. The default does not apply to a record its submitter filed as environmental, food or other (BioSample package `Pathogen.env`, or `epi_type` `environmental/other`): there `feces` or `wound` is often an animal's, so it goes to review as `no_host_env_package`. `resolve` writes the package to the `package` column.

**Any animal is a host.** The lexicon's word list names the hosts seen most often (`chicken`, `Bos taurus`, `dog`). A host it does not name is recognised by the names NCBI Taxonomy gives to animals: a scientific name of a chordate at any rank (`Larus dominicanus`, `Macaca`, `Cervidae`), the species name of a mollusc, a crustacean or another animal with a common name (`Magallana gigas`, `Musca domestica`), or the common name of any animal (`gull`, `sika deer`). The value must be a name or end with one (`European rook`, `white stork nestling`), never just contain one, and plants and bacteria are not placed. The names ship with the package (`src/onehealth/data/animal_names.tsv.gz`), so `classify` still works offline.

**Meat is food, the animal needs a host.** How meat is sold or processed means Food, whatever the species: `ground turkey`, `chicken breast`, `comminuted chicken`, `raw intact chicken`, `smoked salmon`. With no host, a meat species in the text is Food too (`fish`, `poultry`, `chicken liver`), as above. It is the animal with a host, with a gut site (`chicken cecal content`), or alive (`live chicken`, `beef cattle`); `duck`, `pig`, `cow` and `cattle` always mean the animal. A carcass is Food with no host and the animal with one. A swab of a surface or of food premises is Environmental even when its text names the food (`environmental swab sponge non food contact surface`, `poultry processing plant`), and a `food handler` is Human.

Bare `clinical` / `clinical isolate` is a patient sample, as in the NCBI `Pathogen.cl` package ("clinical or host-associated"). It is not an environmental swab, and it follows the host like any specimen.

NCBI `source_type` (`human`, `animal`, `food`, `environmental`) is the class of a row that has one; `classify_reason` is `source_type` where the free text could not settle it and `source_type_over_text` where the text said something else. `resolve` writes it in a column of its own. The `tissue` and the NCBI material type (`whole organism`, `pure culture`) are only text, and the material type is never treated as source evidence.

KlebNET `Animal - food producing` is Animal, not Food. `animal feed` is not retail Food and has no L1 here, so it goes to review (`animal-feed`). A taxon-only `isolation_source` (`Escherichia coli`, `Bos taurus`) or an unknown marker (`unknown`, `missing: lab stock`) is ignored. Taxon-only is a guess from word shape, so a lexicon rule written for source text (one that does not read `host`) always wins over it: add `bedpan washer` to such a rule and it is classified. The host or `source_type` then decides, and with nothing else the row goes to review (`taxon_only_source` / `unknown_source`). Bare `sink` or `faucet` is not a hospital fixture; `hospital sink` / `hospital faucet` are. A fixture without its setting (`Rm 14 Drain`, `Sink swab`) goes to review (`fixture-no-setting`), even with host *Homo sapiens*.

Matching is token-based, so `wound drainage` is Human and `hospital drain` is healthcare-environment. The word `drain` does not match `drainage`.

A complete user override TSV always wins. An override row with a missing/empty/invalid `l1` or `paper_rollup`, or two disagreeing rows for the same ID, is queued for review. The lexicon is not applied silently in those cases. Class names are case-insensitive. An override on a `GCF_` accession also applies to its `GCA_` twin, and an unversioned accession covers every version.

Override columns: `id` (or BioSample/assembly/SRA/PDT), `l1`, optional `healthcare_environment`, `paper_rollup`, `note`.

**Counting.** Summary and plot counts are per input row. If the same isolate was listed twice (for example by BioSample and by PDT), the summary's `dataset` section reports `rows_sharing_a_biosample`. Remove duplicates before reporting isolate counts. A shared cluster or plasmid (`mixed=yes`) is a lead, not evidence of transmission. `Unclassified` never makes a cluster mixed.

## Lexicon

`onehealth` ships one lexicon, `src/onehealth/data/lexicon.yaml`, now version 8, and every classified row records `lexicon_version`. When the lexicon changes its version goes up and the changelog says what changed; earlier versions are in the repository history, not in the package. Pass `--lexicon` to use a YAML of your own. It is recorded under its own `version`; one with no version is `custom-<digest>`, and an edited copy that keeps the shipped version is `8-modified-<digest>`, so a table never claims the shipped lexicon when something else made it. The file header explains `priority`, `exact` vs `phrases`, `host_dependent: true` (specimen rules), and `l1: Unclassified` (send to review). This is a small maintained table, not a wrapper of an external R notebook.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

Besides small hand-built cases, the tests run on real public data: 307 *Klebsiella* isolates from NCBI Pathogen Detection with NCBI's own metadata, SNP calls and AMRFinderPlus results, and 166,158 Pathogen Detection records from 26 organism groups that carry NCBI's `source_type`. Provenance and the rebuild scripts are in `tests/data/ncbi_pd_klebsiella/` and `tests/data/ncbi_pd_source_types/`.

GitHub Actions runs the whole suite on Python 3.10 to 3.13 for every pull request (`.github/workflows/tests.yml`).

## Licence

MIT; see [LICENSE](LICENSE). The NCBI data under `tests/data/` and the names in `src/onehealth/data/animal_names.tsv.gz` come from NCBI, which places no restrictions on the use or distribution of its molecular and taxonomy databases ([policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)).
