# report

Runs the population analyses that One Health genomic studies publish, by source, and writes them as one self-contained HTML report with every figure and every number behind it.

```bash
onehealth report work/classified.tsv -o work/report \
    --features amrfinder.tsv --distances snp_dists.tsv
```

The input is the output of [`classify`](classify.md). Each analysis runs when its data is there:

| You give | You get |
|---|---|
| The classified table | Sampling by source |
| A lineage column, or `--lineages` | Lineage diversity, accumulation, sharing, enrichment and overlap |
| `--features` | Resistance, virulence and plasmid features by source |
| `--distances` | Genomic linkage within and between sources |

`report` reads the output of established tools and does no typing, gene calling or SNP calling itself.

## Options

| Option | What it does |
|---|---|
| `-o`, `--out` | Output directory (required). |
| `--features` | Gene or feature calls, one row per isolate and feature; see [Features](#features). |
| `--distances` | Pairwise SNP distances from snp-dists; see [Distances](#distances). |
| `--lineages` | A typing table from Kleborate, mlst or PopPUNK; see [Lineages](#lineages). |
| `--lineage-column` | The lineage column of the input table, if it is not found by name. |
| `--threshold` | SNPs within which two isolates are close relatives (default 21). Choose one for your organism and pipeline; see [The SNP threshold](#the-snp-threshold). |
| `--group-by` | `source` (the default) or `paper_rollup`; see [Source groups](#source-groups). |
| `--stratify-by` | Shuffle source labels only within the same `year`, `country`, or `year,country` in every permutation test. |
| `--dereplicate` | Run the lineage and feature analyses on one isolate per genomic cluster per source. Needs `--distances`. |
| `--permutations` | Permutations for every test (default 9,999). `0` skips the tests, for a quick look. |
| `--seed` | Random seed, 0 or more (default 1). The same seed gives the same p-values. |
| `--title` | Report title (default "One Health population report"). |
| `--note` | A line under the title, for example a description of the dataset. |

## Output

```text
work/report/
├── report.html      one file: figures, tables, data quality, methods, references
├── figures/         every figure as SVG and PNG
└── tables/          every number behind every figure, as TSV
```

`report.html` needs no other files and no internet connection, so it can be emailed or archived as it is. It follows the reader's light or dark setting, and hovering over a bar or cell shows its numbers.

The header records the `onehealth` version, the date and every parameter (grouping, threshold, permutations, seed, dereplication, strata), so a report can be reproduced.

## Source groups

With `--group-by source`, isolates are grouped as Human, Animal, Environmental, Food, Healthcare environment and Laboratory. Healthcare environment is split out of Environmental by the `healthcare_environment` column, because hospital sinks and drains are a reservoir of their own in most studies.

With `--group-by paper_rollup`, the groups are the three classes of the paper rollup: Human, Animal and Environmental.

Either way:

- Unclassified isolates are counted but left out of every comparison between sources.
- Rows that share a BioSample are one isolate. The first row's values are kept, and every ID of every row is used to match lineages, features and distances.
- Groups keep the same order and colour in every figure.

## Sections

| # | Section | Needs | Tables |
|---|---|---|---|
| 1 | Sampling by source | | `01_sampling.tsv` |
| 2 | Lineage diversity | lineages | `02_diversity.tsv` |
| 3 | Lineage accumulation | lineages | `03_rarefaction.tsv` |
| 4 | Common lineages and where they were found | lineages | `04_lineage_sources.tsv`, `04_lineage_enrichment.tsv` |
| 5 | Lineage overlap between sources | lineages | `05_lineage_overlap.tsv` |
| 6 | Resistance, virulence and plasmid features by source | `--features` | `06_feature_prevalence.tsv`, `06_feature_tests.tsv` |
| 7 | How close are isolates from the same and different sources? | `--distances` | `07_distances.tsv` |
| 8 | Sensitivity to the SNP threshold | `--distances` | `08_threshold_sweep.tsv` |
| 9 | Which sources share close relatives? | `--distances` | `09_linked_isolates.tsv` |
| 10 | Is cross-source sharing more or less than chance? | `--distances`, permutations | `10_observed_expected.tsv` |
| 11 | Clusters that span sources | `--distances` | `11_cross_source_pairs.tsv`, `11_clusters.tsv` |

Each section states its finding in one sentence, then shows the figure, a caption with the method, the table, and links to the TSV files. A section without its data is left out, and the Data quality section says why.

### Sampling

Isolates per collection year, by source. Uneven sampling between sources limits every comparison that follows, so it comes first. Without collection dates the table counts isolates per source.

### Lineages

- **Diversity:** Simpson's index of diversity per source, the chance that two isolates from it belong to different lineages, with approximate 95% confidence intervals: D ± 2 standard errors, from the exact sampling variance of the index (Simpson 1949). The large-sample variance usually quoted (Grundmann et al. 2001) is exactly zero when lineages are equally common, so a source of ten isolates in ten lineages had the interval 1.00–1.00; the two agree for large samples.
- **Accumulation:** rarefaction curves, the expected number of lineages in random subsamples of each source. Compare sources at the same number of isolates, because raw lineage counts favour the larger sample. A curve that is still rising means more sampling would find more lineages.
- **Common lineages:** the 20 most common lineages split by source, and which lineages are enriched in or depleted from a source compared with chance.
- **Overlap:** an UpSet-style chart of how many lineages each combination of sources shares.

Sharing a lineage does not mean recent transmission; the genomic linkage sections test that.

### Features

The share of genotyped isolates from each source that carry each feature, with Wilson 95% intervals, and a permutation chi-squared test of whether prevalence depends on source. The heatmap shows the 40 most common features that some, but not all, isolates carry; the tables list every feature.

Isolates from one clone are not independent, so the test overstates significance when clones dominate. Rerun with `--dereplicate`, or use a GWAS method that models population structure, such as pyseer.

### Genomic linkage

Isolates within `--threshold` SNPs of each other are joined into clusters by single linkage.

- **Distances:** pairwise SNP distances up to 100 SNPs, for pairs from the same source and from different sources.
- **Threshold sweep:** sharing events at 0, 1, 2, 5, 10, 20, 50 and 100 SNPs, and at your threshold. A conclusion that holds across thresholds is more robust.
- **Linked isolates:** for each pair of sources, the share of isolates from one with at least one close relative from the other.
- **Observed and expected:** sharing events between each pair of sources compared with the range expected when source labels are shuffled over the same genomes.
- **Clusters that span sources:** clusters with isolates from two or more sources, each drawn as a minimum spanning tree labelled with SNP distances (up to 12 are drawn; the tables list all), and how many years apart the linked isolates were collected.

A sharing event is one source pair connected inside one cluster, counted once per cluster. Pairs from the same host individual are not counted, when the table has a `host_id`, `host_subject_id`, `patient_id`, `subject_id`, `individual_id` or `animal_id` column. Events are also given per 100 isolates of the two sources.

Shared clusters are candidates for follow-up, not proof of the direction or route of transmission.

## Statistics

- **Permutation tests** shuffle source labels and keep everything else fixed, including the genomes, so sampling effort is the same in every permutation. P-values are two-sided, computed as (b+1)/(m+1) so they are never zero, and corrected for multiple testing with Benjamini–Hochberg. Results are called significant at q < 0.05.
- **`--stratify-by year,country`** shuffles labels only among isolates from the same collection year and country. Without it, sources sampled in different years or places also look less linked than chance; the caption of section 10 says so. If most isolates are in strata with a single source, the tests have little power, and Data quality says so.
- **`--dereplicate`** joins isolates of the same source that are within `--threshold` SNPs and keeps the one with the lowest identifier. Isolates from different sources, and unclassified isolates, are never merged. Lineage and feature analyses use the representatives; sampling and linkage use every isolate.
- **`--permutations 0`** skips every test. The report then gives no p- or q-values, calls nothing significant and leaves out section 10; the sharing events are written to `08_sharing_events.tsv` instead.

The Methods section of the report describes what was run in the words of a methods paragraph, and cites each method.

Every statistic here is checked against R, the vegan package or scipy, or against the exact answer, each time the tests run. [The statistics, and how each was checked](../guides/statistics.md) gives the comparison.

### The SNP threshold

The default of 21 SNPs separated hospital clusters of *Klebsiella pneumoniae* in Europe (David et al. 2019); 22 was optimal in a Norwegian One Health study (Hetland et al. 2025). The right threshold depends on the organism, the SNP-calling pipeline and the time scale you care about, so choose it and say why. The threshold sweep shows how much your conclusions depend on it.

## Input files

Names in the features, distances and lineages files are matched to isolates by any of their IDs, in any case: the query ID, BioSample, assembly, SRA or PDT accession, or the `strain`, `isolate`, `sample_name` or `name` column. An assembly file name such as `GCA_000005845.2_ASM584v2_genomic.fna` matches its assembly, and a sequence file matches the isolate it is named after: `ISO1.fasta`, `assemblies/ISO1.fna.gz` and shovill's `ISO1/contigs.fa` all match `ISO1`. The most specific of these forms that names an isolate decides. A name that matches no isolate, or more than one, is not used, and Data quality reports how many.

### Features

A long table with one row per isolate and feature. AMRFinderPlus output run with `--name` fits as it is.

| Column | Taken from the first of |
|---|---|
| Isolate | `id`, `name`, `sample`, `sample_id`, `isolate`, `strain`, `genome`, `assembly`, `biosample`, `query_id` |
| Feature | `feature`, `gene`, `gene_symbol`, `element_symbol` (spaces or underscores) |
| Class, optional | `class`, `feature_class`, `group`, `drug_class`, `category` |

Every isolate named in the file counts as genotyped. To record a genotyped isolate with no hits, give it one row with an empty feature or `-`. Isolates not in the file are not genotyped and are left out of the prevalence, rather than counted as negative.

### Distances

Pairwise SNP distances from snp-dists, in its default square format or its long `-m` format (name, name, distance), tab- or comma-separated. A pair missing from a long table has no distance: it is not counted and never links. A square matrix must list its rows in the same order as its columns; one with rows missing or reordered stops with an error rather than putting distances on the wrong pairs. It must also be the full, symmetric matrix: one that holds a single triangle, or a negative distance, stops with an error, because the other triangle would be read as every pair being 0 SNPs apart.

The distances must be SNP counts. If they are all fractions between 0 and 1, as Mash or ANI distances are, every pair is within any SNP threshold; the report still runs, and Data quality says that the linkage sections mean nothing for that file.

### Lineages

A lineage column in the input table is found by name: the first of `st`, `mlst_st`, `sequence_type`, `mlst`, `lineage`, `sublineage`, `sl`, `poppunk_cluster`, `cluster_id` or `cluster`. Use `--lineage-column` for another name.

Or pass a typing table with `--lineages`: Kleborate output, a PopPUNK `clusters.csv`, mlst output, or any table with an isolate column and a lineage column, tab- or comma-separated.

- **mlst** output is read with or without its header row. The isolate is the file name in its first column, and the lineage is the ST; `-` (no ST assigned) gives no lineage.
- **No type assigned.** A value that says so is not a lineage: `NF` and `ND` (SRST2), `Novel` (ARIBA), `untypeable`, `not typed` and the like, in any case and with a trailing `*` or `?`. In a sequence-type column (`st`, `mlst_st`, `sequence_type`, `mlst`, or Kleborate's `…__ST`), so is `0`, which Kleborate writes when the alleles match no ST. A cluster column keeps its cluster `0`. These isolates count as untyped, and Data quality says how many there are. Left as lineages, every untyped isolate would belong to one lineage found in every source.
- **Any other table** needs a header row. Its lineage column is the first of `st`, `mlst_st`, `sequence_type`, `lineage`, `sublineage`, `sl`, `poppunk_cluster`, `cluster` or a Kleborate `…__ST` column; its isolate column is the first of `id`, `strain`, `name`, `taxon`, `sample`, `sample_id`, `isolate`, `genome`, `assembly`, `biosample`, or else its first column.

An isolate given two different lineages, under one name or under two of its IDs, stops with an error.

### Other columns

| Column | Used for |
|---|---|
| `collection_date` | Sampling by year, years between linked isolates, and `--stratify-by year`. `2019`, `2019-05` and `2019-05-14` all work, and so do the INSDC forms with a month name, `14-May-2019` and `May-2019`. A numeric date that ends in the year, such as `14/05/2019`, counts as the middle of that year, because its day and month cannot be told apart. Data quality says how many dates could not be read at all. |
| `geo_loc_name`, `country` or `geo_loc` | `--stratify-by country`. The country is the part before the colon (`Norway: Oslo`), and historical names are mapped to current ones, so `Czech Republic` and `Czechia` are one country. |
| `host_id`, `patient_id`, ... | Leaving out pairs from the same host individual (above). |

## References

The report cites these where they are used:

- Hunter PR, Gaston MA (1988) Numerical index of the discriminatory ability of typing systems: an application of Simpson's index of diversity. *J Clin Microbiol* 26:2465–2466. <https://doi.org/10.1128/jcm.26.11.2465-2466.1988>
- Grundmann H, Hori S, Tanner G (2001) Determining confidence intervals when measuring genetic diversity and the discriminatory abilities of typing methods for microorganisms. *J Clin Microbiol* 39:4190–4192. <https://doi.org/10.1128/JCM.39.11.4190-4192.2001>
- Simpson EH (1949) Measurement of diversity. *Nature* 163:688. <https://doi.org/10.1038/163688a0>
- Hurlbert SH (1971) The nonconcept of species diversity: a critique and alternative parameters. *Ecology* 52:577–586. <https://doi.org/10.2307/1934145>
- Wilson EB (1927) Probable inference, the law of succession, and statistical inference. *J Am Stat Assoc* 22:209–212. <https://doi.org/10.1080/01621459.1927.10502953>
- Phipson B, Smyth GK (2010) Permutation p-values should never be zero. *Stat Appl Genet Mol Biol* 9:Article39. <https://doi.org/10.2202/1544-6115.1585>
- David S et al. (2019) Epidemic of carbapenem-resistant *Klebsiella pneumoniae* in Europe is driven by nosocomial spread. *Nat Microbiol* 4:1919–1929. <https://doi.org/10.1038/s41564-019-0492-8>
- Muloi DM et al. (2022) Population genomics of *Escherichia coli* in livestock-keeping households across a rapidly developing urban landscape. *Nat Microbiol* 7:581–589. <https://doi.org/10.1038/s41564-022-01079-y>
- Thorpe HA et al. (2022) A large-scale genomic snapshot of *Klebsiella* spp. isolates in Northern Italy reveals limited transmission between clinical and non-clinical settings. *Nat Microbiol* 7:2054–2067. <https://doi.org/10.1038/s41564-022-01263-0>
- Hetland MAK et al. (2025) A genome-wide One Health study of *Klebsiella pneumoniae* in Norway reveals overlapping populations but few recent transmission events across reservoirs. *Genome Med* 17:42. <https://doi.org/10.1186/s13073-025-01466-0>
- Vezina B, Morampalli BR, Nguyen HA, Gomez-Simmonds A, Peleg AY, Macesic N (2025) The rise and global spread of IMP carbapenemases (1996–2023): a genomic epidemiology study. *Nat Commun* 17:183. <https://doi.org/10.1038/s41467-025-66874-7>
