# One Health analysis: literature review, assessment and roadmap

**Date:** 2026-10-17
**Question asked:** Can `onehealth` help with the kind of analyses recent One Health genomic studies publish, and what would make it a standardised tool that runs those analyses and presents the results clearly?

**Short answer:** Yes, for the part every study currently does by hand. Recent studies all start by harmonising messy source metadata, then run the same small set of population analyses: lineage diversity and sharing, resistance by source, and genomic linkage between sources. `onehealth` already covered the first step. This round adds the second as `onehealth report`, which produces a self-contained HTML report, publication figures and every underlying table. Genome processing (assembly, typing, SNP calling) stays with established tools that `onehealth` reads from.

---

## Contents

1. [What recent One Health genomic studies do](#1-what-recent-one-health-genomic-studies-do)
2. [The common analysis set](#2-the-common-analysis-set)
3. [Where onehealth fits](#3-where-onehealth-fits)
4. [What was built: onehealth report](#4-what-was-built-onehealth-report)
5. [Roadmap](#5-roadmap)
6. [Limitations to keep in mind](#6-limitations-to-keep-in-mind)
7. [References](#7-references)

---

## 1. What recent One Health genomic studies do

I read the methods of the larger recent studies in detail and skimmed smaller ones. The table lists what each did. Numbers are the authors'.

| Study | Organism, setting | Sources | Analyses and thresholds |
|---|---|---|---|
| Thorpe et al. 2022, Pavia [L1] | 3,482 *Klebsiella* genomes, 15 species, 17 months, one city | Hospital, community, livestock, pets, water, soil, plants, invertebrates | PopPUNK strain clusters; Kleborate resistance and virulence scores by source; permutation tests (10,000) with Benjamini–Hochberg for species/lineage–source association; lineage accumulation curves; transmission events at 0, 1, 2, 5, 10 and 20 SNPs, counted once per source pair per cluster and normalised by isolates in the pair. |
| Hetland et al. 2025, Norway [L2] | 3,255 *K. pneumoniae* species complex genomes, 2001–2020 | Human infection and carriage, animals, marine | Simpson's diversity of sublineages per niche (0.99 / 0.98 / 0.96); 12.5% of sublineages in ≥2 niches holding 48.5% of genomes; ≤22-SNP threshold chosen with cutpointR; 130 within-niche vs 9 between-niche sharing events; chi-squared, Kolmogorov–Smirnov and Bonferroni-corrected binomial tests; GWAS with pyseer; dated phylogenies. |
| Winkler et al. 2026, Norway [L3] | 578 KpSC isolates, 1,415 closed plasmids | Human, animal, marine | Plasmid clustering (Pling); 36% of plasmid clusters span more than one niche; animal isolates enriched for virulence plasmids. |
| Muloi et al. 2022, Nairobi [L4] | 1,338 *E. coli* genomes, 99 households | Humans, livestock, wildlife | Sharing pairs at ≤10 cgMLST alleles (cutpointR); clusters counted as n−1 links; observed vs expected sharing given sampling effort (multinomial resampling); resistome sharing by Jaccard similarity. |
| Jacob et al. 2024, Vellore [L5] | 146 *K. pneumoniae* genomes | Clinical, livestock, hospital sewage | Clonal, AMR-gene, virulence and plasmid overlap between settings; found limited overlap. |
| Rodrigues et al. 2022, Europe [L6] | KpSC from chicken meat and salads, 5 countries | Food | Prevalence (60% of chicken meat by qPCR); identical genotypes within the food sector. |
| Rahman et al. 2026, Bangladesh [L7] | ESBL *E. coli*, WHO Tricycle protocol | Hospital wastewater, livestock effluent, aquaculture, river | Niche clustering; 0-SNP hospital–river pairs consistent with effluent spread. |
| Xing et al. 2026, preprint [L8] | Strain-resolved metagenomics plus isolates | Farmers, animals, environment | Strain-sharing networks across compartments. |

Two older anchors are still the reference points:

- **SNP thresholds.** About 21 SNPs separated hospital clusters in the European *K. pneumoniae* survey [L9].
- **Interpretation.** Shared lineages across humans and livestock did not mean recent transmission in the UK *E. coli* study [L10].

**The consistent findings:** sources overlap at the lineage level, but recent transmission is mostly within a source, and cross-source links are rare. The analyses that support this are standard, and almost the same across papers.

---

## 2. The common analysis set

| Analysis | Why studies run it | How it is done | In `onehealth report` |
|---|---|---|---|
| Sampling by source and time | Every comparison depends on sampling balance | Counts by source × year | ✓ |
| Lineage diversity per source | Are some reservoirs more diverse? | Simpson's index with CI [M1, M2] | ✓ |
| Lineage accumulation | Diversity depends on sample size | Rarefaction or accumulation curves [M3] | ✓ |
| Lineage overlap | How much of the population is shared? | % of lineages in ≥2 sources and % of isolates they hold [L2] | ✓ (UpSet-style plot) |
| Lineage–source association | Which lineages are source-specific? | Permutation of source labels, BH [L1] | ✓ |
| Resistance, virulence and plasmid features by source | Where are clinically important genes? | Prevalence and chi-squared tests [L2]; Kleborate scores [L1, T1] | ✓ (permutation chi-squared, Wilson CIs) |
| Pairwise distances within vs between sources | Are same-source isolates closer? | Distance distributions [L2] | ✓ |
| Threshold sensitivity | Thresholds are organism- and pipeline-specific | Sweep 0–20 SNPs [L1] | ✓ (0–100 SNPs plus the chosen value) |
| Strain-sharing events between sources | Is there recent transmission? | Single-linkage clusters at a threshold, one event per source pair per cluster [L1, L2] | ✓ |
| Share of isolates with a close relative in another source | Intuitive effect size | e.g. "4.0% of human infection isolates linked to pigs" [L2] | ✓ |
| Observed vs expected sharing | Adjust for sampling effort | Resampling or permutation null [L4] | ✓ |
| Cross-source clusters | What to follow up | Networks or trees per cluster | ✓ (minimum spanning trees) |
| Time between linked isolates | Recent vs historical links | Years apart [L2] | ✓ (median, per-pair table) |
| Threshold calibration | Choose a defensible threshold | cutpointR on epidemiologically linked pairs [L2, L4] | Roadmap |
| Resistome and plasmid sharing | Genes move on plasmids | Jaccard similarity [L4]; plasmid clusters [L3] | Roadmap |
| Clone-collapsed (dereplicated) tests | Outbreak clusters inflate per-isolate tests | One representative per cluster of closely related genomes [L11] | ✓ (`--dereplicate`) |
| Structure-aware association | Lineages can still share features after dereplication | GWAS (pyseer) [T3] | Out of scope; documented caveat |
| Dated phylogenies | When did lineages cross? | BactDate etc. [L2] | Out of scope |

---

## 3. Where onehealth fits

A One Health genomic study is a chain of steps:

```text
reads ─► assemblies ─► typing (Kleborate / mlst / AMRFinderPlus) ─► distances (snippy or SKA2 + snp-dists)
                                                                         │
NCBI IDs or local metadata ─► onehealth resolve ─► onehealth classify ─► onehealth report ─► HTML, figures, tables
```

**What onehealth does that others don't.** Genome tools (Kleborate [T1], AMRFinderPlus [T2], SKA2 [T4], PopPUNK [M5]) and visualisation platforms (Microreact [T5]) are mature. No standard tool turns heterogeneous isolation-source text into One Health categories reproducibly. Every study above did that step by hand, which is why their source categories differ. `onehealth` makes that step explicit, versioned and auditable, and now carries the result through to the standard population analyses.

**What it should not try to do.** Assembly, typing, SNP calling, phylogenetics and GWAS. These are solved by the tools above. `onehealth` reads their outputs:

- AMRFinderPlus with `--name`
- Kleborate v2 or v3
- mlst
- PopPUNK `clusters.csv`
- snp-dists square or long format

---

## 4. What was built: `onehealth report`

### Usage

```bash
onehealth classify metadata.tsv -o classified.tsv
onehealth report classified.tsv -o report/ \
    --features amrfinder.tsv \      # optional: isolate, gene, class
    --distances snp_dists.tsv \     # optional: snp-dists output
    --lineages kleborate.tsv \      # optional: if ST is not already a column
    --threshold 21                  # SNPs; choose for your organism and pipeline
```

To try it on synthetic data:

```bash
python examples/demo/make_demo.py work/demo
onehealth classify work/demo/metadata.tsv -o work/demo/classified.tsv
onehealth report work/demo/classified.tsv -o work/demo/report \
    --features work/demo/amrfinder.tsv --distances work/demo/snp_dists.tsv \
    --note "Synthetic demonstration data"
```

### Outputs

- `report.html` is a single self-contained file. It has a stat-tile summary, one section per analysis (a one-sentence finding, the figure, a caption, the data table and download links), a data-quality section, an auto-written methods paragraph with numbered references, and hover tooltips on every mark.
- `figures/` holds every figure as SVG and 220-dpi PNG for papers.
- `tables/` holds every number behind every figure as TSV.

Sections appear only when their input is present. A plain classified table gives sampling and lineage sections. Adding `--features` gives the feature heatmap, and adding `--distances` gives the five linkage sections.

### Statistical choices and why

| Choice | Reason |
|---|---|
| One record per isolate (rows sharing a BioSample collapsed) | Avoids pseudo-replication; the earlier summary counted duplicates. |
| Source groups = L1, with healthcare environment split out | Hospital environment vs patients is a key One Health contrast (Kizny Gordon et al. 2017, in REVIEW.md). `--group-by paper_rollup` gives the three-way rollup. |
| Simpson's index with Grundmann CI | Standard for typing data [M1, M2]; reported by [L2]. |
| Exact rarefaction (Hurlbert) | Compares richness at equal sample size [M3]. |
| Permutation tests, p = (b+1)/(m+1), BH correction | Valid for small and sparse tables; never-zero p-values [M4]; same approach as [L1]. |
| Wilson 95% intervals for prevalence | Good coverage at small n and near 0 or 100% [M6]. |
| Events counted once per source pair per cluster, same-host pairs excluded, normalised per 100 isolates of the pair | Exactly the definition in [L1], reused by [L2]. |
| Default threshold 21 SNPs, always with a sweep | 21–22 SNPs is the published *K. pneumoniae* value [L9, L2]. Other organisms differ (≤10 cgMLST alleles for *E. coli* in [L4]), so the sweep is always shown. |
| Permutation null for cross-source sharing | Holds sampling effort fixed, in the spirit of [L4]. |
| Fixed random seed (default 1), recorded in the report | Results reproduce exactly. |

### Design choices

- **Colour follows the source everywhere.** Human is always blue and Animal always orange, from a categorical palette validated for colour-vision deficiency in light and dark mode.
- **Every chart has a data table.** This keeps the report accessible and gives a lighter colour that fails contrast a text fallback.
- **Dark mode is designed, not inverted.** The inline figures are re-themed through CSS variables, while the downloadable SVG and PNG files stay in the light publication style.
- **Fonts are fixed.** Figures use DejaVu Sans, which ships with matplotlib, so they render identically on every machine.
- **Reported language is deliberately cautious.** The report says "close relatives" and "sharing events", never "transmission from A to B", because genomic proximity does not give direction [L10].

---

## 5. Roadmap

Ordered by value to a study, with a rough size. Items 1–4 make the existing analyses more defensible. Items 5–8 widen what the tool covers.

| # | Item | Why | Size |
|---|---|---|---|
| 0 | **Done: classification benchmark.** `onehealth validate` compares classes with curated labels. On 4,351 BioSamples curated for [L11], coverage is 84.5% and accuracy 99.2%. Environmental recall is the weak point (13.5%): hospital-environment strings the lexicon does not know ("Shower 3", "bedpan washer"), environmental swabs recorded with a human host, and records with no source at all. | Measure before changing the lexicon, so every later change shows its effect. | S |
| 1 | **Threshold calibration.** Given known linked pairs (same patient, household or ward), choose the SNP threshold that best separates linked from unlinked pairs, as cutpointR did in [L2, L4]. | Replaces a borrowed default with a data-driven one. | S |
| 2 | ~~**Clone-collapsed association tests.**~~ **Done:** `--dereplicate` keeps one isolate per genomic cluster per source, as in [L11]. On the demonstration data, intrinsic blaSHV-1 went from q = 0.012 to q = 0.21, while the simulated source-associated genes stayed significant. | Isolates from one clone are not independent; per-isolate tests overstate significance. pyseer [T3] stays the strict option. | S |
| 3 | ~~**Time- and place-stratified permutation null.**~~ **Done:** `--stratify-by year,country` restricts every permutation test to isolates from the same year and country. On 307 NCBI Pathogen Detection isolates, expected human–environmental events fell from 8.45 to 5.69 (3 observed). | The unstratified null ignores sampling time and place, so separated sampling also looks like a barrier. | S |
| 4 | **Run file.** `onehealth run study.yaml` captures inputs, threshold, grouping, seed and lexicon version in one file stored with the report. | A standardised, reviewable record of each analysis. | S |
| 5 | **Resistome and plasmid sharing.** Jaccard similarity of AMR gene sets [L4], and sharing of plasmid clusters from MOB-suite or Pling [L3] between sources. | Genes often move on plasmids even when strains do not [L1, L3]. | M |
| 6 | **Microreact export.** Write metadata plus colours and an optional tree in Microreact format [T5]. | Interactive tree, map and timeline exploration without new code. | S |
| 7 | **Workflow wrapper.** A Snakemake or Nextflow recipe from assemblies to report: Kleborate or mlst, AMRFinderPlus, SKA2 or snippy, snp-dists, then onehealth. | One command from genomes to report; the genome tools remain separate. | M |
| 8 | **Metadata standards.** Map outputs to NCBI `source_type` vocabulary and PHA4GE-style contextual fields [T6], and read the WHO Tricycle categories [L7]. | Results comparable across studies. | M |
| 9 | **Organism presets.** Default thresholds and lineage columns per organism (*K. pneumoniae*, *E. coli*, *Salmonella*), each with its citation. | Fewer wrong defaults. | S |
| 10 | **Source attribution.** Model the probable source of human isolates. | Requires training data and validation; a research project, not a feature. | L |

---

## 6. Limitations to keep in mind

- **Shared lineage is not transmission.** Close genomic relatives are evidence of a connection, not of its direction, route or timing [L10]. The report words its findings accordingly.
- **Sampling drives results.** Unequal sampling across sources, years and places shapes every comparison. The sampling section comes first for that reason.
- **Clonal structure inflates per-isolate tests.** On the demonstration data, intrinsic genes present in about 97% of isolates came out "significant", because absences cluster within clones. Run with `--dereplicate` to test one isolate per genomic cluster; lineage-level structure remains, which only a GWAS models.
- **Thresholds are not portable.** SNP counts depend on the reference, masking and calling pipeline. Always read the sensitivity sweep.
- **Feature absence must be explicit.** Only isolates listed in the features file count as genotyped. AMRFinderPlus lists hits only, so an isolate with no hits must appear with an empty feature row, or it will be treated as not genotyped.
- **MDR is not derived automatically.** Counting drug classes from raw gene lists would include intrinsic genes (for *K. pneumoniae*, blaSHV, fosA and oqxAB), making almost every isolate "MDR". Supply MDR from Kleborate's class counts as a feature instead [T1, M7].
- **The demonstration data is synthetic.** It shows how the report looks and behaves, not real biology.

---

## 7. References

Literature was retrieved from PubMed.

**One Health studies**

- **L1** Thorpe HA et al. (2022) A large-scale genomic snapshot of *Klebsiella* spp. isolates in Northern Italy reveals limited transmission between clinical and non-clinical settings. *Nat Microbiol* 7:2054–2067. <https://doi.org/10.1038/s41564-022-01263-0>
- **L2** Hetland MAK et al. (2025) A genome-wide One Health study of *Klebsiella pneumoniae* in Norway reveals overlapping populations but few recent transmission events across reservoirs. *Genome Med* 17:42. <https://doi.org/10.1186/s13073-025-01466-0>
- **L3** Winkler MA et al. (2026) A One Health study of *Klebsiella pneumoniae* species complex plasmids shows a highly diverse and ecologically adaptable plasmidome. *Microb Genom* 12(2). <https://doi.org/10.1099/mgen.0.001629>
- **L4** Muloi DM et al. (2022) Population genomics of *Escherichia coli* in livestock-keeping households across a rapidly developing urban landscape. *Nat Microbiol* 7:581–589. <https://doi.org/10.1038/s41564-022-01079-y>
- **L5** Jacob JJ et al. (2024) Limited evidence of spillover of antimicrobial-resistant *Klebsiella pneumoniae* from animal/environmental reservoirs to humans in Vellore, India. *J Epidemiol Glob Health* 14:1668–1677. <https://doi.org/10.1007/s44197-024-00323-4>
- **L6** Rodrigues C et al. (2022) High prevalence of *Klebsiella pneumoniae* in European food products: a multicentric study comparing culture and molecular detection methods. *Microbiol Spectr* 10:e0237621. <https://doi.org/10.1128/spectrum.02376-21>
- **L7** Rahman A et al. (2026) Environmental dissemination of multidrug-resistant *Escherichia coli* across One Health interfaces in Mymensingh, Bangladesh. *npj Antimicrob Resist* 4. <https://doi.org/10.1038/s44259-026-00202-x>
- **L8** Xing J et al. (2026, preprint) Longitudinal cross-species transmission of microbiomes and resistomes across farmers, animals and environment. *medRxiv*. <https://doi.org/10.64898/2026.05.06.26352545>
- **L9** David S et al. (2019) Epidemic of carbapenem-resistant *Klebsiella pneumoniae* in Europe is driven by nosocomial spread. *Nat Microbiol* 4:1919–1929. <https://doi.org/10.1038/s41564-019-0492-8>
- **L10** Ludden C et al. (2019) One Health genomic surveillance of *Escherichia coli* demonstrates distinct lineages and mobile genetic elements in isolates from humans versus livestock. *mBio* 10:e02693-18. <https://doi.org/10.1128/mBio.02693-18>
- **L11** Vezina B, Morampalli BR, Nguyen HA, Gomez-Simmonds A, Peleg AY, Macesic N (2025) The rise and global spread of IMP carbapenemases (1996–2023): a genomic epidemiology study. *Nat Commun* 17:183. <https://doi.org/10.1038/s41467-025-66874-7>

**Methods**

- **M1** Hunter PR, Gaston MA (1988) Numerical index of the discriminatory ability of typing systems: an application of Simpson's index of diversity. *J Clin Microbiol* 26:2465–2466. <https://doi.org/10.1128/jcm.26.11.2465-2466.1988>
- **M2** Grundmann H, Hori S, Tanner G (2001) Determining confidence intervals when measuring genetic diversity and the discriminatory abilities of typing methods for microorganisms. *J Clin Microbiol* 39:4190–4192. <https://doi.org/10.1128/JCM.39.11.4190-4192.2001>
- **M3** Hurlbert SH (1971) The nonconcept of species diversity: a critique and alternative parameters. *Ecology* 52:577–586. <https://doi.org/10.2307/1934145>
- **M4** Phipson B, Smyth GK (2010) Permutation p-values should never be zero. *Stat Appl Genet Mol Biol* 9:Article39. <https://doi.org/10.2202/1544-6115.1585>
- **M5** Lees JA et al. (2019) Fast and flexible bacterial genomic epidemiology with PopPUNK. *Genome Res* 29:304–316. <https://doi.org/10.1101/gr.241455.118>
- **M6** Wilson EB (1927) Probable inference, the law of succession, and statistical inference. *J Am Stat Assoc* 22:209–212. <https://doi.org/10.1080/01621459.1927.10502953>
- **M7** Magiorakos AP et al. (2012) Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria: an international expert proposal for interim standard definitions for acquired resistance. *Clin Microbiol Infect* 18:268–281. <https://doi.org/10.1111/j.1469-0691.2011.03570.x>

**Tools and standards**

- **T1** Lam MMC et al. (2021) A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex (Kleborate). *Nat Commun* 12:4188. <https://doi.org/10.1038/s41467-021-24448-3>
- **T2** Feldgarden M et al. (2021) AMRFinderPlus and the Reference Gene Catalog. *Sci Rep* 11:12728. <https://doi.org/10.1038/s41598-021-91456-0>
- **T3** Lees JA et al. (2018) pyseer: a comprehensive tool for microbial pangenome-wide association studies. *Bioinformatics* 34:4310–4312. <https://doi.org/10.1093/bioinformatics/bty539>
- **T4** Derelle R et al. (2024) Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis (SKA2). *Genome Res* 34:1661–1673. <https://doi.org/10.1101/gr.279449.124>
- **T5** Argimón S et al. (2016) Microreact: visualizing and sharing data for genomic epidemiology and phylogeography. *Microb Genom* 2:e000093. <https://doi.org/10.1099/mgen.0.000093>
- **T6** Griffiths EJ et al. (2022) Future-proofing and maximizing the utility of metadata: the PHA4GE SARS-CoV-2 contextual data specification package. *GigaScience* 11:giac003. <https://doi.org/10.1093/gigascience/giac003>
