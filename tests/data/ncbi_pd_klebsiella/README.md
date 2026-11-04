# NCBI Pathogen Detection *Klebsiella* test fixture

Real public data used by `tests/test_real_data.py`. Nothing here was simulated or hand-written.

| | |
|---|---|
| Source | NCBI Pathogen Detection, *Klebsiella* release `PDG000000012.2554` (17 October 2026) |
| URL | <https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/PDG000000012.2554/> |
| Rebuild | `python tests/data/ncbi_pd_klebsiella/build_fixture.py` (standard library only, about 2 minutes) |
| Size | 307 isolates in 22 NCBI SNP clusters, 2,676 within-cluster pairs |

## Files

- **metadata.tsv**: BioSample metadata exactly as NCBI publishes it (`host`, `isolation_source`, `source_type`, `epi_type`, `collection_date`, `geo_loc_name`, and others), plus each isolate's NCBI SNP cluster (`PDS_acc`). `NULL` is written as an empty cell.
- **features.tsv**: NCBI's AMRFinderPlus results from the `AMR_genotypes` and `virulence_genotypes` columns, one row per isolate and gene. Allele status suffixes such as `=COMPLETE` are removed. NCBI wraps these list cells in quotes; the quotes are removed too (before 4 November 2026 the first and last gene of each list kept a stray quote, so 23 genes were counted under two names). An isolate genotyped by NCBI with no hits has one row with an empty feature.
- **snp_dists.tsv**: pairwise SNP differences within each cluster. They are counted from NCBI's per-cluster SNP calls (`SNP_trees/<PDS>.tar.gz`, `variation.vcf.gz`), using sites where both isolates have a PASS call.

## Selection

The build script takes SNP clusters of 4–30 isolates spanning at least two collection years, in order of cluster accession. It keeps the first 16 that contain both NCBI `epi_type` "clinical" and "environmental/other" isolates, and the first 6 that are clinical only.

## Known differences from NCBI

- **Distances.** NCBI defines clusters with its own internal SNP distance. Counting PASS sites in the published VCFs reproduces 21 of the 22 clusters exactly at NCBI's 50-SNP linkage. In `PDS000036208.8`, one isolate is 58 SNPs from its nearest cluster member by this count. The fixture records the counts as computed; `test_real_data.py` checks this agreement.
- **Missing pairs.** Pairs in different clusters are not listed. NCBI clusters are built by single linkage at 50 SNPs, so such pairs cannot fall within a smaller threshold.
- **Country names.** Old records use historical names, for example "Czech Republic" (37 isolates) and "Czechia" (11).

## Terms of use

NCBI states that it places no restrictions on the use or distribution of the data in its molecular databases, while noting that original submitters may hold rights. See the [NCBI website and data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/). Cite NCBI Pathogen Detection and the original BioProjects when using these records beyond testing.
