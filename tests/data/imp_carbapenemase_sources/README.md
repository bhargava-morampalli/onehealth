# IMP carbapenemase source-classification benchmark

Real public data used by `tests/test_validation.py` to measure the classifier against manually curated source labels. Nothing here was simulated or hand-written.

| | |
|---|---|
| Truth labels | Vezina B, Morampalli BR, Nguyen HA, Gomez-Simmonds A, Peleg AY, Macesic N (2025) The rise and global spread of IMP carbapenemases (1996–2023): a genomic epidemiology study. *Nat Commun* 17:183. <https://doi.org/10.1038/s41467-025-66874-7>. Supplementary Data, Table S1, column `isolation_type_nenad` |
| Classifier input | NCBI BioSample records fetched from NCBI Datasets on 17 October 2026 |
| Rebuild | `python tests/data/imp_carbapenemase_sources/build_fixture.py` (about 1 minute; run with the project environment) |
| Size | 4,351 BioSamples, one row each |
| Licence | The article and its supplementary data are CC BY 4.0. NCBI places no restrictions on its molecular databases ([policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)). |

## Columns

| Column | Meaning |
|---|---|
| `biosample` | BioSample accession |
| `host`, `isolation_source`, `host_disease`, `sample_type` | What NCBI publishes today, mapped with onehealth's own `attributes_from_biosample`, so the classifier sees exactly what `onehealth resolve` would give it |
| `package` | NCBI BioSample package (e.g. `Pathogen.cl.1.0`); not used by the classifier |
| `curated_label` | The paper's curated category, as published |
| `truth` | The curated category mapped to the paper's three source classes (below) |
| `truth_healthcare` | `yes` when the curated category was a clinical (hospital) environment |
| `paper_isolation_source` | The isolation source as listed in Table S1, for reference |

## Mapping to the paper's source classes

| Curated category | `truth` |
|---|---|
| human, clinical | Human |
| clinical environment, clinical/clinical environment, environmental/other | Environmental |
| animal | Animal |
| drain ("drain drainage") | Unclassified; these 9 are the paper's unclassified samples |

The 69 genomes in Table S1 without a public BioSample are excluded. One BioSample had two curated labels (human and clinical); both map to Human.

## Baseline (lexicon v2)

Measured with `onehealth classify` then `onehealth validate`:

| Metric | Value |
|---|---|
| Coverage (classified, not sent to review) | 84.5% |
| Accuracy (of classified isolates) | 99.2% |
| Human recall / precision | 92.0% / 99.1% |
| Environmental recall / precision | 13.5% / 100% |
| Animal recall / precision | 74.5% / 86.4% |

Most of the gap is missing information, not wrong calls. 555 of the 674 review cases have no isolation source and no host in their current BioSample record; the curators had other context, such as NCBI's isolation type. EBI-submitted records (SAME) rarely carry an isolation source. The paper's `isolation_source` column sometimes holds values that the current BioSample does not, so it is kept only for reference.
