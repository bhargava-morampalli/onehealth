# NCBI `source_type` benchmark

Real public data used by `tests/test_source_type_benchmark.py` to measure the classifier against a label that `onehealth` did not make. Nothing here was simulated or hand-written.

| | |
|---|---|
| Source | NCBI Pathogen Detection metadata tables of ten organism groups, read on 14 November 2026. The releases are listed at the top of `records.tsv` and in `build_fixture.py`. |
| Label | `source_type`, NCBI's own One Health category (`human`, `animal`, `food`, `environmental`), filled in by the submitter |
| Rebuild | `python tests/data/ncbi_pd_source_types/build_fixture.py` (standard library only, under a minute) |
| Size | 73,094 records in 6,144 rows |
| Licence | NCBI places no restrictions on its molecular databases ([policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)). |

## How it is used

`classify` reads `source_type` when a record has one, so on these records it would simply repeat the label. The test hides it: each record is classified from `isolation_source`, `host` and `host_disease` alone, and the L1 class is compared with the label. That measures how well the lexicon reads source text, which is all `classify` has for the two records in three that carry no `source_type`.

## Columns

| Column | Meaning |
|---|---|
| `split` | `development` or `heldout` (below) |
| `organism_group` | The Pathogen Detection organism group |
| `isolation_source`, `host`, `host_disease` | As NCBI publishes them; `NULL` is written as an empty cell |
| `ncbi_source_type` | The label. Named so that `classify` does not read it. |
| `records` | How many records have exactly these values |

A record with neither an `isolation_source` nor a `host` is left out: no text classifier can place it.

## The two splits

The tables are large (Salmonella is 900 MB), so the build script reads evenly spaced 3 MB parts of each with HTTP range requests.

| Split | Records | What it is |
|---|---|---|
| `development` | 33,372 | Eight parts of each of six groups: Salmonella, Escherichia coli and Shigella, Listeria, Campylobacter, Klebsiella, Acinetobacter. Lexicon v5 was written while looking at these records. |
| `heldout` | 39,722 | Eight other parts of the same six tables, half-way between the first eight, and four groups not in the development split: Staphylococcus aureus, Vibrio parahaemolyticus, Bacillus cereus group, Cronobacter. Nobody looked at these before v5 was finished, and v5 was not changed afterwards. |

The held-out records are different records, but many use the same wording as the development records, because the same laboratories submitted both. That is how real data repeats, and it is why the four unseen organism groups are there too.

## Baseline

Coverage is the share of records classified and not sent to review. Accuracy is the share of classified records whose L1 class matches the label. v6 is the lexicon that ships, and the tests hold it to its columns. The v4 and v5 columns are what the lexicon scored when each was current.

| Held out | v4 | v5 | v6 |
|---|---|---|---|
| Coverage | 76.9% | 88.7% | 89.1% |
| Accuracy | 80.9% | 93.9% | 93.9% |
| Human precision / recall | 98.0% / 98.6% | 98.0% / 98.6% | 98.0% / 98.6% |
| Animal precision / recall | 45.4% / 89.0% | 73.7% / 88.3% | 74.3% / 91.7% |
| Food precision / recall | 92.6% / 30.0% | 97.2% / 70.6% | 97.3% / 70.6% |
| Environmental precision / recall | 99.6% / 43.8% | 99.0% / 77.6% | 99.0% / 77.6% |

| Development | v4 | v5 | v6 |
|---|---|---|---|
| Coverage | 83.9% | 93.9% | 94.3% |
| Accuracy | 78.5% | 94.9% | 95.0% |
| Human precision / recall | 99.4% / 98.2% | 99.4% / 98.2% | 99.4% / 98.2% |
| Animal precision / recall | 39.2% / 85.8% | 73.6% / 87.4% | 74.6% / 92.3% |
| Food precision / recall | 85.9% / 30.8% | 98.0% / 87.0% | 98.0% / 87.0% |
| Environmental precision / recall | 98.9% / 59.1% | 98.8% / 79.3% | 99.4% / 79.3% |

v5 was written from the development split and never saw the held-out split. v6's rules are general (any animal NCBI Taxonomy names) and were not written from either split, but by then the held-out figures for v5 had been read once, so the held-out split is no longer unseen in the strict sense.

## What the label is and is not

`source_type` is the submitter's choice, not a curated truth, and some disagreements are about definitions:

- **Farm samples.** NCBI submitters label manure, poultry-house swabs and barn samples `environmental`; `onehealth` calls a farm setting Animal.
- **Slaughter samples.** The US Food Safety and Inspection Service labels caecal contents taken at slaughter `food` (`animal-cattle-steer`, host bovine); `onehealth` calls them Animal.
- **A species on its own.** `chicken` or `turkey` with no host is labelled `food` far more often than `animal` in these records. `onehealth` still calls it Animal, because the word names the animal and the lexicon cannot know the submitter was a food laboratory.

These three account for most of what is left of the Animal precision gap.
