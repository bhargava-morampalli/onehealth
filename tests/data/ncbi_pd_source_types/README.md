# NCBI `source_type` benchmark

Real public data used by `tests/test_source_type_benchmark.py` to measure the classifier against a label that `onehealth` did not make. Nothing here was simulated or hand-written.

| | |
|---|---|
| Source | NCBI Pathogen Detection metadata tables of 18 organism groups, read on 14 and 15 November 2026. The releases are listed at the top of `records.tsv` and in `build_fixture.py`. |
| Label | `source_type`, NCBI's own One Health category (`human`, `animal`, `food`, `environmental`), filled in by the submitter |
| Rebuild | `python tests/data/ncbi_pd_source_types/build_fixture.py` (standard library only, under a minute) |
| Size | 139,554 records in 11,833 rows |
| Licence | NCBI places no restrictions on its molecular databases ([policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)). |

## How it is used

In normal use `classify` takes a record's `source_type` as its class, so on these records it would simply repeat the label. The tests hide it: each record is classified from `isolation_source`, `host` and `host_disease` alone, and the L1 class is compared with the label. That measures how well the lexicon reads source text, which is all `classify` has for the two records in three that carry no `source_type`.

The same records also back the figures in `documentation/guides/definitions.md` on how NCBI submitters record manure, farms and a bare `chicken` (`tests/test_ncbi_conventions.py`).

## Columns

| Column | Meaning |
|---|---|
| `split` | `development`, `heldout`, `validation` or `final` (below) |
| `organism_group` | The Pathogen Detection organism group |
| `isolation_source`, `host`, `host_disease` | As NCBI publishes them; `NULL` is written as an empty cell |
| `ncbi_source_type` | The label. Named so that `classify` does not read it. |
| `records` | How many records have exactly these values |

A record with neither an `isolation_source` nor a `host` is left out: no text classifier can place it.

## The four splits

The tables are large (Salmonella is 900 MB), so the build script reads evenly spaced parts of each with HTTP range requests. The byte ranges of the splits do not overlap, so no record is in two splits.

| Split | Records | What it is | What it was used for |
|---|---|---|---|
| `development` | 33,372 | Eight parts of each of six groups: Salmonella, Escherichia coli and Shigella, Listeria, Campylobacter, Klebsiella, Acinetobacter | Lexicon v5 and v7 were written while looking at these. |
| `heldout` | 39,722 | Other parts of the same six tables, and Staphylococcus aureus, Vibrio parahaemolyticus, Bacillus cereus group, Cronobacter | Unseen by v5. v7 was checked against it and three rules were narrowed. |
| `validation` | 32,208 | Other parts of the same six tables, and Clostridium perfringens, Vibrio cholerae, Yersinia enterocolitica | Meant to be read once after v7. It showed that a farm with an animal host is the animal in NCBI's records, and that rule was corrected. |
| `final` | 34,252 | The last parts of the same six tables, and Vibrio vulnificus, Listeria innocua, Staphylococcus pseudintermedius, Clostridium botulinum, Aeromonas | Read once, after v7 was finished. v7 was not changed afterwards. |

The records of one split often use the same wording as another's, because the same laboratories submitted both. That is how real data repeats, and it is why each split adds organism groups the earlier ones do not have.

**When the lexicon changes again**, the `final` split is no longer unseen for the new version once it has been looked at. Add a split to `build_fixture.py` (other organism groups, or a newer release of these) and report the new version on that.

## Baseline

Coverage is the share of records classified and not sent to review. Accuracy is the share of classified records whose L1 class matches the label. v7 is the lexicon that ships, and the tests hold it to these figures.

| v7 | Development | Held out | Validation | Final |
|---|---|---|---|---|
| Coverage | 95.8% | 90.3% | 96.7% | 95.5% |
| Accuracy | 97.5% | 96.1% | 96.1% | 96.5% |
| Human precision / recall | 99.4% / 98.3% | 98.1% / 98.7% | 99.0% / 99.1% | 99.2% / 98.5% |
| Animal precision / recall | 85.9% / 92.5% | 84.2% / 90.8% | 86.5% / 93.6% | 91.7% / 97.0% |
| Food precision / recall | 99.1% / 89.0% | 98.3% / 73.9% | 99.2% / 86.5% | 98.6% / 83.3% |
| Environmental precision / recall | 99.5% / 92.9% | 98.2% / 87.3% | 99.1% / 86.4% | 97.7% / 90.2% |

Earlier versions on the held-out split, as each scored when it was current: v4 76.9% coverage and 80.9% accuracy, v5 88.7% and 93.9%, v6 89.1% and 93.9%. On the final split v6 scores 94.6% and 95.3%.

## What the label is and is not

`source_type` is the submitter's choice, not a curated truth.

- **Submitters do not all agree.** Manure is `environmental` in 646 records and `animal` in 101; a bare `chicken` with no host is `food` in 810 and `animal` in 56. `onehealth` goes with the majority (`documentation/guides/definitions.md`), so the minority count as wrong calls here.
- **Slaughter samples.** The US Food Safety and Inspection Service labels caecal contents taken at slaughter `food` (`animal-cattle-steer`, host bovine); read as free text they are Animal.
- **Some labels are mistakes**: livers of a bald eagle and two bats labelled `food`.
