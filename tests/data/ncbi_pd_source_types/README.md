# NCBI `source_type` benchmark

Real public data used by `tests/test_source_type_benchmark.py` to measure the classifier against a label that `onehealth` did not make. Nothing here was simulated or hand-written.

| | |
|---|---|
| Source | NCBI Pathogen Detection metadata tables of 26 organism groups, read on 14, 15 and 16 November 2026. The releases are listed at the top of `records.tsv` and in `build_fixture.py`. |
| Label | `source_type`, NCBI's own One Health category (`human`, `animal`, `food`, `environmental`), filled in by the submitter |
| Rebuild | `python tests/data/ncbi_pd_source_types/build_fixture.py` (standard library only, about two minutes) |
| Size | 166,158 records in 13,696 rows |
| Licence | NCBI places no restrictions on its molecular databases ([policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)). |

## How it is used

In normal use `classify` takes a record's `source_type` as its class, so on these records it would simply repeat the label. The tests hide it: each record is classified from `isolation_source`, `host` and `host_disease` alone, and the L1 class is compared with the label. That measures how well the lexicon reads source text, which is all `classify` has for the two records in three that carry no `source_type`.

The same records also back the figures in `documentation/guides/definitions.md` on how NCBI submitters record manure, farms and meat species (`tests/test_ncbi_conventions.py`). Those figures use the first four splits, 139,554 records; the `fresh` split is kept apart to test them.

## Columns

| Column | Meaning |
|---|---|
| `split` | `development`, `heldout`, `validation`, `final` or `fresh` (below) |
| `organism_group` | The Pathogen Detection organism group |
| `isolation_source`, `host`, `host_disease` | As NCBI publishes them; `NULL` is written as an empty cell |
| `ncbi_source_type` | The label. Named so that `classify` does not read it. |
| `records` | How many records have exactly these values |

A record with neither an `isolation_source` nor a `host` is left out: no text classifier can place it.

## The five splits

The tables are large (Salmonella is 900 MB), so the build script reads evenly spaced parts of each with HTTP range requests. The byte ranges of the splits do not overlap, so no record is in two splits.

| Split | Records | What it is | What it was used for |
|---|---|---|---|
| `development` | 33,372 | Eight parts of each of six groups: Salmonella, Escherichia coli and Shigella, Listeria, Campylobacter, Klebsiella, Acinetobacter | Lexicon v5 and v7 were written while looking at these. |
| `heldout` | 39,722 | Other parts of the same six tables, and Staphylococcus aureus, Vibrio parahaemolyticus, Bacillus cereus group, Cronobacter | Unseen by v5. v7 was checked against it and three rules were narrowed. |
| `validation` | 32,208 | Other parts of the same six tables, and Clostridium perfringens, Vibrio cholerae, Yersinia enterocolitica | Meant to be read once after v7. It showed that a farm with an animal host is the animal in NCBI's records, and that rule was corrected. |
| `final` | 34,252 | Other parts of the same six tables, and Vibrio vulnificus, Listeria innocua, Staphylococcus pseudintermedius, Clostridium botulinum, Aeromonas | Read once, after v7 was finished; v7 was not changed afterwards. Lexicon v8 was written with all four splits above in view. |
| `fresh` | 26,604 | The parts still unread of four of the six tables (Listeria and Acinetobacter had none left), and Vibrio alginolyticus, Enterococcus faecalis, Enterococcus faecium, Citrobacter freundii, Aeromonas veronii, Pseudomonas aeruginosa, Enterobacter hormaechei, Serratia | Drawn after v8 was finished, and read once. One value was corrected afterwards (below). |

The records of one split often use the same wording as another's, because the same laboratories submitted both. That is how real data repeats, and it is why each split adds organism groups the earlier ones do not have.

The eight organism groups of the `fresh` split were chosen by name before any of their records were read, and the rebuild reproduced the four earlier splits byte for byte.

**When the lexicon changes again**, the `fresh` split is no longer unseen for the new version once it has been looked at. Add a split to `build_fixture.py` (other organism groups, or a newer release of these) and report the new version on that. Salmonella and Escherichia coli still have unread parts; the other four shared tables do not.

## Baseline

Coverage is the share of records classified and not sent to review. Accuracy is the share of classified records whose L1 class matches the label.

**The figure to quote for v8** is from the `fresh` split as it was first read, with v8 finished and nothing yet changed:

| Fresh split | v7 | v8, as first read |
|---|---|---|
| Coverage | 97.1% | 97.8% |
| Accuracy | 95.4% | 96.7% |
| Human precision / recall | 98.8% / 95.8% | 98.8% / 95.9% |
| Animal precision / recall | 90.7% / 91.2% | 97.4% / 92.0% |
| Food precision / recall | 98.7% / 91.7% | 98.2% / 97.4% |
| Environmental precision / recall | 81.4% / 84.2% | 81.4% / 85.1% |

That reading showed 296 records described only as `healthcare facility`, all human to their submitter. v7 had taken the phrase for a hospital fixture with no labelled record to go by. It was corrected: the value is now read like `hospital`. Nothing else was changed. The lexicon that ships therefore scores higher on this split than the table says (last column below), and that figure is not a fair test.

The lexicon that ships, on every split. The tests hold it to these figures:

| v8 | Development | Held out | Validation | Final | Fresh |
|---|---|---|---|---|---|
| Coverage | 98.0% | 97.1% | 98.8% | 98.3% | 97.8% |
| Accuracy | 98.5% | 97.9% | 98.1% | 98.3% | 97.9% |
| Human precision / recall | 99.4% / 99.4% | 99.4% / 99.4% | 99.1% / 99.2% | 99.4% / 98.8% | 98.8% / 98.6% |
| Animal precision / recall | 93.3% / 95.0% | 92.4% / 95.1% | 94.6% / 97.4% | 97.2% / 97.8% | 97.4% / 92.0% |
| Food precision / recall | 99.0% / 94.8% | 98.4% / 92.9% | 99.2% / 96.0% | 98.5% / 94.4% | 98.2% / 97.4% |
| Environmental precision / recall | 99.5% / 95.6% | 98.3% / 91.9% | 99.1% / 89.4% | 97.7% / 95.1% | 93.9% / 85.1% |

Earlier versions on the held-out split, as each scored when it was current: v4 76.9% coverage and 80.9% accuracy, v5 88.7% and 93.9%, v6 89.1% and 93.9%, v7 90.3% and 96.1%. On the final split, read once for v7: v6 94.6% and 95.3%, v7 95.5% and 96.5%.

## What the label is and is not

`source_type` is the submitter's choice, not a curated truth.

- **Submitters do not all agree.** Manure is `environmental` in 646 records and `animal` in 101; a bare `chicken` with no host is `food` in 810 and `animal` in 56; a bare `hospital` is `human` in all 328. `onehealth` goes with the majority (`documentation/guides/definitions.md`), so the minority count as wrong calls here.
- **Slaughter samples.** The US Food Safety and Inspection Service labels caecal contents taken at slaughter `food` (`animal-cattle-steer`, host bovine); read as free text they are Animal.
- **Some labels are mistakes**: livers of a bald eagle and two bats labelled `food`.
