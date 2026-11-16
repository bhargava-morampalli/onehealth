# Reference values for the statistics

`onehealth report` computes its statistics with its own code. This directory holds the values that R and the vegan package give for the same inputs, so that `tests/test_stats_reference.py` can compare the two every time the tests run. The documentation page "The statistics, and how each was checked" (`documentation/guides/statistics.md`) describes the result.

| File | What it is |
|---|---|
| `cases.tsv` | 88 inputs: counts from the real NCBI *Klebsiella* fixture (`tests/data/ncbi_pd_klebsiella`), and edge cases written by hand. |
| `make_cases.py` | Wrote `cases.tsv`. |
| `reference.R` | Reads `cases.tsv` and computes the reference values with base R and vegan. Uses nothing from `onehealth`. |
| `reference.tsv` | The 567 values `reference.R` wrote, with R 4.3.3 and vegan 2.6-4. |

## What R computes

| Quantity | R function |
|---|---|
| Simpson's index of diversity | `vegan::simpson.unb` |
| Sampling variance of the index | the variance of the index over 2,000,000 multinomial samples (`rmultinom`), with its standard error |
| Rarefaction | `vegan::rarefy` |
| Wilson interval | `prop.test(correct = FALSE)` |
| Benjamini–Hochberg q-values | `p.adjust(method = "BH")` |
| Chi-squared statistic | `chisq.test(correct = FALSE)` |
| Permutation p-value of the chi-squared test | `chisq.test(simulate.p.value = TRUE, B = 1e6)` |
| Lineage enrichment p-value | exact: a sum of hypergeometric probabilities (`dhyper`) |

## Rebuilding

```bash
Rscript tests/data/stats_reference/reference.R
```

It takes about 20 seconds and needs R with the vegan package (`apt install r-base-core r-cran-vegan` on Debian or Ubuntu). The simulations use a fixed seed, so the same version of R writes the same file. The tests do not need R: they read `reference.tsv`.

`cases.tsv` is fixed. `make_cases.py` takes its counts from the *Klebsiella* fixture as the classifier places those records, so running it again after the lexicon changes can change the counts. That would not make the comparison wrong, since R and `onehealth` always get the same inputs, but `reference.R` must then be run again, and the numbers quoted in the documentation updated (a test checks them).

## Two more checks that need no file

`tests/test_stats_reference.py` also compares the tables `report` makes from the real *Klebsiella* data with scipy, and compares every permutation p-value with the exact one on a dataset of ten isolates, where all 4,200 ways of relabelling the isolates can be listed.
