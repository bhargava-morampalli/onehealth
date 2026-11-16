# The statistics, and how each was checked

`report` computes every number with `onehealth`'s own code, so that installing it needs nothing heavier than numpy. That raises a fair question: is the code right?

To answer it, every statistic in a report is compared with an independent implementation each time the test suite runs: **R and the vegan package**, **scipy**, or, for the permutation tests, **the exact answer from listing every possible relabelling**. They agree. This page says what was compared and how closely.

## Summary

| Statistic | In the report | Checked against | Agreement |
|---|---|---|---|
| Simpson's index of diversity (Hunter & Gaston 1988) | Lineages: diversity | `vegan::simpson.unb` | 12 decimal places |
| Its sampling variance, behind the 95% interval (Simpson 1949) | Lineages: diversity | the index of 2,000,000 samples drawn in R | within 0.2% |
| Rarefaction (Hurlbert 1971) | Lineages: accumulation | `vegan::rarefy`; scipy's hypergeometric distribution | 9 decimal places |
| Wilson 95% interval (Wilson 1927) | Features: prevalence | R `prop.test(correct = FALSE)`; scipy `binomtest` | 12 decimal places |
| Benjamini–Hochberg q-values | every test | R `p.adjust(method = "BH")`; scipy `false_discovery_control` | 12 decimal places |
| Chi-squared statistic | Features: test | R `chisq.test(correct = FALSE)`; scipy `chi2_contingency` | 10 significant digits |
| Permutation p-value of that test | Features: test | R `chisq.test(simulate.p.value = TRUE)` with 1,000,000 random tables; every relabelling | within random error |
| Lineage enrichment p-value | Lineages: common lineages | the hypergeometric distribution (R `dhyper`), which is exact; every relabelling | within random error |
| Clusters by single linkage | Genomic linkage | scipy `connected_components` | the same clusters |
| Sharing events and their p-value | Genomic linkage: observed and expected | every relabelling, with the events counted from their definition | the same events; p within random error |
| Shuffling within strata (`--stratify-by`) | every test | every relabelling within year | each one equally often |

"Within random error" is explained [below](#permutation-tests).

## The comparison with R

`tests/data/stats_reference/` holds three files:

- **`cases.tsv`**, 88 inputs. Most are counts from the real NCBI *Klebsiella* data that ships with the tests: isolates per lineage in each source, and isolates with each resistance gene in each source. The rest cover the edges: one isolate, lineages all equally common, a gene in every isolate, none of 40, a p-value that was not tested.
- **`reference.R`**, which reads those inputs and computes 567 values with base R and vegan. It uses nothing from `onehealth`.
- **`reference.tsv`**, the values it wrote (R 4.3.3, vegan 2.6-4).

The tests give the same inputs to `onehealth` and compare. A few of the values, from the real data:

| Statistic | onehealth | R |
|---|---|---|
| Simpson's index, 231 human isolates in 22 lineages | 0.941013 | 0.941013 |
| Expected lineages in 58 of those isolates | 19.2841 | 19.2841 |
| Wilson 95% interval, `oqxB` in 151 of 231 human isolates | 59.0261–71.2071% | 59.0261–71.2071% |
| Chi-squared, `aac(6')-Ib-cr5` by source | 15.0862 | 15.0862 |
| Its permutation p-value | 0.0012 | 0.0012 |
| Enrichment of lineage `PDS000045679.7` in Human (27 isolates, 24.9 expected), p | 0.2247 | 0.2219 (exact) |

The largest difference in any value that does not involve chance is about 0.00000000001 (rarefaction). To compute the reference again:

```bash
Rscript tests/data/stats_reference/reference.R
```

## Permutation tests

A permutation test asks how often shuffled source labels give a result at least as extreme as the real one. `onehealth` tries 9,999 random shuffles by default, so its p-value is an estimate, with a standard error of at most 0.005. Two things are checked.

**That the estimate is of the right quantity.** Where the exact p-value is known, the estimate must be within 4.5 standard errors of it:

- For lineage enrichment with no strata, the number of isolates of a lineage that fall in a source follows a hypergeometric distribution, so the exact p-value is a sum of its probabilities. Checked for 107 lineage and source pairs.
- For the feature test, R draws 1,000,000 random tables with the same totals, which gives the p-value to about three decimal places. Checked for 47 genes.
- For every test, with and without `--stratify-by`, on a dataset of ten isolates: there are 4,200 ways to hand three source labels back to ten isolates (900 within two years), and all of them can be listed. The statistic is recomputed for each in a few plain lines of test code, not with `onehealth`'s functions. That gives the exact p-value for enrichment, for the feature test and for sharing events.

**That the shuffle is fair.** In 90,000 shuffles of those ten isolates, each possible relabelling comes up equally often, within and without strata.

## What rests on a definition, not a package

Two quantities have no package to compare with, because the papers `onehealth` follows define them:

- **A sharing event** (Thorpe et al. 2022; Hetland et al. 2025): one pair of sources connected inside one cluster by at least one pair of isolates within the threshold, not counting two isolates of one host. The test counts events from this sentence, with scipy's clusters, and `onehealth` must give the same counts.
- **The 95% interval of Simpson's index** is D ± 2 standard errors (Grundmann et al. 2001), with the exact variance in place of the large-sample one. The variance is checked; the choice of ± 2 standard errors is a convention, and the interval is approximate.

## One change this check made

The Wilson interval used 1.96 for the 97.5% point of the normal distribution. R and scipy use 1.959964. `onehealth` now does too, so that it agrees with them to the last digit. Interval bounds change in the fifth significant digit: for 51 of 150 isolates, 26.9032–41.8959% became 26.9033–41.8957%. Nothing printed in the HTML report changes.

Nothing else needed correcting.

## References

The methods are cited in full on the [`report`](../commands/report.md#references) page.

- Oksanen J et al. (2022) vegan: Community Ecology Package. R package version 2.6-4. <https://CRAN.R-project.org/package=vegan>
- R Core Team (2024) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna. <https://www.R-project.org/>
- Virtanen P et al. (2020) SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods* 17:261–272. <https://doi.org/10.1038/s41592-019-0686-2>
