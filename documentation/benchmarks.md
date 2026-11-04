# Benchmarks

How well does `classify` do on real data? Two public datasets ship with `onehealth` in `tests/data/`, and the test suite fails if a change makes the results below worse.

!!! note
    The v3 additions to the lexicon were found in the review queues of these two datasets. The v3 numbers show the additions work as intended; they are not a measure of accuracy on data the lexicon has not seen. Run [`validate`](commands/validate.md) on labels from your own study before relying on the classes.

## IMP carbapenemase isolates

4,351 public isolates from a genomic epidemiology study of IMP carbapenemases, whose authors curated the source of every isolate by hand (Vezina et al. 2025, *Nat Commun* 17:183, Supplementary Table S1). `onehealth` classifies the BioSample records NCBI published on 17 October 2026, and its paper rollup is compared with the curated class.

```bash
onehealth classify tests/data/imp_carbapenemase_sources/samples.tsv -o work/imp/classified.tsv
onehealth validate work/imp/classified.tsv \
    --truth tests/data/imp_carbapenemase_sources/samples.tsv -o work/imp/validation.tsv
```

| | Lexicon v2 | Lexicon v3 |
|---|---|---|
| Coverage | 84.5% | 84.8% |
| Accuracy | 99.2% | 99.5% |
| Agreement | 83.6% | 84.2% |
| Human recall / precision | 92.0% / 99.1% | 92.2% / 99.4% |
| Environmental recall / precision | 13.5% / 100% | 17.6% / 100% |
| Animal recall / precision | 74.5% / 86.4% | 74.5% / 86.4% |
| Sent to review | 674 | 661 |

The confusion matrix with v3:

| Curated \ onehealth | Human | Environmental | Animal | Unclassified |
|---|---|---|---|---|
| Human | 3,547 | 0 | 5 | 295 |
| Environmental | 12 | 78 | 1 | 353 |
| Animal | 0 | 0 | 38 | 13 |
| Unclassified | 9 | 0 | 0 | 0 |

### What the numbers mean

- **Most of the gap is missing information.** 555 of the 661 rows sent to review have no isolation source and no host in their BioSample record; the curators had other context. 87 have a source no rule matches, 17 a host that contradicts the source, and 2 a host the lexicon does not recognise.
- **Environmental recall is low** because most environmental isolates in this study have no isolation source in their BioSample. None of the 78 called Environmental is wrong.
- **12 environmental isolates are called Human.** Nine are room drains and swabs without a setting word (`Rm 14 Drain`, `Sink swab`) recorded with host *Homo sapiens*: no rule matches the source, so the host decides. Three are stool samples with host *Homo sapiens*.
- **9 drain samples** the curators could not classify are called Human, also from host *Homo sapiens*.

## NCBI Pathogen Detection *Klebsiella*

307 isolates in 22 SNP clusters from NCBI Pathogen Detection (release `PDG000000012.2554`, 17 October 2026), with BioSample metadata exactly as NCBI publishes it. There are no curated labels, but NCBI's own `epi_type` (clinical, or environmental/other) is a coarse check.

```bash
onehealth classify tests/data/ncbi_pd_klebsiella/metadata.tsv -o work/pd/classified.tsv
```

| | Lexicon v2 | Lexicon v3 |
|---|---|---|
| Classified | 253 of 307 | 263 of 307 |

With v3:

| NCBI `epi_type` | Isolates | Human | Environmental | Animal | Review |
|---|---|---|---|---|---|
| clinical | 234 | 233 | 0 | 0 | 1 |
| environmental/other | 35 | 1 | 28 | 1 | 5 |
| none | 38 | 0 | 0 | 0 | 38 |

No clinical isolate is called Environmental or Food. The Animal isolate is a faecal swab from a pig, which NCBI files under environmental/other. The Human one is `feces` with no host, the `+no_host` default. The 38 isolates without an `epi_type` have no isolation source or host either, so they go to review as `missing_source`.

The same dataset, with its AMRFinderPlus genes and SNP distances, is what the test suite runs `report` on.
