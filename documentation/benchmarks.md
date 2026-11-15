# Benchmarks

How well does `classify` do on real data? Three public datasets ship with `onehealth` in `tests/data/`, and the test suite fails if a change makes the results below worse.

!!! note
    The first two datasets are *Klebsiella* and other hospital organisms, and the lexicon was tuned on them. Their numbers show it works on the data it was built from; they are not a measure of accuracy on data it has not seen. The third dataset has a [held-out part](#ncbi-source_type-held-out) for that. Run [`validate`](commands/validate.md) on labels from your own study before relying on the classes.

All figures are for the lexicon that ships, version 6. What each earlier version scored is in the [changelog](changelog.md).

## IMP carbapenemase isolates

4,351 public isolates from a genomic epidemiology study of IMP carbapenemases, whose authors curated the source of every isolate by hand (Vezina et al. 2025, *Nat Commun* 17:183, Supplementary Table S1). `onehealth` classifies the BioSample records NCBI published on 17 October 2026, and its paper rollup is compared with the curated class.

```bash
onehealth classify tests/data/imp_carbapenemase_sources/samples.tsv -o work/imp/classified.tsv
onehealth validate work/imp/classified.tsv \
    --truth tests/data/imp_carbapenemase_sources/samples.tsv -o work/imp/validation.tsv
```

| | |
|---|---|
| Coverage | 84.8% |
| Accuracy | 99.8% |
| Agreement | 84.7% |
| Human recall / precision | 92.1% / 99.9% |
| Environmental recall / precision | 17.8% / 100% |
| Animal recall / precision | 100% / 89.5% |
| Sent to review | 667 |

The confusion matrix:

| Curated \ onehealth | Human | Environmental | Animal | Unclassified |
|---|---|---|---|---|
| Human | 3,545 | 0 | 5 | 297 |
| Environmental | 3 | 79 | 1 | 361 |
| Animal | 0 | 0 | 51 | 0 |
| Unclassified | 0 | 0 | 0 | 9 |

### What the numbers mean

- **Most of the gap is missing information.** 555 of the 667 rows sent to review have no isolation source and no host in their BioSample record; the curators had other context. 75 name a fixture without its setting (`Shower 3`, `Rm 14 Drain`), 18 have a source no rule matches, 17 a host that contradicts the source, and 2 a host the lexicon does not recognise.
- **Every animal isolate is found.** 13 of the 51 have a host that no word list named (11 American crows and 2 silver gulls, *Chroicocephalus novaehollandiae*). They went to review until v6, which recognises any animal host by its name in NCBI Taxonomy.
- **Environmental recall is low** because most environmental isolates in this study have no isolation source in their BioSample. None of the 79 called Environmental is wrong.
- **3 environmental isolates are called Human.** All three are stool samples with host *Homo sapiens*.
- **A fixture without its setting goes to review.** 9 room drains and swabs (`Rm 14 Drain`, `Sink swab`) have host *Homo sapiens*, and the curators say they are environmental; 9 more drain samples the curators could not classify. All 18 go to review and not to Human, and two human samples go with them (`Drain Swab`, `drain inserting site`), because a fixture word alone cannot tell them apart.

## NCBI Pathogen Detection *Klebsiella*

307 isolates in 22 SNP clusters from NCBI Pathogen Detection (release `PDG000000012.2554`, 17 October 2026), with BioSample metadata exactly as NCBI publishes it. There are no curated labels, but NCBI's own `epi_type` (clinical, or environmental/other) is a coarse check.

```bash
onehealth classify tests/data/ncbi_pd_klebsiella/metadata.tsv -o work/pd/classified.tsv
```

262 of the 307 isolates are classified:

| NCBI `epi_type` | Isolates | Human | Environmental | Animal | Review |
|---|---|---|---|---|---|
| clinical | 234 | 233 | 0 | 0 | 1 |
| environmental/other | 35 | 0 | 28 | 1 | 6 |
| none | 38 | 0 | 0 | 0 | 38 |

No clinical isolate is called Environmental or Food. The Animal isolate is a faecal swab from a pig, which NCBI files under environmental/other. The 38 isolates without an `epi_type` have no isolation source or host either, so they go to review as `missing_source`.

One isolate is `feces` with no host, filed by its submitter as environmental/other. The `+no_host` default would call it Human; it goes to review as `no_host_env_package`; see [`classify`](commands/classify.md#specimens-follow-the-host).

The same dataset, with its AMRFinderPlus genes and SNP distances, is what the test suite runs `report` on.

## NCBI `source_type`, held out

The two datasets above are hospital organisms, and most of their isolates are from patients. This one asks how the lexicon does on food, animal and environmental sources, on records it was not written from.

Some submitters to NCBI Pathogen Detection fill in `source_type`, NCBI's own One Health category (`human`, `animal`, `food`, `environmental`). `classify` reads it when it is there, so here it is hidden: each record is classified from `isolation_source` and `host` alone, and the L1 class is compared with the label. That is the situation of the two NCBI records in three that have no `source_type`.

The fixture holds 73,094 such records from ten organism groups, in two parts (`tests/data/ncbi_pd_source_types/`):

- **Development**, 33,372 records from six groups (Salmonella, *E. coli* and Shigella, Listeria, Campylobacter, Klebsiella, Acinetobacter). Lexicon v5 was written while looking at these.
- **Held out**, 39,722 records: other parts of the same six tables, and four groups not in the development part (*Staphylococcus aureus*, *Vibrio parahaemolyticus*, *Bacillus cereus* group, Cronobacter). Nobody looked at these before v5 was finished, and v5 was not changed afterwards.

| Held out | v4 | v5 | Now (v6) |
|---|---|---|---|
| Coverage | 76.9% | 88.7% | 89.1% |
| Accuracy | 80.9% | 93.9% | 93.9% |
| Human precision / recall | 98.0% / 98.6% | 98.0% / 98.6% | 98.0% / 98.6% |
| Animal precision / recall | 45.4% / 89.0% | 73.7% / 88.3% | 74.3% / 91.7% |
| Food precision / recall | 92.6% / 30.0% | 97.2% / 70.6% | 97.3% / 70.6% |
| Environmental precision / recall | 99.6% / 43.8% | 99.0% / 77.6% | 99.0% / 77.6% |

On the development part, v6 classifies 94.3% and gets 95.0% of those right (v5: 93.9% and 94.9%; v4: 83.9% and 78.5%). The v4 and v5 columns are what the lexicon scored when each was current; the tests hold v6 to its own column.

v5 was written from the development part and never saw the held-out part. v6 was not written from either: its rules are general (any animal NCBI Taxonomy names, and the word `human` in a specimen), and were checked on the development part and on NCBI records with no label. By then the held-out figures for v5 had been read once, as a total and a list of its commonest errors, so the held-out part is no longer unseen in the strict sense.

The confusion matrix, held out:

| NCBI `source_type` \ onehealth | Human | Animal | Food | Environmental | Review |
|---|---|---|---|---|---|
| human | 12,596 | 52 | 0 | 0 | 128 |
| animal | 105 | 4,482 | 61 | 1 | 238 |
| food | 41 | 1,157 | 9,535 | 69 | 2,710 |
| environmental | 117 | 342 | 204 | 6,634 | 1,250 |

### What the numbers mean

- **v4 called retail meat Animal.** `raw intact chicken`, `ground turkey` and `chicken breast` matched the animal word and nothing else, so fewer than half of v4's Animal calls agreed with NCBI. v5 reads how the meat is sold or processed as Food.
- **v4 called food-plant swabs Food.** `environmental swab sponge non food contact surface` contains the word food. v5 reads a swab of a surface or of food premises as Environmental.
- **v5 left animals with an unlisted host in review.** v6 places any animal NCBI Taxonomy names, which moves 173 held-out records from review to Animal, 166 of them labelled animal by NCBI. It also sends 12 to review that v5 called Food from a word in the source: with the host recognised, host and source now disagree.
- **A Human or Environmental call is reliable**, 98% and 99% of the time, in every version.
- **Animal is still the weakest class.** Most of the 1,157 food records called Animal are a species on its own (`chicken`, `turkey`), which the lexicon reads as the animal; chicken organs (`chicken liver`), which can be giblets or a necropsy sample; and caecal samples taken at slaughter, which the US Food Safety and Inspection Service labels food. Most of the 342 environmental records called Animal are manure and farm buildings, which NCBI submitters label environmental and `onehealth` calls a farm setting. Whether those are errors depends on your definitions; `source_type` is the submitter's choice, not a curated truth: among the food records v6 newly calls Animal are livers of a bald eagle, two bats and a raccoon.
- **Accuracy differs by organism.** Held out, v6 is right for 99% of classified Cronobacter and Acinetobacter records, 97% of Salmonella, 94% of Campylobacter, 92% of Listeria and *E. coli*, and 88% of *Vibrio parahaemolyticus*, where seafood is both the animal and the food.
