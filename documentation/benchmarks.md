# Benchmarks

How well does `classify` do on real data? Three public datasets ship with `onehealth` in `tests/data/`, and the test suite fails if a change makes the results below worse.

!!! note
    The first two datasets are *Klebsiella* and other hospital organisms, and the lexicon was tuned on them. Their numbers show it works on the data it was built from; they are not a measure of accuracy on data it has not seen. The third dataset has a [part set aside](#ncbi-source_type-held-out) for that. Run [`validate`](commands/validate.md) on labels from your own study before relying on the classes.

All figures are for the lexicon that ships, version 8. What each earlier version scored is in the [changelog](changelog.md).

## IMP carbapenemase isolates

4,351 public isolates from a genomic epidemiology study of IMP carbapenemases, whose authors curated the source of every isolate by hand (Vezina et al. 2025, *Nat Commun* 17:183, Supplementary Table S1). `onehealth` classifies the BioSample records NCBI published on 17 October 2026, and its paper rollup is compared with the curated class.

```bash
onehealth classify tests/data/imp_carbapenemase_sources/samples.tsv -o work/imp/classified.tsv
onehealth validate work/imp/classified.tsv \
    --truth tests/data/imp_carbapenemase_sources/samples.tsv -o work/imp/validation.tsv
```

| | |
|---|---|
| Coverage | 84.9% |
| Accuracy | 99.8% |
| Agreement | 84.7% |
| Human recall / precision | 92.2% / 99.9% |
| Environmental recall / precision | 18.0% / 100% |
| Animal recall / precision | 100% / 91.1% |
| Sent to review | 665 |

The confusion matrix:

| Curated \ onehealth | Human | Environmental | Animal | Unclassified |
|---|---|---|---|---|
| Human | 3,546 | 0 | 5 | 296 |
| Environmental | 4 | 80 | 0 | 360 |
| Animal | 0 | 0 | 51 | 0 |
| Unclassified | 0 | 0 | 0 | 9 |

### What the numbers mean

- **Most of the gap is missing information.** 555 of the 665 rows sent to review have no isolation source and no host in their BioSample record; the curators had other context. 74 name a fixture without its setting (`Shower 3`, `Rm 14 Drain`), 18 have a host that contradicts the source, 16 a source no rule matches, and 2 a host the lexicon does not recognise.
- **Every animal isolate is found.** 13 of the 51 have a host that no word list named (11 American crows and 2 silver gulls, *Chroicocephalus novaehollandiae*). They went to review until v6, which recognises any animal host by its name in NCBI Taxonomy.
- **Environmental recall is low** because most environmental isolates in this study have no isolation source in their BioSample. None of the 80 called Environmental is wrong. One of them is `manure` with a pig host, which the curators call environmental, as NCBI's records do; until v7 it was Animal.
- **4 environmental isolates are called Human.** Three are stool samples with host *Homo sapiens*. The fourth, since v8, is described only as `hospital`: the curators place it in the hospital environment, and in NCBI's records that value is a patient's in 328 records of 328, so `onehealth` calls it Human.
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

Some submitters to NCBI Pathogen Detection fill in `source_type`, NCBI's own One Health category (`human`, `animal`, `food`, `environmental`). In normal use that field is the class ([Definitions](guides/definitions.md#1-a-records-own-source_type-is-its-class)), so here it is hidden: each record is classified from `isolation_source` and `host` alone, and the L1 class is compared with the label. That is the situation of the two NCBI records in three that have no `source_type`.

The fixture holds 166,158 such records from 26 organism groups, in five parts that share no record (`tests/data/ncbi_pd_source_types/`):

| Part | Records | What it was used for |
|---|---|---|
| Development | 33,372 | Lexicon v5 and v7 were written while looking at these. |
| Held out | 39,722 | Unseen by v5. v7 was checked against it, and three rules were narrowed. |
| Validation | 32,208 | Meant to be read once after v7. It showed that one v7 rule went against NCBI's records (a farm with an animal host is the animal), and the rule was corrected. |
| Final | 34,252 | Read once, after v7 was finished: the figure to quote for v7. Lexicon v8 was written with all four parts above in view. |
| **Fresh** | **26,604** | **Drawn after v8 was finished, and read once. These are the figures to quote.** |

Each part after the first adds organism groups the earlier parts do not have. The fresh part adds *Vibrio alginolyticus*, *Enterococcus faecalis*, *Enterococcus faecium*, *Citrobacter freundii*, *Aeromonas veronii*, *Pseudomonas aeruginosa*, *Enterobacter hormaechei* and Serratia, chosen by name before any of their records were seen.

| Fresh part | Before (v7) | Now (v8), as first read |
|---|---|---|
| Coverage | 97.1% | 97.8% |
| Accuracy | 95.4% | 96.7% |
| Human precision / recall | 98.8% / 95.8% | 98.8% / 95.9% |
| Animal precision / recall | 90.7% / 91.2% | 97.4% / 92.0% |
| Food precision / recall | 98.7% / 91.7% | 98.2% / 97.4% |
| Environmental precision / recall | 81.4% / 84.2% | 81.4% / 85.1% |

The confusion matrix with v8, fresh part, as first read:

| NCBI `source_type` \ onehealth | Human | Animal | Food | Environmental | Review |
|---|---|---|---|---|---|
| human | 10,515 | 2 | 0 | 345 | 81 |
| animal | 119 | 5,239 | 90 | 63 | 186 |
| food | 1 | 106 | 7,595 | 7 | 87 |
| environmental | 9 | 30 | 47 | 1,821 | 234 |

27 human records called Laboratory are not shown.

**One correction was made after this reading.** 296 of the 345 human records called Environmental have one description, `healthcare facility`. Version 7 had taken that for a hospital fixture, with no labelled record to go by; these 296 are all human to their submitter. It is now read like `hospital`: the place a patient isolate came from. With that one change the lexicon that ships scores 97.8% coverage and 97.9% accuracy on this part (Environmental precision 93.9%, Human recall 98.6%), but the part has now been seen, so the table above is the fair test. Nothing else was changed after the reading.

The other parts, with v8. They are not independent of v8, which was written with them in view:

| | Development | Held out | Validation | Final |
|---|---|---|---|---|
| Coverage | 98.0% | 97.1% | 98.8% | 98.3% |
| Accuracy | 98.5% | 97.9% | 98.1% | 98.3% |
| Animal precision / recall | 93.3% / 95.0% | 92.4% / 95.1% | 94.6% / 97.4% | 97.2% / 97.8% |
| Food precision / recall | 99.0% / 94.8% | 98.4% / 92.9% | 99.2% / 96.0% | 98.5% / 94.4% |
| Environmental precision / recall | 99.5% / 95.6% | 98.3% / 91.9% | 99.1% / 89.4% | 97.7% / 95.1% |

### How the lexicon got here

On the held-out part, the one every version since v4 can be compared on:

| Held out | v4 | v5 | v6 | v7 | v8 |
|---|---|---|---|---|---|
| Coverage | 76.9% | 88.7% | 89.1% | 90.3% | 97.1% |
| Accuracy | 80.9% | 93.9% | 93.9% | 96.1% | 97.9% |
| Animal precision | 45.4% | 73.7% | 74.3% | 84.2% | 92.4% |
| Environmental recall | 43.8% | 77.6% | 77.6% | 87.3% | 91.9% |

- **v4 called retail meat Animal** (`raw intact chicken`, `ground turkey`) and food-plant swabs Food. v5 fixed both.
- **v5 left animals with an unlisted host in review.** v6 places any animal NCBI Taxonomy names.
- **v6 defined farm settings and a bare `chicken` its own way.** v7 follows NCBI: manure and a farm are environmental, a farm with a host is the host, and `chicken` on its own with no host is food ([Definitions](guides/definitions.md)).
- **v7 still called `fish`, `poultry` or `chicken liver` with no host Animal, and left `oyster` and `potato` in review.** v8 treats every meat species as it treated `chicken`, and names the foods: Food recall on this part went from 73.9% to 92.9%.

### What the numbers mean

- **A Human, Animal or Food call is reliable**: 99%, 97% and 98% of them match NCBI's label on the fresh part.
- **Animal was the weakest class, and is not any more.** Until v8, most wrong calls were food records called Animal by an animal word (`chicken liver`, `fish`, `poultry`). Those are Food now, and Animal precision went from 90.7% to 97.4%. What remains are the slaughter classes of the US Food Safety and Inspection Service with a host (`animal-cattle-steer`, host bovine), which NCBI labels food in a minority of records.
- **Environmental is now the weakest**, because of the one description above. Without it, 94% of Environmental calls match.
- **Animal recall is 92%.** 119 animal records are called Human: faeces with no host and no word for the animal, which no text can tell from a patient's. 186 go to review, most of them an animal named in `isolation_source` with no host (`Larus cachinnans`, `gull`, `fox`). A species name alone in that field is not taken for a description of a sample.
- **`source_type` is the submitter's choice, not a curated truth.** Some of what counts as a wrong call here is a wrong label: the held-out part has livers of a bald eagle, two bats and a raccoon labelled food.
- **Accuracy differs by organism.** On the fresh part v8 was right for 100% of classified Serratia, *A. veronii* and *E. hormaechei* records, 99% of Campylobacter, 98% of Salmonella, Klebsiella and *E. faecium*, 92% of *E. coli* (98% after the correction), 91% of *E. faecalis* and 88% of *C. freundii*.
- **These figures are for free text alone.** In normal use a record with a `source_type` takes it as its class, so the lexicon's errors apply only to records without one.
