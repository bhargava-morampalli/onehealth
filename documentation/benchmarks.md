# Benchmarks

How well does `classify` do on real data? Three public datasets ship with `onehealth` in `tests/data/`, and the test suite fails if a change makes the results below worse.

!!! note
    The first two datasets are *Klebsiella* and other hospital organisms, and the lexicon was tuned on them. Their numbers show it works on the data it was built from; they are not a measure of accuracy on data it has not seen. The third dataset has a [held-out part](#ncbi-source_type-held-out) for that. Run [`validate`](commands/validate.md) on labels from your own study before relying on the classes.

All figures are for the lexicon that ships, version 7. What each earlier version scored is in the [changelog](changelog.md).

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
| Sent to review | 666 |

The confusion matrix:

| Curated \ onehealth | Human | Environmental | Animal | Unclassified |
|---|---|---|---|---|
| Human | 3,546 | 0 | 5 | 296 |
| Environmental | 3 | 80 | 0 | 361 |
| Animal | 0 | 0 | 51 | 0 |
| Unclassified | 0 | 0 | 0 | 9 |

### What the numbers mean

- **Most of the gap is missing information.** 555 of the 666 rows sent to review have no isolation source and no host in their BioSample record; the curators had other context. 74 name a fixture without its setting (`Shower 3`, `Rm 14 Drain`), 18 have a host that contradicts the source, 17 a source no rule matches, and 2 a host the lexicon does not recognise.
- **Every animal isolate is found.** 13 of the 51 have a host that no word list named (11 American crows and 2 silver gulls, *Chroicocephalus novaehollandiae*). They went to review until v6, which recognises any animal host by its name in NCBI Taxonomy.
- **Environmental recall is low** because most environmental isolates in this study have no isolation source in their BioSample. None of the 80 called Environmental is wrong. One of them is `manure` with a pig host, which the curators call environmental, as NCBI's records do; until v7 it was Animal.
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

Some submitters to NCBI Pathogen Detection fill in `source_type`, NCBI's own One Health category (`human`, `animal`, `food`, `environmental`). In normal use that field is the class ([Definitions](guides/definitions.md#1-a-records-own-source_type-is-its-class)), so here it is hidden: each record is classified from `isolation_source` and `host` alone, and the L1 class is compared with the label. That is the situation of the two NCBI records in three that have no `source_type`.

The fixture holds 139,554 such records from 18 organism groups, in four parts that share no record (`tests/data/ncbi_pd_source_types/`):

| Part | Records | What it was used for |
|---|---|---|
| Development | 33,372 | Lexicon v5 and v7 were written while looking at these. |
| Held out | 39,722 | Unseen by v5. v7 was checked against it, and three rules were narrowed. |
| Validation | 32,208 | Meant to be read once after v7. It showed that one v7 rule went against NCBI's records (a farm with an animal host is the animal), and the rule was corrected. |
| **Final** | **34,252** | **Read once, after v7 was finished. v7 was not changed afterwards. These are the figures to quote.** |

Each part after the first adds organism groups the earlier parts do not have; the final part adds *Vibrio vulnificus*, *Listeria innocua*, *Staphylococcus pseudintermedius*, *Clostridium botulinum* and Aeromonas.

| Final part | Before (v6) | Now (v7) |
|---|---|---|
| Coverage | 94.6% | 95.5% |
| Accuracy | 95.3% | 96.5% |
| Human precision / recall | 99.2% / 98.5% | 99.2% / 98.5% |
| Animal precision / recall | 88.2% / 97.5% | 91.7% / 97.0% |
| Food precision / recall | 98.5% / 81.3% | 98.6% / 83.3% |
| Environmental precision / recall | 97.4% / 78.7% | 97.7% / 90.2% |

The confusion matrix with v7, final part:

| NCBI `source_type` \ onehealth | Human | Animal | Food | Environmental | Review |
|---|---|---|---|---|---|
| human | 9,336 | 1 | 0 | 59 | 78 |
| animal | 73 | 9,058 | 105 | 0 | 99 |
| food | 0 | 723 | 8,849 | 44 | 1,010 |
| environmental | 6 | 96 | 25 | 4,346 | 344 |

The other parts, with v7. They are higher or lower than the final part because their organisms differ, and they are not independent of v7:

| | Development | Held out | Validation |
|---|---|---|---|
| Coverage | 95.8% | 90.3% | 96.7% |
| Accuracy | 97.5% | 96.1% | 96.1% |
| Animal precision / recall | 85.9% / 92.5% | 84.2% / 90.8% | 86.5% / 93.6% |
| Food precision / recall | 99.1% / 89.0% | 98.3% / 73.9% | 99.2% / 86.5% |
| Environmental precision / recall | 99.5% / 92.9% | 98.2% / 87.3% | 99.1% / 86.4% |

### How the lexicon got here

On the held-out part, the one every version since v4 can be compared on:

| Held out | v4 | v5 | v6 | v7 |
|---|---|---|---|---|
| Coverage | 76.9% | 88.7% | 89.1% | 90.3% |
| Accuracy | 80.9% | 93.9% | 93.9% | 96.1% |
| Animal precision | 45.4% | 73.7% | 74.3% | 84.2% |
| Environmental recall | 43.8% | 77.6% | 77.6% | 87.3% |

- **v4 called retail meat Animal** (`raw intact chicken`, `ground turkey`) and food-plant swabs Food. v5 fixed both.
- **v5 left animals with an unlisted host in review.** v6 places any animal NCBI Taxonomy names.
- **v6 defined farm settings and a bare `chicken` its own way.** v7 follows NCBI: manure and a farm are environmental, a farm with a host is the host, and `chicken` on its own with no host is food ([Definitions](guides/definitions.md)).

### What the numbers mean

- **A Human, Food or Environmental call is reliable**: 99%, 99% and 98% of them match NCBI's label on the final part.
- **Animal is the weakest class**, at 92%. 721 of the 1,132 wrong calls on the final part are food records called Animal by an animal word: a species named with other words that do not mark it as meat, and the slaughter classes of the US Food Safety and Inspection Service (`animal-cattle-steer`), which NCBI labels food.
- **Food recall is 83%** because 1,010 food records go to review: a food no rule names. They are not wrong calls, and one [decision](commands/classify.md#decisions) per food settles them.
- **`source_type` is the submitter's choice, not a curated truth.** Some of what counts as a wrong call here is a wrong label: the held-out part has livers of a bald eagle, two bats and a raccoon labelled food.
- **Accuracy differs by organism.** On the final part v7 is right for 100% of classified Klebsiella and *S. pseudintermedius* records, 98% of Listeria, 97% of Salmonella, 95% of Campylobacter and *L. innocua*, and 93% of *E. coli*.
- **These figures are for free text alone.** In normal use a record with a `source_type` takes it as its class, so the lexicon's errors apply only to records without one.
