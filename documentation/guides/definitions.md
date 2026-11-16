# What the classes mean, and where onehealth follows NCBI

!!! warning "Read this before you compare sources"
    A One Health class is a definition, not a measurement. Whether manure is an animal sample or an environmental one changes how many isolates are in each group, and so changes every comparison between groups. `onehealth` does not invent its own definitions where NCBI's records show a convention: **it follows NCBI**. This page says exactly where, with the numbers behind each choice, and how to choose differently for your own study.

There are three such conventions. All three were adopted in lexicon version 7 (November 2026); results made with an earlier version use different definitions and are not comparable class for class.

| # | Convention | In short |
|---|---|---|
| 1 | [A record's own `source_type` is its class](#1-a-records-own-source_type-is-its-class) | Where the submitter told NCBI the One Health category, that category is used as it is. |
| 2 | [Farm settings are environmental, and a farm with a host is the host](#2-farm-settings) | Manure, litter, a barn and a farm are Environmental. Until version 7 they were Animal. |
| 3 | [`chicken` or `turkey` on its own is meat](#3-a-species-that-is-also-the-name-of-its-meat) | With no host these two words are Food. Until version 7 they were Animal. |

The numbers on this page come from 139,554 NCBI Pathogen Detection records whose submitters filled in NCBI's `source_type`, from 18 organism groups (`tests/data/ncbi_pd_source_types/`). The test suite checks every figure quoted here against that data (`tests/test_ncbi_conventions.py`).

## The order of evidence

`classify` takes the class of a row from the best evidence the row has, in this order:

1. **A curator's override** for that isolate ([overrides](../commands/classify.md#overrides)).
2. **A curator's decision** for its source text ([decisions](../commands/classify.md#decisions)).
3. **NCBI's `source_type`**, where the record has one that names a class (convention 1).
4. **The lexicon**, reading the free text of `isolation_source`, `host` and the other fields.

So everything below can be overruled for your study, per isolate or per text, without changing `onehealth`.

## 1. A record's own `source_type` is its class

NCBI's `source_type` is, in NCBI's words, a "controlled vocabulary describing the isolation_source. Choose the best fit term: Human, Animal, Food, Environmental, Other". It belongs to the BioSample package for One Health enteric pathogens, and it is chosen by the submitter, who knows what was sampled.

**What onehealth does.** When the `source_type` column of a row is `human`, `animal`, `food` or `environmental`, that is the row's class. The lexicon reads only free text, so it is not allowed to disagree with the record's own category, and it is not allowed to send such a row to review.

| The free text... | Class | `classify_reason` | `lexicon_suggestion` |
|---|---|---|---|
| is silent, or agrees | NCBI's | the lexicon rule that agrees, as before | the same class |
| cannot settle it: host and source conflict, or a rule sends it to review | NCBI's | `source_type` | the classes the text left open, e.g. `Animal|Food` |
| says something else | NCBI's | `source_type_over_text` | what the lexicon would have said |

Two reasons, so that you can count each and look at the rows: in a sample of 159,045 NCBI records, 51,288 (32%) have a `source_type` that names a class; 456 of them are `source_type` and 141 are `source_type_over_text`.

**What changed.** Until version 7, `source_type` was one more text field. A word in `isolation_source` could outrank it: `raw milk` with `source_type` environmental was Food. A host that pointed elsewhere sent the row to review: `animal-cattle-steer` with host bovine and `source_type` food. Both now take NCBI's category.

**What to keep in mind.**

- **Only a third of NCBI records have it.** For the rest, the lexicon decides, and conventions 2 and 3 apply.
- **Any other value declares nothing**: `other`, `not provided`, an empty cell. So does a `tissue` or NCBI `sample_type` column, which are free text.
- **Submitters make mistakes**, and `onehealth` now repeats them: the liver of a bald eagle labelled `food`, a deer's intestine labelled `human`. The reason `source_type_over_text` is where to look for them.
- **The healthcare flag is still the lexicon's.** NCBI has no such category. `hospital sink` with `source_type` environmental is Environmental with `healthcare_environment` yes.
- **`resolve` writes `source_type` in a column of its own**, as well as inside `sample_type`, so that `classify` can tell the submitter's category from a tissue or a material type.

## 2. Farm settings

Until version 7, `manure`, `litter`, `barn` and `farm` were Animal. NCBI submitters record them differently, and `onehealth` now follows NCBI.

**What animals leave behind is environmental.** In the benchmark records:

| Text has | Records | Environmental | Animal |
|---|---|---|---|
| `manure` | 747 | 646 (86.5%) | 101 |
| `manure`, no host given | 609 | 567 (93.1%) | 42 |
| `litter` | 77 | 77 | 0 |
| `bedding` | 68 | 68 | 0 |

`onehealth` calls manure, litter and bedding **Environmental**, with or without a host: `dairy cattle manure` with host *Bos taurus* is Environmental, and the host names the animals it came from. Be aware that with a host NCBI's records are divided (79 environmental, 59 animal); `onehealth` goes with the majority.

**A place is the place, unless a host says what was sampled there.**

| Text has | Records | Environmental | Animal | Other |
|---|---|---|---|---|
| `farm`, no host given | 187 | 141 (75.4%) | 16 | 30 food |
| `farm`, with a host | 641 | 125 | 507 (79.1%) | 9 human |
| `barn` | 19 | 19 | 0 | |
| `poultry shed` | 69 | 69 | 0 | |
| `drag swab` | 152 | 152 | 0 | |

So a farm, a barn, a coop or a poultry shed with **no host is Environmental**: the sample is of the place. With **a host, the host decides**: `pig farm` with host *Sus scrofa* is Animal, and `farm` with host *Homo sapiens* is Human. Swabs of a poultry house (`drag swab`, `boot swab`) are Environmental.

A specimen still follows its host, as everywhere: `feces, calf barn` with host *Bos taurus* is Animal, and `farm worker stool` is Human.

**If your study defines these differently**, say so with [decisions](../commands/classify.md#decisions): one line for `manure` makes every manure sample Animal.

## 3. A species that is also the name of its meat

`chicken` names a bird and its meat. NCBI asks for a `host` when a sample comes from an animal ("the natural host to the organism from which the sample was obtained"), so a species word in `isolation_source` with no host is usually not the animal. In the benchmark records with **no host**:

| `isolation_source` is exactly | Records | Food | Animal |
|---|---|---|---|
| `chicken` | 867 | 810 (93.4%) | 56 |
| `turkey` | 256 | 255 (99.6%) | 1 |

**What onehealth does.** With no host, `chicken` and `turkey` on their own are **Food**, with the reason `meat-species-no-host` so that you can count them. With a host, the host decides as before: `chicken` with host *Gallus gallus* is Animal.

Only these two words are treated this way, because only they have this much evidence. `duck`, `fish`, `pig` and `cow` on their own are still Animal, and so is `chicken` with another word that does not make it a cut of meat (`live chicken`, `chicken cecum`).

**There is no agreed standard for this word.** The standards for sample metadata (NCBI's packages, the PHA4GE specification, the FoodOn ontology) do not say what a bare `chicken` means; they ask for it not to be left bare: a `host` for an animal, a food product term for meat, a `source_type`. `onehealth` goes with what NCBI's records show. A veterinary study whose `chicken` is a bird should say so with a [decision](../commands/classify.md#decisions): one line, `chicken` with `l1` Animal.

## What is still onehealth's own judgement

NCBI has no convention for these, so they are choices, and each can be counted in the output:

- **`healthcare_environment`.** A hospital fixture is flagged only when the text says where it is (`hospital sink`, `ICU drain`), or is a sink or shower drain. A floor drain is not flagged: in NCBI records it is a food plant's.
- **A specimen with no host is Human**, with the reason ending in `+no_host`. Not when the record was filed as environmental, food or other (`no_host_env_package`), and not when the same text names a food (`specimen_or_food`); those go to review.
- **A host that is a plant goes to review** (`host_unrecognised`). Only animals and people are placed.
- **`Laboratory`** is a class of its own, and rolls up to Environmental with Food.
- **A fixture named without its setting** (`Sink swab`, `Rm 14 Drain`) goes to review.

## How well this works

On 34,252 benchmark records that were set aside until version 7 was finished, `classify` reads the free text alone (with `source_type` hidden) and matches NCBI's category for 96.5% of the records it classifies; see [Benchmarks](../benchmarks.md#ncbi-source_type-held-out).
