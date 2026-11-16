# What the classes mean, and where onehealth follows NCBI

!!! warning "Read this before you compare sources"
    A One Health class is a definition, not a measurement. Whether manure is an animal sample or an environmental one changes how many isolates are in each group, and so changes every comparison between groups. `onehealth` does not invent its own definitions where NCBI's records show a convention: **it follows NCBI**. This page says exactly where, with the numbers behind each choice, and how to choose differently for your own study.

There are four such conventions. The first three were adopted in lexicon version 7 (November 2026), and version 8 widened the third from `chicken` to every meat species; results made with an earlier version use different definitions and are not comparable class for class. The fourth is older, and was measured against NCBI's records afterwards and kept.

| # | Convention | In short |
|---|---|---|
| 1 | [A record's own `source_type` is its class](#1-a-records-own-source_type-is-its-class) | Where the submitter told NCBI the One Health category, that category is used as it is. |
| 2 | [Farm settings are environmental, and a farm with a host is the host](#2-farm-settings) | Manure, litter, a barn and a farm are Environmental. Until version 7 they were Animal. |
| 3 | [A meat species with no host is meat](#3-a-species-that-is-also-the-name-of-its-meat) | With no host, `chicken`, `turkey`, `fish`, `poultry`, `salmon`, `beef` or `lamb` in the text is Food. Until version 7 they were Animal; version 8 applies it to the word anywhere in the text. |
| 4 | [A specimen with no host is a person's](#4-a-specimen-with-no-host-is-a-persons) | `blood` or `stool` with no host named is Human. Unchanged; NCBI's records agree for 98.0%. |

The numbers on this page come from 139,554 NCBI Pathogen Detection records whose submitters filled in NCBI's `source_type`, from 18 organism groups (`tests/data/ncbi_pd_source_types/`; a further 26,604 records are kept apart to test the result). The test suite checks every figure quoted here against that data (`tests/test_ncbi_conventions.py`).

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

Version 7 made these two whole values Food. The same holds when the word is part of a longer text, and for the other species sold as meat. Of the benchmark records with no host that this rule decides, by the word it finds:

| Word in the text | Records | Food |
|---|---|---|
| `chicken` (`chicken liver`, `packaged chicken`) | 2,194 | 2,099 (95.7%) |
| `turkey` | 313 | 312 (99.7%) |
| `beef` | 224 | 219 (97.8%) |
| `poultry` | 171 | 165 (96.5%) |
| `fish` | 162 | 161 (99.4%) |
| `salmon` | 56 | 55 (98.2%) |
| `animal-cattle-...`, a product class of the US meat inspection service | 53 | 51 (96.2%) |
| `lamb`, `bison` | 17 | 17 |
| All | 3,190 | 3,079 (96.5%) |

**What onehealth does.** With no host, one of these words in `isolation_source` makes the row **Food**, with the reason `meat-species-no-host` so that you can count them. So `fish`, `poultry` and `chicken liver` are Food; until version 8 they were Animal. With a host, the host decides as before: `chicken` with host *Gallus gallus* is Animal, and so is `animal-cattle-steer` with host bovine.

**What still means the animal**, with no host:

- anything that says what was sampled: a specimen (`chicken feces`), a gut site (`chicken cecal content`, `chicken gut`), `live chicken`, `beef cattle`;
- a place or what animals leave there, which is Environmental (convention 2): `poultry litter`, `fish farm`;
- the words without this evidence: `duck`, `pig`, `cow`, `cattle`, `swine`, `bovine`, `wild bird`. A bare `cattle` with no host is an animal in all 46 NCBI records.

**There is no agreed standard for these words.** The standards for sample metadata (NCBI's packages, the PHA4GE specification, the FoodOn ontology) do not say what a bare `chicken` means; they ask for it not to be left bare: a `host` for an animal, a food product term for meat, a `source_type`. `onehealth` goes with what NCBI's records show. A veterinary study whose `chicken` is a bird should say so with a [decision](../commands/classify.md#decisions): one line, `chicken` with `l1` Animal.

## 4. A specimen with no host is a person's

`blood`, `stool`, `urine` and `wound` name a specimen and say nothing about whose it is. Many clinical records give no host, because to the submitter it goes without saying. So when the text is a specimen, no host is given and nothing in the text names an animal, `onehealth` calls the row **Human**, with a reason that ends in `+no_host` so that you can count these rows.

This default is older than version 7 and was `onehealth`'s own judgement. In November 2026 it was measured against NCBI's categories, like the conventions above. In the benchmark records that `classify` places this way from the text alone:

| | Records | Human | Animal | Environmental | Food |
|---|---|---|---|---|---|
| A specimen, no host | 15,290 | 14,982 (98.0%) | 271 | 21 | 16 |

So the default stands. NCBI's submitters mean a person by a host-less specimen in 49 records of 50.

**Where the default is not applied.**

- **The text names an animal.** `pig feces` always was Animal. Since version 8 any name NCBI Taxonomy gives to an animal counts, so `stool (Odocoileus hemionus)`, `pigeon stool` and `giant panda feces` are Animal too (reason `specimens+animal-name-in-specimen`); version 7 called them Human. 79 of the 81 benchmark records placed this way are animal in NCBI.
- **The text names a food.** `pork blood sausage` and `oyster homogenized tissue dilution` go to review as `specimen_or_food`. The second is shellfish from a monitoring programme, environmental or food in NCBI; version 7 called it Human.
- **The record was filed as environmental, food or other** (the BioSample package `Pathogen.env`, or `epi_type` `environmental/other`). A host-less specimen there goes to review as `no_host_env_package`. In a sample of 159,045 NCBI records there are 266 of these; of the 45 that also have a `source_type`, 34 are animal, 10 human and 1 environmental, so no default would be right.

**What the 2% are.** Almost all are the faeces or tissue of an animal, submitted with no host and no word for the animal (`feces`, `stool`, `internal tissue liver`). No text classifier can tell these from a patient's; a `source_type` or a `host` can.

## What is still onehealth's own judgement

NCBI has no convention for these, so they are choices, and each can be counted in the output:

- **`healthcare_environment`.** A hospital fixture is flagged only when the text says where it is (`hospital sink`, `ICU drain`), or is a sink or shower drain. A floor drain is not flagged: in NCBI records it is a food plant's. The whole value `hospital` or `healthcare facility` is not a fixture: it is where a patient isolate came from, and is Human (`hospital`: 328 NCBI records of 328).
- **A host that is a plant goes to review** (`host_unrecognised`). Only animals and people are placed.
- **`Laboratory`** is a class of its own, and rolls up to Environmental with Food.
- **A fixture named without its setting** (`Sink swab`, `Rm 14 Drain`) goes to review.

## How well this works

On 26,604 benchmark records that were drawn after version 8 was finished, `classify` reads the free text alone (with `source_type` hidden), classifies 97.8% of them and matches NCBI's category for 96.7% of those; see [Benchmarks](../benchmarks.md#ncbi-source_type-held-out).
