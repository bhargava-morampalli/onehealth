# Changelog

What changed in `onehealth`, newest first, in terms of what you see as a user. It starts with the work after the October 2026 code review and the first version of `onehealth report`. Each entry links to its pull request, where the full reasoning and checks are.

## 16 November 2026: lint and type checks in CI, and a monthly refresh of the animal names ([#19](https://github.com/bhargava-morampalli/onehealth/pull/19))

No result changes: the lexicon is still version 8, and the table of animal names holds the same names.

### Lint and type checks

Nothing checked the code for mistakes the tests would not reach, such as an undefined name in a branch no test runs.

- **ruff** checks for undefined names, unused imports and variables, and likely bugs (bugbear), with the rules in `pyproject.toml`. Style is not checked. It found three inner functions that read a loop's variable, which were safe where they are but break when code is moved; they now take the value as an argument.
- **mypy** checks the types in `src/onehealth`. It found about twenty places where the annotations did not say what the code does (rows that hold numbers and text, an array that started as an empty list, an `Optional` that a branch already rules out). None was a bug at run time; each is now annotated or written so that the type is plain.
- **CI** runs both in a job of their own, on Python 3.12: the numpy that `uv.lock` resolves, and its type stubs, differ between Python versions. Both pass on 3.10 to 3.13 as well.
- `ruff`, `mypy` and `types-PyYAML` are in the `dev` extra, pinned by `uv.lock`.

### Keeping the animal names current

The table of animal names was a snapshot of NCBI Taxonomy from 14 November 2026, with no way to see what a newer one would change.

- **The table records the date of its NCBI dump**, and `classify` prints it after the lexicon version: `lexicon v8 (seed v8), animal names from NCBI Taxonomy of 2026-11-14`. The table was rebuilt from the same dump to add the date: its names are identical, and the same dump now always gives the same file, byte for byte.
- **`tools/refresh_animal_names.py`** rebuilds the table from a newer dump and reports what changes: the names added and removed, and every record of the three benchmark datasets whose class changes, with whether the change agrees with NCBI's label. If no name changed, it changes nothing.
- **A monthly workflow** (`.github/workflows/taxonomy.yml`) runs it, runs the tests with the new table, and opens a pull request with the report. A new table changes classes only through a reviewed pull request and a new version. Opening the pull request needs a repository setting, described in the workflow and the lexicon guide.
- Tried on NCBI Taxonomy of 16 November 2026: 24 new scientific names, and no class change in 170,816 benchmark records.

### Also

- `REVIEW.md`, `ROADMAP.md` and the working notes in `docs/` are removed. They described the tool as of October 2026 (lexicon v2); their findings were fixed, and what changed since is in this changelog.
- The package metadata names its author without an email address.
- **Known limitation, documented:** an animal named alone in `isolation_source` with no host (`gull`, `fox`) goes to review, not to Animal. The `classify` page says why and how to settle such rows with a decision.

## 16 November 2026: install without a clone; the linkage test is about a hundred times faster on close relatives ([#18](https://github.com/bhargava-morampalli/onehealth/pull/18))

No result changes: every number in a report is what it was.

### Install without a clone

The only install the documentation gave was an editable one from a clone of the repository.

- **`pip install "git+https://github.com/bhargava-morampalli/onehealth.git"`** installs the newest version; add `@` and a tag (`@v0.2.0`) for a fixed one, to name in a study's methods. `uv tool install` and `pipx install` take the same address. The Install page and the README now lead with this; the clone is for changing the code and running the tests.
- **The package is checked as a user gets it.** A new job in the test workflow builds the wheel, installs it in an empty environment, and runs `classify` and `report` from another directory, so a data file left out of the package fails the build.

### The linkage test on thousands of close relatives

With 4,000 isolates that are mostly close relatives, there are millions of pairs within the SNP threshold, and every one of them was read again for each of the 9,999 permutations.

| 4,000 isolates, 3 million pairs within 100 SNPs | Before | Now |
|---|---|---|
| 9,999 permutations | 9 minutes or more | 3 to 6 seconds |
| Everything else in the linkage sections | 19 seconds | 5 seconds |

- **Large clusters are counted with matrix products**, a batch of permutations at a time, where many of a cluster's pairs are close. Other pairs are still read one by one.
- **Clusters at every threshold of the sweep come from one pass** over the pairs, in order of distance; there were ten passes.
- **The tables are built with sorts**, where there was one search per pair.
- **The results are the same**, to the last digit, for the same seed. Events, expected ranges, p-values, clusters and tables were compared between the old and new code on 48 random datasets (with unclassified isolates, strata, missing distances and shuffled rows) and on the 4,000-isolate one, and are identical. In that last run, on a busy machine, the old code took 23.6 minutes and the new 10.5 seconds. The exhaustive check of the linkage test against every possible relabelling now runs for both ways of counting.

## 16 November 2026: a `resolve` run can be finished, and a failed lookup no longer looks like a record with no source ([#17](https://github.com/bhargava-morampalli/onehealth/pull/17))

No class changes: the lexicon is still version 8.

### A lookup that did not happen is not "no source"

When NCBI gave an error for an ID, or a run was stopped, `resolve` sent the ID to review with an empty `isolation_source` and `host`. `classify` then gave the row the reason `missing_source`, the same as a BioSample whose submitter gave neither. A study run during an NCBI outage would count those isolates as records without a source, and nothing said otherwise.

- **`classify` keeps the two apart.** A row `resolve` could not get a record for is in review as `unresolved:<reason>`: `unresolved:ncbi_error`, `unresolved:interrupted`, `unresolved:biosample_not_found`. `missing_source` now means only what it says: the record was fetched and gives no source and no host.
- **`classify` warns**, with the count by reason, and says how many were never looked up and can be tried again.
- **`resolve` warns** at the end of a run when NCBI gave no usable answer for some IDs. Before, that was visible only in the review file.
- **`report`** says under Data quality how many isolates without a class were never resolved.
- A row with a source is classified as before, whatever `resolve` said, and an override still wins.

### `resolve --resume` finishes an earlier run

A rerun looked every ID up again. For a list of tens of thousands, or PDTs from a 900 MB table, an interruption near the end cost the whole run.

- **`--resume`** reads the earlier output, keeps every ID that run settled, and looks up only the rest: IDs that were `interrupted`, got an `ncbi_error`, had no network, or are new. An ID NCBI answered about (no such record, several BioSamples) is not asked again. The result is the table one uninterrupted run would give.
- **A PDT that was not found** is looked up again only when the command now names `--pd-organism`. When every PDT is settled, the Pathogen Detection table is not downloaded again.
- **A killed run keeps its work.** While it looks IDs up, `resolve` saves the table every 30 seconds, with the IDs still to do marked `interrupted`. Before, only Ctrl-C saved anything; a power cut or a job that ran out of time left nothing.
- **Every table is written whole or not at all**: to a temporary file that then takes the place of the old one. This holds for every command. A run that dies while writing leaves the earlier file as it was.
- Ctrl-C and the NCBI-error warning now say to run the same command with `--resume`.

## 16 November 2026: 0.3.0, a meat species with no host is food; animals named in the text (lexicon v8) ([#16](https://github.com/bhargava-morampalli/onehealth/pull/16))

**This version changes what some classes mean again**, in the same direction as 0.2.0: where NCBI's records show a convention, `onehealth` follows it. Results made with 0.2.x (lexicon v7) are not comparable class for class. The Definitions page has the figures.

On 26,604 NCBI records drawn after v8 was finished and read once, with `source_type` hidden, v8 classifies 97.8% and is right for 96.7% of those (v7: 97.1% and 95.4%). Animal precision goes from 90.7% to 97.4% and Food recall from 91.7% to 97.4%.

### 1. A meat species with no host is Food, anywhere in the text

0.2.0 made the whole values `chicken` and `turkey` Food when there is no host. NCBI's records show the same for the word inside a longer text, and for the other species sold as meat. Of 3,190 benchmark records with no host that the rule now decides, 3,079 (96.5%) are food.

- **`fish`, `poultry`, `salmon`, `beef`, `lamb`, `chicken liver`, `packaged chicken`, `marine fish`** with no host are Food. They were Animal, or in review (`beef`).
- **With a host, the host decides**, as before: `chicken liver` with host *Gallus gallus* is Animal.
- **What says what was sampled still wins**: a specimen (`chicken feces`), a gut site (`chicken cecal content`, `chicken gut`), `live chicken`, `beef cattle`, a farm or litter.
- **`duck`, `pig`, `cow`, `cattle`, `swine`, `bovine` still mean the animal.** A bare `cattle` with no host is an animal in all 46 NCBI records.

### 2. An animal named in a specimen's own text

A specimen with no host is Human by default. 0.2.1 measured that default and found one kind of record it gets wrong that a text can settle: `stool (Odocoileus hemionus)`, `pigeon stool`, `giant panda feces`. The word lists do not know these animals.

- **Any name NCBI Taxonomy gives to an animal now counts**, with the reason `specimens+animal-name-in-specimen`. 79 of the 81 benchmark records placed this way are animal in NCBI; 0.2.1 called them Human.
- **Ordinary words are animals' names too**, so a single word counts only when it names a vertebrate and is not an everyday word: `pigeon` and `deer` count; `tip`, `white` and `sponge` (a moth, a butterfly, an invertebrate) do not, nor `sole`, `ray` or `monitor`. A name of two or more words always counts. In 5,611 patient specimens of a separate NCBI sample, none is mistaken for an animal's.
- **Shellfish tissue is not a patient.** `oyster homogenized tissue dilution` was Human, because `tissue` is a specimen. `oyster` is now a food, so it goes to review as `specimen_or_food`.

### 3. More hosts are recognised

- **A name at the end of a host value**: `European rook`, `free-range white-tailed deer`, `forest musk deer`, and `white stork nestling` (a word for the animal's age may follow). A name in the middle is still not used, nor one in a list (`bologna, hot dogs, sausage-like products`).
- **Oysters, shrimps, crabs and flies by their species name**: `Magallana gigas`, `Litopenaeus vannamei`, `Musca domestica`. The name table now holds the species names of molluscs, of crustaceans and of any other animal NCBI gives a common name, and synonyms. It grew from 0.7 to 1.4 MB.
- **`pets`** as a host is Animal.

### 4. Other changes in lexicon v8

- **A carcass with an animal host is the animal.** `carcass swab` with host porcine went to review as a conflict. In NCBI's records it is animal in 211 of 211. With no host a carcass is still Food.
- **An animal's gut is the animal**: `cecal content`, `cecum`, `intestine`, `gut`, `dropping`. With a human host the host decides (`intestines of patients`).
- **A sponge, a wipe, dust or a bare `swab`** is Environmental (`sponge`: 5,532 NCBI records of 5,669). A swab with a host is the host's: `swab` with host *Homo sapiens* is Human.
- **Foods that were in review**: `oyster`, `clam`, `mussel`, `prawn`, `crab`, `tuna`, `grain`, `cereal`, `oatmeal`, `rice`, `potato`, `seed`, `almond`, `peach`, `guacamole`, `queso fresco`, `hot dog`, `pig ears` and more, and anything `frozen`, `cooked` or `dried`. Each is food in at least 95% of at least twenty NCBI records.
- **`hospital` as the whole value is Human**: the place a patient isolate came from, in 328 NCBI records of 328. So is `healthcare facility`. `hospital sink` is still a fixture. One record of the IMP benchmark is described only as `hospital`, and its curators place it in the hospital environment; `onehealth` now disagrees with them on that row.
- **Body fluids and the placenta** are specimens: `cerebral spinal fluid`, `synovial fluid`, `placenta`.
- One rule key is new for a lexicon of your own: `within`.

### How it was measured

The benchmark grew to 166,158 records from 26 organism groups, in five parts that share no record. v8 was written with the first four in view. The fifth, the fresh part, was drawn when v8 was finished: the unread parts of the Salmonella, *E. coli*, Campylobacter and Klebsiella tables, and eight organism groups chosen by name before any of their records were seen. The rebuild reproduced the four earlier parts byte for byte.

| Fresh part, 26,604 records | v7 | v8, as first read |
|---|---|---|
| Coverage | 97.1% | 97.8% |
| Accuracy | 95.4% | 96.7% |
| Animal precision / recall | 90.7% / 91.2% | 97.4% / 92.0% |
| Food precision / recall | 98.7% / 91.7% | 98.2% / 97.4% |
| Human precision / recall | 98.8% / 95.8% | 98.8% / 95.9% |
| Environmental precision / recall | 81.4% / 84.2% | 81.4% / 85.1% |

**One correction was made after the reading**, and it is the reason for the low Environmental precision: 296 records described only as `healthcare facility`, all human to their submitter. 0.2.0 had added that phrase as a hospital fixture with no labelled record to go by. It is now read like `hospital`. With it, the lexicon that ships scores 97.9% accuracy on the fresh part, but that part has now been seen, so the table above is the fair test. Nothing else was changed after the reading.

What v8 still gets wrong there: 119 animal records called Human (faeces with no host and no word for the animal), and 186 animal records sent to review, most of them an animal named in `isolation_source` with no host (`Larus cachinnans`, `gull`, `fox`).

On the four earlier parts, which v8 was written with, coverage goes from 90.3–96.7% to 97.1–98.8% and accuracy from 96.1–97.5% to 97.9–98.5%. On the IMP benchmark, coverage is 84.9% and accuracy 99.8% as before; one row changes (`hospital`, above). No row of the NCBI *Klebsiella*, demo or example data changes class.

### Also

- The package version is 0.3.0, and `lexicon_version` is 8.

## 16 November 2026: 0.2.1, the statistics checked against R and scipy ([#15](https://github.com/bhargava-morampalli/onehealth/pull/15))

No class changes: the lexicon is still version 7. One number in the report's tables changes in its fifth digit.

### The statistics are checked against independent implementations

`report` computes its statistics with its own code. Each one is now compared with an independent implementation every time the tests run, and they agree. The new page "The statistics, and how each was checked" (`documentation/guides/statistics.md`) shows the comparison.

- **R and vegan.** `tests/data/stats_reference/` holds 88 inputs, most of them counts from the real NCBI *Klebsiella* data, and the 567 values that R 4.3.3 and vegan 2.6-4 give for them. Simpson's index (`vegan::simpson.unb`), the Wilson interval (`prop.test`) and Benjamini–Hochberg q-values (`p.adjust`) agree to 12 decimal places, rarefaction (`vegan::rarefy`) to 9, and the chi-squared statistic (`chisq.test`) to 10 significant digits. The variance behind the interval of Simpson's index is within 0.2% of the variance over 2,000,000 samples drawn in R.
- **Permutation p-values against exact ones.** Lineage enrichment is compared with the hypergeometric distribution for 107 lineage and source pairs, and the feature test with R's p-value from 1,000,000 random tables for 47 genes. On a dataset of ten isolates, every one of the 4,200 possible relabellings is listed (900 within strata), which gives the exact p-value for enrichment, the feature test and sharing events, with and without `--stratify-by`.
- **scipy.** The tables made from the real *Klebsiella* data are recomputed with scipy, and the single-linkage clusters at 0, 5, 21 and 50 SNPs are the connected components scipy finds.
- **The tests were tested.** Eight errors were put into the statistics code one at a time (a one-sided test for a two-sided one, strata ignored, same-host pairs counted, the large-sample variance, and others). Each made the new tests fail.
- scipy is needed to run the tests, not to use `onehealth`.

### One correction

- **The Wilson interval uses the exact 97.5% point of the normal distribution**, 1.959964, as R and scipy do. It used 1.96. Interval bounds in `06_feature_prevalence.tsv` change in the fifth significant digit: for 51 of 150 isolates, 26.9032–41.8959% became 26.9033–41.8957%. Nothing printed in the HTML report changes. Nothing else needed correcting.

### A specimen with no host: the Human default was measured, and kept

`blood` or `stool` with no host named is classified Human. That was `onehealth`'s own judgement. It was measured against NCBI's `source_type`, like the three conventions of 0.2.0: of 15,235 benchmark records placed this way, NCBI's submitters call 14,735 (96.7%) human. The default stands, and the Definitions page now gives it as a fourth convention with these figures.

The measurement also showed two kinds of record that are wrongly Human and that a text could settle. They are documented, and are to be corrected in the next lexicon version:

- an animal named in a way the word list does not know, in a specimen's own text: `stool (Odocoileus hemionus)`, `pigeon stool` (78 benchmark records);
- tissue of shellfish and fish from monitoring programmes: `oyster homogenized tissue dilution` (120 records).

## 16 November 2026: 0.2.0, the classes follow NCBI (lexicon v7) ([#14](https://github.com/bhargava-morampalli/onehealth/pull/14))

**This version changes what some classes mean.** Where NCBI's records show a convention, `onehealth` now follows it, in three places. Results made with 0.1.0 (lexicon v6 and earlier) use the old definitions and are not comparable class for class; `v0.1.0` is the tag to check out to reproduce them. The new Definitions page in the documentation (`documentation/guides/definitions.md`) gives the figures behind each convention and how to choose differently for a study.

### 1. A record's own `source_type` is its class

NCBI's `source_type` is its controlled One Health category (`human`, `animal`, `food`, `environmental`), chosen by the submitter. Until now it was one more text field: a word in `isolation_source` could outrank it, and a host that pointed elsewhere sent the row to review. Now, where it names a class, that is the row's class.

- `classify_reason` is `source_type` where the free text could not settle it (`animal-cattle-steer`, host bovine, `source_type` food: review before, Food now), and `source_type_over_text` where the text said something else (`raw milk` with `source_type` environmental: Food before, Environmental now). `lexicon_suggestion` keeps what the lexicon said.
- In a sample of 159,045 NCBI records, 51,288 have a `source_type` that names a class; 456 of them are now decided as `source_type` and 141 as `source_type_over_text`. None goes to review.
- `resolve` writes `source_type` in a new column of its own. `sample_type` held it mixed with the tissue and NCBI's material type, so `classify` could not tell which it was.
- Overrides and decisions still win. Any other value (`other`, `not provided`) declares nothing.

### 2. Farm settings are Environmental, and a farm with a host is the host

Until now `manure`, `litter`, `barn` and `farm` were Animal. In 139,554 NCBI records with a `source_type`:

- **Manure, litter and bedding are Environmental**, with or without a host. NCBI submitters record manure as environmental in 646 records of 747, litter in 77 of 77.
- **A farm, a barn, a coop or a poultry shed with no host is Environmental** (a farm: 141 records of 187). **With a host, the host decides**: `pig farm` with host *Sus scrofa* is Animal (507 of 641 farm records with a host are animal), and `farm` with host *Homo sapiens* is Human. The reason is then `farm-environment+<host rule>`.
- Swabs of a poultry house (`drag swab`, `boot swab`) are Environmental.
- `farm worker stool` is no longer Animal: `farm` is not an animal word any more.

### 3. `chicken` or `turkey` on its own, with no host, is Food

They were Animal. NCBI asks for a host when a sample is from an animal, and a bare `chicken` with no host is food in 810 records of 867, a bare `turkey` in 255 of 256. The reason is `meat-species-no-host`. With a host the host decides as before, and `duck`, `fish`, `pig` and `live chicken` are still the animal. A veterinary study whose `chicken` is a bird says so with one line in a `--decisions` file.

### Other changes in lexicon v7

- **A floor drain is not assumed to be in a hospital.** `floor drain` had the healthcare flag; in NCBI records it is a food plant's (`floor drain in produce prep room`). It is Environmental without the flag. `hospital floor drain` and `sink drain` keep it.
- **A specimen word in a food is not a patient.** `pork blood sausage` was Human, because `blood` is a specimen and a specimen with no host defaults to Human. With no host and a food named in the same text, the row now goes to review as `specimen_or_food`.
- **Fixtures and equipment of food premises** are Environmental even when the text names the food: `drain-raw meat`, `deli case`, `cutting boards`, `utensils`, `cooler`, `floor`.
- **Water** in any phrase is Environmental (`pond water`, `tap water`), with `creek`, `stream`, `runoff`, `lagoon` and `pond`. `bottled water` is Food, and `water buffalo` is Animal.
- **More body sites and fluids** are specimens: `respiratory`, `rectal`, `perirectal`, `urinary`, `vaginal`, `throat`, `aspirate`, `secretion`, `ascites`, `body fluid`. Organs are not, since in NCBI records a `lung` or `intestine` is an animal necropsy sample.
- **Hospital equipment** is healthcare environment: `bed rails`, `ECG monitors`, `nebulisers`, `syringe drivers`, `infusion stands`, `stethoscopes`, `hospital bed`, `patient room`. Not `bronchoscope` or `duodenoscope`: in the IMP benchmark those are patient samples taken through the scope.
- **`human ...` as a whole value** is Human: `Human Listeriosis`, `human infection`, `Human, Other`.
- **Accents no longer split words.** `pâté` was read as the letters `p` and `t`; it is the word `pate`.
- Three rule keys are new for a lexicon of your own: `no_host`, `setting_of` and `yields_to_host`.

### How it was measured

The benchmark of NCBI records with a `source_type` grew from 73,094 records to 139,554, from 18 organism groups, in four parts that share no record. With `source_type` hidden:

| Final part, 34,252 records | v6 | v7 |
|---|---|---|
| Coverage | 94.6% | 95.5% |
| Accuracy | 95.3% | 96.5% |
| Animal precision | 88.2% | 91.7% |
| Environmental recall | 78.7% | 90.2% |
| Food recall | 81.3% | 83.3% |
| Animal recall | 97.5% | 97.0% |

The final part was read once, after v7 was finished, and v7 was not changed afterwards. The road there was less tidy and is recorded in the benchmark's README: v7 was checked against the held-out part, where three of its additions proved too broad and were narrowed (`canal` also matched "alimentary canal"); a validation part meant to be read once then showed that a farm with an animal host is the animal in NCBI's records, not a conflict, and that rule was corrected; so a fourth part was drawn for the final figures.

On the IMP benchmark, coverage is 84.9% and accuracy 99.8% (v6: 84.8% and 99.8%). Two rows change: `rectal infection` is now Human, and `manure` with a pig host is Environmental, which is what the curators say. Animal precision goes from 89.5% to 91.1%. No row of the NCBI *Klebsiella*, demo or example data changes class.

### Also

- The package version is 0.2.0, and `lexicon_version` is 7.
- The tag `v0.1.0` marks the last version with the old definitions.

## 15 November 2026: any animal host, the review queue by source text, and a licence ([#13](https://github.com/bhargava-morampalli/onehealth/pull/13))

Three things from the list the last review left open, and one it proposed that turned out not to be worth building.

### Any animal is a host (lexicon v6)

The lexicon named the hosts seen most often (`chicken`, `Bos taurus`, `dog`), so a gull, a macaque or a water buffalo went to review as a host it could not place. v6 recognises any animal host by the names NCBI Taxonomy gives it:

- the scientific name of a chordate at any rank (`Larus dominicanus`, `Macaca`, `Cervidae`), a species with a subspecies or breed (`Bos primigenius taurus`), or a genus with `sp.`;
- the common name of any animal, singular or plural (`gull`, `sika deer`, `mud crab`).

Only a whole name counts, never a word inside a longer value, so `sea water` is not the fish genus *Sea* and `white stork nestling` stays in review. Scientific names are limited to chordates because other animals have genera called *Sepsis*, *Cancer* and *Colon*. Plants and bacteria are not placed. The names (82,711 scientific, 35,780 common) ship with the package, so `classify` still needs no network; `tools/build_animal_names.py` rebuilds them from the NCBI Taxonomy dump.

It is a rule in the lexicon, not code around it: a rule can say `taxonomy: animal` or `taxonomy: human`. The two new rules, `host-taxonomy` and `host-taxonomy-human`, come after the word-list rules of the same priority, so a host those know keeps its reason.

| | v5 | v6 |
|---|---|---|
| IMP benchmark: animal isolates found | 38 of 51 | 51 of 51 |
| IMP benchmark: coverage / accuracy | 84.5% / 99.8% | 84.8% / 99.8% |
| NCBI `source_type`, held out: Animal recall | 88.3% | 91.7% |
| NCBI `source_type`, held out: coverage / accuracy | 88.7% / 93.9% | 89.1% / 93.9% |

The 13 IMP isolates are 11 American crows and 2 silver gulls. In a sample of 159,045 NCBI records, 495 more rows are classified, and 43 that v5 called Environmental or Food from a word in the source now go to review, because the host is recognised as an animal and disagrees. No row of the NCBI *Klebsiella*, demo or example data changes.

### `human` in a specimen's own text

`pig feces` was Animal because the text names the animal, but `human feces` was Human only by the no-host default, and since the last release it went to review in a record filed as environmental, food or other. The new rule `human-word` lets the word decide, with the reason `specimens+human-word`. It is `specimen_only`, a new rule key: on its own the word decides nothing, so `human contact surface` is still a surface. In the NCBI sample, 340 of the 583 rows called Human by the no-host default now have this reason, so `+no_host` counts only the calls that really are a default, and 9 rows come back from review.

### The review queue by source text, and decisions

`classify` now also writes `<out>.review_groups.tsv`: each distinct source text in the review queue once, with the number of rows that carry it, largest first. Texts are compared by their words, so `Bed Rails`, `bed_rails` and `bed-rails` are one. In the NCBI sample the 3,794 review rows that say anything hold 823 distinct texts, and the 100 most common cover 69% of them.

The table has empty `l1`, `healthcare_environment` and `note` columns. Fill in `l1` for the texts you can place and pass the table back with the new `--decisions` option: each decision applies to every row with that text.

- A decision wins over the lexicon, for rows it classifies as well as rows in review, so a food-safety study can decide that a bare `chicken` is Food. An override for one isolate still wins over a decision for its text.
- A mistake stops the run: a class that is not a class, two decisions for one text, or a decision with no text.
- Decided rows have `classify_reason` `decision` and `override_applied` `yes`.
- `--review-groups` chooses where the table is written.

### Licence

The code is now under the MIT licence (`LICENSE`), and the package metadata says so.

### Proposed, measured, and not built: reading `IFSAC_category` and `ontological_term`

The last review suggested reading these two NCBI fields as evidence, expecting them to settle records with no `source_type`. They do not occur there: of 107,757 sampled records without a usable `source_type`, 158 have an IFSAC category and 152 an ontological term. Submitters who fill them in fill in `source_type` too, and `classify` already reads that. Where `source_type` is present, `classify` agrees with it for 98.0% of records. So nothing was added.

### Also

- `lexicon_version` is 6. A lexicon of your own without the new rules behaves as before.
- The README no longer says the documentation site will be published once the repository is public: it is.

## 14 November 2026: lexicon v5, one lexicon, and fixes from a review of the whole tool ([#12](https://github.com/bhargava-morampalli/onehealth/pull/12))

A review ran `classify` on 159,045 NCBI Pathogen Detection records from six organism groups and compared it with NCBI's own `source_type`, read every source file, and tried the commands on awkward input. Everything it found that has one right answer is fixed here, each with a test. What it found that is a decision is listed at the end.

### One lexicon

`onehealth` now ships one lexicon, `src/onehealth/data/lexicon.yaml`, in place of a file per version. v2, v3 and v4 shipped next to the current one so that `--lexicon v3` could reproduce earlier results; each one added tests, benchmark columns and lookup code, and nobody outside the project had results to reproduce. They are in the repository history. `--lexicon v4` now stops with an error that says so, and `--lexicon` takes a YAML of your own, which can be an old version's file. Every table still records the `lexicon_version` that made it.

### Lexicon v5: meat is food, a swab of the premises is not

The lexicon is now version 5, and classified tables record `lexicon_version` 5.

- **Retail meat was called Animal.** `raw intact chicken`, `ground turkey`, `chicken breast`, `comminuted chicken` and `smoked salmon` matched the animal word and no food rule. They are Food in v5. A species on its own (`chicken`) is still the animal, and so is an organ (`chicken liver`), which can be a necropsy sample; `chicken liver (giblets)` is Food. A cut of meat with an animal host goes to review as `host_source_conflict`, as `chicken meat` with a chicken host always did.
- **Foods that no rule knew** are Food: `pork`, `sausage`, `ham`, `shrimp`, `peanut butter`, `ice cream`, `avocado`, `enoki mushroom` and others. They are the foods NCBI records name most often. Two-word foods used to go to review as `taxon_only_source`, taken for species names.
- **Swabs of food premises were called Food**, because their text contains the word: `environmental swab sponge non food contact surface`, `food processing environment`. So were `poultry processing plant` and `salmon processing plant` called Animal. They are Environmental in v5.
- **`food handler` was Food.** It is Human. **`food animal` was Food.** It is Animal.
- **`fecal`, `faecal` and `stools`** are specimens, like `feces`, `faeces` and `stool`.

Effect on NCBI records whose submitters filled in `source_type`, with that field hidden from the classifier (the new benchmark in `tests/data/ncbi_pd_source_types/`). The held-out records were not looked at until v5 was finished, and v5 was not changed afterwards:

| Held out, 39,722 records | v4 | v5 |
|---|---|---|
| Coverage | 76.9% | 88.7% |
| Accuracy | 80.9% | 93.9% |
| Animal precision | 45.4% | 73.7% |
| Food recall | 30.0% | 70.6% |
| Environmental recall | 43.8% | 77.6% |
| Human precision | 98.0% | 98.0% |

On the IMP benchmark v5 changes one row, a pork product that no v4 rule matched: 84.5% coverage and 99.8% accuracy, as v4. No row of the NCBI *Klebsiella*, demo or example data changes class between v4 and v5. The v4 figures were measured before v4 was removed.

Animal is still the weakest class: a bare `chicken` or `turkey` is read as the animal and is often meat. The Benchmarks page (`documentation/benchmarks.md`) has the confusion matrix and what the label does and does not mean.

### `classify`

- **A specimen with no host is no longer called Human when its record was filed as not clinical.** `feces` or `wound` with no host defaulted to Human whatever else the record said. If the record carries the BioSample package `Pathogen.env` ("Pathogen: environmental/food/other"), or `epi_type` `environmental/other` in a Pathogen Detection table, the row now goes to review as `no_host_env_package`. In the NCBI sample, 199 of the 586 host-less specimens with a known package were filed this way. Any other package, or no package column, keeps the default and its `+no_host` reason. It is a change in `classify`, not in the lexicon; it moves one row of the NCBI *Klebsiella* data (`feces`, Brazil) from Human to review.
- **`lexicon_version` tells an edited lexicon from the shipped one.** A copy of the shipped lexicon with a rule changed was recorded as the shipped version, and a file with no `version` was too. An edited copy that keeps the shipped version is now `5-modified-<digest>`, and a file with no version is `custom-<digest>` (the first 8 hex digits of its SHA-256). A file with its own version name is recorded under that name, as before.

### `resolve`

- **New column `package`**, the BioSample package of each record (`Pathogen.cl.1.0`, `Pathogen.env.1.0`, ...), after `sample_type`. `classify` uses it for the rule above. A `--metadata` table's `attribute_package` or `epi_type` fills it too.
- **Ctrl-C no longer loses the run.** Output was written only at the end, so stopping a long lookup wrote nothing. What NCBI has already answered is now written, the IDs not reached are marked `interrupted`, and `resolve` exits with status 130.
- **A reply of the wrong shape sends one ID to review.** A reply that was valid JSON but not shaped as NCBI documents it (text where a list of attributes should be) stopped the whole run with a traceback. It is now `ncbi_error:unexpected reply ...` for that ID.
- **A BioSample read from E-utilities is the one asked for.** The search covered all fields and took the first hit, which can be another sample that mentions the accession. It now searches the accession field first, and uses a record only if it carries the accession.
- **`NCBI_API_KEY`** is read when `--api-key` is not given, which keeps the key out of the shell history. Requests name `onehealth` and its version as the client.

### Input files

- **NCBI Pathogen Detection metadata is read as downloaded.** Its header starts with `#`, so it was dropped as a comment and the first isolate became the header; `classify` then stopped with an error that showed that isolate's values as column names. The same goes for NCBI's `assembly_summary.txt`. A `#` line is the header when it is the last one above the first row, is tab-separated, and names more columns `onehealth` reads than the line below it; a commented-out row or a commented-out copy of the header is still a comment.
- **A column name that appears twice is an error.** The second column silently replaced the first, so with two `host` columns, blood from a patient was classified by the second host.
- **Output cells that a spreadsheet would run** (starting with `=`, or with `+`, `-` or `@` and containing `(` or `|`) are written with a leading apostrophe. The text comes from public BioSample records. `-`, `-5` and `+ve` are not touched.

### `report`

- **Simpson confidence intervals no longer collapse to a point.** The large-sample variance behind the interval is exactly zero when lineages are equally common, so a source of ten isolates in ten lineages had D = 1.00 with the interval 1.00–1.00, and the report reads intervals that do not overlap as a difference. The interval now uses the exact sampling variance of the index (Simpson 1949): 0.91–1.00 for that source. For large samples the two agree. D itself is unchanged; `ci_low` and `ci_high` in `02_diversity.tsv` are a little wider for small sources.
- **"Not typed" is not a lineage.** `NF` and `ND` (SRST2), `Novel` (ARIBA), `untypeable` and the like were kept as lineage names, and so was the `0` that Kleborate writes in an ST column. Every untyped isolate then belonged to one lineage, which was reported as shared between sources. These isolates now count as untyped, and Data quality says how many. A cluster column keeps its cluster `0`.
- **A distance matrix that is one triangle, or has a negative distance, stops with an error.** Only the upper triangle is used, so a lower-triangle file was read as every pair being 0 SNPs apart.
- **Mash or ANI distances are flagged.** If every distance is a fraction between 0 and 1, every pair is within any SNP threshold. The report still runs, and Data quality says the linkage sections mean nothing for that file.
- **More date formats.** The INSDC forms with a month name (`15-Mar-2019`, `Mar-2019`) are read; they used to be dropped from the counts by year and the year strata without a word. A numeric date that ends in the year (`15/03/2019`) counts as the middle of that year. Data quality reports dates that still cannot be read.

### Tests

- GitHub Actions now runs the whole suite on Python 3.10 to 3.13 for every pull request. Until now only `tests/test_docs.py` ran there.
- `tests/data/ncbi_pd_source_types/` is the new benchmark: 73,094 real NCBI records in a development and a held-out split, with the script that rebuilds it.

### Found, and left for a decision

- **No licence.** The repository has no `LICENSE` file, so nobody else may legally reuse the code. Which licence is the author's choice. A `CITATION.cff` would tell users how to cite it.
- **A species on its own.** `chicken` or `turkey` with no host is Animal, and in NCBI food-safety records it is usually meat. Calling it Food would be wrong for veterinary data.
- **Farm settings.** Manure, barns and poultry houses are Animal here and `environmental` in NCBI's `source_type`.
- **Hosts are a word list.** A gull, macaque or water buffalo host goes to review as unrecognised. Classifying hosts by NCBI taxonomy would fix that for every species at once.
- **A specimen with no host and no package** is still Human by default (`+no_host`).
- **`resolve` cannot resume.** An interrupted run keeps its answers but a rerun looks every ID up again.

## 7 November 2026: fixes from a bug review ([#11](https://github.com/bhargava-morampalli/onehealth/pull/11))

A review found the cases below; each fix has a test. None of them changes a single row of the IMP, NCBI Pathogen Detection, demo or example data under any lexicon version, or any report table made from them. They change how `classify` works rather than the lexicon, though, so rows like the ones described can classify differently under v2, v3 and v4 alike.

- **An unknown or species-name `isolation_source` no longer hides other fields.** A rule reads its fields in order and stops at the first match. When `isolation_source` was `Bos taurus` or `unknown tissue`, that match was then thrown away, and the same rule never looked at `source_type`. `Bos taurus` with `source_type` animal is now Animal instead of going to review.
- **`unknown animal` is Animal.** In `isolation_source` it was read as a missing value and sent to review. The species is unknown, but the class is not; the same goes for `unknown food`. `unknown` with only general words (`unknown source`, `unknown environment`) still means nothing, so unknown never becomes Environmental.
- **A specimen that another field contradicts goes to review.** `urine` with host *Homo sapiens* and `source_type` food was called Human; it is now `host_source_conflict`. Words in the specimen's own text do not count: `stool sample from food handler` is still Human.
- **Overrides and truth labels can be keyed by strain** in a table with no BioSample, assembly or ID column. An overrides file can use `strain` or `isolate`; a truth file can also use `sample`, `sample_name` or `name`. `classify` now warns about override IDs that match no row, and `validate` says how many truth IDs matched no row. An override that matched nothing used to be ignored without a word.
- **`validate` counts an isolate once** when the truth file names two of its assemblies, and stops if their labels differ.
- **`report --dereplicate` no longer merges unclassified isolates.** It had inflated the number of collapsed clusters given in Methods.
- **SNP matrices:** a square matrix with rows missing, or in a different order from its columns, now stops with an error. A truncated one used to be read as a long table with a made-up isolate.
- **Lineages:** an isolate given two different lineages under two of its IDs stops with an error; before, one of them silently won. Spaces around column names in a CSV lineage table are ignored.
- The README points to lexicon v4, and the package builds without a warning on `uv` 0.10 and later.

## 5 November 2026: lexicon v4, fixtures without a setting go to review ([#10](https://github.com/bhargava-morampalli/onehealth/pull/10))

The seed lexicon is now v4, and classified tables record `lexicon_version` 4. v3 still ships, and `classify --lexicon v3` reproduces earlier results exactly.

A fixture named without its setting, such as `Rm 14 Drain`, `Sink swab` or `Shower 3`, says neither whose sample it is nor where the fixture is. Without a host it already went to review. With host *Homo sapiens*, which submitters often record for hospital fixture swabs, v3 called it Human. v4 sends it to review, with the reason `fixture-no-setting`.

- A setting word, a specimen or an animal host still decides. `sink drain` and `hospital sink` are healthcare environment, and a bare drain with a pig host is Animal, as other farm settings are.
- `drain tip`, the tip of a surgical drain, is now a specimen, so the host decides, as for `catheter tip`.
- A rule that sends rows to review no longer counts as a class when `classify` checks the host against the source. Without this, a bare drain with an animal host would have gone to review as a conflict. No v2 or v3 result changes: their only such rule, `animal-feed`, has the highest priority.

Effect on the IMP benchmark:

| | v3 | v4 |
|---|---|---|
| Coverage | 84.8% | 84.5% |
| Accuracy | 99.5% | 99.8% |
| Environmental isolates called Human | 12 | 3 |
| Review queue | 661 | 681 |

The 9 room drains and swabs that the curators call environmental, and the 9 drain samples they could not classify, now go to review instead of Human. Two human samples go with them (`Drain Swab`, `drain inserting site`). The NCBI Pathogen Detection and demo data are unchanged. As with v3, the change was found by looking at the benchmark, so this is not an independent test.

## 4 November 2026: input files that were misread without a word ([#9](https://github.com/bhargava-morampalli/onehealth/pull/9))

A comma-separated file was read as a table with one column. `classify` then found no host or isolation source and sent every row to review as `missing_source`, and `resolve --metadata` matched none of its rows, so every ID was looked up at NCBI instead (or sent to review with `--no-network`). Both commands finished as if nothing were wrong.

- **A comma- or semicolon-separated table** now stops with an error that says so, in every command that reads a table: metadata, overrides, truth labels, features. Semicolons are what Excel writes in many European locales. A comma inside a tab-separated table is still just text.
- **An ID list with a comma or semicolon** stops with an error naming the line. No accession contains one, so the line is a CSV row or several IDs on one line; before, it was looked up as one ID and came back `unknown_id`.
- **`classify` on a table with none of the columns it reads** stops with an error. The usual cause is an ID list given to `classify` instead of `resolve`, whose first ID was read as the column name and lost. With `--overrides` and an ID column it is a warning, and the overrides still apply.
- **`resolve --metadata` without an ID column** stops with an error, instead of matching nothing.
- The README said Excel's UTF-8 export works. That export is CSV, so it never did; the README and the Input files page of the documentation now say what does.

`report --lineages` reads mlst's default output, which has no header row; before, it stopped with "no lineage column". The README already listed mlst, but only output with a header row worked.

- The isolate is mlst's FILE column, and file names now match isolates in the lineages, features and distances files: `ISO1.fasta`, `assemblies/ISO1.fna.gz` and shovill's `ISO1/contigs.fa` all find `ISO1`. The most specific form that names an isolate decides. If it names two, the name is left unmatched rather than tried in a shorter form.
- An ST of `-` (none assigned) gives no lineage. mlst's `--csv` output works too.
- A name listed twice in a lineage table with two different lineages now stops with an error; before, the last one silently won.

Reports on existing inputs are unchanged: the NCBI Pathogen Detection and demo reports are identical before and after.

## 4 November 2026: documentation site ([#8](https://github.com/bhargava-morampalli/onehealth/pull/8))

A documentation site in `documentation/`, built with MkDocs: install, a tutorial on the example data, a page for each command with every option and output column, guides to input files and the lexicon, the benchmarks, and this changelog. A GitHub Actions workflow builds it on every pull request and fails on a broken link. It publishes the site to GitHub Pages from `main` once the repository is public, and does nothing while it is private. A test fails if a command gains an option its page does not mention.

To build it locally: `pip install -e ".[docs]"`, then `mkdocs serve`.

## 4 November 2026: lexicon v3 ([#7](https://github.com/bhargava-morampalli/onehealth/pull/7))

The seed lexicon is now v3, and classified tables record `lexicon_version` 3. v2 still ships: `classify --lexicon v2` reproduces earlier results exactly. The additions are general terms found in the review queues of the IMP benchmark and the NCBI Pathogen Detection test data. Each is marked `v3` in the file.

- **Hospital environment:** `bedpan`, `hospital waste water` (two words), and plurals of the drain phrases (`sink drains`, `floor drains`, ...). Matching is by whole word, so `drains` never matched `drain`.
- **Other environment:** `waste water`, `wastewater treatment`, `WWTP`, `influent`, `environmental swab(s)`, and `water` as a whole value.
- **Patient samples (the host decides):** any `catheter`, not only `catheter tip`, and diagnoses ending in `tract infection`. A source that is only `human` or `patient` is Human, and is a conflict with an animal host.
- **Laboratory:** `passaging`.

Effect:

| | v2 | v3 |
|---|---|---|
| IMP benchmark coverage | 84.5% | 84.8% |
| IMP benchmark accuracy | 99.2% | 99.5% |
| IMP review queue | 674 | 661 |
| NCBI Pathogen Detection rows classified | 253 of 307 | 263 of 307 |

Every newly classified row agrees with its curated or NCBI label, and no row moves from a correct class to a wrong one. Both datasets were used to find the additions, so this shows they work as intended; it is not a measure of accuracy on new data. As a cross-check, the additions found in one dataset changed nothing wrongly in the other. One change makes coverage look worse: 12 `environmental swab` samples recorded with host *Homo sapiens* were called Human by v2, which the curation says is wrong. v3 sends them to review as a host/source conflict.

Left as they are, on purpose:

- A bare fixture without a setting word (`Shower 3`, `Sink swab`) still matches no rule, like bare `sink` and `faucet`. A wrong hospital flag would quietly move isolates into the report's healthcare-environment group. Counting them as hospital fixtures would lift IMP coverage to 86.0%. Such a row goes to review, unless its host decides: 9 room drains and sink swabs in the IMP benchmark have host *Homo sapiens* and are called Human. (Corrected 4 November 2026: this entry first said they all go to review.)
- Strings tied to one study (`ensuite 7/8`, `arjo bathroom`) and typos are not added.

## 4 November 2026: PDT lookups that work ([#6](https://github.com/bhargava-morampalli/onehealth/pull/6))

`resolve` could find a PDT only when NCBI had assembled the genome itself and named the assembly after it. Most genomes keep their submitter's name (`ASM80747v2`), so 60 of 65 real Klebsiella PDTs came back `pdt_unresolved`. New option: `--pd-organism Klebsiella` (several groups, comma-separated). It reads that Pathogen Detection organism group's current metadata table from NCBI, which lists every PDT with its BioSample and assembly. The table is streamed, not saved: about 200 MB for Klebsiella and 900 MB for Salmonella. Reading stops once every PDT is found.

- On the same 917 real IDs, all 65 PDTs now resolve (5 before), and every BioSample agrees with the Pathogen Detection fixture. All 307 fixture PDTs resolve to the right BioSample and assembly. The run is also faster, 30 s instead of 83 s, because the per-PDT searches are skipped.
- An outdated or unversioned PDT gets its current version. The `pdt` and `assembly` columns are filled in from the table.
- A PDT the table does not list, or one listed without a BioSample, still gets the old assembly-name lookup. If the table cannot be read, `resolve` says so and uses that lookup for every PDT.
- A misspelt group stops the run before any lookup and suggests the closest names, for example "Did you mean Klebsiella?". Without the option nothing changes, except a hint when PDTs stay unresolved.

## 4 November 2026: edge-case sweep ([#5](https://github.com/bhargava-morampalli/onehealth/pull/5))

Every command was run against about 90 awkward inputs: missing, empty or mis-encoded files, Excel exports, broken lexicons, impossible options, tiny and mismatched tables. Problems found and fixed:

**Errors you can read.** A user mistake (missing file, a directory instead of a file, a file that is not UTF-8, a broken lexicon, `-o` pointing at a file, an image name ending in `.txt`) used to end in a Python traceback. It now prints one line, such as `onehealth classify: error: No such file or directory: metadata.tsv`, and exits 1. A real bug still shows its full traceback.

- Lexicon files are checked rule by rule: invalid YAML, a rule without `id` or `l1`, a non-numeric priority, `rules` that is not a list, unquoted `yes`/`no`. A single phrase written as `phrases: blood` is now one phrase. Before, it was read letter by letter and matched nothing.
- `report --threshold` must be a finite number of SNPs, 0 or more. `nan` crashed; negative or infinite values gave meaningless reports. `--seed` must be 0 or more.
- `resolve` retries a read that times out or a connection dropped mid-reply, then sends just that ID to review as `ncbi_error`. Before, it ended the whole run.

**Files from Excel and other tools are read correctly.**

- **Byte-order mark.** Excel's "UTF-8" export starts with one, and it became part of the first column name. As a result, overrides were silently not applied, local metadata did not join, the first ID in a list was "unknown", and `validate` could not find its ID column.
- **A cell starting with a double quote** (`"wound swab`) used to swallow the rows after it: two rows went in, one came out. Now one line is one row. A fully quoted cell, as Excel and pandas write it, is still unwrapped.
- `NaN`, `#N/A`, `not determined` and `undetermined` count as missing values.
- Classifying a table directly, standard INSDC missing-value terms such as `not collected` and `not provided` were reported as `taxon_only_source`. They are now `unknown_source`.
- Local metadata exported from NCBI Datasets (`Assembly Accession`, `Assembly BioSample Accession`) now joins.

**Small, empty and mismatched inputs.**

- A report on a single isolate crashed (empty legend). So did a report with a threshold beyond the 100-SNP histogram axis (the canvas grew billions of pixels wide).
- `summarize`, `plot`, `validate` and `report` stop with a clear message when given a table that was never classified. Before, they reported everything as Unclassified.
- `validate` stops when the truth file gives one isolate two different labels, directly or under two of its IDs. Before, the last label silently won.
- An overrides file with no ID or no `l1` column is an error. Before, it was ignored, or sent every isolate to review.
- When a report section cannot be drawn (empty features file, no variable features, a distance matrix sharing fewer than two isolates with the table, an empty lineage table), the data-quality section says why. Before, the section silently vanished.

**Test data.** The NCBI Pathogen Detection gene fixture (`tests/data/ncbi_pd_klebsiella/features.tsv`) had stray quotes on the first and last gene of every list. The cause was NCBI's quoted list cells, which the build script never unwrapped. 23 genes were counted under two names. The fixture was rebuilt from the same release. The real-data report's gene table goes from 186 to 168 features, and no other result changes.

## 4 November 2026: NCBI lookups in batches ([#4](https://github.com/bhargava-morampalli/onehealth/pull/4))

`resolve` fetches BioSample and assembly reports from NCBI Datasets 100 at a time instead of one at a time. On 917 real IDs: 1,309 requests and 478 s before, 203 requests and 76 s after, with identical output. Anything a batch does not answer cleanly is still looked up on its own, so results cannot change. Genome reports are paged (20 per page by default), so batches follow every page.

## 4 November 2026: `--permutations 0` skips the tests ([#3](https://github.com/bhargava-morampalli/onehealth/pull/3))

`report --permutations 0` crashed in the linkage section, and gave p = 1.0 everywhere else as if a test had run. It now skips the permutation tests for a quick look: p- and q-values are blank, the report says the tests were skipped, and the observed-vs-expected section is omitted. Negative counts are refused.

## 3 November 2026: faster report, much less memory ([#2](https://github.com/bhargava-morampalli/onehealth/pull/2))

- **Gene-by-source test:** about 10 minutes before, 13 s after, at 10,000 isolates × 1,000 genes. It uses float matrix products and many permutations per product, with identical p-values.
- **Reading a 5,000-isolate SNP matrix:** 30 s and 2.7 GB before, 2.5 s and 0.24 GB after.
- **Linkage and dereplication:** only pairs within reach of a threshold are kept. 10,000 isolates now run in 0.93 GB, where the old code would have needed an estimated 11 GB.

Every table and figure is byte-for-byte the same as before.

## 3 November 2026: fixes from a code review ([#1](https://github.com/bhargava-morampalli/onehealth/pull/1))

- **Tables keyed by name.** A table that names isolates by `strain`, `sample` or `Name` gave a report with 0 isolates, and `summarize` counted all rows as one isolate. Isolates are now found by BioSample, ID, assembly, SRA/PDT, strain or sample name.
- **Species-name guess.** The guess that a source is "only a species name" overruled the lexicon, so rules for phrases like `bedpan washer` were silently ignored. A rule for source text now wins.
- **Bad NCBI replies.** One garbled reply ended a whole `resolve` run. Now it sends just that ID to review.
- **API key.** `--api-key` now raises the request rate to NCBI's 10 per second, and a rejected key is kept out of the output table.
- **Diversity headline.** It no longer reads "D = nan".
- **Help text.** Help and README install steps are corrected.

## Known, not yet fixed

- **Setting words apart from the fixture.** `washroom sink in hospital intensive care unit` goes to review, because rules match whole phrases (`hospital sink`), not words spread through a sentence. That needs a change to how rules match, not more lexicon entries.
