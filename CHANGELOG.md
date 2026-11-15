# Changelog

What changed in `onehealth`, newest first, in terms of what you see as a user. It starts with the work after the October 2026 code review and the first version of `onehealth report`. Each entry links to its pull request, where the full reasoning and checks are.

## Unreleased: any animal host, the review queue by source text, and a licence (branch `host-taxonomy-ncbi-terms-grouped-review`)

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
