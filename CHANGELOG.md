# Changelog

What changed in `onehealth`, newest first, in terms of what you see as a user. It starts with the work after the October 2026 code review and the first version of `onehealth report`. Each entry links to its pull request, where the full reasoning and checks are.

## Unreleased: fixes from a bug review (branch `bug-fixes-nov6`)

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
