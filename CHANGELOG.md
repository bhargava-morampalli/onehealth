# Changelog

What changed in `onehealth`, newest first, in terms of what you see as a user. It starts with the work after the October 2026 code review and the first version of `onehealth report`. Each entry links to its pull request, where the full reasoning and checks are.

## Unreleased: edge-case sweep (branch `edge-case-sweep`)

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

**Internal notes are no longer published.** `ROADMAP.md`, `REVIEW.md` and `docs/` are working notes. They are now in `.gitignore` and removed from the repository (they remain in its history). Links to them from the README, lexicon and code are reworded.

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

- **PDT lookups mostly fail.** 60 of 65 real Pathogen Detection PDTs come back `pdt_unresolved`, because `resolve` searches for an assembly named after the PDT, and most genomes keep their submitter's name. The fix is to look PDTs up in Pathogen Detection's own metadata.
- **Lexicon gaps.** Lexicon additions found by the IMP benchmark, such as `Shower 3`, `bedpan washer` and `WC sink drains`, are not made yet. They change results, so they should come with a new lexicon version.
