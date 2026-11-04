# Code and method review: `onehealth` v0.1.0

**Date:** 2026-10-17 (later changes are in [CHANGELOG.md](CHANGELOG.md))
**Scope:** every file under `src/`, `tests/`, `examples/`, the README, and the existing notes in `docs/`.
**Result:** 8 groups of scientific or methodological problems and 15 software bugs were found and fixed. 5 design questions are left for you to decide (section 6). The test suite grew from 73 to 133 tests, and all pass.

---

## Contents

1. [Summary: the five most important findings](#1-summary-the-five-most-important-findings)
2. [How the review was done](#2-how-the-review-was-done)
3. [Scientific and methodological fixes](#3-scientific-and-methodological-fixes)
4. [Software bugs fixed](#4-software-bugs-fixed)
5. [Simplifications](#5-simplifications)
6. [Decisions left to you](#6-decisions-left-to-you-not-changed)
7. [Known limitations still present](#7-known-limitations-still-present)
8. [Verification: before vs after](#8-verification-before-vs-after)
9. [What changed for you as a user](#9-what-changed-for-you-as-a-user)
10. [Files changed](#10-files-changed)
11. [References](#11-references)

---

## 1. Summary: the five most important findings

| # | Finding | Effect on results | Status |
|---|---|---|---|
| 1 | Specimen words (`blood`, `urine`, `faeces`, `tissue`, `swab`) were treated as **Human**. | `pig faeces`, `dog urine`, `chicken tissue` were silently called **Human**. A real animal host plus a specimen word was always sent to review. | Fixed (S1) |
| 2 | Bare `clinical` / `clinical isolate` was mapped to **healthcare environment**. | Human patient isolates were counted as **Environmental**. This inflates the exact Human–Environmental overlap a One Health study looks for. | Fixed (S2) |
| 3 | The local-metadata join was **not reproducible**. | The same input gave different `isolation_source` values on different runs: 5 of 8 hash seeds gave one answer and 3 gave another. | Fixed (B1) |
| 4 | NCBI's own One Health field (`source_type`) was shadowed, and the material type (`whole organism`) counted as source evidence. | Records with `source_type = human` could never be classified. | Fixed (S3) |
| 5 | The summary mixed classified and unclassified rows, and `Unclassified` made clusters look "mixed". | Overlap counts suggested cross-compartment sharing that was really "unknown source". | Fixed (S8) |

---

## 2. How the review was done

1. **Read everything.** I read all source, tests, fixtures and README, plus the two notes in `docs/` and the git history.
2. **Probed with realistic inputs.** I ran about 40 NCBI-style `isolation_source` / `host` combinations through the classifier, along with edge-case inputs for each command.
3. **Reproduced every suspected bug** before changing code. Each one is now a regression test.
4. **Checked scientific choices against the literature** and official NCBI and INSDC documentation. Sources are in [section 11](#11-references).
5. **Compared old and new behaviour** on the same inputs ([section 8](#8-verification-before-vs-after)).
6. **Preferred the simpler design** when two designs gave the same result. `classify.py` went from 286 to 183 lines.

---

## 3. Scientific and methodological fixes

### S1. Specimens now follow the host

**Problem.** The rule `human-specimens` (priority 80) gave L1 = Human to any text containing `blood`, `urine`, `faeces`, `stool`, `tissue`, `wound`, `swab` and similar words. Animal words such as `pig` or `chicken` had priority 40, so the specimen word always won.

| Input | Before | After |
|---|---|---|
| `pig feces` | Human | **Animal** |
| `chicken tissue` | Human | **Animal** |
| `dog urine` | Human | **Animal** |
| `feces` + host `Bos taurus` | review (conflict) | **Animal** |
| `blood` + host `Gallus gallus` | review (conflict) | **Animal** |
| `plant tissue` + host `Zea mays` | Human | **review** (`host_unrecognised`) |
| `pig bite wound` + host `Homo sapiens` | review | review (`host_source_conflict`), unchanged |

**Why this is wrong.** Anatomical sites and specimen types are species-neutral. The cross-species anatomy ontology Uberon models `blood` or `urine` as one concept shared by humans and animals [R1]. In NCBI BioSample, `host` names the organism the sample came from, and `isolation_source` describes the material [R2]. So the host, not the specimen word, decides Human vs Animal. The plant case matters for *Klebsiella*: *K. variicola* is a well-known plant endophyte from maize, rice, banana and sugar cane [R3].

A code comment in the old `_host_source_conflict` described exactly this intended behaviour ("Specimen from an animal host is Animal"), but the code never implemented it.

**Fix.** Lexicon rules can now be marked `host_dependent: true`. A host-dependent rule takes its L1 from other evidence:

- a recognised animal host, or an animal word in the text, gives **Animal**;
- a human host gives **Human**;
- both at once gives **review** (`host_source_conflict`);
- a host the lexicon cannot place, such as a plant, gives **review** (`host_unrecognised`);
- no host at all falls back to the rule's default, **Human**. The reason ends in `+no_host`, so these calls can be counted (see decision D3).

Priorities alone cannot express this. Food must beat animal (`chicken meat` is Food), environment must beat animal (`pig farm wastewater` is Environmental), and specimen must beat environment (`wound irrigation` is Human). Requiring animal to beat specimen as well creates a cycle, so a small amount of code is required.

Files: `classify.py` (`_classify_specimen`), `lexicon.py`, `v2.yaml` (`specimens`, `wound-drainage`, `clinical-sample`).

### S2. Bare "clinical" is a patient sample, not a hospital surface

**Problem.** `clinical`, `clinical isolate` and `clinical isolates` mapped to Environmental with `healthcare_environment = yes`. This applied even when the host was *Homo sapiens*.

| Input | Before | After |
|---|---|---|
| `clinical` | Environmental, hc = yes | **Human**, hc = no |
| `clinical isolate` + `Homo sapiens` | Environmental, hc = yes | **Human** |
| `clinical` + `Bos taurus` | review | **Animal** (veterinary clinical) |
| `clinical environment` | Environmental, hc = yes | unchanged |

**Why this is wrong.** In NCBI, "clinical" means a host-associated sample. The BioSample package for such isolates is titled *"Pathogen: clinical or host-associated"*, and `host` is mandatory in it [R4]. The Pavia study that this project's test notes cite uses "clinical" as a setting next to, and distinct from, "community, animal and environmental" settings [R5]. Treating clinical isolates as hospital-environment swabs would move human patient isolates into the Environmental rollup. That creates the Human–Environmental overlap the study is trying to measure.

**Fix.** The old `clinical-exact` rule became `clinical-sample`, a host-dependent rule (S1). It also accepts `clinical sample(s)` and `clinical specimen(s)`. `clinical environment` and the hospital-fixture rules are unchanged.

### S3. NCBI `source_type` is now used, and NCBI `sample_type` is no longer treated as source evidence

**Problem.** The `sample_type` column was filled from NCBI `sample_type`, then `source_type`, then `tissue`, taking the first non-empty value.

- NCBI `sample_type` is the **material**: "cell culture, mixed culture, tissue sample, whole organism, single cell, metagenomic assembly" [R2]. It says nothing about the One Health source, but it shadowed the next field.
- NCBI `source_type` is a **controlled One Health category**: `human`, `animal`, `food`, `environmental`, `other` [R2]. It was hidden behind `sample_type` and was not recognised by the lexicon anyway.
- When `isolation_source` was blank, any non-blank `sample_type` made the record "unknown_source". So even `source_type = human` could never classify.

| Input | Before | After |
|---|---|---|
| `source_type = human` | review | **Human** |
| `source_type = environmental` | review | **Environmental** |
| `sample_type = whole organism` only | review `unknown_source` | review `missing_source` (correct label) |
| `isolation_source = unknown` + `source_type = food` | review | **Food** |
| host `Homo sapiens` + `source_type = food` | Food (host ignored) | **review** (`host_source_conflict`) |

**Fix.** The alias order is now `source_type`, then `tissue`, then `sample_type`, in both `constants.py` and `ncbi.py`. The lexicon recognises the `source_type` words. The host-vs-source conflict check now compares the host with all source fields, not only `isolation_source`.

### S4. Uninformative `isolation_source` no longer blocks a usable host

**Problem.** A taxon name in `isolation_source` (`Klebsiella pneumoniae`) or an unknown marker (`unknown`) sent the row to review even when `host = Homo sapiens` was present. A blank `isolation_source` with the same host was classified Human. Two equally empty sources were treated differently.

**Why.** A taxon name or "unknown" carries no source information, so it should behave like a blank field. Misplaced values like these are a documented, common problem in BioSample metadata [R6].

**Fix.** These values are dropped as evidence, and the host or `source_type` decides. With nothing else, the row still goes to review as `taxon_only_source` or `unknown_source`. `Bos taurus` typed into `isolation_source` with no host stays review, as your README specifies.

The INSDC missing-value terms (`restricted access`, `missing: lab stock`, `missing: control sample` and the rest [R7]) are now also recognised as unknown.

### S5. Lexicon words that caused silent misclassification (lexicon v1 to v2)

Each change fixes a concrete wrong answer.

| v1 entry | Wrong result | v2 change | Evidence |
|---|---|---|---|
| `air` (any position) | `chicken air sac` became Environmental | Exact `air`, plus `air sample`, `ambient air`, `indoor air`, `outdoor air` | Air sacs are the classic colibacillosis lesion site in poultry [R8, R9] |
| `marine` | `marine fish` became Environmental | Replaced by `sea water`, `marine water` | Marine is an adjective for animals too |
| `dairy` | `dairy farm` became Food | `dairy product(s)` | IFSAC classes dairy *products* as food [R10] |
| `milk` (Food, 75) beat everything animal | `mastitis milk` became Food | New rule `animal-mastitis` (88) on `isolation_source`, `sample_type` and `host_disease` | *Klebsiella* and *E. coli* are major causes of bovine clinical mastitis [R11] |
| no `breast milk` | `breast milk` became Food | Added as a specimen (Human by default) | Human specimen |
| no `meat` / `retail` | `retail chicken` became Animal | `meat` and `retail` added; three phrases they cover removed | IFSAC raw-meat categories [R10] |
| `ventilator` (any position) | `ventilator-associated pneumonia` would be a hospital fixture | Exact `ventilator`, plus `ventilator tubing` / `ventilator circuit` | Pneumonia is a patient infection [R12] |
| `milk tank` in `animal-sites` | Never fired, because `milk` (75) always beat it (65) | Removed | Dead rule |
| `sputum` listed twice; `hospital sink` / `hospital surface` in two rules | Harmless duplicates | Removed | Tidiness |
| `paper_rollup:` and `l1:` blocks in the YAML | Looked configurable but were never read | Removed; the header says the rollup lives in code | Misleading |

New entries also cover common host binomials (`Meleagris gallopavo`, `Anas platyrhynchos`, `Oryctolagus cuniculus`), so a turkey, duck or rabbit host is recognised under S1.

The file is now `v2.yaml`, and every output row records `lexicon_version = 2`. Old and new outputs can therefore be told apart.

### S6. "wound drainage" and "sink drain" in one string no longer picks Environmental

When two rules tied at the top priority with different L1, the old `_allowed_host_env_split` chose the healthcare-environment rule. The only pair this could ever affect was `wound-drainage` vs `hospital-drain`, both at priority 95. That is an arbitrary tie-break, and it contradicts the README rule that conflicting evidence goes to review. The function was removed, so the tie is now a `rule_conflict` and goes to review.

### S7. Review rows no longer claim "not healthcare environment"

Every `Unclassified` row had `healthcare_environment = no`. That is an assertion, not an absence of data. The summary's `healthcare_environment` count therefore counted every unknown row as "no". Unclassified rows now have an empty value.

### S8. Summary statistics were misleading

| Problem | Fix |
|---|---|
| A row with `classify_status = review` was counted as `ok_empty_source` when `resolve_status` was `ok`. | Only classify-`ok` rows are `ok_*`. Review rows are counted under their reason, such as `missing_source`. |
| A cluster with `Human` and `Unclassified` isolates was `mixed = yes`. | `Unclassified` is not a source and is excluded from `l1_set` and `mixed`. |
| `count` meant "isolates" in every section except overlap, where it meant "distinct L1 classes". | `count` now means isolates in every section. The classes are listed in `l1_set`. |
| The same isolate listed twice (by BioSample and by PDT, as in `examples/`) was counted twice, with no warning. | A new `dataset` section reports `rows`, `unique_isolates` and `rows_sharing_a_biosample` (see decision D5). |

---

## 4. Software bugs fixed

| ID | Bug | How it showed up | Fix |
|---|---|---|---|
| B1 | **Non-deterministic local join.** Candidate keys were tried in set order, which depends on Python's hash seed. Later rows also silently overwrote earlier rows that shared a key. | `GCA_000005845.1` received `.2`'s metadata on 3 of 8 runs. | The exact accession is tried first, then related keys in sorted order. Two *different* local rows for one key now give `local_metadata_conflict` (review). |
| B2 | `resolve` took output columns from the **first row only**. | If the first ID was unresolvable, `cluster_id` and `plasmid_id` disappeared from the whole file. | Columns are the union over all rows. |
| B3 | `read_tsv` **crashed** on a row with more cells than the header. A comment line above the header became the header. | `AttributeError: 'list' object has no attribute 'strip'`. | Comment lines are removed first. Trailing empty cells are ignored. Extra non-empty cells raise a clear error with the line number. |
| B4 | An ID-list TSV keyed by `biosample` or `assembly` was not recognised. | The header word `biosample` was read as an input ID. | All ID column aliases are recognised. Otherwise the first column is used, without the header. |
| B5 | Overrides keyed on `GCF_…` did not apply to the `GCA_…` twin, unlike the local join. | A curator's override was silently ignored. | The record's twin and versionless keys are matched. An override on one version does not apply to another version. |
| B6 | Override `paper_rollup` was not validated, and `l1` was case-sensitive. | `paper_rollup = Hospital` was accepted. `animal` was rejected. | Values are validated (`override_invalid_paper_rollup`) and matched case-insensitively. |
| B7 | An unversioned PDT ID (`PDT000002407`) never tried the assembly-name wildcard. | Always `pdt_unresolved`. | The wildcard search runs whenever the exact name misses. |
| B8 | `resolve` blanked any NCBI value **starting with** "unknown". | `unknown river` was erased from the resolved table. | Only whole-value null tokens are blanked. Partial text is kept for the curator, and classify still treats it as unknown. |
| B9 | No retry on NCBI HTTP 429 or 5xx. | A transient rate-limit error became a random review row, and reruns differed. | Up to 3 retries with back-off. 4xx client errors are not retried. |
| B10 | The `ncbi_error:` reason contained the **full request URL**, including `api_key=` and `email=`, and the raw response body. | The API key was written into a TSV you might share. Tabs or newlines from the body could break the table. | The query string is stripped and whitespace collapsed. |
| B11 | `exact` lexicon entries compared raw strings. | `Clinical isolate.` with a trailing period did not match. | Exact entries are compared as word tokens. |
| B12 | A misspelt lexicon field such as `isolation source` was silently ignored. | The rule never fired. | Unknown field names now raise an error when the lexicon loads. |
| B13 | `plot --by` changed only a panel title. Bar order followed input order. The y-axis had fractional ticks for counts. | Figures were not comparable between datasets. | `--by` removed. Fixed category order with zero-count bars. Integer ticks. |
| B14 | No email warning for live NCBI use. | NCBI asks E-utilities users for `tool` and `email` so it can contact them before blocking an IP [R13]. | A warning is printed when `--email` is missing. |
| B15 | A lexicon rule with `l1: Unclassified` produced `status = ok` with `l1 = Unclassified`. | The output was internally inconsistent. | Such rules now mean "send to review". Animal feed uses this. |

---

## 5. Simplifications

| Before | After |
|---|---|
| Four pre-checks in `classify.py` (missing, unknown, taxon-only, animal feed), an allowed-split exception, and a separate `_from_hits` suggestion path. | One evidence filter, one specimen function and one conflict check. The suggestion comes from the same code path, so it can no longer disagree with the classification. |
| The animal-feed phrase list was hard-coded in Python. | It is a lexicon rule (`animal-feed`, `l1: Unclassified`). Vocabulary now lives in data, not code. |
| `best_hits` de-duplicated hits that could never be duplicates. | Removed. |
| Unused `FIELD_ALIASES["any"]`, unused `optional_group_column`, and a copy of `RESOLVE_FIELDS` pasted into `cli.py`. | Removed or reused. |
| `_index_local` called `metadata_view` 9 times per row. | Called once. |
| The overlap and evidence functions in `summarize.py` each built rows by hand. | One `_row()` helper. |

`classify.py` went from 286 to 183 lines, and `summarize.py` from 120 to 103. `ncbi.py` grew from 313 to 350 lines, for the retry logic, URL scrubbing and their explanations.

---

## 6. Decisions left to you (not changed)

These are modelling choices, not errors. Each is defensible in some study designs, so I flagged them rather than changing them.

### D1. Laboratory (and Food) roll up to "Environmental"

The README table defines this, so it is kept. My concern: the One Health definition from the One Health High-Level Expert Panel covers people, animals, plants and ecosystems [R14]. A laboratory, type or reference strain is not a sample from any of those compartments at the time it was sequenced. Counting lab strains as Environmental can inflate the Environmental group and create false Human–Environmental overlap. **Suggestion:** roll Laboratory to `Unclassified`, or exclude it from overlap analyses. It is a one-line change in `derive_paper_rollup`.

### D2. Unqualified drains count as healthcare environment

`sink drain`, `floor drain`, `drain water` and `p-trap` get `healthcare_environment = yes` with no hospital word. The README says bare `sink` is *not* a hospital fixture, so this is inconsistent. In *Klebsiella* hospital studies the assumption is usually right, because hospital drains, sinks and faucets are the main reservoirs [R15]. But floor drains are also the classic *Listeria* sample in food plants [R16]. The paper rollup is Environmental either way, so only the healthcare flag is affected. **Suggestion:** if your data mixes settings, require a hospital word, as is already done for sinks.

### D3. A specimen with no host defaults to Human

This keeps your tested behaviour: `wound drainage` and `blood` alone are Human. It is an inference, not an observation. Every such call now has a reason ending in `+no_host`, so the summary shows how many there are. **Alternative:** send them to review. That is safer, but it adds curation work.

### D4. A host binomial in `isolation_source` goes to review

`Bos taurus` in `isolation_source` with no host stays review, as your README specifies. A human curator would usually read it as a bovine sample. Kept as specified.

### D5. Counts are per input row, not per isolate

Duplicates are reported in the `dataset` section but not removed, so row-level outputs still match the input one-to-one. If you report isolate counts in a paper, de-duplicate by BioSample first.

---

## 7. Known limitations still present

- **The lexicon is a seed.** Plurals do not match (`product` vs `products`). Many specimen names are missing, such as `synovial fluid` or `nasopharyngeal aspirate`. These go to review, which is safe but adds curation work.
- **The taxon-name guard is a heuristic.** Any 2–3-word phrase it does not know can be flagged as `taxon_only_source`. The outcome is still review, but the label can be wrong. *Updated November 2026:* a lexicon rule written for source text now overrides the guess, so adding a phrase such as `bedpan washer` to your own lexicon works; words ending in -ing, -tion, -sion, -ment or -ness and missing-value terms (`not collected`) are no longer flagged ([CHANGELOG.md](CHANGELOG.md)).
- **Priority silently resolves mixed strings.** For example, `soil contaminated with faeces` with no host becomes Human, because the specimen rule outranks soil. The README promises review for conflicts, but cross-priority mixes are resolved by rank.
- **`farm` counts as an animal word.** So `farm worker stool` with no host becomes Animal. With host `Homo sapiens` it is correctly a conflict.
- **Shared cluster or plasmid is not transmission.** `mixed = yes` only says that one cluster ID contains isolates from several sources. Sequence types and PopPUNK clusters are lineage-level groupings [R17]. Transmission claims need SNP thresholds; about 21 SNPs separated hospital clusters in the European *K. pneumoniae* survey [R18]. One Health studies with shared STs have still found the human and livestock populations distinct [R19]. The Pavia study found only occasional transmission between clinical and non-clinical settings [R5].
- **Header lines starting with `#`** are treated as comments. NCBI `assembly_summary.txt` uses `#` before its header, so strip it first.
- **`docs/` is git-ignored but still tracked.** Adding a folder to `.gitignore` does not untrack files already committed. If you meant to drop them, run `git rm -r --cached docs`. The notes there also describe bugs fixed in commit `37d97b0`. They cite the Pavia paper as "Nature Microbiology 13, 1520–1531"; the correct reference is *Nat Microbiol* 7, 2054–2067 (2022) [R5].

---

## 8. Verification: before vs after

**Tests.** 73 before, 133 after, all passing (`uv run pytest`). The new tests are in `test_classify_host_logic.py`, `test_lexicon.py` and `test_io.py`, plus additions to the resolve, override, summary and CLI tests.

**Determinism.** The same local join ran under 8 different hash seeds. Before the fix: 5 runs said `urine` and 3 said `river water`. After: 8 of 8 say `urine`.

**End-to-end.** The README example commands and the test fixtures run cleanly offline.

**Classification on identical inputs** (original commit `5c0b5c4` vs this version):

| isolation_source | host / other | Before | After |
|---|---|---|---|
| pig feces | | Human | **Animal** |
| chicken tissue | | Human | **Animal** |
| dog urine | | Human | **Animal** |
| feces | Bos taurus | review | **Animal** |
| blood | Gallus gallus | review | **Animal** |
| urine | chicken | review | **Animal** |
| plant tissue | Zea mays | Human | **review** |
| clinical | | Environmental (hc yes) | **Human** |
| clinical isolate | Homo sapiens | Environmental (hc yes) | **Human** |
| clinical | Bos taurus | review | **Animal** |
| Klebsiella pneumoniae | Homo sapiens | review | **Human** |
| mastitis milk | | Food | **Animal** |
| milk | Bos taurus, host_disease mastitis | review | **Animal** |
| breast milk | | Food | **Human** |
| chicken air sac | | Environmental | **Animal** |
| marine fish | | Environmental | **Animal** |
| dairy farm | | Food | **Animal** |
| retail chicken | | Animal | **Food** |
| *(blank)* | source_type human | review | **Human** |
| *(blank)* | source_type environmental | review | **Environmental** |
| unknown | source_type food | review | **Food** |
| *(blank)* | Homo sapiens, source_type food | Food | **review** |
| wound drainage, sink drain | | Environmental | **review** |
| Clinical isolate. | | review | **Human** |

**Unchanged controls:** `wound drainage` (Human), `hospital drain` (Environmental, hc yes), `hospital sink` + *Homo sapiens* (Environmental, hc yes), `clinical environment`, `unknown`, bare `sink`, `Animal - food producing`, `animal feed`, `Escherichia coli`, `Bos taurus`, `cloacal swab` + *Gallus gallus*, `retail meat`, `type strain`, `river water`.

---

## 9. What changed for you as a user

- **Reason strings were renamed.** `human-specimens` is now `specimens+<decider>`, for example `specimens+human-host`, `specimens+animal-hosts` or `specimens+no_host`. `clinical-exact` is now `clinical-sample+…`, and `animal_feed` is now `animal-feed`. New reasons are `host_unrecognised` and `local_metadata_conflict`, plus `override_invalid_paper_rollup`. Update any script that filters on the old names.
- **`healthcare_environment`** is empty, not `no`, for Unclassified rows.
- **`lexicon_version`** is now `2`. Results made with v1 cannot be reproduced with the current code; use git commit `5c0b5c4` for that.
- **`onehealth plot --by`** no longer exists. The figure always shows both panels, as it did before.
- **The summary** has a new `dataset` section, and `count` in the overlap sections now means isolates.
- **Live `resolve`** prints a warning without `--email`, and retries rate-limited requests.

---

## 10. Files changed

| File | Change |
|---|---|
| `src/onehealth/classify.py` | Rewritten around S1, S4, S6 and S7; simpler. |
| `src/onehealth/data/lexicon/v2.yaml` | Renamed from `v1.yaml`. Changes in S2, S3 and S5; documented header. |
| `src/onehealth/lexicon.py` | `host_dependent` flag, token-based exact match, field validation. |
| `src/onehealth/constants.py` | Alias order, INSDC term, lexicon version 2, `Unclassified` in rollup values. |
| `src/onehealth/text.py` | INSDC `missing:` terms; vocabulary for the taxon guard. |
| `src/onehealth/io_tsv.py` | B3, B4; assembly twin keys for overrides (B5). |
| `src/onehealth/resolve.py` | B1 (deterministic join and conflict detection). |
| `src/onehealth/ncbi.py` | B7, B8, B9, B10; `source_type` first (S3). |
| `src/onehealth/overrides.py` | B5, B6. |
| `src/onehealth/summarize.py` | S8. |
| `src/onehealth/plot.py` | B13. |
| `src/onehealth/cli.py` | B2, B13, B14. |
| `README.md` | Describes the corrected behaviour. |
| `tests/` | 3 new files, 5 updated; 133 tests. |

---

## 11. References

Literature was retrieved from PubMed; NCBI and INSDC documents from their official sites.

- **R1** Mungall CJ et al. (2012) Uberon, an integrative multi-species anatomy ontology. *Genome Biol* 13:R5. <https://doi.org/10.1186/gb-2012-13-1-r5>
- **R2** NCBI BioSample attribute definitions (`host`, `isolation_source`, `sample_type`, `source_type`, `tissue`). <https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/>
- **R3** Rosenblueth M et al. (2004) *Klebsiella variicola*, a novel species with clinical and plant-associated isolates. *Syst Appl Microbiol* 27:27–35. <https://doi.org/10.1078/0723-2020-00261>
- **R4** NCBI BioSample package *Pathogen: clinical or host-associated* (Pathogen.cl.1.0). <https://www.ncbi.nlm.nih.gov/biosample/docs/packages/Pathogen.cl.1.0/>
- **R5** Thorpe HA et al. (2022) A large-scale genomic snapshot of *Klebsiella* spp. isolates in Northern Italy reveals limited transmission between clinical and non-clinical settings. *Nat Microbiol* 7:2054–2067. <https://doi.org/10.1038/s41564-022-01263-0>
- **R6** Gonçalves RS, Musen MA (2019) The variable quality of metadata about biological samples used in biomedical experiments. *Sci Data* 6:190021. <https://doi.org/10.1038/sdata.2019.21>
- **R7** INSDC missing value reporting. <https://www.insdc.org/technical-specifications/missing-value-reporting/>
- **R8** Kathayat D et al. (2021) Avian pathogenic *Escherichia coli* (APEC): an overview of virulence and pathogenesis factors, zoonotic potential, and control strategies. *Pathogens* 10:467. <https://doi.org/10.3390/pathogens10040467>
- **R9** Huff GR et al. (2000) Turkey osteomyelitis complex (air sac inoculation and air sacculitis scoring). *Poult Sci* 79:1050–1056. <https://doi.org/10.1093/ps/79.7.1050>
- **R10** Richardson LC et al. (2017) An updated scheme for categorizing foods implicated in foodborne disease outbreaks: a tri-agency collaboration (IFSAC). *Foodborne Pathog Dis* 14:701–710. <https://doi.org/10.1089/fpd.2017.2324>
- **R11** Nery Garcia BL et al. (2024) ESBL-producing and other antimicrobial-resistant Gram-negative pathogens isolated from bovine mastitis: a One Health perspective. *Antibiotics* 13:391. <https://doi.org/10.3390/antibiotics13050391>
- **R12** Paczosa MK, Mecsas J (2016) *Klebsiella pneumoniae*: going on the offense with a strong defense. *Microbiol Mol Biol Rev* 80:629–661. <https://doi.org/10.1128/MMBR.00078-15>
- **R13** NCBI. A General Introduction to the E-utilities (usage guidelines: 3 requests/s without an API key, 10 with one; `tool` and `email`). <https://www.ncbi.nlm.nih.gov/books/NBK25497/> and <https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities>
- **R14** One Health High-Level Expert Panel; Adisasmito WB et al. (2022) One Health: a new definition for a sustainable and healthy future. *PLoS Pathog* 18:e1010537. <https://doi.org/10.1371/journal.ppat.1010537>
- **R15** Kizny Gordon AE et al. (2017) The hospital water environment as a reservoir for carbapenem-resistant organisms causing hospital-acquired infections: a systematic review. *Clin Infect Dis* 64:1435–1444. <https://doi.org/10.1093/cid/cix132>
- **R16** Chowdhury B, Anand S (2023) Environmental persistence of *Listeria monocytogenes* and its implications in dairy processing plants. *Compr Rev Food Sci Food Saf* 22:4573–4599. <https://doi.org/10.1111/1541-4337.13234>
- **R17** Lees JA et al. (2019) Fast and flexible bacterial genomic epidemiology with PopPUNK. *Genome Res* 29:304–316. <https://doi.org/10.1101/gr.241455.118>
- **R18** David S et al. (2019) Epidemic of carbapenem-resistant *Klebsiella pneumoniae* in Europe is driven by nosocomial spread. *Nat Microbiol* 4:1919–1929. <https://doi.org/10.1038/s41564-019-0492-8>
- **R19** Ludden C et al. (2019) One Health genomic surveillance of *Escherichia coli* demonstrates distinct lineages and mobile genetic elements in isolates from humans versus livestock. *mBio* 10:e02693-18. <https://doi.org/10.1128/mBio.02693-18>
