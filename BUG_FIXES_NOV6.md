# Bug-Fix Audit Log — 6 November 2026 (`bug-fixes-nov6`)

This document records the bugs and edge cases identified during the technical audit of `onehealth` (`main` at `035912f`), the root cause of each issue, the minimal fix implemented, and the regression tests added to verify the fix.

---

## 1. Taxon-Only or Unknown `isolation_source` Shadowing Later Rule Fields

- **Component:** `src/onehealth/classify.py` (`_classify_with_lexicon`), `src/onehealth/lexicon.py` (`match_rules`)
- **Problem Description:**
  `match_rules` iterates through each rule's `fields` tuple in order (for example, `animal-hosts` reads `[host, isolation_source, sample_type]` and `specimens` reads `[isolation_source, sample_type]`) and stops at the first field that matches (`break`), producing at most one `RuleHit` per rule. Previously, `_classify_with_lexicon` ran `match_rules(view, lexicon)` on the raw view and afterwards stripped hits where `hit.field == "isolation_source"` if `isolation_source` was an unknown marker (`is_unknown_value(source)`) or a bare scientific name (`_taxon_only(source, hits)`).
  When `isolation_source` matched a multi-field rule and was then discarded, that rule never evaluated its remaining fields (`sample_type`), violating the documented behavior that an unknown or taxon-only `isolation_source` is treated *"as if the field were blank"*.
  - *Reproduction before fix:*
    - `classify_record({"isolation_source": "Bos taurus", "source_type": "animal"}, lexicon)` returned `Unclassified` (`taxon_only_source`) instead of `Animal` (`animal-hosts`), because `"Bos taurus"` matched `animal-hosts` on `isolation_source`, shadowed `sample_type="animal"`, and was then dropped as taxon-only.
    - `classify_record({"isolation_source": "unknown tissue", "sample_type": "blood"}, lexicon)` returned `Unclassified` (`unknown_source`) instead of `Human` (`specimens+no_host`), because `"unknown tissue"` matched `specimens` on `isolation_source`, shadowed `sample_type="blood"`, and was then dropped by `is_unknown_value`.
- **Solution Implemented:**
  Blank any field in `RULE_FIELDS` that satisfies `is_unknown_value` before calling `match_rules`, and if `isolation_source` is detected as `_taxon_only`, re-run `match_rules` with `isolation_source` set to `""`.
- **How It Was Implemented:**
  - In `src/onehealth/classify.py` (`_classify_with_lexicon`), constructed `cleaned = {field: "" if is_unknown_value(view.get(field)) else view.get(field, "") for field in RULE_FIELDS}` prior to calling `match_rules(cleaned, lexicon)`.
  - When `_taxon_only(cleaned["isolation_source"], hits)` is true, re-evaluated `hits = match_rules({**cleaned, "isolation_source": ""}, lexicon)`.
  - `_no_evidence_reason(view)` continues to inspect the original `view` so review reason codes (`unknown_source`, `taxon_only_source`, `missing_source`) are unchanged when no other field matches.
  - Added `test_taxon_or_unknown_in_source_does_not_shadow_sample_type` in `tests/test_classify_host_logic.py`.

---

## 2. Unknown Markers Not Filtered on `host`, `sample_type`, or `host_disease`

- **Component:** `src/onehealth/classify.py` (`_classify_with_lexicon`)
- **Problem Description:**
  `is_unknown_value` recognises partial unknown markers such as `"unknown animal"`, `"unknown food"`, `"unknown environment"`, and INSDC `"missing: <reason>"` strings. However, `_classify_with_lexicon` only checked `is_unknown_value` on `isolation_source`, leaving `host`, `sample_type`, and `host_disease` unfiltered.
  - *Reproduction before fix:*
    - `classify_record({"host": "unknown animal"}, lexicon)` classified the isolate as `Animal` (`classify_status="ok"`).
    - `classify_record({"sample_type": "unknown food"}, lexicon)` classified the isolate as `Food` (`classify_status="ok"`).
- **Solution Implemented:**
  Apply `is_unknown_value` across all four classification fields (`RULE_FIELDS = ("isolation_source", "host", "sample_type", "host_disease")`) before matching lexicon rules.
- **How It Was Implemented:**
  - Handled by the same `cleaned` dictionary comprehension in `src/onehealth/classify.py` (`_classify_with_lexicon`) described in Issue #1.
  - Added `test_unknown_markers_in_host_or_sample_type_stay_unclassified` in `tests/test_classify_host_logic.py`.

---

## 3. Specimen Rules (`host_dependent: true`) Bypassing Non-Host Source Conflicts

- **Component:** `src/onehealth/classify.py` (`_classify_specimen`)
- **Problem Description:**
  When the winning rule is `host_dependent: true` (`wound-drainage`, `clinical-sample`, `specimens`), `_classify_with_lexicon` delegates directly to `_classify_specimen` and returns without calling `_host_source_conflict(hits)`. Inside `_classify_specimen`, only hits with `hit.rule.l1 in {"Human", "Animal"}` were inspected (`deciders`), silently ignoring lower-priority `Food`, `Environmental`, or `Laboratory` hits in other fields (such as `source_type="food"` or `source_type="environmental"`).
  - *Reproduction before fix:*
    - `classify_record({"host": "Homo sapiens", "source_type": "food"}, lexicon)` returned `review` (`host_source_conflict`), but adding a specimen (`{"isolation_source": "urine", "host": "Homo sapiens", "source_type": "food"}`) silently classified the row as `Human` (`ok`).
    - `classify_record({"isolation_source": "blood", "source_type": "environmental"}, lexicon)` classified the row as `Human` (`specimens+no_host`) instead of `host_source_conflict`.
- **Solution Implemented:**
  In `_classify_specimen`, inspect non-`host_dependent` hits (excluding `Unclassified` review rules and sub-phrases of a matched specimen phrase on the same field, such as `"milk"` inside `"breast milk"`). If any remaining hit assigns a non-host L1 (`Food`, `Environmental`, or `Laboratory`), return `review` with `classify_reason="host_source_conflict"` and the competing L1 classes in `lexicon_suggestion`.
- **How It Was Implemented:**
  - Updated `_classify_specimen` in `src/onehealth/classify.py` to filter `other_hits` not subsumed by `specimen_hits` on the same field (`contains_phrase(tokenize(s.matched), tokenize(hit.matched))`) and check `non_host_l1 = {hit.rule.l1 for hit in other_hits if hit.rule.l1 not in {"Human", "Animal"}}`.
  - Added `test_specimen_with_conflicting_non_host_source_type_goes_to_review` in `tests/test_classify_host_logic.py`.

---

## 4. Incomplete Support for `strain` / `sample`-Keyed Tables in `overrides` and `validate`

- **Component:** `src/onehealth/io_tsv.py` (`row_match_ids`), `src/onehealth/overrides.py` (`_OVERRIDE_ID_COLUMNS`), `src/onehealth/validate.py` (`_TRUTH_ID_COLUMNS`), `src/onehealth/cli.py` (`_check_classify_columns`)
- **Problem Description:**
  PR #1 taught `isolate_key` and `report` / `summarize` to identify isolates by `_ISOLATE_NAME_ALIASES = ("strain", "isolate", "sample", "sample_name", "name")` when a table has no NCBI accession or `id` column. However, `row_match_ids` only checked accession and `id` columns, returning an empty set `{}` for strain- or sample-keyed rows. As a result, `classify --overrides` and `validate` could not match rows in tables keyed by `strain` or `sample`.
  - *Reproduction before fix:*
    - `lookup_override({"strain": "KP1", "isolation_source": "blood"}, overrides)` returned `None`.
    - `validate([{"strain": "KP1", "paper_rollup": "Human"}], {"KP1": "Human"})` matched `0` isolates (`unmatched=["KP1"]`).
- **Solution Implemented:**
  Fall back to `{isolate_key(row)}` in `row_match_ids(row)` when a row has no accession or `id` column, and accept strain/isolate columns in `_OVERRIDE_ID_COLUMNS`, `_TRUTH_ID_COLUMNS`, and `_check_classify_columns`.
- **How It Was Implemented:**
  - In `src/onehealth/io_tsv.py` (`row_match_ids`), added fallback `if not matched and (key := isolate_key(row)): return {key}` so rows with accession/ID columns are unaffected while strain/sample-only rows yield their `isolate_key`.
  - Added `"strain"` and `"isolate"` to `_OVERRIDE_ID_COLUMNS` in `src/onehealth/overrides.py`, `_ISOLATE_NAME_ALIASES` to `_TRUTH_ID_COLUMNS` in `src/onehealth/validate.py`, and `_ISOLATE_NAME_ALIASES` to the override ID check in `src/onehealth/cli.py` (`_check_classify_columns`).
  - Added `test_overrides_apply_to_tables_keyed_by_strain_or_sample` in `tests/test_overrides.py` and `test_validate_matches_strain_tables_and_deduplicates_shared_biosamples` in `tests/test_validation.py`.

---

## 5. Single-Linkage Clusters in `linkage()` Including `same_host` Pairs

- **Component:** `src/onehealth/analysis.py` (`linkage` -> `edges_at`)
- **Problem Description:**
  In `linkage()`, `edges_at(t)` built connected components via `_components(len(isolates), pair_a[hit], pair_b[hit])` using all pairs with `pair_d <= t` before filtering out `same_host` pairs (`keep = hit & ~same_host`). Consequently, two isolates from the same individual host (`host_id`) within the SNP threshold formed a 2-isolate cluster in `link.clusters` with zero qualifying edges (`edges = []`), or attached to an existing cluster as a disconnected node with no edge in `fig_clusters`'s minimum spanning tree.
  - *Reproduction before fix:*
    - Two `Human` isolates sharing `host_id="P1"` at 1 SNP distance produced a cluster of 2 isolates with `edges=[]` in `link.clusters`.
- **Solution Implemented:**
  Exclude `same_host` pairs before computing single-linkage connected components in `edges_at(t)`.
- **How It Was Implemented:**
  - In `src/onehealth/analysis.py` (`edges_at`), changed `hit = pair_d <= t; clusters = _components(len(isolates), pair_a[hit], pair_b[hit]); keep = hit & ~same_host` to `keep = (pair_d <= t) & ~same_host; clusters = _components(len(isolates), pair_a[keep], pair_b[keep])`.
  - Added `test_linkage_clusters_exclude_same_host_only_pairs` in `tests/test_analysis.py`.

---

## 6. `dereplicate()` Collapsing `Unclassified` Isolates and Counting Them in `clusters_collapsed`

- **Component:** `src/onehealth/analysis.py` (`dereplicate`)
- **Problem Description:**
  In `dereplicate()`, `Unclassified` isolates have `iso.source == ""`. The linkage condition `link = sources[a] == sources[b]` evaluated to `True` (`"" == ""`) for pairs of `Unclassified` isolates within the SNP threshold, collapsing them together and incrementing `summary.clusters_collapsed` even though `isolates_before` and `isolates_after` only count classified isolates (`len(dataset.classified)`).
  - *Reproduction before fix:*
    - Running `dereplicate` on 1 `Human` isolate and 2 `Unclassified` isolates at 0 SNPs reported `isolates_before=1, isolates_after=1, clusters_collapsed=1` and reduced `derep.unclassified` from 2 to 1.
- **Solution Implemented:**
  Only link pairs where both isolates have a non-empty (classified) source and share the same source.
- **How It Was Implemented:**
  - In `src/onehealth/analysis.py` (`dereplicate`), changed `link = sources[a] == sources[b]` to `link = (sources[a] != "") & (sources[a] == sources[b])`.
  - Added `test_dereplicate_does_not_collapse_unclassified_isolates` in `tests/test_analysis.py`.

---

## 7. Truncated $2 \times 2$ Square Distance Matrix Misread as 3-Column Long Format

- **Component:** `src/onehealth/io_tsv.py` (`read_distance_matrix`)
- **Problem Description:**
  `read_distance_matrix` checked `if names and _first_cells(path, delimiter) == names:` to route square matrices to `_read_square`, falling back to `_read_long` otherwise. A $2 \times 2$ square matrix has 3 columns (`snp-dists\tA\tB`). When truncated to a single row (`A\t0\t7`), `_first_cells` (`['A']`) did not equal `names` (`['A', 'B']`), so `read_distance_matrix` fell back to `_read_long`, treated the square header as a long-format header, and parsed `A\t0\t7` as distance `7.0` between isolate `'A'` and a phantom isolate named `'0'`.
- **Solution Implemented:**
  Route matrices with more than 3 header columns (`len(names) > 2`), or 3-column matrices whose first data row starts with the first column's header name (`first[0] == names[0]`) and whose last header cell is non-numeric (`not _is_float(names[1])`), to `_read_square` so they raise `InputError: distance matrix is not square`.
- **How It Was Implemented:**
  - Added `_is_float` and tightened the `_read_square` dispatch condition in `src/onehealth/io_tsv.py` (`read_distance_matrix`).
  - Added `test_read_distance_matrix_rejects_truncated_2x2_square` in `tests/test_analysis.py`.

---

## 8. Duplicate BioSample Rows in `validate()`, Conflicting Alias Lineages in `apply_lineages()`, and Unstripped CSV Headers in `_table_pairs()`

- **Component:** `src/onehealth/validate.py` (`validate`), `src/onehealth/analysis.py` (`apply_lineages`), `src/onehealth/io_tsv.py` (`_table_pairs`)
- **Problem Description:**
  1. `validate()` deduplicated matched rows by Python object identity (`id(row)`). When a classified table contained two rows sharing a BioSample (for example, two assemblies `GCA_1` and `GCA_2` for `SAMN1`), they were distinct `dict` objects, so `validate()` counted the isolate twice if the truth file keyed on assembly and failed to detect conflicting truth labels across the two assemblies.
  2. `apply_lineages()` built `by_key = {key: _clean(lineages[name]) for name, key in mapping.items()}` without checking whether two distinct external names mapping to the same `iso.key` (such as `SAMN1` and `GCA_000000001.1`) carried conflicting lineages.
  3. `_table_pairs()` used `csv.DictReader` for comma-separated lineage tables without stripping whitespace around header column names (unlike `read_tsv`).
- **Solution Implemented:**
  1. Key `seen` in `validate()` by `isolate_key(row) or str(id(row))` so rows sharing a BioSample are counted once per isolate and checked for truth-label conflicts.
  2. Check for conflicting non-empty lineages mapped to the same `iso.key` in `apply_lineages()` and raise `InputError`.
  3. Strip header and cell whitespace in the `csv.DictReader` branch of `_table_pairs()`.
- **How It Was Implemented:**
  - Updated `seen` key in `src/onehealth/validate.py` (`validate`).
  - Updated `by_key` loop in `src/onehealth/analysis.py` (`apply_lineages`).
  - Stripped keys and values in `src/onehealth/io_tsv.py` (`_table_pairs`).
  - Added regression assertions in `tests/test_validation.py` (`test_validate_matches_strain_tables_and_deduplicates_shared_biosamples`) and `tests/test_analysis.py` (`test_apply_lineages_rejects_conflicting_aliases_and_strips_csv_headers`).

---

## 9. Documentation & Build Metadata Polish

- **`README.md`:** Updated the seeded lexicon reference from `src/onehealth/data/lexicon/v3.yaml` to `src/onehealth/data/lexicon/v4.yaml` (and noted `--lexicon v3` or `v2` for historical reproduction).
- **`CHANGELOG.md`:** Dated the merged PR #10 (`v4` lexicon) entry and recorded the 6 November 2026 bug fixes.
- **`pyproject.toml`:** Widened the `uv_build` upper bound from `<0.10.0` to `<0.15.0` so builds do not warn or fail on newer `uv` versions (`0.10+`).
