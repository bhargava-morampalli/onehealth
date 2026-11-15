# classify

Assigns each row a One Health class from its host and isolation source, using the [lexicon](../guides/lexicon.md). Rows with missing, unknown or conflicting evidence go to a review queue instead of a guess.

```bash
onehealth classify work/resolved.tsv -o work/classified.tsv
```

The input can be the output of [`resolve`](resolve.md) or your own metadata table; see [Input files](../guides/input-files.md).

## Options

| Option | What it does |
|---|---|
| `-o`, `--out` | Output TSV (required). |
| `--overrides` | A curator overrides file. A complete override always wins; see [Overrides](#overrides). |
| `--lexicon` | A lexicon YAML of your own. The default is the one lexicon that ships with onehealth; see [The lexicon](../guides/lexicon.md#version). |
| `--review` | Where to write the review queue (default `<out>.review.tsv`). |

`classify` prints how many rows were classified and which lexicon version was used.

## The classes

**L1** is one of `Human`, `Animal`, `Environmental`, `Food`, `Laboratory`, or `Unclassified` for rows sent to review.

**healthcare_environment** is `yes` for hospital surfaces, sinks, drains, wastewater and other healthcare fixtures, `no` for other classified rows, and empty for `Unclassified` rows. A hospital sink is Environmental with `healthcare_environment = yes`, not Human.

**paper_rollup** collapses L1 into the three classes most One Health papers report:

| L1 | healthcare_environment | paper_rollup |
|---|---|---|
| Human | no | Human |
| Animal | any | Animal |
| Environmental, Food or Laboratory | any | Environmental |
| any | yes | Environmental |
| Unclassified | | Unclassified |

## What the classes are based on

The lexicon reads four fields of each row:

| Field | Taken from the first non-empty of |
|---|---|
| isolation_source | `isolation_source`, `isolation source`, `source`, `sample_source` |
| host | `host` |
| sample_type | `source_type`, `tissue`, `sample_type` |
| host_disease | `host_disease`, `host disease` |

NCBI's `source_type` (`human`, `animal`, `food`, `environmental`) is a One Health category, so it is read before `tissue` and `sample_type`. NCBI's material type (`whole organism`, `pure culture`) is never treated as evidence.

One more column is read, but never by the lexicon: the BioSample package the record was submitted under, from the first non-empty of `package`, `attribute_package`, `biosample_package` and `epi_type`. [`resolve`](resolve.md) writes it as `package`, and NCBI Pathogen Detection tables have `attribute_package` and `epi_type`. It is used in one place only, for a specimen with no host (below).

A table with none of these columns stops with an error, because every row would go to review. The usual cause is an ID list given to `classify` instead of [`resolve`](resolve.md). With `--overrides` and an ID column, it is a warning instead, and only the rows the overrides name are classified.

### Specimens follow the host

Body sites and specimens (blood, urine, faeces, wound, tissue, catheter, "clinical") occur in any host:

- With an animal host or an animal word (`pig faeces`, urine from *Canis lupus familiaris*), the row is **Animal**.
- With a human host, it is **Human**.
- If another field names food, an environment or a laboratory, such as `source_type` food, it goes to review as `host_source_conflict`. Words in the specimen's own text do not count: `stool from a food handler` is still Human.
- With both, it goes to review as `host_source_conflict`.
- With a host the lexicon does not recognise, such as a plant, it goes to review as `host_unrecognised`. Any animal NCBI Taxonomy names is recognised, not only those on the lexicon's word list; see [Hosts by taxonomy](../guides/lexicon.md#hosts-by-taxonomy).
- With no host at all, it is **Human**, and the reason ends in `+no_host` so you can count these calls.
- With no host at all in a record its submitter filed as not clinical, it goes to review as `no_host_env_package`. That is the BioSample package `Pathogen.env` ("Pathogen: environmental/food/other"), or `epi_type` `environmental/other` in a Pathogen Detection table. `feces` or `wound` in such a record is often an animal's. In a sample of 159,045 NCBI Pathogen Detection records, 586 specimens had no host and a known package, and 199 of them (34%) were filed this way. Any other package, or no package column, keeps the Human default.

### Values that are not evidence

- **Missing-value words** in `isolation_source` (`unknown`, `missing`, `not collected`, `NaN`, ...) are ignored, and the host or `source_type` decides. The full list is in [Input files](../guides/input-files.md#missing-values). `unknown` followed by a word that names something still counts, in any field: `unknown animal` is Animal, while `unknown source` and `unknown environment` mean nothing.
- **A bare scientific name** (`Escherichia coli`, `Bos taurus`) in `isolation_source` is ignored the same way: taxonomy belongs in `host`.

With nothing else to go on, the row goes to review.

## Output

The input columns, followed by:

| Column | Contents |
|---|---|
| `l1` | The class |
| `healthcare_environment` | `yes`, `no`, or empty for Unclassified |
| `paper_rollup` | The three-class rollup |
| `classify_status` | `ok` or `review` |
| `classify_reason` | Why (below) |
| `lexicon_version` | The lexicon version that made the call. For a lexicon of your own: its `version`, or `custom-<digest>` if it has none, or `<number>-modified-<digest>` for an edited copy of a shipped one; see [The lexicon](../guides/lexicon.md#version). |
| `override_applied` | `yes` when an override decided |
| `lexicon_suggestion` | For an override, the lexicon's own answer; for a review row, the competing classes, if any |

### classify_reason

| Reason | Status | Meaning |
|---|---|---|
| a rule id, e.g. `hospital-drain` | ok | That lexicon rule matched |
| `<specimen rule>+<host rule>`, e.g. `specimens+human-host` | ok | A specimen, classed by its host |
| `<specimen rule>+no_host` | ok | A specimen with no host, defaulted to Human |
| `override` | ok | An override decided |
| `host_source_conflict` | review | Host and source point to different classes |
| `host_unrecognised` | review | A specimen from a host the lexicon cannot place |
| `no_host_env_package` | review | A specimen with no host, in a record filed as environmental, food or other |
| `rule_conflict:<rules>` | review | Rules of equal priority disagree |
| `animal-feed`, `fixture-no-setting` | review | A rule that sends matches to review: animal feed, or a fixture such as a drain or sink named without its setting |
| `taxon_only_source` | review | `isolation_source` is only a scientific name, and nothing else decides |
| `unknown_source` | review | `isolation_source` is a missing-value word, and nothing else decides |
| `missing_source` | review | No isolation source and no host |
| `unmatched` | review | No rule matched |
| `override_empty`, `override_conflict`, `override_invalid_l1`, `override_invalid_healthcare`, `override_invalid_paper_rollup` | review | A problem with the override for this row |

One exception to `host_source_conflict`: a healthcare-environment source with a human host is classified as healthcare environment, because submitters often record *Homo sapiens* as the host of hospital sink and surface swabs.

## Overrides

An overrides file records a curator's decisions. It is a TSV with these columns:

| Column | Required | Contents |
|---|---|---|
| `id` (or `query_id`, `accession`, `biosample`, `assembly`, `sra`, `pdt`, `strain`, `isolate`) | yes | The isolate |
| `l1` (or `source_l1`, `classification`) | yes | The class |
| `healthcare_environment` | no | `yes` or `no` (default `no`) |
| `paper_rollup` | no | Defaults to the rollup of `l1` and `healthcare_environment` |
| `note` (or `reason`, `comment`) | no | Kept for your records |

- A complete override always wins over the lexicon, and the lexicon's answer is kept in `lexicon_suggestion`.
- Class names are case-insensitive.
- An override on an assembly also covers its GCA/GCF twin, and an unversioned accession covers every version.
- A strain or isolate name matches only a table that has no BioSample, assembly or ID column. `classify` warns about every override ID that matches no row, since it is not applied.
- A row with an empty or invalid `l1`, or two rows that disagree about the same isolate, sends that isolate to review rather than silently using the lexicon.
- A file without an ID or `l1` column is an error.

## Working through the review queue

The review queue (`<out>.review.tsv`) lists every row `classify` would not guess, with the reason.

1. Sort by `classify_reason`. `missing_source` rows have no information to go on: check the BioSample or your lab records.
2. For `unmatched` rows, look for a word the lexicon should know. If it is general, add it to a copy of the lexicon and pass `--lexicon`; see [The lexicon](../guides/lexicon.md).
3. For individual isolates, add a row to an overrides file.
4. Measure the effect with [`validate`](validate.md) before relying on a changed lexicon.
