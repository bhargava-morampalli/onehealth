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
| `--review-groups` | Where to write the review queue by distinct source text (default `<out>.review_groups.tsv`); see [Working through the review queue](#working-through-the-review-queue). |
| `--decisions` | A review-groups table with its `l1` column filled in. Each decision applies to every row with that source text; see [Decisions](#decisions). |

`classify` prints how many rows were classified and which lexicon version was used.

!!! warning "The classes follow NCBI"
    A record's own `source_type` is its class, farm settings are Environmental, a meat species with no host (`chicken`, `fish`, `poultry`) is Food, and a specimen with no host is Human. [Definitions](../guides/definitions.md) gives the numbers behind each of these and says how to choose differently for your study.

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

NCBI's material type (`whole organism`, `pure culture`) is never treated as evidence.

### `source_type` decides

**When the `source_type` column of a row is `human`, `animal`, `food` or `environmental`, that is the row's class.** It is NCBI's controlled One Health category, chosen by the submitter, and the lexicon is not allowed to disagree with it or to send such a row to review. The reason says what happened:

- the lexicon's own rule, when the free text is silent or agrees;
- `source_type`, when the free text could not settle it (a host that conflicts with the source, a rule that sends to review);
- `source_type_over_text`, when the free text said something else. `lexicon_suggestion` holds what the lexicon would have said.

Any other value (`other`, `not provided`, empty) declares nothing, and the lexicon decides. Overrides and decisions still win. See [Definitions](../guides/definitions.md#1-a-records-own-source_type-is-its-class).

One more column is read, but never by the lexicon: the BioSample package the record was submitted under, from the first non-empty of `package`, `attribute_package`, `biosample_package` and `epi_type`. [`resolve`](resolve.md) writes it as `package`, and NCBI Pathogen Detection tables have `attribute_package` and `epi_type`. It is used in one place only, for a specimen with no host (below).

A table with none of these columns stops with an error, because every row would go to review. The usual cause is an ID list given to `classify` instead of [`resolve`](resolve.md). With `--overrides` and an ID column, it is a warning instead, and only the rows the overrides name are classified.

### Specimens follow the host

Body sites and specimens (blood, urine, faeces, wound, tissue, catheter, "clinical") occur in any host:

- With an animal host or an animal word (`pig faeces`, urine from *Canis lupus familiaris*), the row is **Animal**. With no host, any name NCBI Taxonomy gives to an animal counts as the word: `stool (Odocoileus hemionus)`, `pigeon stool`, `giant panda feces`.
- With a human host, or the word `human` in its own text (`human feces`, `Human, blood culture`), it is **Human**. The word counts only for a specimen: `human contact surface` is not a person.
- If another field names food, an environment or a laboratory, it goes to review as `host_source_conflict`, unless that field is `source_type`, which decides (above). Words in the specimen's own text do not count: `stool from a food handler` is still Human.
- With both an animal and a human, it goes to review as `host_source_conflict`.
- With no host, and a food named in the same text (`pork blood sausage`, `oyster homogenized tissue dilution`), it goes to review as `specimen_or_food`: the specimen word is an ingredient as likely as a sample.
- With a host the lexicon does not recognise, such as a plant, it goes to review as `host_unrecognised`. Any animal NCBI Taxonomy names is recognised, not only those on the lexicon's word list; see [Hosts by taxonomy](../guides/lexicon.md#hosts-by-taxonomy).
- With no host at all, it is **Human**, and the reason ends in `+no_host` so you can count these calls. NCBI's submitters mean a person by such a record in 98.0% of 15,290 benchmark records ([Definitions](../guides/definitions.md#4-a-specimen-with-no-host-is-a-persons)).
- With no host at all in a record its submitter filed as not clinical, it goes to review as `no_host_env_package`. That is the BioSample package `Pathogen.env` ("Pathogen: environmental/food/other"), or `epi_type` `environmental/other` in a Pathogen Detection table. `feces` or `wound` in such a record is often an animal's. In a sample of 159,045 NCBI Pathogen Detection records, 586 specimens had no host and a known package, and 199 of them (34%) were filed this way. Any other package, or no package column, keeps the Human default.

### Values that are not evidence

- **Missing-value words** in `isolation_source` (`unknown`, `missing`, `not collected`, `NaN`, ...) are ignored, and the host or `source_type` decides. The full list is in [Input files](../guides/input-files.md#missing-values). `unknown` followed by a word that names something still counts, in any field: `unknown animal` is Animal, while `unknown source` and `unknown environment` mean nothing.
- **A bare scientific name** (`Escherichia coli`, `Bos taurus`) in `isolation_source` is ignored the same way: taxonomy belongs in `host`.

With nothing else to go on, the row goes to review.

!!! note "Known limitation: an animal named alone in `isolation_source`"
    Some submitters put the animal in `isolation_source` and leave `host` empty: `Larus cachinnans`, `gull`, `fox`. Such a row goes to review (`taxon_only_source` or `unmatched`), not to Animal. That is deliberate for now: with no host, a species name there is the animal in some NCBI records and its meat in others (a catfish species sold as fillets), and telling them apart needs a lexicon change measured on records nobody has looked at. On the most recent benchmark part this sent 151 of 5,697 animal records to review. Nothing is called wrongly: settle these rows in your study with one [decision](#decisions) per text, for example `gull` with `l1` Animal.

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
| `override_applied` | `yes` when an override or a decision decided |
| `lexicon_suggestion` | For an override or a decision, the lexicon's own answer; for a review row, the competing classes, if any |

### classify_reason

| Reason | Status | Meaning |
|---|---|---|
| a rule id, e.g. `hospital-drain` | ok | That lexicon rule matched |
| `<specimen rule>+<host rule>`, e.g. `specimens+human-host` | ok | A specimen, classed by its host |
| `<specimen rule>+no_host` | ok | A specimen with no host, defaulted to Human |
| `specimens+animal-name-in-specimen` | ok | A specimen with no host whose text names an animal NCBI Taxonomy knows |
| `farm-environment+<host rule>` | ok | A farm or barn, classed by its host |
| `carcass+<host rule>`, `sampling-material+<host rule>`, `animal-sites+human-host` | ok | A carcass, a bare swab or an intestine, classed by its host |
| `meat-species-no-host` | ok | A meat species (`chicken`, `turkey`, `fish`, `poultry`, `salmon`, `beef`, `lamb`) in the text, with no host: Food |
| `source_type` | ok | NCBI's `source_type` decided where the free text could not |
| `source_type_over_text` | ok | NCBI's `source_type` decided although the free text said something else |
| `override` | ok | An override for this isolate decided |
| `decision` | ok | A decision for this source text decided |
| `host_source_conflict` | review | Host and source point to different classes |
| `host_unrecognised` | review | A specimen from a host the lexicon cannot place |
| `no_host_env_package` | review | A specimen with no host, in a record filed as environmental, food or other |
| `specimen_or_food` | review | A specimen word with no host, in a text that also names a food |
| `rule_conflict:<rules>` | review | Rules of equal priority disagree |
| `animal-feed`, `fixture-no-setting` | review | A rule that sends matches to review: animal feed, or a fixture such as a drain or sink named without its setting |
| `taxon_only_source` | review | `isolation_source` is only a scientific name, and nothing else decides |
| `unknown_source` | review | `isolation_source` is a missing-value word, and nothing else decides |
| `missing_source` | review | No isolation source and no host: the record gives neither |
| `unresolved:<resolve reason>`, e.g. `unresolved:ncbi_error` | review | No isolation source and no host because [`resolve`](resolve.md#stopping-failing-and-finishing-a-run) could not get the record: NCBI gave an error, the run was stopped, or NCBI has no such record. The source is unknown, not missing |
| `unmatched` | review | No rule matched |
| `override_empty`, `override_conflict`, `override_invalid_l1`, `override_invalid_healthcare`, `override_invalid_paper_rollup` | review | A problem with the override for this row |

Four things are not a `host_source_conflict`:

- a healthcare-environment source with a human host is classified as healthcare environment, because submitters often record *Homo sapiens* as the host of hospital sink and surface swabs;
- manure, litter or bedding with an animal host is Environmental: the host names the animals it came from;
- a farm, barn or poultry shed with a host is classified by the host: the host is what was sampled there;
- a carcass, a bare `swab` or an intestine with a host is classified by the host too: `carcass swab` with host porcine is Animal, and `swab` with host *Homo sapiens* is Human.

The last three follow NCBI's records; see [Definitions](../guides/definitions.md#2-farm-settings).

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

## Decisions

An override names one isolate. A decision names a **source text**, and applies to every row that says the same thing: the same words in `isolation_source`, `host`, `sample_type` and `host_disease`. One decision for `bed rails` classifies every bed rail in the table, however it is typed (`Bed Rails`, `bed_rails`, `bed-rails`).

The decisions file is the review-groups table that `classify` writes (below), with its `l1` column filled in:

| Column | Required | Contents |
|---|---|---|
| `isolation_source`, `host`, `sample_type`, `host_disease` | at least one | The text the decision is about |
| `l1` (or `source_l1`, `classification`) | yes | The class. A row with an empty `l1` is not decided yet and is skipped. |
| `healthcare_environment` | no | `yes` or `no` (default `no`) |
| `paper_rollup` | no | Defaults to the rollup of `l1` and `healthcare_environment` |
| `environmental_package` | no | `yes` when the decision is for records filed as environmental, food or other; see below |
| `note` (or `reason`, `comment`) | no | Kept for your records |

- **A decision wins over the lexicon**, for rows the lexicon classifies as well as rows in review. A food-safety study can decide that a bare `chicken` is Food; the lexicon's answer is kept in `lexicon_suggestion`.
- **An override wins over a decision**, since it is about that one isolate.
- **The whole text must match.** `chicken` with host *Gallus gallus* is a different text from `chicken` with no host, so a decision for one does not reach the other. Missing-value words count as empty.
- **A mistake stops the run**: a class that is not a class, two different decisions for one text, or a decision with no text at all (it would apply to every row with no host and no source). A decision reaches many rows at once, so none of these is allowed to pass.
- **`environmental_package`.** A specimen with no host is classified differently in a record filed as environmental, food or other. The table `classify` writes says which of the two a group is, and a decision then applies only to that kind of record. A table you write by hand can leave the column out; its decisions are about the text alone.
- `classify` says how many rows the decisions classified, and warns about decisions that matched no row.

## Working through the review queue

`classify` writes the rows it would not guess twice:

- `<out>.review.tsv` lists every row, with its reason.
- `<out>.review_groups.tsv` lists every **distinct source text** once, with the number of rows that carry it, largest first. In a sample of 159,045 NCBI Pathogen Detection records, the 2,626 review rows that say anything hold 694 distinct texts, and the 100 most common cover 69% of them. This is the shorter table to work through.

| Column | Contents |
|---|---|
| `rows` | How many review rows say this |
| `classify_reason`, `lexicon_suggestion` | Why they are in review, and the classes in competition, if any |
| `isolation_source`, `host`, `sample_type`, `host_disease` | The text, as first typed. A missing-value word is shown as empty. |
| `environmental_package` | `yes` for records filed as environmental, food or other |
| `example_id` | One isolate with this text |
| `l1`, `healthcare_environment`, `note` | Empty, for your decision |

1. Start at the top of the review-groups table: the first lines cover the most rows. Rows with no text at all (`missing_source`) have nothing to decide from; check the BioSample or your lab records for those.
2. Fill in `l1` for the texts you can place, and `healthcare_environment` where it is `yes`. Leave the rest empty.
3. Run `classify` again with `--decisions` pointing at the table you filled in. The same file can be kept and extended from one dataset to the next.
4. For an isolate that is an exception to its text, add a row to an [overrides](#overrides) file.
5. If a word is general, add it to a copy of the lexicon and pass `--lexicon`; see [The lexicon](../guides/lexicon.md).
6. Measure the effect with [`validate`](validate.md) before relying on decisions or a changed lexicon.
