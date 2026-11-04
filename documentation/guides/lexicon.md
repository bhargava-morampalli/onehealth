# The lexicon

`classify` decides each class with a lexicon: a YAML file of rules that map words and phrases in the metadata to a One Health class. The rules are data, not code, so you can read every one, and change them without touching `onehealth`.

## Versions

| Version | Shipped | Use |
|---|---|---|
| v3 | November 2026 | The default |
| v2 | October 2026 | `--lexicon v2`, to reproduce earlier results exactly |

Every classified row records the version that made the call in `lexicon_version`. What changed between versions is in the [changelog](../changelog.md).

## A rule

```yaml
- id: hospital-drain
  fields: [isolation_source, sample_type]
  phrases:
    - hospital drain
    - sink drain
    - sink trap
  l1: Environmental
  healthcare_environment: "yes"
  priority: 95
```

This is part of a rule in v3, shortened.

| Key | Required | Meaning |
|---|---|---|
| `id` | yes | The rule's name. It is written to `classify_reason` when the rule decides. |
| `l1` | yes | `Human`, `Animal`, `Environmental`, `Food`, `Laboratory`, or `Unclassified` to send matches to review. |
| `fields` | no | Which fields to read: `isolation_source`, `host`, `sample_type`, `host_disease` (default `[isolation_source]`). See [Input files](input-files.md#columns-classify-reads) for the columns behind each field. |
| `phrases` | no | Words or phrases that may appear anywhere in the field. |
| `exact` | no | Values that must be the whole field. |
| `healthcare_environment` | no | `"yes"` for a healthcare setting, or `"no"` (the default). |
| `priority` | no | A whole number; higher wins (default 0). |
| `host_dependent` | no | `true` for a body site or specimen that occurs in any host; see below. |
| `note` | no | Why the rule exists, for the next reader. |

A rule with a missing `id` or `l1`, or an unknown class or field, stops with an error naming the rule. So does a phrase YAML reads as true or false, such as an unquoted `yes`, `no`, `on` or `off`: put it in quotes.

## How rules match

- **By whole word.** Values are lower-cased, underscores read as spaces, and split into words at anything that is not a letter or digit. `drain` matches `Sink drain 3` but not `drainage`, and `sink drains` needs its own phrase.
- **`phrases` anywhere, `exact` whole.** `environment` as an exact value matches the value `Environment` but not `hospital environment`.
- **One hit per rule:** the first of its fields that matches.

## How a class is chosen

1. **Missing values** and **bare scientific names** in `isolation_source` are dropped; see [Missing values](input-files.md#missing-values).
2. **The highest priority wins.** Rules of equal priority that disagree about the class or the healthcare flag send the row to review as `rule_conflict:<rules>`.
3. **A rule with `l1: Unclassified`** sends the row to review with its id as the reason, such as `animal-feed`.
4. **A specimen rule** (`host_dependent: true`) lets the host decide: an animal host or animal word makes it Animal, a human host makes it Human, both is `host_source_conflict`, a host no rule recognises is `host_unrecognised`, and no host at all gives the rule's own `l1` with the reason `<rule>+no_host`.
5. **Any other rule** decides, unless the host points to a different class than the source (`host_source_conflict`). The exception is a healthcare-environment source with a human host, which stays a healthcare environment: submitters often record *Homo sapiens* as the host of hospital sink and surface swabs.
6. **No rule** sends the row to review as `missing_source`, `unknown_source`, `taxon_only_source` or `unmatched`; see [`classify`](../commands/classify.md#classify_reason).

The paper rollup is not part of the lexicon. It follows from `l1` and `healthcare_environment`, as in [`classify`](../commands/classify.md#the-classes).

## The rules in v3

| Priority | Rule | Class | Reads |
|---|---|---|---|
| 99 | `animal-feed` | review | isolation source, sample type |
| 95 | `wound-drainage` | specimen (host decides) | isolation source, sample type |
| 95 | `hospital-drain` | Environmental, healthcare | isolation source, sample type |
| 92 | `clinical-environment` | Environmental, healthcare | isolation source, sample type |
| 88 | `animal-food-producing` | Animal | isolation source, sample type, host |
| 88 | `animal-mastitis` | Animal | isolation source, sample type, host disease |
| 85 | `healthcare-fixtures` | Environmental, healthcare | isolation source, sample type |
| 80 | `clinical-sample` | specimen (host decides) | isolation source |
| 80 | `specimens` | specimen (host decides) | isolation source, sample type |
| 75 | `food-products` | Food | isolation source, sample type, host |
| 70 | `laboratory` | Laboratory | isolation source, sample type, host |
| 70 | `environment-outdoor` | Environmental | isolation source, sample type |
| 65 | `animal-sites` | Animal | isolation source, sample type |
| 60 | `environment-generic` | Environmental | isolation source, sample type |
| 40 | `animal-hosts` | Animal | host, isolation source, sample type |
| 20 | `human-host` | Human | host, sample type |
| 20 | `human-source` | Human | isolation source |

The phrases of each rule, and a note on why it is there, are in [`src/onehealth/data/lexicon/v3.yaml`](https://github.com/bhargava-morampalli/onehealth/blob/main/src/onehealth/data/lexicon/v3.yaml). Phrases new in v3 are marked `v3`.

Some values are in no rule on purpose. A bare `sink`, `shower` or `faucet` says nothing about the setting, and a wrong healthcare flag would quietly move isolates into the report's healthcare-environment group. Such a row is decided by its host, or goes to review without one.

## Changing the lexicon

Add words to a copy, not to the shipped file:

1. Copy the shipped lexicon. In a clone of the repository it is `src/onehealth/data/lexicon/v3.yaml`.
2. Give the copy its own `version`, such as `"3-mylab"`, so its results are not mistaken for the shipped v3. A file without a `version` is recorded as the default version.
3. Add phrases to a rule, or a rule of your own. For example, a study that sampled only hospital sinks can classify a bare `sink`:

    ```yaml
    - id: study-sinks
      fields: [isolation_source, sample_type]
      phrases:
        - sink
      l1: Environmental
      healthcare_environment: "yes"
      priority: 85
      note: All sinks in this study were in hospital wards.
    ```

4. Run `classify --lexicon my-lexicon.yaml`, then [`validate`](../commands/validate.md) against labels you trust, and compare with the shipped lexicon on the same labels.

For a word that only one study uses, an [overrides file](../commands/classify.md#overrides) is often simpler than a rule.
