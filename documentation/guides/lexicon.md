# The lexicon

`classify` decides each class with a lexicon: a YAML file of rules that map words and phrases in the metadata to a One Health class. The rules are data, not code, so you can read every one, and change them without touching `onehealth`.

## Versions

| Version | Shipped | Use |
|---|---|---|
| v5 | November 2026 | The default |
| v4 | November 2026 | `--lexicon v4`, to reproduce earlier results exactly |
| v3 | November 2026 | `--lexicon v3`, the same for older results |
| v2 | October 2026 | `--lexicon v2`, the same for older results |

Every classified row records the version that made the call in `lexicon_version`. What changed between versions is in the [changelog](../changelog.md).

A lexicon of your own is recorded under the `version` it carries. A file with no `version` is recorded as `custom-` and the first 8 hex digits of the file's SHA-256, and an edited copy that still carries a shipped number is recorded as, say, `5-modified-1a2b3c4d`. So a table never claims a shipped lexicon that did not make it.

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

This is part of a rule in v5, shortened.

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
3. **A rule with `l1: Unclassified`** sends the row to review with its id as the reason, such as `animal-feed` or `fixture-no-setting`. Such a rule names no class, so it only decides when it is the best match; it never counts as a conflict with the host.
4. **A specimen rule** (`host_dependent: true`) lets the host decide: an animal host or animal word makes it Animal, a human host makes it Human, both is `host_source_conflict`, another field that names food, an environment or a laboratory (such as `source_type`) is `host_source_conflict` too, a host no rule recognises is `host_unrecognised`, and no host at all gives the rule's own `l1` with the reason `<rule>+no_host`. The one exception to that default is a record its submitter filed as environmental, food or other and not as clinical: it goes to review as `no_host_env_package`; see [`classify`](../commands/classify.md#specimens-follow-the-host).
5. **Any other rule** decides, unless the host points to a different class than the source (`host_source_conflict`). The exception is a healthcare-environment source with a human host, which stays a healthcare environment: submitters often record *Homo sapiens* as the host of hospital sink and surface swabs.
6. **No rule** sends the row to review as `missing_source`, `unknown_source`, `taxon_only_source` or `unmatched`; see [`classify`](../commands/classify.md#classify_reason).

The paper rollup is not part of the lexicon. It follows from `l1` and `healthcare_environment`, as in [`classify`](../commands/classify.md#the-classes).

## The rules in v5

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
| 78 | `environmental-sampling` | Environmental | isolation source, sample type |
| 78 | `food-handler` | Human | isolation source, sample type |
| 75 | `food-products` | Food | isolation source, sample type, host |
| 75 | `food-items` | Food | isolation source, sample type |
| 70 | `laboratory` | Laboratory | isolation source, sample type, host |
| 70 | `environment-outdoor` | Environmental | isolation source, sample type |
| 65 | `animal-sites` | Animal | isolation source, sample type |
| 60 | `environment-generic` | Environmental | isolation source, sample type |
| 40 | `animal-hosts` | Animal | host, isolation source, sample type |
| 30 | `fixture-no-setting` | review | isolation source, sample type |
| 20 | `human-host` | Human | host, sample type |
| 20 | `human-source` | Human | isolation source |

The phrases of each rule, and a note on why it is there, are in [`src/onehealth/data/lexicon/v5.yaml`](https://github.com/bhargava-morampalli/onehealth/blob/main/src/onehealth/data/lexicon/v5.yaml). Changes in v3, v4 and v5 are marked `v3`, `v4` and `v5`.

**Meat is food, the animal is not.** A species word on its own (`chicken`, `turkey`, `pig`) means the animal. How meat is sold or processed means food, whatever the species: `ground turkey`, `chicken breast`, `comminuted chicken`, `raw intact chicken`, `smoked salmon`. v4 called all of these Animal. `food-items` lists such words and the foods that NCBI records name most often; because it reads source text only, two-word foods such as `peanut butter` are not taken for species names. Three kinds of word are left out on purpose:

- **A bare body part** (`breast`, `liver`, `thigh`) is as often a patient's, so cuts are listed with their animal (`chicken breast`).
- **An organ** (`chicken liver`) is retail giblets or a necropsy sample, so it stays Animal unless the text says `giblets`.
- **Words that NCBI records use both ways:** bare `beef` (a quarter are beef cattle), `egg`, `seafood`, `fish`, `frozen` and `raw`.

A cut of meat with an animal host (`ground chicken`, host *Gallus gallus*) goes to review as `host_source_conflict`, as `chicken meat` with that host always did.

**A swab of the premises is not the food.** `environmental-sampling` makes a swab or sponge of a surface, and a sample of the place where food is made, Environmental even when the text names the food or the animal: `environmental swab sponge non food contact surface`, `food processing environment`, `poultry processing plant`. v4 called these Food or Animal. It sits above the food rules and below the specimen rules. `food-handler` is a person, so `food handler` is Human, not Food.

A fixture on its own, such as `Rm 14 Drain`, `Sink swab` or `Shower 3`, says neither whose sample it is nor where the fixture is. Calling it a hospital fixture could be wrong, and a wrong healthcare flag would quietly move isolates into the report's healthcare-environment group. Calling it Human because the host is *Homo sapiens* is often wrong too, since submitters record that host for hospital fixture swabs. So `fixture-no-setting` sends it to review. Any rule above it still decides: `sink drain` and `hospital sink` are healthcare environment, `drain tip` is a specimen, and an animal host makes it Animal, as other farm settings do.

## Changing the lexicon

Add words to a copy, not to the shipped file:

1. Copy the shipped lexicon. In a clone of the repository it is `src/onehealth/data/lexicon/v5.yaml`.
2. Give the copy its own `version`, such as `"5-mylab"`, so its results are easy to tell apart. If you forget, they are still not mistaken for the shipped v5: an edited copy that keeps the number `5` is recorded as `5-modified-<digest>`, and a file without a `version` as `custom-<digest>`.
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
