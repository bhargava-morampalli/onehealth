# The lexicon

`classify` decides each class with a lexicon: a YAML file of rules that map words and phrases in the metadata to a One Health class. The rules are data, not code, so you can read every one, and change them without touching `onehealth`.

## Version

`onehealth` ships one lexicon, [`src/onehealth/data/lexicon.yaml`](https://github.com/bhargava-morampalli/onehealth/blob/main/src/onehealth/data/lexicon.yaml). It is version **8**, and every classified row records the version that made the call in `lexicon_version`, so a table always says which rules it came from.

When the lexicon changes, its version goes up, and the [changelog](../changelog.md) says what changed and what it did to the [benchmarks](../benchmarks.md). Earlier versions are not shipped. They are in the repository history, so results made with one can be reproduced by running the `onehealth` commit that made them, or by passing that version's YAML file to `--lexicon`.

A lexicon of your own is recorded under the `version` it carries. A file with no `version` is recorded as `custom-` and the first 8 hex digits of the file's SHA-256, and an edited copy that still says version 7 is recorded as, say, `8-modified-1a2b3c4d`. So a table never claims the shipped lexicon when something else made it.

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

This is part of a rule, shortened.

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
| `taxonomy` | no | `animal` or `human`: the rule also matches a field whose value is a name NCBI Taxonomy gives to an animal, or to people; see [Hosts by taxonomy](#hosts-by-taxonomy). |
| `within` | no | With `taxonomy: animal`: `true` when the name may be anywhere in the field. For a rule that reads source text (`stool (Odocoileus hemionus)`). |
| `specimen_only` | no | `true` for a rule that only says whose specimen it is (`human`). It decides Human or Animal for a specimen rule on the same row, and is ignored otherwise. |
| `no_host` | no | `true` for a rule that applies only to a row with no host (`chicken` with no host is meat). |
| `setting_of` | no | `animal` for what animals leave behind (manure): an animal host is not a contradiction, and the row stays the rule's class. |
| `yields_to_host` | no | `true` for a place (a farm) or a thing that says nothing of whose it is (a carcass, a swab): with a human or animal host, the host decides. |
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
4. **A specimen rule** (`host_dependent: true`) lets the host decide: an animal host or animal word makes it Animal, a human host makes it Human, both is `host_source_conflict`, another field that names food, an environment or a laboratory is `host_source_conflict` too, a host no rule recognises is `host_unrecognised`, and no host at all gives the rule's own `l1` with the reason `<rule>+no_host`, unless the same text names a food (`specimen_or_food`). The one exception to that default is a record its submitter filed as environmental, food or other and not as clinical: it goes to review as `no_host_env_package`; see [`classify`](../commands/classify.md#specimens-follow-the-host).
5. **Any other rule** decides, unless the host points to a different class than the source (`host_source_conflict`). Three exceptions: a healthcare-environment source with a human host stays a healthcare environment, since submitters often record *Homo sapiens* as the host of hospital sink and surface swabs; a `setting_of: animal` rule (manure) with an animal host keeps its class; and a `yields_to_host` rule (a farm, a carcass, a bare swab) takes the class of the host.
6. **No rule** sends the row to review as `missing_source`, `unknown_source`, `taxon_only_source` or `unmatched`; see [`classify`](../commands/classify.md#classify_reason).
7. **Above all of this, a record's own `source_type` decides** when it names a class. That is not a lexicon rule but part of `classify`; see [Definitions](definitions.md#1-a-records-own-source_type-is-its-class).

The paper rollup is not part of the lexicon. It follows from `l1` and `healthcare_environment`, as in [`classify`](../commands/classify.md#the-classes).

## The rules in v8

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
| 79 | `animal-sites` | Animal; with a human host, Human | isolation source, sample type |
| 79 | `animal-alive` | Animal | isolation source, sample type |
| 78 | `environmental-sampling` | Environmental | isolation source, sample type |
| 78 | `food-handler` | Human | isolation source, sample type |
| 77 | `farmed-food` | Food | isolation source, sample type |
| 76 | `farm-environment` | Environmental; with a host, the host's class | isolation source, sample type |
| 75 | `food-products` | Food | isolation source, sample type, host |
| 75 | `carcass` | Food; with a host, the host's class | isolation source, sample type |
| 75 | `food-items` | Food | isolation source, sample type |
| 70 | `laboratory` | Laboratory | isolation source, sample type, host |
| 70 | `environment-outdoor` | Environmental | isolation source, sample type |
| 65 | `animal-waste` | Environmental | isolation source, sample type |
| 60 | `environment-generic` | Environmental | isolation source, sample type |
| 50 | `sampling-material` | Environmental; with a host, the host's class | isolation source, sample type |
| 45 | `meat-species-no-host` | Food, with no host only | isolation source |
| 40 | `animal-hosts` | Animal | host, isolation source, sample type |
| 40 | `host-taxonomy` | Animal, by NCBI Taxonomy | host |
| 30 | `fixture-no-setting` | review | isolation source, sample type |
| 20 | `human-host` | Human | host, sample type |
| 20 | `host-taxonomy-human` | Human, by NCBI Taxonomy | host |
| 20 | `human-word` | Human, for a specimen only | isolation source |
| 20 | `animal-name-in-specimen` | Animal, for a specimen with no host only, by NCBI Taxonomy | isolation source |
| 20 | `human-source` | Human | isolation source |

The phrases of each rule, and a note on why it is there, are in [`src/onehealth/data/lexicon.yaml`](https://github.com/bhargava-morampalli/onehealth/blob/main/src/onehealth/data/lexicon.yaml). What each version added is marked `v3` to `v8` there.

### Hosts by taxonomy

`animal-hosts` is a list of the hosts seen most often (`chicken`, `bos taurus`, `dog`). No list covers every species a pathogen is isolated from, so a gull, a macaque or a water buffalo used to go to review as a host the lexicon could not place. `host-taxonomy` recognises any animal host by the names NCBI Taxonomy gives it:

- **the scientific name of a chordate** (a vertebrate or one of its relatives) at any rank: `Larus dominicanus`, `Macaca`, `Cervidae`, `Aves`;
- **the species name of a mollusc, a crustacean, or any other animal that also has a common name** (v8): `Magallana gigas` (an oyster), `Litopenaeus vannamei` (a shrimp), `Musca domestica`;
- **a species with a subspecies or breed**, or a genus with `sp.`: `Bos primigenius taurus`, `Meleagris sp.`;
- **the common name of any animal**, singular or plural: `gull`, `sika deer`, `rhesus macaque`, `mud crab`.

The host must be a name, or end with one (v8): `European rook`, `free-range white-tailed deer` and `white stork nestling` are recognised, because English puts the noun last and a word for the animal's age may follow it. A name in the middle of a value is not used, nor one in a list (`bologna, hot dogs, sausage-like products`), so `sea water` is not the fish genus *Sea*. A genus on its own counts only for chordates: other animals have genera called *Sepsis*, *Cancer* and *Colon*, and a host field can say `sepsis` or `cancer`. Plants, fungi and bacteria are not placed.

`host-taxonomy-human` does the same for the genus *Homo* however it is written (`Homo`, `Homo sapiens sapiens`). Both rules come after the word-list rule of the same priority, so a host that rule knows keeps its reason (`animal-hosts`, `human-host`).

**In a specimen's own text** (v8). `animal-name-in-specimen` reads `isolation_source` with `within: true`: `stool (Odocoileus hemionus)`, `pigeon stool` and `giant panda feces` are Animal, where the no-host default made them Human. It is stricter than the host rule, because ordinary words are animals' names too. A name of two or more words always counts. A single word counts only when it is the common name of a vertebrate and not also an everyday word: `pigeon` and `deer` count; `tip`, `white` and `sponge` (a moth, a butterfly and an invertebrate) do not, nor do `sole`, `ray` or `monitor`. The rule is `specimen_only` and `no_host`: it decides whose specimen it is and nothing else, so a species name alone in `isolation_source` is still not a description of a sample, and goes to review as `taxon_only_source`.

The names are in `src/onehealth/data/animal_names.tsv.gz`: about 204,000 scientific names and 36,000 common names from an NCBI Taxonomy dump, first of 14 November 2026 and refreshed since (below). `classify` prints the date of the dump after the lexicon version (`lexicon v8 (seed v8), animal names from NCBI Taxonomy of 2026-11-14`), so a methods section can name both. A lexicon of your own gets the same behaviour by giving a rule `taxonomy: animal`, and loses it by leaving the rule out.

#### Keeping the names current

NCBI adds and renames species every day, so the table falls behind. A newer table can change classes: a host that names a newly listed animal is placed, where it went to review before. So the table is never replaced without review, and `classify` prints the date of the table it used.

- **Once a month**, the workflow `.github/workflows/taxonomy.yml` rebuilds the table from the newest NCBI Taxonomy with `tools/refresh_animal_names.py`. If no name changed, it stops there.
- **If names changed**, it runs the tests with the new table and opens a pull request. The pull request lists the names added and removed, every record of the benchmark datasets whose class changes (and whether the change agrees with NCBI's label), and what to do before merging: check the changes, update counts that tests quote, and raise the version if any class changes. The pull request adds its own section to the changelog. When the tests pass and no class changes, merging it is all there is to do.
- **By hand**, `python tools/refresh_animal_names.py` does the same on your machine, and `python tools/build_animal_names.py` rebuilds the table with no report. Both take a dump you already have (`new_taxdump.tar.gz`) instead of downloading one. The same dump always gives the same file, byte for byte.

The pull request needs the repository setting that lets GitHub Actions create pull requests (Settings, Actions, General, Workflow permissions). Without it, the workflow pushes the branch and says so.

### Whose specimen

A specimen follows its host, and an animal word in the specimen's own text counts as a host: `pig feces` is Animal. `human-word` does the same for people: `human feces` and `Human, blood culture` are Human, with the reason `specimens+human-word`. Before v6 they were Human only by the no-host default, and went to review in a record filed as environmental, food or other.

The rule is `specimen_only`, because the word decides nothing on its own: `human contact surface` is a surface and `human listeriosis` names no sample. `patient` is not listed, since a veterinary record has patients too.

**Farm settings and meat species follow NCBI.** `animal-waste`, `farm-environment` and `meat-species-no-host` came in v7 and change what manure, a farm and `chicken` with no host mean; v8 widened the last from the whole values `chicken` and `turkey` to `chicken`, `turkey`, `poultry`, `fish`, `salmon`, `beef` and `lamb` anywhere in the text. [Definitions](definitions.md) has the numbers from NCBI's records behind each.

**Meat is food, the animal needs a host.** How meat is sold or processed means food, whatever the species: `ground turkey`, `chicken breast`, `comminuted chicken`, `raw intact chicken`, `smoked salmon`. Before v5 all of these were Animal. `food-items` lists such words and the foods that NCBI records name most often; because it reads source text only, two-word foods such as `peanut butter` are not taken for species names. With no host, a meat species in the text is food too (`meat-species-no-host`): `fish`, `poultry`, `chicken liver`. What still means the animal:

- **A host**: `chicken liver` with host *Gallus gallus* is a necropsy sample.
- **A gut site, alive, or a herd** (`animal-sites`, `animal-alive`): `chicken cecal content`, `chicken gut`, `live chicken`, `beef cattle`.
- **A word without the evidence**: `duck`, `pig`, `cow`, `cattle`, `swine`, `bovine`.

A bare body part (`breast`, `liver`, `thigh`) is as often a patient's, so cuts are listed only with their animal (`chicken breast`). `egg`, `seafood` and `raw` are not listed: NCBI records use them both ways.

**A carcass is meat, or the animal it was** (v8). `carcass rinse` with no host is Food. `carcass swab` with host porcine is Animal, as in all 211 NCBI records with a carcass and an animal host; until v8 it went to review.

A cut of meat or milk with an animal host (`ground chicken`, host *Gallus gallus*; `milk`, host bovine) goes to review as `host_source_conflict`: NCBI's records are divided.

**A sponge, a wipe, dust or a bare `swab`** (v8, `sampling-material`) is a sample of a surface: Environmental, below the food rules, so `turkey carcass sponge` is still the carcass. A swab of a host is the host: `swab` with host *Homo sapiens* is Human.

**A swab of the premises is not the food.** `environmental-sampling` makes a swab or sponge of a surface, and a sample of the place where food is made, Environmental even when the text names the food or the animal: `environmental swab sponge non food contact surface`, `food processing environment`, `poultry processing plant`. Before v5 these were Food or Animal. It sits above the food rules and below the specimen rules. `food-handler` is a person, so `food handler` is Human, not Food.

A fixture on its own, such as `Rm 14 Drain`, `Sink swab` or `Shower 3`, says neither whose sample it is nor where the fixture is. Calling it a hospital fixture could be wrong, and a wrong healthcare flag would quietly move isolates into the report's healthcare-environment group. Calling it Human because the host is *Homo sapiens* is often wrong too, since submitters record that host for hospital fixture swabs. So `fixture-no-setting` sends it to review. Any rule above it still decides: `sink drain` and `hospital sink` are healthcare environment, `drain tip` is a specimen, and an animal host makes it Animal.

## Changing the lexicon

Add words to a copy, not to the shipped file:

1. Copy the shipped lexicon. In a clone of the repository it is `src/onehealth/data/lexicon.yaml`.
2. Give the copy its own `version`, such as `"7-mylab"`, so its results are easy to tell apart. If you forget, they are still not mistaken for the shipped lexicon: an edited copy that keeps the version `7` is recorded as `7-modified-<digest>`, and a file without a `version` as `custom-<digest>`.
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
