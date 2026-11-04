# validate

Compares the classes from [`classify`](classify.md) with labels you trust, such as a hand-curated subset of your own data. Run it before relying on `onehealth` for a new organism, source or lexicon.

```bash
onehealth validate work/classified.tsv --truth curated.tsv -o work/validation.tsv
```

## Options

| Option | What it does |
|---|---|
| `--truth` | A TSV with an isolate ID column and a label column (required). |
| `--truth-column` | The label column in `--truth` (default `truth`). |
| `--level` | Which class to compare: `paper_rollup` (the default) or `l1`. |
| `-o`, `--out` | Summary TSV (required). The other tables are written next to it. |

## The truth file

- **ID column:** the first of `query_id`, `id`, `accession`, `isolate_id`, `sample_id`, `biosample`, `biosample_acc`, `assembly` or `target_acc`. Each ID is matched against every ID on a classified row, so a BioSample label finds a row that was resolved from an assembly. An assembly also matches its GCA/GCF twin and its unversioned form.
- **Labels** are class names of the chosen `--level`, in any case. `Environment` is read as `Environmental`. `Unclassified` is allowed, for isolates the curators could not place either.
- Rows with an empty ID or label are skipped.

`validate` stops with an error rather than guess when:

- a label is not a class of `--level`. Map your labels first, for example `clinical` to `Human`;
- one isolate has two different labels, in one row each or under two of its IDs. Keeping either label would bias the accuracy without anyone noticing.

## The three numbers

| Metric | Meaning |
|---|---|
| coverage | Of isolates with a real class in the truth, the share `onehealth` classified. The rest went to the review queue. |
| accuracy | Of the isolates `onehealth` classified, the share whose class matches the truth. |
| agreement | Of all isolates, the share where `onehealth` and the truth agree, counting `Unclassified` on both sides as agreement. |

Coverage and accuracy are the pair to watch. Review is the safe answer for unclear evidence, so a low coverage with a high accuracy means "the metadata is thin", not "the classifier is wrong".

## Output

`validate` prints the three numbers and the confusion matrix, and writes four tables. For `-o work/validation.tsv`:

| File | Contents |
|---|---|
| `validation.tsv` | `isolates`, `coverage`, `accuracy`, `agreement` and `unmatched_truth_ids` (truth IDs not found in the classified table), each with its meaning |
| `validation.by_class.tsv` | Per class: isolates in the truth, isolates predicted, correct, recall and precision |
| `validation.confusion.tsv` | Truth classes as rows, `onehealth` classes as columns |
| `validation.disagreements.tsv` | Every isolate where they differ, with its host, isolation source, sample type and `classify_reason` |

Recall is the share of isolates truly in a class that `onehealth` put there; precision is the share of isolates `onehealth` put in a class that are truly in it.

## Using the disagreements

Sort `validation.disagreements.tsv` by `classify_reason`:

- **Predicted `Unclassified`:** the row went to review. Most often the metadata has no isolation source or host. If a word the lexicon should know is missing, see [The lexicon](../guides/lexicon.md).
- **Predicted a different class:** look at the host and isolation source. Either the lexicon rule is too broad, or the truth label used information the metadata does not have.

After changing the lexicon, run `classify` and `validate` again on the same truth set. [Benchmarks](../benchmarks.md) shows the results on two public datasets.
