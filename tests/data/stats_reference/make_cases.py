"""Write cases.tsv: the inputs on which onehealth's statistics are compared with R.

    python tests/data/stats_reference/make_cases.py

Most cases are counts taken from the real NCBI Klebsiella fixture
(tests/data/ncbi_pd_klebsiella): isolates per lineage in each source, isolates with
each resistance gene in each source. A few are written by hand to cover the edges
(one isolate, every lineage equally common, a gene in every isolate).

The file is fixed once written. reference.R reads it and writes reference.tsv, and
tests/test_stats_reference.py feeds the same inputs to onehealth and compares. So the
comparison does not depend on how the classifier places a record today.

One line per case: kind, name, values.

* counts       isolates per lineage in one source: "5,3,2"
* proportion   isolates with a feature, isolates genotyped: "7,40"
* pvalues      p-values to adjust; NA is a test that was not run: "0.01,NA,0.2"
* features     one dataset: isolates per source, then for each feature the isolates
               with it in each source: "sizes=120,60,20;geneA=30,2,1;geneB=..."
* lineages     one dataset: for each lineage the isolates in each source:
               "ST1=30,2,1;ST2=..."
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from onehealth import analysis as A  # noqa: E402
from onehealth.classify import classify_rows  # noqa: E402
from onehealth.io_tsv import read_features, read_tsv  # noqa: E402
from onehealth.lexicon import load_lexicon  # noqa: E402

DATA = ROOT / "tests" / "data" / "ncbi_pd_klebsiella"
OUT = Path(__file__).resolve().parent / "cases.tsv"


def main() -> int:
    classified = classify_rows(read_tsv(DATA / "metadata.tsv"), load_lexicon())
    dataset = A.build_dataset(classified, lineage_column="PDS_acc")
    typed = [iso for iso in dataset.classified if iso.lineage]
    groups = [g for g in dataset.groups if any(iso.source == g for iso in typed)]
    lines = ["# Inputs for the comparison with R. Made by make_cases.py; see README.md.",
             f"# Sources, in the order of every per-source list below: {', '.join(groups)}",
             "kind\tcase\tvalues"]

    def add(kind: str, case: str, values: str) -> None:
        lines.append(f"{kind}\t{case}\t{values}")

    # Lineages per source: Simpson's index, its interval and rarefaction.
    for group in groups:
        counts = sorted(Counter(iso.lineage for iso in typed if iso.source == group).values(), reverse=True)
        add("counts", f"klebsiella_{group.lower().replace(' ', '_')}", ",".join(map(str, counts)))
    add("counts", "two_singletons", "1,1")
    add("counts", "one_lineage", "12")
    add("counts", "ten_singletons", ",".join(["1"] * 10))
    add("counts", "four_equal", "5,5,5,5")
    add("counts", "one_dominant", "96,1,1,1,1")
    add("counts", "small_uneven", "5,3,2")
    add("counts", "large", "400,300,150,80,40,20,5,3,1,1")

    # Lineage by source: enrichment.
    table: dict[str, Counter] = {}
    for iso in typed:
        table.setdefault(iso.lineage, Counter())[iso.source] += 1
    ordered = sorted(table, key=lambda lin: (-sum(table[lin].values()), lin))
    add("lineages", "klebsiella", ";".join(f"{lin}={','.join(str(table[lin][g]) for g in groups)}" for lin in ordered))
    add("lineages", "restricted", "A=20,0,0;B=0,20,0;C=5,5,10;D=2,1,1;E=10,10,10")

    # Features per source: prevalence, its interval and the test of independence.
    presence, _ = read_features(DATA / "features.tsv")
    mapping, _ = A.match_names(list(presence), dataset)
    by_key: dict[str, set[str]] = {}
    for name, key in mapping.items():
        by_key.setdefault(key, set()).update(presence[name])
    isolates = [iso for iso in dataset.classified if iso.key in by_key]
    fgroups = [g for g in dataset.groups if any(iso.source == g for iso in isolates)]
    assert fgroups == groups, (fgroups, groups)
    sizes = [sum(iso.source == g for iso in isolates) for g in groups]
    counts: dict[str, list[int]] = {}
    for iso in isolates:
        for feature in by_key[iso.key]:
            counts.setdefault(feature, [0] * len(groups))[groups.index(iso.source)] += 1
    total = sum(sizes)
    testable = sorted((f for f, c in counts.items() if 0 < sum(c) < total), key=lambda f: (-sum(counts[f]), f))
    # Every fourth gene, so common and rare ones are both there, and the file stays small.
    chosen = testable[::4]
    add("features", "klebsiella",
        f"sizes={','.join(map(str, sizes))};" + ";".join(f"{f}={','.join(map(str, counts[f]))}" for f in chosen))
    add("features", "two_sources", "sizes=30,12;all_in_one=10,0;even=10,4;rare=1,0;nearly_all=29,12;opposite=2,11")
    for feature in chosen[:12]:
        for g, x, n in zip(groups, counts[feature], sizes):
            add("proportion", f"{feature}|{g}", f"{x},{n}")
    for case, values in (("none_of_40", "0,40"), ("all_of_40", "40,40"), ("one_of_1", "1,1"), ("none_of_1", "0,1"),
                         ("one_of_3", "1,3"), ("half_of_1000", "500,1000"), ("three_of_10000", "3,10000")):
        add("proportion", case, values)

    add("pvalues", "mixed", "0.0001,0.0004,0.0019,0.0095,0.0201,0.0278,0.0298,0.0344,0.0459,0.3240,0.4262,0.5719,0.6528,0.7590,1")
    add("pvalues", "ties_and_ones", "0.05,0.05,0.05,1,1,0.0001,0.5")
    add("pvalues", "with_untested", "0.01,NA,0.04,0.03,NA,0.5")
    add("pvalues", "single", "0.03")
    add("pvalues", "unsorted", "0.9,0.0002,0.6,0.0450,0.0100,0.2")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(lines) - 3} cases; sources {groups}; {len(chosen)} features, {len(ordered)} lineages", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
