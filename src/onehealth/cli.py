"""onehealth resolve | classify | summarize | plot | validate | report"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from onehealth import __version__
from onehealth.classify import classify_fieldnames, classify_rows
from onehealth.constants import (
    HOST_ALIASES,
    HOST_DISEASE_ALIASES,
    ISOLATION_SOURCE_ALIASES,
    RESOLVE_FIELDS,
    SAMPLE_TYPE_ALIASES,
    SEED_LEXICON_VERSION,
)
from onehealth.errors import InputError
from onehealth.io_tsv import (
    ID_LIST_ALIASES,
    ISOLATE_NAME_ALIASES,
    first_present,
    read_id_list,
    read_tsv,
    require_columns,
    write_tsv,
)
from onehealth.lexicon import load_lexicon
from onehealth.ncbi import NCBIClient
from onehealth.overrides import load_overrides, unmatched_overrides
from onehealth.plot import plot_counts
from onehealth.resolve import load_local_metadata, resolve_ids
from onehealth.summarize import summarize, summary_fieldnames


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        return 2
    try:
        return args.handler(args)
    except (InputError, OSError) as exc:
        # A problem with the user's files or options: one line, no traceback. Anything else is a bug and keeps it.
        print(f"onehealth {args.command}: error: {_describe(exc)}", file=sys.stderr)
        return 1


def _describe(exc: Exception) -> str:
    if isinstance(exc, OSError) and exc.strerror and exc.filename:
        return f"{exc.strerror}: {exc.filename}"
    return str(exc)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="onehealth",
        description=(
            "Resolve mixed NCBI isolate IDs and classify isolation sources "
            "into One Health L1 classes. Missing, empty, or conflicting "
            "evidence goes to a review queue."
        ),
    )
    parser.add_argument("--version", action="version", version=f"onehealth {__version__}")
    sub = parser.add_subparsers(dest="command")

    resolve_p = sub.add_parser("resolve", help="Resolve GCA/GCF, BioSample, SRA, or PDT IDs")
    resolve_p.add_argument("ids", type=Path, help="One ID per line, or a TSV with an id column")
    resolve_p.add_argument("-o", "--out", type=Path, required=True, help="Resolved metadata TSV")
    resolve_p.add_argument("--metadata", type=Path, help="Local metadata TSV to join (no download)")
    resolve_p.add_argument("--no-network", action="store_true", help="Do not call NCBI")
    resolve_p.add_argument("--email", default="", help="Contact email for NCBI E-utilities")
    resolve_p.add_argument("--api-key", default="", help="Optional NCBI API key")
    resolve_p.add_argument("--review", type=Path, help="Review-queue TSV (default: <out>.review.tsv)")
    resolve_p.add_argument(
        "--pd-organism",
        type=_names,
        default=(),
        metavar="GROUP",
        help="Pathogen Detection organism group(s) of your PDT IDs, comma-separated (e.g. Klebsiella, "
        "Escherichia_coli_Shigella). Their current metadata table is read from NCBI to find each PDT's BioSample.",
    )
    resolve_p.set_defaults(handler=cmd_resolve)

    classify_p = sub.add_parser("classify", help="Assign L1, healthcare-environment, and paper rollup")
    classify_p.add_argument("input", type=Path, help="Resolved or user metadata TSV")
    classify_p.add_argument("-o", "--out", type=Path, required=True, help="Classified TSV")
    classify_p.add_argument("--overrides", type=Path, help="Curator override TSV; complete rows always win")
    classify_p.add_argument("--lexicon", type=Path,
                            help=f"Lexicon YAML, or a version that ships with onehealth such as v2 (default: seeded v{SEED_LEXICON_VERSION})")
    classify_p.add_argument("--review", type=Path, help="Review-queue TSV (default: <out>.review.tsv)")
    classify_p.set_defaults(handler=cmd_classify)

    summarize_p = sub.add_parser(
        "summarize",
        help="Count classes and source evidence; overlap only if cluster/plasmid columns exist",
    )
    summarize_p.add_argument("input", type=Path, help="Classified TSV")
    summarize_p.add_argument("-o", "--out", type=Path, required=True, help="Summary TSV")
    summarize_p.set_defaults(handler=cmd_summarize)

    plot_p = sub.add_parser("plot", help="Bar charts of L1 and paper-rollup counts")
    plot_p.add_argument("input", type=Path, help="Classified TSV")
    plot_p.add_argument("-o", "--out", type=Path, required=True, help="PNG path")
    plot_p.set_defaults(handler=cmd_plot)

    validate_p = sub.add_parser(
        "validate",
        help="Compare onehealth classes with curated labels (coverage, accuracy, confusion matrix)",
    )
    validate_p.add_argument("input", type=Path, help="Classified TSV (from onehealth classify)")
    validate_p.add_argument("--truth", type=Path, required=True, help="TSV with an isolate ID column and a label column")
    validate_p.add_argument("--truth-column", default="truth", help="Label column in --truth (default: truth)")
    validate_p.add_argument("--level", choices=("paper_rollup", "l1"), default="paper_rollup",
                            help="Which onehealth class to compare (default: paper_rollup)")
    validate_p.add_argument("-o", "--out", type=Path, required=True,
                            help="Summary TSV; per-class, confusion and disagreement tables are written next to it")
    validate_p.set_defaults(handler=cmd_validate)

    report_p = sub.add_parser(
        "report",
        help="Standard One Health analyses with an HTML report, figures and tables",
        description=(
            "Lineage diversity and sharing, resistance/virulence by source, and genomic linkage "
            "between sources, written as report.html plus figures/ (SVG, PNG) and tables/ (TSV)."
        ),
    )
    report_p.add_argument("input", type=Path, help="Classified TSV (from onehealth classify)")
    report_p.add_argument("-o", "--out", type=Path, required=True, help="Output directory")
    report_p.add_argument("--features", type=Path, help="Long table: isolate, feature[, class] (AMRFinderPlus --name output fits)")
    report_p.add_argument("--distances", type=Path, help="Pairwise SNP distances from snp-dists (square or -m long format)")
    report_p.add_argument("--lineages", type=Path, help="Typing table with an ST/lineage/cluster column (Kleborate, mlst, PopPUNK)")
    report_p.add_argument("--lineage-column", help="Lineage column in the input (default: first of st, mlst, lineage, cluster ...)")
    report_p.add_argument(
        "--threshold",
        type=_threshold,
        default=21,
        help="SNP threshold for close relatives (default 21, K. pneumoniae; choose one for your organism and pipeline)",
    )
    report_p.add_argument("--group-by", choices=("source", "paper_rollup"), default="source")
    report_p.add_argument(
        "--stratify-by",
        type=_strata,
        default=(),
        metavar="FIELDS",
        help="Shuffle source labels only within strata in every permutation test: year, country, or year,country",
    )
    report_p.add_argument(
        "--dereplicate",
        action="store_true",
        help="Lineage and feature analyses use one isolate per genomic cluster per source (within --threshold SNPs); "
        "needs --distances",
    )
    report_p.add_argument("--permutations", type=_permutations, default=9999, metavar="N",
                          help="Permutations for every test (default 9999; 0 skips the tests for a quick look)")
    report_p.add_argument("--seed", type=_seed, default=1, help="Random seed, 0 or more (default 1)")
    report_p.add_argument("--title", default="One Health population report")
    report_p.add_argument("--note", default="", help="Banner text under the title, e.g. the dataset description")
    report_p.set_defaults(handler=cmd_report)

    return parser


def cmd_validate(args: argparse.Namespace) -> int:
    from onehealth.validate import read_truth, validate

    rows = read_tsv(args.input)
    require_columns(rows, args.input, [args.level], _CLASSIFY_FIRST)
    result = validate(rows, read_truth(args.truth, args.truth_column, args.level), args.level)
    summary = [
        {"metric": "isolates", "value": result.isolates, "meaning": "isolates with a truth label found in the input"},
        {"metric": "coverage", "value": f"{result.coverage:.4f}", "meaning": "share of truly classifiable isolates that onehealth classified"},
        {"metric": "accuracy", "value": f"{result.accuracy:.4f}", "meaning": "share of classified isolates whose class matches the truth"},
        {"metric": "agreement", "value": f"{result.agreement:.4f}", "meaning": "share of all isolates where class and truth agree (Unclassified included)"},
        {"metric": "unmatched_truth_ids", "value": len(result.unmatched), "meaning": "truth IDs not found in the input"},
    ]
    write_tsv(args.out, summary, ["metric", "value", "meaning"])
    by_class = [{k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()} for row in result.per_class()]
    write_tsv(_sidecar(args.out, "by_class"), by_class, ["class", "truth", "predicted", "correct", "recall", "precision"])
    classes = result.classes
    confusion = [{"truth": t, **{p: result.count(t, p) for p in classes}} for t in classes]
    write_tsv(_sidecar(args.out, "confusion"), confusion, ["truth", *classes])
    write_tsv(_sidecar(args.out, "disagreements"), result.disagreements,
              ["id", "truth", "predicted", "host", "isolation_source", "sample_type", "classify_reason"])

    print(f"validate: {result.isolates} isolates; coverage {_percent(result.coverage)}, accuracy "
          f"{_percent(result.accuracy)}, agreement {_percent(result.agreement)}", file=sys.stderr)
    if result.unmatched:
        print(f"validate: {len(result.unmatched)} truth ID{'s' if len(result.unmatched) > 1 else ''} matched no row "
              f"(e.g. {', '.join(result.unmatched[:3])})", file=sys.stderr)
    print("truth \\ onehealth  " + "  ".join(f"{p:>13}" for p in classes), file=sys.stderr)
    for row in confusion:
        print(f"{row['truth']:<19}" + "  ".join(f"{row[p]:>13}" for p in classes), file=sys.stderr)
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    from onehealth.report import build_report  # matplotlib/numpy load only when needed

    if args.dereplicate and not args.distances:
        print("onehealth report: --dereplicate needs --distances", file=sys.stderr)
        return 2

    rows = read_tsv(args.input)
    require_columns(rows, args.input, ["l1", "paper_rollup"] if args.group_by == "paper_rollup" else ["l1"], _CLASSIFY_FIRST)
    path = build_report(
        rows,
        args.out,
        features_path=args.features,
        distances_path=args.distances,
        lineages_path=args.lineages,
        lineage_column=args.lineage_column,
        threshold=args.threshold,
        by=args.group_by,
        permutations=args.permutations,
        seed=args.seed,
        title=args.title,
        note=args.note,
        dereplicate=args.dereplicate,
        stratify=args.stratify_by,
    )
    print(f"wrote {path}", file=sys.stderr)
    return 0


def cmd_resolve(args: argparse.Namespace) -> int:
    ids = read_id_list(args.ids)
    local = load_local_metadata(args.metadata)
    client = None
    allow_network = not args.no_network
    if allow_network:
        if not args.email:
            print(
                "warning: NCBI asks E-utilities users to send --email (and tool) so it can "
                "contact you before blocking an IP",
                file=sys.stderr,
            )
        client = NCBIClient(email=args.email, api_key=args.api_key, log=lambda note: print(f"resolve: {note}", file=sys.stderr))
    elif args.pd_organism:
        print("warning: --pd-organism needs NCBI, so it has no effect with --no-network", file=sys.stderr)
    rows = resolve_ids(ids, local_rows=local, client=client, allow_network=allow_network, pd_groups=list(args.pd_organism))
    # Union over all rows, so optional cluster/plasmid columns survive even when
    # the first row has none.
    extra = sorted({key for row in rows for key in row} - set(RESOLVE_FIELDS))
    fieldnames = [*RESOLVE_FIELDS, *extra]
    write_tsv(args.out, rows, fieldnames)
    _write_review(args.review or _sidecar(args.out, "review"), rows, fieldnames, status_field="resolve_status")
    _print_counts("resolve", rows, "resolve_status")
    unresolved = sum(1 for row in rows if row.get("resolve_reason") == "pdt_unresolved")
    if unresolved and allow_network and not args.pd_organism:
        print(f"resolve: {unresolved} PDT IDs were not found by assembly name; pass --pd-organism (for example "
              "Klebsiella) to look them up in NCBI Pathogen Detection's own table", file=sys.stderr)
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    rows = read_tsv(args.input)
    lexicon = load_lexicon(args.lexicon)
    overrides = load_overrides(args.overrides)
    _check_classify_columns(rows, args.input, has_overrides=bool(overrides))
    classified = classify_rows(rows, lexicon, overrides)
    fieldnames = classify_fieldnames(classified)
    write_tsv(args.out, classified, fieldnames)
    _write_review(args.review or _sidecar(args.out, "review"), classified, fieldnames, status_field="classify_status")
    _print_counts("classify", classified, "classify_status")
    print(f"lexicon v{lexicon.version} (seed v{SEED_LEXICON_VERSION})", file=sys.stderr)
    missed = unmatched_overrides(rows, overrides)
    if missed:
        # Never applied, and nothing in the output would show it.
        print(f"warning: {len(missed)} override ID{'s' if len(missed) > 1 else ''} in {args.overrides} matched no row "
              f"and {'were' if len(missed) > 1 else 'was'} not applied (e.g. {', '.join(missed[:3])}). A strain or "
              "sample name only matches a table with no BioSample, assembly or ID column.", file=sys.stderr)
    return 0


def cmd_summarize(args: argparse.Namespace) -> int:
    rows = read_tsv(args.input)
    require_columns(rows, args.input, ["l1"], _CLASSIFY_FIRST)
    summary = summarize(rows)
    write_tsv(args.out, summary, summary_fieldnames())
    print(f"wrote {args.out} ({len(summary)} rows)", file=sys.stderr)
    return 0


def cmd_plot(args: argparse.Namespace) -> int:
    rows = read_tsv(args.input)
    require_columns(rows, args.input, ["l1", "paper_rollup"], _CLASSIFY_FIRST)
    plot_counts(rows, args.out)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


_CLASSIFY_FIRST = "Is this the output of onehealth classify?"

# Every column the lexicon can read (io_tsv.metadata_view).
_EVIDENCE_COLUMNS = (*ISOLATION_SOURCE_ALIASES, *HOST_ALIASES, *SAMPLE_TYPE_ALIASES, *HOST_DISEASE_ALIASES)


def _check_classify_columns(rows: list[dict[str, str]], path: Path, has_overrides: bool) -> None:
    """A table with none of the columns classify reads would send every row to review.

    The usual cause is an ID list given to classify instead of resolve; its first ID is
    then read as the header. Overrides can still classify the rows they name, but only
    if the table has an ID column to match them on.
    """
    if not rows or first_present(rows[0], _EVIDENCE_COLUMNS) is not None:
        return
    found = ", ".join(list(rows[0])[:8]) + (", ..." if len(rows[0]) > 8 else "")
    problem = (f"{path}: no column classify can read (isolation_source, host, source_type, ...; found: {found})")
    if has_overrides and first_present(rows[0], (*ID_LIST_ALIASES, *ISOLATE_NAME_ALIASES)) is not None:
        print(f"warning: {problem}, so only rows with an override can be classified", file=sys.stderr)
        return
    raise InputError(f"{problem}. For a list of IDs, run onehealth resolve first.")


def _percent(value: float) -> str:
    return "n/a" if math.isnan(value) else f"{value:.1%}"


def _names(text: str) -> tuple[str, ...]:
    return tuple(dict.fromkeys(part.strip() for part in text.split(",") if part.strip()))


def _threshold(text: str) -> float:
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {text!r}") from None
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be a number of SNPs, 0 or more")
    return value


def _seed(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a whole number: {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError("must be 0 or more")
    return value


def _permutations(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a whole number: {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError("must be 0 or more (0 skips the permutation tests)")
    return value


def _strata(text: str) -> tuple[str, ...]:
    fields = tuple(dict.fromkeys(part.strip().lower() for part in text.split(",") if part.strip()))
    unknown = [f for f in fields if f not in ("year", "country")]
    if unknown:
        raise argparse.ArgumentTypeError(f"unknown strata {unknown}; use year, country or year,country")
    return fields


def _sidecar(path: Path, suffix: str) -> Path:
    return path.with_name(f"{path.stem}.{suffix}{path.suffix}")


def _write_review(path: Path, rows: list[dict[str, str]], fieldnames: list[str], status_field: str) -> None:
    review = [row for row in rows if row.get(status_field) == "review"]
    write_tsv(path, review, fieldnames)


def _print_counts(stage: str, rows: list[dict[str, str]], status_field: str) -> None:
    ok = sum(1 for row in rows if row.get(status_field) == "ok")
    review = sum(1 for row in rows if row.get(status_field) == "review")
    print(f"{stage}: {len(rows)} rows ({ok} ok, {review} review)", file=sys.stderr)
