"""onehealth resolve | classify | summarize | plot"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from onehealth import __version__
from onehealth.classify import classify_fieldnames, classify_rows
from onehealth.constants import SEED_LEXICON_VERSION
from onehealth.io_tsv import read_id_list, read_tsv, write_tsv
from onehealth.lexicon import load_lexicon
from onehealth.ncbi import NCBIClient
from onehealth.overrides import load_overrides
from onehealth.plot import plot_counts
from onehealth.resolve import load_local_metadata, resolve_ids
from onehealth.summarize import summarize, summary_fieldnames


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        return 2
    return args.handler(args)


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
    resolve_p.set_defaults(handler=cmd_resolve)

    classify_p = sub.add_parser("classify", help="Assign L1, healthcare-environment, and paper rollup")
    classify_p.add_argument("input", type=Path, help="Resolved or user metadata TSV")
    classify_p.add_argument("-o", "--out", type=Path, required=True, help="Classified TSV")
    classify_p.add_argument("--overrides", type=Path, help="Curator override TSV; complete rows always win")
    classify_p.add_argument("--lexicon", type=Path, help="Lexicon YAML (default: seeded v1)")
    classify_p.add_argument("--review", type=Path, help="Review-queue TSV (default: <out>.review.tsv)")
    classify_p.set_defaults(handler=cmd_classify)

    summarize_p = sub.add_parser("summarize", help="Count classes; overlap only if cluster/plasmid columns exist")
    summarize_p.add_argument("input", type=Path, help="Classified TSV")
    summarize_p.add_argument("-o", "--out", type=Path, required=True, help="Summary TSV")
    summarize_p.set_defaults(handler=cmd_summarize)

    plot_p = sub.add_parser("plot", help="Bar charts of L1 and paper-rollup counts")
    plot_p.add_argument("input", type=Path, help="Classified TSV")
    plot_p.add_argument("-o", "--out", type=Path, required=True, help="PNG path")
    plot_p.add_argument("--by", choices=("l1", "paper_rollup"), default="l1")
    plot_p.set_defaults(handler=cmd_plot)

    return parser


def cmd_resolve(args: argparse.Namespace) -> int:
    ids = read_id_list(args.ids)
    local = load_local_metadata(args.metadata)
    client = None
    allow_network = not args.no_network
    if allow_network:
        client = NCBIClient(email=args.email, api_key=args.api_key)
    rows = resolve_ids(ids, local_rows=local, client=client, allow_network=allow_network)
    # Preserve optional cluster/plasmid columns from local metadata when present.
    fieldnames = list(rows[0].keys()) if rows else [
        "query_id",
        "id_type",
        "biosample",
        "assembly",
        "sra",
        "pdt",
        "organism",
        "host",
        "isolation_source",
        "host_disease",
        "sample_type",
        "geo_loc_name",
        "collection_date",
        "strain",
        "isolate",
        "resolve_status",
        "resolve_reason",
    ]
    write_tsv(args.out, rows, fieldnames)
    _write_review(args.review or _sidecar(args.out, "review"), rows, fieldnames, status_field="resolve_status")
    _print_counts("resolve", rows, "resolve_status")
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    rows = read_tsv(args.input)
    lexicon = load_lexicon(args.lexicon)
    overrides = load_overrides(args.overrides)
    classified = classify_rows(rows, lexicon, overrides)
    fieldnames = classify_fieldnames(classified)
    write_tsv(args.out, classified, fieldnames)
    _write_review(args.review or _sidecar(args.out, "review"), classified, fieldnames, status_field="classify_status")
    _print_counts("classify", classified, "classify_status")
    print(f"lexicon v{lexicon.version} (seed v{SEED_LEXICON_VERSION})", file=sys.stderr)
    return 0


def cmd_summarize(args: argparse.Namespace) -> int:
    rows = read_tsv(args.input)
    summary = summarize(rows)
    write_tsv(args.out, summary, summary_fieldnames())
    print(f"wrote {args.out} ({len(summary)} rows)", file=sys.stderr)
    return 0


def cmd_plot(args: argparse.Namespace) -> int:
    rows = read_tsv(args.input)
    plot_counts(rows, args.out, by=args.by)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


def _sidecar(path: Path, suffix: str) -> Path:
    return path.with_name(f"{path.stem}.{suffix}{path.suffix}")


def _write_review(path: Path, rows: list[dict[str, str]], fieldnames: list[str], status_field: str) -> None:
    review = [row for row in rows if row.get(status_field) == "review"]
    write_tsv(path, review, fieldnames)


def _print_counts(stage: str, rows: list[dict[str, str]], status_field: str) -> None:
    ok = sum(1 for row in rows if row.get(status_field) == "ok")
    review = sum(1 for row in rows if row.get(status_field) == "review")
    print(f"{stage}: {len(rows)} rows ({ok} ok, {review} review)", file=sys.stderr)
