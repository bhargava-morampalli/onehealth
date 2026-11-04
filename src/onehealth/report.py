"""`onehealth report`: run the standard analyses and write a self-contained HTML
report, publication figures (SVG + PNG) and every underlying table (TSV).

Design rules (palette, marks, accessibility) follow a validated categorical
palette: sources keep a fixed colour everywhere, every chart has a data table,
and the HTML adapts to light and dark mode.
"""

from __future__ import annotations

import html
import io
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from onehealth import __version__  # noqa: E402
from onehealth import analysis as A  # noqa: E402
from onehealth.io_tsv import read_distance_matrix, read_features, read_lineages, write_tsv  # noqa: E402

# --- theme ---------------------------------------------------------------------------
# token: (light, dark). Light hexes are unique so inline SVG can be re-themed by CSS.
TOKENS = {
    "surface": ("#fcfcfb", "#1a1a19"),
    "page": ("#f9f9f7", "#0d0d0d"),
    "ink": ("#0b0b0b", "#f0efec"),
    "ink-2": ("#52514e", "#c3c2b7"),
    "muted": ("#898781", "#898781"),
    "grid": ("#e1e0d9", "#2c2c2a"),
    "axis": ("#c3c2b7", "#383835"),
    "accent": ("#6250d6", "#9085e9"),
    "accent-wash": ("#d9d4f6", "#38345c"),
    "absent": ("#eceae4", "#2a2a28"),
    "s-human": ("#2a78d6", "#3987e5"),
    "s-animal": ("#eb6834", "#d95926"),
    "s-environmental": ("#1baf7a", "#199e70"),
    "s-food": ("#eda100", "#c98500"),
    "s-healthcare": ("#e87ba4", "#d55181"),
    "s-laboratory": ("#008300", "#018401"),
    "s-unclassified": ("#a3a29c", "#6b6a65"),
    # binned sequential (blue): zero, then five steps
    "q0": ("#f0efec", "#262624"),
    "q1": ("#cde2fb", "#16467f"),
    "q2": ("#9ec5f4", "#1c5cab"),
    "q3": ("#5598e7", "#2c73cc"),
    "q4": ("#256abf", "#5598e7"),
    "q5": ("#104281", "#9ec5f4"),
}
LIGHT = {name: pair[0] for name, pair in TOKENS.items()}
SOURCE_TOKEN = {
    "Human": "s-human",
    "Animal": "s-animal",
    "Environmental": "s-environmental",
    "Food": "s-food",
    "Healthcare environment": "s-healthcare",
    "Laboratory": "s-laboratory",
    "Unclassified": "s-unclassified",
}
SOURCE_MARKER = {"Human": "o", "Animal": "s", "Environmental": "^", "Food": "D", "Healthcare environment": "v", "Laboratory": "P"}
SEQ = ["q1", "q2", "q3", "q4", "q5"]


def colour(source: str) -> str:
    return LIGHT[SOURCE_TOKEN.get(source, "s-unclassified")]


plt.rcParams.update(
    {
        # DejaVu Sans ships with matplotlib: identical figures on every machine, full glyph set.
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "svg.fonttype": "path",
        "svg.hashsalt": "onehealth",
        "figure.facecolor": LIGHT["surface"],
        "axes.facecolor": LIGHT["surface"],
        "savefig.facecolor": LIGHT["surface"],
        "axes.edgecolor": LIGHT["axis"],
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": LIGHT["ink-2"],
        "axes.labelsize": 10.5,
        "axes.titlesize": 10,
        "axes.titlecolor": LIGHT["ink"],
        "text.color": LIGHT["ink"],
        "xtick.color": LIGHT["axis"],
        "ytick.color": LIGHT["axis"],
        "xtick.labelcolor": LIGHT["ink-2"],
        "ytick.labelcolor": LIGHT["ink-2"],
        "grid.color": LIGHT["grid"],
        "grid.linewidth": 0.8,
        "legend.frameon": False,
        "legend.fontsize": 10,
        "lines.linewidth": 1.6,  # ~2px
        "lines.solid_capstyle": "round",
    }
)


class Figure:
    """A matplotlib figure plus per-mark tooltip text."""

    _counter = 0  # tooltip ids must be unique across the whole HTML page

    def __init__(self, width: float, height: float, **kwargs):
        self.fig, self.ax = plt.subplots(figsize=(width, height), **kwargs)
        self.tips: dict[str, str] = {}

    @classmethod
    def wrap(cls, fig, ax) -> "Figure":
        self = cls.__new__(cls)
        self.fig, self.ax, self.tips = fig, ax, {}
        return self

    def tip(self, artist, text: str) -> None:
        Figure._counter += 1
        gid = f"tip{Figure._counter}"
        artist.set_gid(gid)
        self.tips[gid] = text


def _hgrid(ax) -> None:
    ax.grid(axis="y")
    ax.set_axisbelow(True)


def _vgrid(ax) -> None:
    ax.grid(axis="x")
    ax.set_axisbelow(True)


def _legend_sources(ax, sources, loc="upper left", anchor=(0, 1.02), ncol=None, kind="patch") -> None:
    if kind == "patch":
        handles = [Patch(facecolor=colour(s), label=s) for s in sources]
    else:  # legends mirror the mark: a short line key for line charts
        handles = [Line2D([], [], color=colour(s), linewidth=2, label=s) for s in sources]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=anchor, ncol=ncol or len(sources),
              handlelength=1.0, handleheight=0.9, columnspacing=1.2, borderaxespad=0.2)


def _gap() -> dict:
    """2px surface gap between touching fills."""
    return {"edgecolor": LIGHT["surface"], "linewidth": 1.2}


def pct(value: float, digits: int = 1) -> str:
    return "n/a" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.{digits}f}%"


def fmt_p(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    return "<0.001" if value < 0.001 else f"{value:.3f}"


# --- figures ---------------------------------------------------------------------------


def fig_sampling(dataset: A.Dataset) -> Figure:
    years = A.counts_by_year(dataset)
    sources = [g for g in dataset.groups if any(iso.source == g for iso in dataset.classified)]
    if len(years) >= 2:
        f = Figure(8.2, 3.2)
        bottom = np.zeros(len(years))
        x = np.arange(len(years))
        for s in sources:
            values = np.array([row[s] for row in years], dtype=float)
            bars = f.ax.bar(x, values, width=0.62, bottom=bottom, color=colour(s), **_gap())
            for bar, row, v in zip(bars, years, values):
                if v:
                    f.tip(bar, f"{row['year']} · {s}: {int(v)} isolates")
            bottom += values
        f.ax.set_xticks(x, [str(row["year"]) for row in years])
        f.ax.set_ylabel("Isolates")
        _hgrid(f.ax)
        _legend_sources(f.ax, sources)
    else:
        rows = A.source_counts(dataset)
        f = Figure(8.2, 0.45 * len(rows) + 0.9)
        y = np.arange(len(rows))[::-1]
        for yy, row in zip(y, rows):
            bar = f.ax.barh(yy, row["isolates"], height=0.6, color=colour(row["source"]))[0]
            f.tip(bar, f"{row['source']}: {row['isolates']} isolates ({pct(row['percent'])})")
            f.ax.text(row["isolates"], yy, f"  {row['isolates']}", va="center", color=LIGHT["ink-2"], fontsize=9)
        f.ax.set_yticks(y, [row["source"] for row in rows])
        f.ax.set_xlabel("Isolates")
        _vgrid(f.ax)
    f.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    return f


def fig_diversity(rows: list[dict]) -> Figure:
    f = Figure(8.2, 0.42 * len(rows) + 1.0)
    y = np.arange(len(rows))[::-1]
    for yy, row in zip(y, rows):
        f.ax.plot([row["ci_low"], row["ci_high"]], [yy, yy], color=LIGHT["ink-2"], linewidth=1.5, solid_capstyle="butt")
        dot = f.ax.scatter([row["simpson"]], [yy], s=70, color=colour(row["source"]), edgecolors=LIGHT["surface"],
                           linewidths=2, zorder=3)
        f.tip(dot, f"{row['source']}: D = {row['simpson']:.3f} (95% CI {row['ci_low']:.3f}–{row['ci_high']:.3f}); "
                   f"{row['lineages']} lineages in {row['isolates']} isolates")
    f.ax.set_yticks(y, [f"{r['source']}  (n={r['isolates']})" for r in rows])
    # A source with a single typed isolate has no D (NaN); min() would not skip it.
    low = min((r["ci_low"] for r in rows if not math.isnan(r["ci_low"])), default=0.0)
    f.ax.set_xlim(max(0.0, low - 0.05), 1.005)
    f.ax.set_xlabel("Simpson's index of diversity (lineages)")
    _vgrid(f.ax)
    return f


def fig_rarefaction(rows: list[dict], groups: list[str]) -> Figure:
    f = Figure(8.2, 3.4)
    for g in groups:
        pts = [r for r in rows if r["source"] == g]
        if not pts:
            continue
        x = [p["isolates_sampled"] for p in pts]
        y = [p["expected_lineages"] for p in pts]
        line = f.ax.plot(x, y, color=colour(g))[0]
        f.tip(line, f"{g}: {y[-1]:.0f} lineages in {x[-1]} isolates")
        f.ax.scatter([x[-1]], [y[-1]], s=40, color=colour(g), edgecolors=LIGHT["surface"], linewidths=2, zorder=3)
    f.ax.set_xlabel("Isolates sampled")
    f.ax.set_ylabel("Expected lineages")
    _hgrid(f.ax)
    _legend_sources(f.ax, [g for g in groups if any(r["source"] == g for r in rows)], kind="line")
    return f


def fig_lineages(rows: list[dict], groups: list[str], enriched: set[str], top: int = 20) -> Figure:
    rows = rows[:top]
    f = Figure(8.2, 0.3 * len(rows) + 1.1)
    y = np.arange(len(rows))[::-1]
    left = np.zeros(len(rows))
    for g in groups:
        values = np.array([r[g] for r in rows], dtype=float)
        bars = f.ax.barh(y, values, left=left, height=0.66, color=colour(g), **_gap())
        for bar, r, v in zip(bars, rows, values):
            if v:
                f.tip(bar, f"{r['lineage']} · {g}: {int(v)} of {r['isolates']} isolates")
        left += values
    labels = [f"{r['lineage']}{' *' if r['lineage'] in enriched else ''}" for r in rows]
    f.ax.set_yticks(y, labels)
    f.ax.set_xlabel("Isolates")
    _vgrid(f.ax)
    _legend_sources(f.ax, groups)
    return f


def fig_upset(intersections: list[dict], groups: list[str]) -> Figure:
    n = len(intersections)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.2, 3.0 + 0.28 * len(groups)), sharex=True,
                                      gridspec_kw={"height_ratios": [2.2, 0.32 * len(groups) + 0.4], "hspace": 0.06})
    f = Figure.wrap(fig, top)
    x = np.arange(n)
    for xx, row in zip(x, intersections):
        bar = top.bar(xx, row["lineages"], width=0.6, color=LIGHT["ink-2"])[0]
        f.tip(bar, f"{row['lineages']} lineages found only in: {row['sources']}")
        top.text(xx, row["lineages"], str(row["lineages"]), ha="center", va="bottom", fontsize=8.5, color=LIGHT["ink-2"])
    top.set_ylabel("Lineages")
    _hgrid(top)
    top.spines["bottom"].set_visible(False)
    top.tick_params(axis="x", length=0)
    for k, g in enumerate(groups):
        yy = len(groups) - 1 - k
        for xx, row in zip(x, intersections):
            member = g in row["members"]
            bottom.scatter([xx], [yy], s=55, color=colour(g) if member else LIGHT["absent"],
                           edgecolors=LIGHT["surface"], linewidths=1.5, zorder=3)
    for xx, row in zip(x, intersections):
        ys = [len(groups) - 1 - groups.index(g) for g in row["members"]]
        if len(ys) > 1:
            bottom.plot([xx, xx], [min(ys), max(ys)], color=LIGHT["ink-2"], linewidth=1.5, zorder=2)
    bottom.set_yticks(range(len(groups))[::-1], groups)
    bottom.set_ylim(-0.6, len(groups) - 0.4)
    bottom.set_xticks([])
    for side in ("left", "bottom"):
        bottom.spines[side].set_visible(False)
    bottom.tick_params(axis="y", length=0)
    return f


def _bin(value: float, edges: list[float]) -> str:
    if value is None or math.isnan(value) or value <= 0:
        return "q0"
    for k, edge in enumerate(edges):
        if value <= edge:
            return SEQ[k]
    return SEQ[-1]


def _bin_legend(ax, edges: list[float], unit: str, anchor=(0, 1.02)) -> None:
    """Four edges -> six classes: 0, <=e1, e1-e2, e2-e3, e3-e4, >e4."""
    labels = ["0", f"≤{edges[0]:g}{unit}"] + [f"{a:g}–{b:g}{unit}" for a, b in zip(edges, edges[1:])]
    labels.append(f">{edges[-1]:g}{unit}")
    tokens = ["q0", *SEQ[: len(edges) + 1]]
    handles = [Patch(facecolor=LIGHT[t], label=label) for t, label in zip(tokens, labels)]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=anchor, ncol=len(handles), handlelength=1.0,
              columnspacing=1.0, borderaxespad=0.2)


def _heatmap(f: Figure, rows: list[str], cols: list[str], value, tip, edges: list[float], unit: str) -> None:
    ax = f.ax
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            v = value(r, c)
            cell = Rectangle((j, len(rows) - 1 - i), 1, 1, facecolor=LIGHT[_bin(v, edges)], edgecolor=LIGHT["surface"],
                             linewidth=1.5)
            ax.add_patch(cell)
            f.tip(cell, tip(r, c))
    ax.set_xlim(0, len(cols))
    ax.set_ylim(0, len(rows))
    ax.set_xticks(np.arange(len(cols)) + 0.5, cols, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(rows))[::-1] + 0.5, rows)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    _bin_legend(ax, edges, unit)


def fig_features(result: A.FeatureResult, top: int = 40) -> Figure | None:
    tests = [t for t in result.tests if 0 < t["percent_overall"] < 100]
    tests = sorted(tests, key=lambda t: -t["percent_overall"])[:top]  # most common variable features
    tests.sort(key=lambda t: (t["class"], -t["percent_overall"], t["feature"]))
    if not tests:
        return None
    lookup = {(r["feature"], r["source"]): r for r in result.rows}
    labels = {t["feature"]: f"{t['feature']}{'  *' if t['q'] < 0.05 else ''}" for t in tests}
    names = [t["feature"] for t in tests]
    f = Figure(2.2 + 0.85 * len(result.groups), 0.27 * len(names) + 1.6)

    def value(r, c):
        return lookup[(r, c)]["percent"]

    def tip(r, c):
        row = lookup[(r, c)]
        return (f"{r} · {c}: {row['present']}/{row['genotyped']} ({pct(row['percent'])}, "
                f"95% CI {pct(row['ci_low'], 0)}–{pct(row['ci_high'], 0)})")

    _heatmap(f, names, result.groups, value, tip, [10, 25, 50, 75], "%")
    f.ax.set_yticks(np.arange(len(names))[::-1] + 0.5, [labels[n] for n in names])
    classes = [next(t["class"] for t in tests if t["feature"] == n) for n in names]
    if any(classes):
        sec = f.ax.secondary_yaxis("right")
        sec.set_yticks(np.arange(len(names))[::-1] + 0.5, [c.capitalize() if c else "" for c in classes])
        sec.tick_params(length=0, labelcolor=LIGHT["muted"], labelsize=9)
        sec.spines["right"].set_visible(False)
    return f


def fig_distance_histogram(link: A.LinkageResult, limit: int = 100) -> Figure:
    f = Figure(8.2, 3.2)
    bins = np.arange(0, limit + 2, 2)
    for values, label, token in ((link.close_within, "Same source", "muted"),
                                 (link.close_between, "Different sources", "accent")):
        close = values[values <= limit]
        counts, _ = np.histogram(close, bins=bins)
        line = f.ax.step(bins[:-1], counts, where="post", color=LIGHT[token], label=f"{label} ({len(close)} pairs ≤{limit} SNPs)")[0]
        f.tip(line, f"{label}: {len(close)} pairs within {limit} SNPs ({int((values <= link.threshold).sum())} "
                    f"within the {link.threshold:g}-SNP threshold)")
    f.ax.axvline(link.threshold, color=LIGHT["ink-2"], linewidth=1)
    f.ax.text(link.threshold, f.ax.get_ylim()[1], f"  threshold {link.threshold:g} SNPs", va="top", fontsize=8.5,
              color=LIGHT["ink-2"])
    f.ax.set_xlabel("Pairwise SNP distance")
    f.ax.set_ylabel("Isolate pairs")
    f.ax.set_xlim(0, limit)
    _hgrid(f.ax)
    f.ax.legend(loc="upper right")
    return f


def fig_sweep(link: A.LinkageResult) -> Figure:
    f = Figure(8.2, 3.0)
    x = np.arange(len(link.sweep))
    for key, label, token in (("within_events", "Within a source", "muted"), ("between_events", "Between sources", "accent")):
        y = [row[key] for row in link.sweep]
        line = f.ax.plot(x, y, color=LIGHT[token], marker="o", markersize=6, markeredgecolor=LIGHT["surface"],
                         markeredgewidth=1.5, label=label)[0]
        f.tip(line, label + ": " + ", ".join(f"{row['threshold']:g} SNPs → {row[key]}" for row in link.sweep))
    chosen = next(i for i, row in enumerate(link.sweep) if row["threshold"] == link.threshold)
    f.ax.axvline(chosen, color=LIGHT["ink-2"], linewidth=1)
    f.ax.set_xticks(x, [f"{row['threshold']:g}" for row in link.sweep])
    f.ax.set_xlabel("SNP threshold")
    f.ax.set_ylabel("Sharing events")
    _hgrid(f.ax)
    f.ax.legend(loc="upper left")
    return f


def fig_linked(link: A.LinkageResult) -> Figure:
    lookup = {(r["from_source"], r["to_source"]): r for r in link.linked_fraction}
    groups = link.groups
    f = Figure(2.4 + 0.9 * len(groups), 0.5 * len(groups) + 1.9)

    def value(r, c):
        return lookup[(r, c)]["percent"]

    def tip(r, c):
        row = lookup[(r, c)]
        return f"{row['with_close_relative']} of {row['isolates']} {r} isolates ({pct(row['percent'])}) have a {c} isolate within {link.threshold:g} SNPs"

    _heatmap(f, groups, groups, value, tip, [1, 5, 10, 25], "%")
    f.ax.set_ylabel("Isolates from…")
    f.ax.set_xlabel("…with a close relative from")
    return f


def fig_observed_expected(link: A.LinkageResult) -> Figure | None:
    rows = [r for r in link.events if r["kind"] == "between" and (r["events"] or r["expected_high"] > 0)]
    if not rows:
        return None
    rows.sort(key=lambda r: -(r["events"] - r["expected"]))
    f = Figure(8.2, 0.36 * len(rows) + 1.2)
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows):
        band = f.ax.barh(yy, r["expected_high"] - r["expected_low"], left=r["expected_low"], height=0.5,
                         color=LIGHT["accent-wash"])[0]
        f.tip(band, f"{r['source_a']} ↔ {r['source_b']}: expected {r['expected']:.1f} "
                    f"(95% range {r['expected_low']:.0f}–{r['expected_high']:.0f}) if sharing ignored source")
        dot = f.ax.scatter([r["events"]], [yy], s=60, color=LIGHT["accent"], edgecolors=LIGHT["surface"], linewidths=2,
                           zorder=3)
        f.tip(dot, f"{r['source_a']} ↔ {r['source_b']}: {r['events']} observed events (q = {fmt_p(r['q'])})")
    f.ax.set_yticks(y, [f"{r['source_a']} ↔ {r['source_b']}{'  *' if r['q'] < 0.05 else ''}" for r in rows])
    f.ax.set_xlim(left=-0.8)
    f.ax.set_xlabel(f"Cross-source sharing events at ≤{link.threshold:g} SNPs")
    _vgrid(f.ax)
    handles = [Patch(facecolor=LIGHT["accent-wash"], label="Expected range (95%) if source did not matter"),
               Line2D([], [], marker="o", linestyle="", markersize=7, markerfacecolor=LIGHT["accent"],
                      markeredgecolor=LIGHT["surface"], label="Observed")]
    f.ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, borderaxespad=0.2)
    return f


def fig_clusters(link: A.LinkageResult, limit: int = 12) -> Figure | None:
    clusters = [c for c in link.clusters if c["n_sources"] >= 2][:limit]
    if not clusters:
        return None
    cols = min(4, len(clusters))
    rows = math.ceil(len(clusters) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(1.8 * cols + 0.4, 1.9 * rows + 0.5), squeeze=False)
    f = Figure.wrap(fig, axes[0][0])
    for ax in axes.flat:
        ax.set_axis_off()
    present = set()
    for ax, c in zip(axes.flat, clusters):
        n = len(c["members"])
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
        xs, ys = np.cos(angles), np.sin(angles)
        if n == 2:
            xs, ys = np.array([-0.6, 0.6]), np.array([0.0, 0.0])
        for i, j, d in _spanning_tree(n, c["edges"]):
            ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color=LIGHT["axis"], linewidth=1.4, zorder=1)
            ax.text((xs[i] + xs[j]) / 2, (ys[i] + ys[j]) / 2, f"{d:g}", fontsize=8.5, color=LIGHT["ink-2"], ha="center",
                    va="center", zorder=2, bbox={"boxstyle": "round,pad=0.15", "facecolor": LIGHT["surface"], "edgecolor": "none"})
        for k, (key, source) in enumerate(zip(c["members"], c["member_sources"])):
            label = source or "Unclassified"
            present.add(label)
            dot = ax.scatter([xs[k]], [ys[k]], s=90, marker=SOURCE_MARKER.get(label, "o"), color=colour(label),
                             edgecolors=LIGHT["surface"], linewidths=2, zorder=3)
            f.tip(dot, f"{key} · {label}")
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.45, 1.6)
        lineage = c["lineages"] or "untyped"
        ax.set_title(f"{lineage} · {n} isolates", fontsize=9.5, color=LIGHT["ink-2"], pad=2)
    order = [g for g in [*link.groups, "Unclassified"] if g in present]
    handles = [Line2D([], [], marker=SOURCE_MARKER.get(s, "o"), linestyle="", markersize=7, markerfacecolor=colour(s),
                      markeredgecolor=LIGHT["surface"], label=s) for s in order]
    fig.legend(handles=handles, loc="upper left", ncol=len(handles), frameon=False, bbox_to_anchor=(0.01, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return f


def _spanning_tree(n: int, edges: list[tuple[int, int, float]]) -> list[tuple[int, int, float]]:
    """Minimum spanning tree (Kruskal) of a cluster's within-threshold edges, for a legible drawing."""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree = []
    for i, j, d in sorted(edges, key=lambda e: e[2]):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
            tree.append((i, j, d))
    return tree


# --- rendering -----------------------------------------------------------------------------

_HEX_TO_VAR = sorted(((light.lower(), name) for name, (light, _) in TOKENS.items()), key=lambda item: -len(item[0]))


def _render(fig: Figure, stem: str, figures_dir: Path) -> str:
    """Write <stem>.svg/.png (light theme) and return inline SVG re-themed with CSS variables."""
    meta = {"Date": None, "Creator": f"onehealth {__version__}"}
    fig.fig.savefig(figures_dir / f"{stem}.svg", format="svg", bbox_inches="tight", metadata=meta)
    fig.fig.savefig(figures_dir / f"{stem}.png", format="png", dpi=220, bbox_inches="tight")
    buffer = io.StringIO()
    fig.fig.savefig(buffer, format="svg", bbox_inches="tight", metadata=meta)
    plt.close(fig.fig)
    svg = buffer.getvalue()
    svg = svg[svg.index("<svg") :]
    # Natural size (pt -> px), shrinking on narrow screens; text then matches the page scale.
    width = re.search(r'<svg[^>]*?\swidth="([\d.]+)pt"', svg)
    natural = f"width:{float(width.group(1)) * 96 / 72:.0f}px;" if width else ""
    svg = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"\s+height="[^"]*"', rf'\1 style="{natural}max-width:100%;height:auto"', svg, count=1)
    svg = svg.replace("<svg ", f'<svg role="img" aria-labelledby="{stem}-title" ', 1)
    for hex_value, name in _HEX_TO_VAR:
        svg = re.sub(re.escape(hex_value), f"var(--{name})", svg, flags=re.IGNORECASE)
    for gid, text in fig.tips.items():
        svg = svg.replace(f'<g id="{gid}">', f'<g id="{gid}" class="tip" tabindex="0" data-tip="{html.escape(text, quote=True)}">')
    return svg


def _table_html(rows: list[dict], columns: list[tuple[str, str]], limit: int = 400) -> str:
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = []
    for row in rows[:limit]:
        cells = []
        for key, _ in columns:
            value = row.get(key, "")
            if isinstance(value, float):
                text = _fmt(value)
            elif isinstance(value, (list, tuple, set)):
                text = "; ".join(map(str, value))
            else:
                text = str(value)
            cells.append(f"<td>{html.escape(text)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    more = f'<p class="muted small">Showing {limit} of {len(rows)} rows; the TSV has all of them.</p>' if len(rows) > limit else ""
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>{more}'


def _fmt(value: float) -> str:
    if math.isnan(value):
        return ""
    if value and abs(value) < 0.01:
        return f"{value:.2g}"
    return f"{value:,.3f}".rstrip("0").rstrip(".")


def _write_table(tables_dir: Path, stem: str, rows: list[dict]) -> str:
    clean = [{k: v for k, v in row.items() if not isinstance(v, (list, tuple, set))} for row in rows]
    columns = list(clean[0]) if clean else []
    formatted = [
        {k: (f"{v:.6g}" if isinstance(v, float) and not math.isnan(v) else "" if isinstance(v, float) else v) for k, v in row.items()}
        for row in clean
    ]
    write_tsv(tables_dir / f"{stem}.tsv", formatted, columns)
    return f"tables/{stem}.tsv"


# --- report ---------------------------------------------------------------------------------


def build_report(
    rows: list[dict[str, str]],
    out_dir: Path,
    *,
    features_path: Path | None = None,
    distances_path: Path | None = None,
    lineages_path: Path | None = None,
    lineage_column: str | None = None,
    threshold: float = 21,
    by: str = "source",
    permutations: int = 9999,
    seed: int = 1,
    title: str = "One Health population report",
    note: str = "",
    dereplicate: bool = False,
    stratify: tuple[str, ...] = (),
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir, tables_dir = out_dir / "figures", out_dir / "tables"
    figures_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)

    dataset = A.build_dataset(rows, by=by, lineage_column=lineage_column)
    quality: list[str] = []
    if lineages_path:
        unmatched = A.apply_lineages(dataset, read_lineages(lineages_path))
        if unmatched:
            quality.append(f"{len(unmatched)} names in the lineage table matched no isolate (e.g. {', '.join(unmatched[:3])}).")
    sections: list[str] = []
    kpis: list[tuple[str, str]] = []
    refs = References()
    methods: list[str] = []

    classified = dataset.classified
    groups_present = [g for g in dataset.groups if any(iso.source == g for iso in classified)]
    kpis.append(("Isolates", f"{len(dataset.isolates):,}"))
    kpis.append(("With a source class", f"{len(classified):,} ({100 * len(classified) / max(1, len(dataset.isolates)):.0f}%)"))
    kpis.append(("Source groups", str(len(groups_present))))

    names, matrix = read_distance_matrix(distances_path) if distances_path else ([], None)
    # Lineage and feature analyses can run on one representative per genomic cluster per source.
    population, derep = dataset, None
    if dereplicate:
        if matrix is None:
            raise ValueError("--dereplicate needs --distances")
        population, derep = A.dereplicate(dataset, names, matrix, threshold)
        kpis.append(("After dereplication", f"{derep.isolates_after:,} of {derep.isolates_before:,}"))
        if derep.not_in_matrix:
            quality.append(f"{derep.not_in_matrix} classified isolates have no distances and were kept as their own "
                           "representatives during dereplication.")
    derep_note = " Counts use one representative per genomic cluster per source (dereplicated)." if derep else ""
    strata_text = " and ".join({"year": "collection year", "country": "country"}[f] for f in stratify)
    tested = permutations > 0  # --permutations 0 skips every permutation test
    if not tested:
        quality.append("Permutation tests were skipped (--permutations 0), so no p- or q-values are reported and "
                       "nothing is called significant.")

    # 1. Sampling
    counts = A.source_counts(dataset)
    if counts:
        biggest = max(counts, key=lambda r: r["isolates"])
        finding = (f"{len(classified)} isolates have a source class across {len(groups_present)} groups. "
                   f"{biggest['source']} is the largest group ({biggest['isolates']}, {pct(biggest['percent'], 0)}).")
        table = A.counts_by_year(dataset) or counts
        stem = "01_sampling"
        sections.append(_section(
            stem, "Sampling by source", finding, _render(fig_sampling(dataset), stem, figures_dir),
            "Bars count isolates per collection year, coloured by source. Unclassified isolates are not shown. "
            "Uneven sampling between sources limits every comparison below.",
            _table_html(table, [(k, k.replace("_", " ").capitalize()) for k in table[0]]),
            _write_table(tables_dir, stem, table)))
        methods.append(
            f"Isolates were grouped by One Health source ({'L1 class with healthcare environment split out' if by == 'source' else 'paper rollup'}); "
            f"rows sharing a BioSample were counted once ({dataset.duplicates} duplicate rows collapsed).")
    if derep:
        methods.append(
            f"Lineage and feature analyses used one representative per genomic cluster per source: isolates from the "
            f"same source within {threshold:g} SNPs were joined by single linkage and the isolate with the lowest "
            f"identifier was kept ({derep.isolates_before} to {derep.isolates_after} isolates; "
            f"{derep.clusters_collapsed} clusters collapsed), following the dereplication in Vezina et al. "
            f"{refs.cite('vezina')}. Isolates from different sources were never merged, and linkage analyses used all isolates.")
    if stratify and tested:
        strata = A.strata_summary(population.isolates, stratify)
        methods.append(
            f"All permutation tests shuffled source labels only among isolates with the same {strata_text} "
            f"({strata['strata']} strata). {strata['isolates_in_single_source_strata']} of {strata['isolates']} "
            "classified isolates were in strata with a single source, so their labels could not change.")
        if strata["isolates"] and strata["isolates_in_single_source_strata"] > strata["isolates"] / 2:
            quality.append(f"Most classified isolates are in {strata_text} strata that contain a single source, so the "
                           "stratified permutation tests have little power.")

    # 2. Lineages
    typed = [iso for iso in population.classified if iso.lineage]
    if population.lineage_column and typed:
        div = A.diversity(population)
        rare = A.rarefaction_curves(population)
        lineage_rows = A.lineage_sources(population)
        summary = A.lineage_sharing_summary(population)
        enrichment = A.lineage_enrichment(population, permutations=permutations, seed=seed, strata=stratify)
        significant = [r for r in enrichment if r["q"] < 0.05]
        enriched = {r["lineage"] for r in significant if r["direction"] == "enriched"}
        share = 100 * summary["shared_lineages"] / max(1, summary["lineages"])
        held = 100 * summary["isolates_in_shared"] / max(1, summary["typed_isolates"])
        kpis.append(("Lineages", f"{summary['lineages']:,}"))
        kpis.append(("Lineages in ≥2 sources", f"{summary['shared_lineages']} ({share:.0f}%)"))
        r_hg, r_gr, r_hu = (refs.cite(k) for k in ("hunter", "grundmann", "hurlbert"))
        r_th = refs.cite("thorpe") if tested else ""  # cited here only for the enrichment test
        r_he = refs.cite("hetland")
        r_ps = refs.cite("phipson") if tested else ""
        stem = "02_diversity"
        sections.append(_section(
            stem, "Lineage diversity",
            _diversity_finding(div),
            _render(fig_diversity(div), stem, figures_dir),
            f"Simpson's index D is the chance that two isolates drawn from a source belong to different lineages "
            f"{r_hg}; lines are 95% confidence intervals {r_gr}.{derep_note}",
            _table_html(div, [("source", "Source"), ("isolates", "Isolates"), ("lineages", "Lineages"),
                              ("simpson", "Simpson D"), ("ci_low", "95% CI low"), ("ci_high", "95% CI high")]),
            _write_table(tables_dir, stem, div)))
        stem = "03_rarefaction"
        sections.append(_section(
            stem, "Lineage accumulation",
            "Curves that are still rising have not reached saturation, so more sampling would find more lineages.",
            _render(fig_rarefaction(rare, groups_present), stem, figures_dir),
            f"Expected number of lineages in random subsamples of each source (rarefaction {r_hu}). Compare sources at "
            f"the same x value; raw lineage counts favour the larger sample.{derep_note}",
            _table_html(rare, [("source", "Source"), ("isolates_sampled", "Isolates sampled"),
                               ("expected_lineages", "Expected lineages")]),
            _write_table(tables_dir, stem, rare)))
        stem = "04_lineages"
        sections.append(_section(
            stem, "Common lineages and where they were found",
            f"{summary['shared_lineages']} of {summary['lineages']} lineages ({share:.0f}%) were found in two or more "
            f"sources; they hold {held:.0f}% of typed isolates. "
            + ("Enrichment was not tested (--permutations 0)." if not tested else
               f"{len(significant)} lineage–source combinations differ from chance (q < 0.05; * on the axis marks enrichment)."
               if significant else "No lineage was significantly enriched in a source."),
            _render(fig_lineages(lineage_rows, groups_present, enriched), stem, figures_dir),
            "The 20 most common lineages, split by source."
            + (f" Enrichment uses {permutations:,} permutations of source labels with Benjamini–Hochberg correction "
               f"{r_th}{r_ps}." if tested else "")
            + derep_note,
            _table_html(enrichment or lineage_rows,
                        [("lineage", "Lineage"), ("source", "Source"), ("observed", "Observed"), ("expected", "Expected"),
                         ("ratio", "Obs/Exp"), ("p", "p"), ("q", "q (BH)"), ("direction", "Direction")]
                        if enrichment else [("lineage", "Lineage"), ("isolates", "Isolates"), ("sources", "Sources")]),
            _write_table(tables_dir, "04_lineage_sources", lineage_rows),
            extra_downloads=[_write_table(tables_dir, "04_lineage_enrichment", enrichment)] if enrichment else []))
        inter = A.lineage_intersections(population)
        stem = "05_lineage_overlap"
        sections.append(_section(
            stem, "Lineage overlap between sources",
            "Each column is a combination of sources; the bar counts lineages found in exactly that combination.",
            _render(fig_upset(inter, groups_present), stem, figures_dir),
            f"An UpSet-style view of lineage sharing, the overlap measure reported by Hetland et al. {r_he}. "
            f"Sharing a lineage does not imply recent transmission; see the genomic linkage sections.{derep_note}",
            _table_html(inter, [("sources", "Sources"), ("lineages", "Lineages")]),
            _write_table(tables_dir, stem, inter)))
        methods.append(
            f"Lineage diversity per source was measured with Simpson's index {r_hg} and 95% confidence intervals {r_gr}, "
            f"and compared with rarefaction curves {r_hu}. "
            + (f"Lineage–source associations were tested by permuting source labels ({permutations:,} permutations, "
               f"seed {seed}), with two-sided p-values computed as (b+1)/(m+1) {r_ps} and Benjamini–Hochberg "
               f"correction, following Thorpe et al. {r_th}." if tested else
               "Lineage–source associations were not tested (--permutations 0)."))
    else:
        quality.append("No lineage column was found (ST, lineage, cluster ...), so lineage analyses were skipped. "
                       "Add one to the table or pass --lineages.")

    # 3. Features
    if features_path:
        presence, classes = read_features(features_path)
        result = A.feature_prevalence(population, presence, classes, permutations=permutations, seed=seed,
                                      strata=stratify)
        if result.unmatched:
            quality.append(f"{len(result.unmatched)} isolates in the features file matched no classified isolate "
                           f"(e.g. {', '.join(result.unmatched[:3])}).")
        fig = fig_features(result) if result.rows else None
        if fig:
            sig = [t for t in result.tests if t["q"] < 0.05]
            r_wi = refs.cite("wilson")
            r_ch = refs.cite("hetland") if tested else ""
            stem = "06_features"
            sections.append(_section(
                stem, "Resistance, virulence and plasmid features by source",
                (f"{len(sig)} of {len([t for t in result.tests if not math.isnan(t['p'])])} variable features differ in "
                 "prevalence between sources (q < 0.05, marked *). " if tested else
                 "Differences in prevalence between sources were not tested (--permutations 0). ")
                + "Features carried by every isolate are omitted from the chart.",
                _render(fig, stem, figures_dir),
                f"Cell shade is the percentage of genotyped isolates from that source carrying the feature; hover or "
                f"open the table for counts and Wilson 95% intervals {r_wi}."
                + ((f" Tests are permutation chi-squared tests of independence from source {r_ch}, BH-corrected."
                    f"{derep_note} "
                    + ("Dereplication removes outbreak clusters, but lineages can still share features; a GWAS that "
                       "models population structure (pyseer) is the stricter check." if derep else
                       "Isolates from one clone are not independent, so these tests overstate significance when clones "
                       "dominate; rerun with --dereplicate, or use a GWAS that models population structure (pyseer)."))
                   if tested else derep_note),
                _table_html(result.tests, [("feature", "Feature"), ("class", "Class"),
                                           ("isolates_with_feature", "Isolates"), ("percent_overall", "% overall"),
                                           ("chi2", "Chi-squared"), ("p", "p"), ("q", "q (BH)")]),
                _write_table(tables_dir, "06_feature_prevalence", result.rows),
                extra_downloads=[_write_table(tables_dir, "06_feature_tests", result.tests)]))
            methods.append(
                f"Feature prevalence per source is reported with Wilson 95% intervals {r_wi}; association with source "
                + (f"was tested with a permutation chi-squared test ({permutations:,} permutations) and BH correction, "
                   f"as in Hetland et al. {r_ch}." if tested else "was not tested (--permutations 0).")
                + " Only isolates present in the features file count as genotyped "
                f"({', '.join(f'{g}: {n}' for g, n in result.genotyped.items())}).")

    # 4. Genomic linkage
    if distances_path:
        link = A.linkage(dataset, names, matrix, threshold=threshold, permutations=permutations, seed=seed,
                         strata=stratify)
        if link.unmatched:
            quality.append(f"{len(link.unmatched)} names in the distance matrix matched no isolate "
                           f"(e.g. {', '.join(link.unmatched[:3])}).")
        between = [r for r in link.events if r["kind"] == "between"]
        within = [r for r in link.events if r["kind"] == "within"]
        n_between = sum(r["events"] for r in between)
        n_within = sum(r["events"] for r in within)
        kpis.append((f"Cross-source links (≤{threshold:g} SNPs)", f"{n_between}"))
        r_th, r_he, r_da = (refs.cite(k) for k in ("thorpe", "hetland", "david"))
        r_mu, r_ps = (refs.cite("muloi"), refs.cite("phipson")) if tested else ("", "")
        close_rows = [
            {"pairs": label, "all_pairs": total, f"within_{threshold:g}_snps": int((values <= threshold).sum()),
             "within_100_snps": int((values <= 100).sum())}
            for label, total, values in (("Same source", link.pairs_within, link.close_within),
                                         ("Different sources", link.pairs_between, link.close_between))
        ]
        stem = "07_distances"
        sections.append(_section(
            stem, "How close are isolates from the same and different sources?",
            f"{int((link.close_within <= threshold).sum())} same-source pairs and "
            f"{int((link.close_between <= threshold).sum())} cross-source pairs are within {threshold:g} SNPs.",
            _render(fig_distance_histogram(link), stem, figures_dir),
            f"Pairwise SNP distances up to 100 SNPs. The vertical line is the linkage threshold. Thresholds are "
            f"organism- and pipeline-specific: 21 SNPs separated hospital clusters of K. pneumoniae {r_da} and 22 was "
            f"optimal in a Norwegian One Health study {r_he}.",
            _table_html(close_rows, [(k, k.replace("_", " ").capitalize()) for k in close_rows[0]]),
            _write_table(tables_dir, stem, close_rows)))
        stem = "08_threshold_sweep"
        sections.append(_section(
            stem, "Sensitivity to the SNP threshold",
            f"At {threshold:g} SNPs there are {n_within} within-source and {n_between} between-source sharing events. "
            "A conclusion that holds across thresholds is more robust.",
            _render(fig_sweep(link), stem, figures_dir),
            f"Events are counted once per source pair per cluster {r_th}. Counts can fall at large thresholds because "
            "clusters merge.",
            _table_html(link.sweep, [("threshold", "SNP threshold"), ("within_events", "Within-source events"),
                                     ("between_events", "Between-source events")]),
            _write_table(tables_dir, stem, link.sweep),
            # Normally in the observed-vs-expected section, which needs permutations.
            extra_downloads=[] if tested else [_write_table(tables_dir, "08_sharing_events", link.events)]))
        stem = "09_linked_isolates"
        sections.append(_section(
            stem, "Which sources share close relatives?",
            "Each cell is the share of isolates from the row source with at least one close relative from the column source.",
            _render(fig_linked(link), stem, figures_dir),
            f"Close relative = within {threshold:g} SNPs, excluding isolates from the same host individual {r_he}. "
            "The diagonal is within-source sharing.",
            _table_html(link.linked_fraction, [("from_source", "From"), ("to_source", "To"), ("isolates", "Isolates"),
                                               ("with_close_relative", "With close relative"), ("percent", "%")]),
            _write_table(tables_dir, stem, link.linked_fraction)))
        oe = fig_observed_expected(link) if tested else None
        if oe:
            stem = "10_observed_expected"
            sig = [r for r in between if r["q"] < 0.05]
            lower = [r for r in sig if r["events"] < r["expected"]]
            higher = [r for r in sig if r["events"] > r["expected"]]
            if sig:
                finding = (f"{len(sig)} source pairs differ from chance (q < 0.05): {len(lower)} share fewer close "
                           f"relatives than expected and {len(higher)} share more.")
                if lower and not higher:
                    finding += " That pattern fits transmission that happens mostly within sources."
            else:
                finding = "No source pair differs from the permutation expectation."
            sections.append(_section(
                stem, "Is cross-source sharing more or less than chance?",
                finding,
                _render(oe, stem, figures_dir),
                f"Shaded bands show the 95% range of events when source labels are shuffled over the same genomes "
                f"({permutations:,} permutations), which keeps sampling effort fixed {r_mu}{r_th}. * marks q < 0.05 {r_ps}. "
                + (f"Labels were shuffled only among isolates with the same {strata_text}, so differences in when or "
                   "where sources were sampled do not drive the expectation." if stratify else
                   "Caveat: shuffling ignores when and where samples were taken, so sources sampled at different times "
                   "or places also show fewer links than chance; rerun with --stratify-by year,country."),
                _table_html(link.events, [("source_a", "Source A"), ("source_b", "Source B"), ("events", "Events"),
                                          ("events_per_100", "Per 100 isolates"), ("expected", "Expected"),
                                          ("expected_low", "2.5%"), ("expected_high", "97.5%"), ("q", "q (BH)")]),
                _write_table(tables_dir, stem, link.events)))
        cl = fig_clusters(link)
        if cl:
            stem = "11_cross_source_clusters"
            cross = [c for c in link.clusters if c["n_sources"] >= 2]
            gaps = [e["years_apart"] for e in link.cross_edges if not math.isnan(e["years_apart"])]
            timing = f" Linked cross-source isolates were collected a median {np.median(gaps):.1f} years apart." if gaps else ""
            sections.append(_section(
                stem, "Clusters that span sources",
                f"{len(cross)} clusters contain isolates from two or more sources.{timing}",
                _render(cl, stem, figures_dir),
                "Each panel is one cluster; nodes are isolates (shape and colour = source). Lines show the minimum "
                "spanning tree of pairs within the threshold, labelled with SNP distance. Shared clusters are candidates "
                "for follow-up, not proof of direction or route.",
                _table_html(link.cross_edges, [("cluster", "Cluster"), ("isolate_a", "Isolate A"), ("source_a", "Source A"),
                                               ("isolate_b", "Isolate B"), ("source_b", "Source B"),
                                               ("lineage", "Lineage"), ("snps", "SNPs"), ("years_apart", "Years apart")]),
                _write_table(tables_dir, "11_cross_source_pairs", link.cross_edges),
                extra_downloads=[_write_table(tables_dir, "11_clusters", link.clusters)]))
        swept = ", ".join(f"{row['threshold']:g}" for row in link.sweep)
        methods.append(
            f"Pairwise SNP distances were read for {link.matched} isolates. Isolates within {threshold:g} SNPs were "
            f"joined by single linkage; each unique source pair connected inside a cluster counted as one sharing "
            f"event, excluding pairs from the same host individual {r_th}{r_he}. Events were normalised per 100 isolates "
            "of the two sources"
            + (f" and compared with {permutations:,} permutations of source labels over the fixed genetic network "
               f"{r_mu}, with two-sided p-values {r_ps} and BH correction" if tested else "")
            + f". Results were repeated at {swept} SNPs.")

    if dataset.unclassified:
        quality.append(f"{dataset.unclassified} isolates have no source class and are excluded from source comparisons.")
    if dataset.duplicates:
        quality.append(f"{dataset.duplicates} rows repeated an isolate already listed (same BioSample or ID) and were "
                       "collapsed.")
    if dataset.unnamed:
        quality.append(f"{dataset.unnamed} rows have no identifier (BioSample, ID, strain or sample name). Each is "
                       "counted as its own isolate, but cannot be matched to lineages, features or distances.")

    page = _page(title, note, kpis, sections, quality, methods, refs,
                 dict(threshold=threshold, by=by, permutations=permutations, seed=seed, dereplicate=bool(derep),
                      stratify=stratify))
    path = out_dir / "report.html"
    path.write_text(page, encoding="utf-8")
    return path


def _diversity_finding(rows: list[dict]) -> str:
    """Headline for the diversity section. Sources with one typed isolate have no D and are skipped."""
    measured = [r for r in rows if not math.isnan(r["simpson"])]
    if not measured:
        return "No source has two or more typed isolates, so lineage diversity cannot be measured."
    top = max(measured, key=lambda r: r["simpson"])
    if len(measured) == 1:
        return (f"Only {top['source']} has two or more typed isolates (D = {top['simpson']:.2f}), so diversity "
                "cannot be compared between sources.")
    return (f"Lineage diversity was highest in {top['source']} (D = {top['simpson']:.2f}). "
            "Overlapping intervals mean the sources are not clearly different.")


def _section(stem, heading, finding, svg, caption, table, download, extra_downloads=()) -> str:
    links = " · ".join(
        f'<a href="{html.escape(p)}">{html.escape(Path(p).name)}</a>' for p in (download, *extra_downloads)
    )
    figure_links = f'<a href="figures/{stem}.svg">SVG</a> · <a href="figures/{stem}.png">PNG</a>'
    return f"""
<section class="card" id="{stem}">
  <h2 id="{stem}-title">{html.escape(heading)}</h2>
  <p class="finding">{html.escape(finding)}</p>
  <figure>{svg}<figcaption>{html.escape(caption)}</figcaption></figure>
  <details><summary>Data table</summary>{table}</details>
  <p class="downloads">Download: {figure_links} · {links}</p>
</section>"""


class References:
    LIBRARY = {
        "vezina": "Vezina B, Morampalli BR, Nguyen HA, Gomez-Simmonds A, Peleg AY, Macesic N (2025) The rise and global spread of IMP carbapenemases (1996–2023): a genomic epidemiology study. Nat Commun 17:183. https://doi.org/10.1038/s41467-025-66874-7",
        "thorpe": "Thorpe HA et al. (2022) A large-scale genomic snapshot of Klebsiella spp. isolates in Northern Italy reveals limited transmission between clinical and non-clinical settings. Nat Microbiol 7:2054–2067. https://doi.org/10.1038/s41564-022-01263-0",
        "hetland": "Hetland MAK et al. (2025) A genome-wide One Health study of Klebsiella pneumoniae in Norway reveals overlapping populations but few recent transmission events across reservoirs. Genome Med 17:42. https://doi.org/10.1186/s13073-025-01466-0",
        "muloi": "Muloi DM et al. (2022) Population genomics of Escherichia coli in livestock-keeping households across a rapidly developing urban landscape. Nat Microbiol 7:581–589. https://doi.org/10.1038/s41564-022-01079-y",
        "david": "David S et al. (2019) Epidemic of carbapenem-resistant Klebsiella pneumoniae in Europe is driven by nosocomial spread. Nat Microbiol 4:1919–1929. https://doi.org/10.1038/s41564-019-0492-8",
        "hunter": "Hunter PR, Gaston MA (1988) Numerical index of the discriminatory ability of typing systems: an application of Simpson's index of diversity. J Clin Microbiol 26:2465–2466. https://doi.org/10.1128/jcm.26.11.2465-2466.1988",
        "grundmann": "Grundmann H, Hori S, Tanner G (2001) Determining confidence intervals when measuring genetic diversity and the discriminatory abilities of typing methods for microorganisms. J Clin Microbiol 39:4190–4192. https://doi.org/10.1128/JCM.39.11.4190-4192.2001",
        "hurlbert": "Hurlbert SH (1971) The nonconcept of species diversity: a critique and alternative parameters. Ecology 52:577–586. https://doi.org/10.2307/1934145",
        "phipson": "Phipson B, Smyth GK (2010) Permutation p-values should never be zero. Stat Appl Genet Mol Biol 9:Article39. https://doi.org/10.2202/1544-6115.1585",
        "wilson": "Wilson EB (1927) Probable inference, the law of succession, and statistical inference. J Am Stat Assoc 22:209–212. https://doi.org/10.1080/01621459.1927.10502953",
    }

    def __init__(self) -> None:
        self.order: list[str] = []

    def cite(self, key: str) -> str:
        if key not in self.order:
            self.order.append(key)
        return f"[{self.order.index(key) + 1}]"

    def html(self) -> str:
        items = "".join(f"<li>{html.escape(self.LIBRARY[k])}</li>" for k in self.order)
        return f"<ol class='refs'>{items}</ol>"


def _page(title, note, kpis, sections, quality, methods, refs, params) -> str:
    light = "\n".join(f"  --{name}: {pair[0]};" for name, pair in TOKENS.items())
    dark = "\n".join(f"  --{name}: {pair[1]};" for name, pair in TOKENS.items())
    tiles = "".join(f'<div class="tile"><div class="label">{html.escape(k)}</div><div class="value">{html.escape(v)}</div></div>' for k, v in kpis)
    toc = "".join(
        f'<a href="#{m.group(1)}">{html.escape(html.unescape(m.group(2)))}</a>'
        for s in sections for m in [re.search(r'id="(\d\d_[a-z_]+)-title">([^<]+)<', s)] if m
    )
    banner = f'<p class="banner">{html.escape(note)}</p>' if note else ""
    issues = "".join(f"<li>{html.escape(q)}</li>" for q in quality) or "<li>No data-quality issues were detected.</li>"
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    param_text = (f"grouping: {params['by']}; SNP threshold: {params['threshold']:g}; permutations: "
                  f"{params['permutations']:,}; random seed: {params['seed']}; "
                  f"dereplicated: {'yes' if params['dereplicate'] else 'no'}; "
                  f"permutation strata: {', '.join(params['stratify']) or 'none'}")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{
  color-scheme: light;
{light}
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
  color-scheme: dark;
{dark}
  }}
}}
:root[data-theme="dark"] {{
  color-scheme: dark;
{dark}
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--page); color: var(--ink); font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }}
main {{ max-width: 1040px; margin: 0 auto; padding: 32px 16px 64px; }}
header h1 {{ font-size: 28px; line-height: 1.2; margin: 0 0 6px; font-weight: 650; }}
.sub {{ color: var(--ink-2); margin: 0 0 16px; }}
.banner {{ border-left: 3px solid var(--accent); background: var(--surface); padding: 10px 14px; border-radius: 6px; color: var(--ink-2); }}
.kpis {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin: 20px 0; }}
.tile {{ background: var(--surface); border: 1px solid var(--grid); border-radius: 12px; padding: 14px 16px; }}
.tile .label {{ color: var(--ink-2); font-size: 13px; }}
.tile .value {{ font-size: 24px; font-weight: 600; margin-top: 2px; }}
nav.toc {{ display: flex; flex-wrap: wrap; gap: 6px 14px; margin: 8px 0 24px; font-size: 14px; }}
a {{ color: var(--s-human); text-underline-offset: 2px; }}
.card {{ background: var(--surface); border: 1px solid var(--grid); border-radius: 14px; padding: 24px; margin: 18px 0; }}
.card h2 {{ font-size: 19px; margin: 0 0 6px; font-weight: 620; }}
.finding {{ margin: 0 0 14px; font-size: 16px; }}
figure {{ margin: 0; }}
figure svg {{ display: block; max-width: 100%; height: auto; margin: 0 auto; }}
figcaption, .muted {{ color: var(--ink-2); font-size: 13.5px; margin-top: 10px; }}
.small {{ font-size: 12.5px; }}
details {{ margin-top: 12px; }}
summary {{ cursor: pointer; color: var(--ink-2); font-size: 14px; }}
.table-wrap {{ overflow-x: auto; margin-top: 8px; }}
table {{ border-collapse: collapse; font-size: 13px; width: 100%; font-variant-numeric: tabular-nums; }}
th, td {{ text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--grid); white-space: nowrap; }}
th {{ color: var(--ink-2); font-weight: 600; }}
.downloads {{ font-size: 13px; color: var(--ink-2); margin: 12px 0 0; }}
.tip:focus {{ outline: 2px solid var(--accent); outline-offset: 1px; }}
.tip:hover {{ opacity: 0.82; }}
#tooltip {{ position: fixed; pointer-events: none; z-index: 10; max-width: 320px; background: var(--surface);
  color: var(--ink); border: 1px solid var(--grid); border-radius: 8px; padding: 8px 10px; font-size: 13px;
  box-shadow: 0 4px 18px rgba(0,0,0,0.12); display: none; }}
ol.refs {{ font-size: 13px; color: var(--ink-2); padding-left: 20px; }}
ul.quality {{ color: var(--ink-2); font-size: 14px; }}
@media print {{ .card {{ break-inside: avoid; }} #tooltip {{ display: none !important; }} }}
</style>
</head>
<body>
<main>
<header>
  <h1>{html.escape(title)}</h1>
  <p class="sub">onehealth {__version__} · generated {stamp} · {html.escape(param_text)}</p>
  {banner}
</header>
<div class="kpis">{tiles}</div>
<nav class="toc">{toc}<a href="#methods">Methods</a></nav>
{''.join(sections)}
<section class="card" id="quality">
  <h2>Data quality</h2>
  <ul class="quality">{issues}</ul>
</section>
<section class="card" id="methods">
  <h2>Methods</h2>
  <p>{' '.join(methods)}</p>
  <p class="muted">Every figure is also saved as SVG and PNG under figures/, and every number under tables/.
  Shared lineages or close genetic relatedness are evidence of a connection, not of its direction or route.</p>
  <h2>References</h2>
  {refs.html()}
</section>
</main>
<div id="tooltip" role="tooltip"></div>
<script>
(() => {{
  const tip = document.getElementById('tooltip');
  const show = (el, x, y) => {{
    tip.textContent = el.getAttribute('data-tip');
    tip.style.display = 'block';
    const w = tip.offsetWidth, h = tip.offsetHeight;
    tip.style.left = Math.min(x + 14, window.innerWidth - w - 8) + 'px';
    tip.style.top = Math.max(8, y - h - 12) + 'px';
  }};
  document.addEventListener('pointermove', (e) => {{
    const el = e.target.closest && e.target.closest('.tip');
    if (el) show(el, e.clientX, e.clientY); else tip.style.display = 'none';
  }});
  document.addEventListener('focusin', (e) => {{
    const el = e.target.closest && e.target.closest('.tip');
    if (!el) return;
    const r = el.getBoundingClientRect();
    show(el, r.left + r.width / 2, r.top);
  }});
  document.addEventListener('focusout', () => {{ tip.style.display = 'none'; }});
}})();
</script>
</body>
</html>
"""
