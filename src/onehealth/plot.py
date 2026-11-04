"""Simple count plots. No cluster calling."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from onehealth.constants import L1_VALUES, PAPER_ROLLUP_VALUES  # noqa: E402
from onehealth.errors import InputError  # noqa: E402


def plot_counts(rows: list[dict[str, str]], dest: Path) -> None:
    """Two panels: L1 and paper rollup. Categories keep a fixed order, zeros included."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    kind = dest.suffix.lstrip(".").lower()
    if kind and kind not in fig.canvas.get_supported_filetypes():
        plt.close(fig)
        raise InputError(f"{dest}: cannot write .{kind} images; use .png, .svg or .pdf")
    dest.parent.mkdir(parents=True, exist_ok=True)
    _bar(axes[0], rows, "l1", L1_VALUES, "L1 source class")
    _bar(axes[1], rows, "paper_rollup", PAPER_ROLLUP_VALUES, "Paper rollup")
    fig.tight_layout()
    fig.savefig(dest, dpi=120)
    plt.close(fig)


def _bar(axis, rows: list[dict[str, str]], field: str, order: tuple[str, ...], title: str) -> None:
    counter = Counter((row.get(field) or "Unclassified") for row in rows)
    labels = list(order) + sorted(set(counter) - set(order))
    values = [counter.get(label, 0) for label in labels]
    axis.bar(range(len(labels)), values, color="#3b6d8c")
    axis.set_xticks(range(len(labels)), labels, rotation=30, ha="right")
    axis.set_ylabel("Rows")
    axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    axis.set_title(title)
    if not rows:
        axis.text(0.5, 0.5, "No rows", ha="center", va="center", transform=axis.transAxes)
