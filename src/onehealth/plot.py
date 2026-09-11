"""Simple count plots. No cluster calling."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def plot_counts(rows: list[dict[str, str]], dest: Path, by: str = "l1") -> None:
    if by not in {"l1", "paper_rollup"}:
        raise ValueError("plot --by must be l1 or paper_rollup")
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    _bar(axes[0], rows, "l1", "L1 source class")
    _bar(axes[1], rows, "paper_rollup", "Paper rollup")
    if by == "paper_rollup":
        axes[1].set_title("Paper rollup (requested)")
    fig.tight_layout()
    fig.savefig(dest, dpi=120)
    plt.close(fig)


def _bar(axis, rows: list[dict[str, str]], field: str, title: str) -> None:
    counter = Counter((row.get(field) or "Unclassified") for row in rows)
    labels = list(counter.keys())
    values = [counter[label] for label in labels]
    axis.bar(range(len(labels)), values, color="#3b6d8c")
    axis.set_xticks(range(len(labels)), labels, rotation=30, ha="right")
    axis.set_ylabel("Isolates")
    axis.set_title(title)
    if not labels:
        axis.text(0.5, 0.5, "No rows", ha="center", va="center", transform=axis.transAxes)
