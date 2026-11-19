"""Shared matplotlib styling (dataviz reference palette, recessive axes)."""

import matplotlib.pyplot as plt

# Categorical slots of the dataviz reference palette, in fixed order.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK = "#2b2b2b"
MUTED = "#6b6b6b"
GRID = "#e6e6e3"


def axes(xlabel, ylabel, figsize=(5, 3.4)):
    """A single axis with thin, recessive spines and grid."""
    fig, ax = plt.subplots(figsize=figsize)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=INK, fontsize=10)
    ax.set_ylabel(ylabel, color=INK, fontsize=10)
    return fig, ax


def bars(ax, x, heights, color=SERIES[0], labels=True, fmt="{:.3f}"):
    """Thin bars with a value label on top of each."""
    b = ax.bar(x, heights, width=0.6, color=color, edgecolor="white", linewidth=1)
    if labels:
        for xi, h in zip(x, heights):
            ax.annotate(fmt.format(h), (xi, h), xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8, color=INK)
    return b
