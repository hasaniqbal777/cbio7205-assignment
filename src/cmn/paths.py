"""Shared paths and figure saving."""

from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
FIGS = SRC / "figs"
REPORT_FIGS = SRC.parent / "report" / "fig"
DATA = SRC.parent / "data"


def save_fig(fig, name, **kwargs):
    """Save a matplotlib figure to src/figs/ and a copy to report/fig/ (used by LaTeX)."""
    kwargs.setdefault("bbox_inches", "tight")
    kwargs.setdefault("dpi", 300)
    for d in (FIGS, REPORT_FIGS):
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / name, **kwargs)
