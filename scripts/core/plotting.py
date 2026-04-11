"""Shared plotting utilities — matching output/ analysis style.

Style:
- Arial 10pt base, 11pt titles
- All four spines visible
- Ticks outward
- Colorblind-friendly palette
- DPI 200
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Color palette ─────────────────────────────────────────────────────────────
PALETTE = {
    "blue":   "#2166AC",
    "red":    "#B2182B",
    "green":  "#4CAF50",
    "yellow": "#CCBB44",
    "cyan":   "#66CCEE",
    "purple": "#AA3377",
    "grey":   "#BBBBBB",
    "dark":   "#332288",
    "light_grey": "#DDDDDD",
    "orange": "#FF9800",
}

# Column widths (inches)
COL1 = 7       # single panel
COL1_5 = 10    # 1.5 column
COL2 = 14      # double column

# ── rcParams ─────────────────────────────────────────────────────────────────
STYLE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "regular",
    "axes.labelsize": 10,
    "axes.labelweight": "regular",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "legend.title_fontsize": 9,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.pad": 3,
    "ytick.major.pad": 3,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.grid": False,
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "legend.frameon": False,
    "legend.borderpad": 0.3,
    "legend.handlelength": 1.2,
    "legend.handletextpad": 0.4,
    "lines.linewidth": 1.5,
    "lines.markersize": 6,
    "patch.linewidth": 0.6,
    "boxplot.flierprops.markersize": 3,
    "boxplot.flierprops.markeredgewidth": 0.5,
    "boxplot.boxprops.linewidth": 0.8,
    "boxplot.whiskerprops.linewidth": 0.8,
    "boxplot.capprops.linewidth": 0.8,
    "boxplot.medianprops.linewidth": 1.0,
}


def apply_style(cfg):
    """Apply shared analysis style globally."""
    plt.rcdefaults()
    plt.rcParams.update(STYLE_RC)
    plot_cfg = cfg.get("plot", {})
    dpi = plot_cfg.get("dpi", 200)
    plt.rcParams.update({"figure.dpi": dpi, "savefig.dpi": dpi})


def style_axes(ax):
    """Finalize axes styling."""
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.tick_params(axis="both", direction="out", width=0.8, length=3.5, pad=3)


def get_recipe_colors(cfg):
    """Condition color mapping."""
    default_map = {
        "R1": PALETTE["blue"],
        "R2": PALETTE["red"],
    }
    colors = dict(cfg.get("plot", {}).get("recipe_colors", default_map))
    conditions = cfg.get("experiment", {}).get("conditions", [])
    fallback = [
        PALETTE["blue"],
        PALETTE["red"],
        PALETTE["green"],
        PALETTE["purple"],
        PALETTE["cyan"],
        PALETTE["yellow"],
        PALETTE["dark"],
    ]
    for idx, cond in enumerate(conditions):
        colors.setdefault(cond, fallback[idx % len(fallback)])
    return colors


def get_color(name):
    """Get named color from palette."""
    return PALETTE.get(name, name)


def get_plot_format(cfg=None):
    """Return normalized plot format extension without leading dot."""
    if not cfg:
        return "png"
    fmt = str(cfg.get("plot", {}).get("format", "png")).strip().lower()
    return fmt.lstrip(".") or "png"


def plot_filename(cfg, filename):
    """Return a filename with the configured plot extension."""
    fmt = get_plot_format(cfg)
    path = Path(filename)
    if path.suffix:
        return str(path.with_suffix(f".{fmt}"))
    return f"{path}.{fmt}"


def plot_path(cfg, directory, filename):
    """Return an output path under directory with the configured plot extension."""
    return os.path.join(directory, plot_filename(cfg, filename))


def save_figure(fig, path, cfg=None):
    """Save figure at publication quality."""
    dpi = 200
    fmt = get_plot_format(cfg)
    output_path = plot_filename(cfg, os.fspath(path))
    if cfg:
        dpi = cfg.get("plot", {}).get("dpi", 200)
    fig.savefig(output_path, format=fmt, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none", pad_inches=0.05)
    plt.close(fig)
    return output_path
