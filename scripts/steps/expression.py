"""Step 2a: Detected Genes — polar bar chart of detected genes per sample by expression level."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core.config_data import (
    get_filtered_expr,
    get_norm_method,
    get_sample_display_name,
    get_sample_manifest,
    get_sample_names,
)
from ..core.config_runtime import ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, save_figure, plot_filename

log = setup_logger("sample_analysis.detected_genes")

# Expression bins: (label_template, min, max, color)
# {u} is replaced by the normalization unit (CPM or FPKM) at runtime
_BIN_DEFS = [
    ("0-1 {u} (low)", 0, 1, "#d1e5f0"),
    ("1-5 {u}", 1, 5, "#fddbc7"),
    ("5-10 {u}", 5, 10, "#f4a582"),
    ("10-100 {u}", 10, 100, "#d6604d"),
    (">=100 {u}", 100, np.inf, "#b2182b"),
]


def run(cfg):
    log.info("── Step 2a: Expression ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    out_dir = ensure_output_dir(cfg, "02_sample_analysis/expression")
    norm_label = get_norm_method(cfg).upper()
    bins = [(t.format(u=norm_label), lo, hi, c) for t, lo, hi, c in _BIN_DEFS]

    # Load filtered expression
    cpm = _load_expr(cfg, samples)
    all_samples = samples["all"]
    n_samples = len(all_samples)

    # Count genes in each expression bin per sample
    bin_counts = {}
    for label, lo, hi, _ in bins:
        counts = []
        for s in all_samples:
            if hi == np.inf:
                counts.append((cpm[s] >= lo).sum())
            else:
                counts.append(((cpm[s] >= lo) & (cpm[s] < hi)).sum())
        bin_counts[label] = np.array(counts)

    # Total detected genes (expr >= 1) per sample
    detected = np.array([(cpm[s] >= 1).sum() for s in all_samples])

    # Save TSV
    tsv_data = {"sample": all_samples, "total_genes": [cpm.shape[0]] * n_samples,
                "detected_expr1": detected}
    for label, _, _, _ in bins:
        tsv_data[label] = bin_counts[label]
    pd.DataFrame(tsv_data).to_csv(
        os.path.join(out_dir, "detected_genes.tsv"), sep="\t", index=False)

    # ── Polar bar chart ──
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": "polar"})

    # Start from top, go clockwise
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    angles = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)
    bar_width = 2 * np.pi / n_samples * 0.8

    # Stacked bars from inside out
    bottoms = np.zeros(n_samples)
    for label, _, _, color in bins:
        vals = bin_counts[label]
        ax.bar(angles, vals, width=bar_width, bottom=bottoms,
               color=color, edgecolor="white", linewidth=0.5, label=label)
        bottoms += vals

    # Labels: sample name outside, total detected at bar top
    ax.set_xticks(angles)
    ax.set_xticklabels(
        [get_sample_display_name(cfg, s) for s in all_samples],
        fontsize=9,
    )

    for i, angle in enumerate(angles):
        ax.text(angle, bottoms[i] + 100, str(detected[i]),
                ha="center", va="bottom", fontsize=9, color="#555")

    ax.set_title("Detected Genes per Sample by Expression Level",
                 fontsize=12, pad=20)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.15, -0.08),
              fontsize=8, title="Expression interval", title_fontsize=9,
              frameon=True, fancybox=True)

    # Clean up radial ticks
    ax.set_yticks([2000, 4000, 6000])
    ax.set_yticklabels(["2000", "4000", "6000"], fontsize=9, color="#999")
    ax.set_rlabel_position(0)

    detected_genes_plot = plot_filename(cfg, "detected_genes_polar.png")
    save_figure(fig, os.path.join(out_dir, "detected_genes_polar.png"), cfg)

    write_readme(out_dir, "2a", "Expression", {
        "detected_genes.tsv": f"Per-sample gene detection counts across expression-level bins ({norm_label} intervals).",
        detected_genes_plot: "Polar bar chart of detected genes per sample, stacked by expression level.",
    })

    manifest = get_sample_manifest(cfg)
    note = ""
    if manifest is not None:
        note = (
            "\n- Exploration mode: repeated physical samples retained from `sample_manifest.tsv`\n"
        )

    doc = f"""\
### Parameters

- **Detection threshold**: {norm_label} >= 1
- **Total genes in matrix**: {cpm.shape[0]:,}
- **Samples**: {n_samples}
{note}

### Results

- Detected genes ({norm_label} >= 1): {detected.min():,}–{detected.max():,} per sample

{img(cfg, "02_sample_analysis/expression", "detected_genes_polar.png", "Detected genes per sample by expression level")}
"""
    update_section(cfg, "2a", "Expression", doc)

    log.info(f"  Detected genes: {detected.min()}-{detected.max()} ({norm_label}>=1)")
    log.info(f"  Saved {detected_genes_plot}")


def _load_expr(cfg, samples):
    """Load filtered normalized expression (CPM or FPKM)."""
    return get_filtered_expr(cfg, log2=False)[samples["all"]]
