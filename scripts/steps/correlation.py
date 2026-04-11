"""Step 2b: Correlation Analysis — sample-sample correlation heatmap."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core.config_data import (
    get_filtered_expr,
    get_sample_display_name,
    get_sample_manifest,
    get_sample_names,
)
from ..core.config_runtime import get_conditions, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, save_figure, plot_filename, COL1_5

log = setup_logger("sample_analysis.correlation")


def run(cfg):
    corr_cfg = cfg["sample_analysis"]["correlation"]
    if not corr_cfg["enabled"]:
        log.info("Correlation analysis disabled, skipping")
        return

    log.info("── Step 2b: Correlation Analysis ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    out_dir = ensure_output_dir(cfg, "02_sample_analysis/correlation")
    log2cpm = _load_log2expr(cfg, samples)
    method = corr_cfg["method"]

    corr_matrix = log2cpm.corr(method=method)
    corr_matrix.to_csv(os.path.join(out_dir, "correlation_matrix.tsv"), sep="\t")
    log.info(f"  {method.capitalize()} correlation computed")

    # ── Heatmap (upper triangle, sample names on diagonal) ──
    annot_fontsize = corr_cfg.get("annot_fontsize", 8)
    n = len(corr_matrix)
    sample_names = list(corr_matrix.columns)
    vals = corr_matrix.values

    fig, ax = plt.subplots(figsize=(7, 6))
    vmin = vals[vals < 1].min() - 0.005 if (vals < 1).any() else 0.9

    # Draw upper triangle cells (i <= j)
    for i in range(n):
        for j in range(i, n):
            # Color mapping
            norm_val = (vals[i, j] - vmin) / (1 - vmin)
            color = plt.cm.RdYlBu_r(norm_val)

            rect = plt.Rectangle((j, i), 1, 1, facecolor=color,
                                  edgecolor="white", linewidth=0.5)
            ax.add_patch(rect)

            if i == j:
                # Diagonal: sample name (white background)
                rect.set_facecolor("white")
                ax.text(j + 0.5, i + 0.5, get_sample_display_name(cfg, sample_names[i]),
                        ha="center", va="center", fontsize=9, fontweight="bold")
            else:
                # Upper triangle: correlation value
                ax.text(j + 0.5, i + 0.5, f"{vals[i, j]:.3f}",
                        ha="center", va="center", fontsize=annot_fontsize)

    ax.set_xlim(0, n)
    ax.set_ylim(n, 0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    norm_label = cfg.get("normalization", {}).get("method", "cpm").upper()
    ax.set_title(f"Sample Correlation (log2{norm_label}, {method.capitalize()})", fontsize=10, pad=10)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r",
                                norm=plt.Normalize(vmin=vmin, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=20, pad=0.02)
    cb.set_label(f"{method.capitalize()} r", fontsize=10)
    cb.ax.tick_params(labelsize=9)

    correlation_heatmap = plot_filename(cfg, "correlation_heatmap.png")
    save_figure(fig, os.path.join(out_dir, "correlation_heatmap.png"), cfg)
    log.info(f"  Saved {correlation_heatmap}")

    # Summary
    conditions = get_conditions(cfg)
    lines = [f"Correlation Summary ({method})", "=" * 40]

    def _stats_text(values):
        values = np.asarray(values, dtype=float)
        if values.size == 0:
            return "min=NA mean=NA max=NA"
        return f"min={values.min():.4f} mean={values.mean():.4f} max={values.max():.4f}"

    for cond in conditions:
        sub = corr_matrix.loc[samples[cond], samples[cond]]
        vals = sub.values[np.triu_indices_from(sub.values, k=1)]
        lines.append(f"  {cond}: {_stats_text(vals)}")
    between = corr_matrix.loc[samples[conditions[0]], samples[conditions[1]]].values.flatten()
    lines.append(f"  Between: {_stats_text(between)}")
    with open(os.path.join(out_dir, "correlation_summary.txt"), "w") as f:
        f.write("\n".join(lines))

    write_readme(out_dir, "2b", "Correlation Analysis", {
        "correlation_matrix.tsv": f"Sample-by-sample {method} correlation matrix on log2CPM values.",
        correlation_heatmap: f"Upper-triangle heatmap of pairwise {method} correlations with annotated values.",
        "correlation_summary.txt": f"Within-condition and between-condition {method} correlation statistics.",
    })

    manifest = get_sample_manifest(cfg)
    note = ""
    if manifest is not None:
        note = "\n- Exploration mode: repeated physical samples retained from `sample_manifest.tsv`\n"

    doc = f"""\
### Parameters

- **Method**: {method.capitalize()}
- **Input**: log2(CPM + 1) filtered expression matrix
{note}

### Results

- Between-condition correlation: {_stats_text(between)}

{img(cfg, "02_sample_analysis/correlation", "correlation_heatmap.png", f"Sample-sample {method} correlation heatmap")}
"""
    update_section(cfg, "2b", "Correlation Analysis", doc)


def _load_log2expr(cfg, samples):
    return get_filtered_expr(cfg, log2=True)[samples["all"]]
