"""Step 2c: PCA — Principal Component Analysis on filtered log2CPM."""

import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

from ..core.config_data import (
    get_filtered_expr,
    get_sample_condition,
    get_sample_display_name,
    get_sample_manifest,
    get_sample_names,
)
from ..core.config_runtime import get_conditions, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, save_figure, get_recipe_colors, plot_filename

log = setup_logger("sample_analysis.pca")


def run(cfg):
    pca_cfg = cfg["sample_analysis"]["pca"]
    if not pca_cfg["enabled"]:
        log.info("PCA disabled, skipping")
        return

    log.info("── Step 2c: PCA ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    out_dir = ensure_output_dir(cfg, "02_sample_analysis/pca")
    colors = get_recipe_colors(cfg)
    conditions = get_conditions(cfg)
    strain = cfg["project"].get("strain", "")
    manifest = get_sample_manifest(cfg)

    log2cpm = _load_log2expr(cfg, samples)
    n_comp = min(pca_cfg["n_components"], log2cpm.shape[1])

    pca = PCA(n_components=n_comp)
    scores = pca.fit_transform(log2cpm.T)
    var_exp = pca.explained_variance_ratio_

    pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(n_comp)],
        "variance_explained": var_exp,
        "cumulative": np.cumsum(var_exp),
    }).to_csv(os.path.join(out_dir, "pca_variance.tsv"), sep="\t", index=False)

    pd.DataFrame(
        scores, index=samples["all"],
        columns=[f"PC{i+1}" for i in range(n_comp)],
    ).to_csv(os.path.join(out_dir, "pca_scores.tsv"), sep="\t")

    if n_comp >= 2:
        log.info(f"  PC1={var_exp[0]:.1%}, PC2={var_exp[1]:.1%}")
    else:
        log.info(f"  PC1={var_exp[0]:.1%}")

    plot_pairs = _valid_plot_pairs(n_comp, pca_cfg.get("plot_pairs", [[1, 2], [1, 3]]))
    if not plot_pairs:
        log.warning("  No valid PCA plot pairs for the available component count; skipping PCA scatter plots")

    plot_files = []
    for pc_x, pc_y in plot_pairs:
        fig, ax = plt.subplots(figsize=(8, 6))

        for idx, sample in enumerate(samples["all"]):
            cond = get_sample_condition(cfg, sample)
            ax.scatter(
                scores[idx, pc_x - 1], scores[idx, pc_y - 1],
                c=colors.get(cond, "#999"),
                marker="o",
                s=100, edgecolors="black", linewidths=0.5, zorder=3,
            )
            ax.annotate(
                get_sample_display_name(cfg, sample),
                (scores[idx, pc_x - 1], scores[idx, pc_y - 1]),
                fontsize=9, ha="left", va="bottom",
                xytext=(4, 4), textcoords="offset points",
            )

        # Legend
        for cond in conditions:
            ax.scatter([], [], c=colors.get(cond, "#999"), marker="o",
                       s=60, edgecolors="black", linewidths=0.5, label=cond)
        ax.legend(fontsize=8, loc="best")

        ax.set_xlabel(f"PC{pc_x} ({var_exp[pc_x-1]*100:.1f}%)")
        ax.set_ylabel(f"PC{pc_y} ({var_exp[pc_y-1]*100:.1f}%)")
        title = f"PCA \u2014 {' vs '.join(conditions)} ({strain})"
        if manifest is not None:
            title = f"PCA \u2014 {' vs '.join(conditions)} ({strain}, exploratory repeats)"
        ax.set_title(title)

        # Reference lines at origin
        ax.axhline(0, ls="--", c="grey", lw=1.0, alpha=0.5)
        ax.axvline(0, ls="--", c="grey", lw=1.0, alpha=0.5)

        fig.tight_layout()
        fname = f"pca_pc{pc_x}_pc{pc_y}.png"
        plot_name = plot_filename(cfg, fname)
        save_figure(fig, os.path.join(out_dir, fname), cfg)
        plot_files.append((pc_x, pc_y, plot_name))
        log.info(f"  Saved {plot_name}")

    file_descs = {
        "pca_variance.tsv": "Variance explained and cumulative variance for each principal component.",
        "pca_scores.tsv": "Sample scores (coordinates) in PC space.",
    }
    for pc_x, pc_y, plot_name in plot_files:
        file_descs[plot_name] = f"Scatter plot of samples on PC{pc_x} vs PC{pc_y}, colored by condition."
    write_readme(out_dir, "2c", "PCA", file_descs)

    note = ""
    if manifest is not None:
        note = "\n- Exploration mode: repeated physical samples retained from `sample_manifest.tsv`\n"

    pc_lines = []
    for i in range(min(n_comp, 3)):
        pc_lines.append(f"- PC{i+1}: {var_exp[i]:.1%} variance explained")
    if n_comp >= 3:
        pc_lines.append(f"- Cumulative (PC1-3): {sum(var_exp[:3]):.1%}")
    else:
        pc_lines.append(f"- Cumulative (PC1-{n_comp}): {sum(var_exp[:n_comp]):.1%}")

    plot_sections = "\n\n".join(
        img(cfg, "02_sample_analysis/pca", f"pca_pc{pc_x}_pc{pc_y}.png", f"PCA — PC{pc_x} vs PC{pc_y}")
        for pc_x, pc_y, _ in plot_files
    )

    doc = f"""\
### Parameters

- **Components**: {n_comp}
- **Input**: log2(CPM + 1) filtered expression matrix
{note}

### Results

{chr(10).join(pc_lines)}
"""
    if plot_sections:
        doc = f"{doc.rstrip()}\n\n{plot_sections}\n"
    update_section(cfg, "2c", "PCA", doc)


def _load_log2expr(cfg, samples):
    return get_filtered_expr(cfg, log2=True)[samples["all"]]


def _valid_plot_pairs(n_comp, plot_pairs):
    """Return unique PCA axis pairs that fit within the available components."""
    valid = []
    seen = set()
    for pc_x, pc_y in plot_pairs:
        pair = (int(pc_x), int(pc_y))
        if pair[0] < 1 or pair[1] < 1 or pair[0] == pair[1]:
            continue
        if pair[0] > n_comp or pair[1] > n_comp:
            continue
        if pair in seen:
            continue
        seen.add(pair)
        valid.append(pair)
    return valid
