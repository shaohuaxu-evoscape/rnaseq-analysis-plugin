"""Step 2d: Sample Dendrogram — hierarchical clustering of 10 samples.

Distance metric: 1 - Pearson r, Ward linkage.
"""

import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
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

log = setup_logger("sample_analysis.dendrogram")


def run(cfg):
    log.info("── Step 2d: Sample Dendrogram ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    out_dir = ensure_output_dir(cfg, "02_sample_analysis/dendrogram")
    colors = get_recipe_colors(cfg)
    conditions = get_conditions(cfg)
    strain = cfg["project"].get("strain", "")
    manifest = get_sample_manifest(cfg)

    log2cpm = _load_log2expr(cfg, samples)

    # Pearson correlation → distance
    corr = log2cpm.corr(method="pearson")
    dist = 1 - corr.values
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2
    condensed = squareform(dist)
    link = linkage(condensed, method="ward")

    sample_names = corr.columns.tolist()

    # Map sample → recipe color
    label_colors = {}
    for s in sample_names:
        cond = get_sample_condition(cfg, s)
        label_colors[s] = colors.get(cond, "black")

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    display_labels = [get_sample_display_name(cfg, s) for s in sample_names]
    label_map = dict(zip(display_labels, sample_names))
    dendrogram(link, labels=display_labels, ax=ax,
               leaf_font_size=9, leaf_rotation=45)

    # Color tick labels by recipe
    for lbl in ax.get_xticklabels():
        sample_name = label_map.get(lbl.get_text(), lbl.get_text())
        lbl.set_color(label_colors.get(sample_name, "black"))

    ax.set_ylabel("Distance (1 \u2212 Pearson r, Ward linkage)")
    title = f"Sample Dendrogram \u2014 {' vs '.join(conditions)} ({strain})"
    if manifest is not None:
        title = f"Sample Dendrogram \u2014 {' vs '.join(conditions)} ({strain}, exploratory repeats)"
    ax.set_title(title)
    fig.tight_layout()

    dendrogram_plot = plot_filename(cfg, "dendrogram.png")
    save_figure(fig, os.path.join(out_dir, "dendrogram.png"), cfg)
    log.info(f"  Saved {dendrogram_plot}")

    # Save linkage info
    linkage_df = pd.DataFrame(link, columns=["cluster1", "cluster2", "distance", "n_samples"])
    linkage_df.to_csv(os.path.join(out_dir, "linkage.tsv"), sep="\t", index=False)

    write_readme(out_dir, "2d", "Sample Dendrogram", {
        dendrogram_plot: "Hierarchical clustering dendrogram (1 - Pearson r distance, Ward linkage).",
        "linkage.tsv": "Linkage matrix from hierarchical clustering (cluster pairs, distances, sizes).",
    })

    note = ""
    if manifest is not None:
        note = "\n- Exploration mode: repeated physical samples retained from `sample_manifest.tsv`\n"

    doc = f"""\
### Parameters

- **Distance metric**: 1 - Pearson r
- **Linkage method**: Ward
- **Samples**: {len(sample_names)}
{note}

### Results

- Strain: {strain}

{img(cfg, "02_sample_analysis/dendrogram", "dendrogram.png", "Hierarchical clustering dendrogram (Ward linkage)")}
"""
    update_section(cfg, "2d", "Sample Dendrogram", doc)


def _load_log2expr(cfg, samples):
    return get_filtered_expr(cfg, log2=True)[samples["all"]]
