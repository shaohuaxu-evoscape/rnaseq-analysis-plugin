"""Step 5a: Gene Clustering — K-means on top variable genes."""

import os
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt

from ..core.config_data import get_filtered_expr, get_sample_manifest, get_sample_names
from ..core.config_runtime import ensure_output_dir, get_conditions, get_timepoints, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import (apply_style, style_axes, save_figure,
                        get_recipe_colors, plot_filename, COL1_5)

log = setup_logger("cluster_analysis.gene_clustering")


def run(cfg):
    clust_cfg = cfg["sample_analysis"]["clustering"]
    if not clust_cfg["enabled"]:
        log.info("Clustering disabled, skipping")
        return

    log.info("── Step 5a: Gene Clustering ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    manifest = get_sample_manifest(cfg)
    out_dir = ensure_output_dir(cfg, "05_cluster_analysis/gene_clustering")
    conditions = get_conditions(cfg)
    timepoints = get_timepoints(cfg)
    colors = get_recipe_colors(cfg)

    log2cpm = _load_log2expr(cfg, samples)
    top_n = clust_cfg["top_n_genes"]
    gene_var = log2cpm.var(axis=1).sort_values(ascending=False)
    data = log2cpm.loc[gene_var.head(top_n).index].copy()
    log.info(f"  Top {top_n} variable genes selected")

    # Z-score normalize per gene
    z_data = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)
    z_data = z_data.dropna()  # drop zero-variance genes

    # K selection
    k_min, k_max = clust_cfg["k_range"]
    k_range = range(k_min, k_max + 1)
    rs = clust_cfg["random_state"]
    inertias, silhouettes = [], []
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=rs, n_init=20)
        labels = km.fit_predict(z_data.values)
        inertias.append(km.inertia_)
        silhouettes.append(silhouette_score(z_data.values, labels))

    # Final K
    fixed_k = clust_cfg.get("fixed_k")
    if fixed_k is None:
        fixed_k = list(k_range)[np.argmax(silhouettes)]
    log.info(f"  K={fixed_k}")

    # K selection plot (reference style: 7x4, dual-axis, selected K line, combined legend)
    fig, ax1 = plt.subplots(figsize=(7, 4))
    ks = list(k_range)

    # Inertia (left axis)
    ln1 = ax1.plot(ks, inertias, "s-", color="#1f77b4", ms=5, lw=1.5, label="Inertia (elbow)")
    ax1.set_ylabel("Inertia (within-cluster SS)", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")

    # Silhouette (right axis)
    ax2 = ax1.twinx()
    ln2 = ax2.plot(ks, silhouettes, "o-", color="#d62728", ms=5, lw=1.5, label="Silhouette score")
    ax2.set_ylabel("Silhouette score", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")

    # Selected K vertical line
    ax1.axvline(fixed_k, ls="--", color="#333", lw=1.0, alpha=0.7)

    # Combined legend
    lns = ln1 + ln2 + [plt.Line2D([0], [0], ls="--", color="#333", lw=1.0)]
    labs = ["Inertia (elbow)", "Silhouette score", f"Selected K={fixed_k}"]
    ax1.legend(lns, labs, fontsize=8, loc="center right")

    ax1.set_xlabel("Number of clusters (K)")
    ax1.set_xticks(ks)
    ax1.set_title("K Selection \u2014 Elbow + Silhouette")
    # Keep all four spines for this plot
    for spine in ["top", "right", "left", "bottom"]:
        ax1.spines[spine].set_visible(True)
        ax1.spines[spine].set_linewidth(0.5)
        ax2.spines[spine].set_visible(True)
        ax2.spines[spine].set_linewidth(0.5)
    fig.tight_layout()
    k_selection_plot = plot_filename(cfg, "k_selection.png")
    save_figure(fig, os.path.join(out_dir, "k_selection.png"), cfg)

    # Final clustering

    km = KMeans(n_clusters=fixed_k, random_state=rs, n_init=30)
    labels = km.fit_predict(z_data.values)
    z_data_c = z_data.copy()
    z_data_c["cluster"] = labels + 1  # 1-indexed

    # Save assignments
    pd.DataFrame({"gene": z_data_c.index, "cluster": z_data_c["cluster"]}).to_csv(
        os.path.join(out_dir, "cluster_assignments.tsv"), sep="\t", index=False)

    # ── Cluster profiles ──
    cond1, cond2 = get_conditions(cfg)[:2]
    cond1_samples = samples[cond1]
    cond2_samples = samples[cond2]
    n_cols = 4
    n_rows = int(np.ceil(fixed_k / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 3.8 * n_rows))
    if fixed_k == 1:
        axes = np.array([[axes]])
    axes = np.atleast_2d(axes)

    for ci in range(fixed_k):
        r, col = divmod(ci, n_cols)
        ax = axes[r, col]
        cluster_id = ci + 1
        mask = z_data_c["cluster"] == cluster_id
        genes = z_data_c.loc[mask]
        n_genes = mask.sum()

        cond1_vals = _collapse_condition_timepoints(
            genes, manifest, cond1, cond1_samples, timepoints
        )
        cond2_vals = _collapse_condition_timepoints(
            genes, manifest, cond2, cond2_samples, timepoints
        )

        # Divergence score and trend arrows
        div = cond2_vals.mean() - cond1_vals.mean()
        cond1_trend = "\u2191" if cond1_vals[:, -1].mean() > cond1_vals[:, 0].mean() else "\u2193"
        cond2_trend = "\u2191" if cond2_vals[:, -1].mean() > cond2_vals[:, 0].mean() else "\u2193"

        # Individual gene traces (semi-transparent)
        for j in range(n_genes):
            ax.plot(timepoints, cond1_vals[j], color=colors[cond1], alpha=0.03, lw=1.0)
            ax.plot(timepoints, cond2_vals[j], color=colors[cond2], alpha=0.03, lw=1.0)

        # Mean ± SD shading
        for vals, color, label in [
            (cond1_vals, colors[cond1], cond1),
            (cond2_vals, colors[cond2], cond2),
        ]:
            mean = vals.mean(axis=0)
            std = vals.std(axis=0)
            ax.fill_between(timepoints, mean - std, mean + std, alpha=0.12, color=color)
            ax.plot(timepoints, mean, "o-", color=color, ms=5, lw=2, label=label, zorder=5)

        # Title: cluster name bold, trend line with colored condition labels
        ax.set_title(f"Cluster {cluster_id} (n={n_genes})", fontsize=10, fontweight="bold")
        # Colored trend annotation below title
        ax.text(0.5, 1.0, f"{cond1}{cond1_trend}", color=colors[cond1], fontsize=8, fontweight="bold",
                ha="right", va="top", transform=ax.transAxes)
        ax.text(0.5, 1.0, f"  {cond2}{cond2_trend}  div={div:+.2f}", color=colors[cond2], fontsize=8,
                fontweight="bold", ha="left", va="top", transform=ax.transAxes)

        ax.set_xlabel("Time (h)", fontsize=8)
        ax.axhline(0, ls="--", c="grey", lw=1.0)
        ax.set_xticks(timepoints)
        style_axes(ax)

        if ci == 0:
            ax.legend(fontsize=8)

    for ci in range(fixed_k, n_rows * n_cols):
        r, col = divmod(ci, n_cols)
        axes[r, col].set_visible(False)

    fig.supylabel("Z-score", fontsize=11)
    fig.suptitle(f"Time-series Gene Clusters (top {top_n} variable genes, K={fixed_k})",
                 fontsize=13, y=1.01)
    fig.tight_layout()
    cluster_profiles_plot = plot_filename(cfg, "cluster_profiles.png")
    save_figure(fig, os.path.join(out_dir, "cluster_profiles.png"), cfg)
    log.info(f"  Saved {cluster_profiles_plot}")

    # ── Heatmap (with cluster color sidebar, hierarchical ordering within clusters) ──
    all_samples = cond1_samples + cond2_samples
    data_only = z_data_c[all_samples]
    cluster_col = z_data_c["cluster"]

    ordered_idx = []
    for c in sorted(cluster_col.unique()):
        genes_in_c = data_only.loc[cluster_col == c]
        if len(genes_in_c) > 1:
            dist = pdist(genes_in_c.values, metric="euclidean")
            link = linkage(dist, method="ward")
            order = leaves_list(link)
            ordered_idx.extend(genes_in_c.index[order].tolist())
        else:
            ordered_idx.extend(genes_in_c.index.tolist())

    ordered_data = data_only.loc[ordered_idx]
    ordered_clusters = cluster_col.loc[ordered_idx]

    vmax = 2.5
    fig, (ax_cbar, ax_heat) = plt.subplots(
        1, 2, figsize=(COL1_5, max(4, len(ordered_data) * 0.012 + 2)),
        gridspec_kw={"width_ratios": [0.3, 10]},
    )

    # Cluster color sidebar
    n_clusters = ordered_clusters.max()
    cmap_cluster = plt.colormaps.get_cmap("tab10").resampled(n_clusters)
    cluster_colors = cmap_cluster(ordered_clusters.values - 1)
    ax_cbar.imshow(cluster_colors.reshape(-1, 1, 4), aspect="auto", interpolation="nearest")
    ax_cbar.set_xticks([])
    ax_cbar.set_yticks([])
    ax_cbar.set_ylabel(f"Genes (n={len(ordered_data)})")

    # Heatmap
    im = ax_heat.imshow(ordered_data.values, aspect="auto", cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, interpolation="nearest")
    col_labels = [s.replace("-", "\n") for s in ordered_data.columns]
    ax_heat.set_xticks(range(len(col_labels)))
    ax_heat.set_xticklabels(col_labels, fontsize=9)
    ax_heat.set_yticks([])
    ax_heat.set_title("Z-scored Expression by Cluster")

    cb = fig.colorbar(im, ax=ax_heat, shrink=0.3, aspect=25, label="Z-score", pad=0.02)
    cb.ax.tick_params(labelsize=5)

    # Cluster legend below z-score colorbar
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=cmap_cluster(c - 1), label=f"C{c}")
                      for c in sorted(ordered_clusters.unique())]
    cb.ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.8),
                 fontsize=8, frameon=False, handlelength=1.0, handleheight=1.0,
                 ncol=1, title="Cluster", title_fontsize=9)
    cluster_heatmap_plot = plot_filename(cfg, "cluster_heatmap.png")
    save_figure(fig, os.path.join(out_dir, "cluster_heatmap.png"), cfg)
    log.info(f"  Saved {cluster_heatmap_plot}")

    write_readme(out_dir, "5a", "Gene Clustering", {
        k_selection_plot: "Elbow + silhouette plot for K selection with selected K highlighted",
        "cluster_assignments.tsv": "Per-gene cluster assignments (gene ID and 1-indexed cluster number)",
        cluster_profiles_plot: "Per-cluster Z-score expression profiles with condition-level mean +/- SD and individual gene traces",
        cluster_heatmap_plot: "Z-scored expression heatmap with cluster color sidebar and hierarchical ordering within clusters",
    })

    best_sil = max(silhouettes)
    doc = (
        "### Parameters\n\n"
        f"- Top **{top_n}** variable genes selected from filtered log2CPM\n"
        f"- K range tested: {k_min}--{k_max}; final K = **{fixed_k}**\n"
        f"- Best silhouette score: {best_sil:.3f}\n"
        f"- Z-score normalization per gene; K-means with n_init=30, random_state={rs}\n\n"
        "### Results\n\n"
        f"- **{len(z_data)}** genes clustered into **{fixed_k}** clusters "
        f"(after dropping zero-variance genes)\n"
        f"- Cluster sizes: {', '.join(f'C{c+1}={int((labels==c).sum())}' for c in range(fixed_k))}\n\n"
        f"{img(cfg, '05_cluster_analysis/gene_clustering', 'cluster_profiles.png', 'Cluster expression profiles with condition-level mean and gene traces')}\n\n"
        f"{img(cfg, '05_cluster_analysis/gene_clustering', 'cluster_heatmap.png', 'Z-scored expression heatmap by cluster')}\n"
    )
    update_section(cfg, "5a", "Gene Clustering", doc)

    log.info("  Gene clustering complete")


def _load_log2expr(cfg, samples):
    return get_filtered_expr(cfg, log2=True)[samples["all"]]


def _collapse_condition_timepoints(data, manifest, condition, sample_ids, timepoints):
    """Collapse repeated samples to one mean profile per timepoint."""
    if manifest is None:
        return data[sample_ids].values

    rows = manifest.loc[manifest["sample_id"].isin(sample_ids)].copy()
    values = []
    for tp in timepoints:
        tp_samples = rows.loc[rows["timepoint"].astype(int) == int(tp), "sample_id"].tolist()
        if not tp_samples:
            raise ValueError(
                f"No samples found for {condition} at {tp}h in cluster plotting"
            )
        values.append(data[tp_samples].mean(axis=1).to_numpy())
    return np.column_stack(values)
