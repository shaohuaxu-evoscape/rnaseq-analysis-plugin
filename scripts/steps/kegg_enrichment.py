"""Step 3c: KEGG Pathway Enrichment — Fisher's exact test on DE genes.

Uses gene→UniProt→KEGG mapping chain and bulk pathway-gene query.
Gene names are matched by the prefix configured in project.gene_id_prefix.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from ..core.config_runtime import ensure_output_dir, get_conditions, write_readme
from ..core.de_helpers import load_de_data, direction_column_name
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.net import fetch_cached as _fetch_cached
from ..core.plotting import apply_style, style_axes, save_figure, plot_filename
from ..core.stats import bh_fdr

log = setup_logger("differential_analysis.kegg_enrichment")


def _parse_kegg_tsv(text):
    rows = []
    for line in text.strip().split("\n"):
        if "\t" in line:
            rows.append(line.split("\t", 1))
    return rows


def _build_id_mapping(cache_dir, organism, taxid, gene_id_prefix):
    """Build gene → KEGG gene mapping via UniProt cross-reference.

    Chain: gene name → UniProt accession → KEGG gene ID
    """
    # KEGG ↔ UniProt
    text = _fetch_cached(
        f"https://rest.kegg.jp/conv/uniprot/{organism}",
        cache_dir, f"kegg_uniprot_conv.tsv"
    )
    uniprot2kegg = {}
    for row in _parse_kegg_tsv(text):
        kegg_id = row[0].strip()
        uniprot_acc = row[1].strip().replace("up:", "")
        uniprot2kegg[uniprot_acc] = kegg_id

    # UniProt gene names
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           f"?query=organism_id:{taxid}&fields=accession,gene_names&format=tsv")
    text = _fetch_cached(url, cache_dir, "uniprot_genes.tsv")
    gene_to_uniprot = {}
    for line in text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].strip()
        for gn in parts[1].strip().split():
            if gene_id_prefix and gn.startswith(gene_id_prefix):
                gene_to_uniprot[gn] = acc
                break

    # Chain: gene → UniProt → KEGG
    gene_to_kegg = {}
    for gene, uniprot_acc in gene_to_uniprot.items():
        if uniprot_acc in uniprot2kegg:
            gene_to_kegg[gene] = uniprot2kegg[uniprot_acc]

    log.info(f"  KEGG genes: {len(uniprot2kegg)}")
    log.info(f"  UniProt entries with prefix '{gene_id_prefix}': {len(gene_to_uniprot)}")
    log.info(f"  Gene → KEGG mapped: {len(gene_to_kegg)}")
    return gene_to_kegg


def _build_pathway_sets(cache_dir, organism, gene_to_kegg):
    """Build pathway gene sets using bulk KEGG query."""
    # Pathway names
    text = _fetch_cached(
        f"https://rest.kegg.jp/list/pathway/{organism}",
        cache_dir, "kegg_pathway_list.tsv"
    )
    pathway_names = {}
    for row in _parse_kegg_tsv(text):
        pid = row[0].replace("path:", "")
        name = row[1].split(" - ")[0].strip()
        pathway_names[pid] = name

    # Bulk pathway → gene link (single API call!)
    text = _fetch_cached(
        f"https://rest.kegg.jp/link/{organism}/pathway",
        cache_dir, "kegg_pathway_genes.tsv"
    )
    kegg_to_gene = {v: k for k, v in gene_to_kegg.items()}
    pathway_genes = {}
    for row in _parse_kegg_tsv(text):
        pid = row[0].replace("path:", "")
        kegg_gene = row[1].strip()
        gene = kegg_to_gene.get(kegg_gene)
        if gene:
            pathway_genes.setdefault(pid, set()).add(gene)

    log.info(f"  Pathways loaded: {len(pathway_names)}")
    log.info(f"  Pathways with mapped genes: {len(pathway_genes)}")
    return pathway_names, pathway_genes


def run(cfg):
    """Run KEGG pathway enrichment analysis."""
    kegg_cfg = cfg["differential_analysis"]["kegg_enrichment"]
    if not kegg_cfg["enabled"]:
        log.info("KEGG enrichment disabled, skipping")
        return

    log.info("── Step 3c: KEGG Enrichment ──")
    apply_style(cfg)

    out_dir = ensure_output_dir(cfg, "03_differential_analysis/kegg_enrichment")
    cache_dir = ensure_output_dir(cfg, "03_differential_analysis/kegg_enrichment/cache")
    organism = kegg_cfg["organism"]

    # Load DE data
    de_genes, background_genes, de_direction = load_de_data(cfg)
    log.info(f"  DE genes: {len(de_genes)}, Background: {len(background_genes)}")
    cond2 = get_conditions(cfg)[1]
    up_col = direction_column_name("up", cond2)
    down_col = direction_column_name("down", cond2)

    # Build mappings
    taxid = cfg.get("differential_analysis", {}).get("go_enrichment", {}).get("organism_taxid") \
            or cfg["project"].get("organism_taxid", 0)
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")
    gene_to_kegg = _build_id_mapping(cache_dir, organism, taxid, gene_id_prefix)
    pathway_names, pathway_genes = _build_pathway_sets(cache_dir, organism, gene_to_kegg)

    # Fisher's exact test
    min_genes = kegg_cfg["min_pathway_genes"]
    max_genes = kegg_cfg["max_pathway_genes"]
    sig_thresh = kegg_cfg["sig_threshold"]

    de_set = set(de_genes)
    bg_set = set(background_genes)
    N = len(bg_set)
    K = len(de_set & bg_set)

    results = []
    for pid, pgenes in pathway_genes.items():
        pgenes_in_bg = pgenes & bg_set
        n = len(pgenes_in_bg)
        if n < min_genes or n > max_genes:
            continue

        overlap = de_set & pgenes_in_bg
        k = len(overlap)
        if k == 0:
            continue

        n_up = sum(1 for g in overlap if de_direction.get(g) == "up")
        n_down = k - n_up

        cell_d = N - K - n + k
        if cell_d < 0:
            continue
        table = [[k, K - k], [n - k, cell_d]]
        _, pval = fisher_exact(table, alternative="greater")

        results.append({
            "pathway_id": pid,
            "pathway_name": pathway_names.get(pid, pid),
            "DE_in_pathway": k, up_col: n_up, down_col: n_down,
            "pathway_size": n, "gene_ratio": k / n,
            "bg_size": N, "DE_total": K, "pvalue": pval,
            "overlap_genes": ";".join(sorted(overlap)),
        })

    if not results:
        log.warning("  No KEGG enrichment results")
        return

    df = pd.DataFrame(results).sort_values("pvalue").reset_index(drop=True)
    df["fdr"] = bh_fdr(df["pvalue"].values)
    df = df.sort_values("fdr")
    df.to_csv(os.path.join(out_dir, "kegg_enrichment.tsv"), sep="\t", index=False)

    n_sig_p = (df["pvalue"] < sig_thresh).sum()
    n_sig_fdr = (df["fdr"] < sig_thresh).sum()
    log.info(f"  Pathways tested: {len(df)}, sig p<{sig_thresh}: {n_sig_p}, sig FDR<{sig_thresh}: {n_sig_fdr}")

    # Plot
    top_n = kegg_cfg["top_n_plot"]
    plot_df = df[df["pvalue"] < sig_thresh].head(top_n)
    if not plot_df.empty:
        _plot_enrichment(plot_df, out_dir, cfg, sig_thresh)

    write_readme(out_dir, "3c", "KEGG Enrichment", {
        "kegg_enrichment.tsv": "KEGG pathway enrichment results with Fisher's exact test p-values, FDR, and overlap genes",
        plot_filename(cfg, "kegg_enrichment.png"): "Dot plot of top enriched KEGG pathways by gene ratio and significance",
        "cache/kegg_uniprot_conv.tsv": "Cached KEGG-to-UniProt ID conversion table",
        "cache/uniprot_genes.tsv": "Cached UniProt accession-to-gene-name mapping",
        "cache/kegg_pathway_list.tsv": "Cached KEGG pathway names for the organism",
        "cache/kegg_pathway_genes.tsv": "Cached KEGG pathway-to-gene membership links",
    })

    subdir = "03_differential_analysis/kegg_enrichment"
    doc = (
        f"### Parameters\n"
        f"- Organism: {organism}\n"
        f"- Min pathway genes: {min_genes}\n"
        f"- Max pathway genes: {max_genes}\n"
        f"- Significance threshold: {sig_thresh}\n"
        f"- Top N plot: {top_n}\n\n"
        f"### Results\n"
        f"- DE genes: {len(de_genes)}, Background: {len(background_genes)}\n"
        f"- Pathways tested: {len(df)}\n"
        f"- Significant (p < {sig_thresh}): {n_sig_p}\n"
        f"- Significant (FDR < {sig_thresh}): {n_sig_fdr}\n\n"
        f"{img(cfg, subdir, 'kegg_enrichment.png', 'KEGG pathway enrichment dot plot')}\n"
    )
    update_section(cfg, "3c", "KEGG Enrichment", doc)

    log.info("  KEGG enrichment complete")


def _plot_enrichment(df, out_dir, cfg, sig_thresh):
    """Combined bar + dot plot for enriched KEGG pathways, matching reference style."""
    from ..core.plotting import COL2
    df = df.sort_values("pvalue", ascending=False).reset_index(drop=True)
    n = len(df)

    neg_log_p = -np.log10(df["pvalue"].values.clip(min=1e-50))
    norm = mcolors.Normalize(vmin=max(0, neg_log_p.min() - 0.2),
                             vmax=neg_log_p.max() + 0.2)
    cmap = plt.cm.RdPu

    sizes = df["DE_in_pathway"].values.astype(float)
    s_min, s_max = 60, 600
    def scale_size(c):
        if sizes.max() == sizes.min():
            return (s_min + s_max) / 2
        return s_min + (c - sizes.min()) / (sizes.max() - sizes.min()) * (s_max - s_min)

    fig_h = max(4, n * 0.5 + 2)
    fig = plt.figure(figsize=(COL2, fig_h))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.06], wspace=0.12,
                          height_ratios=[1, 1], hspace=0.05)
    ax_bar = fig.add_subplot(gs[:, 0])
    ax_dot = fig.add_subplot(gs[:, 1])

    y_pos = np.arange(n)

    # Left: bar chart
    colors = cmap(norm(neg_log_p))
    ax_bar.barh(y_pos, df["DE_in_pathway"].values, color=colors,
                edgecolor="black", linewidth=0.5, height=0.6)
    for i, (_, row) in enumerate(df.iterrows()):
        marker = " **" if row["fdr"] < sig_thresh else " *"
        ax_bar.text(row["DE_in_pathway"] + 0.2, i,
                    f"{int(row['DE_in_pathway'])}{marker}",
                    va="center", fontsize=8)
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(df["pathway_name"].values, fontsize=9)
    ax_bar.set_xlabel("Number of DE Genes", fontsize=10)
    ax_bar.set_title("Gene Count", fontsize=11)
    ax_bar.set_ylim(-0.6, n - 0.4)
    ax_bar.grid(axis="x", alpha=0.3, ls="--")
    style_axes(ax_bar)

    # Right: dot plot
    scaled = np.array([scale_size(c) for c in sizes])
    sc = ax_dot.scatter(df["gene_ratio"].values, y_pos, s=scaled,
                        c=neg_log_p, cmap=cmap, norm=norm,
                        edgecolors="black", linewidths=0.5, zorder=3)
    xmin, xmax = df["gene_ratio"].min(), df["gene_ratio"].max()
    x_pad = max((xmax - xmin) * 0.25, 0.05)
    ax_dot.set_xlim(xmin - x_pad, xmax + x_pad)
    ax_dot.set_ylim(-0.8, n - 0.2)
    ax_dot.set_yticks(y_pos)
    ax_dot.set_yticklabels([])
    ax_dot.set_xlabel("Gene Ratio (DE / Pathway)", fontsize=10)
    ax_dot.set_title("Gene Ratio", fontsize=11)
    ax_dot.grid(axis="x", alpha=0.3, ls="--")
    style_axes(ax_dot)

    # Colorbar
    ax_cbar = fig.add_subplot(gs[0, 2])
    cbar = fig.colorbar(sc, cax=ax_cbar)
    cbar.set_label("$-$log$_{10}$(p)", fontsize=9)

    # Size legend
    ax_leg = fig.add_subplot(gs[1, 2])
    ax_leg.set_axis_off()
    legend_vals = sorted(set([int(sizes.min()), int(sizes.max())]))
    legend_vals = [max(1, v) for v in legend_vals]
    handles = []
    for v in legend_vals:
        handles.append(ax_leg.scatter([], [], s=scale_size(v), c="grey",
                                      edgecolors="black", linewidths=0.5, label=str(v)))
    ax_leg.legend(handles=handles, title="Count", loc="upper center",
                  fontsize=9, title_fontsize=10, framealpha=0.9,
                  borderpad=1.2, labelspacing=2.5, handletextpad=1.5)

    conds = get_conditions(cfg)
    fig.suptitle(f"KEGG Pathway Enrichment — {conds[0]} vs {conds[1]}\n"
                 f"Fisher's exact test (* p<{sig_thresh}, ** FDR<{sig_thresh})",
                 fontsize=12, y=1.02)

    kegg_enrichment_plot = plot_filename(cfg, "kegg_enrichment.png")
    save_figure(fig, os.path.join(out_dir, "kegg_enrichment.png"), cfg)
    log.info(f"  Saved {kegg_enrichment_plot}")
