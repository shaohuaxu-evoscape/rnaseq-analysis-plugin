"""Step 3b: GO Enrichment — Fisher's exact test on DE genes.

Two-step UniProt mapping: accession→gene name + accession→GO terms.
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

log = setup_logger("differential_analysis.go_enrichment")

ONTOLOGY_NAMES = {"BP": "Biological Process", "MF": "Molecular Function", "CC": "Cellular Component"}


def _build_go_mapping(cache_dir, taxid, gene_id_prefix):
    """Two-step mapping: UniProt accession→gene name, then accession→GO terms."""
    # Step 1: accession → gene name (filtered by gene_id_prefix)
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           f"?query=organism_id:{taxid}&fields=accession,gene_names&format=tsv")
    text = _fetch_cached(url, cache_dir, f"uniprot_genes.tsv")
    acc_to_gene = {}
    for line in text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].strip()
        for gn in parts[1].strip().split():
            if gene_id_prefix and gn.startswith(gene_id_prefix):
                acc_to_gene[acc] = gn
                break

    # Step 2: accession → GO terms (BP, MF, CC)
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           f"?query=organism_id:{taxid}+AND+go:*&fields=accession,go_id,go_p,go_f,go_c&format=tsv")
    text = _fetch_cached(url, cache_dir, f"uniprot_go.tsv")

    gene_go = {}   # gene name → set of GO IDs
    go_terms = {}  # GO ID → (ontology, term_name)

    for line in text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        acc = parts[0].strip()
        gene = acc_to_gene.get(acc)
        if not gene:
            continue
        if gene not in gene_go:
            gene_go[gene] = set()

        for col_idx, ontology in [(2, "BP"), (3, "MF"), (4, "CC")]:
            if col_idx >= len(parts) or not parts[col_idx].strip():
                continue
            for entry in parts[col_idx].strip().split("; "):
                if "[GO:" not in entry:
                    continue
                bracket_start = entry.rfind("[GO:")
                go_id = entry[bracket_start + 1:entry.rfind("]")]
                term_name = entry[:bracket_start].strip()
                gene_go[gene].add(go_id)
                go_terms[go_id] = (ontology, term_name)

    log.info(f"  Genes with GO annotation: {len(gene_go)}")
    log.info(f"  Unique GO terms: {len(go_terms)}")

    # Build term → gene set
    term_genes = {}
    for gene, terms in gene_go.items():
        for go_id in terms:
            term_genes.setdefault(go_id, set()).add(gene)

    return gene_go, go_terms, term_genes


def run(cfg):
    """Run GO enrichment analysis."""
    go_cfg = cfg["differential_analysis"]["go_enrichment"]
    if not go_cfg["enabled"]:
        log.info("GO enrichment disabled, skipping")
        return

    log.info("── Step 3b: GO Enrichment ──")
    apply_style(cfg)

    out_dir = ensure_output_dir(cfg, "03_differential_analysis/go_enrichment")
    cache_dir = ensure_output_dir(cfg, "03_differential_analysis/go_enrichment/cache")
    # Load DE data
    de_genes, background_genes, de_direction = load_de_data(cfg)
    log.info(f"  DE genes: {len(de_genes)}, Background: {len(background_genes)}")
    cond2 = get_conditions(cfg)[1]
    up_col = direction_column_name("up", cond2)
    down_col = direction_column_name("down", cond2)

    # Build GO mapping
    taxid = go_cfg.get("organism_taxid") or cfg["project"].get("organism_taxid", 0)
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")
    gene_go, go_terms, term_genes = _build_go_mapping(cache_dir, taxid, gene_id_prefix)

    # Run enrichment
    min_genes = go_cfg["min_term_genes"]
    max_genes = go_cfg["max_term_genes"]
    sig_thresh = go_cfg["sig_threshold"]

    de_set = set(de_genes)
    bg_set = set(background_genes)
    N = len(bg_set)
    K = len(de_set & bg_set)

    results = []
    for go_id, tgenes in term_genes.items():
        tgenes_in_bg = tgenes & bg_set
        n = len(tgenes_in_bg)
        if n < min_genes or n > max_genes:
            continue
        overlap = de_set & tgenes_in_bg
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

        ontology, term_name = go_terms[go_id]
        results.append({
            "go_id": go_id, "term_name": term_name, "ontology": ontology,
            "DE_in_term": k, up_col: n_up, down_col: n_down,
            "term_size": n, "gene_ratio": k / n,
            "bg_size": N, "DE_total": K, "pvalue": pval,
            "overlap_genes": ";".join(sorted(overlap)),
        })

    if not results:
        log.warning("  No GO enrichment results")
        return

    df = pd.DataFrame(results)

    # BH FDR per ontology
    dfs = []
    for ont in df["ontology"].unique():
        sub = df[df["ontology"] == ont].sort_values("pvalue").reset_index(drop=True)
        sub["fdr"] = bh_fdr(sub["pvalue"].values)
        dfs.append(sub)
    df = pd.concat(dfs, ignore_index=True).sort_values("pvalue")

    df.to_csv(os.path.join(out_dir, "go_enrichment.tsv"), sep="\t", index=False)

    n_sig_p = (df["pvalue"] < sig_thresh).sum()
    n_sig_fdr = (df["fdr"] < sig_thresh).sum()
    log.info(f"  Results: {len(df)} terms, {n_sig_p} sig (p<{sig_thresh}), {n_sig_fdr} sig (FDR<{sig_thresh})")

    # Plot: combined bubble plot (BP + CC + MF in one figure)
    top_n = go_cfg["top_n_plot"]
    _plot_combined_bubble(df, go_cfg["ontologies"], top_n, sig_thresh, out_dir, cfg)

    write_readme(out_dir, "3b", "GO Enrichment", {
        "go_enrichment.tsv": "GO enrichment results with Fisher's exact test p-values, FDR, and overlap genes per term",
        plot_filename(cfg, "go_enrichment.png"): "Combined bubble plot of top enriched GO terms across BP, MF, and CC ontologies",
        "cache/uniprot_genes.tsv": "Cached UniProt accession-to-gene-name mapping",
        "cache/uniprot_go.tsv": "Cached UniProt accession-to-GO-term mapping",
    })

    subdir = "03_differential_analysis/go_enrichment"
    doc = (
        f"### Parameters\n"
        f"- Organism taxid: {taxid}\n"
        f"- Min term genes: {min_genes}\n"
        f"- Max term genes: {max_genes}\n"
        f"- Significance threshold: {sig_thresh}\n"
        f"- Top N plot: {top_n}\n\n"
        f"### Results\n"
        f"- DE genes: {len(de_genes)}, Background: {len(background_genes)}\n"
        f"- Terms tested: {len(df)}\n"
        f"- Significant (p < {sig_thresh}): {n_sig_p}\n"
        f"- Significant (FDR < {sig_thresh}): {n_sig_fdr}\n\n"
        f"{img(cfg, subdir, 'go_enrichment.png', 'GO enrichment bubble plot')}\n"
    )
    update_section(cfg, "3b", "GO Enrichment", doc)

    log.info("  GO enrichment complete")


def _plot_combined_bubble(df, ontologies, top_n, sig_thresh, out_dir, cfg):
    """Combined bar + dot plot: BP/MF/CC stacked vertically, matching reference style."""
    from ..core.plotting import COL2

    panels = []
    for ont in ontologies:
        ont_df = df[(df["ontology"] == ont) & (df["pvalue"] < sig_thresh)].head(top_n)
        if not ont_df.empty:
            panels.append((ont, ont_df.sort_values("pvalue", ascending=False).reset_index(drop=True)))

    if not panels:
        return

    # Global scales across all ontologies
    all_neg_log = np.concatenate([
        -np.log10(sig["pvalue"].values.clip(min=1e-50)) for _, sig in panels
    ])
    global_norm = mcolors.Normalize(vmin=all_neg_log.min() - 0.1, vmax=all_neg_log.max() + 0.1)
    cmap = plt.cm.RdPu

    all_counts = np.concatenate([sig["DE_in_term"].values for _, sig in panels])
    cnt_min, cnt_max = all_counts.min(), all_counts.max()
    s_min, s_max = 50, 500

    def scale_size(c):
        if cnt_max == cnt_min:
            return (s_min + s_max) / 2
        return s_min + (c - cnt_min) / (cnt_max - cnt_min) * (s_max - s_min)

    # Layout: N rows × 3 cols (bar, dot, colorbar/legend)
    heights = [len(sig) for _, sig in panels]
    total_terms = sum(heights)
    fig_h = max(6, total_terms * 0.38 + len(panels) * 0.3 + 1.5)
    fig = plt.figure(figsize=(COL2, fig_h))
    gs = fig.add_gridspec(len(panels), 3, width_ratios=[1, 1, 0.06],
                          height_ratios=heights, wspace=0.10, hspace=0.15)

    last_sc = None
    for row_idx, (ont, sig) in enumerate(panels):
        n = len(sig)
        ax_bar = fig.add_subplot(gs[row_idx, 0])
        ax_dot = fig.add_subplot(gs[row_idx, 1])
        y_pos = np.arange(n)
        neg_log = -np.log10(sig["pvalue"].values.clip(min=1e-50))

        # Left: bar chart
        colors = cmap(global_norm(neg_log))
        ax_bar.barh(y_pos, sig["DE_in_term"].values, color=colors,
                    edgecolor="black", linewidth=0.5, height=0.65)
        for i, (_, row) in enumerate(sig.iterrows()):
            marker = " **" if row["fdr"] < sig_thresh else " *"
            ax_bar.text(row["DE_in_term"] + 0.15, i,
                        f"{int(row['DE_in_term'])}{marker}",
                        va="center", fontsize=8)
        ax_bar.set_yticks(y_pos)
        ax_bar.set_yticklabels(sig["term_name"].values, fontsize=8)
        ax_bar.set_ylim(-0.6, n - 0.4)
        ax_bar.grid(axis="x", alpha=0.3, ls="--")
        ax_bar.set_ylabel(ONTOLOGY_NAMES[ont], fontsize=10, fontweight="bold")
        if row_idx == len(panels) - 1:
            ax_bar.set_xlabel("Number of DE Genes", fontsize=10)
        style_axes(ax_bar)

        # Right: dot plot
        sizes = sig["DE_in_term"].values.astype(float)
        scaled = np.array([scale_size(c) for c in sizes])
        sc = ax_dot.scatter(sig["gene_ratio"].values, y_pos, s=scaled, c=neg_log,
                            cmap=cmap, norm=global_norm, edgecolors="black",
                            linewidths=0.5, zorder=3)
        last_sc = sc
        xmin, xmax = sig["gene_ratio"].min(), sig["gene_ratio"].max()
        x_pad = (xmax - xmin) * 0.25 if xmax > xmin else 0.05
        ax_dot.set_xlim(xmin - x_pad, xmax + x_pad)
        ax_dot.set_ylim(-0.6, n - 0.4)
        ax_dot.set_yticks(y_pos)
        ax_dot.set_yticklabels([])
        ax_dot.grid(axis="x", alpha=0.3, ls="--")
        if row_idx == len(panels) - 1:
            ax_dot.set_xlabel("Gene Ratio (DE / Term)", fontsize=10)
        if row_idx == 0:
            ax_bar.set_title("Gene Count", fontsize=11)
            ax_dot.set_title("Gene Ratio", fontsize=11)
        style_axes(ax_dot)

    # Colorbar
    if last_sc is not None:
        ax_cbar = fig.add_subplot(gs[0, 2])
        cbar = fig.colorbar(last_sc, cax=ax_cbar)
        cbar.set_label("$-$log$_{10}$(p)", fontsize=9)

    # Size legend
    ax_leg = fig.add_subplot(gs[-1, 2])
    ax_leg.set_axis_off()
    legend_vals = sorted(set([int(cnt_min), int(cnt_max)]))
    legend_vals = [max(1, v) for v in legend_vals]
    handles = []
    for v in legend_vals:
        handles.append(ax_leg.scatter([], [], s=scale_size(v), c="grey",
                                      edgecolors="black", linewidths=0.5, label=str(v)))
    ax_leg.legend(handles=handles, title="Count", loc="upper center",
                  fontsize=9, title_fontsize=10, framealpha=0.9,
                  borderpad=1.0, labelspacing=2.0, handletextpad=1.2)

    conds = get_conditions(cfg)
    fig.suptitle(f"GO Enrichment: All Ontologies — {conds[0]} vs {conds[1]}\n"
                 f"Fisher's exact test (* p<{sig_thresh}, ** FDR<{sig_thresh})",
                 fontsize=12, y=1.005)

    go_enrichment_plot = plot_filename(cfg, "go_enrichment.png")
    save_figure(fig, os.path.join(out_dir, "go_enrichment.png"), cfg)
    log.info(f"  Saved {go_enrichment_plot}")
