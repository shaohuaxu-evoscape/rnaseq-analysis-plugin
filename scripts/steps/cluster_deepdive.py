"""Step 5b: Cluster Deep-Dive — per-cluster enrichment analysis."""

import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from ..core.config_data import get_filtered_expr, get_sample_manifest, get_sample_names
from ..core.config_runtime import ensure_output_dir, get_conditions, get_path, get_timepoints, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, get_recipe_colors, plot_filename, PALETTE
from ..core.stats import bh_fdr

log = setup_logger("cluster_analysis.cluster_deepdive")

# Default clusters to analyze (overridden by config cluster_deepdive.target_clusters)
_DEFAULT_TARGET_CLUSTERS = [3, 4, 6, 7]


def _load_go_mapping(cfg):
    """Load cached GO mapping from step 3b."""
    out_root = get_path(cfg, "output_dir")
    cache_dir = os.path.join(out_root, "03_differential_analysis", "go_enrichment", "cache")

    # accession → gene name
    acc_df = pd.read_csv(os.path.join(cache_dir, "uniprot_genes.tsv"), sep="\t")
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")
    acc_to_gene = {}
    for acc, names in zip(acc_df.iloc[:, 0].astype(str), acc_df.iloc[:, 1].astype(str)):
        acc = acc.strip()
        for gn in names.strip().split():
            if gene_id_prefix and gn.startswith(gene_id_prefix):
                acc_to_gene[acc] = gn
                break

    # accession → GO terms
    go_df = pd.read_csv(os.path.join(cache_dir, "uniprot_go.tsv"), sep="\t")
    gene_go = {}
    go_terms = {}

    # Map column names to ontology codes
    ontology_cols = {}
    for col in go_df.columns:
        col_lower = col.lower()
        if "biological process" in col_lower:
            ontology_cols[col] = "BP"
        elif "molecular function" in col_lower:
            ontology_cols[col] = "MF"
        elif "cellular component" in col_lower:
            ontology_cols[col] = "CC"

    for rec in go_df.to_dict("records"):
        acc = str(rec[go_df.columns[0]]).strip()
        gene = acc_to_gene.get(acc)
        if not gene:
            continue
        if gene not in gene_go:
            gene_go[gene] = set()

        for col_name, ontology in ontology_cols.items():
            val = str(rec[col_name]).strip() if pd.notna(rec[col_name]) else ""
            if not val or val == "nan":
                continue
            for entry in val.split("; "):
                if "[GO:" not in entry:
                    continue
                bracket_start = entry.rfind("[GO:")
                go_id = entry[bracket_start + 1:entry.rfind("]")]
                term_name = entry[:bracket_start].strip()
                gene_go[gene].add(go_id)
                go_terms[go_id] = (ontology, term_name)

    term_genes = {}
    for gene, terms in gene_go.items():
        for go_id in terms:
            term_genes.setdefault(go_id, set()).add(gene)

    return gene_go, go_terms, term_genes


def _load_kegg_mapping(cfg):
    """Load cached KEGG mapping from step 3c."""
    out_root = get_path(cfg, "output_dir")
    cache_dir = os.path.join(out_root, "03_differential_analysis", "kegg_enrichment", "cache")

    # KEGG gene → UniProt accession (KEGG cache has no header)
    conv_df = pd.read_csv(os.path.join(cache_dir, "kegg_uniprot_conv.tsv"), sep="\t", header=None)
    conv_df.columns = ["kegg_gene", "uniprot"]

    # UniProt accession → gene name
    acc_df = pd.read_csv(os.path.join(cache_dir, "uniprot_genes.tsv"), sep="\t")
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")
    acc_to_gene = {}
    for acc, names in zip(acc_df.iloc[:, 0].astype(str), acc_df.iloc[:, 1].astype(str)):
        acc = acc.strip()
        for gn in names.strip().split():
            if gene_id_prefix and gn.startswith(gene_id_prefix):
                acc_to_gene[acc] = gn
                break

    # KEGG gene → gene name
    kegg_to_gene = {}
    for kegg_gene, uniprot in zip(conv_df["kegg_gene"].astype(str), conv_df["uniprot"].astype(str)):
        up_acc = uniprot.replace("up:", "")
        gene = acc_to_gene.get(up_acc)
        if gene:
            kegg_to_gene[kegg_gene] = gene

    # Pathway → gene set (KEGG cache has no header)
    pw_df = pd.read_csv(os.path.join(cache_dir, "kegg_pathway_genes.tsv"), sep="\t", header=None)
    pw_df.columns = ["pathway", "kegg_gene"]

    pathway_genes = {}
    for pw_raw, kg in zip(pw_df["pathway"].astype(str), pw_df["kegg_gene"].astype(str)):
        pw = pw_raw.replace("path:", "")
        gene = kegg_to_gene.get(kg)
        if gene:
            pathway_genes.setdefault(pw, set()).add(gene)

    # Pathway names — KEGG cache has no header row
    pathway_names = {}
    pw_list = pd.read_csv(os.path.join(cache_dir, "kegg_pathway_list.tsv"), sep="\t", header=None)
    pw_list.columns = ["pathway", "name"]
    for pw_id, name in zip(pw_list["pathway"].astype(str), pw_list["name"].astype(str)):
        pw_id = pw_id.strip().replace("path:", "")
        name = name.strip().split(" - ")[0].strip()
        pathway_names[pw_id] = name

    return pathway_genes, pathway_names


def _run_enrichment(test_genes, background_genes, term_genes, min_size=5, max_size=500):
    """Fisher's exact test for enrichment of test_genes in each term."""
    test_set = set(test_genes)
    bg_set = set(background_genes)
    N = len(bg_set)
    K = len(test_set & bg_set)

    results = []
    for term_id, tgenes in term_genes.items():
        tgenes_in_bg = tgenes & bg_set
        n = len(tgenes_in_bg)
        if n < min_size or n > max_size:
            continue
        overlap = test_set & tgenes_in_bg
        k = len(overlap)
        if k == 0:
            continue

        cell_d = N - K - n + k
        if cell_d < 0:
            continue
        table = [[k, K - k], [n - k, cell_d]]
        _, pval = fisher_exact(table, alternative="greater")
        fold = (k / K) / (n / N) if K > 0 and N > 0 else 0

        results.append({
            "term_id": term_id, "overlap": k, "term_size": n,
            "cluster_in_bg": K, "bg_size": N,
            "gene_ratio": k / n, "fold_enrichment": fold,
            "pvalue": pval, "overlap_genes": ";".join(sorted(overlap)),
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results).sort_values("pvalue")
    df["fdr"] = bh_fdr(df["pvalue"].values)
    return df


def run(cfg):
    adv_cfg = cfg.get("cluster_analysis", {}).get("cluster_deepdive", {})
    if adv_cfg.get("enabled") is False:
        log.info("Cluster deep-dive disabled, skipping")
        return

    log.info("── Step 5b: Cluster Deep-Dive ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    manifest = get_sample_manifest(cfg)
    out_dir = ensure_output_dir(cfg, "05_cluster_analysis/cluster_deepdive")
    conditions = get_conditions(cfg)
    timepoints = get_timepoints(cfg)
    colors = get_recipe_colors(cfg)

    # Load data
    log2cpm = _load_log2expr(cfg, samples)
    assign = pd.read_csv(
        os.path.join(get_path(cfg, "output_dir"),
                     "05_cluster_analysis", "gene_clustering", "cluster_assignments.tsv"),
        sep="\t")
    de_all = pd.read_csv(
        os.path.join(get_path(cfg, "output_dir"),
                     "03_differential_analysis", "de_screening", "de_results_all.tsv"),
        sep="\t", index_col=0)
    de_genes_df = pd.read_csv(
        os.path.join(get_path(cfg, "output_dir"),
                     "03_differential_analysis", "de_screening", "de_genes.tsv"),
        sep="\t", index_col=0)

    background = list(log2cpm.index)

    # Load annotation mappings
    log.info("  Loading GO mapping...")
    gene_go, go_terms, go_term_genes = _load_go_mapping(cfg)
    log.info(f"  GO: {len(go_terms)} terms, {len(gene_go)} annotated genes")

    log.info("  Loading KEGG mapping...")
    kegg_pathway_genes, kegg_pathway_names = _load_kegg_mapping(cfg)
    log.info(f"  KEGG: {len(kegg_pathway_names)} pathways")

    # Z-score for profile plots
    data = log2cpm.loc[assign["gene"]]
    z = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

    cond1, cond2 = get_conditions(cfg)[:2]
    cond1_samples = samples[cond1]
    cond2_samples = samples[cond2]

    # ── Per-cluster analysis ──
    cluster_summaries = []

    available_clusters = set(assign["cluster"].unique())
    configured_targets = adv_cfg.get("target_clusters")
    if configured_targets:
        target_clusters = sorted(set(configured_targets) & available_clusters)
        if not target_clusters:
            log.warning(f"  Configured target_clusters {configured_targets} not found in data; falling back to all clusters")
            target_clusters = sorted(available_clusters)
        elif len(target_clusters) < len(configured_targets):
            missing = sorted(set(configured_targets) - available_clusters)
            log.warning(f"  target_clusters {missing} not found in data, skipped")
    else:
        target_clusters = sorted(available_clusters)

    for cid in target_clusters:
        log.info(f"\n  ── Cluster {cid} ──")
        genes = list(assign.loc[assign["cluster"] == cid, "gene"])
        n_total = len(genes)

        # DE overlap
        de_in_cluster = [g for g in genes if g in de_genes_df.index]
        n_de = len(de_in_cluster)
        n_up = sum(1 for g in de_in_cluster if de_genes_df.loc[g, "direction"] == "up")
        n_down = n_de - n_up
        log.info(f"    Genes: {n_total}, DE: {n_de} ({n_up} up, {n_down} down)")

        # Z-score profile stats
        zc = z.loc[genes]
        cond1_mean = _collapse_cluster_profile(
            zc, manifest, cond1, cond1_samples, timepoints
        )
        cond2_mean = _collapse_cluster_profile(
            zc, manifest, cond2, cond2_samples, timepoints
        )
        div = cond2_mean.mean() - cond1_mean.mean()

        # GO enrichment
        go_results = _run_enrichment(genes, background, go_term_genes)
        if not go_results.empty:
            go_results["term_name"] = go_results["term_id"].map(
                lambda x: go_terms.get(x, ("", x))[1])
            go_results["ontology"] = go_results["term_id"].map(
                lambda x: go_terms.get(x, ("", ""))[0])
            n_go_sig = (go_results["fdr"] < 0.05).sum()
            log.info(f"    GO: {len(go_results)} terms tested, {n_go_sig} FDR<0.05")
        else:
            n_go_sig = 0

        # KEGG enrichment
        kegg_results = _run_enrichment(genes, background, kegg_pathway_genes)
        if not kegg_results.empty:
            kegg_results["pathway_name"] = kegg_results["term_id"].map(
                lambda x: kegg_pathway_names.get(x, x))
            n_kegg_sig = (kegg_results["fdr"] < 0.05).sum()
            log.info(f"    KEGG: {len(kegg_results)} pathways tested, {n_kegg_sig} FDR<0.05")
        else:
            n_kegg_sig = 0

        # Save per-cluster results
        cdir = ensure_output_dir(cfg, f"05_cluster_analysis/cluster_deepdive/C{cid}")
        if not go_results.empty:
            go_results.to_csv(os.path.join(cdir, "go_enrichment.tsv"), sep="\t", index=False)
        if not kegg_results.empty:
            kegg_results.to_csv(os.path.join(cdir, "kegg_enrichment.tsv"), sep="\t", index=False)

        # Save DE gene details for this cluster
        if de_in_cluster:
            de_detail = de_all.loc[de_in_cluster].copy()
            de_detail.to_csv(os.path.join(cdir, "de_genes.tsv"), sep="\t")

        cluster_summaries.append({
            "cluster": cid, "n_genes": n_total, "n_de": n_de,
            "n_up": n_up, "n_down": n_down, "div": div,
            "n_go_sig": n_go_sig, "n_kegg_sig": n_kegg_sig,
            "go_results": go_results, "kegg_results": kegg_results,
            "cond1_mean": cond1_mean, "cond2_mean": cond2_mean,
            "genes": genes, "de_in_cluster": de_in_cluster,
        })

    # Save summary table
    summary_df = pd.DataFrame([{k: v for k, v in s.items()
                                 if k not in ("go_results", "kegg_results", "cond1_mean",
                                              "cond2_mean", "genes", "de_in_cluster")}
                                for s in cluster_summaries])
    summary_df.to_csv(os.path.join(out_dir, "cluster_summary.tsv"), sep="\t", index=False)

    # ── Visualizations ──

    # Fig 1: Overview panel — profile + DE composition per cluster
    _plot_overview(cluster_summaries, timepoints, colors, out_dir, cfg)

    # Fig 2: GO enrichment comparison across clusters
    _plot_go_comparison(cluster_summaries, out_dir, cfg)

    # Fig 3: KEGG enrichment comparison across clusters
    _plot_kegg_comparison(cluster_summaries, out_dir, cfg)

    write_readme(out_dir, "5b", "Cluster Deep-Dive", {
        "cluster_summary.tsv": "Per-cluster summary: gene count, DE overlap (up/down), divergence score, and number of significant GO/KEGG terms",
        plot_filename(cfg, "fig_overview.png"): "Overview panel with Z-score expression profiles (top) and DE composition bar charts (bottom) per cluster",
        plot_filename(cfg, "fig_go_comparison.png"): "Top GO enrichment terms per cluster as horizontal bar charts with significance markers",
        plot_filename(cfg, "fig_kegg_comparison.png"): "Top KEGG pathway enrichment per cluster as horizontal bar charts with significance markers",
        "C*/go_enrichment.tsv": "Per-cluster GO enrichment results (Fisher's exact test, FDR-corrected) with overlap genes",
        "C*/kegg_enrichment.tsv": "Per-cluster KEGG pathway enrichment results (Fisher's exact test, FDR-corrected) with overlap genes",
        "C*/de_genes.tsv": "DE gene details for genes within each cluster",
    })

    cluster_lines = "\n".join(
        f"- **C{s['cluster']}**: {s['n_genes']} genes, "
        f"{s['n_de']} DE ({s['n_up']} up, {s['n_down']} down), "
        f"{s['n_go_sig']} GO sig, {s['n_kegg_sig']} KEGG sig"
        for s in cluster_summaries
    )
    doc = (
        "### Method\n\n"
        f"Fisher's exact test enrichment (GO + KEGG) for clusters "
        f"{', '.join(f'C{c}' for c in target_clusters)} against filtered gene background. "
        "BH-FDR correction applied.\n\n"
        "### Results\n\n"
        f"{cluster_lines}\n\n"
        f"{img(cfg, '05_cluster_analysis/cluster_deepdive', 'fig_overview.png', 'Cluster profiles and DE composition')}\n\n"
        f"{img(cfg, '05_cluster_analysis/cluster_deepdive', 'fig_go_comparison.png', 'GO enrichment comparison across clusters')}\n\n"
        f"{img(cfg, '05_cluster_analysis/cluster_deepdive', 'fig_kegg_comparison.png', 'KEGG pathway enrichment comparison across clusters')}\n"
    )
    update_section(cfg, "5b", "Cluster Deep-Dive", doc)

    log.info("\n  Cluster deep-dive complete")


def _plot_overview(summaries, timepoints, colors, out_dir, cfg):
    """Overview: 2-row panel. Top: Z-score profiles. Bottom: DE bar chart."""
    n = len(summaries)
    cond1, cond2 = get_conditions(cfg)[:2]
    fig, axes = plt.subplots(2, n, figsize=(3.2 * n, 5.5),
                             gridspec_kw={"height_ratios": [1.3, 1]})

    for i, s in enumerate(summaries):
        cid = s["cluster"]

        # Top: profile
        ax = axes[0, i]
        ax.fill_between(timepoints,
                         s["cond1_mean"] - 0.3, s["cond1_mean"] + 0.3,
                         alpha=0.15, color=colors[cond1])
        ax.fill_between(timepoints,
                         s["cond2_mean"] - 0.3, s["cond2_mean"] + 0.3,
                         alpha=0.15, color=colors[cond2])
        ax.plot(timepoints, s["cond1_mean"], "o-", color=colors[cond1],
                ms=5, lw=1.5, label=cond1)
        ax.plot(timepoints, s["cond2_mean"], "s-", color=colors[cond2],
                ms=5, lw=1.5, label=cond2)
        ax.axhline(0, ls="--", color="grey", lw=1.0)
        ax.set_title(f"C{cid} (n={s['n_genes']})\ndiv={s['div']:+.2f}",
                     fontsize=8, fontweight="bold")
        ax.set_xticks(timepoints)
        ax.set_xticklabels([str(t) for t in timepoints], fontsize=9)
        if i == 0:
            ax.set_ylabel("Z-score")
            ax.legend(fontsize=8, loc="best")
        style_axes(ax)

        # Bottom: DE bar
        ax2 = axes[1, i]
        bars = ax2.bar(["Up", "Down", "Non-DE"],
                       [s["n_up"], s["n_down"], s["n_genes"] - s["n_de"]],
                       color=[PALETTE["red"], PALETTE["blue"], PALETTE["grey"]],
                       edgecolor="black", linewidth=0.3)
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax2.text(bar.get_x() + bar.get_width() / 2, h + 2, str(int(h)),
                         ha="center", va="bottom", fontsize=9)
        ax2.set_title(f"DE: {s['n_de']}/{s['n_genes']}")
        ax2.set_ylabel("Count" if i == 0 else "")
        style_axes(ax2)

    fig.suptitle("Cluster Deep-Dive: Expression Profiles & DE Composition",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    overview_plot = plot_filename(cfg, "fig_overview.png")
    save_figure(fig, os.path.join(out_dir, "fig_overview.png"), cfg)
    log.info(f"  Saved {overview_plot}")


def _plot_go_comparison(summaries, out_dir, cfg):
    """Top GO terms per cluster, horizontal bar chart (same style as KEGG)."""
    top_n = 10
    all_rows = []
    for s in summaries:
        go = s["go_results"]
        if go.empty:
            continue
        top = go.head(top_n).copy()
        top["cluster"] = f"C{s['cluster']}"
        all_rows.append(top)

    if not all_rows:
        log.info("  No GO results to plot")
        return

    combined = pd.concat(all_rows, ignore_index=True)

    clusters = [f"C{s['cluster']}" for s in summaries]
    n_clusters = len(clusters)

    fig, axes = plt.subplots(n_clusters, 1, figsize=(7, 2.2 * n_clusters))
    if n_clusters == 1:
        axes = [axes]

    cmap = plt.colormaps.get_cmap("tab10").resampled(n_clusters)

    for i, cname in enumerate(clusters):
        ax = axes[i]
        sub = combined[combined["cluster"] == cname].sort_values("pvalue").head(top_n)
        if sub.empty:
            ax.set_visible(False)
            continue

        sub = sub.sort_values("pvalue", ascending=False)
        names = sub["term_name"].values
        nlp = -np.log10(np.clip(sub["pvalue"].values, 1e-50, None))
        p_sig = sub["pvalue"].values < 0.05
        fdr_sig = sub["fdr"].values < 0.05

        bar_colors = [cmap(i) if fs else (*cmap(i)[:3], 0.5) if ps
                      else (*cmap(i)[:3], 0.15)
                      for ps, fs in zip(p_sig, fdr_sig)]
        y = np.arange(len(names))
        ax.barh(y, nlp, color=bar_colors, edgecolor="black", linewidth=0.3, height=0.7)

        # Overlap count + significance markers
        for j, (_, row) in enumerate(sub.iterrows()):
            marker = "**" if row["fdr"] < 0.05 else "*" if row["pvalue"] < 0.05 else ""
            ax.text(nlp[j] + 0.05, j, f"{row['overlap']}{marker}",
                    va="center", fontsize=9)

        n_p = int(p_sig.sum())
        n_fdr = int(fdr_sig.sum())
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("-log10(p)" if i == n_clusters - 1 else "")
        ax.set_title(f"{cname} — GO Enrichment  Top {len(sub)}  "
                     f"({n_fdr} FDR<0.05, {n_p} p<0.05)",
                     fontsize=8, fontweight="bold", color=cmap(i))
        style_axes(ax)

    fig.suptitle("GO Enrichment by Cluster  (* p<0.05  ** FDR<0.05)",
                 fontsize=10, y=1.01)
    fig.tight_layout()
    go_comparison_plot = plot_filename(cfg, "fig_go_comparison.png")
    save_figure(fig, os.path.join(out_dir, "fig_go_comparison.png"), cfg)
    log.info(f"  Saved {go_comparison_plot}")


def _plot_kegg_comparison(summaries, out_dir, cfg):
    """Top KEGG pathways per cluster, horizontal bar chart."""
    top_n = 10
    all_rows = []
    for s in summaries:
        kegg = s["kegg_results"]
        if kegg.empty:
            continue
        top = kegg.head(top_n).copy()
        top["cluster"] = f"C{s['cluster']}"
        all_rows.append(top)

    if not all_rows:
        log.info("  No KEGG results to plot")
        return

    combined = pd.concat(all_rows, ignore_index=True)

    clusters = [f"C{s['cluster']}" for s in summaries]
    n_clusters = len(clusters)

    # One subplot per cluster, stacked vertically
    fig, axes = plt.subplots(n_clusters, 1, figsize=(7, 2.2 * n_clusters))
    if n_clusters == 1:
        axes = [axes]

    cmap = plt.colormaps.get_cmap("tab10").resampled(n_clusters)

    for i, cname in enumerate(clusters):
        ax = axes[i]
        sub = combined[combined["cluster"] == cname].sort_values("pvalue").head(top_n)
        if sub.empty:
            ax.set_visible(False)
            continue

        sub = sub.sort_values("pvalue", ascending=False)
        names = sub["pathway_name"].values
        nlp = -np.log10(np.clip(sub["pvalue"].values, 1e-50, None))
        p_sig = sub["pvalue"].values < 0.05
        fdr_sig = sub["fdr"].values < 0.05

        bar_colors = [cmap(i) if fs else (*cmap(i)[:3], 0.5) if ps
                      else (*cmap(i)[:3], 0.15)
                      for ps, fs in zip(p_sig, fdr_sig)]
        y = np.arange(len(names))
        ax.barh(y, nlp, color=bar_colors, edgecolor="black", linewidth=0.3, height=0.7)

        # Overlap count + significance markers
        for j, (_, row) in enumerate(sub.iterrows()):
            marker = "**" if row["fdr"] < 0.05 else "*" if row["pvalue"] < 0.05 else ""
            ax.text(nlp[j] + 0.05, j, f"{row['overlap']}{marker}",
                    va="center", fontsize=9)

        n_p = int(p_sig.sum())
        n_fdr = int(fdr_sig.sum())
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("-log10(p)" if i == n_clusters - 1 else "")
        ax.set_title(f"{cname} — KEGG Enrichment  Top {len(sub)}  "
                     f"({n_fdr} FDR<0.05, {n_p} p<0.05)",
                     fontsize=8, fontweight="bold", color=cmap(i))
        style_axes(ax)

    fig.suptitle("KEGG Pathway Enrichment by Cluster  (* p<0.05  ** FDR<0.05)",
                 fontsize=10, y=1.01)
    fig.tight_layout()
    kegg_comparison_plot = plot_filename(cfg, "fig_kegg_comparison.png")
    save_figure(fig, os.path.join(out_dir, "fig_kegg_comparison.png"), cfg)
    log.info(f"  Saved {kegg_comparison_plot}")



def _load_log2expr(cfg, samples):
    return get_filtered_expr(cfg, log2=True)[samples["all"]]


def _collapse_cluster_profile(data, manifest, condition, sample_ids, timepoints):
    """Return one mean value per timepoint for a cluster profile."""
    if manifest is None:
        return data[sample_ids].mean(axis=0).to_numpy()

    rows = manifest.loc[manifest["sample_id"].isin(sample_ids)].copy()
    means = []
    for tp in timepoints:
        tp_samples = rows.loc[rows["timepoint"].astype(int) == int(tp), "sample_id"].tolist()
        if not tp_samples:
            raise ValueError(
                f"No samples found for {condition} at {tp}h in cluster deep-dive"
            )
        means.append(data[tp_samples].mean(axis=1).mean())
    return np.asarray(means, dtype=float)
