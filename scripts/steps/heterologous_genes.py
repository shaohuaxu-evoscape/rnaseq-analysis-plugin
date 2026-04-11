"""Step 4c: Heterologous Gene Analysis — target-condition expression of 8 engineered genes.

Two figures:
  fig1: 4x2 timecourse (CPM) for all 8 heterologous genes
  fig2: 2x3 expression level vs log2FC scatter (mean + per timepoint)

Plus summary TSV with per-gene stats and per-timepoint CPM/FC values.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

from ..core.config_data import get_filtered_expr, get_norm_method, get_sample_names, require_sample_manifest
from ..core.config_runtime import get_conditions, get_path, get_timepoints, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, get_recipe_colors, plot_filename, COL2

log = setup_logger("advanced_analysis.heterologous_genes")

# Configure via analysis_case.yaml: advanced_analysis.heterologous_genes.gene_list
# Each entry: {name: "GeneName", function: "description", color: "#hex"}
_DEFAULT_GENE_LIST = []

# Fusion rules: if the fusion gene exists, drop the parts to avoid double-counting
# Each entry: {fusion: "FusedAB", parts: ["GeneA", "GeneB"]}
_DEFAULT_FUSION_RULES = []


def _load_gene_config(adv_cfg):
    """Load gene list, function map, color map, and fusion rules from config or defaults."""
    gene_list = adv_cfg.get("gene_list", _DEFAULT_GENE_LIST)
    fusion_rules = adv_cfg.get("fusion_rules", _DEFAULT_FUSION_RULES)

    genes = [g["name"] for g in gene_list]
    gene_functions = {g["name"]: g.get("function", "") for g in gene_list}
    gene_colors = {g["name"]: g.get("color", "#666666") for g in gene_list}

    return genes, gene_functions, gene_colors, fusion_rules


def run(cfg):
    adv_cfg = cfg.get("advanced_analysis", {}).get("heterologous_genes", {})
    if adv_cfg.get("enabled") is False:
        log.info("Heterologous gene analysis disabled, skipping")
        return

    log.info("── Step 4c: Heterologous Gene Analysis ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    timepoints = get_timepoints(cfg)
    out_dir = ensure_output_dir(cfg, "04_advanced_analysis/heterologous_genes")
    rc = get_recipe_colors(cfg)
    norm_label = get_norm_method(cfg).upper()  # "CPM" or "FPKM"

    # Load gene config from analysis_case.yaml or defaults
    GENES, GENE_FUNCTIONS, GENE_COLORS, fusion_rules = _load_gene_config(adv_cfg)

    # Load filtered expression and DE results
    expr = _load_expr(cfg)
    de_results = _load_de_results(cfg)
    cond1, cond2 = get_conditions(cfg)[:2]
    sample_map = _sample_map_by_timepoint(cfg, timepoints)

    # Resolve gene list: apply fusion rules to avoid double-counting
    active_genes = list(GENES)
    for rule in fusion_rules:
        fusion_id = rule["fusion"]
        parts = set(rule.get("parts", []))
        has_fusion = fusion_id in expr.index
        has_parts = parts & set(expr.index)
        if has_fusion:
            active_genes = [g for g in active_genes if g not in parts]
        else:
            active_genes = [g for g in active_genes if g != fusion_id]

    # Compute per-gene stats
    gene_data = {}
    for gene in active_genes:
        if gene not in expr.index:
            log.warning(f"  Gene {gene} not found in count matrix, skipping")
            continue
        cond1_vals = [_mean_expr_at_tp(expr, gene, sample_map[cond1][tp]) for tp in timepoints]
        cond2_vals = [_mean_expr_at_tp(expr, gene, sample_map[cond2][tp]) for tp in timepoints]
        fcs = [np.log2((v2 + 1) / (v1 + 1)) for v1, v2 in zip(cond1_vals, cond2_vals)]
        mean_expr = np.mean(cond1_vals + cond2_vals)
        mean_fc = np.mean(fcs)
        pval = float(de_results.loc[gene, "pvalue"]) if gene in de_results.index else np.nan
        is_de = bool(de_results.loc[gene, "is_de"]) if gene in de_results.index else False
        gene_data[gene] = {
            "cond1": cond1_vals, "cond2": cond2_vals, "fcs": fcs,
            "mean_expr": mean_expr, "mean_fc": mean_fc, "p_value": pval, "is_de": is_de,
        }

    if not gene_data:
        log.warning("  No heterologous genes found, skipping")
        return

    found = list(gene_data.keys())
    log.info(f"  Found {len(found)}/{len(active_genes)} heterologous genes")
    for gene in found:
        d = gene_data[gene]
        sig = "*" if d["is_de"] else ""
        log.info(f"    {gene:>10}: {norm_label}={d['mean_expr']:>7.1f}, "
                 f"FC={d['mean_fc']:+.3f}, p={d['p_value']:.4f}{sig}")

    # Save summary TSV
    _save_summary(gene_data, timepoints, out_dir, GENE_FUNCTIONS)

    # Figures
    _plot_timecourse(gene_data, timepoints, rc, out_dir, cfg, GENE_FUNCTIONS)
    _plot_expression_vs_fc(gene_data, timepoints, expr, out_dir, cfg, GENE_COLORS)

    write_readme(out_dir, "4c", "Heterologous Gene Analysis", {
        "heterologous_genes_summary.tsv": f"Per-gene summary with function, mean {norm_label}, mean log2FC, p-value, significance, and per-timepoint condition values",
        plot_filename(cfg, "fig_heterologous_timecourse.png"): f"4x2 timecourse panels showing {norm_label} for each heterologous gene with condition comparison",
        plot_filename(cfg, "fig_heterologous_expr_vs_fc.png"): f"2x3 scatter panels of expression level (log10 {norm_label}) vs log2FC per timepoint and mean",
    })

    n_sig = sum(1 for g in gene_data.values() if g["is_de"])
    doc = (
        "### Parameters\n\n"
        f"- **{len(active_genes)}** heterologous genes queried ({len(found)} found in count matrix)\n"
        f"- {norm_label} aggregated within each target condition and timepoint by sample mean\n"
        f"- log2FC = log2(({cond2}+1)/({cond1}+1)) using those timepoint means\n"
        "- Significance reuses the global DE result from step 3a\n\n"
        "### Results\n\n"
        f"- **{len(found)}** genes detected, **{n_sig}** significantly different (FDR < 0.05, |FC| > 1.5-fold)\n"
        + "".join(
            f"- {g}: mean {norm_label}={gene_data[g]['mean_expr']:.0f}, "
            f"FC={gene_data[g]['mean_fc']:+.3f}, p={gene_data[g]['p_value']:.4f}\n"
            for g in found
        )
        + f"\n{img(cfg, '04_advanced_analysis/heterologous_genes', 'fig_heterologous_timecourse.png', f'Heterologous gene {norm_label} timecourse by target condition')}\n\n"
        f"{img(cfg, '04_advanced_analysis/heterologous_genes', 'fig_heterologous_expr_vs_fc.png', 'Expression level vs log2FC scatter')}\n"
    )
    update_section(cfg, "4c", "Heterologous Gene Analysis", doc)

    log.info("  Heterologous gene analysis complete")


def _load_expr(cfg):
    """Load filtered normalized expression (CPM or FPKM)."""
    return get_filtered_expr(cfg, log2=False)


def _load_de_results(cfg):
    """Load DE results indexed by gene."""
    out_root = get_path(cfg, "output_dir")
    path = os.path.join(out_root, "03_differential_analysis", "de_screening", "de_results_all.tsv")
    return pd.read_csv(path, sep="\t", index_col=0)


def _sample_map_by_timepoint(cfg, timepoints):
    """Map target condition and timepoint to contributing sample ids."""
    manifest = require_sample_manifest(
        cfg,
        required_columns={"sample_id", "target_condition", "timepoint", "batch_id"},
    )

    mapping = {}
    for cond in get_conditions(cfg)[:2]:
        mapping[cond] = {}
        sub = manifest[manifest["target_condition"] == cond]
        for tp in timepoints:
            ids = sub.loc[sub["timepoint"] == tp, "sample_id"].tolist()
            if not ids:
                raise ValueError(f"No samples found for {cond} at {tp}h")
            mapping[cond][tp] = ids
    return mapping


def _mean_expr_at_tp(expr_df, gene, sample_ids):
    return float(expr_df.loc[gene, sample_ids].astype(float).mean())


def _save_summary(gene_data, timepoints, out_dir, gene_functions):
    """Save per-gene summary TSV."""
    rows = []
    for gene in gene_data:
        d = gene_data[gene]
        row = {
            "gene": gene,
            "function": gene_functions.get(gene, ""),
            "mean_expr": round(d["mean_expr"], 1),
            "mean_log2FC": round(d["mean_fc"], 3),
            "p_value": round(d["p_value"], 4),
            "significant": d["is_de"],
            "direction_in_condition2": ("UP" if d["mean_fc"] > 0.1
                                        else ("DOWN" if d["mean_fc"] < -0.1
                                              else "neutral")),
        }
        for i, tp in enumerate(timepoints):
            row[f"condition1_{tp}h"] = round(d["cond1"][i], 1)
            row[f"condition2_{tp}h"] = round(d["cond2"][i], 1)
            row[f"log2FC_{tp}h"] = round(d["fcs"][i], 3)
        rows.append(row)

    df = pd.DataFrame(rows)
    path = os.path.join(out_dir, "heterologous_genes_summary.tsv")
    df.to_csv(path, sep="\t", index=False)
    log.info(f"  Saved heterologous_genes_summary.tsv ({len(rows)} genes)")


def _plot_timecourse(gene_data, timepoints, rc, out_dir, cfg, gene_functions):
    """4x2 timecourse figure: expression per gene, condition1 vs condition2."""
    found = list(gene_data.keys())
    nrows = (len(found) + 1) // 2
    fig, axes = plt.subplots(nrows, 2, figsize=(COL2, nrows * 3.5))
    if nrows == 1:
        axes = axes.reshape(1, -1)
    cond1, cond2 = get_conditions(cfg)[:2]
    fig.suptitle(f"Heterologous Gene Expression — {cond1} vs {cond2}",
                 fontsize=14, fontweight="bold", y=0.98)

    for idx, gene in enumerate(found):
        ax = axes[idx // 2, idx % 2]
        d = gene_data[gene]

        ax.plot(timepoints, d["cond1"], "o-", color=rc.get(cond1, "#2166AC"),
                ms=6, lw=1.5, label=cond1)
        ax.plot(timepoints, d["cond2"], "s-", color=rc.get(cond2, "#B2182B"),
                ms=6, lw=1.5, label=cond2)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel(get_norm_method(cfg).upper())
        ax.set_xticks(timepoints)

        sig = "DE" if d["is_de"] else f'p={d["p_value"]:.2f}'
        direction = ("UP" if d["mean_fc"] > 0.1
                     else ("DOWN" if d["mean_fc"] < -0.1 else "neutral"))
        ax.set_title(f"{gene} — {gene_functions.get(gene, '')}")
        ax.text(0.97, 0.95,
                f'mean FC={d["mean_fc"]:+.2f}\n{sig}\n{direction} in {cond2}',
                transform=ax.transAxes, fontsize=8, va="top", ha="right",
                bbox=dict(boxstyle="round", facecolor="#F5F5F5", alpha=0.8))
        ax.legend(fontsize=8, loc="upper left")
        style_axes(ax)

    # Hide unused axes
    for idx in range(len(found), nrows * 2):
        axes[idx // 2, idx % 2].set_visible(False)

    fig.tight_layout()
    heterologous_timecourse_plot = plot_filename(cfg, "fig_heterologous_timecourse.png")
    save_figure(fig, os.path.join(out_dir, "fig_heterologous_timecourse.png"), cfg)
    log.info(f"  Saved {heterologous_timecourse_plot}")


def _plot_expression_vs_fc(gene_data, timepoints, expr, out_dir, cfg, gene_colors):
    """2x3 scatter: expression level vs log2FC per timepoint + mean."""
    found = list(gene_data.keys())
    fc_thresh = cfg.get("differential_analysis", {}).get("de_screening", {}).get("log2fc_threshold", 0.585)
    norm_label = get_norm_method(cfg).upper()

    cond1, cond2 = get_conditions(cfg)[:2]
    tp_scatter = {}
    for tp in timepoints:
        rows = []
        for gene in found:
            idx = timepoints.index(tp)
            v1 = gene_data[gene]["cond1"][idx]
            v2 = gene_data[gene]["cond2"][idx]
            avg = (v1 + v2) / 2
            fc = np.log2((v2 + 1) / (v1 + 1))
            rows.append({"gene": gene, "mean_expr": avg, "log2FC": fc})
        tp_scatter[tp] = pd.DataFrame(rows)

    # Mean data
    mean_df = pd.DataFrame([
        {"gene": g, "mean_expr": gene_data[g]["mean_expr"],
         "log2FC": gene_data[g]["mean_fc"]}
        for g in found
    ])

    fig, axes = plt.subplots(2, 3, figsize=(COL2, 8))
    fig.suptitle(f"Expression Level vs log2FC ({cond2}/{cond1}) — Heterologous Genes",
                 fontsize=14, fontweight="bold", y=0.98)

    _scatter_panel(axes[0, 0], mean_df, "Mean (all timepoints)", norm_label=norm_label, show_ylabel=True, fc_thresh=fc_thresh, gene_colors=gene_colors)
    for i, tp in enumerate(timepoints):
        row, col = (i + 1) // 3, (i + 1) % 3
        _scatter_panel(axes[row, col], tp_scatter[tp], f"{tp}h",
                       norm_label=norm_label, show_ylabel=(col == 0), fc_thresh=fc_thresh, gene_colors=gene_colors)

    fig.tight_layout()
    heterologous_expr_fc_plot = plot_filename(cfg, "fig_heterologous_expr_vs_fc.png")
    save_figure(fig, os.path.join(out_dir, "fig_heterologous_expr_vs_fc.png"), cfg)
    log.info(f"  Saved {heterologous_expr_fc_plot}")


def _scatter_panel(ax, df, title, norm_label="FPKM", show_ylabel=True, fc_thresh=0.585, gene_colors=None):
    """Single scatter panel: log10(expression) vs log2FC."""
    _gc = gene_colors or {}
    for _, row in df.iterrows():
        ax.scatter(np.log10(row["mean_expr"]), row["log2FC"],
                   c=[_gc.get(row["gene"], "#666")], s=40, zorder=3,
                   edgecolors="white", linewidths=0.5)
        offset_y = 0.05 if row["log2FC"] >= 0 else -0.12
        ax.text(np.log10(row["mean_expr"]) + 0.05, row["log2FC"] + offset_y,
                row["gene"], fontsize=8,
                va="bottom" if row["log2FC"] >= 0 else "top")

    ax.axhline(0, color="black", ls="-", lw=1.0)
    ax.axhline(fc_thresh, color="gray", ls=":", lw=0.8, alpha=0.4)
    ax.axhline(-fc_thresh, color="gray", ls=":", lw=0.8, alpha=0.4)

    rho, pval = spearmanr(df["mean_expr"], df["log2FC"])
    ax.text(0.05, 0.95, f"rho={rho:.2f}, p={pval:.3f}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round", facecolor="#FFF8E1", alpha=0.8))

    ax.set_xlabel(f"log10(Mean {norm_label})")
    if show_ylabel:
        ax.set_ylabel("log2FC")
    ax.set_title(title)
    ax.set_xlim(0.8, 4.2)
    ax.set_ylim(-1.8, 1.2)
    style_axes(ax)
