"""Step 3a: Differential Expression Screening — DESeq2 Wald test.

Uses pyDESeq2 on raw filtered counts with design ~ condition + timepoint
(+ batch when batches are crossed with conditions).  Falls back to the
legacy OLS model when pydeseq2 is not installed.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core.config_data import (
    get_filtered_counts,
    get_filtered_expr,
    get_norm_method,
    get_sample_names,
    require_sample_manifest,
)
from ..core.config_runtime import get_conditions, get_timepoints, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, plot_filename, COL1, PALETTE

log = setup_logger("differential_analysis.de_screening")


# ══════════════════════════════════════════════════════════════════════════════
# Main entry point
# ══════════════════════════════════════════════════════════════════════════════

def run(cfg):
    de_cfg = cfg["differential_analysis"]["de_screening"]
    if not de_cfg["enabled"]:
        log.info("DE screening disabled, skipping")
        return

    log.info("── Step 3a: Differential Expression ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    conditions = get_conditions(cfg)
    timepoints = get_timepoints(cfg)
    out_dir = ensure_output_dir(cfg, "03_differential_analysis/de_screening")
    fdr_thresh = de_cfg["fdr_threshold"]
    fc_thresh = de_cfg["log2fc_threshold"]
    cond1, cond2 = conditions[0], conditions[1]

    # Build sample metadata
    metadata = _build_metadata(cfg, cond1, cond2, timepoints, samples)

    # Choose DE method
    method = de_cfg.get("method", "deseq2_wald")
    if method in ("deseq2_wald", "deseq2"):
        df = _run_deseq2_wald(cfg, metadata, samples, cond1, cond2, timepoints, de_cfg)
    else:
        df = _run_ols_fallback(cfg, metadata, samples, cond1, cond2, timepoints)

    # Apply thresholds
    df["is_de"] = (df["fdr"] < fdr_thresh) & (df["mean_log2fc"].abs() > fc_thresh)
    df["direction"] = "ns"
    df.loc[df["is_de"] & (df["mean_log2fc"] > 0), "direction"] = "up"
    df.loc[df["is_de"] & (df["mean_log2fc"] < 0), "direction"] = "down"

    n_up = (df["direction"] == "up").sum()
    n_down = (df["direction"] == "down").sum()
    n_de = df["is_de"].sum()
    log.info(f"  DE genes: {n_de} ({n_up} up, {n_down} down in {cond2})")

    # Save results
    df.to_csv(os.path.join(out_dir, "de_results_all.tsv"), sep="\t")
    de_genes = df[df["is_de"]].copy()
    de_genes.to_csv(os.path.join(out_dir, "de_genes.tsv"), sep="\t")

    cfg.setdefault("_data", {})
    cfg["_data"]["de_results"] = df
    cfg["_data"]["de_genes"] = de_genes

    # Determine method label for reports
    formula = _get_formula(metadata)
    method_label = "DESeq2 Wald test" if method in ("deseq2_wald", "deseq2") else "OLS fixed-effects"

    with open(os.path.join(out_dir, "de_summary.txt"), "w") as f:
        f.write("\n".join([
            "Differential Expression Summary", "=" * 40,
            f"Method: {method_label}",
            f"Formula: {formula}",
            f"FDR < {fdr_thresh}, |log2FC| > {fc_thresh}",
            f"Tested: {len(df)}", f"DE: {n_de} (up={n_up}, down={n_down})",
        ]))

    # ── Volcano plot ──
    _plot_volcano(df, fdr_thresh, fc_thresh, cond1, cond2, out_dir, cfg)

    # ── MA plot ──
    _plot_ma(df, cond1, cond2, out_dir, cfg)

    write_readme(out_dir, "3a", "DE Screening", {
        "de_results_all.tsv": "Full DE statistics for all tested genes",
        "de_genes.tsv": "Subset containing only significant DE genes",
        "de_summary.txt": "Summary of DE analysis parameters and gene counts",
        plot_filename(cfg, "volcano_plot.png"): "Volcano plot of log2FC vs -log10(FDR)",
        plot_filename(cfg, "ma_plot.png"): "MA plot of mean expression vs log2FC",
    })

    subdir = "03_differential_analysis/de_screening"
    doc = (
        f"### Parameters\n"
        f"- Method: {method_label}\n"
        f"- Formula: `{formula}`\n"
        f"- FDR threshold: {fdr_thresh}\n"
        f"- |log2FC| threshold: {fc_thresh}\n"
        f"- Conditions: {cond1} vs {cond2}\n"
        f"- Timepoints: {', '.join(str(t) for t in timepoints)}\n\n"
        f"### Results\n"
        f"- Genes tested: {len(df)}\n"
        f"- DE genes: {n_de} ({n_up} up, {n_down} down in {cond2})\n\n"
        f"{img(cfg, subdir, 'volcano_plot.png', 'Volcano plot')}\n\n"
        f"{img(cfg, subdir, 'ma_plot.png', 'MA plot')}\n"
    )
    update_section(cfg, "3a", "DE Screening", doc)
    log.info("  DE screening complete")


# ══════════════════════════════════════════════════════════════════════════════
# Metadata & design detection
# ══════════════════════════════════════════════════════════════════════════════

def _build_metadata(cfg, cond1, cond2, timepoints, samples):
    """Build sample metadata DataFrame for DE model."""
    manifest = require_sample_manifest(
        cfg,
        required_columns={"sample_id", "target_condition", "timepoint", "batch_id"},
    )

    # Validate completeness
    counts = manifest.groupby(["target_condition", "timepoint"]).size().reset_index()
    for cond in [cond1, cond2]:
        sub = counts[counts["target_condition"] == cond]
        missing = sorted(set(timepoints) - set(sub["timepoint"].tolist()))
        if missing:
            raise ValueError(
                f"Differential analysis step 3a is missing samples for {cond} "
                f"at timepoints: {missing}"
            )
    if len(timepoints) < 2:
        raise ValueError("Differential analysis step 3a requires at least 2 shared timepoints.")

    metadata = (
        manifest.set_index("sample_id")
        .loc[samples["all"], ["target_condition", "timepoint", "batch_id"]]
        .copy()
    )
    metadata["condition"] = metadata["target_condition"]
    metadata["timepoint"] = metadata["timepoint"].astype(str)
    metadata["batch_id"] = metadata["batch_id"].astype(str)
    return metadata


def _is_crossed_design(metadata):
    """Check if any batch contributes samples to both conditions."""
    batches_per_cond = metadata.groupby("condition")["batch_id"].apply(set)
    if len(batches_per_cond) < 2:
        return False
    cond_list = list(batches_per_cond.index)
    return bool(batches_per_cond[cond_list[0]] & batches_per_cond[cond_list[1]])


def _get_formula(metadata):
    """Return the appropriate DESeq2 formula based on design."""
    if _is_crossed_design(metadata):
        return "~ batch_id + condition + timepoint"
    return "~ condition + timepoint"


# ══════════════════════════════════════════════════════════════════════════════
# DESeq2 Wald test
# ══════════════════════════════════════════════════════════════════════════════

def _run_deseq2_wald(cfg, metadata, samples, cond1, cond2, timepoints, de_cfg):
    """Run DESeq2 Wald test on raw filtered counts."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        log.warning("  pydeseq2 not installed, falling back to OLS")
        return _run_ols_fallback(cfg, metadata, samples, cond1, cond2, timepoints)

    # Load raw filtered counts (genes x samples, integer)
    raw_counts = get_filtered_counts(cfg)[samples["all"]]

    # Prepare metadata for pyDESeq2 (samples x factors)
    clinical = metadata[["condition", "timepoint"]].copy()
    if _is_crossed_design(metadata):
        clinical["batch_id"] = metadata["batch_id"]
        design_factors = ["batch_id", "condition", "timepoint"]
        log.info("  Design: crossed (batch in model)")
    else:
        design_factors = ["condition", "timepoint"]
        n_batches = metadata.groupby("condition")["batch_id"].nunique()
        log.info(f"  Design: nested (batch = biological replicate, not in model)")
        log.info(f"  Replicates: {dict(n_batches)}")

    formula = _get_formula(metadata)
    log.info(f"  Formula: {formula}")
    log.info(f"  Running DESeq2 Wald test on {raw_counts.shape[0]} genes, "
             f"{raw_counts.shape[1]} samples...")

    # pyDESeq2 expects counts as samples x genes
    # Ensure reference level for condition factor (cond1 = baseline)
    clinical["condition"] = pd.Categorical(
        clinical["condition"], categories=[cond1, cond2]
    )

    dds = DeseqDataSet(
        counts=raw_counts.T,
        metadata=clinical,
        design=formula,
        refit_cooks=True,
    )
    dds.deseq2()

    # Extract Wald test results for condition effect
    stat = DeseqStats(dds, contrast=["condition", cond2, cond1])
    stat.summary()

    results_df = stat.results_df.copy()
    log.info(f"  DESeq2 complete: {len(results_df)} genes tested")

    # Build output DataFrame with standard column names
    df = pd.DataFrame(index=results_df.index)
    df.index.name = "gene"
    df["mean_log2fc"] = results_df["log2FoldChange"].fillna(0.0)
    df["mean_expr"] = np.log2(results_df["baseMean"].fillna(0.0) + 1)
    df["stat"] = results_df["stat"].fillna(0.0)
    df["pvalue"] = results_df["pvalue"].fillna(1.0)
    df["fdr"] = results_df["padj"].fillna(1.0)

    # Add descriptive per-timepoint fold changes (observed means on log2-normalized data)
    log2expr = get_filtered_expr(cfg, log2=True)[samples["all"]]
    for tp in timepoints:
        fc_series = _timepoint_mean_log2fc(log2expr, metadata, cond1, cond2, tp)
        if isinstance(fc_series, pd.Series):
            df[f"log2fc_{tp}h"] = fc_series.reindex(df.index).fillna(0.0)
        else:
            df[f"log2fc_{tp}h"] = fc_series

    return df


# ══════════════════════════════════════════════════════════════════════════════
# OLS fallback (legacy)
# ══════════════════════════════════════════════════════════════════════════════

def _run_ols_fallback(cfg, metadata, samples, cond1, cond2, timepoints):
    """Legacy OLS fixed-effects model on log2-normalized data."""
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests

    log.info("  Running OLS fallback...")
    log2expr = get_filtered_expr(cfg, log2=True)[samples["all"]]

    # Build design matrix
    design = pd.DataFrame(index=metadata.index)
    design["Intercept"] = 1.0
    design["condition_cond2"] = (metadata["condition"] == cond2).astype(float)

    tp_cats = sorted(metadata["timepoint"].unique())
    for tp in tp_cats[1:]:
        design[f"timepoint_{tp}"] = (metadata["timepoint"] == tp).astype(float)

    if _is_crossed_design(metadata):
        batch_cats = sorted(metadata["batch_id"].unique())
        for b in batch_cats[1:]:
            design[f"batch_{b}"] = (metadata["batch_id"] == b).astype(float)

    coef_idx = design.columns.get_loc("condition_cond2")

    results = []
    for gene in log2expr.index:
        y = log2expr.loc[gene, metadata.index].astype(float).values
        fit = sm.OLS(y, design.values).fit()
        mean_log2fc = float(fit.params[coef_idx])
        p_val = float(fit.pvalues[coef_idx])
        if np.isnan(p_val):
            p_val = 1.0

        row = {
            "gene": gene,
            "mean_log2fc": mean_log2fc,
            "mean_expr": float(np.mean(y)),
            "stat": float(fit.tvalues[coef_idx]),
            "pvalue": p_val,
        }
        for tp in timepoints:
            row[f"log2fc_{tp}h"] = _timepoint_mean_log2fc(
                log2expr, metadata, cond1, cond2, tp
            )
        results.append(row)

    df = pd.DataFrame(results).set_index("gene")
    _, fdr_vals, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["fdr"] = fdr_vals
    return df


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

def _timepoint_mean_log2fc(log2expr, metadata, cond1, cond2, timepoint):
    """Observed per-gene mean difference at a timepoint: mean(cond2) - mean(cond1).

    Returns a Series indexed by gene when called for DataFrame assignment,
    or a float scalar when called for a single gene vector.
    """
    tp_str = str(int(timepoint))
    tp_mask = metadata["timepoint"] == tp_str
    c1_mask = metadata["condition"] == cond1
    c2_mask = metadata["condition"] == cond2

    c1_samples = metadata.index[tp_mask & c1_mask]
    c2_samples = metadata.index[tp_mask & c2_mask]

    if len(c1_samples) == 0 or len(c2_samples) == 0:
        return 0.0

    mean1 = log2expr[c1_samples].mean(axis=1)
    mean2 = log2expr[c2_samples].mean(axis=1)
    return mean2 - mean1


# ══════════════════════════════════════════════════════════════════════════════
# Plots
# ══════════════════════════════════════════════════════════════════════════════

def _plot_volcano(df, fdr_thresh, fc_thresh, cond1, cond2, out_dir, cfg):
    ns = df[~df["is_de"]]
    up = df[df["direction"] == "up"]
    down = df[df["direction"] == "down"]
    n_up, n_down = len(up), len(down)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))
    ax.scatter(ns["mean_log2fc"], -np.log10(ns["fdr"]),
               c=PALETTE["light_grey"], s=3, alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(up["mean_log2fc"], -np.log10(up["fdr"]),
               c=PALETTE["red"], s=4, alpha=0.7, linewidths=0, label=f"Up ({n_up})")
    ax.scatter(down["mean_log2fc"], -np.log10(down["fdr"]),
               c=PALETTE["blue"], s=4, alpha=0.7, linewidths=0, label=f"Down ({n_down})")

    ax.axhline(-np.log10(fdr_thresh), color="#999", linewidth=0.8, linestyle="--")
    ax.axvline(fc_thresh, color="#999", linewidth=0.8, linestyle="--")
    ax.axvline(-fc_thresh, color="#999", linewidth=0.8, linestyle="--")

    ax.set_xlabel("log$_2$FC")
    ax.set_ylabel("$-$log$_{10}$(FDR)")
    ax.set_title(f"{cond2} vs {cond1}")
    ax.legend(markerscale=1.5)
    style_axes(ax)
    save_figure(fig, os.path.join(out_dir, "volcano_plot.png"), cfg)
    log.info(f"  Saved {plot_filename(cfg, 'volcano_plot.png')}")


def _plot_ma(df, cond1, cond2, out_dir, cfg):
    ns = df[~df["is_de"]]
    up = df[df["direction"] == "up"]
    down = df[df["direction"] == "down"]
    n_up, n_down = len(up), len(down)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))
    ax.scatter(ns["mean_expr"], ns["mean_log2fc"],
               c=PALETTE["light_grey"], s=3, alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(up["mean_expr"], up["mean_log2fc"],
               c=PALETTE["red"], s=4, alpha=0.7, linewidths=0, label=f"Up ({n_up})")
    ax.scatter(down["mean_expr"], down["mean_log2fc"],
               c=PALETTE["blue"], s=4, alpha=0.7, linewidths=0, label=f"Down ({n_down})")
    ax.axhline(0, color="#999", linewidth=0.8)
    ax.set_xlabel("log$_2$(baseMean + 1)")
    ax.set_ylabel("log$_2$FC")
    ax.set_title("MA plot")
    ax.legend(markerscale=1.5)
    style_axes(ax)
    save_figure(fig, os.path.join(out_dir, "ma_plot.png"), cfg)
    log.info(f"  Saved {plot_filename(cfg, 'ma_plot.png')}")
