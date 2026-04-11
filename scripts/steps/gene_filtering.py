"""Step 1a: Gene Filtering — remove low-expression genes.

Applies expression-based filter (CPM or FPKM) per condition and exports
filtered count/normalized matrices.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core.config_data import (
    get_sample_condition,
    get_sample_display_name,
    get_sample_manifest,
    get_sample_names,
)
from ..core.config_runtime import get_conditions, get_path, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, get_recipe_colors, plot_filename, COL2

log = setup_logger("normalization.gene_filtering")


def _parse_gene_lengths_from_gtf(gtf_path):
    """Parse gene lengths from GTF by merging exon intervals per gene.

    Returns a Series: gene_id → effective length (bp).
    """
    gene_exons = {}  # gene_id → list of (start, end)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            start, end = int(parts[3]), int(parts[4])
            attrs = parts[8]
            m = re.search(r'gene_id "([^"]+)"', attrs)
            if not m:
                continue
            gid = m.group(1)
            gene_exons.setdefault(gid, []).append((start, end))

    lengths = {}
    for gid, exons in gene_exons.items():
        # Merge overlapping exon intervals
        exons.sort()
        merged = [exons[0]]
        for s, e in exons[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        lengths[gid] = sum(e - s + 1 for s, e in merged)

    return pd.Series(lengths, name="gene_length")


def _normalize(raw, method, gene_lengths=None):
    """Compute CPM, FPKM, or TPM from raw counts.

    Parameters
    ----------
    raw : DataFrame
        Raw count matrix (genes × samples).
    method : str
        "cpm", "fpkm", or "tpm".
    gene_lengths : Series, optional
        Gene lengths in bp, required for FPKM and TPM.

    Returns
    -------
    DataFrame
        Normalized expression matrix.
    """
    if method in ("fpkm", "tpm"):
        if gene_lengths is None:
            raise ValueError(f"gene_lengths required for {method.upper()} normalization")
        gl = gene_lengths.reindex(raw.index)
        missing = gl.isna().sum()
        if missing > 0:
            log.warning(f"  {missing} genes missing from GTF — excluded from {method.upper()}")
            gl = gl.dropna()
            raw = raw.loc[gl.index]
        if method == "fpkm":
            lib_sizes = raw.sum(axis=0)
            # FPKM = reads * 1e9 / (library_size * gene_length)
            norm = raw.div(lib_sizes, axis=1).multiply(1e9).div(gl, axis=0)
        else:
            # TPM = (reads / gene_length) / sum(reads / gene_length) * 1e6
            rate = raw.div(gl, axis=0)
            norm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    else:
        lib_sizes = raw.sum(axis=0)
        norm = raw.div(lib_sizes, axis=1) * 1e6
    return norm


def run(cfg):
    """Run gene filtering step.

    Requires cfg['_data'] populated by normalization step.
    """
    log.info("── Step 1a: Gene Filtering ──")
    apply_style(cfg)

    samples = get_sample_names(cfg)
    out_dir = ensure_output_dir(cfg, "01_normalization/gene_filtering")
    norm_cfg = cfg["normalization"]
    filt_cfg = norm_cfg["gene_filtering"]
    method = norm_cfg.get("method", "cpm").lower()
    # Support both old "cpm_min" and new "min_expr" keys
    min_expr = filt_cfg.get("min_expr", filt_cfg.get("cpm_min", 1.0))
    min_samples = filt_cfg["min_samples"]
    method_label = method.upper()

    # Get data from previous step
    data = cfg.get("_data", {})
    raw = data.get("raw_counts")
    if raw is None:
        counts_path = get_path(cfg, "gene_counts")
        raw = pd.read_csv(counts_path, sep="\t", index_col=0)
        raw = raw[samples["all"]].copy()

    total_before = raw.shape[0]
    # Capture original library sizes before any gene subsetting
    lib_sizes_original = raw.sum(axis=0)

    # Load gene lengths if FPKM or TPM
    gene_lengths = None
    if method in ("fpkm", "tpm"):
        gtf_path = get_path(cfg, "reference_gtf")
        log.info(f"  Loading gene lengths from GTF: {gtf_path}")
        gene_lengths = _parse_gene_lengths_from_gtf(gtf_path)
        log.info(f"  Gene lengths loaded: {len(gene_lengths)} genes")

    # Normalize for filtering
    expr = _normalize(raw, method, gene_lengths)
    # Update raw to match (may have dropped genes missing from GTF)
    if method in ("fpkm", "tpm"):
        raw = raw.loc[expr.index]
        total_before = raw.shape[0]

    # Filter: gene must pass expression threshold in >= min_samples per at least one condition
    conditions = get_conditions(cfg)
    keep_mask = pd.Series(False, index=raw.index)
    for cond in conditions:
        cond_samples = samples[cond]
        cond_pass = (expr[cond_samples] >= min_expr).sum(axis=1) >= min_samples
        keep_mask = keep_mask | cond_pass

    raw_filt = raw.loc[keep_mask].copy()
    total_after = raw_filt.shape[0]
    removed = total_before - total_after

    log.info(f"  Filter: {method_label} >= {min_expr} in >= {min_samples} samples/condition")
    log.info(f"  Before: {total_before} genes → After: {total_after} genes ({removed} removed)")

    # Re-normalize filtered data
    expr_filt = _normalize(raw_filt, method, gene_lengths)
    pc = norm_cfg.get("pseudocount", 1)
    log2_filt = np.log2(expr_filt + pc)

    # Save
    raw_filt.to_csv(os.path.join(out_dir, "filtered_counts.tsv"), sep="\t")
    expr_filt.to_csv(os.path.join(out_dir, f"filtered_{method}.tsv"), sep="\t")
    log2_filt.to_csv(os.path.join(out_dir, f"filtered_log2{method}.tsv"), sep="\t")

    # Save gene lengths if FPKM
    if gene_lengths is not None:
        gl_filt = gene_lengths.reindex(raw_filt.index).dropna().astype(int)
        gl_filt.to_csv(os.path.join(out_dir, "gene_lengths.tsv"), sep="\t", header=True)

    # Update in-memory data (use canonical keys for downstream compatibility)
    if "_data" not in cfg:
        cfg["_data"] = {}
    cfg["_data"]["filtered_counts"] = raw_filt
    cfg["_data"][f"filtered_{method}"] = expr_filt
    cfg["_data"][f"filtered_log2{method}"] = log2_filt

    # Summary
    summary_lines = [
        f"Gene Filtering Summary",
        f"=" * 40,
        f"Method: {method_label}",
        f"Threshold: {method_label} >= {min_expr} in >= {min_samples} samples per condition",
        f"Genes before filtering: {total_before}",
        f"Genes after filtering:  {total_after}",
        f"Genes removed:          {removed}",
        f"Retention rate:          {total_after/total_before:.1%}",
    ]
    summary_text = "\n".join(summary_lines)
    with open(os.path.join(out_dir, "filtering_summary.txt"), "w") as f:
        f.write(summary_text)

    colors = get_recipe_colors(cfg)
    bar_colors = []
    for s in samples["all"]:
        cond = get_sample_condition(cfg, s)
        bar_colors.append(colors.get(cond, "#999"))

    lib_sizes = lib_sizes_original
    lib_sizes_filt = raw_filt.sum(axis=0)

    # Plot 1: gene_filtering.png — gene count bar + expression violin
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 5), gridspec_kw={"width_ratios": [2, 5]})

    # Left: before/after gene count bar chart
    ax = axes[0]
    ax.bar([0, 0.5], [total_before, total_after],
           color=["#BBBBBB", "#2166AC"], width=0.15, edgecolor="none")
    ax.set_xticks([0, 0.5])
    ax.set_xticklabels(["Before", "After"])
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Gene filtering\n({method_label} >= {min_expr}, >= {min_samples} samples/condition)")
    for x, v in zip([0, 0.5], [total_before, total_after]):
        ax.text(x, v + 50, f"{v:,}", ha="center", va="bottom", fontsize=9)
    style_axes(ax)

    # Right: violin plot of expression distribution per sample (after filtering)
    ax = axes[1]
    data_violin = [log2_filt[s].values for s in samples["all"]]
    parts = ax.violinplot(data_violin, positions=range(len(samples["all"])),
                          showmedians=True, showextrema=False, widths=0.7)
    for i, pc_body in enumerate(parts["bodies"]):
        cond = get_sample_condition(cfg, samples["all"][i])
        pc_body.set_facecolor(colors.get(cond, "#999999"))
        pc_body.set_alpha(0.7)
        pc_body.set_edgecolor("none")
    parts["cmedians"].set_color("white")
    parts["cmedians"].set_linewidth(1.0)
    ax.set_xticks(range(len(samples["all"])))
    ax.set_xticklabels(
        [get_sample_display_name(cfg, s) for s in samples["all"]],
        rotation=45,
        ha="right",
    )
    ax.set_ylabel(f"log2({method_label} + 1)")
    ax.set_title("Expression distribution (filtered)")
    style_axes(ax)

    gene_filtering_plot = plot_filename(cfg, "gene_filtering.png")
    save_figure(fig, os.path.join(out_dir, "gene_filtering.png"), cfg)
    log.info(f"  Saved {gene_filtering_plot}")

    # TSV: library size before/after filtering per sample
    lib_size_df = pd.DataFrame({
        "sample": samples["all"],
        "condition": [get_sample_condition(cfg, s) for s in samples["all"]],
        "display_label": [get_sample_display_name(cfg, s) for s in samples["all"]],
        "before_M": (lib_sizes.values / 1e6).round(3),
        "after_M": (lib_sizes_filt.values / 1e6).round(3),
        "change_M": ((lib_sizes_filt.values - lib_sizes.values) / 1e6).round(3),
        "change_pct": ((lib_sizes_filt.values - lib_sizes.values) / lib_sizes.values * 100).round(2),
    })
    lib_size_df.to_csv(os.path.join(out_dir, "library_size.tsv"), sep="\t", index=False)

    file_descs = {
        "filtered_counts.tsv": "Raw read counts after removing low-expression genes (genes x samples).",
        f"filtered_{method}.tsv": f"{method_label}-normalized expression values from filtered counts.",
        f"filtered_log2{method}.tsv": f"log2({method_label} + 1) transformed values; primary input for downstream analysis.",
        "library_size.tsv": "Per-sample library size (M reads) before and after gene filtering with change metrics.",
        "filtering_summary.txt": f"Filtering statistics: {method_label} threshold, gene counts before/after, retention rate.",
        gene_filtering_plot: f"Left: gene count bar chart (before vs after). Right: log2({method_label}+1) violin plot per sample.",
    }
    if gene_lengths is not None:
        file_descs["gene_lengths.tsv"] = "Gene lengths (bp) from GTF used for FPKM calculation."
    write_readme(out_dir, "1a", "Gene Filtering", file_descs)

    manifest = get_sample_manifest(cfg)
    note = ""
    if manifest is not None:
        note = (
            "\n- Exploration mode: repeated samples retained from `sample_manifest.tsv`\n"
            "- Filtering and visualization operate on physical samples; no DE inference is implied\n"
        )

    doc = f"""\
### Parameters

- **Normalization method**: {method_label}
- **Expression threshold**: >= {min_expr}
- **Min samples per condition**: {min_samples}
{note}

### Results

- Genes before filtering: {total_before:,}
- Genes after filtering: {total_after:,}
- Genes removed: {removed:,}
- Retention rate: {total_after/total_before:.1%}

{img(cfg, "01_normalization/gene_filtering", "gene_filtering.png", "Gene filtering: before/after counts and per-sample expression distributions")}
"""
    update_section(cfg, "1a", "Gene Filtering", doc)

    log.info("  Saved library_size.tsv")
