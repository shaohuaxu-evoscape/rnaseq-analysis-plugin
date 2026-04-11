"""Grid search for preprocessing and gene filtering parameter optimization.

Usage:
    python -m scripts.grid_search --mode step1 -c configs/analysis_case.yaml
    python -m scripts.grid_search --mode step0 -c configs/analysis_case.yaml
    python -m scripts.grid_search --mode step1 -c configs/analysis_case.yaml --collect-only
"""

import argparse
import copy
import itertools
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from .prepare_analysis_case import prepare_analysis_case, _set_nested
from .runner import run_pipeline
from .logger import setup_logger
from .plotting import apply_style, style_axes, save_figure

log = setup_logger("grid_search")

PIPELINE_ROOT = Path(__file__).resolve().parent.parent

# ── Parameter grids ──────────────────────────────────────────────────────────

STEP1_GRID = {
    "normalization.method": ["cpm", "tpm"],
    "normalization.gene_filtering.min_expr": [0.5, 1.0, 2.0, 5.0],
    "normalization.gene_filtering.min_samples": [1, 2, 3],
}

STEP0_GRID = {
    "preprocessing.fastp.unqualified_percent_limit": [30, 40, 50],
    "preprocessing.fastp.length_required": [36, 50, 75],
}


# ── Config helpers ───────────────────────────────────────────────────────────

def _get_nested(cfg, dotted_key, default=None):
    """Get a value from a nested dict using dot-separated key path."""
    keys = dotted_key.split(".")
    d = cfg
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def _combo_to_dirname(combo_dict):
    """Convert parameter combo to a filesystem-safe directory name.

    Example: {'normalization.method': 'tpm', ...} -> 'method_tpm__min_expr_1.0__min_samples_2'
    """
    parts = []
    for key, val in combo_dict.items():
        short_key = key.split(".")[-1]
        parts.append(f"{short_key}_{val}")
    return "__".join(parts)


def _make_grid_config(base_config_path, combo_dict, run_name, output_dir=None):
    """Create a modified analysis_case YAML for one grid search combo.

    Returns the path to the generated temporary config file.
    """
    with open(base_config_path) as f:
        cfg = yaml.safe_load(f)

    cfg = copy.deepcopy(cfg)
    cfg["run_name"] = run_name

    if output_dir is not None:
        cfg["output_dir"] = output_dir

    for key, value in combo_dict.items():
        _set_nested(cfg, key, value)

    out_root = PIPELINE_ROOT / "results" / run_name
    out_root.mkdir(parents=True, exist_ok=True)
    temp_path = out_root / "analysis_case_grid.yaml"
    with open(temp_path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)
    return str(temp_path)


# ── Metric collection ────────────────────────────────────────────────────────

def _find_pair_dir(output_root):
    """Find the condition pair subdirectory (e.g., condition1-condition2/)."""
    output_root = Path(output_root)
    if not output_root.exists():
        return None
    for d in output_root.iterdir():
        if d.is_dir() and "-" in d.name and not d.name.startswith("."):
            # Skip summary and other non-pair directories
            if d.name in ("summary", "documents"):
                continue
            return d
    return None


def _collect_metrics(output_root):
    """Parse pipeline output files and return a dict of metrics."""
    metrics = {}
    pair_dir = _find_pair_dir(output_root)
    if pair_dir is None:
        log.warning(f"  No pair directory found in {output_root}")
        return metrics

    # 1. Filtering summary
    filt_summary = pair_dir / "01_normalization" / "gene_filtering" / "filtering_summary.txt"
    if filt_summary.exists():
        text = filt_summary.read_text()
        m = re.search(r"Genes before filtering:\s+(\d+)", text)
        if m:
            metrics["genes_before_filter"] = int(m.group(1))
        m = re.search(r"Genes after filtering:\s+(\d+)", text)
        if m:
            metrics["genes_after_filter"] = int(m.group(1))
        m = re.search(r"Retention rate:\s+([\d.]+)%", text)
        if m:
            metrics["retention_rate"] = float(m.group(1))

    # 2. Correlation summary
    corr_summary = pair_dir / "02_sample_analysis" / "correlation" / "correlation_summary.txt"
    if corr_summary.exists():
        text = corr_summary.read_text()
        within_means = re.findall(r"condition\d+: .* mean=([\d.]+)", text)
        if within_means:
            metrics["within_corr_mean"] = np.mean([float(v) for v in within_means])
        m = re.search(r"Between: .* mean=([\d.]+)", text)
        if m:
            metrics["between_corr_mean"] = float(m.group(1))

    # 3. PCA variance
    pca_file = pair_dir / "02_sample_analysis" / "pca" / "pca_variance.tsv"
    if pca_file.exists():
        pca_df = pd.read_csv(pca_file, sep="\t")
        if len(pca_df) >= 2:
            metrics["pc1_variance"] = float(pca_df.iloc[0]["variance_explained"])
            metrics["pc2_variance"] = float(pca_df.iloc[1]["variance_explained"])
            metrics["pc1_pc2_cumulative"] = float(pca_df.iloc[1]["cumulative"])

    # 4. DE summary
    de_summary = pair_dir / "03_differential_analysis" / "de_screening" / "de_summary.txt"
    if de_summary.exists():
        text = de_summary.read_text()
        m = re.search(r"Tested:\s*(\d+)", text)
        if m:
            metrics["genes_tested"] = int(m.group(1))
        m = re.search(r"DE:\s*(\d+)\s*\(up=(\d+),\s*down=(\d+)\)", text)
        if m:
            metrics["n_de_genes"] = int(m.group(1))
            metrics["n_up"] = int(m.group(2))
            metrics["n_down"] = int(m.group(3))

    return metrics


def _is_combo_complete(output_root):
    """Check if a combo has already been run (DE summary exists)."""
    pair_dir = _find_pair_dir(Path(output_root))
    if pair_dir is None:
        return False
    de_summary = pair_dir / "03_differential_analysis" / "de_screening" / "de_summary.txt"
    return de_summary.exists()


# ── Step 1 grid search ───────────────────────────────────────────────────────

def run_step1_grid(base_config_path, grid=None, steps="1a-3a", collect_only=False):
    """Run gene filtering parameter grid search."""
    grid = grid or STEP1_GRID
    param_names = list(grid.keys())
    combos = list(itertools.product(*grid.values()))
    n_total = len(combos)

    log.info("=" * 60)
    log.info(f"Step 1 Grid Search: {n_total} combinations")
    log.info(f"Parameters: {', '.join(param_names)}")
    log.info("=" * 60)

    results = []
    for i, combo_values in enumerate(combos, 1):
        combo_dict = dict(zip(param_names, combo_values))
        combo_name = _combo_to_dirname(combo_dict)
        run_name = f"grid_search/step1_filtering/{combo_name}"
        output_root = PIPELINE_ROOT / "results" / run_name

        log.info(f"\n[{i}/{n_total}] {combo_name}")

        # Check if already complete
        if _is_combo_complete(output_root):
            log.info("  Already complete, collecting metrics only")
            metrics = _collect_metrics(output_root)
            metrics.update(combo_dict)
            metrics["combo_name"] = combo_name
            metrics["status"] = "ok"
            results.append(metrics)
            continue

        if collect_only:
            log.info("  Not complete, skipping (--collect-only)")
            continue

        try:
            temp_config = _make_grid_config(base_config_path, combo_dict, run_name)
            prepared = prepare_analysis_case(temp_config)
            result = run_pipeline(
                config_path=prepared["analysis_config"],
                steps=steps,
                dry_run=False,
            )
            n_err = sum(
                1 for r in result.get("results", {}).values()
                if r["status"] == "error"
            )
            metrics = _collect_metrics(output_root)
            metrics.update(combo_dict)
            metrics["combo_name"] = combo_name
            metrics["status"] = "ok" if n_err == 0 else "error"
        except Exception as e:
            log.error(f"  FAILED: {e}")
            metrics = {**combo_dict, "combo_name": combo_name, "status": "error", "error": str(e)}

        results.append(metrics)

    if not results:
        log.warning("No results collected")
        return pd.DataFrame()

    df = pd.DataFrame(results)
    summary_dir = PIPELINE_ROOT / "results" / "grid_search" / "step1_filtering" / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(summary_dir / "grid_search_results.tsv", sep="\t", index=False)
    log.info(f"\nResults saved to {summary_dir / 'grid_search_results.tsv'}")

    ok_df = df[df["status"] == "ok"]
    if len(ok_df) > 0:
        _plot_step1_results(ok_df, summary_dir)

    return df


# ── Step 0 grid search ───────────────────────────────────────────────────────

def _reuse_hisat2_index(cfg, combo_shared_dir):
    """Symlink existing HISAT2 index into the combo's shared dir.

    This avoids rebuilding the index for every grid search combo, since the
    index depends only on the reference genome — not on fastp parameters.
    The alignment step auto-discovers the index under shared_output_dir.
    """
    ref_genome = cfg.get("paths", {}).get("reference_genome", "")
    ref_name = Path(ref_genome).name
    for suffix in [".gz", ".fa", ".fasta"]:
        if ref_name.endswith(suffix):
            ref_name = ref_name[: -len(suffix)]

    # Search for an existing index directory under results/shared/
    shared_root = PIPELINE_ROOT / "results" / "shared"
    if not shared_root.exists():
        log.warning("  No existing HISAT2 index found; include step 0a")
        return

    for ht2 in shared_root.rglob(f"{ref_name}.1.ht2"):
        source_index_dir = ht2.parent
        target_index_dir = Path(combo_shared_dir) / "01_preprocessing" / "hisat2_index"
        if target_index_dir.exists() or target_index_dir.is_symlink():
            log.info(f"  HISAT2 index already linked: {target_index_dir}")
            return
        target_index_dir.parent.mkdir(parents=True, exist_ok=True)
        target_index_dir.symlink_to(source_index_dir)
        log.info(f"  Symlinked HISAT2 index: {source_index_dir} → {target_index_dir}")
        return

    log.warning("  No existing HISAT2 index found; include step 0a")


def run_step0_grid(base_config_path, grid=None, steps="0b-3a", collect_only=False):
    """Run fastp preprocessing parameter grid search."""
    grid = grid or STEP0_GRID
    param_names = list(grid.keys())
    combos = list(itertools.product(*grid.values()))
    n_total = len(combos)

    log.info("=" * 60)
    log.info(f"Step 0 Grid Search: {n_total} combinations")
    log.info(f"Parameters: {', '.join(param_names)}")
    log.info("=" * 60)

    results = []
    for i, combo_values in enumerate(combos, 1):
        combo_dict = dict(zip(param_names, combo_values))
        combo_name = _combo_to_dirname(combo_dict)
        run_name = f"grid_search/step0_preprocessing/{combo_name}"
        output_root = PIPELINE_ROOT / "results" / run_name

        log.info(f"\n[{i}/{n_total}] {combo_name}")

        if _is_combo_complete(output_root):
            log.info("  Already complete, collecting metrics only")
            metrics = _collect_metrics(output_root)
            metrics.update(combo_dict)
            metrics["combo_name"] = combo_name
            metrics["status"] = "ok"
            results.append(metrics)
            continue

        if collect_only:
            log.info("  Not complete, skipping (--collect-only)")
            continue

        try:
            # Modify config: remove gene_counts to force preprocessing, set unique shared dir
            with open(base_config_path) as f:
                cfg = yaml.safe_load(f)
            cfg = copy.deepcopy(cfg)
            cfg["run_name"] = run_name
            # Force preprocessing to use combo-specific shared directory
            cfg["shared_batch_name"] = run_name

            # Apply fastp params
            for key, value in combo_dict.items():
                _set_nested(cfg, key, value)

            # Force preprocessing by removing explicit gene_counts
            for batch_id, batch_cfg in cfg.get("batches", {}).items():
                if batch_cfg and "gene_counts" in batch_cfg:
                    del batch_cfg["gene_counts"]

            # Reuse existing HISAT2 index from the default batch to avoid rebuilding
            combo_shared_dir = PIPELINE_ROOT / "results" / "shared" / run_name
            _reuse_hisat2_index(cfg, combo_shared_dir)

            # Write temp config
            output_root.mkdir(parents=True, exist_ok=True)
            temp_config = output_root / "analysis_case_grid.yaml"
            with open(temp_config, "w") as f:
                yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)

            prepared = prepare_analysis_case(str(temp_config))
            result = run_pipeline(
                config_path=prepared["analysis_config"],
                steps=steps,
                dry_run=False,
                allow_hidden=steps.startswith("0"),
            )
            n_err = sum(
                1 for r in result.get("results", {}).values()
                if r["status"] == "error"
            )
            metrics = _collect_metrics(output_root)
            metrics.update(combo_dict)
            metrics["combo_name"] = combo_name
            metrics["status"] = "ok" if n_err == 0 else "error"
        except Exception as e:
            log.error(f"  FAILED: {e}")
            metrics = {**combo_dict, "combo_name": combo_name, "status": "error", "error": str(e)}

        results.append(metrics)

    if not results:
        log.warning("No results collected")
        return pd.DataFrame()

    df = pd.DataFrame(results)
    summary_dir = PIPELINE_ROOT / "results" / "grid_search" / "step0_preprocessing" / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(summary_dir / "grid_search_results.tsv", sep="\t", index=False)
    log.info(f"\nResults saved to {summary_dir / 'grid_search_results.tsv'}")

    ok_df = df[df["status"] == "ok"]
    if len(ok_df) > 0:
        _plot_step0_results(ok_df, summary_dir)

    return df


# ── Visualization ────────────────────────────────────────────────────────────

METRIC_LABELS = {
    "genes_after_filter": "Genes After Filter",
    "retention_rate": "Retention Rate (%)",
    "n_de_genes": "DE Genes",
    "within_corr_mean": "Within-Condition Corr",
    "between_corr_mean": "Between-Condition Corr",
    "pc1_pc2_cumulative": "PC1+PC2 Variance",
}


def _plot_step1_results(df, summary_dir):
    """Generate Step 1 grid search comparison plots."""
    log.info("Generating Step 1 plots...")

    # 1. Metric heatmaps: for each method, min_expr x min_samples
    methods = df["normalization.method"].unique()
    metrics_to_plot = ["genes_after_filter", "n_de_genes", "within_corr_mean",
                       "pc1_pc2_cumulative", "retention_rate", "between_corr_mean"]
    metrics_to_plot = [m for m in metrics_to_plot if m in df.columns]

    n_metrics = len(metrics_to_plot)
    n_methods = len(methods)
    if n_metrics > 0 and n_methods > 0:
        fig, axes = plt.subplots(n_methods, n_metrics,
                                 figsize=(3.5 * n_metrics, 3.5 * n_methods),
                                 squeeze=False)
        for row, method in enumerate(methods):
            mdf = df[df["normalization.method"] == method]
            for col, metric in enumerate(metrics_to_plot):
                ax = axes[row, col]
                if metric not in mdf.columns or mdf[metric].isna().all():
                    ax.set_visible(False)
                    continue
                pivot = mdf.pivot_table(
                    index="normalization.gene_filtering.min_expr",
                    columns="normalization.gene_filtering.min_samples",
                    values=metric,
                    aggfunc="first",
                )
                im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto")
                ax.set_xticks(range(len(pivot.columns)))
                ax.set_xticklabels(pivot.columns)
                ax.set_yticks(range(len(pivot.index)))
                ax.set_yticklabels(pivot.index)
                ax.set_xlabel("min_samples")
                ax.set_ylabel("min_expr")
                label = METRIC_LABELS.get(metric, metric)
                ax.set_title(f"{method.upper()}: {label}", fontsize=10)

                # Annotate cells
                for yi in range(pivot.shape[0]):
                    for xi in range(pivot.shape[1]):
                        val = pivot.values[yi, xi]
                        if pd.notna(val):
                            fmt = ".0f" if val > 10 else ".2f" if val < 1 else ".1f"
                            ax.text(xi, yi, f"{val:{fmt}}", ha="center", va="center",
                                    fontsize=8, color="white" if val > pivot.values.mean() else "black")
                fig.colorbar(im, ax=ax, shrink=0.8)

        fig.suptitle("Step 1 Grid Search: Gene Filtering Parameters", fontsize=14, y=1.02)
        fig.tight_layout()
        save_figure(fig, str(summary_dir / "metric_heatmaps.png"))
        log.info("  Saved metric_heatmaps.png")

    # 2. Pareto frontier: genes_retained vs n_de_genes
    if "genes_after_filter" in df.columns and "n_de_genes" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        for method in methods:
            mdf = df[df["normalization.method"] == method]
            color = "#4477AA" if method == "cpm" else "#EE6677"
            sizes = 50
            if "within_corr_mean" in mdf.columns:
                corr_vals = mdf["within_corr_mean"].fillna(0.9)
                sizes = ((corr_vals - corr_vals.min()) / (corr_vals.max() - corr_vals.min() + 1e-9) * 100 + 30)
            ax.scatter(mdf["genes_after_filter"], mdf["n_de_genes"],
                       s=sizes, c=color, alpha=0.7, edgecolors="white",
                       label=method.upper(), linewidths=0.5)
            # Label points
            for _, row in mdf.iterrows():
                expr = row["normalization.gene_filtering.min_expr"]
                samp = int(row["normalization.gene_filtering.min_samples"])
                ax.annotate(f"e{expr}/s{samp}", (row["genes_after_filter"], row["n_de_genes"]),
                            fontsize=6, alpha=0.7, textcoords="offset points", xytext=(3, 3))

        # Mark Pareto-optimal points
        pareto_mask = _pareto_mask(df[["genes_after_filter", "n_de_genes"]].values)
        pareto_df = df[pareto_mask]
        ax.scatter(pareto_df["genes_after_filter"], pareto_df["n_de_genes"],
                   s=150, facecolors="none", edgecolors="black", linewidths=2,
                   label="Pareto optimal", zorder=5)

        ax.set_xlabel("Genes After Filter")
        ax.set_ylabel("DE Genes")
        ax.set_title("Gene Retention vs DE Discovery Trade-off")
        ax.legend()
        style_axes(ax)
        fig.tight_layout()
        save_figure(fig, str(summary_dir / "pareto_frontier.png"))
        log.info("  Saved pareto_frontier.png")

    # 3. Parameter sensitivity
    _plot_sensitivity(df, summary_dir,
                      param_keys=["normalization.gene_filtering.min_expr",
                                  "normalization.gene_filtering.min_samples"],
                      group_key="normalization.method",
                      filename="parameter_sensitivity.png",
                      title="Step 1: Parameter Sensitivity")


def _plot_step0_results(df, summary_dir):
    """Generate Step 0 grid search comparison plots."""
    log.info("Generating Step 0 plots...")

    param_keys = [k for k in STEP0_GRID.keys() if k in df.columns]
    metrics_to_plot = ["genes_after_filter", "n_de_genes", "retention_rate",
                       "within_corr_mean", "pc1_pc2_cumulative"]
    metrics_to_plot = [m for m in metrics_to_plot if m in df.columns]

    if not param_keys or not metrics_to_plot:
        log.warning("  Insufficient data for Step 0 plots")
        return

    # Heatmap: unqualified_percent_limit x length_required
    n_metrics = len(metrics_to_plot)
    fig, axes = plt.subplots(1, n_metrics, figsize=(3.5 * n_metrics, 3.5), squeeze=False)
    for col, metric in enumerate(metrics_to_plot):
        ax = axes[0, col]
        if metric not in df.columns or df[metric].isna().all():
            ax.set_visible(False)
            continue

        upl_key = "preprocessing.fastp.unqualified_percent_limit"
        lr_key = "preprocessing.fastp.length_required"
        if upl_key in df.columns and lr_key in df.columns:
            pivot = df.pivot_table(index=upl_key, columns=lr_key, values=metric, aggfunc="first")
            im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto")
            ax.set_xticks(range(len(pivot.columns)))
            ax.set_xticklabels(pivot.columns)
            ax.set_yticks(range(len(pivot.index)))
            ax.set_yticklabels(pivot.index)
            ax.set_xlabel("length_required")
            ax.set_ylabel("unqualified_percent_limit")
            label = METRIC_LABELS.get(metric, metric)
            ax.set_title(label, fontsize=10)

            for yi in range(pivot.shape[0]):
                for xi in range(pivot.shape[1]):
                    val = pivot.values[yi, xi]
                    if pd.notna(val):
                        fmt = ".0f" if val > 10 else ".2f" if val < 1 else ".1f"
                        ax.text(xi, yi, f"{val:{fmt}}", ha="center", va="center",
                                fontsize=8, color="white" if val > pivot.values.mean() else "black")
            fig.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle("Step 0 Grid Search: Fastp Parameters", fontsize=14, y=1.02)
    fig.tight_layout()
    save_figure(fig, str(summary_dir / "fastp_sensitivity.png"))
    log.info("  Saved fastp_sensitivity.png")


def _plot_sensitivity(df, summary_dir, param_keys, group_key, filename, title):
    """Plot marginal sensitivity: how each metric changes with one parameter."""
    metrics_to_plot = ["genes_after_filter", "n_de_genes", "within_corr_mean", "pc1_pc2_cumulative"]
    metrics_to_plot = [m for m in metrics_to_plot if m in df.columns]

    if not metrics_to_plot or not param_keys:
        return

    n_params = len(param_keys)
    n_metrics = len(metrics_to_plot)
    fig, axes = plt.subplots(n_metrics, n_params,
                             figsize=(4 * n_params, 3 * n_metrics),
                             squeeze=False)

    groups = df[group_key].unique() if group_key in df.columns else [None]
    colors = ["#4477AA", "#EE6677", "#228833", "#CCBB44"]

    for col, param in enumerate(param_keys):
        if param not in df.columns:
            continue
        short_param = param.split(".")[-1]
        for row, metric in enumerate(metrics_to_plot):
            ax = axes[row, col]
            for gi, group in enumerate(groups):
                if group is not None:
                    gdf = df[df[group_key] == group]
                    label = str(group).upper()
                else:
                    gdf = df
                    label = None
                # Average over other parameters
                agg = gdf.groupby(param)[metric].mean()
                ax.plot(agg.index, agg.values, "o-", color=colors[gi % len(colors)],
                        label=label, markersize=5)
            ax.set_xlabel(short_param)
            label = METRIC_LABELS.get(metric, metric)
            ax.set_ylabel(label)
            if row == 0:
                ax.set_title(short_param)
            if col == 0 and groups[0] is not None:
                ax.legend(fontsize=8)
            style_axes(ax)

    fig.suptitle(title, fontsize=14, y=1.02)
    fig.tight_layout()
    save_figure(fig, str(summary_dir / filename))
    log.info(f"  Saved {filename}")


def _pareto_mask(points):
    """Find Pareto-optimal points (maximize both dimensions)."""
    n = len(points)
    is_pareto = np.ones(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # j dominates i if j >= i on both and strictly > on at least one
            if (points[j] >= points[i]).all() and (points[j] > points[i]).any():
                is_pareto[i] = False
                break
    return is_pareto


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Grid search for preprocessing and gene filtering parameters"
    )
    parser.add_argument(
        "--config", "-c", required=True,
        help="Path to base analysis_case.yaml",
    )
    parser.add_argument(
        "--mode", required=True, choices=["step0", "step1"],
        help="Which grid search to run: step0 (fastp) or step1 (gene filtering)",
    )
    parser.add_argument(
        "--steps", "-s", default=None,
        help="Pipeline steps to run per combo (default: 0b-3a for step0, 1a-3a for step1)",
    )
    parser.add_argument(
        "--collect-only", action="store_true",
        help="Only collect metrics from existing runs, do not run new combos",
    )
    args = parser.parse_args()

    if args.mode == "step1":
        steps = args.steps or "1a-3a"
        run_step1_grid(args.config, steps=steps, collect_only=args.collect_only)
    elif args.mode == "step0":
        steps = args.steps or "0b-3a"
        run_step0_grid(args.config, steps=steps, collect_only=args.collect_only)


if __name__ == "__main__":
    main()
