"""Step 4b: Temporal Causality — multi-module lag-1 cross-correlation analysis.

Two-panel figure:
  Left:  Multi-module temporal dynamics (mean log2FC per timepoint)
  Right: Lag-1 cross-correlation (forward vs reverse) bar chart

Requires organism-specific gene modules defined in config under
advanced_analysis.temporal_causality.modules.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from ..core.config_data import get_filtered_expr, get_sample_names, require_sample_manifest
from ..core.config_runtime import get_conditions, get_timepoints, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, plot_filename

log = setup_logger("advanced_analysis.temporal_causality")

# Default colors for modules
_MODULE_COLORS = {
    "torc1": "#984EA3",
    "sterol": "#377EB8",
    "fa_degradation": "#FF7F00",
    "mva": "#E78AC3",
}

_MODULE_SPECS = [
    ("TORC1", "TORC1 net activity", "torc1"),
    ("Sterol", "Sterol pathway (12 ERG)", "sterol"),
    ("FA-deg", "FA degradation (19 genes)", "fa_degradation"),
    ("MVA", "MVA pathway (20 genes)", "mva"),
]


def _mean_fc(log2cpm, gene_list, r1_sample_groups, r2_sample_groups):
    """Compute mean log2FC (R2-R1) per timepoint for a gene set.

    Parameters
    ----------
    r1_sample_groups, r2_sample_groups : list of list[str]
        Per-timepoint sample ID lists, aligned by index.
    """
    ids = [g for g in gene_list if g in log2cpm.index]
    fcs = []
    for r1_ids, r2_ids in zip(r1_sample_groups, r2_sample_groups):
        r1_mean = log2cpm.loc[ids, r1_ids].mean(axis=1)
        r2_mean = log2cpm.loc[ids, r2_ids].mean(axis=1)
        fcs.append(float((r2_mean - r1_mean).mean()))
    return np.array(fcs)


def _lag1_r(x, y):
    """Lag-1 Pearson correlation: x[t] vs y[t+1]. Returns (r, p-value)."""
    r, p = pearsonr(x[:-1], y[1:])
    return r, p


def run(cfg):
    adv_cfg = cfg.get("advanced_analysis", {}).get("temporal_causality", {})
    if adv_cfg.get("enabled") is False:
        log.info("Temporal causality disabled, skipping")
        return

    log.info("── Step 4b: Temporal Causality ──")

    # Load gene modules from config
    modules_cfg = adv_cfg.get("modules", {})
    if not modules_cfg:
        log.warning("  No gene modules configured in temporal_causality.modules — skipping")
        log.info("  Configure gene lists to enable this step (see references/rnaseq-config-reference.md)")
        return

    apply_style(cfg)

    samples = get_sample_names(cfg)
    timepoints = get_timepoints(cfg)
    out_dir = ensure_output_dir(cfg, "04_advanced_analysis/temporal_causality")

    cond1, cond2 = get_conditions(cfg)[:2]

    # Load expression data
    log2cpm = _load_log2expr(cfg, samples)

    # Build timepoint-aligned sample groups
    r1_groups, r2_groups = _sample_groups_by_timepoint(cfg, cond1, cond2, timepoints)

    # Compute TORC1 net activity: activators_mean_FC - neg_regulators_mean_FC
    torc1_act = modules_cfg.get("torc1_activators", [])
    torc1_neg = modules_cfg.get("torc1_neg_regs", [])
    if torc1_act and torc1_neg:
        act_fc = _mean_fc(log2cpm, torc1_act, r1_groups, r2_groups)
        neg_fc = _mean_fc(log2cpm, torc1_neg, r1_groups, r2_groups)
        torc1 = act_fc - neg_fc
    else:
        torc1 = None

    # Compute module mean FCs
    sterol_genes = modules_cfg.get("sterol", [])
    sterol = _mean_fc(log2cpm, sterol_genes, r1_groups, r2_groups) if sterol_genes else None
    fa_deg_genes = modules_cfg.get("fa_degradation", [])
    fa_deg = _mean_fc(log2cpm, fa_deg_genes, r1_groups, r2_groups) if fa_deg_genes else None
    mva_genes = modules_cfg.get("mva", [])
    mva = _mean_fc(log2cpm, mva_genes, r1_groups, r2_groups) if mva_genes else None

    # Log computed modules
    named_modules = {}
    if torc1 is not None:
        named_modules["TORC1"] = (torc1, _MODULE_COLORS.get("torc1", "#984EA3"))
        log.info(f"  TORC1:  {' -> '.join(f'{v:+.2f}' for v in torc1)}")
    if sterol is not None:
        named_modules["Sterol"] = (sterol, _MODULE_COLORS.get("sterol", "#377EB8"))
        log.info(f"  Sterol: {' -> '.join(f'{v:+.2f}' for v in sterol)}")
    if fa_deg is not None:
        named_modules["FA-deg"] = (fa_deg, _MODULE_COLORS.get("fa_degradation", "#FF7F00"))
        log.info(f"  FA-deg: {' -> '.join(f'{v:+.2f}' for v in fa_deg)}")
    if mva is not None:
        named_modules["MVA"] = (mva, _MODULE_COLORS.get("mva", "#E78AC3"))
        log.info(f"  MVA:    {' -> '.join(f'{v:+.2f}' for v in mva)}")

    if len(named_modules) < 2:
        log.warning("  Need at least 2 modules for lag-1 analysis — skipping")
        return

    # Lag-1 correlations (forward and reverse) for all module pairs
    module_names = list(named_modules.keys())
    pairs = []
    for i, a_name in enumerate(module_names):
        for b_name in module_names[i + 1:]:
            a_vals, a_color = named_modules[a_name]
            b_vals, _ = named_modules[b_name]
            pairs.append((a_name, b_name, a_vals, b_vals, a_color))

    n_pairs = len(timepoints) - 1
    if n_pairs <= 4:
        log.warning(f"  Lag-1 correlations based on only {n_pairs} data pairs — "
                    "results are exploratory and should be interpreted with caution")

    lag_results = []
    for a_name, b_name, a_vals, b_vals, color in pairs:
        fwd_r, fwd_p = _lag1_r(a_vals, b_vals)
        rev_r, rev_p = _lag1_r(b_vals, a_vals)
        lag_results.append({
            "A": a_name, "B": b_name,
            "forward": fwd_r, "forward_pvalue": fwd_p,
            "reverse": rev_r, "reverse_pvalue": rev_p,
            "n_pairs": n_pairs,
            "color": color,
        })
        log.info(f"  {a_name}->{b_name}: fwd={fwd_r:+.2f} (p={fwd_p:.3f}), rev={rev_r:+.2f} (p={rev_p:.3f})")

    # Save data
    module_df = pd.DataFrame({"timepoint": timepoints})
    if torc1 is not None:
        module_df["TORC1_net_activity"] = torc1
    if sterol is not None:
        module_df["Sterol_12ERG"] = sterol
    if fa_deg is not None:
        module_df["FA_deg_no_ADH"] = fa_deg
    if mva is not None:
        module_df["MVA_pathway"] = mva
    module_df.to_csv(os.path.join(out_dir, "module_fc_timeseries.tsv"),
                     sep="\t", index=False)

    lag_df = pd.DataFrame(lag_results)
    lag_df.to_csv(os.path.join(out_dir, "lag1_correlations.tsv"),
                  sep="\t", index=False)

    # Plot
    _plot_two_panel(timepoints, named_modules, lag_results, out_dir, cfg)

    write_readme(out_dir, "4b", "Temporal Causality", {
        "module_fc_timeseries.tsv": "Mean log2FC per timepoint for TORC1 net activity, Sterol (12 ERG), FA degradation (19 genes), and MVA pathway (20 genes)",
        "lag1_correlations.tsv": "Lag-1 Pearson cross-correlations (forward and reverse) between all module pairs",
        plot_filename(cfg, "fig_temporal_causality.png"): "Two-panel figure: four-module temporal dynamics and lag-1 cross-correlation bar chart",
    })

    trajectory_lines = "\n".join(
        f"- {module_name} trajectory: {' -> '.join(f'{v:+.2f}' for v in values)}"
        for module_name, (values, _) in named_modules.items()
    )

    doc = (
        "### Method\n\n"
        "TORC1 net activity = mean FC of 7 activators - mean FC of 9 negative regulators. "
        "Lag-1 Pearson cross-correlation tests temporal precedence between "
        "TORC1, Sterol (12 ERG), FA degradation (19 genes), and MVA (20 genes).\n\n"
        "### Results\n\n"
        f"{trajectory_lines}\n"
        f"- {len(lag_results)} pairwise lag-1 correlations computed (forward vs reverse)\n\n"
        f"{img(cfg, '04_advanced_analysis/temporal_causality', 'fig_temporal_causality.png', 'Four-module temporal dynamics and lag-1 cross-correlation')}\n"
    )
    update_section(cfg, "4b", "Temporal Causality", doc)

    log.info("  Temporal causality complete")


def _plot_two_panel(timepoints, named_modules, lag_results, out_dir, cfg):
    """Two-panel figure: dynamics + lag-1 bar chart."""
    cond1, cond2 = get_conditions(cfg)[:2]
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("TORC1 / Sterol / FA Degradation / MVA — Temporal Causality",
                 fontsize=14, fontweight="bold", y=0.98)
    fig.subplots_adjust(wspace=0.30)

    # ── Left: Four-module time series ──
    ax = ax_left

    plotted_modules = []
    for module_name, label, color_key in _MODULE_SPECS:
        values_color = named_modules.get(module_name)
        if values_color is None:
            continue
        vals, color = values_color
        plotted_modules.append((module_name, vals, color))
        ax.plot(timepoints, vals, "o-", color=color, ms=6, lw=1.5, label=label)
        for i, t in enumerate(timepoints):
            ax.text(t, vals[i] + 0.03, f"{vals[i]:+.2f}", ha="center", va="bottom",
                    fontsize=8, color=color,
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

    # Lag arrows: TORC1 -> Sterol when both modules are present.
    torc1_values = named_modules.get("TORC1")
    sterol_values = named_modules.get("Sterol")
    if torc1_values is not None and sterol_values is not None:
        torc1, torc1_color = torc1_values
        sterol, _ = sterol_values
        for i in range(len(timepoints) - 1):
            ax.annotate("", xy=(timepoints[i+1], sterol[i+1]),
                        xytext=(timepoints[i], torc1[i]),
                        arrowprops=dict(arrowstyle="-|>", color=torc1_color, lw=0.8,
                                        alpha=0.2, connectionstyle="arc3,rad=0.25"))

    ax.axhline(0, ls="--", c="grey", lw=1.0)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(f"Mean log2FC ({cond2} - {cond1})")
    ax.set_title("Four-Module Temporal Dynamics", fontsize=11, fontweight="bold")
    ax.set_xticks(timepoints)
    ax.legend(loc="lower right", fontsize=8)
    style_axes(ax)

    # ── Right: Lag-1 bar chart ──
    ax = ax_right

    pair_names = [f"{r['A']}\n-> {r['B']}" for r in lag_results]
    pair_colors = [r["color"] for r in lag_results]
    fwd_vals = [r["forward"] for r in lag_results]
    rev_vals = [r["reverse"] for r in lag_results]

    x = np.arange(len(pair_names))
    w = 0.35

    bars_fwd = ax.bar(x - w/2, fwd_vals, w, color=pair_colors,
                       edgecolor="white", linewidth=0.4,
                       label="Forward: A[t] -> B[t+1]", alpha=0.85)
    bars_rev = ax.bar(x + w/2, rev_vals, w, color=pair_colors,
                       edgecolor="white", linewidth=0.4,
                       label="Reverse: B[t] -> A[t+1]", alpha=0.35, hatch="//")

    # Value labels
    for bar, val in zip(bars_fwd, fwd_vals):
        y = bar.get_height()
        va = "bottom" if y >= 0 else "top"
        offset = 0.03 if y >= 0 else -0.03
        ax.text(bar.get_x() + bar.get_width()/2, y + offset,
                f"{val:+.2f}", ha="center", va=va, fontsize=9, fontweight="bold")
    for bar, val in zip(bars_rev, rev_vals):
        y = bar.get_height()
        va = "bottom" if y >= 0 else "top"
        offset = 0.03 if y >= 0 else -0.03
        ax.text(bar.get_x() + bar.get_width()/2, y + offset,
                f"{val:+.2f}", ha="center", va=va, fontsize=8, color="#666")

    # Verdict annotations
    for i in range(len(lag_results)):
        fwd, rev = fwd_vals[i], rev_vals[i]
        if fwd > rev + 0.3:
            verdict, color = "A leads B", "#2E7D32"
        elif rev > fwd + 0.3:
            verdict, color = "B leads A", "#C62828"
        elif abs(fwd) < 0.3 and abs(rev) < 0.3:
            verdict, color = "No relation", "#C62828"
        else:
            verdict, color = "Unclear", "#666"
        y_max = max(abs(fwd), abs(rev)) + 0.15
        ax.text(x[i], y_max, verdict, ha="center", va="bottom", fontsize=9,
                fontweight="bold", color=color)

    ax.axhline(0, color="black", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(pair_names, fontsize=9)
    ax.set_ylabel("Lag-1 Pearson r")
    ax.set_title("Lag-1 Cross-Correlation (Forward vs Reverse)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(axis="y", alpha=0.2, ls="--")
    ax.set_ylim(-1.05, 1.25)
    style_axes(ax)

    fig.tight_layout()
    temporal_causality_plot = plot_filename(cfg, "fig_temporal_causality.png")
    save_figure(fig, os.path.join(out_dir, "fig_temporal_causality.png"), cfg)
    log.info(f"  Saved {temporal_causality_plot}")


def _sample_groups_by_timepoint(cfg, cond1, cond2, timepoints):
    """Return per-timepoint sample ID lists for each condition."""
    manifest = require_sample_manifest(
        cfg,
        required_columns={"sample_id", "target_condition", "timepoint", "batch_id"},
    )
    r1_groups, r2_groups = [], []
    for tp in timepoints:
        for cond, groups in [(cond1, r1_groups), (cond2, r2_groups)]:
            sub = manifest[(manifest["target_condition"] == cond) &
                           (manifest["timepoint"].astype(int) == int(tp))]
            ids = sub["sample_id"].tolist()
            if not ids:
                raise ValueError(f"No samples for {cond} at {tp}h in temporal causality")
            groups.append(ids)
    return r1_groups, r2_groups


def _load_log2expr(cfg, samples):
    return get_filtered_expr(cfg, log2=True)[samples["all"]]
