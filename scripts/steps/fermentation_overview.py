"""Step 4d: Fermentation Overview — condition1 vs condition2 fermentation metrics.

2x2 figure combining growth, product, substrate metabolism, and specific
substrate consumption rate.  Plus summary TSV with all metrics at all timepoints.

Data source: configured via advanced_analysis.fermentation_overview in analysis_case.yaml.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core.config_runtime import ensure_output_dir, get_timepoints, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, get_recipe_colors, plot_filename, COL2

log = setup_logger("advanced_analysis.fermentation_overview")

# Configure column mappings via analysis_case.yaml:
#   advanced_analysis.fermentation_overview.metrics
_DEFAULT_METRICS = {
    "biomass_dw": "Dry_Weight",
    "biomass_reactor_weight": "Reactor_Weight",
    "product": "Product",
    "product_ester": "Esterification_ratio",
    "substrate": "Substrate",
    "feed": "Feed_amount",
}

_DEFAULT_NUM_COLS = []


def run(cfg):
    adv_cfg = cfg.get("advanced_analysis", {}).get("fermentation_overview", {})
    if adv_cfg.get("enabled") is False:
        log.info("Fermentation overview disabled, skipping")
        return

    log.info("── Step 4d: Fermentation Overview ──")
    apply_style(cfg)

    out_dir = ensure_output_dir(cfg, "04_advanced_analysis/fermentation_overview")
    rc = get_recipe_colors(cfg)

    # Load configurable column mappings (fall back to built-in defaults)
    metrics = adv_cfg.get("metrics", _DEFAULT_METRICS)
    num_cols = adv_cfg.get("numeric_columns", _DEFAULT_NUM_COLS)
    reactor_id_col = adv_cfg.get("reactor_id_column", "Reactor_ID")
    time_col = adv_cfg.get("time_column", "Run_time (t)")
    reactor_ids = adv_cfg.get("reactor_ids", ["R01", "R02"])
    rna_timepoints = get_timepoints(cfg)

    # Load fermentation data
    pipeline_root = cfg["_pipeline_root"]
    excel_path = adv_cfg.get("data_file", "inputs/H2L-013_检测结果.xlsx")
    if not os.path.isabs(excel_path):
        excel_path = os.path.join(pipeline_root, excel_path)

    sheet = adv_cfg.get("sheet_name", "检测数据")
    df = pd.read_excel(excel_path, sheet_name=sheet)
    df.columns = [c.replace("\n", "").strip() for c in df.columns]

    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Derived: Biomass (g) = DW (mg/g) x Reactor_Weight (g) / 1000
    dw_col = metrics.get("biomass_dw", "Dry_Weight (mg/g)")
    rw_col = metrics.get("biomass_reactor_weight", "Reactor_Weight (Pre_Sampling,g)")
    if dw_col in df.columns and rw_col in df.columns:
        df["Biomass (g)"] = df[dw_col] * df[rw_col] / 1000

    # Split by reactor
    r1 = df[df[reactor_id_col] == reactor_ids[0]].set_index(time_col).sort_index()
    r2 = df[df[reactor_id_col] == reactor_ids[1]].set_index(time_col).sort_index()
    tps = sorted(set(r1.index) & set(r2.index))

    # N/D residual sugar -> 0 after first RNA timepoint
    rs_col = metrics.get("substrate", "Res_Sugar (RS, g/L)")
    first_tp = min(rna_timepoints) if rna_timepoints else 18
    for rx in [r1, r2]:
        if rs_col in rx.columns:
            mask = rx[rs_col].isna() & (rx.index >= first_tp)
            rx.loc[mask, rs_col] = 0.0

    log.info(f"  Loaded {len(tps)} shared timepoints from {os.path.basename(excel_path)}")

    # Save summary TSV
    _save_summary(r1, r2, tps, out_dir, rna_timepoints, metrics)

    # Plot 2x2 overview
    _plot_overview(r1, r2, tps, rc, out_dir, cfg, rna_timepoints, metrics, reactor_ids)

    write_readme(out_dir, "4d", "Fermentation Overview", {
        "fermentation_overview_summary.tsv": "Per-timepoint summary of R1 vs R2 metrics: OD, dry weight, biomass, residual sugar, sugar feed, AXT, esterification, and R2/R1 ratios",
        plot_filename(cfg, "fig_fermentation_overview.png"): "2x2 overview figure: biomass, astaxanthin + esterification, residual sugar + cumulative feed, and specific sugar consumption rate",
    })

    r1_label, r2_label = reactor_ids[0], reactor_ids[1]
    doc = (
        "### Parameters\n\n"
        f"- Data source: `{os.path.basename(excel_path)}` (sheet '{sheet}')\n"
        f"- **{len(tps)}** shared timepoints between {r1_label} and {r2_label}\n"
        f"- RNA sampling timepoints: {', '.join(str(t) + 'h' for t in rna_timepoints)}\n\n"
        "### Results\n\n"
        f"- 2x2 overview: biomass, product + esterification, substrate + feed, sSCR\n\n"
        f"{img(cfg, '04_advanced_analysis/fermentation_overview', 'fig_fermentation_overview.png', f'Fermentation overview {r1_label} vs {r2_label}')}\n"
    )
    update_section(cfg, "4d", "Fermentation Overview", doc)

    log.info("  Fermentation overview complete")


def _valid_tps(r1, r2, col, tps):
    return [t for t in tps
            if pd.notna(r1.loc[t, col]) and pd.notna(r2.loc[t, col])]


def _get_feed_data(rx, tps, feed_col="Sugar_Feed_amount (g)"):
    feed = rx[feed_col].copy()
    dt = pd.Series(rx.index, index=rx.index).diff()
    d_feed = feed.diff().fillna(0)
    rate = np.where(dt > 0, d_feed / dt, 0)
    return feed, rate, dt


def _save_summary(r1, r2, tps, out_dir, rna_timepoints, metrics):
    rows = []
    for t in tps:
        row = {"Time_h": t, "RNA_tp": "*" if t in rna_timepoints else ""}
        for pfx, rx in [("R1", r1), ("R2", r2)]:
            for col, key in [
                ("OD", "OD"),
                ("Dry_Weight (mg/g)", "DW"),
                ("Biomass (g)", "Biomass_g"),
                ("Res_Sugar (RS, g/L)", "RS"),
                ("Sugar_Feed_amount (g)", "Sugar_Feed"),
                ("Total_AXT (%)", "AXT"),
                ("Esterification_ratio (%)", "Ester"),
            ]:
                val = rx.loc[t, col] if (col in rx.columns and t in rx.index) else np.nan
                row[f"{pfx}_{key}"] = val if pd.notna(val) else np.nan
        # Ratios
        for m in ["OD", "DW", "Biomass_g", "RS", "AXT"]:
            v1, v2 = row.get(f"R1_{m}", np.nan), row.get(f"R2_{m}", np.nan)
            if pd.notna(v1) and pd.notna(v2) and v1 > 0:
                row[f"{m}_R2/R1"] = v2 / v1
        rows.append(row)

    summary = pd.DataFrame(rows)
    path = os.path.join(out_dir, "fermentation_overview_summary.tsv")
    summary.to_csv(path, sep="\t", index=False, float_format="%.4f")
    log.info(f"  Saved fermentation_overview_summary.tsv ({len(rows)} timepoints)")


def _plot_overview(r1, r2, tps, rc, out_dir, cfg, rna_timepoints, metrics, reactor_ids):
    """2x2 figure: Biomass, Product+Ester, Substrate+Feed, sSCR."""
    r1_label, r2_label = reactor_ids[0], reactor_ids[1]
    c_r1 = rc.get(r1_label, rc.get("condition1", "#2166AC"))
    c_r2 = rc.get(r2_label, rc.get("condition2", "#B2182B"))

    rs_col = metrics.get("substrate", "Res_Sugar (RS, g/L)")
    prod_col = metrics.get("product", "Total_AXT (%)")
    ester_col = metrics.get("product_ester", "Esterification_ratio (%)")
    feed_col = metrics.get("feed", "Sugar_Feed_amount (g)")
    rw_col = metrics.get("biomass_reactor_weight", "Reactor_Weight (Pre_Sampling,g)")

    rs_tps = _valid_tps(r1, r2, rs_col, tps)
    bm_tps = _valid_tps(r1, r2, "Biomass (g)", tps)
    prod_tps = _valid_tps(r1, r2, prod_col, tps)

    r1_feed, _, _ = _get_feed_data(r1, tps, feed_col)
    r2_feed, _, _ = _get_feed_data(r2, tps, feed_col)

    fig, axes = plt.subplots(2, 2, figsize=(COL2, 9))
    fig.suptitle(f"Fermentation Overview — {r1_label} vs {r2_label}",
                 fontsize=14, fontweight="bold", y=0.98)

    first_tp = min(rna_timepoints) if rna_timepoints else 18
    XTICKS = sorted(set([0] + rna_timepoints + [max(tps)]))

    def _add_rna_lines(ax):
        for t in rna_timepoints:
            ax.axvline(t, color="gray", alpha=0.15, ls="--", lw=0.8)

    # ── Top-left: Biomass ──
    ax = axes[0, 0]
    ax.plot(bm_tps, [r1.loc[t, "Biomass (g)"] for t in bm_tps],
            "o-", color=c_r1, ms=6, lw=1.5, label=r1_label)
    ax.plot(bm_tps, [r2.loc[t, "Biomass (g)"] for t in bm_tps],
            "s-", color=c_r2, ms=6, lw=1.5, label=r2_label)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Biomass (g)")
    ax.set_title("Biomass")
    ax.legend(fontsize=8)
    _add_rna_lines(ax)
    ax.set_xlim(-2, max(tps) + 5)
    ax.set_xticks(XTICKS)
    style_axes(ax)

    # ── Top-right: AXT + Esterification ──
    ax1 = axes[0, 1]
    ln1 = ax1.plot(prod_tps, [r1.loc[t, prod_col] for t in prod_tps],
                   "o-", color=c_r1, ms=6, lw=1.5, label=f"{r1_label} {prod_col}")
    ln2 = ax1.plot(prod_tps, [r2.loc[t, prod_col] for t in prod_tps],
                   "s-", color=c_r2, ms=6, lw=1.5, label=f"{r2_label} {prod_col}")
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel(prod_col)
    ax1.set_title("Product & Esterification")

    est_tps = _valid_tps(r1, r2, ester_col, tps)
    ax2 = ax1.twinx()
    ln3 = ax2.plot(est_tps,
                   [r1.loc[t, ester_col] for t in est_tps],
                   "o--", color=c_r1, ms=5, lw=1.0, alpha=0.5, label=f"{r1_label} Ester%")
    ln4 = ax2.plot(est_tps,
                   [r2.loc[t, ester_col] for t in est_tps],
                   "s--", color=c_r2, ms=5, lw=1.0, alpha=0.5, label=f"{r2_label} Ester%")
    ax2.set_ylabel("Esterification (%)", fontsize=9)
    ax2.axhline(50, color="gray", ls=":", lw=0.8, alpha=0.4)
    lns = ln1 + ln2 + ln3 + ln4
    ax1.legend(lns, [l.get_label() for l in lns], loc="upper left", fontsize=8)
    _add_rna_lines(ax1)
    ax1.set_xlim(-2, max(tps) + 5)
    ax1.set_xticks(XTICKS)
    style_axes(ax1)

    # ── Bottom-left: Residual Sugar + Cumulative Feed ──
    ax1 = axes[1, 0]
    ax1.plot(rs_tps, [r1.loc[t, rs_col] for t in rs_tps],
             "o-", color=c_r1, ms=6, lw=1.5, label=f"{r1_label} Substrate")
    ax1.plot(rs_tps, [r2.loc[t, rs_col] for t in rs_tps],
             "s-", color=c_r2, ms=6, lw=1.5, label=f"{r2_label} Substrate")
    ax1.axhline(0, color="gray", ls="-", alpha=0.2)
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel(rs_col)
    ax1.set_title("Substrate & Cumulative Feed")
    ax1.set_ylim(-2, 60)

    ax2 = ax1.twinx()
    ax2.plot(tps, [r1_feed.loc[t] for t in tps],
             "o--", color=c_r1, ms=4, lw=1.0, alpha=0.4, label=f"{r1_label} Feed")
    ax2.plot(tps, [r2_feed.loc[t] for t in tps],
             "s--", color=c_r2, ms=4, lw=1.0, alpha=0.4, label=f"{r2_label} Feed")
    ax2.set_ylabel(f"Cumul. {feed_col}", fontsize=9)

    lns1, lab1 = ax1.get_legend_handles_labels()
    lns2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(lns1 + lns2, lab1 + lab2, loc="upper left", fontsize=8)
    _add_rna_lines(ax1)
    ax1.set_xlim(-2, max(tps) + 5)
    ax1.set_xticks(XTICKS)
    style_axes(ax1)

    # ── Bottom-right: sSCR ──
    ax = axes[1, 1]
    for rx, rx_feed, color, marker, label in [
        (r1, r1_feed, c_r1, "o", r1_label),
        (r2, r2_feed, c_r2, "s", r2_label),
    ]:
        feed_vals = [rx_feed.loc[t] for t in tps]
        mid_t, sscr = [], []
        for i in range(1, len(tps)):
            t0, t1 = tps[i - 1], tps[i]
            dt = t1 - t0
            if dt <= 0 or t1 <= first_tp:
                continue
            bm0 = rx.loc[t0, "Biomass (g)"]
            bm1 = rx.loc[t1, "Biomass (g)"]
            if pd.notna(bm0) and pd.notna(bm1) and (bm0 + bm1) > 0:
                feed_rate = (feed_vals[i] - feed_vals[i - 1]) / dt
                rs0 = rx.loc[t0, rs_col]
                rs1 = rx.loc[t1, rs_col]
                rw0 = rx.loc[t0, rw_col]
                rw1 = rx.loc[t1, rw_col]
                rs0 = rs0 if pd.notna(rs0) else 0
                rs1 = rs1 if pd.notna(rs1) else 0
                d_total_rs = (rs1 * rw1 / 1000 - rs0 * rw0 / 1000) / dt
                consumption_rate = feed_rate - d_total_rs
                avg_bm = (bm0 + bm1) / 2
                mid_t.append((t0 + t1) / 2)
                sscr.append(consumption_rate / avg_bm)
        ax.plot(mid_t, sscr, f"{marker}-", color=color, ms=6, lw=1.5, label=label)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("sSCR (g substrate / h / g biomass)")
    ax.set_title("Specific Substrate Consumption Rate")
    ax.legend(fontsize=8)
    _add_rna_lines(ax)
    ax.set_xlim(-2, max(tps) + 5)
    ax.set_xticks(XTICKS)
    style_axes(ax)

    fig.tight_layout()
    fermentation_overview_plot = plot_filename(cfg, "fig_fermentation_overview.png")
    save_figure(fig, os.path.join(out_dir, "fig_fermentation_overview.png"), cfg)
    log.info(f"  Saved {fermentation_overview_plot}")
