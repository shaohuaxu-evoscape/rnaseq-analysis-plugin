"""Step 4a: Transporter & Aminopeptidase Analysis.

Comprehensive transporter census by substrate type, nitrogen transporter
analysis, and aminopeptidase differential analysis.

Generates:
  - fig_transporter.png      — N-transporter category trends + key gene CPM
  - fig_aminopeptidase.png   — aminopeptidase vs endopeptidase + DE bar
  - transporter_census.tsv   — per-gene stats
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ..core.config_data import get_filtered_expr, get_sample_names
from ..core.config_runtime import get_conditions, get_path, get_timepoints, ensure_output_dir, write_readme
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, get_recipe_colors, plot_filename, COL2

log = setup_logger("advanced_analysis.transporter_aminopeptidase")

# ── Constants ────────────────────────────────────────────────────────────────

C_UP = "#D73027"
C_DN = "#4575B4"
C_NS = "#999999"

# ── Organism-specific annotations ──
# Configure via analysis_case.yaml advanced_analysis.transporter_aminopeptidase.*
# Each list/dict below is empty by default; provide your organism's data in config.

# Aminopeptidase info: [{gene, name, substrate, family, de, dir}]
_DEFAULT_AMINO_INFO = []

# Transporter manual overrides: {gene_id: substrate_category}
_DEFAULT_TRANSPORTER_OVERRIDES = {}

_DEFAULT_SUBSTRATE_CATEGORIES = {
    "Peptide/oligopeptide": r"peptide transport|oligopeptide transport|oligopeptide transmembrane|dipeptide transmembrane",
    "Amino acid": r"amino acid transmembrane transporter|amino acid transport",
    "Ammonium": r"ammonium transmembrane transporter|ammonium transport|ammonium channel|methylammonium",
    "Urea": r"urea transmembrane transporter|urea transport|urea channel",
    "Sugar/hexose": r"sugar transmembrane transporter|hexose transmembrane|glucose transmembrane|sugar transport|monosaccharide transport",
    "Lipid/sterol": r"lipid transport|sterol transport|fatty acid transport|phospholipid transport",
    "Ion (metal)": r"iron ion transmembrane|zinc ion transmembrane|copper ion transmembrane|metal ion transmembrane|manganese.*transport|iron.*transport|siderophore",
    "Ion (other)": r"potassium ion transmembrane|sodium ion transmembrane|calcium ion transmembrane|chloride transmembrane|phosphate transmembrane|sulfate transmembrane|inorganic phosphate",
    "Nucleotide/nucleoside": r"nucleotide transmembrane|nucleoside transmembrane|purine.*transport",
    "Vitamin/cofactor": r"vitamin transmembrane|thiamine transmembrane|biotin transmembrane|folate transmembrane|riboflavin transmembrane",
    "Drug/xenobiotic": r"drug transmembrane|multidrug|xenobiotic transmembrane",
    "ABC transporter": r"ABC-type|ATPase-coupled transmembrane|ABC transporter",
    "Vesicle/intracellular": r"vesicle-mediated transport|intracellular protein transport|endosomal transport|vacuolar transport|Golgi.*transport",
    "Mitochondrial": r"mitochondrial transport|mitochondrial transmembrane",
    "ER/secretory": r"endoplasmic reticulum.*transport|ER to Golgi|protein insertion into ER",
}


def _load_annotations(ta_cfg):
    """Load organism-specific annotations from config or defaults."""
    return {
        "amino_info": ta_cfg.get("amino_info", _DEFAULT_AMINO_INFO),
        "transporter_overrides": ta_cfg.get("transporter_overrides", _DEFAULT_TRANSPORTER_OVERRIDES),
        "substrate_categories": ta_cfg.get("substrate_categories", _DEFAULT_SUBSTRATE_CATEGORIES),
        "endo_de_info": ta_cfg.get("endo_de_info", _DEFAULT_ENDO_DE_INFO),
    }


# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════

def _load_go_annotations(cfg):
    """Load GO annotations and UniProt gene name mapping from pipeline results."""
    out_root = get_path(cfg, "output_dir")
    cache_dir = os.path.join(out_root, "03_differential_analysis", "go_enrichment", "cache")
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")

    go = pd.read_csv(os.path.join(cache_dir, "uniprot_go.tsv"), sep="\t")
    go["all_go"] = (
        go["Gene Ontology (biological process)"].fillna("") + ";" +
        go["Gene Ontology (molecular function)"].fillna("")
    )

    uni = pd.read_csv(os.path.join(cache_dir, "uniprot_genes.tsv"), sep="\t")

    prefix_pattern = re.compile(rf"({re.escape(gene_id_prefix)}\w+)") if gene_id_prefix else None

    def extract_gene(names):
        if pd.isna(names):
            return None
        if prefix_pattern:
            m = prefix_pattern.search(str(names))
            return m.group(1) if m else None
        return None

    uni["gene"] = uni["Gene Names"].apply(extract_gene)
    entry2gene = uni.dropna(subset=["gene"]).set_index("Entry")["gene"].to_dict()
    return go, entry2gene


def _load_de_results(cfg):
    """Load DE results table from pipeline results."""
    out_root = get_path(cfg, "output_dir")
    de_file = os.path.join(out_root, "03_differential_analysis", "de_screening", "de_results_all.tsv")
    de = pd.read_csv(de_file, sep="\t")
    # Normalize column names to match expected format
    col_map = {
        "mean_log2fc": "mean_log2FC", "fdr": "FDR", "is_de": "is_DE",
    }
    for t in get_timepoints(cfg):
        col_map[f"log2fc_{t}h"] = f"log2FC_{t}h"
    de.rename(columns={k: v for k, v in col_map.items() if k in de.columns},
              inplace=True)
    cond2 = get_conditions(cfg)[1]
    up_label = f"up_in_{cond2}"
    down_label = f"down_in_{cond2}"
    dir_map = {"up": up_label, "down": down_label, "ns": "ns"}
    de["direction"] = de["direction"].map(dir_map).fillna("ns")
    return de


def _load_de_indexed(cfg):
    """Load DE results indexed by gene."""
    de = _load_de_results(cfg)
    return de.set_index("gene")


def _load_log2expr(cfg):
    """Load log2-normalized expression matrix from pipeline results."""
    return get_filtered_expr(cfg, log2=True)


# ══════════════════════════════════════════════════════════════════════════════
# Transporter classification
# ══════════════════════════════════════════════════════════════════════════════

def _find_genes_by_go(go_df, entry2gene, pattern):
    mask = go_df["all_go"].str.contains(pattern, case=False, na=False)
    entries = go_df[mask]["Entry"].unique()
    return sorted(set(entry2gene[e] for e in entries if e in entry2gene))


def _classify_transporters(go_df, entry2gene, transporter_overrides, substrate_categories):
    all_transport = _find_genes_by_go(
        go_df, entry2gene,
        r"transmembrane transporter|transporter activity|transport|channel activity"
    )
    for g in transporter_overrides:
        if g not in all_transport:
            all_transport.append(g)
    all_transport = sorted(set(all_transport))

    gene_categories = {}
    for gene in all_transport:
        yali_entries = [e for e, y in entry2gene.items() if y == gene]
        if not yali_entries:
            continue
        gene_go = go_df[go_df["Entry"].isin(yali_entries)]["all_go"].str.cat(sep=";")

        assigned = []
        for cat_name, pattern in substrate_categories.items():
            if re.search(pattern, gene_go, re.IGNORECASE):
                assigned.append(cat_name)

        if not assigned:
            if re.search(r"transmembrane transporter activity", gene_go, re.I):
                assigned = ["Other transmembrane"]
            else:
                assigned = ["Other transport"]

        if gene in transporter_overrides:
            assigned = [transporter_overrides[gene]]

        gene_categories[gene] = assigned

    return gene_categories


def _compute_category_stats(gene_categories, de_all, fc_cols):
    gene_primary = {}
    for gene, cats in gene_categories.items():
        for cat in cats:
            if cat not in ("Other transmembrane", "Other transport",
                           "Vesicle/intracellular", "ER/secretory"):
                gene_primary[gene] = cat
                break
        else:
            gene_primary[gene] = cats[0]

    rows = []
    for gene, cat in gene_primary.items():
        if gene in de_all.index:
            r = de_all.loc[gene]
            rows.append({
                "gene": gene, "category": cat,
                "mean_log2FC": r["mean_log2FC"],
                "is_DE": r["is_DE"],
                "direction": r["direction"] if r["is_DE"] else "ns",
                "FDR": r["FDR"],
                **{col: r[col] for col in fc_cols},
            })
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1: N-source Transporters (1×3)
# ══════════════════════════════════════════════════════════════════════════════

def _plot_transporter(stats_df, log2cpm, out_dir, cfg, timepoints, fc_cols):
    """N-source transporter figure: category trends + peptide DE + AA DE."""
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")
    samples = get_sample_names(cfg)
    cond1, cond2 = get_conditions(cfg)[:2]
    up_label = f"up_in_{cond2}"
    down_label = f"down_in_{cond2}"

    n_cats = ["Peptide/oligopeptide", "Amino acid", "Ammonium", "Urea"]
    ndf = stats_df[stats_df["category"].isin(n_cats)].copy()

    cat_colors = {
        "Peptide/oligopeptide": "#E41A1C", "Amino acid": "#377EB8",
        "Ammonium": "#4DAF4A", "Urea": "#FF7F00",
    }

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(COL2 + 2.5, 5.5))

    # ── Panel A: N-transporter mean FC by category over time ──
    for cat in n_cats:
        sub = ndf[ndf["category"] == cat]
        if len(sub) == 0:
            continue
        mean_fc = sub[fc_cols].mean().values
        n = len(sub)
        n_de = sub["is_DE"].sum()
        ax_a.plot(timepoints, mean_fc, "o-", color=cat_colors[cat], ms=6, lw=1.5,
                  label=f"{cat} (n={n}, {n_de} DE)")
        if n > 1:
            sem = sub[fc_cols].sem().values
            ax_a.fill_between(timepoints, mean_fc - sem, mean_fc + sem,
                              color=cat_colors[cat], alpha=0.10)
    ax_a.axhline(0, color="black", ls="--", lw=1.0, alpha=0.4)
    ax_a.set_xlabel("Time (h)")
    ax_a.set_ylabel(f"Mean log2FC ({cond2}/{cond1})")
    ax_a.set_title("N-Transporter FC by Category")
    ax_a.legend(fontsize=8, loc="lower left")
    ax_a.set_xticks(timepoints)
    style_axes(ax_a)

    # ── Panel B: DE peptide transporter log2FC temporal ──
    # Configurable via analysis_case.yaml advanced_analysis.transporter_aminopeptidase.de_peptide_genes
    ta_cfg = cfg.get("advanced_analysis", {}).get("transporter_aminopeptidase", {})
    de_peptide_cfg = ta_cfg.get("de_peptide_genes", [])
    de_peptide = [(g["gene"], g["name"], g.get("fmt", "-"), g.get("color", "#E41A1C"))
                  for g in de_peptide_cfg]
    for gene, name, fmt, color in de_peptide:
        row = stats_df[stats_df["gene"] == gene]
        if row.empty:
            continue
        r = row.iloc[0]
        fcs = [r[c] for c in fc_cols]
        fc_mean = r["mean_log2FC"]
        ax_b.plot(timepoints, fcs, fmt, color=color, marker="o", ms=6, lw=1.5,
                  label=f"{name} ({fc_mean:+.1f})")
    ax_b.axhline(0, color="black", ls="--", lw=1.0, alpha=0.5)
    ax_b.set_xlabel("Time (h)")
    ax_b.set_ylabel(f"log2FC ({cond2}/{cond1})")
    ax_b.set_title("DE Peptide Transporter FC")
    ax_b.legend(fontsize=8, loc="best")

    ax_b.set_xticks(timepoints)
    style_axes(ax_b)

    # ── Panel C: DE amino acid transporter log2FC temporal ──
    aa_de = stats_df[(stats_df["category"] == "Amino acid") & (stats_df["is_DE"])].copy()
    aa_up = aa_de[aa_de["direction"] == up_label].sort_values("mean_log2FC", ascending=False)
    aa_dn = aa_de[aa_de["direction"] == down_label].sort_values("mean_log2FC", ascending=True)

    # Gene name mapping (configurable via config)
    AA_NAMES = ta_cfg.get("aa_gene_names", {})

    # Use a red gradient for up, blue gradient for down
    up_colors = ["#D32F2F", "#E57373", "#FFAB91"]
    dn_colors = ["#1565C0", "#42A5F5", "#64B5F6", "#90CAF9", "#BBDEFB"]

    for i, (_, r) in enumerate(aa_up.iterrows()):
        fcs = [r[c] for c in fc_cols]
        short = r["gene"].replace(gene_id_prefix, "")
        name = AA_NAMES.get(r["gene"], "")
        c = up_colors[i] if i < len(up_colors) else up_colors[-1]
        ax_c.plot(timepoints, fcs, "-", color=c, marker="o", ms=6, lw=1.5,
                  label=f"{short} {name} ({r['mean_log2FC']:+.1f})")

    for i, (_, r) in enumerate(aa_dn.iterrows()):
        fcs = [r[c] for c in fc_cols]
        short = r["gene"].replace(gene_id_prefix, "")
        name = AA_NAMES.get(r["gene"], "")
        c = dn_colors[i] if i < len(dn_colors) else dn_colors[-1]
        ax_c.plot(timepoints, fcs, "-", color=c, marker="o", ms=6, lw=1.5,
                  label=f"{short} {name} ({r['mean_log2FC']:+.1f})")

    ax_c.axhline(0, color="black", ls="--", lw=1.0, alpha=0.5)
    ax_c.set_xlabel("Time (h)")
    ax_c.set_ylabel(f"log2FC ({cond2}/{cond1})")
    ax_c.set_title(f"DE Amino Acid Transporter FC ({len(aa_up)}\u2191 {len(aa_dn)}\u2193)")
    ax_c.legend(fontsize=7, loc="best", ncol=2)
    ax_c.set_xticks(timepoints)
    style_axes(ax_c)

    fig.tight_layout()
    transporter_plot = plot_filename(cfg, "fig_transporter.png")
    save_figure(fig, os.path.join(out_dir, "fig_transporter.png"), cfg)
    log.info(f"  Saved {transporter_plot}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2: Aminopeptidase / Endopeptidase (1×3)
# ══════════════════════════════════════════════════════════════════════════════

# DE endopeptidase annotations: {gene_id: (name, family)}
_DEFAULT_ENDO_DE_INFO = {}


def _plot_aminopeptidase(de_all_flat, go_df, entry2gene, out_dir, cfg, timepoints, fc_cols, amino_info, endo_de_info):
    """N-metabolism enzyme figure: trend + amino DE bar + endo DE bar."""
    gene_id_prefix = cfg["project"].get("gene_id_prefix", "")
    cond1, cond2 = get_conditions(cfg)[:2]
    up_label = f"up_in_{cond2}"
    down_label = f"down_in_{cond2}"
    amino_genes = _find_genes_by_go(go_df, entry2gene, r"aminopeptidase")
    endo_genes = _find_genes_by_go(go_df, entry2gene, r"endopeptidase")
    amino_data = de_all_flat[de_all_flat["gene"].isin(amino_genes)]
    endo_data = de_all_flat[de_all_flat["gene"].isin(endo_genes)]

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(COL2 + 2.5, 5.5))

    # ── Panel A: Aminopeptidase vs endopeptidase mean FC ──
    amino_mean = [amino_data[f"log2FC_{t}h"].mean() for t in timepoints]
    endo_mean = [endo_data[f"log2FC_{t}h"].mean() for t in timepoints]

    ax_a.plot(timepoints, amino_mean, "o-", color="#FF9800", ms=6, lw=1.5,
              label=f"Aminopeptidase (n={len(amino_data)})")
    ax_a.plot(timepoints, endo_mean, "o-", color="#7B1FA2", ms=6, lw=1.5,
              label=f"Endopeptidase (n={len(endo_data)})")
    ax_a.axhline(0, color="black", ls="--", lw=1.0, alpha=0.4)
    ax_a.fill_between(timepoints, amino_mean, endo_mean, alpha=0.08, color="#FF9800")
    ax_a.set_xlabel("Time (h)")
    ax_a.set_ylabel(f"Mean log2FC ({cond2}/{cond1})")
    ax_a.set_title("Aminopeptidase vs Endopeptidase FC Trend")
    ax_a.legend(fontsize=8, loc="lower left")
    ax_a.set_xticks(timepoints)
    style_axes(ax_a)

    # ── Panel B: DE aminopeptidase bar by enzyme family ──
    amino_info = [dict(g) for g in amino_info]
    for g in amino_info:
        row = de_all_flat[de_all_flat["gene"] == g["gene"]]
        if not row.empty:
            g["mean_fc"] = row.iloc[0]["mean_log2FC"]
        else:
            g["mean_fc"] = 0.0

    de_amino = sorted([g for g in amino_info if g["de"]],
                      key=lambda x: x["mean_fc"], reverse=True)
    y = np.arange(len(de_amino))
    family_colors = {"Metallo": "#FF9800", "Cysteine": "#9C27B0",
                     "Serine": "#009688", "Aspartyl": "#E53935"}
    bar_colors = [family_colors.get(g["family"], "#999") for g in de_amino]

    ax_b.barh(y, [g["mean_fc"] for g in de_amino], color=bar_colors,
              edgecolor="white", lw=0.5, height=0.55)
    ax_b.axvline(0, color="black", lw=1.0)
    ax_b.set_yticks(y)
    ylabels = [f"{g['gene'].replace(gene_id_prefix, '')}  {g['name']}" for g in de_amino]
    ax_b.set_yticklabels(ylabels, fontsize=8)
    ax_b.invert_yaxis()
    ax_b.set_xlabel(f"Mean log2FC ({cond2}/{cond1})")
    ax_b.set_title(f"DE Aminopeptidases ({len(de_amino)})")
    for i, g in enumerate(de_amino):
        color = C_UP if g["dir"] == "up" else C_DN
        ax_b.get_yticklabels()[i].set_color(color)
    legend_elements = [Patch(facecolor="#FF9800", label="Metallo"),
                       Patch(facecolor="#9C27B0", label="Cysteine"),
                       Patch(facecolor="#009688", label="Serine")]
    ax_b.legend(handles=legend_elements, fontsize=8, loc="lower right")
    style_axes(ax_b)

    # ── Panel C: DE endopeptidase bar by enzyme family ──
    endo_de = endo_data[endo_data["is_DE"]].copy()
    endo_de = endo_de.sort_values("mean_log2FC", ascending=False)

    y_e = np.arange(len(endo_de))
    endo_bar_colors = []
    endo_labels = []
    for _, r in endo_de.iterrows():
        gene = r["gene"]
        name, fam = endo_de_info.get(gene, ("Unknown", "Other"))
        endo_bar_colors.append(family_colors.get(fam, "#999"))
        short = gene.replace(gene_id_prefix, "")
        endo_labels.append(f"{short}  {name}")

    ax_c.barh(y_e, endo_de["mean_log2FC"].values, color=endo_bar_colors,
              edgecolor="white", lw=0.5, height=0.55)
    ax_c.axvline(0, color="black", lw=1.0)
    ax_c.set_yticks(y_e)
    ax_c.set_yticklabels(endo_labels, fontsize=8)
    ax_c.invert_yaxis()
    ax_c.set_xlabel(f"Mean log2FC ({cond2}/{cond1})")
    n_up = (endo_de["direction"] == up_label).sum()
    n_dn = (endo_de["direction"] == down_label).sum()
    ax_c.set_title(f"DE Endopeptidases ({n_up}\u2191 {n_dn}\u2193)")
    for i, (_, r) in enumerate(endo_de.iterrows()):
        color = C_UP if r["direction"] == up_label else C_DN
        ax_c.get_yticklabels()[i].set_color(color)
    legend_elements_e = [Patch(facecolor="#FF9800", label="Metallo"),
                         Patch(facecolor="#E53935", label="Aspartyl"),
                         Patch(facecolor="#009688", label="Serine"),
                         Patch(facecolor="#9C27B0", label="Cysteine")]
    ax_c.legend(handles=legend_elements_e, fontsize=8, loc="lower right")
    style_axes(ax_c)

    fig.tight_layout()
    aminopeptidase_plot = plot_filename(cfg, "fig_aminopeptidase.png")
    save_figure(fig, os.path.join(out_dir, "fig_aminopeptidase.png"), cfg)
    log.info(f"  Saved {aminopeptidase_plot}")


# ══════════════════════════════════════════════════════════════════════════════
# Main entry point
# ══════════════════════════════════════════════════════════════════════════════

def run(cfg):
    ta_cfg = cfg.get("advanced_analysis", {}).get("transporter_aminopeptidase", {})
    if not ta_cfg.get("enabled", True):
        log.info("Transporter/aminopeptidase analysis disabled, skipping")
        return

    log.info("── Step 4a: Transporter & Aminopeptidase Analysis ──")
    apply_style(cfg)

    out_dir = ensure_output_dir(cfg, "04_advanced_analysis/aminopeptidase")

    # Load organism-specific annotations from config or defaults
    annot = _load_annotations(ta_cfg)

    # Load data
    log.info("  Loading GO annotations...")
    go_df, entry2gene = _load_go_annotations(cfg)

    log.info("  Loading DE results...")
    de_all_flat = _load_de_results(cfg)
    de_all_idx = de_all_flat.set_index("gene")

    log.info("  Loading log2CPM...")
    log2cpm = _load_log2expr(cfg)

    # Classify transporters
    log.info("  Classifying transporters...")
    gene_categories = _classify_transporters(
        go_df, entry2gene, annot["transporter_overrides"], annot["substrate_categories"])
    log.info(f"  Total transporter genes: {len(gene_categories)}")

    # Derive timepoints and FC column names from config
    timepoints = get_timepoints(cfg)
    fc_cols = [f"log2FC_{t}h" for t in timepoints]

    stats_df = _compute_category_stats(gene_categories, de_all_idx, fc_cols)
    log.info(f"  With expression data: {len(stats_df)}, DE: {stats_df['is_DE'].sum()}")

    # Save census TSV
    tsv_path = os.path.join(out_dir, "transporter_census.tsv")
    stats_df.to_csv(tsv_path, sep="\t", index=False)
    log.info(f"  Saved {tsv_path}")

    # Generate figures
    log.info("  Generating transporter figure...")
    _plot_transporter(stats_df, log2cpm, out_dir, cfg, timepoints, fc_cols)

    log.info("  Generating aminopeptidase figure...")
    _plot_aminopeptidase(de_all_flat, go_df, entry2gene, out_dir, cfg, timepoints, fc_cols,
                         annot["amino_info"], annot["endo_de_info"])

    write_readme(out_dir, "4a", "Transporter & Aminopeptidase", {
        "transporter_census.tsv": "Per-gene transporter stats with substrate category, mean log2FC, DE status, and per-timepoint fold changes",
        plot_filename(cfg, "fig_transporter.png"): "N-transporter category FC trends, DE peptide transporter temporal FC, and DE amino acid transporter FC",
        plot_filename(cfg, "fig_aminopeptidase.png"): "Aminopeptidase vs endopeptidase mean FC trend, DE aminopeptidase bar by enzyme family, and DE endopeptidase bar",
    })

    n_de = stats_df["is_DE"].sum()
    n_cats = stats_df["category"].nunique()
    doc = (
        "### Method\n\n"
        "Transporter census by substrate category from GO annotations, "
        "with manual overrides for nitrogen transporters. "
        "Aminopeptidase/endopeptidase analysis using UniProt family classification.\n\n"
        "### Results\n\n"
        f"- **{len(stats_df)}** transporter genes classified into **{n_cats}** substrate categories\n"
        f"- **{n_de}** DE transporters (FDR < 0.05, |mean log2FC| > 0.585)\n"
        f"- **{len(annot['amino_info'])}** aminopeptidases curated ({sum(1 for a in annot['amino_info'] if a['de'])} DE)\n\n"
        f"{img(cfg, '04_advanced_analysis/aminopeptidase', 'fig_transporter.png', 'N-source transporter FC trends and key gene profiles')}\n\n"
        f"{img(cfg, '04_advanced_analysis/aminopeptidase', 'fig_aminopeptidase.png', 'Aminopeptidase vs endopeptidase analysis')}\n"
    )
    update_section(cfg, "4a", "Transporter & Aminopeptidase", doc)

    log.info("  Transporter & aminopeptidase analysis complete")
