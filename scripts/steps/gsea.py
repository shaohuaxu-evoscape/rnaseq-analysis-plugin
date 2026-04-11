"""Step 3d: GSEA — Gene Set Enrichment Analysis using gseapy.

Uses custom KEGG gene sets mapped to organism-specific gene IDs.
Falls back to saving ranking file for external analysis.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core.config_runtime import ensure_output_dir, get_path, write_readme
from ..core.de_helpers import load_de_results
from ..core.docgen import update_section, img
from ..core.logger import setup_logger
from ..core.plotting import apply_style, save_figure

log = setup_logger("differential_analysis.gsea")


def run(cfg):
    """Run GSEA using gseapy prerank with custom organism-specific gene sets."""
    gsea_cfg = cfg["differential_analysis"]["gsea"]
    if not gsea_cfg["enabled"]:
        log.info("GSEA disabled, skipping")
        return

    log.info("── Step 3d: GSEA ──")
    apply_style(cfg)

    out_dir = ensure_output_dir(cfg, "03_differential_analysis/gsea")
    # Load DE results for ranking
    de_results = load_de_results(cfg)
    if de_results is None:
        log.warning("  No DE results found, skipping GSEA")
        return

    # Build ranked gene list
    metric = gsea_cfg["ranking_metric"]
    if metric == "mean_log2fc":
        ranking = de_results["mean_log2fc"].sort_values(ascending=False)
    elif metric == "signed_neg_log10p":
        sign = np.sign(de_results["mean_log2fc"])
        ranking = (sign * (-np.log10(de_results["pvalue"].clip(1e-300)))).sort_values(ascending=False)
    else:
        ranking = de_results["mean_log2fc"].sort_values(ascending=False)

    ranking.to_csv(os.path.join(out_dir, "gsea_ranking.tsv"), sep="\t", header=True)
    log.info(f"  Ranked {len(ranking)} genes by {metric}")

    # Load custom KEGG gene sets from KEGG enrichment cache
    gene_sets = _load_kegg_gene_sets(cfg, cfg["project"].get("gene_id_prefix", ""))
    if not gene_sets:
        log.warning("  No KEGG gene sets available. Run KEGG enrichment (3c) first.")
        log.info("  Ranking file saved for external GSEA analysis")
        return

    log.info(f"  Loaded {len(gene_sets)} KEGG gene sets")

    # Try gseapy prerank
    try:
        import gseapy as gp

        pre_res = gp.prerank(
            rnk=ranking,
            gene_sets=gene_sets,
            min_size=gsea_cfg.get("min_size", 10),
            max_size=gsea_cfg.get("max_size", 500),
            permutation_num=gsea_cfg.get("permutations", 1000),
            outdir=os.path.join(out_dir, "gseapy_output"),
            seed=42,
            verbose=False,
        )

        if pre_res.res2d is not None and not pre_res.res2d.empty:
            pre_res.res2d.to_csv(
                os.path.join(out_dir, "gsea_results.tsv"), sep="\t", index=False
            )
            n_sig = (pre_res.res2d["FDR q-val"] < 0.25).sum() if "FDR q-val" in pre_res.res2d.columns else 0
            log.info(f"  GSEA results: {len(pre_res.res2d)} gene sets, {n_sig} sig (FDR<0.25)")
        else:
            log.info("  GSEA returned no results")

    except ImportError:
        log.warning("  gseapy not installed — saving ranking only")
        log.info("  Install with: pip install gseapy")

    except Exception as e:
        log.warning(f"  GSEA failed: {e}")
        log.info("  Ranking file saved for external GSEA analysis")

    write_readme(out_dir, "3d", "GSEA", {
        "gsea_ranking.tsv": "Ranked gene list used as GSEA input (sorted by ranking metric)",
        "gsea_results.tsv": "GSEA results with NES, p-values, and FDR q-values per gene set",
        "gseapy_output/": "Directory of gseapy prerank output files including enrichment plots",
    })

    doc = (
        f"### Parameters\n"
        f"- Ranking metric: {metric}\n"
        f"- Gene sets: KEGG pathways ({len(gene_sets)} sets)\n"
        f"- Min size: {gsea_cfg.get('min_size', 10)}\n"
        f"- Max size: {gsea_cfg.get('max_size', 500)}\n"
        f"- Permutations: {gsea_cfg.get('permutations', 1000)}\n\n"
        f"### Results\n"
        f"- Genes ranked: {len(ranking)}\n"
    )
    update_section(cfg, "3d", "GSEA", doc)

    log.info("  GSEA step complete")


def _load_kegg_gene_sets(cfg, gene_id_prefix):
    """Load full KEGG pathway gene sets from cached pathway-gene links.

    Uses the same cache files as step 3c to get ALL genes per pathway,
    not just the DE-overlap subset from enrichment results.
    """
    cache_dir = os.path.join(
        get_path(cfg, "output_dir"),
        "03_differential_analysis", "kegg_enrichment", "cache"
    )
    genes_file = os.path.join(cache_dir, "kegg_pathway_genes.tsv")
    conv_file = os.path.join(cache_dir, "kegg_uniprot_conv.tsv")
    names_file = os.path.join(cache_dir, "kegg_pathway_list.tsv")
    acc_file = os.path.join(cache_dir, "uniprot_genes.tsv")

    if not all(os.path.exists(f) for f in [genes_file, conv_file, acc_file]):
        return {}

    # Build KEGG gene → gene name mapping (same chain as step 3c)
    acc_to_gene = {}
    with open(acc_file) as f:
        for line in f.read().strip().split("\n")[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            acc = parts[0].strip()
            for gn in parts[1].strip().split():
                if gene_id_prefix and gn.startswith(gene_id_prefix):
                    acc_to_gene[acc] = gn
                    break

    uniprot2kegg = {}
    with open(conv_file) as f:
        for line in f.read().strip().split("\n"):
            if "\t" not in line:
                continue
            kegg_id, up_id = line.split("\t", 1)
            uniprot2kegg[up_id.strip().replace("up:", "")] = kegg_id.strip()

    kegg_to_gene = {}
    for acc, gene in acc_to_gene.items():
        kegg_id = uniprot2kegg.get(acc)
        if kegg_id:
            kegg_to_gene[kegg_id] = gene

    # Pathway names
    pathway_names = {}
    if os.path.exists(names_file):
        with open(names_file) as f:
            for line in f.read().strip().split("\n"):
                if "\t" not in line:
                    continue
                pid, name = line.split("\t", 1)
                pid = pid.replace("path:", "")
                pathway_names[pid] = name.split(" - ")[0].strip()

    # Build full pathway → gene set
    gene_sets = {}
    with open(genes_file) as f:
        for line in f.read().strip().split("\n"):
            if "\t" not in line:
                continue
            pid, kegg_gene = line.split("\t", 1)
            pid = pid.replace("path:", "")
            gene = kegg_to_gene.get(kegg_gene.strip())
            if gene:
                name = pathway_names.get(pid, pid)
                gene_sets.setdefault(name, [])
                if gene not in gene_sets[name]:
                    gene_sets[name].append(gene)

    return gene_sets
