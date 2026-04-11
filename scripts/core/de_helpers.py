"""Shared DE data loading helpers used by enrichment and advanced analysis steps."""

import os
import pandas as pd

from .config_data import get_filtered_expr
from .config_runtime import get_path


def load_de_data(cfg):
    """Load DE genes, background, and direction info.

    Returns
    -------
    tuple
        (de_genes: list, background_genes: list, direction: dict)
    """
    data = cfg.get("_data", {})

    de_results = data.get("de_results")
    if de_results is not None:
        de_genes = list(de_results[de_results["is_de"]].index)
        direction = de_results.loc[de_results["is_de"], "direction"].to_dict()
    else:
        out_root = get_path(cfg, "output_dir")
        path = os.path.join(out_root, "03_differential_analysis", "de_screening", "de_genes.tsv")
        de_df = pd.read_csv(path, sep="\t", index_col=0)
        de_genes = list(de_df.index)
        direction = de_df["direction"].to_dict()

    log2expr = get_filtered_expr(cfg, log2=True)
    background = list(log2expr.index)

    return de_genes, background, direction


def load_de_results(cfg):
    """Load full DE results table (all genes).

    Returns
    -------
    pd.DataFrame or None
    """
    data = cfg.get("_data", {})
    de_results = data.get("de_results")
    if de_results is not None:
        return de_results
    out_root = get_path(cfg, "output_dir")
    path = os.path.join(out_root, "03_differential_analysis", "de_screening", "de_results_all.tsv")
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t", index_col=0)
    return None


def direction_column_name(prefix, condition):
    """Sanitize condition name for use in column headers."""
    safe = "".join(ch if ch.isalnum() else "_" for ch in str(condition)).strip("_")
    return f"{prefix}_in_{safe}"
