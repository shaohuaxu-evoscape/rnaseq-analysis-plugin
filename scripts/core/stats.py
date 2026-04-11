"""Shared statistical utilities."""

import numpy as np
from statsmodels.stats.multitest import multipletests


def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvalues)
    if n == 0:
        return np.array([])
    _, fdr, _, _ = multipletests(pvalues, method="fdr_bh")
    return fdr
