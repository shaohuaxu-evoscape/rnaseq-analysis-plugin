# DE Screening

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Differential expression screening using DESeq2 Wald test on raw counts. Produces genome-wide DE statistics, significant gene lists, volcano and MA plots.

## Step Info

- **ID**: 3a
- **Module**: Differential Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 3a`

## Configuration

```yaml
differential_analysis:
  de_screening:
    enabled: true
    method: "deseq2_wald"          # DESeq2 Wald test; "ols" for legacy fallback
    fdr_threshold: 0.05
    log2fc_threshold: 0.585        # = log2(1.5), i.e., 1.5-fold change
    fdr_method: "bh"               # Benjamini-Hochberg
```

### DESeq2 Model (auto-detected)

- **Crossed design** (one batch contains both conditions): `~ batch_id + condition + timepoint`
- **Nested design** (each batch = one condition): `~ condition + timepoint`

No manual configuration needed; the step auto-detects the design from the sample manifest.

## Outputs

Output location: `results/<run_name>/<pair_name>/03_differential_analysis/de_screening/`

| File | Description |
|------|-------------|
| `de_results_all.tsv` | Full DE statistics for all tested genes (log2FC, FDR, direction, per-timepoint diffs) |
| `de_genes.tsv` | Significant DE gene subset (passing FDR and FC thresholds) |
| `de_summary.txt` | Analysis parameters and gene counts (total, up, down) |
| `volcano_plot.png` | Mean log2FC vs -log10(FDR) |
| `ma_plot.png` | Mean expression vs mean log2FC |

## Interpretation

- **de_results_all.tsv columns**: gene_id, mean_log2fc, fdr, direction (up/down/ns), per-timepoint log2FC values
- **Volcano plot**: Genes in top-left and top-right quadrants are significantly DE; horizontal dashed line = FDR threshold, vertical lines = FC threshold
- **MA plot**: Shows whether fold changes are expression-dependent (they should not be)
- **de_summary.txt**: Quick check for total DE count and up/down balance

## In-Memory Data

| Key | Description |
|-----|-------------|
| `de_results` | Full DE results DataFrame (all genes) |
| `de_genes` | Significant DE genes DataFrame |

## Common Issues

See [Troubleshooting](rnaseq-troubleshooting.md):
- §1: Singular design matrix
- §2: Missing timepoints
