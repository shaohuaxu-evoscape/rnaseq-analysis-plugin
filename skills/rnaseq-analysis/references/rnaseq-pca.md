# PCA

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Principal Component Analysis on filtered log2-normalized expression. Visualizes sample grouping, condition separation, and potential batch effects.

## Step Info

- **ID**: 2c
- **Module**: Sample Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 2c`

## Configuration

```yaml
sample_analysis:
  pca:
    enabled: true
    n_components: 5
    plot_pairs: [[1, 2], [1, 3]]  # PC pairs to plot
```

- **n_components**: Number of principal components to compute
- **plot_pairs**: Each `[x, y]` pair generates a scatter plot of PC_x vs PC_y

## Outputs

Output location: `results/<run_name>/<pair_name>/02_sample_analysis/pca/`

| File | Description |
|------|-------------|
| `pca_variance.tsv` | Per-PC variance explained and cumulative |
| `pca_scores.tsv` | Sample scores in PC space |
| `pca_pc{x}_pc{y}.png` | Scatter plots for each configured PC pair |

## Interpretation

- PC1 should capture the most variance; check if it separates conditions or timepoints
- If PC1 separates batches instead of conditions, batch effects may be significant
- Samples from the same condition/timepoint should cluster together
- High variance on PC1 (> 30%) is typical for well-separated conditions
