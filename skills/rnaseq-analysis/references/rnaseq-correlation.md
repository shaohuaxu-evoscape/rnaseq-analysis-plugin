# Correlation Analysis

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Sample-to-sample correlation heatmap. Evaluates biological reproducibility and condition separation.

## Step Info

- **ID**: 2b
- **Module**: Sample Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 2b`

## Configuration

```yaml
sample_analysis:
  correlation:
    enabled: true
    method: "pearson"              # "pearson" or "spearman"
    plot_upper_triangle: true
    annot_fontsize: 8
```

## Outputs

Output location: `results/<run_name>/<pair_name>/02_sample_analysis/correlation/`

| File | Description |
|------|-------------|
| `correlation_matrix.tsv` | Sample-by-sample correlation matrix |
| `correlation_heatmap.png` | Upper-triangle heatmap with annotated values |
| `correlation_summary.txt` | Within-condition and between-condition statistics |

## Interpretation

- Within-condition correlations should be high (> 0.95 for Pearson on log2 data)
- Between-condition correlations should be lower than within-condition
- `correlation_summary.txt` reports mean within-group and between-group correlations
- Low within-condition correlation may indicate batch effects or outlier samples
