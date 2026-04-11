# Expression Summary

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Polar bar chart of detected genes per sample by expression level. Provides a quick visual overview of sample quality and expression distribution.

## Step Info

- **ID**: 2a
- **Module**: Sample Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 2a`

## Configuration

```yaml
sample_analysis:
  expression_summary:
    enabled: true
    stats: ["mean", "median", "cv"]
```

## Outputs

Output location: `results/<run_name>/<pair_name>/02_sample_analysis/expression/`

| File | Description |
|------|-------------|
| `detected_genes.tsv` | Per-sample detection counts across expression bins |
| `detected_genes_polar.png` | Polar bar chart stacked by expression level |

## Interpretation

- Each spoke represents one sample; stacked bars show gene counts in expression bins
- Samples with similar total detected genes and similar bin distributions indicate good consistency
- Outlier samples (significantly fewer detected genes) may indicate quality issues
