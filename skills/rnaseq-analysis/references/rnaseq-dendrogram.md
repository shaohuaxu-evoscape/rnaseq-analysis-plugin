# Dendrogram

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Hierarchical clustering dendrogram of samples. Distance metric: 1 - Pearson r, Ward linkage.

## Step Info

- **ID**: 2d
- **Module**: Sample Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 2d`

## Configuration

No step-specific configuration. Uses default distance metric (1 - Pearson r) and Ward linkage.

## Outputs

Output location: `results/<run_name>/<pair_name>/02_sample_analysis/dendrogram/`

| File | Description |
|------|-------------|
| `dendrogram.png` | Hierarchical clustering dendrogram |
| `linkage.tsv` | Linkage matrix (cluster pairs, distances, sizes) |

## Interpretation

- Samples from the same condition should cluster together in a subtree
- Large distances between condition clusters indicate strong differential expression
- If samples from different conditions are intermixed, condition effects may be weak or confounded
- Compare with PCA (step 2c) results for consistency
