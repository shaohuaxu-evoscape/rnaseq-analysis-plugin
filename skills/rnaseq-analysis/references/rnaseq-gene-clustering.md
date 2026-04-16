# Gene Clustering

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

K-means clustering on top variable genes. Identifies co-expression patterns across conditions and timepoints. Includes automatic K selection via elbow method and silhouette analysis.

## Step Info

- **ID**: 5a
- **Module**: Cluster Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 5a`

## Configuration

```yaml
cluster_analysis:
  gene_clustering:
    enabled: true
    top_n_genes: 2000              # Top N genes by variance for clustering
    k_range: [4, 12]               # K-means search range [min, max]
    fixed_k: 7                     # Fixed K (overrides auto-selection if set)
    random_state: 42
```

- **top_n_genes**: Higher values include more genes but may dilute signal
- **k_range**: Range to search for optimal K (silhouette + elbow)
- **fixed_k**: Set to skip auto-selection and use a specific K

## Outputs

Output location: `results/<run_name>/<pair_name>/05_cluster_analysis/gene_clustering/`

| File | Description |
|------|-------------|
| `k_selection.png` | Elbow plot (inertia) + silhouette score for K range |
| `cluster_assignments.tsv` | Gene ID and 1-indexed cluster number |
| `cluster_profiles.png` | Per-cluster Z-score expression profiles with condition comparison |
| `cluster_heatmap.png` | Z-scored heatmap with cluster sidebar annotation |

## Interpretation

- **K selection plot**: Choose K where silhouette score peaks or elbow flattens
- **Cluster profiles**: Each cluster shows a distinct expression pattern; compare condition1 vs condition2 profiles
- **Heatmap**: Visual overview of all clusters; clear block structure indicates good separation
- Clusters with opposing profiles between conditions contain the most biologically interesting genes
- Use `cluster_assignments.tsv` as input for Cluster Deep-Dive (step 5b) enrichment analysis
