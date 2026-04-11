# Cluster Deep-Dive

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Per-cluster enrichment analysis. Runs GO and KEGG enrichment on each gene cluster independently, produces overview figures comparing clusters, and identifies DE gene composition within each cluster.

## Step Info

- **ID**: 5b
- **Module**: Cluster Analysis
- **Depends on**: 5a (Gene Clustering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 5b`

## Configuration

```yaml
cluster_analysis:
  cluster_deepdive:
    enabled: true
    target_clusters: null          # null = all clusters; or list of cluster IDs e.g., [1, 3, 5]
```

- **target_clusters**: Analyze all clusters (null) or specify a subset

## Outputs

Output location: `results/<run_name>/<pair_name>/05_cluster_analysis/cluster_deepdive/`

| File | Description |
|------|-------------|
| `cluster_summary.tsv` | Per-cluster summary (gene count, DE overlap, divergence, GO/KEGG sig counts) |
| `fig_overview.png` | Z-score profiles + DE composition bar charts |
| `fig_go_comparison.png` | Top GO terms per cluster (bar charts) |
| `fig_kegg_comparison.png` | Top KEGG pathways per cluster (bar charts) |
| `C{id}/go_enrichment.tsv` | Per-cluster GO enrichment results |
| `C{id}/kegg_enrichment.tsv` | Per-cluster KEGG enrichment results |
| `C{id}/de_genes.tsv` | DE genes within each cluster |

## Interpretation

- **cluster_summary.tsv**: Quick comparison of cluster sizes, DE gene overlap, and enrichment richness
- **Overview figure**: Top panel shows expression profiles, bottom panel shows DE up/down/ns composition per cluster
- **GO/KEGG comparison figures**: Side-by-side enrichment across clusters reveals functional specialization
- Clusters with high DE overlap and strong enrichment are the most biologically informative
- Per-cluster `C{id}/` directories contain detailed enrichment data for further investigation

## Dependencies

This step reuses cached GO and KEGG annotation data from steps 3b and 3c. If those steps have not been run, the enrichment analysis within the deep-dive will be skipped.
