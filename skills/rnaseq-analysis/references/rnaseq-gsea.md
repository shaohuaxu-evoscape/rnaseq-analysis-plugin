# GSEA

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Gene Set Enrichment Analysis using gseapy. Uses custom KEGG gene sets mapped to organism-specific IDs (not human gene sets). Falls back to saving the ranking file for external analysis if gseapy is unavailable.

## Step Info

- **ID**: 3d
- **Module**: Differential Analysis
- **Depends on**: 3a (DE Screening)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 3d`

## Configuration

```yaml
differential_analysis:
  gsea:
    enabled: true
    ranking_metric: "mean_log2fc"  # Column from de_results_all.tsv used for ranking
    min_size: 10                   # Minimum gene set size
    max_size: 500                  # Maximum gene set size
    permutations: 1000             # Number of permutations
```

- **ranking_metric**: `"mean_log2fc"` (default) or `"signed_neg_log10p"` for signed significance ranking
- Uses KEGG pathway gene sets from step 3c cache (not external GMT files)

## Outputs

Output location: `results/<run_name>/<pair_name>/03_differential_analysis/gsea/`

| File | Description |
|------|-------------|
| `gsea_ranking.tsv` | Ranked gene list (gene_id, ranking_score) |
| `gsea_results.tsv` | GSEA results (pathway, NES, p-value, FDR q-value, leading edge genes) |
| `gseapy_output/` | gseapy output directory with enrichment plots |

## Interpretation

- **NES (Normalized Enrichment Score)**: Positive = enriched in upregulated genes, negative = enriched in downregulated genes
- **FDR q-value**: Multiple testing corrected significance
- Focus on pathways with |NES| > 1.5 and FDR < 0.25 (standard GSEA threshold)
- GSEA complements Fisher's exact test (steps 3b/3c) by detecting coordinated but subtle expression changes
- Leading edge genes are the core subset driving the enrichment signal
