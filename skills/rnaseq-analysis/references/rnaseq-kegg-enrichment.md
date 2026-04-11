# KEGG Enrichment

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

KEGG pathway enrichment analysis using Fisher's exact test on DE genes. Uses a gene-to-UniProt-to-KEGG mapping chain and bulk pathway-gene query.

## Step Info

- **ID**: 3c
- **Module**: Differential Analysis
- **Depends on**: 3a (DE Screening)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 3c`

## Configuration

```yaml
differential_analysis:
  kegg_enrichment:
    enabled: true
    organism: ""                   # KEGG organism code -- MUST FILL IN (e.g., "yli", "sce")
    sig_threshold: 0.05
    min_pathway_genes: 5
    max_pathway_genes: 500
    top_n_plot: 20
```

- **organism**: Required. KEGG organism code (e.g., `"yli"` for Y. lipolytica, `"sce"` for S. cerevisiae)
- **min/max_pathway_genes**: Filter pathways by gene set size

## Outputs

Output location: `results/<run_name>/<pair_name>/03_differential_analysis/kegg_enrichment/`

| File | Description |
|------|-------------|
| `kegg_enrichment.tsv` | KEGG pathway enrichment results (pathway, p-value, FDR, overlap genes) |
| `kegg_enrichment.png` | Dot plot of top enriched pathways |
| `cache/kegg_uniprot_conv.tsv` | KEGG-to-UniProt ID conversion |
| `cache/uniprot_yli_genes.tsv` | UniProt accession-to-gene mapping |
| `cache/kegg_pathway_list.tsv` | KEGG pathway names |
| `cache/kegg_pathway_genes.tsv` | KEGG pathway-to-gene links |

## Interpretation

- **Dot plot**: X-axis = gene ratio, Y-axis = pathway names, size = gene count, color = -log10(FDR)
- Pathways related to the experimental perturbation (e.g., nitrogen metabolism, amino acid biosynthesis) should be enriched
- Cross-reference with GO enrichment (step 3b) for consistency
- Check overlap gene lists to identify key genes driving each pathway enrichment

## Common Issues

See [Troubleshooting](rnaseq-troubleshooting.md) §6 for network request failures. KEGG data is cached in `cache/` after first successful download.
