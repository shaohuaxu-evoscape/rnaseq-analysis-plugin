# Heterologous Genes

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Expression analysis of engineered heterologous/transgenes. Generates timecourse expression panels and expression-vs-fold-change scatter plots. Supports fusion gene rules where a fused construct replaces individual parts.

## Step Info

- **ID**: 4c
- **Module**: Advanced Analysis
- **Depends on**: 1a (Gene Filtering)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 4c`

## Configuration

```yaml
advanced_analysis:
  heterologous_genes:
    enabled: true                  # Default: false (disabled unless configured)
    gene_list:
      - name: "GeneA"
        function: "Description of gene function"
        color: "#D32F2F"
      - name: "GeneB"
        function: "Another description"
        color: "#1976D2"
    fusion_rules:
      - fusion: "FusedAB"
        parts: ["GeneA", "GeneB"]  # Skip parts when fusion ID is present
```

- **gene_list**: Required for this step. Each entry has name (matching gene ID), function description, and plot color
- **fusion_rules**: Optional. When the fusion gene ID is found, its parts are skipped from individual analysis

## Outputs

Output location: `results/<run_name>/<pair_name>/04_advanced_analysis/heterologous_genes/`

| File | Description |
|------|-------------|
| `heterologous_genes_summary.tsv` | Per-gene summary (function, mean expression, mean log2FC, p-value, per-timepoint values) |
| `fig_heterologous_timecourse.png` | Multi-panel timecourse (expression per condition across timepoints) |
| `fig_heterologous_expr_vs_fc.png` | Scatter plots: expression level vs log2 fold change |

## Interpretation

- **Timecourse panels**: Each gene gets its own panel showing both conditions over time
- **Expression vs FC scatter**: Genes in the upper-right quadrant are highly expressed AND upregulated
- Low or zero expression of a heterologous gene may indicate silencing, loss, or poor codon optimization
