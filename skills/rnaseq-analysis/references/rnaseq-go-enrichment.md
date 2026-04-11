# GO Enrichment

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Gene Ontology enrichment analysis using Fisher's exact test on DE genes. Tests Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) ontologies. Uses UniProt two-step mapping: accession-to-gene + accession-to-GO.

## Step Info

- **ID**: 3b
- **Module**: Differential Analysis
- **Depends on**: 3a (DE Screening)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 3b`

## Configuration

```yaml
differential_analysis:
  go_enrichment:
    enabled: true
    organism_taxid: 0              # NCBI taxonomy ID -- MUST FILL IN
    ontologies: ["BP", "MF", "CC"]
    sig_threshold: 0.05
    min_term_genes: 5
    max_term_genes: 500
    top_n_plot: 20
```

- **organism_taxid**: Required. Examples: 4952 (Y. lipolytica), 559292 (S. cerevisiae S288C)
- **min/max_term_genes**: Filter GO terms by gene set size to avoid overly broad or specific terms
- **top_n_plot**: Number of top terms to show in the bubble plot

Also requires `project.gene_id_prefix` to be set (e.g., `"YALI1_"`) for UniProt-to-gene matching.

## Outputs

Output location: `results/<run_name>/<pair_name>/03_differential_analysis/go_enrichment/`

| File | Description |
|------|-------------|
| `go_enrichment.tsv` | GO enrichment results (term, ontology, p-value, FDR, overlap gene list) |
| `go_enrichment.png` | Combined bubble plot across BP/MF/CC |
| `cache/uniprot_yli_genes.tsv` | Cached UniProt accession-to-gene mapping |
| `cache/uniprot_yli_go.tsv` | Cached UniProt accession-to-GO mapping |

## Interpretation

- **Bubble plot**: X-axis = gene ratio, Y-axis = GO terms, size = gene count, color = -log10(FDR)
- Focus on BP terms for biological insight, MF for molecular mechanisms, CC for subcellular localization
- Check the overlap gene list in the TSV for genes driving each enrichment
- Compare up-regulated vs. down-regulated enrichments (direction column)

## Common Issues

See [Troubleshooting](rnaseq-troubleshooting.md) §6 for network request failures. UniProt/GO data is cached in `cache/` after first successful download.
