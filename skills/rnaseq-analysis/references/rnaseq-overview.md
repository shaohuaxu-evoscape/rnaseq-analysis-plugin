# RNA-seq Pipeline Overview

## Core Concepts

- **Analysis Case**: A pairwise DE comparison defined by `configs/analysis_case.yaml`
- **Module**: A group of related analysis steps (6 modules, numbered 0-5)
- **Step**: An atomic analysis unit within a module, exports `run(cfg)`
- **Hidden Preprocessing**: Steps 0a-0d auto-trigger when gene_counts.tsv is missing. When `remote.enabled` is true, these run on a remote server via MCP remote-linux — see `rnaseq-remote-preprocessing.md`
- **Config**: `configs/analysis_case.yaml` is the sole user-editable entry point

## Pipeline Architecture

```
analysis_case.yaml -> prepare_analysis_case -> runner -> step modules
                                                         |
Module 0 (hidden): HISAT2 Index -> Fastp -> HISAT2 Align -> HTSeq Count
Module 1: Gene Filtering
Module 2: Expression | Correlation | PCA | Dendrogram
Module 3: DE Screening -> GO Enrichment | KEGG Enrichment | GSEA
Module 4: Transporter | Temporal Causality | Heterologous Genes | Fermentation
Module 5: Gene Clustering -> Cluster Deep-Dive
```

## Plugin Scripts

```bash
# Run full pipeline
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b

# Run a single step
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_de.py -c configs/analysis_case.yaml

# List available steps
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps

# Project setup
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_setup.py
```

## Pipeline Steps

**Step selection syntax**: `all` | `1a-2d` | `1a,3b,5a` | `normalization` | `1-3`

| ID | Step | Description | Depends On |
|----|------|-------------|------------|
| 1a | Gene Filtering | CPM/TPM normalization + low-expression filtering | -- |
| 2a | Expression Summary | Per-sample detected genes by expression bin | 1a |
| 2b | Correlation | Sample-to-sample correlation heatmap | 1a |
| 2c | PCA | Principal component analysis | 1a |
| 2d | Dendrogram | Hierarchical clustering dendrogram | 1a |
| 3a | DE Screening | DESeq2 Wald test, volcano/MA plots | 1a |
| 3b | GO Enrichment | Gene Ontology (BP/MF/CC) enrichment | 3a |
| 3c | KEGG Enrichment | KEGG pathway enrichment | 3a |
| 3d | GSEA | Gene Set Enrichment Analysis | 3a |
| 4a | Transporter & Aminopeptidase | Substrate classification and expression trends | 3a |
| 4b | Temporal Causality | Time-series cross-correlation | 3a |
| 4c | Heterologous Genes | Engineered transgene expression analysis | 1a |
| 4d | Fermentation Overview | Fermentation metrics summary | -- |
| 5a | Gene Clustering | K-means clustering on top variable genes | 1a |
| 5b | Cluster Deep-Dive | Per-cluster enrichment and expression | 5a |
