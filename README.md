# rnaseq-analysis-plugin

Claude Code plugin for RNA-seq time-series differential analysis. Config-driven pipeline with DESeq2 Wald test, GO/KEGG enrichment, GSEA, and K-means clustering.

## Install

```bash
claude /plugin add evoscapebio/rnaseq-analysis-plugin
```

## Quick Start

After installing the plugin, start a new RNA-seq project:

```
/rnaseq-analysis
```

The skill auto-detects your project state:
- **No config found** — launches an interactive setup wizard
- **Config exists, no results** — suggests pipeline execution commands
- **Results exist** — shows overview and offers interpretation

## What's Included

### Skill

The `rnaseq-analysis` skill provides intelligent routing and 21 reference documents covering:

- **3 Operations**: Setup Wizard, Run Pipeline, Grid Search
- **14 Analysis Steps**: Gene Filtering, PCA, DE Screening, GO/KEGG Enrichment, GSEA, Clustering, and more
- **3 Guides**: Config Reference, Step Development, Troubleshooting

### Scripts

18 standalone CLI scripts under `scripts/`:

| Script | Description |
|--------|-------------|
| `rnaseq_run.py` | Run pipeline steps (full, partial, dry-run, resume) |
| `rnaseq_setup.py` | Interactive project setup |
| `rnaseq_grid_search.py` | Parameter optimization |
| `rnaseq_filter.py` | Gene Filtering (step 1a) |
| `rnaseq_de.py` | DE Screening (step 3a) |
| `rnaseq_go.py` | GO Enrichment (step 3b) |
| `rnaseq_kegg.py` | KEGG Enrichment (step 3c) |
| `rnaseq_gsea.py` | GSEA (step 3d) |
| `rnaseq_clustering.py` | Gene Clustering (step 5a) |
| ... | 9 more step scripts |

### Pipeline Steps

| Module | Steps |
|--------|-------|
| Normalization | Gene Filtering |
| Sample Analysis | Expression, Correlation, PCA, Dendrogram |
| Differential Analysis | DE Screening, GO Enrichment, KEGG Enrichment, GSEA |
| Advanced Analysis | Transporter, Temporal Causality, Heterologous Genes, Fermentation |
| Cluster Analysis | Gene Clustering, Cluster Deep-Dive |

## Project Setup

After installing the plugin, create a new analysis project:

```bash
mkdir my-rnaseq-project && cd my-rnaseq-project
```

Then invoke `/rnaseq-analysis`. The skill automatically:
1. Scaffolds the directory structure (`configs/`, `inputs/`, `results/`)
2. Copies the config template
3. Launches an interactive wizard to fill in your project details

## Requirements

Python >= 3.10 with:

```
numpy pandas scipy scikit-learn matplotlib statsmodels PyYAML openpyxl pydeseq2
```

Optional: `gseapy` for GSEA analysis.

## License

MIT
