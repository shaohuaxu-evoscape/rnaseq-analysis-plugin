---
name: rnaseq-analysis
version: 1.0.0
description: >
  RNA-seq time-series differential analysis pipeline. DESeq2 Wald test, GO/KEGG
  enrichment, GSEA, clustering. Supports interactive project setup, execution,
  result interpretation, step development, and troubleshooting.
  TRIGGER when: user mentions RNA-seq, transcriptomics, differential expression,
  DE analysis, gene enrichment, GO/KEGG enrichment, GSEA, gene clustering,
  initializing or setting up an analysis project, running the pipeline,
  interpreting DE/enrichment results, adding analysis steps, or asking about
  gene counts, volcano plots, PCA, or sample correlation.
argument-hint: "Leave empty for auto mode, or describe your need: run pipeline, interpret results, add step, etc."
---

# RNA-seq Analysis Pipeline (v1)

**CRITICAL — Before executing any operation, use Read to load the linked reference file.**

User request: $ARGUMENTS

## Project State Detection

### analysis_case.yaml
```yaml
!`cat configs/analysis_case.yaml 2>/dev/null || echo "NOT_FOUND"`
```

### Existing Results
```
!`ls results/*/pipeline_run_manifest.json 2>/dev/null | head -5 || echo "NO_RESULTS"`
```

### Local Data Scan (project folder only)
```
!`echo "=== Reference genomes ===" && find inputs/ref -name "*.fa" -o -name "*.fasta" 2>/dev/null | head -5; echo "=== Reference annotations ===" && find inputs/ref -name "*.gtf" -o -name "*.gff" 2>/dev/null | head -5; echo "=== Gene counts ===" && find inputs -maxdepth 2 -name "gene_counts*.tsv" -o -name "gene_counts*.csv" 2>/dev/null | head -5; find results/shared -name "gene_counts.tsv" 2>/dev/null | head -3; echo "=== FASTQ dirs ===" && ls -d inputs/fastq*/ inputs/raw*/ 2>/dev/null | head -5 || echo "none"`
```

> **Scope rule:** This scan is strictly limited to the current project directory. Do NOT use paths from external context (CLAUDE.md, project-level knowledge, sibling directories, or any source outside this project folder).

### Remote Preprocessing Config
```
!`python3 -c "import yaml; c=yaml.safe_load(open('configs/analysis_case.yaml')); r=c.get('remote',{}); print(f'enabled={r.get(\"enabled\",False)} host={r.get(\"host\",\"\")} data_dir={r.get(\"data_dir\",\"\")}') if r else print('NOT_CONFIGURED')" 2>/dev/null || echo "NOT_CONFIGURED"`
```

---

## Plugin Scripts

This plugin provides per-step CLI scripts at `${CLAUDE_PLUGIN_ROOT}/scripts/`.

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

---

## Behavior Routing

Based on the detection results above, choose the corresponding behavior:

### Case A: New Project (analysis_case.yaml is NOT_FOUND)

**First, scaffold the project directory structure:**

```bash
mkdir -p configs inputs/ref results
cp ${CLAUDE_PLUGIN_ROOT}/configs/analysis_case.template.yaml configs/analysis_case.yaml
```

Then enter the interactive setup wizard. Read [`references/rnaseq-setup-wizard.md`](references/rnaseq-setup-wizard.md) and follow the 8-round protocol to fill in the config.

### Case B: Configured but Not Run (yaml exists, no results)

**Check if remote preprocessing is needed:** If `remote.enabled = true` AND no gene_counts.tsv exists locally:

1. **Primary (fast)**: Run the SSH preprocessing script:
   ```bash
   bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
   ```
2. **Fallback (interactive debug)**: If the script fails or the user wants step-by-step control, read [`references/rnaseq-remote-preprocessing.md`](references/rnaseq-remote-preprocessing.md) and execute via MCP remote-linux.

After gene_counts.tsv is pulled back, continue with local analysis.

**Otherwise** (local mode or gene_counts already available), show config summary and suggest:
```
Configuration is ready. Suggestions:
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --dry-run
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b
```
Read [`references/rnaseq-run-pipeline.md`](references/rnaseq-run-pipeline.md) for execution details.

### Case C: Results Exist

Show run status, then act based on user request:
- No arguments -> show results overview (DE count, step status)
- "interpret" / "results" -> read the relevant step reference for interpretation guidance
- "run X" -> read [`references/rnaseq-run-pipeline.md`](references/rnaseq-run-pipeline.md)
- "add step" -> read [`references/rnaseq-extend-step.md`](references/rnaseq-extend-step.md)

### Case D: User Has Explicit Request ($ARGUMENTS is non-empty)

Handle the request directly by reading the relevant reference file below.

---

## Core Concepts

- **Analysis Case**: A pairwise DE comparison defined by `configs/analysis_case.yaml`
- **Module**: A group of related analysis steps (6 modules, numbered 0-5)
- **Step**: An atomic analysis unit within a module, exports `run(cfg)`
- **Hidden Preprocessing**: Steps 0a-0d auto-trigger when gene_counts.tsv is missing. When `remote.enabled` is true, these run on a remote server via MCP remote-linux — see [`Remote Preprocessing`](references/rnaseq-remote-preprocessing.md)
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

---

## Operations

| Operation | Description |
|-----------|-------------|
| [`Setup Wizard`](references/rnaseq-setup-wizard.md) | Interactive project setup — generates analysis_case.yaml |
| [`Run Pipeline`](references/rnaseq-run-pipeline.md) | Execute analysis steps (full, partial, dry-run, resume) |
| [`Grid Search`](references/rnaseq-grid-search.md) | Parameter optimization for filtering and preprocessing |
| [`Remote Preprocessing`](references/rnaseq-remote-preprocessing.md) | Run steps 0a-0d on remote server via MCP remote-linux |

---

## Pipeline Steps

### Normalization

| Step | ID | Description | Depends On |
|------|----|-------------|------------|
| [`Gene Filtering`](references/rnaseq-gene-filtering.md) | 1a | CPM/TPM normalization + low-expression filtering | -- |

### Sample Analysis

| Step | ID | Description | Depends On |
|------|----|-------------|------------|
| [`Expression Summary`](references/rnaseq-expression-summary.md) | 2a | Per-sample detected genes by expression bin | 1a |
| [`Correlation`](references/rnaseq-correlation.md) | 2b | Sample-to-sample correlation heatmap | 1a |
| [`PCA`](references/rnaseq-pca.md) | 2c | Principal component analysis | 1a |
| [`Dendrogram`](references/rnaseq-dendrogram.md) | 2d | Hierarchical clustering dendrogram | 1a |

### Differential Analysis

| Step | ID | Description | Depends On |
|------|----|-------------|------------|
| [`DE Screening`](references/rnaseq-de-screening.md) | 3a | DESeq2 Wald test, volcano/MA plots | 1a |
| [`GO Enrichment`](references/rnaseq-go-enrichment.md) | 3b | Gene Ontology (BP/MF/CC) enrichment | 3a |
| [`KEGG Enrichment`](references/rnaseq-kegg-enrichment.md) | 3c | KEGG pathway enrichment | 3a |
| [`GSEA`](references/rnaseq-gsea.md) | 3d | Gene Set Enrichment Analysis | 3a |

### Advanced Analysis

| Step | ID | Description | Depends On |
|------|----|-------------|------------|
| [`Transporter & Aminopeptidase`](references/rnaseq-transporter-aminopeptidase.md) | 4a | Substrate classification and expression trends | 3a |
| [`Temporal Causality`](references/rnaseq-temporal-causality.md) | 4b | Time-series cross-correlation | 3a |
| [`Heterologous Genes`](references/rnaseq-heterologous-genes.md) | 4c | Engineered transgene expression analysis | 1a |
| [`Fermentation Overview`](references/rnaseq-fermentation-overview.md) | 4d | Fermentation metrics summary | -- |

### Cluster Analysis

| Step | ID | Description | Depends On |
|------|----|-------------|------------|
| [`Gene Clustering`](references/rnaseq-gene-clustering.md) | 5a | K-means clustering on top variable genes | 1a |
| [`Cluster Deep-Dive`](references/rnaseq-cluster-deepdive.md) | 5b | Per-cluster enrichment and expression | 5a |

**Step selection syntax**: `all` | `1a-2d` | `1a,3b,5a` | `normalization` | `1-3`

---

## Guides

| Guide | Description |
|-------|-------------|
| [`Config Reference`](references/rnaseq-config-reference.md) | Per-field documentation for analysis_case.yaml |
| [`Extend: Add Step`](references/rnaseq-extend-step.md) | Step development guide, templates, API reference |
| [`Remote Deployment`](references/rnaseq-remote-deployment.md) | Deploy shared pipeline on remote server |
| [`Troubleshooting`](references/rnaseq-troubleshooting.md) | Common errors and fixes |
