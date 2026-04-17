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

## Behavior Routing

Based on the detection results above, choose the corresponding behavior:

### Case A: New Project (analysis_case.yaml is NOT_FOUND)

**First, determine the project folder name:**
- Extract a short identifier from $ARGUMENTS (e.g. "run_test" from "创建1个run_test的转录组分析项目")
- If no name can be parsed, use `rna_project` as the default
- Create and scaffold inside that subfolder:

```bash
PROJECT_DIR=<parsed_name_or_rna_project>
mkdir -p ${PROJECT_DIR}/configs ${PROJECT_DIR}/inputs/ref ${PROJECT_DIR}/results
cp ${CLAUDE_PLUGIN_ROOT}/configs/analysis_case.template.yaml ${PROJECT_DIR}/configs/analysis_case.yaml
```

Then enter the interactive setup wizard. Read [`references/wizard-1-basics.md`](references/wizard-1-basics.md) and follow the staged protocol. All subsequent file operations (config read/write, data scan, results) use `${PROJECT_DIR}/` as the base. Each stage file directs to the next — load only the current stage.

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

## Reference Index

Step references (read on demand):
- `1a` Gene Filtering → `references/rnaseq-gene-filtering.md`
- `2a` Expression / `2b` Correlation / `2c` PCA / `2d` Dendrogram → `references/rnaseq-{name}.md`
- `3a` DE Screening / `3b` GO / `3c` KEGG / `3d` GSEA → `references/rnaseq-{name}.md`
- `4a` Transporter / `4b` Temporal Causality / `4c` Heterologous / `4d` Fermentation → `references/rnaseq-{name}.md`
- `5a` Gene Clustering / `5b` Cluster Deep-Dive → `references/rnaseq-{name}.md`

Operation references:
- Setup Wizard → `references/wizard-1-basics.md` (staged: 1-basics → 2-data → 3-remote → 4-analysis → 5-confirm)
- Run Pipeline → `references/rnaseq-run-pipeline.md`
- Grid Search → `references/rnaseq-grid-search.md`
- Remote Preprocessing → `references/rnaseq-remote-preprocessing.md`
- Config Reference → `references/rnaseq-config-reference.md`
- Extend Step → `references/rnaseq-extend-step.md`
- Troubleshooting → `references/rnaseq-troubleshooting.md`
- Overview (concepts, architecture) → `references/rnaseq-overview.md`
