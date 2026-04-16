---
name: "rnaseq-analysis"
description: "Use when the user wants RNA-seq or transcriptomics workflow help in this repo: project setup, analysis-case config authoring, pipeline execution, DESeq2-based differential expression, GO/KEGG enrichment, GSEA, clustering, result interpretation, remote preprocessing, grid search, or extending pipeline steps."
---

# RNA-seq Analysis Skill

## When to use
- Set up a new RNA-seq analysis project in this repository.
- Create or edit `configs/analysis_case.yaml`.
- Run or debug the pipeline under `scripts/`.
- Interpret outputs such as DE tables, volcano plots, PCA, correlation, enrichment, or clustering.
- Add or modify pipeline steps.

## Core workflow
1. Detect project state inside the current workspace:
   - Check whether `configs/analysis_case.yaml` exists.
   - Check whether `results/*/pipeline_run_manifest.json` exists.
   - Check for `inputs/ref/*`, `results/shared/**/gene_counts.tsv`, and FASTQ directories under `inputs/`.
2. Choose the path that matches the state:
   - No config: scaffold `configs/`, `inputs/ref/`, and `results/`, then follow the setup wizard reference.
   - Config exists but no results: prepare or validate inputs, then run the pipeline.
   - Results exist: summarize run status and read the relevant step reference before interpreting or extending anything.
3. Prefer the Python entrypoints in this repo:
   - `python scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps ...`
   - `python scripts/rnaseq_setup.py`
   - `python scripts/rnaseq_grid_search.py --mode ... -c configs/analysis_case.yaml`
4. Before explaining a step, interpreting a result, or editing a workflow, read the matching file in `references/`.

## Entry points
- Full pipeline: `scripts/rnaseq_run.py`
- Setup helper: `scripts/rnaseq_setup.py`
- Grid search: `scripts/rnaseq_grid_search.py`
- Individual step runners: `scripts/rnaseq_*.py`
- Core execution logic: `scripts/core/`
- Step implementations: `scripts/steps/`

## Reference routing
- Setup or config authoring: read `references/rnaseq-setup-wizard.md` and `references/rnaseq-config-reference.md`
- Running or resuming the pipeline: read `references/rnaseq-run-pipeline.md`
- Remote preprocessing or deployment: read `references/rnaseq-remote-preprocessing.md` or `references/rnaseq-remote-deployment.md`
- Parameter sweeps: read `references/rnaseq-grid-search.md`
- Troubleshooting: read `references/rnaseq-troubleshooting.md`
- Step-specific interpretation:
  - `1a` gene filtering: `references/rnaseq-gene-filtering.md`
  - `2a-2d` sample analysis: expression, correlation, PCA, dendrogram references
  - `3a-3d` differential analysis: DE, GO, KEGG, GSEA references
  - `4a-4d` advanced analysis: transporter, temporal, heterologous, fermentation references
  - `5a-5b` clustering: gene clustering and cluster deep-dive references
- Adding or changing a step: read `references/rnaseq-extend-step.md`

## Operating rules
- Stay inside the current project directory; do not assume sibling projects or external paths.
- Treat `configs/analysis_case.yaml` as the sole user-edited analysis entry point.
- Prefer `python` CLI entrypoints over ad hoc imports when executing the pipeline.
- When remote preprocessing is enabled, use `scripts/rnaseq_remote_preprocess.sh` first and fall back to the remote preprocessing reference only if needed.
- When interpreting outputs, cite concrete files under `results/` rather than paraphrasing from memory.

## Validation
- After changing pipeline code or configs, run the narrowest relevant test first, then broader validation if needed.
- For repository-wide validation, use `pytest -q`.

