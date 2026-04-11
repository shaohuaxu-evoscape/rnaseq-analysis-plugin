# Setup Wizard

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Interactive project setup that generates `configs/analysis_case.yaml` from user input. Automatically triggered when no config file exists (Behavior Routing Case A).

## Pre-conditions

- `configs/analysis_case.yaml` does NOT exist
- Template available at `configs/analysis_case.template.yaml`

## Principle

1-2 questions per round, with sensible defaults. Auto-detect what can be inferred from local files to minimize user input.

## Round 1: Project Basics

Ask: project name, organism, strain (optional).

## Round 2: Data Sources

Scan local `inputs/` and `results/shared/`:
- Found `.fa` + `.gtf` -> auto-fill reference genome paths
- Found `gene_counts.tsv` -> auto-fill count matrix path
- None found -> ask for FASTQ directory path

## Round 3: Condition Definition

If gene_counts.tsv is available, read the header to list available conditions:

```python
# Column format: {condition}-{timepoint}, e.g., R1-18, R2-48
header = open(gene_counts_path).readline().strip().split('\t')[1:]
conditions = sorted(set(col.rsplit('-', 1)[0] for col in header))
timepoints = sorted(set(int(col.rsplit('-', 1)[1]) for col in header))
```

Let the user choose which two conditions to compare.

## Round 4: Timepoints

Show detected timepoints, ask whether to use all or select a subset.

## Round 5: Enrichment Analysis

Ask: NCBI Taxonomy ID and KEGG organism code. Can be skipped if unsure.

## Round 6: Advanced Analysis (Optional)

Ask: whether to enable heterologous gene analysis, fermentation overview, transporter analysis. Default: all disabled.

## Round 7: Confirm & Generate

Show config preview -> user confirms -> Write generates `configs/analysis_case.yaml` -> suggest next commands.

When generating the config, populate values into the structure from `configs/analysis_case.template.yaml`.

## Post-Setup

After generating the config, suggest:

```
Configuration saved. Next steps:
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --dry-run   # preview
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b             # full run
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps                                             # list steps
```

See [Run Pipeline](rnaseq-run-pipeline.md) for execution details.
