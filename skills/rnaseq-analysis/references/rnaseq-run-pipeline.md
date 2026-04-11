# Run Pipeline

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Execute the RNA-seq analysis pipeline. Requires a valid `configs/analysis_case.yaml`.

## Command

```bash
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml [OPTIONS]
```

## Options

| Flag | Description |
|------|-------------|
| `-c, --config <path>` | Path to analysis_case.yaml (required) |
| `-s, --steps <spec>` | Step selection (default: `1a-2d`) |
| `-n, --dry-run` | Preview execution plan without running |
| `--resume` | Resume from pipeline_run_manifest.json checkpoint |
| `-l, --list-steps` | List all available steps and exit |

## Step Selection Syntax

| Format | Example | Description |
|--------|---------|-------------|
| Range | `1a-2d` | All steps from 1a to 2d inclusive |
| List | `1a,3b,5a` | Specific steps (comma-separated) |
| All | `all` | Every public step |
| Module range | `1-3` | All steps in modules 1 through 3 |
| Module name | `normalization` | All steps in the named module |

## Common Patterns

```bash
# Full pipeline
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps all

# Normalization + sample QC
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-2d

# Differential analysis + enrichment
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 3a-3d

# Single step (dependencies auto-resolved)
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 3b

# Dry run (preview only)
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --dry-run

# Resume interrupted run
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --resume

# List available steps
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps
```

## Auto Dependency Resolution

Requesting a step automatically includes all its dependencies. For example:
- `--steps 3b` (GO Enrichment) auto-adds `1a` (Gene Filtering) + `3a` (DE Screening)
- `--steps 5b` (Cluster Deep-Dive) auto-adds `1a` + `5a` (Gene Clustering)

Already-completed steps are skipped when using `--resume`.

## DESeq2 Model Auto-Detection

The DE screening step (3a) automatically selects the statistical model:

- **Crossed design** (a single batch contains both conditions): `~ batch_id + condition + timepoint`
- **Nested design** (each batch = one condition): `~ condition + timepoint`

No manual configuration needed.

## Resume Behavior

`--resume` reads `pipeline_run_manifest.json` and skips completed steps. Important:
- Config fingerprint (SHA-256) must match -- if `analysis_case.yaml` changed, resume starts fresh
- Changing `--steps` with the same config is fine
- See [Troubleshooting](rnaseq-troubleshooting.md) §7 for resume issues

## Output Structure

```
results/
├── <run_name>/
│   ├── analysis_config.yaml          # Resolved runtime config
│   ├── sample_manifest.tsv           # Sample metadata
│   ├── analysis_gene_counts.tsv      # Merged count matrix
│   ├── pipeline_run_manifest.json    # Checkpoint for --resume
│   ├── documents/                    # Auto-generated reports
│   └── <pair_name>/
│       ├── 01_normalization/
│       ├── 02_sample_analysis/
│       ├── 03_differential_analysis/
│       ├── 04_advanced_analysis/
│       └── 05_cluster_analysis/
└── shared/
    └── <batch_id>/
        └── 01_preprocessing/
```
