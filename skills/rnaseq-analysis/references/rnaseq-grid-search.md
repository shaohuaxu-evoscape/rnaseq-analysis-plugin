# Grid Search

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Systematic parameter optimization for preprocessing and gene filtering.

## Command

```bash
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_grid_search.py --mode <mode> -c configs/analysis_case.yaml [OPTIONS]
```

## Modes

### step0 — Preprocessing Parameters

Sweeps over fastp QC settings to find optimal trimming parameters.

```bash
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_grid_search.py --mode step0 -c configs/analysis_case.yaml
```

**Parameter grid**:

| Parameter | Values |
|-----------|--------|
| `unqualified_percent_limit` | 30, 40, 50 |
| `length_required` | 36, 50, 75 |

Total: 9 combinations.

### step1 — Gene Filtering Parameters

Sweeps over normalization method and filtering thresholds.

```bash
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_grid_search.py --mode step1 -c configs/analysis_case.yaml
```

**Parameter grid**:

| Parameter | Values |
|-----------|--------|
| `normalization.method` | cpm, tpm |
| `gene_filtering.min_expr` | 0.5, 1.0, 2.0, 5.0 |
| `gene_filtering.min_samples` | 1, 2, 3 |

Total: 24 combinations.

## Options

| Flag | Description |
|------|-------------|
| `--collect-only` | Skip execution, only collect and visualize results from existing runs |

## Output Structure

```
results/<run_name>/
└── grid_search/
    ├── step0_preprocessing/
    │   └── <param_combo>/          # One directory per parameter combination
    └── step1_filtering/
        └── <param_combo>/
```

## Workflow

1. Run grid search: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_grid_search.py --mode step1 -c configs/analysis_case.yaml`
2. Review results: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_grid_search.py --mode step1 -c configs/analysis_case.yaml --collect-only`
3. Update `analysis_case.yaml` with the optimal parameters
4. Re-run the full pipeline with the new parameters
