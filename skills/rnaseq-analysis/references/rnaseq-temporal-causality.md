# Temporal Causality

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Time-series cross-correlation analysis across functional gene modules. Computes lag-1 Pearson correlations to infer temporal relationships between metabolic pathways (e.g., TORC1, sterol biosynthesis, fatty acid degradation, mevalonate pathway).

## Step Info

- **ID**: 4b
- **Module**: Advanced Analysis
- **Depends on**: 3a (DE Screening)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 4b`

## Configuration

```yaml
advanced_analysis:
  temporal_causality:
    enabled: true
```

Gene module definitions are configured within the step source code. Override by providing custom gene lists in the config.

## Outputs

Output location: `results/<run_name>/<pair_name>/04_advanced_analysis/temporal_causality/`

| File | Description |
|------|-------------|
| `module_fc_timeseries.tsv` | Mean log2FC per timepoint for each gene module |
| `lag1_correlations.tsv` | Lag-1 Pearson cross-correlations (forward and reverse) |
| `fig_temporal_causality.png` | Two-panel figure: temporal dynamics + lag-1 correlation bars |

## Interpretation

- **Left panel (dynamics)**: Shows how mean fold change evolves across timepoints for each module
- **Right panel (lag-1 bars)**: High positive lag-1 correlation (r > 0.7) suggests module A at time t predicts module B at time t+1
- Asymmetric correlations (forward >> reverse) suggest causal direction
- This analysis requires >= 4 timepoints to be meaningful

## Common Issues

See [Troubleshooting](rnaseq-troubleshooting.md). Configure gene lists to enable this step if default modules are not relevant to your organism.
