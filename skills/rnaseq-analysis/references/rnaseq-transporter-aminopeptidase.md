# Transporter & Aminopeptidase

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Comprehensive transporter census by substrate type, nitrogen transporter expression analysis, and aminopeptidase differential analysis. Classifies transporters using GO annotations and analyzes expression trends across timepoints.

## Step Info

- **ID**: 4a
- **Module**: Advanced Analysis
- **Depends on**: 3a (DE Screening)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 4a`

## Configuration

```yaml
advanced_analysis:
  transporter_aminopeptidase:
    enabled: true
    # Optional organism-specific overrides:
    # amino_info: [...]           # Aminopeptidase annotation list
    # transporter_overrides: {}   # Manual transporter classification overrides
    # substrate_categories: {}    # Substrate category GO regex patterns
    # endo_de_info: {}            # Endopeptidase annotations
```

## Outputs

Output location: `results/<run_name>/<pair_name>/04_advanced_analysis/aminopeptidase/`

| File | Description |
|------|-------------|
| `transporter_census.tsv` | Per-gene transporter stats (category, mean log2FC, DE status, per-timepoint FCs) |
| `fig_transporter.png` | N-transporter category FC trends, DE peptide/amino acid transporters |
| `fig_aminopeptidase.png` | Aminopeptidase vs endopeptidase analysis |

## Interpretation

- **Transporter census**: Shows distribution of transporters by substrate type and their differential expression patterns
- **N-transporter trends**: Time-series fold changes reveal how nitrogen transport shifts across fermentation
- **Aminopeptidase panel**: Compares expression changes between aminopeptidases and endopeptidases to understand proteolytic response
