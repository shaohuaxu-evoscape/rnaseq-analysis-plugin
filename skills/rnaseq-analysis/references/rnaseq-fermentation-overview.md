# Fermentation Overview

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Fermentation condition summary combining growth metrics, product accumulation, sugar metabolism, and specific consumption rates. Reads from an external Excel data file.

## Step Info

- **ID**: 4d
- **Module**: Advanced Analysis
- **Depends on**: -- (no pipeline dependencies)
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 4d`

## Configuration

```yaml
advanced_analysis:
  fermentation_overview:
    enabled: true                  # Default: false (disabled unless configured)
    data_file: "inputs/fermentation_data.xlsx"
    sheet_name: "Sheet1"
    # Optional column name overrides:
    # reactor_id_column: "Reactor_ID"
    # time_column: "Run_time (t)"
    # reactor_ids: ["R01", "R02"]
    # metrics:
    #   biomass_dw: "Dry_Weight (mg/g)"
    #   biomass_reactor_weight: "Reactor_Weight (Pre_Sampling,g)"
    #   product: "Total_AXT (%)"
    #   product_ester: "Esterification_ratio (%)"
    #   substrate: "Res_Sugar (RS, g/L)"
    #   feed: "Sugar_Feed_amount (g)"
```

- **data_file**: Path to Excel file with fermentation data (relative to project root)
- **sheet_name**: Excel sheet name containing the data
- **reactor_ids**: Map to condition1 and condition2 respectively
- **metrics**: Column name mappings (customize for your Excel file format)

## Outputs

Output location: `results/<run_name>/<pair_name>/04_advanced_analysis/fermentation_overview/`

| File | Description |
|------|-------------|
| `fermentation_overview_summary.tsv` | Per-timepoint metrics (OD, DW, biomass, sugar, product, R2/R1 ratios) |
| `fig_fermentation_overview.png` | 2x2 figure: biomass, product + esterification, residual sugar + feed, specific sugar consumption rate |

## Interpretation

- **Top-left (biomass)**: Growth curves for both conditions; divergence indicates differential growth
- **Top-right (product)**: Astaxanthin accumulation and esterification ratio over time
- **Bottom-left (sugar)**: Residual sugar depletion and feed addition patterns
- **Bottom-right (sSCR)**: Specific sugar consumption rate reveals metabolic efficiency differences
