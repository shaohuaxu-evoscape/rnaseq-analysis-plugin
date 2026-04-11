# Gene Filtering

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Remove low-expression genes. Applies expression-based filter (CPM/TPM/FPKM) per condition and exports filtered count/normalized matrices. This is the foundation step -- all downstream analysis depends on its output.

## Step Info

- **ID**: 1a
- **Module**: Normalization
- **Depends on**: --
- **Run command**: `python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a`

## Configuration

```yaml
normalization:
  method: "tpm"          # "cpm" | "tpm" | "fpkm"
  log_transform: true
  pseudocount: 1         # log2(expression + pseudocount)
  gene_filtering:
    min_expr: 1.0        # Minimum expression threshold (in normalized units)
    min_samples: 2       # Minimum samples per condition passing threshold
```

- **method**: TPM/FPKM require gene lengths from `paths.reference_gtf`; CPM does not
- **gene_filtering logic**: A gene passes if in at least one condition, >= `min_samples` samples have expression >= `min_expr`

## Outputs

Output location: `results/<run_name>/<pair_name>/01_normalization/gene_filtering/`

| File | Description |
|------|-------------|
| `filtered_counts.tsv` | Raw counts after gene filtering (genes x samples) |
| `filtered_{method}.tsv` | Normalized expression matrix |
| `filtered_log2{method}.tsv` | Log2-transformed normalized expression (primary downstream input) |
| `library_size.tsv` | Per-sample library size before/after filtering |
| `filtering_summary.txt` | Gene count statistics (total, passed, removed) |
| `gene_lengths.tsv` | Gene lengths from GTF (if TPM/FPKM used) |
| `gene_filtering.png` | Bar chart (before/after gene counts) + violin plot (expression distribution) |

## Interpretation

- Check `filtering_summary.txt` for the number of genes retained vs. removed
- A reasonable pass rate is 40-70% of genes retained
- If too many genes removed: lower `min_expr` or `min_samples`
- If too few genes removed: raise `min_expr`
- The violin plot should show a clear separation between the filtered-out low-expression peak and the retained distribution

## In-Memory Data

This step caches results for downstream steps:

| Key | Description |
|-----|-------------|
| `filtered_counts` | Raw counts after filtering (genes x samples) |
| `filtered_{method}` | Normalized expression |
| `filtered_log2{method}` | Log2-transformed expression |

## Common Issues

See [Troubleshooting](rnaseq-troubleshooting.md) §3 for TPM/FPKM gene loss warnings.
