# Extend: Add Step

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Complete reference for adding new analysis steps to the RNA-seq pipeline.

## 1. Step File Template

New step files go in `scripts/<module_name>/sNN_<step_name>.py`.

```python
"""Step Xa: Step Display Name -- one-line description of what this step does.

Outputs:
  - result_table.tsv: Main result table
  - summary_plot.png: Summary figure
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config_data import (
    get_filtered_expr,
    get_sample_names,
    get_sample_condition,
    get_sample_display_name,
    require_sample_manifest,
)
from ..config_runtime import (
    ensure_output_dir,
    get_conditions,
    get_timepoints,
    write_readme,
)
from ..de_helpers import load_de_data, load_de_results
from ..docgen import update_section, img
from ..logger import setup_logger
from ..plotting import (
    apply_style,
    style_axes,
    save_figure,
    get_recipe_colors,
    plot_filename,
    PALETTE,
    COL1,
)

log = setup_logger("module_name.step_name")


def run(cfg):
    """Run the analysis step.

    Parameters
    ----------
    cfg : dict
        Pipeline config dict loaded by config_runtime.load_config().
    """
    # Check if enabled (if step has a toggle)
    step_cfg = cfg.get("section_name", {}).get("step_key", {})
    if not step_cfg.get("enabled", True):
        log.info("Step name disabled, skipping")
        return

    log.info("-- Step Xa: Step Name --")
    apply_style(cfg)

    # -- Get config and data --
    samples = get_sample_names(cfg)
    conditions = get_conditions(cfg)
    timepoints = get_timepoints(cfg)
    colors = get_recipe_colors(cfg)
    out_dir = ensure_output_dir(cfg, "0N_module/step_subdir")

    # Load expression data
    log2expr = get_filtered_expr(cfg, log2=True)   # log2(TPM+1), genes x samples
    expr = get_filtered_expr(cfg, log2=False)       # TPM, genes x samples

    # Load DE data (if needed)
    de_genes, background, direction = load_de_data(cfg)
    de_results = load_de_results(cfg)  # Full DE table (all genes)

    # -- Core computation --
    result_df = pd.DataFrame(...)

    # -- Save results --
    result_df.to_csv(os.path.join(out_dir, "result_table.tsv"), sep="\t")
    log.info(f"  Saved result_table.tsv ({len(result_df)} rows)")

    # -- Plot --
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))
    # ... binddata ...
    ax.set_xlabel("X Label")
    ax.set_ylabel("Y Label")
    ax.set_title("Plot Title")
    style_axes(ax)
    plot_name = plot_filename(cfg, "summary_plot.png")
    save_figure(fig, os.path.join(out_dir, "summary_plot.png"), cfg)
    log.info(f"  Saved {plot_name}")

    # -- README --
    write_readme(out_dir, "Xa", "Step Name", {
        "result_table.tsv": "Main result table containing...",
        plot_name: "Summary figure showing...",
    })

    # -- Auto report --
    subdir = "0N_module/step_subdir"
    doc = f"""\
### Parameters

- **Param1**: {step_cfg.get("param1", "default")}
- **Param2**: {step_cfg.get("param2", "default")}

### Results

- Found {len(result_df)} items

{img(cfg, subdir, "summary_plot.png", "Summary figure title")}
"""
    update_section(cfg, "Xa", "Step Name", doc)

    # -- Optional: cache to _data --
    cfg.setdefault("_data", {})["step_result"] = result_df
```

## 2. Step Registration

Register in `configs/step_registry.yaml`:

```yaml
modules:
  module_name:
    display_name: "Module Display Name"
    module_number: N
    public: true
    steps:
      - id: "Na"
        module_path: "module_name.sNN_step_name"
        display_name: "Step Display Name"
        depends_on: ["1a"]     # Steps that must complete first
```

**ID naming rules**:
- Format: `{module_number}{lowercase_letter}`, e.g., `1a`, `3b`, `4d`
- Module number matches the `module_number` field
- Letters increment in order within the module

**module_path rules**:
- Corresponds to `scripts/{module_name}/sNN_{step_name}.py`
- No `scripts.` prefix (runner adds it automatically)

**depends_on rules**:
- List all direct dependency step IDs
- Runner automatically computes transitive closure (e.g., `3b` depends on `3a`, `3a` depends on `1a`, so requesting `3b` auto-adds `1a` and `3a`)
- Omit `depends_on` for steps with no dependencies

**public field**:
- `true`: User can specify via `--steps`
- `false`: Hidden step, internal use only (e.g., preprocessing 0a-0d)

## 3. Config Access API

### 3.1 Runtime Helpers (scripts/config_runtime.py)

```python
from ..config_runtime import (
    load_config,          # load_config(path, overrides=None) -> cfg dict
    get_batch_id,         # get_batch_id(cfg) -> "batch_001"
    get_conditions,       # get_conditions(cfg) -> ["condition1", "condition2"]
    get_timepoints,       # get_timepoints(cfg) -> [18, 24, 48, 66, 78]
    get_pair_name,        # get_pair_name(cfg) -> "condition1-condition2"
    get_path,             # get_path(cfg, "gene_counts") -> absolute path
    ensure_output_dir,    # ensure_output_dir(cfg, "03_differential_analysis/go") -> create and return absolute path
    ensure_shared_dir,    # ensure_shared_dir(cfg, "01_preprocessing/fastp") -> batch shared directory
    write_readme,         # write_readme(out_dir, "3a", "DE Screening", {"file": "desc"})
)
```

### 3.2 Data Helpers (scripts/config_data.py)

```python
from ..config_data import (
    get_sample_names,         # -> {"condition1": [...], "condition2": [...], "all": [...]}
    get_all_sample_names,     # -> includes all conditions (not just target)
    get_sample_condition,     # (cfg, "sample_id") -> "condition1"
    get_sample_display_name,  # (cfg, "sample_id") -> "batch_001|CTRL|18h"
    get_sample_manifest,      # -> DataFrame or None
    require_sample_manifest,  # -> DataFrame (raises if missing)
    get_filtered_expr,        # (cfg, log2=True) -> DataFrame (genes x samples)
    get_filtered_counts,      # (cfg) -> DataFrame (raw integer counts, genes x samples)
    get_norm_method,          # -> "tpm" / "cpm" / "fpkm"
    get_preprocessing_raw_dir,  # -> Path (raw FASTQ directory)
    get_module1_source_conditions,  # -> ["R1", "R2"] or None
)
```

`get_sample_names(cfg)` return structure:
```python
{
    "condition1": ["condition1_batch_001_CTRL-18", "condition1_batch_001_CTRL-24", ...],
    "condition2": ["condition2_batch_001_TREAT-18", "condition2_batch_001_TREAT-24", ...],
    "all": [... concatenation of the above ...]
}
```

`get_filtered_expr(cfg, log2=True)` return:
- Rows: gene IDs
- Columns: sample IDs
- Values: log2(TPM + 1) or raw TPM (depends on log2 param)
- Source: step 1a output, or from `_data` cache

### 3.3 DE Helpers (scripts/de_helpers.py)

```python
from ..de_helpers import (
    load_de_data,      # (cfg) -> (de_genes: list, background: list, direction: dict)
    load_de_results,   # (cfg) -> DataFrame (all genes) or None
)
```

`load_de_data(cfg)` returns:
- `de_genes`: list of significant DE gene IDs
- `background`: list of all filtered gene IDs (for enrichment background)
- `direction`: `{"Gene001": "up", "Gene002": "down", ...}`

### 3.4 Statistics Helpers (scripts/stats.py)

```python
from ..stats import bh_fdr  # bh_fdr(pvalues) -> numpy array of FDR values
```

### 3.5 Network Helpers (scripts/net.py)

```python
from ..net import (
    fetch_with_retry,  # (url, timeout=120, retries=3) -> text
    fetch_cached,      # (url, cache_dir, cache_file) -> text (auto-cached locally)
)
```

## 4. Plotting Conventions

### 4.1 Standard Flow

```python
from ..plotting import apply_style, style_axes, save_figure, plot_filename, get_recipe_colors, PALETTE, COL1

# 1. Call at start of run()
apply_style(cfg)

# 2. Create figure
fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))

# 3. Style axes after bindingdata
style_axes(ax)

# 4. Save (auto-handles format and DPI)
save_figure(fig, os.path.join(out_dir, "plot.png"), cfg)
```

### 4.2 Colors

```python
# Condition colors (from config, with defaults)
colors = get_recipe_colors(cfg)
colors["condition1"]  # "#4477AA" (blue)
colors["condition2"]  # "#EE6677" (red)

# Palette
PALETTE["blue"]    # "#2166AC"
PALETTE["red"]     # "#B2182B"
PALETTE["green"]   # "#4CAF50"
PALETTE["yellow"]  # "#CCBB44"
PALETTE["cyan"]    # "#66CCEE"
PALETTE["purple"]  # "#AA3377"
PALETTE["grey"]    # "#BBBBBB"
PALETTE["dark"]    # "#332288"
PALETTE["light_grey"]  # "#DDDDDD"
```

### 4.3 Column Width Constants

```python
COL1 = 7       # Single panel (inches)
COL1_5 = 10    # 1.5 column
COL2 = 14      # Double panel
```

### 4.4 Style Details

- Font: Arial 10pt (titles 11pt)
- All four spines visible, linewidth 0.8
- Ticks outward, length 3.5
- DPI: from `cfg["plot"]["dpi"]` (default 200)
- Format: from `cfg["plot"]["format"]` (default png)
- Background: white
- Grid: off
- Legend: no border

### 4.5 Filename Handling

```python
# Get filename with correct extension
plot_filename(cfg, "volcano_plot.png")  # -> "volcano_plot.png" or "volcano_plot.pdf"

# Get full output path
from ..plotting import plot_path
plot_path(cfg, out_dir, "volcano_plot.png")
```

## 5. Documentation Conventions

### 5.1 write_readme

Generates README.txt in each step's output directory:

```python
write_readme(out_dir, "3a", "DE Screening", {
    "de_results_all.tsv": "Full DE statistics for all tested genes",
    "de_genes.tsv": "Significant DE genes only",
    "volcano_plot.png": "Volcano plot of mean log2FC vs -log10(FDR)",
})
```

### 5.2 update_section

Updates the auto-generated Markdown analysis report:

```python
from ..docgen import update_section, img

content = f"""\
### Parameters
- **Method**: DESeq2 Wald test
- **FDR threshold**: {fdr_thresh}

### Results
- DE genes: {n_de}

{img(cfg, "03_differential_analysis/de_screening", "volcano_plot.png", "Volcano plot")}
"""
update_section(cfg, "3a", "DE Screening", content)
```

- `update_section(cfg, step_id, step_name, markdown_content)`
  - Auto-inserts/replaces the step's section in the report
  - Sections delimited by HTML comments: `<!-- SECTION:3a:START -->` ... `<!-- SECTION:3a:END -->`
  - Auto-adds heading and timestamp

- `img(cfg, subdir, filename, caption)` -- generates relative image reference to comparison result directory
- `shared_img(cfg, subdir, filename, caption)` -- generates relative image reference to shared result directory

## 6. In-Memory Data Sharing

Steps share intermediate results via `cfg["_data"]` to avoid redundant disk I/O:

```python
# Write (in producer step)
cfg.setdefault("_data", {})["my_result"] = result_df

# Read (in consumer step)
cached = cfg.get("_data", {}).get("my_result")
if cached is None:
    # Fall back to disk
    cached = pd.read_csv(...)
```

**Established key conventions**:

| Key | Producer | DataFrame shape |
|-----|----------|-----------------|
| `raw_counts` | prepare | genes x samples (integer) |
| `filtered_counts` | 1a | genes x samples (integer, filtered) |
| `filtered_{method}` | 1a | genes x samples (normalized) |
| `filtered_log2{method}` | 1a | genes x samples (log2 transformed) |
| `de_results` | 3a | genes x [mean_log2fc, fdr, direction, ...] |
| `de_genes` | 3a | DE genes x [same columns] |

Where `{method}` = `tpm`, `cpm`, or `fpkm`.

## 7. Testing

```bash
# Run all tests
python -m pytest tests/

# Test step parsing (no data required)
python -m pytest tests/test_runner_steps.py -v

# Verify step registration
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps

# Dry run to verify dependencies
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps Xa --dry-run
```

Existing tests:
- `tests/test_runner_steps.py`: step parsing, ranges, dependency closure, deduplication
- `tests/test_analysis_case_preprocessing.py`: analysis case preparation and preprocessing integration
