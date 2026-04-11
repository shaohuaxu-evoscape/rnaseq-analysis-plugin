# Config Reference

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Per-field documentation for `configs/analysis_case.yaml` — the sole user-editable config entry point.

## run_name

- **Type**: string
- **Default**: `"analysis_case_01"`
- **Purpose**: Output directory name (`results/<run_name>/`)
- **Common changes**: Use a descriptive name for each analysis scenario, e.g., `"h2l013_r1_vs_r2"`

## module1_overrides

- **Type**: dict
- **Default**: `{}`
- **Purpose**: Forwarded when a batch needs hidden preprocessing (module 1 fallback)
- **Common changes**: Usually left as empty dict

## batches

- **Type**: dict of batch_id -> batch_config
- **Required**: Yes
- **Purpose**: Define input data sources

Each batch supports two fields:
```yaml
batches:
  batch1:
    gene_counts: "path/to/gene_counts.tsv"   # Pre-computed count matrix
    raw_fastq_dir: "/path/to/raw/fastq"       # Raw FASTQ directory (preprocessing fallback)
```

**Resolution priority**:
1. `gene_counts` path exists -> use directly
2. Default module 1 output path (`results/shared/{batch_name}/01_preprocessing/gene_counts/gene_counts.tsv`) exists -> use
3. `raw_fastq_dir` exists -> automatically run hidden preprocessing (0a-0d) to generate count matrix
4. None of the above -> error

**gene_counts.tsv format**:
- TSV, first column is gene ID, subsequent columns are `{condition}-{timepoint}` (e.g., `R1-18`, `R2-48`)
- Integer count values

## target_conditions

- **Type**: dict of condition_name -> list of {batch, condition}
- **Required**: Yes, must have exactly 2 conditions
- **Purpose**: Define the two comparison groups

```yaml
target_conditions:
  condition1:                    # Internal condition name (used in output dirs, plot legends)
    - batch: batch1              # References a key in batches
      condition: R1              # Column prefix in gene_counts.tsv
  condition2:
    - batch: batch1
      condition: R2
```

- `condition1` / `condition2` are internal names used for output directories and plot legends
- `condition` must match the column prefix in gene_counts.tsv
- Cross-batch merging supported: a target_condition can include members from multiple batches

## project

- **Type**: dict
- **Purpose**: Project metadata, used in report titles and enrichment analysis

```yaml
project:
  name: "My Project"
  organism: "Species name"
  strain: "Strain ID"
  organism_taxid: 4952          # NCBI taxonomy ID (required for GO/KEGG enrichment)
  gene_id_prefix: "YALI1_"     # Gene name prefix matching count matrix gene IDs
```

### project.organism_taxid

- **Type**: int
- **Default**: 0 (must be filled in for enrichment steps)
- **Purpose**: NCBI taxonomy ID used to query UniProt for gene-to-GO and gene-to-KEGG mappings. Also used by the KEGG enrichment step if `go_enrichment.organism_taxid` is not set.
- **Examples**: 4952 (Y. lipolytica), 559292 (S. cerevisiae S288C), 284812 (S. pombe)

### project.gene_id_prefix

- **Type**: string
- **Default**: `""` (must be filled in for enrichment steps)
- **Purpose**: Prefix used to match gene names from UniProt responses against gene IDs in the count matrix. When UniProt returns multiple gene names per entry, only names starting with this prefix are used.
- **Examples**: `"YALI1_"` (Y. lipolytica A316), `"Y"` (S. cerevisiae systematic names like YAL001C)

## experiment

### experiment.timepoints

- **Type**: list of int (optional)
- **Default**: Auto-detect all shared timepoints
- **Purpose**: Specify timepoints for analysis; final set is the intersection of config values and actually available shared timepoints
- **Common changes**: Remove low-quality timepoints, e.g., `[18, 24, 48, 66]` (drop 78h)

## paths

### paths.reference_genome

- **Type**: string (relative or absolute path)
- **Default**: none
- **Purpose**: Reference genome FASTA for HISAT2 indexing (hidden preprocessing)

### paths.reference_gtf

- **Type**: string (relative or absolute path)
- **Default**: none
- **Purpose**: Gene annotation GTF for gene length calculation (TPM/FPKM) and htseq-count

## preprocessing

Parameters for hidden preprocessing (steps 0a-0d). Only used when a batch needs to generate gene_counts.tsv from raw FASTQ.

### preprocessing.threads

- **Type**: int
- **Default**: 8
- **Purpose**: Thread count for fastp/hisat2/samtools; also caps parallel htseq-count workers

### preprocessing.test_n_samples

- **Type**: int
- **Default**: 0
- **Purpose**: 0 = process all samples; N > 0 = process only the first N samples (for testing)

### preprocessing.fastp

```yaml
fastp:
  enabled: true
  cut_mean_quality: 20           # Sliding window mean quality threshold
  qualified_quality_phred: 20    # Base quality threshold
  unqualified_percent_limit: 40  # Max percentage of low-quality bases allowed
  length_required: 50            # Minimum read length
```

### preprocessing.hisat2_index / hisat2_align / htseq_count

Each contains an `enabled: true` toggle. Usually keep defaults.

## normalization

### normalization.method

- **Type**: string (`"cpm"` | `"tpm"` | `"fpkm"`)
- **Default**: `"tpm"`
- **Purpose**: Normalization method
- **Note**: TPM and FPKM require gene lengths from `paths.reference_gtf`; CPM does not

### normalization.log_transform

- **Type**: bool
- **Default**: `true`
- **Purpose**: Whether to apply log2 transformation

### normalization.pseudocount

- **Type**: number
- **Default**: 1
- **Purpose**: Pseudocount for log2 transformation: `log2(expression + pseudocount)`

### normalization.gene_filtering

```yaml
gene_filtering:
  min_expr: 1.0      # Minimum expression threshold (in CPM/TPM/FPKM units)
  min_samples: 2      # Minimum number of samples that must pass the threshold
```

Filtering logic: a gene passes if in **at least one condition**, >= `min_samples` samples have expression >= `min_expr`.

## sample_analysis

### sample_analysis.expression_summary

```yaml
expression_summary:
  enabled: true
  stats: ["mean", "median", "cv"]
```

### sample_analysis.correlation

```yaml
correlation:
  enabled: true
  method: "pearson"              # "pearson" or "spearman"
  plot_upper_triangle: true
  annot_fontsize: 8
```

### sample_analysis.pca

```yaml
pca:
  enabled: true
  n_components: 5
  plot_pairs: [[1, 2], [1, 3]]  # PC pairs to plot
```

### sample_analysis.clustering

```yaml
clustering:
  enabled: true
  top_n_genes: 2000              # Top N genes by variance
  k_range: [4, 12]               # K-means search range
  fixed_k: 7                     # Fixed K value (if not searching)
  random_state: 42
```

## differential_analysis

### differential_analysis.de_screening

```yaml
de_screening:
  enabled: true
  method: "deseq2_wald"          # DESeq2 Wald test on raw counts; "ols" for legacy fallback
  fdr_threshold: 0.05
  log2fc_threshold: 0.585        # = log2(1.5), i.e., 1.5-fold change
  fdr_method: "bh"               # Benjamini-Hochberg
```

**DESeq2 model** (auto-detected):
- Crossed design: `~ batch_id + condition + timepoint`
- Nested design: `~ condition + timepoint`

### differential_analysis.go_enrichment

```yaml
go_enrichment:
  enabled: true
  organism_taxid: 0              # NCBI taxonomy ID -- FILL THIS IN
  ontologies: ["BP", "MF", "CC"]
  sig_threshold: 0.05
  min_term_genes: 5
  max_term_genes: 500
  top_n_plot: 20
```

### differential_analysis.kegg_enrichment

```yaml
kegg_enrichment:
  enabled: true
  organism: ""                   # KEGG organism code -- FILL THIS IN (e.g., "sce", "eco")
  sig_threshold: 0.05
  min_pathway_genes: 5
  max_pathway_genes: 500
  top_n_plot: 20
```

### differential_analysis.gsea

```yaml
gsea:
  enabled: true
  ranking_metric: "mean_log2fc"  # Column used for gene ranking
  min_size: 10
  max_size: 500
  permutations: 1000
```

## advanced_analysis

### advanced_analysis.transporter_aminopeptidase

```yaml
transporter_aminopeptidase:
  enabled: true
  # Optional overrides (organism-specific):
  # amino_info: [...]           # Aminopeptidase annotation list [{gene, name, substrate, family, de, dir}]
  # transporter_overrides: {}   # Manual transporter classification overrides {gene_id: category}
  # substrate_categories: {}    # Substrate category GO regex patterns {category: regex}
  # endo_de_info: {}            # Endopeptidase annotations {gene_id: [name, family]}
```

### advanced_analysis.temporal_causality

- `enabled: true` -- Time-series cross-correlation analysis

### advanced_analysis.heterologous_genes

```yaml
heterologous_genes:
  enabled: true
  # gene_list -- Heterologous gene list, each entry has name/function/color
  # gene_list:
  #   - name: "GeneA"
  #     function: "Description"
  #     color: "#D32F2F"
  # fusion_rules -- Fusion gene rules: skip parts when fusion ID is present
  # fusion_rules:
  #   - fusion: "FusedAB"
  #     parts: ["GeneA", "GeneB"]
```

### advanced_analysis.fermentation_overview

```yaml
fermentation_overview:
  enabled: true
  data_file: "inputs/fermentation_data.xlsx"
  sheet_name: "Sheet1"
  # Optional overrides:
  # reactor_id_column: "Reactor_ID"       # Reactor ID column name in Excel
  # time_column: "Run_time (t)"           # Time column name in Excel
  # reactor_ids: ["R01", "R02"]           # Reactor IDs for the two conditions
  # metrics:                              # Column name mappings for each metric
  #   biomass_dw: "Dry_Weight (mg/g)"
  #   biomass_reactor_weight: "Reactor_Weight (Pre_Sampling,g)"
  #   product: "Total_AXT (%)"
  #   product_ester: "Esterification_ratio (%)"
  #   substrate: "Res_Sugar (RS, g/L)"
  #   feed: "Sugar_Feed_amount (g)"
```

## cluster_analysis

### cluster_analysis.cluster_deepdive

- `enabled: true` -- Per-cluster GO/KEGG enrichment and expression profile analysis

## plot

```yaml
plot:
  dpi: 300                       # Image resolution
  format: "png"                  # Output format (png/pdf/svg)
  style: "nature"                # Style preset
  colormap: "RdPu"               # Heatmap colormap
  recipe_colors:                 # Condition color mapping
    condition1: "#4477AA"
    condition2: "#EE6677"
  font_family: "Arial"
  figsize_default: [3.5, 2.8]   # Default figure size (inches)
```

**recipe_colors**: keys must match `target_conditions` names (condition1, condition2).
