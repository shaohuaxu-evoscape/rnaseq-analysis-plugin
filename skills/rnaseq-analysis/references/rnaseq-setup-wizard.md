# Setup Wizard

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Interactive project setup that scaffolds the directory structure and generates `configs/analysis_case.yaml` from user input. Automatically triggered when no config file exists (Behavior Routing Case A).

## Pre-conditions

- `configs/analysis_case.yaml` does NOT exist
- Current directory is the intended project root

## Round 0: Project Scaffolding (automatic)

Create the directory structure and copy the config template. Do NOT ask — just execute:

```bash
mkdir -p configs inputs/ref results
cp ${CLAUDE_PLUGIN_ROOT}/configs/analysis_case.template.yaml configs/analysis_case.yaml
```

Then inform the user:

```
Project scaffolded:
  configs/analysis_case.yaml   — config file (we will fill this in together)
  inputs/                      — place your data here (gene counts, reference genome, FASTQ)
  results/                     — pipeline outputs will go here
```

## Principle

1-2 questions per round, with sensible defaults. Auto-detect what can be inferred from local files to minimize user input.

## Round 1: Project Basics

Ask: project name, organism, strain (optional).

Update `project.name`, `project.organism`, `project.strain` in the config.

## Round 2: Data Sources

**Scan scope: ONLY within the current project directory.** Do NOT use paths from external context sources (CLAUDE.md, other project directories, sibling folders, or any prior knowledge). If a file is not found within `inputs/` or `results/shared/`, it is considered absent — ask the user.

Scan `inputs/` and `results/shared/` within the project root:
- Found `.fa` / `.fasta` in `inputs/ref/` -> auto-fill `paths.reference_genome`
- Found `.gtf` / `.gff` in `inputs/ref/` -> auto-fill `paths.reference_gtf`
- Found `gene_counts*.tsv` in `inputs/` or `results/shared/` -> auto-fill `batches.batch1.gene_counts`
- Nothing found -> proceed to Round 2b (remote preprocessing)

## Round 2b: Remote Preprocessing (if no local gene_counts)

If no local gene_counts.tsv found, ask:

1. "Do you need to run preprocessing on a remote server?" (yes/no)
2. If yes, ask for:
   - Remote host (hostname or IP)
   - SSH username
   - Remote path to raw FASTQ data directory
   - Remote working directory for outputs
   - Remote path to reference genome (.fa)
   - Remote path to gene annotation (.gtf)

Update the `remote` section in config:
```yaml
remote:
  enabled: true
  host: "bioalgo-ws01"
  user: "shaohua"
  data_dir: "/data/rna/20260313"
  work_dir: "/home/shaohua/rnaseq"
  reference_genome: "/ref/A316.v1.fa"
  reference_gtf: "/ref/A316.v1.gtf"
  threads: 8
```

If the user answers "no", ask for a local gene counts file path or FASTQ directory path instead.

## Round 3: Condition Definition

If gene_counts.tsv is available, read the header to list available conditions:

```python
# Column format: {condition}-{timepoint}, e.g., R1-18, R2-48
header = open(gene_counts_path).readline().strip().split('\t')[1:]
conditions = sorted(set(col.rsplit('-', 1)[0] for col in header))
timepoints = sorted(set(int(col.rsplit('-', 1)[1]) for col in header))
```

Let the user choose which two conditions to compare. Update `target_conditions`.

## Round 4: Timepoints

Show detected timepoints, ask whether to use all or select a subset.

Update `experiment.timepoints` if a subset is chosen (otherwise leave auto-detect).

## Round 5: Normalization

Ask: normalization method (TPM default, CPM if no GTF available). Show current defaults for gene filtering.

Update `normalization.method` if changed.

## Round 6: Enrichment Analysis

Ask: NCBI Taxonomy ID and KEGG organism code. Can be skipped if unsure.

Update `project.organism_taxid`, `project.gene_id_prefix`, `differential_analysis.go_enrichment.organism_taxid`, `differential_analysis.kegg_enrichment.organism`.

## Round 7: Advanced Analysis (Optional)

Ask: whether to enable heterologous gene analysis, fermentation overview, transporter analysis. Default: all disabled.

Update the relevant `advanced_analysis.*.enabled` flags.

## Round 8: Confirm & Generate

Show a summary of all configured values. User confirms -> Write the final `configs/analysis_case.yaml`.

Then suggest next commands:

```
Configuration saved. Next steps:
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --dry-run   # preview
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b             # full run
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps                                             # list steps
```

See [Run Pipeline](rnaseq-run-pipeline.md) for execution details.
