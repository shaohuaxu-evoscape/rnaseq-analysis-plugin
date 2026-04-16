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

### 2a — Reference Genome (local scan)

A genome lives in a subfolder under `inputs/ref/`, containing exactly one `.fa` (or `.fasta`) and one `.gtf` (or `.gff`). Scan by listing subfolders:

```bash
ls -d inputs/ref/*/  2>/dev/null
```

For each subfolder found, confirm it contains a `.fa` and a `.gtf`:
```bash
find inputs/ref/{name}/ -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \)
find inputs/ref/{name}/ -maxdepth 1 \( -name "*.gtf" -o -name "*.gff" \)
```

- **Found exactly one subfolder (with valid `.fa` + `.gtf`)** → confirm with user: *"Found genome: `{folder_name}`. Use this one?"*
  - Yes → record as local genome; auto-fill `paths.reference_genome` and `paths.reference_gtf` to the files inside
  - No → treat as not found (proceed to shared genome selection in Round 2b)
- **Found multiple subfolders** → list them by folder name, ask user to pick one
- **No subfolders found** → record as "no local genome" (handled in Round 2b Step 3)

### 2b — Gene Counts (local scan)

Scan `inputs/` and `results/shared/` for gene counts:

- Found `gene_counts*.tsv` → auto-fill `batches.batch1.gene_counts`, **skip to Round 3**
- Not found → ask the user:

  > **Option A — I already have a gene counts file**
  > Please provide the absolute path to the file.
  > (Set `batches.batch1.gene_counts`, proceed to Round 3)

  > **Option B — I need to run preprocessing from raw FASTQ**
  > (Proceed to Round 2b: Remote Preprocessing Setup)

## Round 2b: Remote Preprocessing Setup

This path applies when the user has no gene counts yet and needs to run Step 0 (HISAT2 → HTSeq) on a remote server.

**Step 1 — Collect experiment design** (conditions and timepoints are not yet in a file header, so ask the user):

Ask:
1. What are the sample condition names? *(e.g., R1, R2, R3)*
2. What timepoints are included? *(e.g., 18, 24, 48, 66, 78)*
3. Which two conditions will be compared for DE analysis? *(e.g., R1 vs R2)*

Update `target_conditions` and `experiment.timepoints` with these values now — they will be confirmed/corrected after Step 0 produces the gene_counts header.

**Step 2 — Collect basic remote server configuration:**

Ask for:
- Remote host (hostname or IP)
- SSH username (the logged-in user — NOT necessarily "shaohua")
- Remote path to raw FASTQ data directory

Do NOT ask for `deploy_dir`, `work_dir`, or genome paths — handled below:
- `deploy_dir` is always `/home/shaohua/evoprojects/rnaseq-analysis-plugin` (shared admin scripts, fixed)
- `work_dir` is always `/home/{user}/{project_name}` (derived from SSH username and project name)
- Genome paths are resolved in Step 3 below

**Step 3 — Reference Genome Selection and Sync:**

Three sub-cases based on the Round 2a scan result:

**Case 1: Local genome found in `inputs/ref/{genome_folder}/`**

Upload the entire genome folder to the remote server:
```bash
# Create ref dir on remote (via mcp__remote-linux__Bash)
mkdir -p {work_dir}/ref

# Upload genome folder from local (via local Bash)
scp -r inputs/ref/{genome_folder}/ {user}@{host}:{work_dir}/ref/
```
Then find the `.fa` and `.gtf` inside the uploaded folder and set in config:
```yaml
remote:
  reference_genome: "{work_dir}/ref/{genome_folder}/{genome.fa}"
  reference_gtf:    "{work_dir}/ref/{genome_folder}/{genome.gtf}"
```

**Case 2: No local genome — use shared genome library**

List available genome folders on the remote server via `mcp__remote-linux__Bash`:
```bash
ls -d /home/shaohua/evoprojects/ref/*/
```

Each subfolder is one genome (contains one `.fa` + one `.gtf`). Present the folder names to the user and ask them to choose. Then copy the chosen folder to the user's directory on the remote server:
```bash
# Run on remote server via mcp__remote-linux__Bash
cp -r /home/shaohua/evoprojects/ref/{chosen_genome}/ {work_dir}/ref/
```
Find the `.fa` and `.gtf` inside the copied folder and set in config:
```yaml
remote:
  reference_genome: "{work_dir}/ref/{chosen_genome}/{genome.fa}"
  reference_gtf:    "{work_dir}/ref/{chosen_genome}/{genome.gtf}"
```

**Case 3: No local genome, user wants a custom path**

Ask the user for the full remote path to the genome folder. These files are already on the remote server and do not need to be copied.
```yaml
remote:
  reference_genome: "/custom/path/{genome_folder}/{genome.fa}"
  reference_gtf:    "/custom/path/{genome_folder}/{genome.gtf}"
```

**Step 4 — Write complete remote config:**

```yaml
remote:
  enabled: true
  host: "bioalgo-ws01"
  user: "alice"                                                      # logged-in user's SSH username
  data_dir: "/data/rna/20260313"                                     # shared raw FASTQ location
  work_dir: "/home/alice/rna_test"                                   # user's output dir (/home/{user}/{project_name})
  deploy_dir: "/home/shaohua/evoprojects/rnaseq-analysis-plugin"    # shared scripts — always shaohua's folder
  reference_genome: "{work_dir}/ref/{genome.fa}"                     # resolved in Step 3
  reference_gtf:    "{work_dir}/ref/{genome.gtf}"                    # resolved in Step 3
  threads: 8
```

After config is saved, remind the user:
```
Next: run Step 0 on the remote server to generate gene_counts.tsv, then re-run the pipeline locally.
  bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
```

Skip Rounds 3–4 for now (conditions/timepoints already set above). Resume from Round 5.

## Round 3: Condition Definition

*(Skip this round if coming from Round 2b — conditions were already collected.)*

Read the gene_counts.tsv header to list available conditions:

```python
# Column format: {condition}-{timepoint}, e.g., R1-18, R2-48
header = open(gene_counts_path).readline().strip().split('\t')[1:]
conditions = sorted(set(col.rsplit('-', 1)[0] for col in header))
timepoints = sorted(set(int(col.rsplit('-', 1)[1]) for col in header))
```

Show detected conditions and let the user choose which two to compare. Update `target_conditions`.

## Round 4: Timepoints

*(Skip this round if coming from Round 2b — timepoints were already collected.)*

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
