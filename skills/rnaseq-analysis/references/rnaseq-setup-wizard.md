# Setup Wizard

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Interactive project setup that scaffolds the directory structure and generates `configs/analysis_case.yaml` from user input. Automatically triggered when no config file exists (Behavior Routing Case A).

## Pre-conditions

- `configs/analysis_case.yaml` does NOT exist
- Current directory is the intended project root

## Interaction Style

**Use a menu-driven format throughout.** Whenever options can be enumerated, present them as a numbered list. The user selects by entering a number. Free-text input is used only when no finite list applies (e.g., project name, file path, custom value).

**Do NOT add any section headers, round labels, or step labels before menus or questions** (e.g., do NOT write "Round 1 —", "Step 1:", "Round 2 — Data Sources", or any similar prefix in output shown to the user). Just show the menu or question directly.

After each selection, echo confirmation before presenting the next question:
```
✓ Selected: [what the user chose]
```

**Menu format** (use consistently — NO numbering prefix in the box title):
```
┌─ Title ───────────────────────────────────────────────────┐
│  [1]  Option A
│  [2]  Option B   ← recommended
│  [3]  Enter custom value
└───────────────────────────────────────────────────────────┘
Enter number:
```

Mark the recommended/default option with `← default` or `← recommended`.

## Project Scaffolding (automatic)

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

## Project Basics

**Project name (free text).** Ask separately before showing any menu:

```
请输入项目名称（自由文本，例如 "RNA Analysis Testcase"）：
```

Wait for the user's reply, then echo:
```
✓ Project name : {name}
```

**Organism (numbered menu).** Only show this menu AFTER receiving the project name. Lead with an introductory sentence:

```
请选择实验物种：

┌─ Organism ────────────────────────────────────────────────┐
│  [1]  Yarrowia lipolytica
│  [2]  Saccharomyces cerevisiae
│  [3]  Escherichia coli
│  [4]  Homo sapiens
│  [5]  Mus musculus
│  [6]  Other (enter name)
└───────────────────────────────────────────────────────────┘
输入编号：
```

**Strain (free text, optional).** Only show after organism is confirmed. Lead with an introductory sentence:

```
请输入菌株名称（选填，例如 "A316"；直接回车跳过）：
```

Echo confirmation after all three are collected:
```
✓ Project  : {name}
✓ Organism : {organism}
✓ Strain   : {strain or "—"}
```

Update `project.name`, `project.organism`, `project.strain` in the config.

## Data Sources

**Scan scope: ONLY within the current project directory.** Do NOT use paths from external context sources (CLAUDE.md, other project directories, sibling folders, or any prior knowledge). If a file is not found within `inputs/` or `results/shared/`, it is considered absent — ask the user.

### Reference Genome (local scan)

A genome lives in a subfolder under `inputs/ref/`, containing exactly one `.fa` (or `.fasta`) and one `.gtf` (or `.gff`). Scan by listing subfolders:

```bash
ls -d inputs/ref/*/  2>/dev/null
```

For each subfolder found, confirm it contains a `.fa` and a `.gtf`:
```bash
find inputs/ref/{name}/ -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \)
find inputs/ref/{name}/ -maxdepth 1 \( -name "*.gtf" -o -name "*.gff" \)
```

**Found exactly one subfolder (with valid `.fa` + `.gtf`)** → present as menu:
```
┌─ Reference Genome ────────────────────────────────────────┐
│  Found genome in inputs/ref/:
│
│  [1]  {folder_name}   ← use this
│  [2]  No, I'll use a different genome
└───────────────────────────────────────────────────────────┘
Enter number:
```

**Found multiple subfolders** → list all as numbered options, user picks one.

**No subfolders found** → skip to Gene Counts scan (genome handled during remote config).

Echo:
```
✓ Genome : {folder_name}  (or "will be selected during remote setup")
```

### Gene Counts (local scan)

Scan `inputs/` and `results/shared/` for gene counts. If found, auto-fill and skip to Condition Definition. If not found, present menu:

```
┌─ Gene Counts ─────────────────────────────────────────────┐
│  No gene_counts.tsv found in inputs/.
│
│  [1]  I already have a gene counts file (provide path)
│  [2]  I need to run preprocessing from raw FASTQ
└───────────────────────────────────────────────────────────┘
Enter number:
```

- **[1]** → ask for absolute path as free text. Set `batches.batch1.gene_counts`. Proceed to Condition Definition.
- **[2]** → proceed to Remote Preprocessing Setup.

Echo:
```
✓ Gene counts : {path or "will be generated by remote preprocessing"}
```

## Remote Preprocessing Setup

This path applies when the user has no gene counts yet and needs to run Step 0 (HISAT2 → HTSeq) on a remote server.

**Remote server connection:**

First, check whether `~/.ssh/config` exists and contains `Host` entries (excluding `Host *`):

```bash
grep '^Host ' ~/.ssh/config 2>/dev/null | grep -v '\*' | awk '{print $2}'
```

**If SSH config entries are found** — present them as a numbered menu:

```
┌─ Remote Server ───────────────────────────────────────────┐
│  以下服务器已在 ~/.ssh/config 中配置：
│
│  [1]  azure
│  [2]  bioalgo-ws01
│  [3]  {other Host entries}
│  [N]  Enter manually
└───────────────────────────────────────────────────────────┘
请选择服务器：
```

- Selecting a configured host → the SSH alias already encapsulates hostname, user, and IdentityFile; treat it as key-based auth automatically. Skip the auth method question and echo:
  ```
  ✓ Remote host : {alias}  (config from ~/.ssh/config)
  ✓ Auth method : Key-based（SSH config 中已定义）
  ```
  Set in config:
  ```yaml
  remote:
    host: "{alias}"
    ssh_auth: "key"
    ssh_key: ""    # resolved from ~/.ssh/config at runtime
    user: ""       # resolved from ~/.ssh/config at runtime
  ```

- Selecting "Enter manually" → ask for host as free text, then proceed to auth method question below.

**If no SSH config entries found** — ask for host as free text:

```
请输入远程服务器地址（SSH alias 或 hostname）：_______________
```

Then ask for the SSH authentication method:

```
┌─ SSH Authentication ──────────────────────────────────────┐
│  [1]  Key-based authentication (public/private key)   ← recommended
│  [2]  Password authentication
└───────────────────────────────────────────────────────────┘
输入编号：
```

- **[1] Key-based** — no username or password needed. Ask only for the private key path (default: `~/.ssh/id_rsa`):
  ```
  请输入私钥路径（直接回车使用默认值 ~/.ssh/id_rsa）：_______________
  ```
  Verify the file exists locally:
  ```bash
  ls {key_path}
  ```
  If not found, warn and ask to re-enter or skip. Set in config:
  ```yaml
  remote:
    ssh_auth: "key"
    ssh_key: "~/.ssh/id_rsa"
    user: ""    # auto-detected from ~/.ssh/config at runtime
  ```

- **[2] Password** — ask for SSH username as free text (password entered at runtime):
  ```
  请输入 SSH username：_______________
  ```
  Set in config:
  ```yaml
  remote:
    ssh_auth: "password"
    user: "{username}"
  ```

Echo confirmation:
```
✓ Remote host : {host}
✓ Auth method : {key-based / password}
```

**Batches, conditions, and timepoints:**

A **comparison group** (condition1 / condition2) can contain samples from multiple batches and multiple condition labels. Each comparison group is defined as a list of `{batch, condition}` pairs.

**Define batches:**

Ask how many batches are involved (free text, default 1):

```
涉及几个批次？（直接回车默认为 1）：_______________
```

For each batch (repeat until all batches are defined), ask for the batch ID (the remote data directory name) as free text:

```
Batch {n} — 批次 ID（例如 20260313）：_______________
```

Internally assign names `batch1`, `batch2`, … in order. Update `batches` in config:
```yaml
batches:
  batch1:
    gene_counts: ""    # to be filled after preprocessing
  batch2:
    gene_counts: ""
```

Echo summary:
```
✓ Batches defined:
    batch1  →  20260313
    batch2  →  20251103
```

**Define comparison groups:**

First, ask the user to name the two comparison groups so they know what they're assigning to:

```
请为两个比较组命名（这将作为分析结果中的组别标签）：

  Condition 1 名称（例如 Control、R1_recipe）：_______________
  Condition 2 名称（例如 Treatment、R2_recipe）：_______________
```

Echo:
```
✓ Condition 1 : {name1}
✓ Condition 2 : {name2}
```

Then go through each batch one by one. For each batch, ask for a single-line input using this delimiter convention:
- **Comma** (`,`) — separates multiple conditions within the same group
- **Semicolon** (`;`) — separates the two groups

Format: `{name1_conditions} ; {name2_conditions}`. A side that is left blank means no conditions from this batch for that group.

```
batch1（20260313）条件分配
  格式：{name1}的条件 ; {name2}的条件   （组内多个用逗号，留空表示该批次不参与）
  示例：R1, R2 ; R3        →  {name1}: R1, R2  /  {name2}: R3
        R1 ;               →  {name1}: R1      /  {name2}: 无
        ; T1, T2           →  {name1}: 无       /  {name2}: T1, T2

  > _______________
```

Parse by splitting on `;` first (yields two parts), then splitting each part on `,`. Strip whitespace from all tokens. If a part is empty, assign no conditions to that group from this batch.

Repeat for each remaining batch. After all batches are processed, echo the full summary:

```
✓ Condition 1（{name1}）:
    batch1 / R1
    batch1 / R2
✓ Condition 2（{name2}）:
    batch1 / R3
    batch2 / T1
```

If a condition label is entered for both groups in the same batch, warn and ask to re-enter that batch's assignments.

**Timepoints:**

Default is auto-detection: leave `experiment.timepoints` unset, and the pipeline will detect all timepoints shared between the two conditions from gene_counts.tsv column names at runtime.

Ask only if the user wants to restrict to a subset:

```
┌─ Timepoints ──────────────────────────────────────────────┐
│  [1]  Auto-detect from data   ← default
│  [2]  Specify a subset manually
└───────────────────────────────────────────────────────────┘
输入编号：
```

- **[1] Auto-detect**: leave `experiment.timepoints` commented out in config.
- **[2] Specify**: ask for timepoints as free text (comma-separated), then set `experiment.timepoints: [...]` in config.

Update `target_conditions` in config.

Show derived paths and ask to confirm:

```
Derived paths (will be set automatically):
  data_dir   = /fold/fermentation-rna-data/data/{batch_name}/
  shared_dir = /fold/home/shaohua/evoprojects/common/shared/{batch_name}/
  work_dir   = /fold/home/{user}/{project_name}/
  deploy_dir = /fold/home/shaohua/evoprojects/rnaseq-analysis-plugin  (fixed)

┌─ Confirm derived paths ───────────────────────────────────┐
│  [1]  Looks correct, continue
│  [2]  Edit manually
└───────────────────────────────────────────────────────────┘
Enter number:
```

After confirming, verify FASTQ data on the remote server via `mcp__remote-linux__Bash`:
```bash
ls /fold/fermentation-rna-data/data/{batch_name}/
```

Two FASTQ naming patterns exist across batches. Auto-detect by checking filenames:

```python
import re

fq_files = [f for f in files if f.endswith('.R1.fq.gz')]

# Try Pattern A: {condition}-{timepoint}.R1.fq.gz
pattern_a = [re.match(r'^(.+)-(\d+)\.R1\.fq\.gz$', f) for f in fq_files]
is_pattern_a = sum(1 for m in pattern_a if m) >= len(fq_files) * 0.8
```

**Pattern A detected** (`{condition}-{timepoint}.{R1|R2}.fq.gz`, e.g. `R1-18.R1.fq.gz`):

Show detection summary:
```
Found in /fold/fermentation-rna-data/data/20260313/:
  Pattern    : A  ({condition}-{timepoint}.{R1|R2}.fq.gz)
  Conditions : R1, R2, R5, R6
  Timepoints : 18, 24, 48, 66, 78
  Total      : 20 samples (paired-end)
```

Cross-check declared conditions → warn if any missing, note extras that will be skipped.

Auto-generate `inputs/sample_manifest.tsv` (sample_id = `{condition}-{timepoint}`):
```
sample_id  condition  timepoint  r1_file           r2_file
R1-18      R1         18         R1-18.R1.fq.gz    R1-18.R2.fq.gz
R1-24      R1         24         R1-24.R1.fq.gz    R1-24.R2.fq.gz
R2-18      R2         18         R2-18.R1.fq.gz    R2-18.R2.fq.gz
...
```
Filter to declared conditions only. Write manifest to `inputs/sample_manifest.tsv`.

---

**Pattern B detected** (`{sample_id}.{R1|R2}.fq.gz`, e.g. `W1.R1.fq.gz`):

Show detection summary:
```
Found in /fold/fermentation-rna-data/data/20251103/:
  Pattern    : B  ({sample_id}.{R1|R2}.fq.gz)
  Sample IDs : W1, W2, W3, Y1, Y2, Y3
```

Present experiment design as a menu:
```
┌─ Experiment Design ───────────────────────────────────────┐
│  Filenames contain only sample IDs (no condition/timepoint).
│  What is the experiment design?
│
│  [1]  Time-series — each sample is a different timepoint
│  [2]  Biological replicates — samples are replicates at the same timepoint
└───────────────────────────────────────────────────────────┘
Enter number:
```

**Sub-case B1: Time-series** — ask user to map each sample ID to condition and timepoint as free text. Recommended sample_id: `{condition}-{timepoint}`:
```
sample_id  condition  timepoint  r1_file      r2_file
W-18       W          18         W1.R1.fq.gz  W1.R2.fq.gz
W-24       W          24         W2.R1.fq.gz  W2.R2.fq.gz
W-48       W          48         W3.R1.fq.gz  W3.R2.fq.gz
Y-18       Y          18         Y1.R1.fq.gz  Y1.R2.fq.gz
...
```

**Sub-case B2: Biological replicates** — ask user to map each sample ID to condition/timepoint. Recommended sample_id: `{condition}-{timepoint}-rep{n}`:
```
sample_id   condition  timepoint  r1_file      r2_file
W-0-rep1    W          0          W1.R1.fq.gz  W1.R2.fq.gz
W-0-rep2    W          0          W2.R1.fq.gz  W2.R2.fq.gz
W-0-rep3    W          0          W3.R1.fq.gz  W3.R2.fq.gz
Y-0-rep1    Y          0          Y1.R1.fq.gz  Y1.R2.fq.gz
...
```
> Note: For replicate designs, time-series step 4b (temporal causality) is not applicable and will be disabled.

Once mapping is provided, write `inputs/sample_manifest.tsv` accordingly.

---

Set `preprocessing.fastq_pattern` in config:
```yaml
preprocessing:
  fastq_pattern: "auto"          # wizard writes "A" or "B" after detection
  sample_manifest: "inputs/sample_manifest.tsv"
```

**Reference Genome Selection and Sync:**

Three sub-cases based on the Reference Genome local scan result:

**Case 1: Local genome found in `inputs/ref/{genome_folder}/`**

Upload the genome folder to the shared reference library on the remote server:
```bash
# Create shared ref dir on remote (via mcp__remote-linux__Bash)
mkdir -p /fold/home/shaohua/evoprojects/common/ref/{genome_folder}

# Upload genome folder from local (via local Bash)
scp -r inputs/ref/{genome_folder}/ {user}@{host}:/fold/home/shaohua/evoprojects/common/ref/
```
Then find the `.fa` and `.gtf` inside the uploaded folder and set in config:
```yaml
remote:
  reference_genome: "/fold/home/shaohua/evoprojects/common/ref/{genome_folder}/{genome.fa}"
  reference_gtf:    "/fold/home/shaohua/evoprojects/common/ref/{genome_folder}/{genome.gtf}"
```

**Case 2: No local genome — use shared genome library**

List available genome folders in the shared reference library via `mcp__remote-linux__Bash`:
```bash
ls -d /fold/home/shaohua/evoprojects/common/ref/*/
```

Present as a numbered menu:
```
┌─ Reference Genome (shared library) ───────────────────────┐
│  Available genomes on remote server:
│
│  [1]  A316
│  [2]  CLIB89
│  [3]  {other folders found}
│  [4]  I have a custom path on the remote server
└───────────────────────────────────────────────────────────┘
Enter number:
```

The genome is referenced in-place from the shared library — no copy needed. Set in config:
```yaml
remote:
  reference_genome: "/fold/home/shaohua/evoprojects/common/ref/{chosen_genome}/{genome.fa}"
  reference_gtf:    "/fold/home/shaohua/evoprojects/common/ref/{chosen_genome}/{genome.gtf}"
```

**Case 3: No local genome, user wants a custom path**

Ask the user for the full remote path to the genome folder as free text:
```yaml
remote:
  reference_genome: "/custom/path/{genome_folder}/{genome.fa}"
  reference_gtf:    "/custom/path/{genome_folder}/{genome.gtf}"
```

**Write complete remote config:**

```yaml
remote:
  enabled: true
  host: "bioalgo-ws01"
  user: "alice"                                                      # logged-in user's SSH username
  batch_name: "20260313"                                             # batch ID; data_dir and shared_dir are auto-derived
  data_dir: "/fold/fermentation-rna-data/data/20260313"             # auto-derived, do not edit manually
  shared_dir: "/fold/home/shaohua/evoprojects/common/shared/20260313"  # auto-derived; Step 0 outputs shared by all users
  work_dir: "/fold/home/alice/rna_test"                             # user's analysis dir (/fold/home/{user}/{project_name})
  deploy_dir: "/fold/home/shaohua/evoprojects/rnaseq-analysis-plugin"  # shared scripts — always shaohua's folder
  reference_genome: "/fold/home/shaohua/evoprojects/common/ref/{genome}/{genome.fa}"  # resolved in Step 3
  reference_gtf:    "/fold/home/shaohua/evoprojects/common/ref/{genome}/{genome.gtf}" # resolved in Step 3
  threads: 8
```

After config is saved, remind the user:
```
Next: run Step 0 on the remote server to generate gene_counts.tsv, then re-run the pipeline locally.
  bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
```

Skip Condition Definition and Timepoints sections for now (conditions/timepoints already set above). Resume from Normalization.

## Condition Definition

*(Skip if coming from Remote Preprocessing Setup — conditions were already collected.)*

Read the gene_counts.tsv header to extract conditions:

```python
# Column format: {condition}-{timepoint}, e.g., R1-18, R2-48
header = open(gene_counts_path).readline().strip().split('\t')[1:]
conditions = sorted(set(col.rsplit('-', 1)[0] for col in header))
timepoints = sorted(set(int(col.rsplit('-', 1)[1]) for col in header))
```

Generate all pairwise combinations and present as a menu:
```
┌─ DE Comparison ───────────────────────────────────────────┐
│  Detected conditions: R1, R2, R5, R6
│  Which two conditions to compare?
│
│  [1]  R1  vs  R2   ← recommended (first two)
│  [2]  R1  vs  R5
│  [3]  R1  vs  R6
│  [4]  R2  vs  R5
│  [5]  R2  vs  R6
│  [6]  R5  vs  R6
└───────────────────────────────────────────────────────────┘
Enter number:
```

Echo:
```
✓ Comparison : {condition1} vs {condition2}
```

Update `target_conditions`.

## Timepoints

*(Skip if coming from Remote Preprocessing Setup — timepoints were already collected.)*

Present detected timepoints as a menu:

```
┌─ Timepoints ──────────────────────────────────────────────┐
│  Detected timepoints: 18, 24, 48, 66, 78
│
│  [1]  Use all: 18, 24, 48, 66, 78   ← default
│  [2]  Select a subset
└───────────────────────────────────────────────────────────┘
Enter number:
```

If **[2]** selected, list each timepoint as a numbered item and ask the user which to EXCLUDE (all are included by default):
```
┌─ Select Timepoints ───────────────────────────────────────┐
│  All timepoints are included by default.
│  Enter the numbers to EXCLUDE (comma-separated):
│
│  [1]  18h
│  [2]  24h
│  [3]  48h
│  [4]  66h
│  [5]  78h
└───────────────────────────────────────────────────────────┘
Exclude numbers (e.g., 3,5), or 0 to keep all:
```

Parse the user's reply as a comma-separated list of numbers to exclude; the remaining timepoints are kept.

Echo:
```
✓ Timepoints : {selected list}
```

Update `experiment.timepoints`.

## Normalization

```
┌─ Normalization Method ────────────────────────────────────┐
│  [1]  TPM   — Transcripts Per Million (requires GTF)   ← default
│  [2]  CPM   — Counts Per Million (no GTF needed)
│  [3]  FPKM  — Fragments Per Kilobase Per Million
└───────────────────────────────────────────────────────────┘
Enter number:
```

Then show current gene filtering defaults and ask to confirm or edit:
```
Gene filtering defaults:
  min_expr    = 1.0   (minimum normalized expression)
  min_samples = 2     (minimum samples above threshold)

┌──────────────────────────────────────────────────────────┐
│  [1]  Keep defaults   ← recommended
│  [2]  Edit values
└──────────────────────────────────────────────────────────┘
Enter number:
```

Echo:
```
✓ Normalization : {method}
✓ Gene filtering: min_expr={x}, min_samples={n}
```

Update `normalization.method` and `normalization.gene_filtering`.

## Enrichment Analysis

Present a menu of common organisms with pre-filled IDs:

```
┌─ Organism for Enrichment ─────────────────────────────────┐
│  Select organism to auto-fill taxid and KEGG code:
│
│  [1]  Yarrowia lipolytica    taxid=4952   KEGG=yli
│  [2]  Saccharomyces cerevisiae  taxid=4932  KEGG=sce
│  [3]  Escherichia coli K-12  taxid=83333  KEGG=eco
│  [4]  Homo sapiens           taxid=9606   KEGG=hsa
│  [5]  Mus musculus           taxid=10090  KEGG=mmu
│  [6]  Enter manually
│  [7]  Skip (fill in later)
└───────────────────────────────────────────────────────────┘
Enter number:
```

If **[6]**, ask for taxid and KEGG code separately as free text.

Echo:
```
✓ Organism taxid : {taxid}
✓ KEGG code      : {kegg}
✓ Gene ID prefix : {prefix}  (e.g., YALI1_ for Y. lipolytica)
```

Update `project.organism_taxid`, `project.gene_id_prefix`, `differential_analysis.go_enrichment.organism_taxid`, `differential_analysis.kegg_enrichment.organism`.

## Advanced Analysis (Optional)

Present all optional modules as a numbered list. Ask the user to enter the numbers they want to enable in a single reply. Default is all disabled (0):

```
┌─ Advanced Analysis ───────────────────────────────────────┐
│  All modules below are disabled by default.
│  Enter numbers to enable (comma-separated), or 0 to skip:
│
│  [1]  Heterologous gene expression (4c)
│  [2]  Fermentation overview (4d)
│  [3]  Transporter & aminopeptidase analysis (4a)
│  [4]  Temporal causality (4b)
└───────────────────────────────────────────────────────────┘
Enable (e.g., 1,3), or 0 to skip all:
```

Parse the user's reply as a comma-separated list of numbers to enable. If "0" or empty, keep all disabled.

For any module toggled on, ask for required parameters:
- **4c (Heterologous genes)**: ask for gene list as free text (name, function, color per gene)
- **4d (Fermentation overview)**: ask for data file path
- **4b (Temporal causality)**: ask for gene module lists (TORC1, sterol, etc.)

Echo:
```
✓ Heterologous genes  : {enabled/disabled}
✓ Fermentation data   : {enabled/disabled}
✓ Transporter analysis: {enabled/disabled}
✓ Temporal causality  : {enabled/disabled}
```

Update the relevant `advanced_analysis.*.enabled` flags.

## Confirm & Generate

Show a full summary table before writing:

```
┌─ Configuration Summary ───────────────────────────────────┐
│
│  Project
│    name       : {name}
│    organism   : {organism}  (taxid={taxid}, KEGG={kegg})
│    strain     : {strain}
│
│  Data
│    gene counts: {path}
│    genome     : {local/remote path}
│
│  Comparison
│    condition1 : {cond1}
│    condition2 : {cond2}
│    timepoints : {list}
│
│  Normalization: {method}  (min_expr={x}, min_samples={n})
│
│  Advanced     : {enabled modules or "none"}
│
│  [1]  Confirm and write configs/analysis_case.yaml   ← 
│  [2]  Go back and edit a specific section
└───────────────────────────────────────────────────────────┘
Enter number:
```

If **[2]**, ask which section to go back to:
```
┌─ Go back to ──────────────────────────────────────────────┐
│  [1]  Project Basics
│  [2]  Data Sources
│  [3]  Conditions
│  [4]  Timepoints
│  [5]  Normalization
│  [6]  Enrichment
│  [7]  Advanced Analysis
└───────────────────────────────────────────────────────────┘
Enter number:
```

On confirmation, write `configs/analysis_case.yaml` and suggest next commands:

```
Configuration saved. Next steps:
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --dry-run   # preview
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b             # full run
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps                                             # list steps
```

See [Run Pipeline](rnaseq-run-pipeline.md) for execution details.
