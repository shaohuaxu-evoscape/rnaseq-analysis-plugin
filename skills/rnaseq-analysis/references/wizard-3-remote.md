# Setup Wizard — Stage 3: Remote Preprocessing Setup

This stage applies when the user has no gene counts yet and needs to run Step 0 (HISAT2 → HTSeq) on a remote server.

## Remote Server Connection

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

- Configured host selected → key-based auth assumed (SSH alias encapsulates all). Set `remote.host`, `remote.ssh_auth: "key"`, `remote.ssh_key: ""`, `remote.user: ""`. Echo `✓ Remote host / Auth method`.
- "Enter manually" → ask for host as free text, then show auth menu below.

**If no SSH config entries found** — ask for host as free text, then show auth menu:

```
┌─ SSH Authentication ──────────────────────────────────────┐
│  [1]  Key-based authentication (public/private key)   ← recommended
│  [2]  Password authentication
└───────────────────────────────────────────────────────────┘
输入编号：
```

- **[1] Key-based** — ask for private key path (default: `~/.ssh/id_rsa`). Verify with `ls {key_path}`; warn if not found. Set `remote.ssh_auth: "key"`, `remote.ssh_key: "{path}"`.
- **[2] Password** — ask for SSH username. Set `remote.ssh_auth: "password"`, `remote.user: "{username}"`.

Echo `✓ Remote host : {host}` / `✓ Auth method : {method}`.

## Batches, Conditions, and Timepoints

**Batches:**

```
┌─ Batch Count ─────────────────────────────────────────────┐
│  涉及几个批次？
│  [1]  1 个批次   ← default  [2]  2 个批次
│  [3]  3 个批次   [4]  更多（手动输入数字）
└───────────────────────────────────────────────────────────┘
输入编号：
```

For each batch ask for its ID (remote data directory name): `Batch {n} — 批次 ID（例如 20260313）：`. Assign internally as `batch1`, `batch2`, …. Echo `✓ Batches defined: batch1 → {id} …`.

**Comparison groups:**

```
请为两个比较组命名（用分号隔开，例如 Control ; Treatment）：_______________
```

For each batch, ask condition assignment on one line (`；` separates groups, `,` separates conditions within a group):

```
batch1（20260313）条件分配
  格式：{name1}的条件 ; {name2}的条件   （组内多个用逗号，留空表示该批次不参与）
  示例：R1, R2 ; R3  |  R1 ;  |  ; T1, T2

  > _______________
```

If a label appears in both groups for the same batch, warn and ask to re-enter. After all batches, echo full summary:
```
✓ Condition 1（{name1}）: batch1/R1, batch1/R2
✓ Condition 2（{name2}）: batch1/R3, batch2/T1
```

**Timepoints:**

```
┌─ Timepoints ──────────────────────────────────────────────┐
│  [1]  Auto-detect from data   ← default
│  [2]  Specify a subset manually
└───────────────────────────────────────────────────────────┘
输入编号：
```

- **[1]**: leave `experiment.timepoints` commented out.
- **[2]**: ask for timepoints (comma-separated); set `experiment.timepoints: [...]`. After gene_counts.tsv is available, validate and silently drop any not found.

Update `target_conditions` in config.

## Derived Paths

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

## FASTQ Pattern Detection

Verify FASTQ data on remote via `mcp__remote-linux__Bash`:
```bash
ls /fold/fermentation-rna-data/data/{batch_name}/
```

Auto-detect pattern (≥80% match):
```python
import re
fq_files = [f for f in files if f.endswith('.R1.fq.gz')]
is_pattern_a = sum(1 for f in fq_files if re.match(r'^(.+)-(\d+)\.R1\.fq\.gz$', f)) >= len(fq_files) * 0.8
```

**Pattern A** (`{condition}-{timepoint}.R1.fq.gz`): show detected conditions/timepoints/total, cross-check declared conditions (warn missing, skip extras). Auto-generate manifest (sample_id = `{condition}-{timepoint}`):
```
sample_id  condition  timepoint  r1_file           r2_file
R1-18      R1         18         R1-18.R1.fq.gz    R1-18.R2.fq.gz
```
Filter to declared conditions. Write to `inputs/sample_manifest.tsv`.

**Pattern B** (`{sample_id}.R1.fq.gz`): show detected sample IDs, then ask design:

```
┌─ Experiment Design ───────────────────────────────────────┐
│  [1]  Time-series — each sample is a different timepoint
│  [2]  Biological replicates — replicates at the same timepoint
└───────────────────────────────────────────────────────────┘
Enter number:
```

- **B1**: user maps sample IDs → condition + timepoint. Recommended sample_id: `{condition}-{timepoint}`.
- **B2**: user maps sample IDs → condition + timepoint + rep. Recommended: `{condition}-{timepoint}-rep{n}`. (Step 4b temporal causality disabled for replicate designs.)

Write `inputs/sample_manifest.tsv`. Set `preprocessing.fastq_pattern: "A"` (or `"B"`) and `sample_manifest: "inputs/sample_manifest.tsv"` in config.

## Reference Genome

**Case 1: Local genome in `inputs/ref/{genome_folder}/`** — upload to remote shared library:
```bash
# remote (via mcp__remote-linux__Bash)
mkdir -p /fold/home/shaohua/evoprojects/common/ref/{genome_folder}
# local
scp -r inputs/ref/{genome_folder}/ {user}@{host}:/fold/home/shaohua/evoprojects/common/ref/
```
Then set `remote.reference_genome` / `remote.reference_gtf` to the uploaded paths (see final config block below).

**Case 2: No local genome** — list shared library via `mcp__remote-linux__Bash` (`ls -d /fold/home/shaohua/evoprojects/common/ref/*/`), present as numbered menu. Referenced in-place, no copy needed.

**Case 3: Custom path** — ask user for full remote genome path as free text.

## Write Remote Config

```yaml
remote:
  enabled: true
  host: "bioalgo-ws01"
  user: "alice"
  batch_name: "20260313"
  data_dir: "/fold/fermentation-rna-data/data/20260313"
  shared_dir: "/fold/home/shaohua/evoprojects/common/shared/20260313"
  work_dir: "/fold/home/alice/rna_test"
  deploy_dir: "/fold/home/shaohua/evoprojects/rnaseq-analysis-plugin"
  reference_genome: "/fold/home/shaohua/evoprojects/common/ref/{genome}/{genome.fa}"
  reference_gtf:    "/fold/home/shaohua/evoprojects/common/ref/{genome}/{genome.gtf}"
  threads: 8
```

Remind the user:
```
Next: run Step 0 on the remote server to generate gene_counts.tsv, then re-run the pipeline locally.
  bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml

The script will pull gene_counts.tsv to:
  results/shared/{batch_name}/01_preprocessing/gene_counts/gene_counts.tsv

Do NOT set batches.batch1.gene_counts in the config — leave it empty so the pipeline
auto-detects from results/shared/. Do NOT copy gene_counts to inputs/.
```

## Post-Pull Prompt

After gene_counts.tsv is confirmed pulled to `results/shared/`, immediately ask:

```
✓ gene_counts.tsv 已拉回本地：results/shared/{batch_name}/01_preprocessing/gene_counts/gene_counts.tsv
  {N} 基因 x {M} 样本

┌─ 下一步 ──────────────────────────────────────────────────┐
│  [1]  立即开始下游分析（继续 wizard-4-analysis.md）
│  [2]  稍后再运行（退出向导）
└───────────────────────────────────────────────────────────┘
输入编号：
```

- **[1]**: Read `references/wizard-4-analysis.md` and continue with Normalization.
- **[2]**: Exit wizard. Remind the user to run:
  ```
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b
  ```

## Next

Read `references/wizard-4-analysis.md` and continue with Normalization.
(Condition Definition and Timepoints are already collected above — skip them.)
