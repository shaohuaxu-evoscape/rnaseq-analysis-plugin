# Remote Preprocessing (Interactive / Debug Mode)

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Step-by-step interactive remote preprocessing via MCP remote-linux. Use this as a **fallback** when the one-shot SSH script fails or when you need to debug individual steps.

**Primary method**: Use the SSH script first:
```bash
bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
```

**This document is for**: Interactive debugging, partial re-runs, or when MCP remote-linux is preferred over the SSH script.

## When to Use

- The SSH script (`rnaseq_remote_preprocess.sh`) failed and you need to diagnose which step broke
- User wants to run a single preprocessing step (e.g., only re-run fastp)
- MCP remote-linux is available and user prefers interactive control
- `remote.enabled: true` is set in the config

## Configuration

Read the `remote` section from `configs/analysis_case.yaml`:

```yaml
remote:
  enabled: true
  host: "bioalgo-ws01"                                # Remote hostname
  user: "alice"                                       # Logged-in user's SSH username
  data_dir: "/data/rna/20260313"                      # Shared raw FASTQ location
  work_dir: "/home/alice/rna_test"                    # User's output directory (/home/{user}/{project})
  deploy_dir: "/home/shaohua/evoprojects/rnaseq-analysis-plugin"  # Shared pipeline scripts (admin, fixed)
  reference_genome: "/ref/A316.v1.fa"
  reference_gtf: "/ref/A316.v1.gtf"
  threads: 8
```

**Key distinction:**
- `deploy_dir` — always `/home/shaohua/evoprojects/rnaseq-analysis-plugin/`. Scripts and step implementations live here. Never changes between users.
- `work_dir` — always `/home/{user}/{project_name}/`. All Step 0 outputs (BAM, counts, QC) land here. Each user gets their own directory.

## Pre-flight Checks

Before running any preprocessing, verify the remote environment. Use `mcp__remote-linux__Bash` for all remote commands.

### 1. Check connectivity

```bash
echo "Connected to $(hostname) as $(whoami)"
```

### 2. Check required tools

```bash
which fastp hisat2 samtools htseq-count
```

If any tool is missing, stop and inform the user.

### 3. Check input data

```bash
# Verify reference files exist
ls -la {reference_genome} {reference_gtf}

# List available FASTQ files
ls {data_dir}/*.R1.fq.gz 2>/dev/null | head -20
```

### 4. Check shared scripts

```bash
ls {deploy_dir}/scripts/rnaseq_run.py
```

If not found, the admin needs to deploy first. See [Remote Deployment](rnaseq-remote-deployment.md).

### 5. Create working directory

```bash
mkdir -p {work_dir}/01_preprocessing/{hisat2_index,clean_fastq,fastp_reports,alignment,gene_counts/sample_counts}
```

## Step 0a: HISAT2 Index

Build the genome index. Skip if index files already exist.

```bash
INDEX_PREFIX="{work_dir}/01_preprocessing/hisat2_index/$(basename {reference_genome} .fa)"

# Check if index already exists
if [ -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "HISAT2 index already exists, skipping"
else
    # Extract splice sites and exons
    hisat2_extract_splice_sites.py {reference_gtf} > ${INDEX_PREFIX}.splice_sites.txt
    hisat2_extract_exons.py {reference_gtf} > ${INDEX_PREFIX}.exons.txt

    # Build index
    hisat2-build -p {threads} \
        --ss ${INDEX_PREFIX}.splice_sites.txt \
        --exon ${INDEX_PREFIX}.exons.txt \
        {reference_genome} ${INDEX_PREFIX}
fi
```

## Step 0b: Fastp QC/Trim

Run fastp for each sample. Discover samples from the data directory.

```bash
# Discover sample IDs
for r1 in {data_dir}/*.R1.fq.gz; do
    SAMPLE=$(basename "$r1" .R1.fq.gz)
    R2="{data_dir}/${SAMPLE}.R2.fq.gz"
    OUT_DIR="{work_dir}/01_preprocessing"

    echo "Processing $SAMPLE ..."
    fastp \
        -i "$r1" -I "$R2" \
        -o "${OUT_DIR}/clean_fastq/${SAMPLE}.R1.clean.fq.gz" \
        -O "${OUT_DIR}/clean_fastq/${SAMPLE}.R2.clean.fq.gz" \
        --json "${OUT_DIR}/fastp_reports/${SAMPLE}.fastp.json" \
        --html "${OUT_DIR}/fastp_reports/${SAMPLE}.fastp.html" \
        --cut_front --cut_tail --cut_window_size 4 \
        --cut_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --thread {threads}
done
```

For large numbers of samples, process in batches and report progress.

## Step 0c: HISAT2 Alignment

Align each sample to the genome index.

```bash
INDEX_PREFIX="{work_dir}/01_preprocessing/hisat2_index/$(basename {reference_genome} .fa)"
OUT_DIR="{work_dir}/01_preprocessing"

for r1 in ${OUT_DIR}/clean_fastq/*.R1.clean.fq.gz; do
    SAMPLE=$(basename "$r1" .R1.clean.fq.gz)
    R2="${OUT_DIR}/clean_fastq/${SAMPLE}.R2.clean.fq.gz"

    echo "Aligning $SAMPLE ..."
    hisat2 -p {threads} --dta -x ${INDEX_PREFIX} \
        -1 "$r1" -2 "$R2" \
        2> "${OUT_DIR}/alignment/${SAMPLE}.hisat2.log" \
    | samtools sort -@ {threads} -o "${OUT_DIR}/alignment/${SAMPLE}.sorted.bam"

    samtools index "${OUT_DIR}/alignment/${SAMPLE}.sorted.bam"
done
```

## Step 0d: HTSeq Count

Quantify gene expression for each sample, then combine into a matrix.

```bash
OUT_DIR="{work_dir}/01_preprocessing"
GTF="{reference_gtf}"

# Count per sample
for bam in ${OUT_DIR}/alignment/*.sorted.bam; do
    SAMPLE=$(basename "$bam" .sorted.bam)
    echo "Counting $SAMPLE ..."
    htseq-count \
        --format=bam --order=pos --mode=union \
        --stranded=no --type=exon --idattr=gene_id \
        "$bam" "$GTF" \
    > "${OUT_DIR}/gene_counts/sample_counts/${SAMPLE}.tsv"
done

# Combine into matrix
echo "Combining counts..."
python3 -c "
import os, pandas as pd
counts_dir = '${OUT_DIR}/gene_counts/sample_counts'
frames = {}
for f in sorted(os.listdir(counts_dir)):
    if not f.endswith('.tsv'): continue
    name = f.replace('.tsv', '')
    df = pd.read_csv(os.path.join(counts_dir, f), sep='\t', header=None, names=['gene', name], index_col=0)
    df = df[~df.index.str.startswith('__')]
    frames[name] = df[name]
combined = pd.DataFrame(frames)
combined.index.name = 'gene_id'
combined.to_csv('${OUT_DIR}/gene_counts/gene_counts.tsv', sep='\t')
print(f'Combined {len(frames)} samples, {len(combined)} genes')
"
```

## Pull Results Back to Local

After all steps complete, use **local Bash** (not remote) to pull the gene_counts.tsv back:

```bash
# Create local directory
mkdir -p inputs/gene_counts

# Pull gene counts (work_dir is the user's project dir on the remote server)
scp {user}@{host}:{work_dir}/01_preprocessing/gene_counts/gene_counts.tsv \
    inputs/gene_counts/

# Optionally pull fastp reports for QC review
mkdir -p inputs/fastp_reports
scp {user}@{host}:{work_dir}/01_preprocessing/fastp_reports/*.json \
    inputs/fastp_reports/
```

Update `batches.batch1.gene_counts` in `configs/analysis_case.yaml` to point to the pulled file:
```yaml
batches:
  batch1:
    gene_counts: "inputs/gene_counts/gene_counts.tsv"
```

## Post-Preprocessing

Verify the gene_counts.tsv:

```bash
head -2 inputs/gene_counts/gene_counts.tsv
wc -l inputs/gene_counts/gene_counts.tsv
```

Then proceed with local analysis:

```
python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b
```

## Error Handling

- **Connection failed**: Check MCP remote-linux configuration in Claude Code settings
- **Tool not found**: Ask user to install the missing tool on the remote server
- **FASTQ not found**: Verify `remote.data_dir` path and file naming convention (`*.R1.fq.gz`)
- **Disk space**: Check `df -h {work_dir}` before starting — preprocessing generates large intermediate files
- **Partial failure**: Each step can be re-run independently. Check which samples completed and resume from there.
