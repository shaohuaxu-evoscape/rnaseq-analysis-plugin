#!/bin/bash
# Remote RNA-seq preprocessing — runs fastp, HISAT2, htseq-count on a remote server.
#
# Usage:
#   bash rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
#
# Reads the `remote` section from the config to determine host, paths, etc.
# Results (gene_counts.tsv) are pulled back to local results/ directory.

set -euo pipefail

# ── Parse arguments ──────────────────────────────────────────────────────────

usage() {
    echo "Usage: $0 -c <config.yaml>"
    echo ""
    echo "Run RNA-seq preprocessing on a remote server."
    echo "Reads remote.host, remote.user, remote.data_dir, etc. from the config."
    exit 1
}

CONFIG=""
while getopts "c:" opt; do
    case $opt in
        c) CONFIG="$OPTARG" ;;
        *) usage ;;
    esac
done

if [ -z "$CONFIG" ]; then
    usage
fi

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

# ── Read config values ───────────────────────────────────────────────────────

read_config() {
    python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print($1)"
}

REMOTE_ENABLED=$(read_config "c.get('remote',{}).get('enabled', False)")
if [ "$REMOTE_ENABLED" != "True" ]; then
    echo "ERROR: remote.enabled is not true in $CONFIG"
    echo "Set remote.enabled: true to use remote preprocessing."
    exit 1
fi

HOST=$(read_config "c['remote']['host']")
USER=$(read_config "c['remote']['user']")
DATA_DIR=$(read_config "c['remote']['data_dir']")
WORK_DIR=$(read_config "c['remote']['work_dir']")
REF_FA=$(read_config "c['remote']['reference_genome']")
REF_GTF=$(read_config "c['remote']['reference_gtf']")
THREADS=$(read_config "c['remote'].get('threads', 8)")
BATCH=$(read_config "list(c['batches'].keys())[0]")

echo "============================================================"
echo "Remote RNA-seq Preprocessing"
echo "============================================================"
echo "  Host:       ${USER}@${HOST}"
echo "  Batch:      ${BATCH}"
echo "  Data dir:   ${DATA_DIR}"
echo "  Work dir:   ${WORK_DIR}"
echo "  Reference:  ${REF_FA}"
echo "  GTF:        ${REF_GTF}"
echo "  Threads:    ${THREADS}"
echo "============================================================"

# ── Pre-flight checks ────────────────────────────────────────────────────────

echo ""
echo "[Pre-flight] Checking remote server..."
ssh "${USER}@${HOST}" "echo 'Connected to \$(hostname) as \$(whoami)'" || {
    echo "ERROR: Cannot connect to ${USER}@${HOST}"
    echo "Check SSH keys and network connectivity."
    exit 1
}

echo "[Pre-flight] Checking required tools..."
ssh "${USER}@${HOST}" "which fastp hisat2 samtools htseq-count" || {
    echo "ERROR: Required tools missing on remote server."
    echo "Install: fastp, hisat2, samtools, htseq-count"
    exit 1
}

echo "[Pre-flight] Checking input data..."
ssh "${USER}@${HOST}" "ls ${DATA_DIR}/*.R1.fq.gz 2>/dev/null | wc -l | xargs -I{} echo '{} FASTQ pairs found'" || {
    echo "ERROR: Cannot access ${DATA_DIR}"
    exit 1
}

# ── Execute preprocessing on remote ──────────────────────────────────────────

echo ""
echo "[Step 0a] HISAT2 Index..."

# Use a script file to avoid heredoc quoting issues
REMOTE_SCRIPT=$(cat << 'SCRIPT_TEMPLATE'
set -e

BATCH="__BATCH__"
DATA_DIR="__DATA_DIR__"
WORK_DIR="__WORK_DIR__"
REF_FA="__REF_FA__"
REF_GTF="__REF_GTF__"
THREADS="__THREADS__"

OUT="${WORK_DIR}/${BATCH}/01_preprocessing"
mkdir -p "${OUT}"/{hisat2_index,clean_fastq,fastp_reports,alignment,gene_counts/sample_counts}

# ── Step 0a: HISAT2 Index ────────────────────────────────────────────────────
INDEX_BASE=$(basename "${REF_FA}" | sed 's/\.gz$//; s/\.fa$//; s/\.fasta$//')
INDEX="${OUT}/hisat2_index/${INDEX_BASE}"

if [ -f "${INDEX}.1.ht2" ]; then
    echo "  HISAT2 index exists, skipping"
else
    echo "  Building HISAT2 index..."
    hisat2_extract_splice_sites.py "${REF_GTF}" > "${INDEX}.splice_sites.txt" 2>/dev/null || true
    hisat2_extract_exons.py "${REF_GTF}" > "${INDEX}.exons.txt" 2>/dev/null || true

    BUILD_ARGS="-p ${THREADS}"
    [ -s "${INDEX}.splice_sites.txt" ] && BUILD_ARGS="${BUILD_ARGS} --ss ${INDEX}.splice_sites.txt"
    [ -s "${INDEX}.exons.txt" ] && BUILD_ARGS="${BUILD_ARGS} --exon ${INDEX}.exons.txt"

    hisat2-build ${BUILD_ARGS} "${REF_FA}" "${INDEX}"
    echo "  HISAT2 index built"
fi

# ── Step 0b: Fastp QC/Trim ──────────────────────────────────────────────────
echo ""
echo "[Step 0b] Fastp QC/Trim..."
N_SAMPLES=0
for r1 in ${DATA_DIR}/*.R1.fq.gz; do
    [ -f "$r1" ] || continue
    SAMPLE=$(basename "$r1" .R1.fq.gz)
    R2="${DATA_DIR}/${SAMPLE}.R2.fq.gz"

    if [ -f "${OUT}/clean_fastq/${SAMPLE}.R1.clean.fq.gz" ]; then
        echo "  ${SAMPLE}: already processed, skipping"
        N_SAMPLES=$((N_SAMPLES + 1))
        continue
    fi

    echo "  ${SAMPLE}: running fastp..."
    fastp \
        -i "$r1" -I "$R2" \
        -o "${OUT}/clean_fastq/${SAMPLE}.R1.clean.fq.gz" \
        -O "${OUT}/clean_fastq/${SAMPLE}.R2.clean.fq.gz" \
        --json "${OUT}/fastp_reports/${SAMPLE}.fastp.json" \
        --html "${OUT}/fastp_reports/${SAMPLE}.fastp.html" \
        --cut_front --cut_tail --cut_window_size 4 \
        --cut_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --thread "${THREADS}" 2>/dev/null
    N_SAMPLES=$((N_SAMPLES + 1))
done
echo "  Fastp complete: ${N_SAMPLES} samples"

# ── Step 0c: HISAT2 Alignment ───────────────────────────────────────────────
echo ""
echo "[Step 0c] HISAT2 Alignment..."
for r1 in ${OUT}/clean_fastq/*.R1.clean.fq.gz; do
    [ -f "$r1" ] || continue
    SAMPLE=$(basename "$r1" .R1.clean.fq.gz)
    R2="${OUT}/clean_fastq/${SAMPLE}.R2.clean.fq.gz"
    BAM="${OUT}/alignment/${SAMPLE}.sorted.bam"

    if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
        echo "  ${SAMPLE}: already aligned, skipping"
        continue
    fi

    echo "  ${SAMPLE}: aligning..."
    hisat2 -p "${THREADS}" --dta -x "${INDEX}" \
        -1 "$r1" -2 "$R2" \
        2> "${OUT}/alignment/${SAMPLE}.hisat2.log" \
    | samtools sort -@ "${THREADS}" -o "$BAM"
    samtools index "$BAM"
done
echo "  Alignment complete"

# ── Step 0d: HTSeq Count ────────────────────────────────────────────────────
echo ""
echo "[Step 0d] HTSeq Count..."
for bam in ${OUT}/alignment/*.sorted.bam; do
    [ -f "$bam" ] || continue
    SAMPLE=$(basename "$bam" .sorted.bam)
    COUNT_FILE="${OUT}/gene_counts/sample_counts/${SAMPLE}.tsv"

    if [ -f "$COUNT_FILE" ]; then
        echo "  ${SAMPLE}: already counted, skipping"
        continue
    fi

    echo "  ${SAMPLE}: counting..."
    htseq-count \
        --format=bam --order=pos --mode=union \
        --stranded=no --type=exon --idattr=gene_id \
        "$bam" "${REF_GTF}" \
    > "$COUNT_FILE"
done

# ── Combine counts ──────────────────────────────────────────────────────────
echo ""
echo "[Combine] Merging sample counts..."
python3 -c "
import os, pandas as pd
d = '${OUT}/gene_counts/sample_counts'
frames = {}
for f in sorted(os.listdir(d)):
    if not f.endswith('.tsv'): continue
    name = f.replace('.tsv', '')
    df = pd.read_csv(os.path.join(d, f), sep='\t', header=None, names=['gene', name], index_col=0)
    df = df[~df.index.str.startswith('__')]
    frames[name] = df[name]
combined = pd.DataFrame(frames)
combined.index.name = 'gene_id'
combined.to_csv('${OUT}/gene_counts/gene_counts.tsv', sep='\t')
print(f'  Combined {len(frames)} samples, {len(combined)} genes')
"

echo ""
echo "Remote preprocessing complete."
SCRIPT_TEMPLATE
)

# Replace placeholders
REMOTE_SCRIPT="${REMOTE_SCRIPT//__BATCH__/$BATCH}"
REMOTE_SCRIPT="${REMOTE_SCRIPT//__DATA_DIR__/$DATA_DIR}"
REMOTE_SCRIPT="${REMOTE_SCRIPT//__WORK_DIR__/$WORK_DIR}"
REMOTE_SCRIPT="${REMOTE_SCRIPT//__REF_FA__/$REF_FA}"
REMOTE_SCRIPT="${REMOTE_SCRIPT//__REF_GTF__/$REF_GTF}"
REMOTE_SCRIPT="${REMOTE_SCRIPT//__THREADS__/$THREADS}"

ssh "${USER}@${HOST}" "bash -s" <<< "$REMOTE_SCRIPT"

# ── Pull results back ────────────────────────────────────────────────────────

echo ""
echo "[Pull] Fetching gene_counts.tsv..."
LOCAL_DIR="results/shared/${BATCH}/01_preprocessing/gene_counts"
mkdir -p "$LOCAL_DIR"
scp "${USER}@${HOST}:${WORK_DIR}/${BATCH}/01_preprocessing/gene_counts/gene_counts.tsv" \
    "${LOCAL_DIR}/gene_counts.tsv"

echo ""
echo "============================================================"
echo "Done. Results at: ${LOCAL_DIR}/gene_counts.tsv"
LINES=$(wc -l < "${LOCAL_DIR}/gene_counts.tsv" | tr -d ' ')
COLS=$(head -1 "${LOCAL_DIR}/gene_counts.tsv" | awk -F'\t' '{print NF-1}')
echo "  ${LINES} genes x ${COLS} samples"
echo "============================================================"
echo ""
echo "Next: python \${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c $CONFIG --steps 1a-5b"
