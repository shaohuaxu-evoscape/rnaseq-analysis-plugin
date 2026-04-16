#!/bin/bash
# Remote preprocessing — SSH to remote server and run Module 0 via the deployed Python scripts.
#
# Usage:
#   bash rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
#
# Prerequisites:
#   - Plugin deployed on remote via rnaseq_remote_setup.sh
#   - remote.enabled: true in config
#   - SSH key access to remote host

set -euo pipefail

usage() {
    echo "Usage: $0 -c <config.yaml>"
    echo ""
    echo "Run RNA-seq preprocessing on a remote server using the deployed plugin."
    echo "Reads remote.host, remote.user, remote.deploy_dir from the config."
    exit 1
}

CONFIG=""
while getopts "c:" opt; do
    case $opt in
        c) CONFIG="$OPTARG" ;;
        *) usage ;;
    esac
done

if [ -z "$CONFIG" ] || [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: ${CONFIG:-<not specified>}"
    usage
fi

# ── Read config ──────────────────────────────────────────────────────────────

read_config() {
    python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print($1)"
}

REMOTE_ENABLED=$(read_config "c.get('remote',{}).get('enabled', False)")
if [ "$REMOTE_ENABLED" != "True" ]; then
    echo "ERROR: remote.enabled is not true in $CONFIG"
    exit 1
fi

HOST=$(read_config "c['remote']['host']")
USER=$(read_config "c['remote']['user']")
DEPLOY_DIR=$(read_config "c['remote'].get('deploy_dir', '')")
CONDA_ENV=$(read_config "c['remote'].get('conda_env', 'rnaseq')")
CONDA_PREFIX=$(read_config "c['remote'].get('conda_prefix', '')")
BATCH=$(read_config "list(c['batches'].keys())[0]")

# Derive Python path: use conda_prefix if set, else fall back to ~/miniconda3
if [ -n "$CONDA_PREFIX" ]; then
    REMOTE_PYTHON="${CONDA_PREFIX}/envs/${CONDA_ENV}/bin/python"
else
    REMOTE_PYTHON="\$HOME/miniconda3/envs/${CONDA_ENV}/bin/python"
fi

if [ -z "$DEPLOY_DIR" ]; then
    echo "ERROR: remote.deploy_dir not set in $CONFIG"
    echo "Set it to the plugin install path on the remote server (e.g., /fold/home/shaohua/evoprojects/rnaseq-analysis-plugin)"
    exit 1
fi

REMOTE_WORKDIR="~/rnaseq-projects/${BATCH}"
REMOTE_LOG="${REMOTE_WORKDIR}/preprocess.log"
REMOTE_DONE="${REMOTE_WORKDIR}/preprocess.done"
REMOTE_FAIL="${REMOTE_WORKDIR}/preprocess.failed"

echo "============================================================"
echo "Remote RNA-seq Preprocessing"
echo "============================================================"
echo "  Host:       ${USER}@${HOST}"
echo "  Deploy dir: ${DEPLOY_DIR}"
echo "  Conda env:  ${CONDA_ENV} (${REMOTE_PYTHON})"
echo "  Batch:      ${BATCH}"
echo "============================================================"

# ── Step 1: Sync config to remote ────────────────────────────────────────────

echo ""
echo "[1/4] Syncing config to remote..."
ssh "${HOST}" "mkdir -p ${REMOTE_WORKDIR}/configs && rm -f ${REMOTE_DONE} ${REMOTE_FAIL}"
scp "$CONFIG" "${HOST}:${REMOTE_WORKDIR}/configs/analysis_case.yaml"
echo "  Config uploaded to ${REMOTE_WORKDIR}/configs/"

# ── Step 2: Launch preprocessing in background via nohup ─────────────────────

echo ""
echo "[2/4] Launching preprocessing on remote (nohup background)..."
ssh "${HOST}" "
  nohup bash -c '
    cd ${REMOTE_WORKDIR} && \
    ${REMOTE_PYTHON} ${DEPLOY_DIR}/scripts/rnaseq_run.py \
      -c configs/analysis_case.yaml --steps preprocessing \
    && touch ${REMOTE_DONE} \
    || touch ${REMOTE_FAIL}
  ' > ${REMOTE_LOG} 2>&1 &
  echo \"  PID \$! started — log: ${REMOTE_LOG}\"
"

echo "  Preprocessing running in background on remote."
echo "  Polling every 60s for completion..."

# ── Step 3: Poll until done or failed ────────────────────────────────────────

ELAPSED=0
while true; do
    sleep 60
    ELAPSED=$((ELAPSED + 60))

    STATUS=$(ssh "${HOST}" "
        if [ -f ${REMOTE_DONE} ]; then echo done;
        elif [ -f ${REMOTE_FAIL} ]; then echo failed;
        else echo running; fi
    " 2>/dev/null || echo "ssh_error")

    TAIL=$(ssh "${HOST}" "tail -3 ${REMOTE_LOG} 2>/dev/null" 2>/dev/null || echo "(log unavailable)")

    echo "  [${ELAPSED}s] ${STATUS} — ${TAIL}"

    if [ "$STATUS" = "done" ]; then
        break
    elif [ "$STATUS" = "failed" ]; then
        echo ""
        echo "ERROR: Remote preprocessing failed. Last log:"
        ssh "${HOST}" "tail -30 ${REMOTE_LOG} 2>/dev/null" || true
        exit 1
    fi
done

echo ""
echo "[3/4] Preprocessing complete. Pulling gene_counts.tsv..."
LOCAL_DIR="results/shared/${BATCH}/01_preprocessing/gene_counts"
mkdir -p "$LOCAL_DIR"

REMOTE_COUNTS="${REMOTE_WORKDIR}/results/shared/${BATCH}/01_preprocessing/gene_counts/gene_counts.tsv"
scp "${HOST}:${REMOTE_COUNTS}" "${LOCAL_DIR}/gene_counts.tsv" || {
    echo "ERROR: Could not pull gene_counts.tsv from remote."
    echo "Check remote path: ${REMOTE_COUNTS}"
    exit 1
}

# ── Step 4: Verify ───────────────────────────────────────────────────────────

echo ""
echo "[4/4] Verifying results..."
LINES=$(wc -l < "${LOCAL_DIR}/gene_counts.tsv" | tr -d ' ')
COLS=$(head -1 "${LOCAL_DIR}/gene_counts.tsv" | awk -F'\t' '{print NF-1}')

echo ""
echo "============================================================"
echo "Done. Results at: ${LOCAL_DIR}/gene_counts.tsv"
echo "  ${LINES} genes x ${COLS} samples"
echo "============================================================"
echo ""
echo "Next: run local analysis:"
echo "  python \${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c $CONFIG --steps 1a-5b"
