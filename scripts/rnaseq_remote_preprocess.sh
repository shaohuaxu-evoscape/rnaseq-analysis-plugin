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
BATCH=$(read_config "list(c['batches'].keys())[0]")

if [ -z "$DEPLOY_DIR" ]; then
    echo "ERROR: remote.deploy_dir not set in $CONFIG"
    echo "Set it to the plugin install path on the remote server (e.g., /home/shaohua/rnaseq-analysis-plugin)"
    exit 1
fi

echo "============================================================"
echo "Remote RNA-seq Preprocessing"
echo "============================================================"
echo "  Host:       ${USER}@${HOST}"
echo "  Deploy dir: ${DEPLOY_DIR}"
echo "  Conda env:  ${CONDA_ENV}"
echo "  Batch:      ${BATCH}"
echo "============================================================"

# ── Step 1: Sync config to remote ────────────────────────────────────────────

echo ""
echo "[1/4] Syncing config to remote..."
ssh "${USER}@${HOST}" "mkdir -p ~/rnaseq-projects/${BATCH}/configs"
scp "$CONFIG" "${USER}@${HOST}:~/rnaseq-projects/${BATCH}/configs/analysis_case.yaml"
echo "  Config uploaded to ~/rnaseq-projects/${BATCH}/configs/"

# ── Step 2: Run preprocessing on remote ──────────────────────────────────────

echo ""
echo "[2/4] Running preprocessing on remote (steps 0a-0d)..."
ssh "${USER}@${HOST}" "cd ~/rnaseq-projects/${BATCH} && conda run --no-banner -n ${CONDA_ENV} python ${DEPLOY_DIR}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps preprocessing" || {
    echo ""
    echo "ERROR: Remote preprocessing failed."
    echo "Debug with MCP remote-linux: read rnaseq-remote-preprocessing.md for step-by-step instructions."
    exit 1
}

# ── Step 3: Pull results back ────────────────────────────────────────────────

echo ""
echo "[3/4] Pulling gene_counts.tsv back to local..."
LOCAL_DIR="results/shared/${BATCH}/01_preprocessing/gene_counts"
mkdir -p "$LOCAL_DIR"

REMOTE_COUNTS="~/rnaseq-projects/${BATCH}/results/shared/${BATCH}/01_preprocessing/gene_counts/gene_counts.tsv"
scp "${USER}@${HOST}:${REMOTE_COUNTS}" "${LOCAL_DIR}/gene_counts.tsv" || {
    echo "ERROR: Could not pull gene_counts.tsv from remote."
    echo "Check if preprocessing completed successfully."
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
