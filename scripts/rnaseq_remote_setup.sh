#!/bin/bash
# Deploy rnaseq-analysis-plugin on a remote server.
#
# Run this script ONCE as shaohua to set up the shared code installation.
# Other users can then call the pipeline via the wrapper script.
#
# Usage:
#   bash rnaseq_remote_setup.sh                          # deploy locally (on the server itself)
#   bash rnaseq_remote_setup.sh -H user@host             # deploy via SSH
#   bash rnaseq_remote_setup.sh -H user@host -d /path    # custom deploy path

set -euo pipefail

REMOTE_HOST=""
DEPLOY_DIR=""
REPO_URL="https://github.com/shaohuaxu-evoscape/rnaseq-analysis-plugin.git"
CONDA_ENV="rnaseq"

usage() {
    echo "Usage: $0 [-H user@host] [-d deploy_dir] [-r repo_url] [-e conda_env]"
    echo ""
    echo "Options:"
    echo "  -H  Remote host (user@host). Omit to deploy on local machine."
    echo "  -d  Deploy directory (default: ~/rnaseq-analysis-plugin)"
    echo "  -r  Git repo URL (default: $REPO_URL)"
    echo "  -e  Conda environment name (default: $CONDA_ENV)"
    exit 1
}

while getopts "H:d:r:e:h" opt; do
    case $opt in
        H) REMOTE_HOST="$OPTARG" ;;
        d) DEPLOY_DIR="$OPTARG" ;;
        r) REPO_URL="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Build the setup commands
SETUP_SCRIPT=$(cat << 'SETUP_EOF'
set -e

DEPLOY_DIR="__DEPLOY_DIR__"
REPO_URL="__REPO_URL__"
CONDA_ENV="__CONDA_ENV__"

# Default deploy dir
if [ -z "$DEPLOY_DIR" ]; then
    DEPLOY_DIR="$HOME/rnaseq-analysis-plugin"
fi

echo "============================================================"
echo "RNA-seq Analysis Plugin — Remote Setup"
echo "============================================================"
echo "  Deploy dir:   $DEPLOY_DIR"
echo "  Repo:         $REPO_URL"
echo "  Conda env:    $CONDA_ENV"
echo "============================================================"

# ── Step 1: Clone or update repo ─────────────────────────────────────────────
echo ""
echo "[1/4] Setting up code..."
if [ -d "$DEPLOY_DIR/.git" ]; then
    echo "  Repo exists, pulling latest..."
    cd "$DEPLOY_DIR" && git pull --quiet
else
    echo "  Cloning repo..."
    git clone "$REPO_URL" "$DEPLOY_DIR"
fi

# ── Step 2: Verify conda environment ─────────────────────────────────────────
echo ""
echo "[2/4] Checking conda environment '$CONDA_ENV'..."
if conda env list 2>/dev/null | grep -q "^${CONDA_ENV} "; then
    echo "  Conda env '$CONDA_ENV' exists"
    # Verify key packages
    conda run -n "$CONDA_ENV" python -c "import pandas, numpy, scipy, yaml; print('  Core packages OK')"
    conda run -n "$CONDA_ENV" python -c "import pydeseq2; print('  pydeseq2 OK')" 2>/dev/null || echo "  WARNING: pydeseq2 not installed (needed for DE analysis)"
else
    echo "  WARNING: Conda env '$CONDA_ENV' not found."
    echo "  Create it with: conda create -n $CONDA_ENV python=3.10 pandas numpy scipy scikit-learn matplotlib statsmodels pyyaml openpyxl"
    echo "  Then: conda activate $CONDA_ENV && pip install pydeseq2"
fi

# ── Step 3: Check bioinformatics tools ───────────────────────────────────────
echo ""
echo "[3/4] Checking bioinformatics tools..."
TOOLS_OK=true
for tool in fastp hisat2 samtools htseq-count; do
    if which "$tool" > /dev/null 2>&1; then
        VERSION=$("$tool" --version 2>&1 | head -1 || echo "unknown")
        echo "  $tool: $VERSION"
    else
        echo "  WARNING: $tool not found in PATH"
        TOOLS_OK=false
    fi
done
if [ "$TOOLS_OK" = false ]; then
    echo ""
    echo "  Install missing tools via conda:"
    echo "    conda install -n $CONDA_ENV -c bioconda fastp hisat2 samtools htseq"
fi

# ── Step 4: Create wrapper script ────────────────────────────────────────────
echo ""
echo "[4/4] Creating wrapper script..."
WRAPPER="$DEPLOY_DIR/scripts/rnaseq-run"
cat > "$WRAPPER" << WRAPPER_EOF
#!/bin/bash
# Wrapper: activate conda env and run the pipeline.
# Usage: /home/shaohua/rnaseq-analysis-plugin/scripts/rnaseq-run -c config.yaml --steps 1a-5b
SCRIPT_DIR="\$(cd "\$(dirname "\$0")" && pwd)"
exec conda run --no-banner -n ${CONDA_ENV} python "\${SCRIPT_DIR}/rnaseq_run.py" "\$@"
WRAPPER_EOF
chmod +x "$WRAPPER"
echo "  Created: $WRAPPER"

# Also create init-project wrapper
INIT_WRAPPER="$DEPLOY_DIR/scripts/rnaseq-init-project"
cat > "$INIT_WRAPPER" << INIT_EOF
#!/bin/bash
# Initialize a new RNA-seq analysis project in the current directory.
# Usage: cd ~/my-project && /home/shaohua/rnaseq-analysis-plugin/scripts/rnaseq-init-project
SCRIPT_DIR="\$(cd "\$(dirname "\$0")" && pwd)"
PLUGIN_DIR="\$(dirname "\$SCRIPT_DIR")"

echo "Initializing RNA-seq analysis project in: \$(pwd)"
mkdir -p configs inputs/ref results
cp "\${PLUGIN_DIR}/configs/analysis_case.template.yaml" configs/analysis_case.yaml
echo ""
echo "Project scaffolded:"
echo "  configs/analysis_case.yaml  — edit this with your project settings"
echo "  inputs/                     — place your data here"
echo "  results/                    — pipeline outputs will go here"
echo ""
echo "Next: edit configs/analysis_case.yaml, then run:"
echo "  \${SCRIPT_DIR}/rnaseq-run -c configs/analysis_case.yaml --steps 1a-5b"
INIT_EOF
chmod +x "$INIT_WRAPPER"
echo "  Created: $INIT_WRAPPER"

# ── Done ─────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "Setup complete!"
echo ""
echo "For other users:"
echo "  1. cd ~/my-project"
echo "  2. $DEPLOY_DIR/scripts/rnaseq-init-project"
echo "  3. Edit configs/analysis_case.yaml"
echo "  4. $DEPLOY_DIR/scripts/rnaseq-run -c configs/analysis_case.yaml --steps 1a-5b"
echo "============================================================"
SETUP_EOF

# Replace placeholders
SETUP_SCRIPT="${SETUP_SCRIPT//__DEPLOY_DIR__/$DEPLOY_DIR}"
SETUP_SCRIPT="${SETUP_SCRIPT//__REPO_URL__/$REPO_URL}"
SETUP_SCRIPT="${SETUP_SCRIPT//__CONDA_ENV__/$CONDA_ENV}"

# Execute
if [ -n "$REMOTE_HOST" ]; then
    echo "Deploying to $REMOTE_HOST..."
    ssh "$REMOTE_HOST" "bash -s" <<< "$SETUP_SCRIPT"
else
    bash -c "$SETUP_SCRIPT"
fi
