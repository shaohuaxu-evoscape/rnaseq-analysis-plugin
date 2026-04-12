# Remote Deployment

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview.

Deploy the RNA-seq analysis pipeline on a remote server for shared use. Code is maintained by one admin (shaohua), other users call it to run QC preprocessing with results saved in their own directories.

## Architecture

```
Remote Server
├── /home/shaohua/rnaseq-analysis-plugin/       ← Shared code (admin maintains)
│   ├── scripts/rnaseq-run                      ← Wrapper (auto-activates conda)
│   ├── scripts/rnaseq-init-project             ← Project scaffolding for new users
│   ├── scripts/core/                           ← Pipeline core modules
│   ├── scripts/steps/                          ← Step implementations
│   └── configs/step_registry.yaml              ← Step definitions
│
├── /home/userA/my-project/                     ← User A's project
│   ├── configs/analysis_case.yaml
│   ├── inputs/ → /data/rna/batch1              ← Symlink to shared data
│   └── results/                                ← User A's outputs
│
└── /home/userB/another-project/                ← User B's project
    ├── configs/analysis_case.yaml
    └── results/                                ← User B's outputs
```

## Admin Setup (One-Time)

Run on the remote server as shaohua, or from local via SSH:

```bash
# From local machine
bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_setup.sh -H shaohua@azure

# Or directly on the remote server
bash rnaseq_remote_setup.sh
```

This will:
1. Clone the plugin repo to `~/rnaseq-analysis-plugin/`
2. Verify conda `rnaseq` environment and Python packages
3. Check bioinformatics tools (fastp, hisat2, samtools, htseq-count)
4. Create wrapper scripts (`rnaseq-run`, `rnaseq-init-project`)

### Requirements

| Component | Requirement |
|-----------|-------------|
| OS | Linux (Ubuntu 20.04+ / CentOS 7+) |
| Python | 3.10+ in conda `rnaseq` environment |
| RAM | >= 16 GB (HISAT2 index building) |
| Disk | >= 50 GB per batch (intermediate BAM files) |
| Tools | fastp >= 0.23, hisat2 >= 2.2, samtools >= 1.15, htseq-count >= 2.0 |
| Network | Git access to GitHub (for clone/update) |

### Install Missing Tools

```bash
conda activate rnaseq
conda install -c bioconda fastp hisat2 samtools htseq
pip install pydeseq2
```

### Update the Deployment

```bash
cd ~/rnaseq-analysis-plugin && git pull
```

Or re-run the setup script (it pulls if the repo already exists).

## User Workflow

### Step 1: Initialize a Project

```bash
cd ~
mkdir my-rnaseq-project && cd my-rnaseq-project
/home/shaohua/rnaseq-analysis-plugin/scripts/rnaseq-init-project
```

This creates `configs/analysis_case.yaml`, `inputs/`, `results/`.

### Step 2: Configure

Edit `configs/analysis_case.yaml`:
- Set batch name, gene counts path (or raw FASTQ path)
- Set target conditions
- Set organism info (taxid, gene prefix)

### Step 3: Run Pipeline

```bash
# Full pipeline
/home/shaohua/rnaseq-analysis-plugin/scripts/rnaseq-run -c configs/analysis_case.yaml --steps 1a-5b

# Preview only
/home/shaohua/rnaseq-analysis-plugin/scripts/rnaseq-run -c configs/analysis_case.yaml --steps 1a-5b --dry-run

# List available steps
/home/shaohua/rnaseq-analysis-plugin/scripts/rnaseq-run --list-steps
```

Results appear in `./results/` under the user's project directory.

## Local Claude Code Integration

From a local machine with the plugin installed, trigger remote preprocessing via SSH:

```bash
bash ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_remote_preprocess.sh -c configs/analysis_case.yaml
```

Config must include:
```yaml
remote:
  enabled: true
  host: "azure"
  user: "shaohua"
  deploy_dir: "/home/shaohua/rnaseq-analysis-plugin"
  conda_env: "rnaseq"
```

The script:
1. Uploads config to remote
2. Runs preprocessing via the deployed Python scripts
3. Pulls `gene_counts.tsv` back to local for downstream analysis

## Troubleshooting

- **conda not found**: Ensure conda is in PATH. Check `~/.bashrc` or `~/.bash_profile`.
- **Permission denied on wrapper scripts**: Run `chmod +x ~/rnaseq-analysis-plugin/scripts/rnaseq-run`
- **Tools not found**: Activate the conda env first: `conda activate rnaseq && which fastp`
- **Disk full**: Check `df -h`. BAM files are large — clean old results with `rm -rf results/shared/*/01_preprocessing/alignment/`
