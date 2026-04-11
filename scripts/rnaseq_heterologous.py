#!/usr/bin/env python3
"""Run Heterologous Gene Analysis (step 4c).

Usage:
    python rnaseq_heterologous.py -c configs/analysis_case.yaml
    python rnaseq_heterologous.py -c configs/analysis_case.yaml --dry-run
"""
import argparse
import os
import sys

# Ensure the plugin scripts package is importable
_plugin_root = os.path.dirname(os.path.abspath(__file__))
_parent = os.path.dirname(_plugin_root)
if _parent not in sys.path:
    sys.path.insert(0, _parent)

from scripts.core.config_runtime import load_config
from scripts.steps.heterologous_genes import run


def main():
    parser = argparse.ArgumentParser(description="Heterologous Gene Analysis (step 4c)")
    parser.add_argument("-c", "--config", required=True, help="Path to analysis_case.yaml")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Print what would run")
    args = parser.parse_args()

    if args.dry_run:
        print("[DRY RUN] Would run step 4c: Heterologous Gene Analysis")
        return

    cfg = load_config(args.config)
    run(cfg)


if __name__ == "__main__":
    main()
