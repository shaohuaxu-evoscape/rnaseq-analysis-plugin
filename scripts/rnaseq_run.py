#!/usr/bin/env python3
"""Run RNA-seq analysis pipeline steps.

Usage:
    python rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b
    python rnaseq_run.py -c configs/analysis_case.yaml --steps 3a-3d --dry-run
    python rnaseq_run.py -c configs/analysis_case.yaml --steps all --resume
    python rnaseq_run.py --list-steps
"""
import os
import sys

# Ensure the plugin scripts package is importable
_plugin_root = os.path.dirname(os.path.abspath(__file__))
_parent = os.path.dirname(_plugin_root)
if _parent not in sys.path:
    sys.path.insert(0, _parent)

from scripts.__main__ import main

if __name__ == "__main__":
    main()
