#!/usr/bin/env python3
"""Parameter optimization via grid search.

Usage:
    python rnaseq_grid_search.py --mode step1 -c configs/analysis_case.yaml
    python rnaseq_grid_search.py --mode step0 -c configs/analysis_case.yaml
    python rnaseq_grid_search.py --mode step1 -c configs/analysis_case.yaml --collect-only
"""
import os
import sys

_plugin_root = os.path.dirname(os.path.abspath(__file__))
_parent = os.path.dirname(_plugin_root)
if _parent not in sys.path:
    sys.path.insert(0, _parent)

from scripts.core.grid_search import main

if __name__ == "__main__":
    main()
