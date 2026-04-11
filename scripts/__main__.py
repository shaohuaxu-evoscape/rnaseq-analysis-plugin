"""CLI entry point: python -m scripts [OPTIONS].

This entrypoint now always follows the analysis-case workflow:
1. resolve or generate per-batch gene_counts.tsv files
2. prepare target-condition inputs for exploratory analysis
3. run the downstream pipeline using the generated analysis config
"""

import argparse
import sys

from .core.prepare_analysis_case import prepare_analysis_case
from .core.runner import run_pipeline, _ensure_registry


def _parse_set_value(raw):
    """Auto-coerce a CLI string value to the appropriate Python type."""
    if raw.lower() in ("true", "yes"):
        return True
    if raw.lower() in ("false", "no"):
        return False
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        pass
    return raw


def _parse_set_args(set_args):
    """Parse --set key=value pairs into a dict of dotted-key → typed value."""
    overrides = {}
    if not set_args:
        return overrides
    for item in set_args:
        if "=" not in item:
            raise ValueError(f"--set requires key=value format, got: {item!r}")
        key, _, raw_value = item.partition("=")
        key = key.strip()
        if not key:
            raise ValueError(f"--set key cannot be empty: {item!r}")
        overrides[key] = _parse_set_value(raw_value.strip())
    return overrides


def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq exploratory analysis pipeline (analysis-case driven)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
  1. Read an analysis-case config passed with -c/--config
  2. Reuse existing batch gene_counts.tsv files or auto-run hidden preprocessing when needed
  3. Build exploratory analysis inputs
  4. Run downstream analysis on the generated config

Recommended:
  python -m scripts -c configs/analysis_case.yaml
  python -m scripts -c configs/analysis_case.yaml --steps 1a-2d
  python -m scripts --list-steps           # show available pipeline steps

Override config parameters:
  python -m scripts -c configs/analysis_case.yaml --set normalization.method=cpm
  python -m scripts -c configs/analysis_case.yaml --set normalization.gene_filtering.min_expr=2.0 --set differential_analysis.de_screening.fdr_threshold=0.01

Note:
  The legacy direct single-matrix entry mode has been removed.
""",
    )
    parser.add_argument(
        "--config", "-c",
        help="Path to analysis-case YAML",
    )
    parser.add_argument(
        "--steps", "-s",
        default="1a-2d",
        help="Pipeline steps to run after hidden preprocessing/analysis-case preparation (default: 1a-2d)",
    )
    parser.add_argument(
        "--set",
        action="append",
        dest="set_overrides",
        metavar="KEY=VALUE",
        help="Override a config parameter using dot-separated key path (repeatable)",
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Prepare analysis-case inputs, then print what downstream steps would run",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from pipeline_run_manifest.json and skip already completed steps for the same resolved config",
    )
    parser.add_argument(
        "--list-steps", "-l",
        action="store_true",
        help="List all available pipeline steps and exit",
    )

    args = parser.parse_args()

    if args.list_steps:
        import re as _re
        _ensure_registry()
        from .core.runner import STEPS, MODULE_DISPLAY_NAMES
        print("\nRNA-seq Pipeline Steps:")
        print("=" * 50)
        current_module = None
        for sid, mod_path, name in STEPS:
            m = _re.match(r"(\d+)", sid)
            if not m:
                continue
            module = m.group(1)
            if module != current_module:
                current_module = module
                module_name = MODULE_DISPLAY_NAMES.get(module, f"Module {module}")
                print(f"\n  Module {module}: {module_name}")
            print(f"    {sid}  {name:30s}  ({mod_path})")
        print()
        sys.exit(0)

    if not args.config:
        parser.error("--config/-c is required")

    try:
        set_overrides = _parse_set_args(args.set_overrides)
    except ValueError as exc:
        parser.error(str(exc))

    prepared = prepare_analysis_case(args.config, case_overrides=set_overrides or None)
    analysis_config = prepared["analysis_config"]

    try:
        result = run_pipeline(
            config_path=analysis_config,
            steps=args.steps,
            dry_run=args.dry_run,
            overrides=None,
            resume=args.resume,
        )
    except ValueError as exc:
        parser.error(str(exc))

    status = result.get("status")
    if status == "complete":
        n_err = sum(
            1 for r in result.get("results", {}).values()
            if r["status"] == "error"
        )
        sys.exit(1 if n_err > 0 else 0)
    elif status == "dry_run":
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
