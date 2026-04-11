"""Runtime config loading and accessor helpers."""

import os
import re
import sys
from datetime import datetime

import pandas as pd
import yaml

from .validation import validate_pipeline_config, validate_sample_manifest


def _find_plugin_root():
    """Find the plugin installation root (where configs/step_registry.yaml lives).

    Walks up from the scripts/core/ directory to find the plugin root.
    Works whether the code is in a plugin install, a symlink, or standalone.
    """
    # Walk up from this file's real location
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.isfile(os.path.join(d, "configs", "step_registry.yaml")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    # Fallback: try cwd
    if os.path.isfile(os.path.join(os.getcwd(), "configs", "step_registry.yaml")):
        return os.getcwd()
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _find_project_root():
    """Discover the user's project root (cwd or nearest dir with pyproject.toml).

    This is the directory where the user runs the pipeline — contains their
    configs/analysis_case.yaml, inputs/, results/, etc.
    """
    d = os.path.abspath(os.getcwd())
    for _ in range(4):
        if os.path.isfile(os.path.join(d, "configs", "analysis_case.yaml")):
            return d
        if os.path.isfile(os.path.join(d, "pyproject.toml")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return os.getcwd()


PLUGIN_ROOT = _find_plugin_root()
PROJECT_ROOT = _find_project_root()


def _parse_columns(gene_counts_path):
    """Parse gene_counts.tsv header -> dict of {condition: [timepoints]}."""
    with open(gene_counts_path, "r") as f:
        header = f.readline().strip().split("\t")

    if len(header) < 2:
        raise ValueError(
            f"gene_counts file appears empty or has no sample columns: {gene_counts_path}"
        )

    pattern = re.compile(r"^(.+)-(\d+)$")
    tp_by_cond = {}
    for col in header[1:]:
        m = pattern.match(col)
        if m:
            cond, tp = m.group(1), int(m.group(2))
            tp_by_cond.setdefault(cond, set()).add(tp)
    return {c: sorted(tps) for c, tps in tp_by_cond.items()}


def load_config(config_path=None, overrides=None):
    """Load pipeline config from YAML file, with optional CLI overrides."""
    if config_path is None:
        config_path = os.path.join(PROJECT_ROOT, "configs", "analysis_case.yaml")
    config_path = os.path.abspath(config_path)

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    cfg = validate_pipeline_config(cfg)
    cfg["paths"].setdefault("raw_fastq_dir", "")
    cfg["paths"].setdefault("gene_counts", "")
    cfg["paths"].setdefault("sample_manifest", "")
    cfg["paths"].setdefault("shared_output_dir", "results/shared/{batch_id}")
    cfg["paths"].setdefault("output_dir", "results/{batch_id}/{pair_name}")

    if overrides:
        if "batch_id" in overrides:
            cfg["project"]["batch_id"] = overrides["batch_id"]
        if "conditions" in overrides:
            cfg["experiment"]["conditions"] = overrides["conditions"]
        paths = cfg["paths"]
        if "raw_fastq_dir" in overrides:
            paths["raw_fastq_dir"] = overrides["raw_fastq_dir"]
        if "gene_counts" in overrides:
            paths["gene_counts"] = overrides["gene_counts"]
        if "sample_manifest" in overrides:
            paths["sample_manifest"] = overrides["sample_manifest"]
        if "shared_output_dir" in overrides:
            paths["shared_output_dir"] = overrides["shared_output_dir"]
        if "output_dir" in overrides:
            paths["output_dir"] = overrides["output_dir"]
        if "module1_source_conditions" in overrides:
            cfg["_module1_source_conditions"] = sorted({
                str(item) for item in overrides["module1_source_conditions"] or []
            })
        if "shared_batch_name" in overrides:
            cfg["_shared_batch_name"] = str(overrides["shared_batch_name"])

    cfg["_project_root"] = PROJECT_ROOT
    cfg["_pipeline_root"] = PROJECT_ROOT

    conditions = cfg["experiment"]["conditions"]
    pair_name = "-".join(conditions[:2]) if len(conditions) >= 2 else conditions[0]
    cfg["_pair_name"] = pair_name

    batch_id = (
        cfg["project"].get("batch_id")
        or cfg.get("run_name")
        or os.path.splitext(os.path.basename(config_path))[0]
    )
    cfg["project"]["batch_id"] = batch_id
    pipeline_root = cfg["_pipeline_root"]
    for key, val in cfg["paths"].items():
        if isinstance(val, str):
            resolved = val.format(batch_id=batch_id, pair_name=pair_name)
            if os.path.isabs(resolved):
                cfg["paths"][key] = resolved
            else:
                cfg["paths"][key] = os.path.join(pipeline_root, resolved)

    manifest_path = cfg["paths"].get("sample_manifest", "")
    if manifest_path and os.path.isfile(manifest_path):
        manifest = validate_sample_manifest(pd.read_csv(manifest_path, sep="\t"))
        if "target_condition" in manifest.columns and "timepoint" in manifest.columns:
            cond_map = {}
            for cond, sub in manifest.groupby("target_condition"):
                cond_map[str(cond)] = sorted(sub["timepoint"].astype(int).unique().tolist())
            cfg["_all_conditions"] = cond_map
            selected = [c for c in get_conditions(cfg) if c in cond_map]
            tp_sets = [set(cond_map[c]) for c in selected]
            shared = set.intersection(*tp_sets) if tp_sets else set()
            cfg["experiment"]["timepoints"] = sorted(shared)
        return validate_pipeline_config(cfg)

    gc_path = cfg["paths"].get("gene_counts", "")
    if gc_path and os.path.isfile(gc_path):
        col_map = _parse_columns(gc_path)
        cfg["_all_conditions"] = col_map
        tp_sets = [set(col_map.get(c, [])) for c in conditions]
        shared = set.intersection(*tp_sets) if tp_sets else set()
        cfg["experiment"]["timepoints"] = sorted(shared)

    return validate_pipeline_config(cfg)


def get_batch_id(cfg):
    """Return validated batch/run identifier."""
    batch_id = str(cfg["project"]["batch_id"]).strip()
    if not batch_id:
        raise ValueError("pipeline config project.batch_id must be set")
    return batch_id


def get_conditions(cfg):
    """Return selected comparison conditions."""
    conditions = cfg.get("experiment", {}).get("conditions")
    validated = validate_pipeline_config(
        {"project": {}, "experiment": {"conditions": conditions}, "paths": {}}
    )
    return validated["experiment"]["conditions"]


def get_timepoints(cfg):
    """Return selected shared timepoints."""
    timepoints = cfg.get("experiment", {}).get("timepoints")
    if timepoints is None:
        raise ValueError("pipeline config experiment.timepoints is not set")
    return [int(item) for item in timepoints]


def get_pair_name(cfg):
    """Return current comparison pair identifier."""
    pair_name = str(cfg.get("_pair_name", "")).strip()
    if not pair_name:
        conditions = get_conditions(cfg)
        pair_name = "-".join(conditions[:2]) if len(conditions) >= 2 else conditions[0]
    return pair_name


def get_path(cfg, key):
    """Return a resolved path from cfg['paths']."""
    paths = cfg.get("paths", {})
    if key not in paths:
        raise KeyError(f"pipeline config paths.{key} is not defined")
    return paths[key]


def ensure_output_dir(cfg, subdir=""):
    """Create output directory under comparison-specific output_dir."""
    out = get_path(cfg, "output_dir")
    if subdir:
        out = os.path.join(out, subdir)
    os.makedirs(out, exist_ok=True)
    return out


def ensure_shared_dir(cfg, subdir=""):
    """Create output directory under shared_output_dir (batch-level, all samples)."""
    out = get_path(cfg, "shared_output_dir")
    if subdir:
        out = os.path.join(out, subdir)
    os.makedirs(out, exist_ok=True)
    return out


def write_readme(out_dir, step_id, step_name, file_descriptions):
    """Write a README.txt describing the contents of an output directory."""
    lines = [
        f"Step {step_id}: {step_name}",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "Files:",
    ]
    for fname, desc in file_descriptions.items():
        lines.append(f"  {fname:40s} — {desc}")
    lines.append("")

    with open(os.path.join(out_dir, "README.txt"), "w") as f:
        f.write("\n".join(lines))
