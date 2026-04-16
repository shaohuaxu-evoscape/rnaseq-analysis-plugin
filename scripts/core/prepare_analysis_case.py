"""Prepare per-batch gene-count inputs for an analysis case.

Primary usage:
    python -m scripts -c configs/analysis_case.yaml --steps 1a-2d

This module remains available as a lower-level preparation helper:
    python -m scripts.prepare_analysis_case -c configs/analysis_case.yaml

The analysis-case file declares:
  - where to find each batch's gene_counts.tsv, if already available
  - optional raw FASTQ directory for batches that need module 1 to run
  - target conditions built from batch/condition pairs

Outputs:
  - analysis_gene_counts.tsv
  - sample_manifest.tsv
  - batch_counts_manifest.tsv
  - generated analysis YAML consumed by the main `python -m scripts` entrypoint
"""

from __future__ import annotations

import argparse
import copy
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import yaml

try:
    from .runner import run_pipeline
    from .validation import validate_analysis_case_config, validate_pipeline_config, validate_sample_manifest
except ImportError:  # pragma: no cover
    from scripts.runner import run_pipeline
    from scripts.validation import validate_analysis_case_config, validate_pipeline_config, validate_sample_manifest


PIPELINE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MODULE1_REL = Path("results") / "shared" / "{shared_batch_name}" / "01_preprocessing" / "gene_counts" / "gene_counts.tsv"
COLUMN_RE = re.compile(r"^(.+)-(\d+)$")


def _set_nested(cfg: dict, dotted_key: str, value) -> None:
    """Set a value in a nested dict using dot-separated key path."""
    keys = dotted_key.split(".")
    d = cfg
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def _resolve_path(value: str | None, base_dir: Path) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        case_relative = base_dir / path
        if case_relative.exists() or str(value).startswith("."):
            path = case_relative.resolve()
        else:
            path = (PIPELINE_ROOT / path).resolve()
    return path


def _project_root_from_case(case_path: Path) -> Path:
    """Infer the user project root from an analysis-case path.

    The common layout is `<project>/configs/analysis_case.yaml`. When the config
    does not live under `configs/`, fall back to its parent directory.
    """
    case_path = Path(case_path).resolve()
    if case_path.parent.name == "configs":
        return case_path.parent.parent
    return case_path.parent


def _load_case(case_path: Path) -> dict:
    with open(case_path, "r") as f:
        cfg = yaml.safe_load(f) or {}
    return validate_analysis_case_config(cfg)


def _requested_conditions_by_batch(target_conditions: OrderedDict[str, list]) -> Dict[str, List[str]]:
    requested = {}
    for members in target_conditions.values():
        for member in members:
            requested.setdefault(member["batch"], set()).add(str(member["condition"]))
    return {batch_id: sorted(items) for batch_id, items in requested.items()}


def _shared_batch_name(batch_id: str, batch_cfg: dict, config_path: Path, case_cfg: dict = None) -> str:
    # Allow explicit override via case-level shared_batch_name (used by grid search)
    if case_cfg and case_cfg.get("shared_batch_name"):
        return case_cfg["shared_batch_name"]
    raw_fastq_dir = _resolve_path(batch_cfg.get("raw_fastq_dir"), config_path.parent)
    if raw_fastq_dir is not None:
        return raw_fastq_dir.name
    return batch_id


def _default_counts_path(shared_batch_name: str, config_path: Path) -> Path:
    project_root = _project_root_from_case(config_path)
    return project_root / str(DEFAULT_MODULE1_REL).format(shared_batch_name=shared_batch_name)


def _available_conditions(path: Path) -> set[str]:
    with open(path, "r") as f:
        header = f.readline().strip().split("\t")
    conditions = set()
    for col in header[1:]:
        match = COLUMN_RE.match(col)
        if match:
            conditions.add(match.group(1))
    return conditions


def _count_matrix_covers_conditions(path: Path, requested_conditions: List[str]) -> bool:
    if not path.is_file():
        return False
    if not requested_conditions:
        return True
    available = _available_conditions(path)
    return set(requested_conditions).issubset(available)


def _ensure_batch_counts(
    yaml_batch_id: str,
    shared_batch_name: str,
    batch_cfg: dict,
    config_path: Path,
    raw_overrides: dict,
    requested_conditions: List[str],
) -> Tuple[Path, str]:
    explicit = _resolve_path(batch_cfg.get("gene_counts"), config_path.parent)
    if explicit and explicit.is_file():
        return explicit, "explicit"

    default_path = _default_counts_path(shared_batch_name, config_path)
    if _count_matrix_covers_conditions(default_path, requested_conditions):
        return default_path, "default_output"

    raw_fastq_dir = _resolve_path(batch_cfg.get("raw_fastq_dir"), config_path.parent)
    if raw_fastq_dir is None:
        raise FileNotFoundError(
            f"Batch '{yaml_batch_id}' has no usable gene_counts.tsv and no raw_fastq_dir for module 1 fallback"
        )

    overrides = dict(raw_overrides)
    overrides.update({
        "batch_id": shared_batch_name,
        "shared_batch_name": shared_batch_name,
        "raw_fastq_dir": str(raw_fastq_dir),
        "shared_output_dir": f"results/shared/{shared_batch_name}",
        "module1_source_conditions": requested_conditions,
    })
    result = run_pipeline(
        config_path=str(config_path),
        steps="preprocessing",
        dry_run=False,
        overrides=overrides,
        allow_hidden=True,
    )
    if result.get("status") != "complete":
        raise RuntimeError(f"Module 1 fallback did not complete for batch '{yaml_batch_id}'")
    n_err = sum(
        1 for item in result.get("results", {}).values()
        if item.get("status") == "error"
    )
    if n_err:
        raise RuntimeError(f"Module 1 fallback failed for batch '{yaml_batch_id}' with {n_err} error(s)")
    if not _count_matrix_covers_conditions(default_path, requested_conditions):
        raise FileNotFoundError(
            f"Batch '{yaml_batch_id}' finished module 1 fallback, but usable counts were not found at {default_path}"
        )
    return default_path, "generated"


def _load_count_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = [str(col) for col in df.columns]
    return df


def _condition_timepoints(df: pd.DataFrame, source_condition: str) -> List[int]:
    vals = []
    for col in df.columns:
        match = COLUMN_RE.match(col)
        if match and match.group(1) == source_condition:
            vals.append(int(match.group(2)))
    return sorted(set(vals))


def _shared_timepoints(
    target_conditions: OrderedDict[str, list],
    batch_tables: Dict[str, pd.DataFrame],
) -> List[int]:
    per_cond = []
    for _, members in target_conditions.items():
        available = None
        for member in members:
            batch_id = member["batch"]
            cond = member["condition"]
            member_tps = set(_condition_timepoints(batch_tables[batch_id], cond))
            if available is None:
                available = member_tps
            else:
                available &= member_tps
        per_cond.append(available or set())
    if not per_cond:
        return []
    return sorted(set.intersection(*per_cond))


def _selected_timepoints(case_cfg: dict, shared_timepoints: List[int]) -> List[int]:
    """Apply optional user-specified timepoint filtering on top of shared timepoints."""
    experiment_cfg = case_cfg.get("experiment", {}) or {}
    requested = experiment_cfg.get("timepoints")
    if not requested:
        return shared_timepoints

    requested_set = {int(tp) for tp in requested}
    selected = [tp for tp in shared_timepoints if tp in requested_set]
    if not selected:
        raise ValueError(
            "Configured experiment.timepoints have no overlap with the shared timepoints "
            f"available in the selected target conditions: {shared_timepoints}"
        )
    return selected


def _select_columns(
    target_conditions: OrderedDict[str, list],
    batch_tables: Dict[str, pd.DataFrame],
    batch_name_map: Dict[str, str],
    shared_timepoints: List[int],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    selected_frames = []
    manifest_rows = []

    for target_condition, members in target_conditions.items():
        for member in members:
            batch_id = member["batch"]
            source_condition = member["condition"]
            shared_batch_name = batch_name_map[batch_id]
            df = batch_tables[batch_id]

            for tp in shared_timepoints:
                source_col = f"{source_condition}-{tp}"
                if source_col not in df.columns:
                    continue
                sample_id = f"{target_condition}_{shared_batch_name}_{source_condition}-{tp}"
                selected_frames.append(df[[source_col]].rename(columns={source_col: sample_id}))
                manifest_rows.append({
                    "sample_id": sample_id,
                    "target_condition": target_condition,
                    "batch_id": shared_batch_name,
                    "source_condition": source_condition,
                    "timepoint": tp,
                    "display_label": f"{shared_batch_name}|{source_condition}|{tp}h",
                })

    if not selected_frames:
        raise ValueError("No columns selected for analysis-case matrix")

    combined = pd.concat(selected_frames, axis=1)
    import logging
    _log = logging.getLogger("prepare_analysis_case")

    n_missing = int(combined.isna().sum().sum())
    if n_missing > 0:
        _log.warning(
            f"  {n_missing} gene/sample cells were absent in one batch and filled with 0"
        )
    filled = combined.fillna(0)
    non_integer_mask = (filled != filled.round())
    n_non_int = int(non_integer_mask.sum().sum())
    if n_non_int > 0:
        _log.warning(
            f"  {n_non_int} non-integer values found in count matrix; rounding to nearest int"
        )
    analysis_counts = filled.round().astype(int)
    manifest = pd.DataFrame(manifest_rows)
    manifest = manifest.sort_values(
        ["target_condition", "timepoint", "batch_id", "source_condition", "sample_id"]
    ).reset_index(drop=True)
    analysis_counts = analysis_counts[manifest["sample_id"].tolist()]
    return analysis_counts, manifest


def _write_analysis_config(
    config_path: Path,
    analysis_config_path: Path,
    analysis_counts_path: Path,
    sample_manifest_path: Path,
    output_root: Path,
    shared_results_root: Path,
    conditions: List[str],
    run_name: str,
    case_cfg: dict | None = None,
) -> None:
    if case_cfg is not None:
        cfg = copy.deepcopy(case_cfg)
    else:
        with open(config_path, "r") as f:
            cfg = copy.deepcopy(yaml.safe_load(f) or {})
    for key in [
        "run_name",
        "output_dir",
        "analysis_config",
        "module1_overrides",
        "batches",
        "target_conditions",
    ]:
        cfg.pop(key, None)
    cfg.setdefault("project", {})
    cfg["project"]["batch_id"] = run_name
    cfg.setdefault("experiment", {})
    cfg["experiment"]["conditions"] = conditions
    cfg.setdefault("paths", {})
    cfg["paths"]["gene_counts"] = str(analysis_counts_path)
    cfg["paths"]["sample_manifest"] = str(sample_manifest_path)
    cfg["paths"]["shared_output_dir"] = str(shared_results_root.resolve())
    cfg["paths"]["output_dir"] = str(output_root / "{pair_name}")
    cfg["documents_dir"] = str(output_root / "documents")

    analysis_config_path.parent.mkdir(parents=True, exist_ok=True)
    cfg = validate_pipeline_config(cfg)
    with open(analysis_config_path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)


def _resolve_analysis_output_paths(case_cfg: dict, case_path: Path) -> tuple[str, Path, Path, Path]:
    """Resolve run name, output root, shared root, and generated config path."""
    case_dir = case_path.parent
    project_root = _project_root_from_case(case_path)
    run_name = str(case_cfg.get("run_name", case_path.stem))
    output_root = _resolve_path(case_cfg.get("output_dir"), case_dir)
    if output_root is None:
        output_root = (project_root / "results" / run_name).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    shared_results_root = (project_root / "results" / "shared").resolve()
    shared_results_root.mkdir(parents=True, exist_ok=True)

    analysis_config_path = _resolve_path(case_cfg.get("analysis_config"), case_dir)
    if analysis_config_path is None:
        analysis_config_path = (output_root / "analysis_config.yaml").resolve()
    return run_name, output_root, shared_results_root, analysis_config_path


def _prepare_batch_tables(
    case_cfg: dict,
    case_path: Path,
    target_conditions: OrderedDict[str, list],
    raw_overrides: dict,
) -> tuple[Dict[str, str], Dict[str, pd.DataFrame], list[dict]]:
    """Resolve counts for each batch and load normalized per-batch tables."""
    requested_by_batch = _requested_conditions_by_batch(target_conditions)
    batch_name_map = {
        batch_id: _shared_batch_name(batch_id, batch_cfg or {}, case_path, case_cfg=case_cfg)
        for batch_id, batch_cfg in case_cfg["batches"].items()
    }
    batch_tables = {}
    batch_manifest_rows = []

    for batch_id, batch_cfg in case_cfg["batches"].items():
        batch_cfg = batch_cfg or {}
        shared_batch_name = batch_name_map[batch_id]
        raw_fastq_dir = _resolve_path(batch_cfg.get("raw_fastq_dir"), case_path.parent)
        counts_path, mode = _ensure_batch_counts(
            yaml_batch_id=batch_id,
            shared_batch_name=shared_batch_name,
            batch_cfg=batch_cfg,
            config_path=case_path,
            raw_overrides=raw_overrides,
            requested_conditions=requested_by_batch.get(batch_id, []),
        )
        batch_tables[batch_id] = _load_count_matrix(counts_path)
        batch_manifest_rows.append({
            "yaml_batch_id": batch_id,
            "batch_id": shared_batch_name,
            "resolution_mode": mode,
            "gene_counts_path": str(counts_path),
            "raw_fastq_dir": str(raw_fastq_dir) if raw_fastq_dir else "",
            "requested_conditions": ",".join(requested_by_batch.get(batch_id, [])),
        })

    return batch_name_map, batch_tables, batch_manifest_rows


def _write_analysis_artifacts(
    output_root: Path,
    analysis_counts: pd.DataFrame,
    sample_manifest: pd.DataFrame,
    batch_manifest_rows: list[dict],
) -> tuple[Path, Path, Path]:
    """Persist generated analysis-case artifacts under the output root."""
    analysis_counts_path = output_root / "analysis_gene_counts.tsv"
    sample_manifest_path = output_root / "sample_manifest.tsv"
    batch_manifest_path = output_root / "batch_counts_manifest.tsv"

    analysis_counts.to_csv(analysis_counts_path, sep="\t", index_label="Gene")
    sample_manifest.to_csv(sample_manifest_path, sep="\t", index=False)
    pd.DataFrame(batch_manifest_rows).to_csv(batch_manifest_path, sep="\t", index=False)
    return analysis_counts_path, sample_manifest_path, batch_manifest_path


def prepare_analysis_case(case_config_path: str | Path, case_overrides: dict | None = None) -> dict:
    """Prepare exploratory inputs from an analysis-case config.

    Parameters
    ----------
    case_config_path : str | Path
        Path to the analysis-case YAML file.
    case_overrides : dict, optional
        Dot-separated key → value overrides applied on top of the loaded
        case config (e.g. ``{"normalization.method": "cpm"}``).
    """
    case_path = Path(case_config_path).resolve()
    case_cfg = _load_case(case_path)
    if case_overrides:
        for key, value in case_overrides.items():
            _set_nested(case_cfg, key, value)
    run_name, output_root, shared_results_root, analysis_config_path = _resolve_analysis_output_paths(case_cfg, case_path)

    target_conditions = OrderedDict(case_cfg["target_conditions"])
    if len(target_conditions) != 2:
        raise ValueError("Exploratory analysis currently supports exactly 2 target conditions")

    raw_overrides = case_cfg.get("module1_overrides", {}) or {}
    batch_name_map, batch_tables, batch_manifest_rows = _prepare_batch_tables(
        case_cfg=case_cfg,
        case_path=case_path,
        target_conditions=target_conditions,
        raw_overrides=raw_overrides,
    )

    shared_timepoints = _shared_timepoints(target_conditions, batch_tables)
    if not shared_timepoints:
        raise ValueError("Target conditions have no shared timepoints")
    shared_timepoints = _selected_timepoints(case_cfg, shared_timepoints)

    analysis_counts, sample_manifest = _select_columns(
        target_conditions=target_conditions,
        batch_tables=batch_tables,
        batch_name_map=batch_name_map,
        shared_timepoints=shared_timepoints,
    )
    sample_manifest = validate_sample_manifest(
        sample_manifest,
        required_columns={"batch_id", "source_condition", "display_label"},
    )

    analysis_counts_path, sample_manifest_path, batch_manifest_path = _write_analysis_artifacts(
        output_root=output_root,
        analysis_counts=analysis_counts,
        sample_manifest=sample_manifest,
        batch_manifest_rows=batch_manifest_rows,
    )

    _write_analysis_config(
        config_path=case_path,
        analysis_config_path=analysis_config_path,
        analysis_counts_path=analysis_counts_path,
        sample_manifest_path=sample_manifest_path,
        output_root=output_root,
        shared_results_root=shared_results_root,
        conditions=list(target_conditions.keys()),
        run_name=run_name,
        case_cfg=case_cfg,
    )

    return {
        "analysis_case": str(case_path),
        "output_root": str(output_root),
        "analysis_gene_counts": str(analysis_counts_path),
        "sample_manifest": str(sample_manifest_path),
        "batch_manifest": str(batch_manifest_path),
        "analysis_config": str(analysis_config_path),
        "conditions": list(target_conditions.keys()),
        "timepoints": shared_timepoints,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare analysis-case inputs for the exploratory workflow"
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        help="Path to analysis-case YAML",
    )
    args = parser.parse_args()

    result = prepare_analysis_case(args.config)

    print("Preparation complete")
    print(f"  analysis counts:  {result['analysis_gene_counts']}")
    print(f"  sample manifest:  {result['sample_manifest']}")
    print(f"  batch manifest:   {result['batch_manifest']}")
    print(f"  analysis config:  {result['analysis_config']}")
    print("Primary entrypoint:")
    print(f"  python -m scripts -c {args.config} --steps 1a-2d")


if __name__ == "__main__":
    main()
