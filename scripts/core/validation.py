"""Validation helpers for analysis-case configs, pipeline configs, and manifests."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd


def _require_mapping(obj, label):
    if not isinstance(obj, Mapping):
        raise ValueError(f"{label} must be a mapping/dictionary")


def _normalize_conditions(conditions, label):
    if not isinstance(conditions, list) or not conditions:
        raise ValueError(f"{label} must be a non-empty list")
    normalized = []
    for item in conditions:
        text = str(item).strip()
        if not text:
            raise ValueError(f"{label} cannot contain empty condition names")
        normalized.append(text)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{label} contains duplicate condition names")
    return normalized


def _normalize_timepoints(timepoints, label):
    if timepoints is None:
        return None
    if not isinstance(timepoints, list):
        raise ValueError(f"{label} must be a list of integers")
    normalized = []
    for item in timepoints:
        try:
            normalized.append(int(item))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} contains a non-integer value: {item!r}") from exc
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{label} contains duplicate timepoints")
    return sorted(normalized)


def validate_analysis_case_config(cfg):
    """Validate and normalize a user-facing analysis-case config."""
    _require_mapping(cfg, "analysis case")
    if "batches" not in cfg or "target_conditions" not in cfg:
        raise ValueError("analysis case must define both 'batches' and 'target_conditions'")

    batches = cfg["batches"]
    target_conditions = cfg["target_conditions"]
    _require_mapping(batches, "analysis case 'batches'")
    _require_mapping(target_conditions, "analysis case 'target_conditions'")
    if not batches:
        raise ValueError("analysis case 'batches' cannot be empty")
    if not target_conditions:
        raise ValueError("analysis case 'target_conditions' cannot be empty")

    for batch_id, batch_cfg in batches.items():
        if not str(batch_id).strip():
            raise ValueError("analysis case batch IDs cannot be empty")
        if batch_cfg is not None and not isinstance(batch_cfg, Mapping):
            raise ValueError(f"batch '{batch_id}' must map to a dictionary or null")

    normalized_targets = {}
    for target_name, members in target_conditions.items():
        target_name = str(target_name).strip()
        if not target_name:
            raise ValueError("target condition names cannot be empty")
        if not isinstance(members, list) or not members:
            raise ValueError(f"target condition '{target_name}' must contain at least one batch member")
        normalized_members = []
        for member in members:
            if not isinstance(member, Mapping):
                raise ValueError(f"target condition '{target_name}' contains a non-dictionary member")
            batch_id = str(member.get("batch", "")).strip()
            condition = str(member.get("condition", "")).strip()
            if not batch_id or not condition:
                raise ValueError(
                    f"target condition '{target_name}' members must define non-empty 'batch' and 'condition'"
                )
            if batch_id not in batches:
                raise ValueError(
                    f"target condition '{target_name}' references unknown batch '{batch_id}'"
                )
            normalized_members.append({"batch": batch_id, "condition": condition})
        normalized_targets[target_name] = normalized_members

    cfg["target_conditions"] = normalized_targets
    cfg.setdefault("experiment", {})
    cfg["experiment"]["timepoints"] = _normalize_timepoints(
        cfg["experiment"].get("timepoints"),
        "analysis case experiment.timepoints",
    )
    return cfg


def validate_sample_manifest(manifest, required_columns=None):
    """Validate and normalize a sample manifest dataframe."""
    if not isinstance(manifest, pd.DataFrame):
        raise ValueError("sample manifest must be a pandas DataFrame")

    manifest = manifest.copy()
    condition_col = None
    if "target_condition" in manifest.columns:
        condition_col = "target_condition"
    elif "merged_condition" in manifest.columns:
        condition_col = "merged_condition"

    base_required = {"sample_id", "timepoint"}
    required = set(required_columns or set())
    if condition_col is None:
        missing = ["target_condition"]
    else:
        missing = []
    missing += sorted((base_required | required) - set(manifest.columns))
    if missing:
        raise ValueError(
            f"sample_manifest missing required columns: {', '.join(missing)}"
        )

    if condition_col != "target_condition":
        manifest = manifest.rename(columns={condition_col: "target_condition"})

    for col in ["sample_id", "target_condition", "batch_id", "source_condition", "display_label"]:
        if col in manifest.columns:
            manifest[col] = manifest[col].astype(str).str.strip()

    if manifest["sample_id"].eq("").any():
        raise ValueError("sample_manifest contains empty sample_id values")
    if manifest["sample_id"].duplicated().any():
        dupes = manifest.loc[manifest["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"sample_manifest contains duplicate sample_id values: {', '.join(dupes[:5])}")
    if manifest["target_condition"].eq("").any():
        raise ValueError("sample_manifest contains empty target_condition values")

    try:
        manifest["timepoint"] = pd.to_numeric(manifest["timepoint"]).astype(int)
    except Exception as exc:
        raise ValueError("sample_manifest contains non-integer timepoint values") from exc

    if "display_label" not in manifest.columns:
        if {"batch_id", "source_condition", "timepoint"}.issubset(manifest.columns):
            manifest["display_label"] = manifest.apply(
                lambda row: f"{row['batch_id']}|{row['source_condition']}|{row['timepoint']}h",
                axis=1,
            )
        else:
            manifest["display_label"] = manifest["sample_id"]

    return manifest


def validate_pipeline_config(cfg):
    """Validate and normalize a runner-consumable pipeline config."""
    _require_mapping(cfg, "pipeline config")
    cfg.setdefault("project", {})
    cfg.setdefault("experiment", {})
    cfg.setdefault("paths", {})
    _require_mapping(cfg["project"], "pipeline config 'project'")
    _require_mapping(cfg["experiment"], "pipeline config 'experiment'")
    _require_mapping(cfg["paths"], "pipeline config 'paths'")

    conditions = cfg["experiment"].get("conditions")
    if conditions is None and "target_conditions" in cfg:
        target_conditions = cfg.get("target_conditions") or {}
        if isinstance(target_conditions, Mapping):
            conditions = list(target_conditions.keys())
            cfg["experiment"]["conditions"] = conditions
    cfg["experiment"]["conditions"] = _normalize_conditions(
        cfg["experiment"].get("conditions"),
        "pipeline config experiment.conditions",
    )
    cfg["experiment"]["timepoints"] = _normalize_timepoints(
        cfg["experiment"].get("timepoints"),
        "pipeline config experiment.timepoints",
    )
    cfg["experiment"].setdefault("sample_pattern", "{condition}-{timepoint}")

    for key, value in list(cfg["paths"].items()):
        if value is None:
            cfg["paths"][key] = ""
        elif not isinstance(value, str):
            raise ValueError(f"pipeline config paths.{key} must be a string path")

    batch_id = str(cfg["project"].get("batch_id", "")).strip()
    if batch_id:
        cfg["project"]["batch_id"] = batch_id
    return cfg
