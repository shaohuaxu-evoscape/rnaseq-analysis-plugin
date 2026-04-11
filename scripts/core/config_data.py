"""Sample manifest, sample naming, and expression I/O helpers."""

import os
from pathlib import Path

import pandas as pd

from .config_runtime import get_batch_id, get_conditions, get_path, get_timepoints
from .validation import validate_sample_manifest


def get_sample_manifest(cfg):
    """Return sample manifest if configured, else None."""
    data = cfg.setdefault("_data", {})
    cached = data.get("sample_manifest")
    if cached is not None:
        return cached

    manifest_path = cfg.get("paths", {}).get("sample_manifest")
    if not manifest_path or not os.path.isfile(manifest_path):
        return None

    manifest = validate_sample_manifest(pd.read_csv(manifest_path, sep="\t"))
    data["sample_manifest"] = manifest
    return manifest


def require_sample_manifest(cfg, required_columns=None):
    """Return a validated sample manifest and enforce extra required columns."""
    manifest = get_sample_manifest(cfg)
    if manifest is None:
        raise ValueError("sample_manifest.tsv is required for this analysis step")
    return validate_sample_manifest(manifest, required_columns=required_columns)


def _sample_manifest_row(cfg, sample_name):
    manifest = get_sample_manifest(cfg)
    if manifest is None:
        return None
    rows = manifest.loc[manifest["sample_id"] == sample_name]
    if rows.empty:
        return None
    return rows.iloc[0]


def get_sample_names(cfg):
    """Sample names for the current comparison pair."""
    manifest = get_sample_manifest(cfg)
    if manifest is not None:
        conditions = get_conditions(cfg)
        samples = {}
        all_samples = []
        sort_cols = [
            col for col in
            ["target_condition", "timepoint", "batch_id", "source_condition", "sample_id"]
            if col in manifest.columns
        ]
        ordered = manifest.sort_values(sort_cols) if sort_cols else manifest
        for cond in conditions:
            cond_samples = ordered.loc[
                ordered["target_condition"] == cond, "sample_id"
            ].tolist()
            samples[cond] = cond_samples
            all_samples.extend(cond_samples)
        samples["all"] = all_samples
        return samples

    conditions = get_conditions(cfg)
    timepoints = get_timepoints(cfg)
    pattern = cfg["experiment"]["sample_pattern"]

    samples = {}
    all_samples = []
    for cond in conditions:
        cond_samples = [
            pattern.format(condition=cond, timepoint=t) for t in timepoints
        ]
        samples[cond] = cond_samples
        all_samples.extend(cond_samples)
    samples["all"] = all_samples
    return samples


def get_all_sample_names(cfg):
    """Sample names for all conditions in the batch."""
    manifest = get_sample_manifest(cfg)
    if manifest is not None:
        samples = {}
        all_samples = []
        sort_cols = [
            col for col in
            ["target_condition", "timepoint", "batch_id", "source_condition", "sample_id"]
            if col in manifest.columns
        ]
        ordered = manifest.sort_values(sort_cols) if sort_cols else manifest
        for cond in ordered["target_condition"].drop_duplicates().tolist():
            cond_samples = ordered.loc[
                ordered["target_condition"] == cond, "sample_id"
            ].tolist()
            samples[cond] = cond_samples
            all_samples.extend(cond_samples)
        samples["all"] = all_samples
        return samples

    all_conds = cfg.get("_all_conditions", {})
    pattern = cfg["experiment"]["sample_pattern"]

    samples = {}
    all_samples = []
    for cond, tps in sorted(all_conds.items()):
        cond_samples = [
            pattern.format(condition=cond, timepoint=t) for t in tps
        ]
        samples[cond] = cond_samples
        all_samples.extend(cond_samples)
    samples["all"] = all_samples
    return samples


def get_sample_condition(cfg, sample_name):
    """Return display grouping condition for a sample."""
    row = _sample_manifest_row(cfg, sample_name)
    if row is not None:
        return str(row["target_condition"])
    return sample_name.rsplit("-", 1)[0]


def get_sample_display_name(cfg, sample_name):
    """Return human-friendly sample label."""
    row = _sample_manifest_row(cfg, sample_name)
    if row is not None:
        return str(row.get("display_label", sample_name))
    return sample_name


def get_module1_source_conditions(cfg):
    """Requested source-condition prefixes for module 1, or None for all."""
    items = cfg.get("_module1_source_conditions")
    if not items:
        return None
    return [str(item) for item in items]


def get_preprocessing_raw_dir(cfg):
    """Resolve the directory that contains raw FASTQ files for module 1."""
    raw_data_dir = Path(get_path(cfg, "raw_fastq_dir"))
    batch_id = get_batch_id(cfg)
    nested = raw_data_dir / batch_id
    return nested if nested.exists() else raw_data_dir


def sample_condition_prefix(sample_id):
    """Return the condition prefix from a sample ID like R1-18."""
    sample_id = str(sample_id)
    return sample_id.split("-", 1)[0]


def filter_sample_ids_by_conditions(sample_ids, conditions):
    """Filter sample IDs by allowed source-condition prefixes."""
    if not conditions:
        return sorted(str(item) for item in sample_ids)
    allowed = {str(item) for item in conditions}
    return sorted(
        str(sample_id) for sample_id in sample_ids
        if sample_condition_prefix(sample_id) in allowed
    )


def get_norm_method(cfg):
    """Return normalization method string ('cpm', 'fpkm', or 'tpm')."""
    return cfg.get("normalization", {}).get("method", "cpm").lower()


def get_filtered_expr(cfg, log2=True):
    """Load filtered normalized expression matrix (in-memory or from disk)."""
    method = get_norm_method(cfg)
    data = cfg.get("_data", {})
    key = f"filtered_log2{method}" if log2 else f"filtered_{method}"
    result = data.get(key)
    if result is not None:
        return result
    fname = f"filtered_log2{method}.tsv" if log2 else f"filtered_{method}.tsv"
    out_root = get_path(cfg, "output_dir")
    path = os.path.join(out_root, "01_normalization", "gene_filtering", fname)
    result = pd.read_csv(path, sep="\t", index_col=0)
    cfg.setdefault("_data", {})[key] = result
    return result


def get_filtered_counts(cfg):
    """Load filtered raw integer counts (in-memory or from disk)."""
    data = cfg.get("_data", {})
    result = data.get("filtered_counts")
    if result is not None:
        return result
    out_root = get_path(cfg, "output_dir")
    path = os.path.join(out_root, "01_normalization", "gene_filtering", "filtered_counts.tsv")
    result = pd.read_csv(path, sep="\t", index_col=0)
    cfg.setdefault("_data", {})["filtered_counts"] = result
    return result
