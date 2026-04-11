"""Step 1b: Fastp Quality Control and Trimming.

Runs fastp on paired-end FASTQ files for adapter trimming,
quality filtering, and generates QC reports.
Requires: fastp
"""

import json
import pathlib
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..core.config_data import get_module1_source_conditions, get_preprocessing_raw_dir
from ..core.config_runtime import ensure_shared_dir as ensure_output_dir
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, plot_filename

log = setup_logger("preprocessing.fastp")


def run(cfg):
    """Run fastp on all samples."""
    raw_cfg = cfg.get("preprocessing", {}).get("fastp", {})
    if not raw_cfg.get("enabled", True):
        log.info("Fastp disabled, skipping")
        return

    log.info("── Step 1b: Fastp QC ──")

    threads = cfg.get("preprocessing", {}).get("threads", 8)
    max_files = cfg.get("preprocessing", {}).get("test_n_samples", 0)
    requested_conditions = get_module1_source_conditions(cfg)

    ensure_output_dir(cfg, "01_preprocessing/fastp")
    clean_dir = pathlib.Path(ensure_output_dir(cfg, "01_preprocessing/clean_fastq"))
    report_dir = pathlib.Path(ensure_output_dir(cfg, "01_preprocessing/fastp_reports"))

    raw_dir = get_preprocessing_raw_dir(cfg)
    r1_files = _requested_r1_files(raw_dir, requested_conditions)
    if not r1_files:
        log.warning(f"  No .R1.fq.gz files found in {raw_dir}")
        return

    if max_files:
        r1_files = r1_files[:max_files]
        log.info(f"  Test mode: {len(r1_files)} samples")
    else:
        log.info(f"  Found {len(r1_files)} samples")

    requested_sample_ids = [r1.name.split(".R1.fq.gz")[0] for r1 in r1_files]
    processed = []
    for r1 in r1_files:
        sample_id = r1.name.split(".R1.fq.gz")[0]
        r2 = raw_dir / f"{sample_id}.R2.fq.gz"
        if not r2.exists():
            log.warning(f"  R2 missing for {sample_id}, skipping")
            continue

        out_r1 = clean_dir / f"{sample_id}.R1.clean.fq.gz"
        out_r2 = clean_dir / f"{sample_id}.R2.clean.fq.gz"
        json_report = report_dir / f"{sample_id}.fastp.json"
        html_report = report_dir / f"{sample_id}.fastp.html"

        # Skip if output exists
        if out_r1.exists() and out_r2.exists():
            log.info(f"  {sample_id}: output exists, skipping")
            processed.append({"sample": sample_id, "json": json_report})
            continue

        log.info(f"  Processing {sample_id}...")
        fastp_cfg = cfg.get("preprocessing", {}).get("fastp", {})
        cmd = [
            "fastp",
            "-i", str(r1), "-I", str(r2),
            "-o", str(out_r1), "-O", str(out_r2),
            "--thread", str(threads),
            "--cut_front", "--cut_tail",
            "--cut_window_size", str(fastp_cfg.get("cut_window_size", 4)),
            "--cut_mean_quality", str(fastp_cfg.get("cut_mean_quality", 20)),
            "--qualified_quality_phred", str(fastp_cfg.get("qualified_quality_phred", 20)),
            "--unqualified_percent_limit", str(fastp_cfg.get("unqualified_percent_limit", 40)),
            "--length_required", str(fastp_cfg.get("length_required", 50)),
            "--html", str(html_report),
            "--json", str(json_report),
            "--report_title", sample_id,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"fastp failed for {sample_id} (rc={result.returncode}):\n{result.stderr}"
            )
        processed.append({"sample": sample_id, "json": json_report})

    # Generate metrics summary + QC plot
    _write_metrics(report_dir)
    _plot_qc(report_dir, cfg)
    log.info(f"  Fastp complete: {len(processed)} samples")

    # Store clean fastq dir for downstream
    data = cfg.setdefault("_data", {})
    data["clean_fastq_dir"] = str(clean_dir)
    data["requested_sample_ids"] = requested_sample_ids


def _requested_r1_files(raw_dir, requested_conditions):
    if not requested_conditions:
        return sorted(raw_dir.glob("*.R1.fq.gz"))

    r1_files = []
    missing = []
    for condition in requested_conditions:
        matches = sorted(raw_dir.glob(f"{condition}-*.R1.fq.gz"))
        if not matches:
            missing.append(condition)
        r1_files.extend(matches)

    if missing:
        raise ValueError(
            "Module 1 could not find FASTQ files for requested source conditions "
            f"{', '.join(missing)} in {raw_dir}"
        )
    unique = {path.name: path for path in r1_files}
    return [unique[name] for name in sorted(unique)]


def _report_items(report_dir):
    items = []
    for json_path in sorted(report_dir.glob("*.fastp.json")):
        items.append({
            "sample": json_path.name.replace(".fastp.json", ""),
            "json": json_path,
        })
    return items


def _write_metrics(report_dir):
    """Write fastp metrics summary TSV."""
    processed = _report_items(report_dir)
    metrics_file = report_dir / "fastp_metrics.tsv"
    header = [
        "Sample", "RawReads(M)", "RawBases(G)", "CleanReads(M)",
        "CleanBases(G)", "ValidBases(%)", "Q30(%)", "GC(%)"
    ]

    with open(metrics_file, "w") as f:
        f.write("\t".join(header) + "\n")
        for item in processed:
            json_path = item["json"]
            if not json_path.exists():
                continue
            try:
                with open(json_path, "r") as jf:
                    data = json.load(jf)
                summary = data.get("summary", {})
                before = summary.get("before_filtering", {})
                after = summary.get("after_filtering", {})

                raw_reads = float(before.get("total_reads", 0))
                raw_bases = float(before.get("total_bases", 0))
                clean_reads = float(after.get("total_reads", 0))
                clean_bases = float(after.get("total_bases", 0))
                q30 = float(after.get("q30_rate", 0))
                gc = float(after.get("gc_content", 0))
                valid_pct = (clean_bases / raw_bases * 100) if raw_bases > 0 else 0

                row = [
                    item["sample"],
                    f"{raw_reads/1e6:.2f}", f"{raw_bases/1e9:.2f}",
                    f"{clean_reads/1e6:.2f}", f"{clean_bases/1e9:.2f}",
                    f"{valid_pct:.2f}", f"{q30*100:.2f}", f"{gc*100:.2f}",
                ]
                f.write("\t".join(row) + "\n")
            except Exception as e:
                log.warning(f"  Failed to parse metrics for {item['sample']}: {e}")

    log.info(f"  Saved fastp_metrics.tsv")


def _plot_qc(report_dir, cfg):
    """Generate QC reads bar plot from fastp reports."""
    apply_style(cfg)
    processed = _report_items(report_dir)
    records = []
    for item in processed:
        json_path = item["json"]
        if not json_path.exists():
            continue
        try:
            with open(json_path, "r") as jf:
                data = json.load(jf)
            summary = data.get("summary", {})
            before = summary.get("before_filtering", {})
            after = summary.get("after_filtering", {})
            records.append({
                "sample": item["sample"],
                "raw_reads": float(before.get("total_reads", 0)),
                "clean_reads": float(after.get("total_reads", 0)),
            })
        except Exception as e:
            log.warning(f"  Failed to parse {item['sample']} for QC plot: {e}")
            continue

    if not records:
        return

    df = pd.DataFrame(records)
    fig, ax = plt.subplots(figsize=(10, 5))
    x = range(len(df))
    width = 0.35
    ax.bar([i - width / 2 for i in x], df["raw_reads"] / 1e6, width,
           label="Raw", color="#aaaaaa")
    ax.bar([i + width / 2 for i in x], df["clean_reads"] / 1e6, width,
           label="Clean", color="#4c72b0")
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=45, ha="right")
    ax.set_ylabel("Reads (millions)")
    ax.set_title("Sequencing Depth — Before/After QC")
    ax.legend()
    style_axes(ax)
    qc_reads_plot = plot_filename(cfg, "qc_reads.png")
    save_figure(fig, str(report_dir / "qc_reads.png"), cfg)
    log.info(f"  Saved {qc_reads_plot}")
