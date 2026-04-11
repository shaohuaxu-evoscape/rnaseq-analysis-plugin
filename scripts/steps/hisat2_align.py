"""Step 1c: HISAT2 Alignment — align clean reads to reference genome.

Runs HISAT2 → samtools sort → samtools index for each sample.
Requires: hisat2, samtools
"""

import pathlib
import re
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..core.config_data import filter_sample_ids_by_conditions, get_module1_source_conditions
from ..core.config_runtime import ensure_shared_dir as ensure_output_dir
from ..core.logger import setup_logger
from ..core.plotting import apply_style, style_axes, save_figure, plot_filename

log = setup_logger("preprocessing.hisat2_align")


def run(cfg):
    """Run HISAT2 alignment on all samples."""
    raw_cfg = cfg.get("preprocessing", {}).get("hisat2_align", {})
    if not raw_cfg.get("enabled", True):
        log.info("HISAT2 alignment disabled, skipping")
        return

    log.info("── Step 1c: HISAT2 Alignment ──")

    threads = cfg.get("preprocessing", {}).get("threads", 8)
    max_files = cfg.get("preprocessing", {}).get("test_n_samples", 0)

    # Get index prefix from previous step or config
    data = cfg.get("_data", {})
    index_prefix = data.get("hisat2_index_prefix")
    if not index_prefix:
        # Try to find it in output dir
        index_dir = pathlib.Path(cfg["paths"]["shared_output_dir"]) / "01_preprocessing" / "hisat2_index"
        ref_name = pathlib.Path(cfg["paths"]["reference_genome"]).name
        for suffix in [".gz", ".fa", ".fasta"]:
            if ref_name.endswith(suffix):
                ref_name = ref_name[:-len(suffix)]
        index_prefix = str(index_dir / ref_name)

    # Get clean fastq dir
    clean_dir = data.get("clean_fastq_dir")
    if not clean_dir:
        clean_dir = str(pathlib.Path(cfg["paths"]["shared_output_dir"]) / "01_preprocessing" / "clean_fastq")
    clean_dir = pathlib.Path(clean_dir)
    requested_sample_ids = set(data.get("requested_sample_ids", []))
    requested_conditions = get_module1_source_conditions(cfg)

    align_dir = pathlib.Path(ensure_output_dir(cfg, "01_preprocessing/alignment"))

    # Find clean R1 files
    r1_files = sorted(clean_dir.glob("*.R1.clean.fq.gz"))
    if requested_sample_ids:
        r1_files = [
            path for path in r1_files
            if path.name.split(".R1.clean.fq.gz")[0] in requested_sample_ids
        ]
    elif requested_conditions:
        allowed_ids = set(filter_sample_ids_by_conditions(
            [path.name.split(".R1.clean.fq.gz")[0] for path in r1_files],
            requested_conditions,
        ))
        r1_files = [
            path for path in r1_files
            if path.name.split(".R1.clean.fq.gz")[0] in allowed_ids
        ]
    if not r1_files:
        log.warning(f"  No clean FASTQ files in {clean_dir}")
        return

    if max_files:
        r1_files = r1_files[:max_files]
        log.info(f"  Test mode: {len(r1_files)} samples")
    else:
        log.info(f"  Found {len(r1_files)} samples to align")

    for r1 in r1_files:
        sample_id = r1.name.split(".R1.clean.fq.gz")[0]
        r2 = clean_dir / f"{sample_id}.R2.clean.fq.gz"
        if not r2.exists():
            log.warning(f"  R2 missing for {sample_id}, skipping")
            continue

        bam_out = align_dir / f"{sample_id}.sorted.bam"
        bai_out = align_dir / f"{sample_id}.sorted.bam.bai"
        summary_file = align_dir / f"{sample_id}.hisat2.log"

        if _bam_is_complete(bam_out, bai_out):
            log.info(f"  {sample_id}: BAM + BAI exist, skipping")
            continue

        if bam_out.exists() or bai_out.exists():
            log.warning(f"  {sample_id}: found incomplete alignment outputs, removing stale files and rerunning")
            for path in [bam_out, bai_out, summary_file]:
                if path.exists():
                    path.unlink()

        log.info(f"  Aligning {sample_id} ({threads} threads)...")

        # HISAT2 → samtools sort → BAM
        hisat2_cmd = [
            "hisat2",
            "-p", str(threads),
            "-x", str(index_prefix),
            "-1", str(r1),
            "-2", str(r2),
            "--summary-file", str(summary_file),
            "--dta",
        ]
        sort_cmd = [
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(bam_out),
            "-",
        ]

        p1 = subprocess.Popen(hisat2_cmd, stdout=subprocess.PIPE)
        try:
            p2 = subprocess.run(sort_cmd, stdin=p1.stdout, check=False)
        finally:
            p1.stdout.close()
            p1.wait()

        if p1.returncode != 0:
            raise RuntimeError(f"HISAT2 failed for {sample_id} (rc={p1.returncode})")
        if p2.returncode != 0:
            raise RuntimeError(f"samtools sort failed for {sample_id} (rc={p2.returncode})")

        # Index BAM
        subprocess.run(["samtools", "index", str(bam_out)], check=True)
        log.info(f"  {sample_id}: aligned + indexed")

    # Store alignment dir
    cfg.setdefault("_data", {})["alignment_dir"] = str(align_dir)

    # Generate alignment summary from HISAT2 logs
    _write_alignment_summary(align_dir, cfg)
    log.info(f"  Alignment complete: {len(r1_files)} samples")


def _bam_is_complete(bam_path, bai_path):
    """Return True when BAM exists, index exists, and BAM passes quick integrity check."""
    if not bam_path.exists() or not bai_path.exists():
        return False
    check = subprocess.run(
        ["samtools", "quickcheck", str(bam_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return check.returncode == 0


def _parse_hisat2_log(log_path):
    """Extract alignment stats from HISAT2 summary log."""
    stats = {}
    text = log_path.read_text()
    m = re.search(r"(\d+) reads; of these:", text)
    if m:
        stats["total_reads"] = int(m.group(1))
    m = re.search(r"(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", text)
    if m:
        stats["unique_mapped"] = int(m.group(1))
    m = re.search(r"(\d+) \([\d.]+%\) aligned concordantly >1 times", text)
    if m:
        stats["multi_mapped"] = int(m.group(1))
    m = re.search(r"([\d.]+)% overall alignment rate", text)
    if m:
        stats["alignment_rate"] = float(m.group(1))
    return stats


def _write_alignment_summary(align_dir, cfg):
    """Parse all HISAT2 logs and generate summary table + plot."""
    apply_style(cfg)
    logs = sorted(align_dir.glob("*.hisat2.log"))
    if not logs:
        return

    records = []
    for lp in logs:
        sample = lp.name.replace(".hisat2.log", "")
        stats = _parse_hisat2_log(lp)
        stats["sample"] = sample
        records.append(stats)

    df = pd.DataFrame(records)
    df.to_csv(align_dir / "alignment_summary.tsv", sep="\t", index=False)
    log.info(f"  Saved alignment_summary.tsv ({len(df)} samples)")

    if "alignment_rate" not in df.columns:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(range(len(df)), df["alignment_rate"], color="#4c72b0")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["sample"], rotation=45, ha="right")
    ax.set_ylabel("Alignment Rate (%)")
    ax.set_title("HISAT2 Alignment Rate")
    ax.set_ylim(0, 100)
    style_axes(ax)
    alignment_rates_plot = plot_filename(cfg, "alignment_rates.png")
    save_figure(fig, str(align_dir / "alignment_rates.png"), cfg)
    log.info(f"  Saved {alignment_rates_plot}")
