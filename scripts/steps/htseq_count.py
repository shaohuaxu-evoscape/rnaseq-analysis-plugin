"""Step 1d: HTSeq-count — count reads per gene feature.

Runs htseq-count on sorted BAM files and combines into a gene count matrix.
Requires: htseq-count
"""

import csv
import gzip
import multiprocessing
import pathlib
import shutil
import subprocess
from typing import Dict, Tuple
import pandas as pd

from ..core.config_data import filter_sample_ids_by_conditions, get_module1_source_conditions
from ..core.config_runtime import ensure_shared_dir as ensure_output_dir
from ..core.logger import setup_logger

log = setup_logger("preprocessing.htseq_count")


def _count_one_bam(args: Tuple[str, str]) -> Tuple[str, Dict[str, int]]:
    """Run htseq-count on a single BAM file and return parsed counts."""
    bam_path_str, gtf_path_str = args
    bam = pathlib.Path(bam_path_str)
    sample_name = bam.name.replace(".sorted.bam", "")

    cmd = [
        "htseq-count",
        "--format=bam",
        "--order=pos",
        "--mode=union",
        "--stranded=no",
        "--type=exon",
        "--idattr=gene_id",
        str(bam),
        gtf_path_str,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    counts = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        gene, count = parts
        if gene.startswith("__"):
            continue
        counts[gene] = int(count)

    return sample_name, counts


def run(cfg):
    """Run htseq-count on all samples and combine."""
    raw_cfg = cfg.get("preprocessing", {}).get("htseq_count", {})
    if not raw_cfg.get("enabled", True):
        log.info("HTSeq-count disabled, skipping")
        return

    log.info("── Step 1d: HTSeq-count ──")

    ref_gtf = pathlib.Path(cfg["paths"]["reference_gtf"])
    max_files = cfg.get("preprocessing", {}).get("test_n_samples", 0)
    threads = cfg.get("preprocessing", {}).get("threads", 8)

    # Get alignment dir
    data = cfg.get("_data", {})
    align_dir = data.get("alignment_dir")
    if not align_dir:
        align_dir = str(pathlib.Path(cfg["paths"]["shared_output_dir"]) / "01_preprocessing" / "alignment")
    align_dir = pathlib.Path(align_dir)
    requested_sample_ids = set(data.get("requested_sample_ids", []))
    requested_conditions = get_module1_source_conditions(cfg)

    counts_dir = pathlib.Path(ensure_output_dir(cfg, "01_preprocessing/gene_counts"))
    sample_counts_dir = counts_dir / "sample_counts"
    sample_counts_dir.mkdir(parents=True, exist_ok=True)
    combined_tsv = counts_dir / "gene_counts.tsv"

    # Find BAM files
    bams = sorted(align_dir.glob("*.sorted.bam"))
    if requested_sample_ids:
        bams = [
            bam for bam in bams
            if bam.name.replace(".sorted.bam", "") in requested_sample_ids
        ]
    elif requested_conditions:
        allowed_ids = set(filter_sample_ids_by_conditions(
            [bam.name.replace(".sorted.bam", "") for bam in bams],
            requested_conditions,
        ))
        bams = [
            bam for bam in bams
            if bam.name.replace(".sorted.bam", "") in allowed_ids
        ]
    if not bams:
        log.warning(f"  No BAM files in {align_dir}")
        return

    if max_files:
        bams = bams[:max_files]

    existing_matrix = _load_existing_matrix(combined_tsv)
    existing_samples = set(existing_matrix.columns) if existing_matrix is not None else set()
    task_bams = []
    for bam in bams:
        sample_name = bam.name.replace(".sorted.bam", "")
        sample_tsv = sample_counts_dir / f"{sample_name}.tsv"
        if sample_tsv.exists() or sample_name in existing_samples:
            log.info(f"  {sample_name}: counts exist, skipping")
            continue
        task_bams.append(bam)

    log.info(
        f"  Counting features for {len(task_bams)} new sample(s); "
        f"{len(existing_samples)} sample(s) already present in combined counts"
    )

    # Decompress GTF if needed
    temp_gtf = None
    gtf_to_use = ref_gtf
    if ref_gtf.suffix == ".gz":
        temp_gtf = counts_dir / ref_gtf.stem
        if not temp_gtf.exists():
            with gzip.open(ref_gtf, "rb") as fin, open(temp_gtf, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        gtf_to_use = temp_gtf

    if task_bams:
        n_workers = min(len(task_bams), max(1, threads), multiprocessing.cpu_count(), 8)
        log.info(
            f"  Running htseq-count with {n_workers} parallel worker(s) "
            f"on {len(task_bams)} new sample(s)..."
        )
        task_args = [(str(bam), str(gtf_to_use)) for bam in task_bams]
        try:
            with multiprocessing.Pool(n_workers) as pool:
                results = pool.map(_count_one_bam, task_args)
        except Exception:
            log.error("  htseq-count pool failed; see traceback above")
            raise
        finally:
            if temp_gtf and temp_gtf.exists():
                temp_gtf.unlink()

        for sample_name, counts in results:
            _write_sample_counts(sample_counts_dir / f"{sample_name}.tsv", counts)
    elif temp_gtf and temp_gtf.exists():
        temp_gtf.unlink()

    sample_counts = _load_sample_counts_dir(sample_counts_dir)
    if existing_matrix is not None:
        for sample_name in existing_matrix.columns:
            if sample_name not in sample_counts:
                sample_counts[sample_name] = existing_matrix[sample_name].astype(int).to_dict()

    if not sample_counts:
        log.warning("  No sample counts available to combine")
        return

    all_genes = set()
    for counts in sample_counts.values():
        all_genes.update(counts.keys())

    sorted_genes = sorted(all_genes)
    sorted_samples = sorted(sample_counts.keys())

    with open(combined_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene"] + sorted_samples)
        for gene in sorted_genes:
            row = [gene] + [sample_counts[sample].get(gene, 0) for sample in sorted_samples]
            writer.writerow(row)

    log.info(f"  Gene counts: {len(sorted_genes)} genes x {len(sorted_samples)} samples")
    log.info(f"  Saved {combined_tsv}")


def _write_sample_counts(path: pathlib.Path, counts: Dict[str, int]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", "Count"])
        for gene in sorted(counts):
            writer.writerow([gene, counts[gene]])


def _load_existing_matrix(path: pathlib.Path):
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", index_col=0)


def _load_sample_counts_dir(sample_counts_dir: pathlib.Path) -> Dict[str, Dict[str, int]]:
    sample_counts = {}
    for sample_path in sorted(sample_counts_dir.glob("*.tsv")):
        sample_name = sample_path.stem
        df = pd.read_csv(sample_path, sep="\t")
        sample_counts[sample_name] = dict(zip(
            df["Gene"].astype(str), df["Count"].astype(int)
        ))
    return sample_counts
