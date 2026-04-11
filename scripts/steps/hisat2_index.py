"""Step 1a: Create HISAT2 Index from reference genome.

Builds HISAT2 index with splice site and exon annotations extracted from GTF.
Requires: hisat2-build, hisat2_extract_splice_sites.py, hisat2_extract_exons.py
"""

import gzip
import pathlib
import shutil
import subprocess

from ..core.config_runtime import ensure_shared_dir as ensure_output_dir
from ..core.logger import setup_logger

log = setup_logger("preprocessing.hisat2_index")


def _run_cmd(cmd, log_file=None):
    """Run shell command, optionally logging to file."""
    log.info(f"  $ {' '.join(str(c) for c in cmd)}")
    if log_file:
        with open(log_file, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        stderr = result.stderr if not log_file else f"see {log_file}"
        raise RuntimeError(f"Command failed (rc={result.returncode}): {stderr}")
    return result


def _ensure_uncompressed(src, dest_dir):
    """Decompress .gz file if needed, return path to uncompressed file."""
    src = pathlib.Path(src)
    if src.suffix == ".gz":
        dest = pathlib.Path(dest_dir) / src.stem
        if not dest.exists():
            log.info(f"  Decompressing {src.name}...")
            with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        return dest
    return src


def run(cfg):
    """Build HISAT2 index."""
    raw_cfg = cfg.get("preprocessing", {}).get("hisat2_index", {})
    if not raw_cfg.get("enabled", True):
        log.info("HISAT2 index disabled, skipping")
        return

    log.info("── Step 1a: HISAT2 Index ──")

    ref_fasta = pathlib.Path(cfg["paths"]["reference_genome"])
    ref_gtf = pathlib.Path(cfg["paths"]["reference_gtf"])
    threads = cfg.get("preprocessing", {}).get("threads", 8)

    out_dir = ensure_output_dir(cfg, "01_preprocessing/hisat2_index")
    out_dir = pathlib.Path(out_dir)

    # Derive index prefix from fasta name
    prefix_name = ref_fasta.name
    for suffix in [".gz", ".fa", ".fasta"]:
        if prefix_name.endswith(suffix):
            prefix_name = prefix_name[:-len(suffix)]
    index_prefix = out_dir / prefix_name

    # Check if index already exists
    if (out_dir / f"{prefix_name}.1.ht2").exists():
        log.info(f"  Index exists at {index_prefix}, skipping")
        cfg.setdefault("_data", {})["hisat2_index_prefix"] = str(index_prefix)
        return

    # Decompress if needed
    temp_fasta = _ensure_uncompressed(ref_fasta, out_dir)
    temp_gtf = _ensure_uncompressed(ref_gtf, out_dir)

    # Extract splice sites
    splice_sites = out_dir / f"{prefix_name}.splice_sites.txt"
    log.info("  Extracting splice sites...")
    with open(splice_sites, "w") as f:
        subprocess.run(
            ["hisat2_extract_splice_sites.py", str(temp_gtf)],
            stdout=f, check=True
        )

    # Extract exons
    exons = out_dir / f"{prefix_name}.exons.txt"
    log.info("  Extracting exons...")
    with open(exons, "w") as f:
        subprocess.run(
            ["hisat2_extract_exons.py", str(temp_gtf)],
            stdout=f, check=True
        )

    # Build index
    log.info(f"  Building HISAT2 index ({threads} threads)...")
    _run_cmd([
        "hisat2-build",
        "-p", str(threads),
        "--ss", str(splice_sites),
        "--exon", str(exons),
        str(temp_fasta),
        str(index_prefix),
    ])

    # Cleanup temp decompressed files
    if ref_fasta.suffix == ".gz" and temp_fasta.exists():
        temp_fasta.unlink()
    if ref_gtf.suffix == ".gz" and temp_gtf.exists():
        temp_gtf.unlink()

    cfg.setdefault("_data", {})["hisat2_index_prefix"] = str(index_prefix)
    log.info(f"  Index built: {index_prefix}")
