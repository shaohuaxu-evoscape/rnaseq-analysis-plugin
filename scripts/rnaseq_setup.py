#!/usr/bin/env python3
"""Interactive project setup — generates configs/analysis_case.yaml.

This script assists with creating a new analysis case configuration.
It scans local data, detects available samples, and guides through
condition definition.

Usage:
    python rnaseq_setup.py
    python rnaseq_setup.py --template configs/analysis_case.template.yaml
"""
import argparse
import os
import sys

# Ensure project root is on path — walk up to find pyproject.toml
_d = os.path.abspath(os.getcwd())
for _ in range(6):
    if os.path.isfile(os.path.join(_d, "pyproject.toml")):
        break
    _d = os.path.dirname(_d)
_project_root = _d
sys.path.insert(0, _project_root)
os.chdir(_project_root)


def scan_local_data():
    """Scan for available data files."""
    info = {}
    # Reference genomes
    for root, _, files in os.walk("inputs/ref"):
        for f in files:
            if f.endswith(".fa") or f.endswith(".fasta"):
                info.setdefault("reference_genome", []).append(os.path.join(root, f))
            if f.endswith(".gtf"):
                info.setdefault("reference_gtf", []).append(os.path.join(root, f))
    # Gene counts
    for root, _, files in os.walk("results/shared"):
        for f in files:
            if f == "gene_counts.tsv":
                info.setdefault("gene_counts", []).append(os.path.join(root, f))
    # FASTQ directories
    for d in ["inputs/fastq", "inputs/raw"]:
        if os.path.isdir(d):
            info.setdefault("fastq_dirs", []).append(d)
    return info


def detect_conditions(gene_counts_path):
    """Read gene count matrix header to detect conditions and timepoints."""
    with open(gene_counts_path) as f:
        header = f.readline().strip().split("\t")[1:]
    conditions = sorted(set(col.rsplit("-", 1)[0] for col in header))
    timepoints = sorted(set(int(col.rsplit("-", 1)[1]) for col in header))
    return conditions, timepoints


def main():
    parser = argparse.ArgumentParser(description="RNA-seq project setup wizard")
    parser.add_argument("--template", default="configs/analysis_case.template.yaml",
                        help="Path to template config")
    args = parser.parse_args()

    print("=" * 50)
    print("RNA-seq Analysis — Project Setup Wizard")
    print("=" * 50)

    data = scan_local_data()

    print("\nLocal data scan:")
    for key, paths in data.items():
        for p in paths:
            print(f"  [{key}] {p}")
    if not data:
        print("  No data files found in inputs/ or results/shared/")

    if data.get("gene_counts"):
        gc_path = data["gene_counts"][0]
        conditions, timepoints = detect_conditions(gc_path)
        print(f"\nDetected conditions: {conditions}")
        print(f"Detected timepoints: {timepoints}")

    print(f"\nTemplate: {args.template}")
    print("\nTo complete setup, use the /rnaseq-analysis skill which will")
    print("guide you through the interactive configuration wizard.")


if __name__ == "__main__":
    main()
