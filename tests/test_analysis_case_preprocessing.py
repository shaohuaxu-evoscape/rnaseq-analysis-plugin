import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from scripts.core import prepare_analysis_case as pac
from scripts.core.config_data import require_sample_manifest
from scripts.core.config_runtime import load_config
from scripts.core.pipeline_state import get_state_path
from scripts.core.plotting import plot_filename, save_figure
from scripts.steps.fastp import _requested_r1_files
from scripts.core.runner import run_pipeline
from scripts.steps.pca import _valid_plot_pairs
from scripts.steps.gene_clustering import _get_clustering_config
from scripts.steps.temporal_causality import _plot_two_panel
from scripts.core.validation import validate_analysis_case_config, validate_pipeline_config, validate_sample_manifest


def _write_counts(path: Path, columns):
    df = pd.DataFrame(
        {
            "Gene": ["YALI1_A", "YALI1_B"],
            **columns,
        }
    )
    df.to_csv(path, sep="\t", index=False)


def _write_counts_rows(path: Path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


class AnalysisCasePreprocessingTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_case(self, raw_fastq_dir: Path, output_root: Path):
        case_path = self.tmp / "analysis_case.yaml"
        case_cfg = {
            "run_name": "analysis_case_test",
            "output_dir": str(output_root),
            "batches": {
                "batch1": {
                    "gene_counts": "",
                    "raw_fastq_dir": str(raw_fastq_dir),
                }
            },
            "target_conditions": {
                "condition1": [{"batch": "batch1", "condition": "R1"}],
                "condition2": [{"batch": "batch1", "condition": "R2"}],
            },
            "project": {
                "name": "DS-002",
                "organism": "Yarrowia lipolytica",
                "strain": "A316",
            },
            "experiment": {"timepoints": [18, 24]},
            "paths": {
                "reference_genome": "inputs/ref/A316.v1.fa",
                "reference_gtf": "inputs/ref/A316.v1.gtf",
            },
        }
        with open(case_path, "w") as f:
            yaml.safe_dump(case_cfg, f, sort_keys=False)
        return case_path

    def test_prepare_analysis_case_uses_shared_batch_name_in_manifests(self):
        raw_dir = self.tmp / "20260313"
        raw_dir.mkdir()
        output_root = self.tmp / "out"
        counts_path = self.tmp / "shared_counts.tsv"
        _write_counts(
            counts_path,
            {
                "R1-18": [10, 20],
                "R1-24": [11, 21],
                "R2-18": [30, 40],
                "R2-24": [31, 41],
            },
        )
        case_path = self._write_case(raw_dir, output_root)

        with mock.patch.object(pac, "_default_counts_path", return_value=counts_path):
            result = pac.prepare_analysis_case(case_path)

        sample_manifest = pd.read_csv(result["sample_manifest"], sep="\t")
        batch_manifest = pd.read_csv(result["batch_manifest"], sep="\t")

        self.assertEqual(set(sample_manifest["batch_id"].astype(str)), {"20260313"})
        self.assertTrue(
            sample_manifest["sample_id"].str.contains("condition1_20260313_R1-18").any()
        )
        self.assertEqual(batch_manifest.loc[0, "yaml_batch_id"], "batch1")
        self.assertEqual(str(batch_manifest.loc[0, "batch_id"]), "20260313")
        self.assertEqual(batch_manifest.loc[0, "requested_conditions"], "R1,R2")
        with open(result["analysis_config"], "r") as f:
            generated_cfg = yaml.safe_load(f)
        self.assertTrue(str(generated_cfg["paths"]["shared_output_dir"]).endswith("results/shared"))

    def test_ensure_batch_counts_reruns_when_shared_counts_missing_requested_condition(self):
        raw_dir = self.tmp / "20260313"
        raw_dir.mkdir()
        case_path = self._write_case(raw_dir, self.tmp / "out")
        counts_path = self.tmp / "shared_counts.tsv"
        _write_counts(
            counts_path,
            {
                "R1-18": [10, 20],
                "R2-18": [30, 40],
            },
        )

        def fake_run_pipeline(config_path, steps, dry_run, overrides, allow_hidden):
            self.assertEqual(steps, "preprocessing")
            self.assertFalse(dry_run)
            self.assertTrue(allow_hidden)
            self.assertEqual(overrides["batch_id"], "20260313")
            self.assertEqual(overrides["shared_batch_name"], "20260313")
            self.assertEqual(overrides["shared_output_dir"], "results/shared/20260313")
            self.assertEqual(overrides["module1_source_conditions"], ["R1", "R5"])
            _write_counts(
                counts_path,
                {
                    "R1-18": [10, 20],
                    "R2-18": [30, 40],
                    "R5-18": [50, 60],
                },
            )
            return {
                "status": "complete",
                "results": {
                    "0a": {"status": "ok"},
                    "0b": {"status": "ok"},
                    "0c": {"status": "ok"},
                    "0d": {"status": "ok"},
                },
            }

        with mock.patch.object(pac, "_default_counts_path", return_value=counts_path):
            with mock.patch.object(pac, "run_pipeline", side_effect=fake_run_pipeline) as mocked:
                out_path, mode = pac._ensure_batch_counts(
                    yaml_batch_id="batch1",
                    shared_batch_name="20260313",
                    batch_cfg={"gene_counts": "", "raw_fastq_dir": str(raw_dir)},
                    config_path=case_path,
                    raw_overrides={},
                    requested_conditions=["R1", "R5"],
                )

        self.assertEqual(out_path, counts_path)
        self.assertEqual(mode, "generated")
        self.assertEqual(mocked.call_count, 1)

    def test_default_output_paths_resolve_under_case_project_root(self):
        project_root = self.tmp / "project"
        case_dir = project_root / "configs"
        case_dir.mkdir(parents=True)
        case_path = case_dir / "analysis_case.yaml"
        case_path.write_text("run_name: my_run\n")

        run_name, output_root, shared_root, analysis_config_path = pac._resolve_analysis_output_paths(
            {"run_name": "my_run"},
            case_path,
        )

        self.assertEqual(run_name, "my_run")
        self.assertEqual(output_root, (project_root / "results" / "my_run").resolve())
        self.assertEqual(shared_root, (project_root / "results" / "shared").resolve())
        self.assertEqual(analysis_config_path, output_root / "analysis_config.yaml")

    def test_default_counts_path_resolves_under_case_project_root(self):
        case_path = self.tmp / "project" / "configs" / "analysis_case.yaml"
        case_path.parent.mkdir(parents=True)
        case_path.write_text("batches: {}\n")

        counts_path = pac._default_counts_path("batch_x", case_path)

        self.assertEqual(
            counts_path,
            (self.tmp / "project" / "results" / "shared" / "batch_x" / "01_preprocessing" / "gene_counts" / "gene_counts.tsv").resolve(),
        )

    def test_requested_r1_files_filters_to_requested_conditions(self):
        raw_dir = self.tmp / "raw"
        raw_dir.mkdir()
        for name in [
            "R1-18.R1.fq.gz",
            "R1-24.R1.fq.gz",
            "R2-18.R1.fq.gz",
            "R5-18.R1.fq.gz",
        ]:
            (raw_dir / name).touch()

        files = _requested_r1_files(raw_dir, ["R1", "R2"])
        self.assertEqual(
            [path.name for path in files],
            ["R1-18.R1.fq.gz", "R1-24.R1.fq.gz", "R2-18.R1.fq.gz"],
        )

        with self.assertRaisesRegex(ValueError, "R6"):
            _requested_r1_files(raw_dir, ["R1", "R6"])

    def test_plot_filename_and_saved_extension_follow_configured_format(self):
        cfg = {"plot": {"format": "pdf", "dpi": 150}}
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])

        out_path = save_figure(fig, self.tmp / "figure.png", cfg)

        self.assertEqual(plot_filename(cfg, "figure.png"), "figure.pdf")
        self.assertTrue(str(out_path).endswith("figure.pdf"))
        self.assertTrue((self.tmp / "figure.pdf").exists())
        self.assertFalse((self.tmp / "figure.png").exists())

    def test_valid_pca_plot_pairs_skip_unavailable_components(self):
        self.assertEqual(
            _valid_plot_pairs(2, [[1, 2], [1, 3], [2, 2], [0, 1], [1, 2]]),
            [(1, 2)],
        )

    def test_get_clustering_config_prefers_cluster_analysis_and_supports_legacy_fallback(self):
        canonical_cfg = {
            "cluster_analysis": {
                "gene_clustering": {"enabled": True, "top_n_genes": 10, "k_range": [2, 3], "random_state": 1}
            },
            "sample_analysis": {
                "clustering": {"enabled": False}
            },
        }
        legacy_cfg = {
            "sample_analysis": {
                "clustering": {"enabled": True, "top_n_genes": 20, "k_range": [3, 4], "random_state": 2}
            }
        }

        self.assertTrue(_get_clustering_config(canonical_cfg)["enabled"])
        self.assertEqual(_get_clustering_config(canonical_cfg)["top_n_genes"], 10)
        self.assertEqual(_get_clustering_config(legacy_cfg)["top_n_genes"], 20)

    def test_temporal_plot_handles_partial_module_set(self):
        cfg = {
            "plot": {"format": "png", "dpi": 72},
            "experiment": {"conditions": ["condition1", "condition2"]},
        }
        out_dir = self.tmp / "temporal"
        out_dir.mkdir()

        named_modules = {
            "Sterol": (np.array([0.1, 0.2, 0.25]), "#377EB8"),
            "MVA": (np.array([-0.2, -0.1, 0.05]), "#E78AC3"),
        }
        lag_results = [
            {
                "A": "Sterol",
                "B": "MVA",
                "forward": 0.7,
                "reverse": -0.1,
                "color": "#377EB8",
            }
        ]

        _plot_two_panel([18, 24, 48], named_modules, lag_results, out_dir, cfg)

        self.assertTrue((out_dir / "fig_temporal_causality.png").exists())

    def test_validate_analysis_case_rejects_unknown_batch_reference(self):
        cfg = {
            "batches": {"batch1": {}},
            "target_conditions": {
                "condition1": [{"batch": "missing_batch", "condition": "R1"}],
            },
        }
        with self.assertRaisesRegex(ValueError, "unknown batch"):
            validate_analysis_case_config(cfg)

    def test_validate_pipeline_config_rejects_duplicate_conditions(self):
        cfg = {
            "project": {"batch_id": "run1"},
            "experiment": {"conditions": ["condition1", "condition1"]},
            "paths": {},
        }
        with self.assertRaisesRegex(ValueError, "duplicate condition"):
            validate_pipeline_config(cfg)

    def test_validate_sample_manifest_rejects_duplicate_sample_ids(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["s1", "s1"],
                "target_condition": ["condition1", "condition2"],
                "timepoint": [18, 24],
            }
        )
        with self.assertRaisesRegex(ValueError, "duplicate sample_id"):
            validate_sample_manifest(manifest)

    def test_require_sample_manifest_normalizes_merged_condition_column(self):
        manifest_path = self.tmp / "sample_manifest.tsv"
        pd.DataFrame(
            {
                "sample_id": ["s1", "s2"],
                "merged_condition": ["condition1", "condition2"],
                "timepoint": [18, 24],
                "batch_id": ["b1", "b1"],
            }
        ).to_csv(manifest_path, sep="\t", index=False)

        cfg = {
            "project": {"batch_id": "run1"},
            "experiment": {"conditions": ["condition1", "condition2"], "timepoints": [18, 24]},
            "paths": {"sample_manifest": str(manifest_path)},
        }

        manifest = require_sample_manifest(cfg, required_columns={"batch_id"})
        self.assertIn("target_condition", manifest.columns)
        self.assertEqual(manifest["target_condition"].tolist(), ["condition1", "condition2"])

    def test_end_to_end_prepare_then_run_1a_and_3a(self):
        counts_path = self.tmp / "counts.tsv"
        output_root = self.tmp / "results"
        _write_counts_rows(
            counts_path,
            [
                {"Gene": "YALI1_A", "R1-18": 100, "R1-24": 120, "R2-18": 300, "R2-24": 360},
                {"Gene": "YALI1_B", "R1-18": 200, "R1-24": 220, "R2-18": 210, "R2-24": 230},
                {"Gene": "YALI1_C", "R1-18": 50, "R1-24": 55, "R2-18": 20, "R2-24": 22},
            ],
        )

        case_path = self.tmp / "analysis_case_e2e.yaml"
        case_cfg = {
            "run_name": "analysis_case_e2e",
            "output_dir": str(output_root),
            "batches": {
                "batch1": {
                    "gene_counts": str(counts_path),
                }
            },
            "target_conditions": {
                "condition1": [{"batch": "batch1", "condition": "R1"}],
                "condition2": [{"batch": "batch1", "condition": "R2"}],
            },
            "project": {
                "name": "DS-002",
                "organism": "Yarrowia lipolytica",
                "strain": "A316",
            },
            "experiment": {"timepoints": [18, 24]},
            "normalization": {
                "method": "cpm",
                "pseudocount": 1,
                "gene_filtering": {"min_expr": 1.0, "min_samples": 1},
            },
            "differential_analysis": {
                "de_screening": {
                    "enabled": True,
                    "fdr_threshold": 0.05,
                    "log2fc_threshold": 0.585,
                }
            },
            "plot": {"dpi": 72, "format": "png"},
        }
        with open(case_path, "w") as f:
            yaml.safe_dump(case_cfg, f, sort_keys=False)

        prepared = pac.prepare_analysis_case(case_path)
        result = run_pipeline(config_path=prepared["analysis_config"], steps="1a,3a")
        resolved_cfg = load_config(prepared["analysis_config"])
        pair_output_dir = Path(resolved_cfg["paths"]["output_dir"])

        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["steps"], ["1a", "3a"])
        self.assertEqual(result["results"]["1a"]["status"], "ok")
        self.assertEqual(result["results"]["3a"]["status"], "ok")
        self.assertEqual(prepared["timepoints"], [18, 24])
        self.assertTrue((pair_output_dir / "01_normalization" / "gene_filtering" / "filtered_counts.tsv").exists())
        self.assertTrue((pair_output_dir / "03_differential_analysis" / "de_screening" / "de_results_all.tsv").exists())
        self.assertTrue((pair_output_dir / "03_differential_analysis" / "de_screening" / "de_genes.tsv").exists())
        self.assertTrue((output_root / "documents" / "condition1-condition2.md").exists())

    def test_resume_skips_completed_steps_using_pipeline_manifest(self):
        counts_path = self.tmp / "counts_resume.tsv"
        output_root = self.tmp / "results_resume"
        _write_counts_rows(
            counts_path,
            [
                {"Gene": "YALI1_A", "R1-18": 80, "R1-24": 100, "R2-18": 240, "R2-24": 300},
                {"Gene": "YALI1_B", "R1-18": 150, "R1-24": 180, "R2-18": 140, "R2-24": 175},
                {"Gene": "YALI1_C", "R1-18": 40, "R1-24": 42, "R2-18": 12, "R2-24": 10},
            ],
        )
        case_path = self.tmp / "analysis_case_resume.yaml"
        case_cfg = {
            "run_name": "analysis_case_resume",
            "output_dir": str(output_root),
            "batches": {"batch1": {"gene_counts": str(counts_path)}},
            "target_conditions": {
                "condition1": [{"batch": "batch1", "condition": "R1"}],
                "condition2": [{"batch": "batch1", "condition": "R2"}],
            },
            "project": {
                "name": "DS-002",
                "organism": "Yarrowia lipolytica",
                "strain": "A316",
            },
            "experiment": {"timepoints": [18, 24]},
            "normalization": {
                "method": "cpm",
                "pseudocount": 1,
                "gene_filtering": {"min_expr": 1.0, "min_samples": 1},
            },
            "differential_analysis": {
                "de_screening": {
                    "enabled": True,
                    "fdr_threshold": 0.05,
                    "log2fc_threshold": 0.585,
                }
            },
            "plot": {"dpi": 72, "format": "png"},
        }
        with open(case_path, "w") as f:
            yaml.safe_dump(case_cfg, f, sort_keys=False)

        prepared = pac.prepare_analysis_case(case_path)
        first = run_pipeline(config_path=prepared["analysis_config"], steps="1a,3a", resume=False)
        second = run_pipeline(config_path=prepared["analysis_config"], steps="1a,3a", resume=True)
        resolved_cfg = load_config(prepared["analysis_config"])
        manifest_path = Path(get_state_path(resolved_cfg))

        self.assertEqual(first["results"]["1a"]["status"], "ok")
        self.assertEqual(first["results"]["3a"]["status"], "ok")
        self.assertEqual(second["results"]["1a"]["status"], "skipped")
        self.assertEqual(second["results"]["3a"]["status"], "skipped")
        self.assertEqual(second["results"]["1a"]["reason"], "resume_checkpoint")
        self.assertTrue(manifest_path.exists())
        manifest = yaml.safe_load(manifest_path.read_text())
        self.assertEqual(manifest["steps"]["1a"]["status"], "ok")
        self.assertEqual(manifest["steps"]["3a"]["status"], "ok")


if __name__ == "__main__":
    unittest.main()
