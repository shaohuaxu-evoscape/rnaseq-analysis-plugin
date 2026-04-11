import unittest

from scripts.core.runner import _parse_steps, _resolve_step_plan


class RunnerStepParsingTests(unittest.TestCase):
    def test_public_all_exposes_only_visible_steps(self):
        self.assertEqual(
            _parse_steps("all"),
            ["1a", "2a", "2b", "2c", "2d", "3a", "3b", "3c", "3d", "4a", "4b", "4c", "4d", "5a", "5b"],
        )

    def test_hidden_preprocessing_only_available_in_internal_mode(self):
        self.assertEqual(_parse_steps("preprocessing", allow_hidden=True), ["0a", "0b", "0c", "0d"])
        with self.assertRaisesRegex(ValueError, "automatic and hidden"):
            _parse_steps("preprocessing")

    def test_module_ranges(self):
        self.assertEqual(_parse_steps("1-2"), ["1a", "2a", "2b", "2c", "2d"])
        self.assertEqual(
            _parse_steps("1-5"),
            ["1a", "2a", "2b", "2c", "2d", "3a", "3b", "3c", "3d", "4a", "4b", "4c", "4d", "5a", "5b"],
        )

    def test_step_range(self):
        self.assertEqual(_parse_steps("3a-3d"), ["3a", "3b", "3c", "3d"])
        self.assertEqual(_parse_steps("4c-5b"), ["4c", "4d", "5a", "5b"])

    def test_individual_steps(self):
        self.assertEqual(_parse_steps("1a,3b,5a"), ["1a", "3b", "5a"])

    def test_module_name(self):
        self.assertEqual(_parse_steps("normalization"), ["1a"])
        self.assertEqual(_parse_steps("cluster_analysis"), ["5a", "5b"])

    def test_invalid_step_rejected(self):
        with self.assertRaises(ValueError):
            _parse_steps("9z")

    def test_deduplication(self):
        self.assertEqual(_parse_steps("1a,1a,1a"), ["1a"])

    def test_dependencies_are_auto_included_in_execution_plan(self):
        plan, auto_added = _resolve_step_plan(["3b"])
        self.assertEqual(plan, ["1a", "3a", "3b"])
        self.assertEqual(auto_added, ["1a", "3a"])

    def test_execution_plan_preserves_registry_order_after_dependency_expansion(self):
        plan, auto_added = _resolve_step_plan(["5b", "2b"])
        self.assertEqual(plan, ["1a", "2b", "5a", "5b"])
        self.assertEqual(auto_added, ["1a", "5a"])


if __name__ == "__main__":
    unittest.main()
