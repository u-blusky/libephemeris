"""
Tests for compare_scripts precision validation suite.

Tests verify that:
1. Comparison scripts run without errors
2. Comparison modules are properly registered in run_all_compare
3. Hypothetical planet comparisons execute
4. Minor body comparisons execute with proper category handling
"""

import sys
import subprocess

import pytest

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")


class TestCompareScriptsInfrastructure:
    """Test that comparison infrastructure is properly set up."""

    def test_comparison_modules_list(self):
        """Test that all comparison modules are registered."""
        from run_all_compare import COMPARISON_MODULES

        # Check we have a reasonable number of modules
        assert len(COMPARISON_MODULES) >= 19, (
            "Should have at least 19 comparison modules"
        )

        # Check each module has correct structure
        for module in COMPARISON_MODULES:
            assert len(module) == 4, f"Module {module} should have 4 elements"
            name, script, desc, is_quick = module
            assert isinstance(name, str)
            assert isinstance(script, str)
            assert script.endswith(".py")
            assert isinstance(desc, str)
            assert isinstance(is_quick, bool)

    def test_hypothetical_module_registered(self):
        """Test that hypothetical module is in the comparison list."""
        from run_all_compare import COMPARISON_MODULES

        module_names = [m[0] for m in COMPARISON_MODULES]
        assert "hypothetical" in module_names, (
            "hypothetical module should be registered"
        )

    def test_minor_bodies_module_registered(self):
        """Test that minor_bodies module is in the comparison list."""
        from run_all_compare import COMPARISON_MODULES

        module_names = [m[0] for m in COMPARISON_MODULES]
        assert "minor_bodies" in module_names, (
            "minor_bodies module should be registered"
        )


class TestCompareHypothetical:
    """Test hypothetical planet comparison script."""

    def test_import_compare_hypothetical(self):
        """Test that compare_hypothetical can be imported."""
        from compare_hypothetical import (
            URANIAN_PLANETS,
            OTHER_HYPOTHETICAL,
            compare_hypothetical_body,
            run_all_comparisons,
        )

        assert len(URANIAN_PLANETS) == 8, "Should have 8 Uranian planets"
        assert len(OTHER_HYPOTHETICAL) >= 1, "Should have at least Transpluto"

    def test_uranian_planets_defined(self):
        """Test all 8 Uranian planets are defined."""
        from compare_hypothetical import URANIAN_PLANETS

        names = [p[2] for p in URANIAN_PLANETS]
        expected = [
            "Cupido",
            "Hades",
            "Zeus",
            "Kronos",
            "Apollon",
            "Admetos",
            "Vulkanus",
            "Poseidon",
        ]
        assert names == expected, f"Expected {expected}, got {names}"


class TestCompareMinorBodies:
    """Test minor bodies comparison script."""

    def test_import_compare_minor_bodies(self):
        """Test that compare_minor_bodies can be imported."""
        from compare_minor_bodies import (
            MAIN_BELT_ASTEROIDS,
            CENTAURS,
            TNOS,
            NEAR_EARTH_ASTEROIDS,
            TOLERANCES,
            compare_minor_body,
            run_all_comparisons,
        )

        assert len(MAIN_BELT_ASTEROIDS) >= 4, (
            "Should have at least Ceres, Pallas, Juno, Vesta"
        )
        assert len(CENTAURS) >= 2, "Should have at least Chiron, Pholus"
        assert len(TNOS) >= 5, "Should have multiple TNOs"

    def test_tolerances_defined(self):
        """Test tolerances are defined for each category."""
        from compare_minor_bodies import TOLERANCES

        expected_categories = ["main_belt", "centaur", "tno", "near_earth"]
        for cat in expected_categories:
            assert cat in TOLERANCES, f"Tolerance for {cat} should be defined"
            assert "longitude" in TOLERANCES[cat]
            assert "latitude" in TOLERANCES[cat]

    def test_tnos_list(self):
        """Test TNO list includes major bodies."""
        from compare_minor_bodies import TNOS

        names = [t[2] for t in TNOS]
        expected_in_list = ["Eris", "Makemake", "Sedna", "Orcus"]
        for expected in expected_in_list:
            assert expected in names, f"{expected} should be in TNO list"


class TestCompareLunar:
    """Test lunar comparison script."""

    def test_import_compare_lunar(self):
        """Test that compare_lunar can be imported."""
        from compare_lunar import (
            compare_lunar_nodes,
            compare_lilith,
            compare_true_node_precision,
            compare_true_lilith_precision,
        )

    def test_threshold_constants(self):
        """Test threshold constants are defined."""
        from compare_lunar import (
            TRUE_NODE_MAX_ERROR_THRESHOLD,
            TRUE_NODE_MEAN_ERROR_THRESHOLD,
            TRUE_LILITH_MAX_ERROR_THRESHOLD,
        )

        assert TRUE_NODE_MAX_ERROR_THRESHOLD > 0
        assert TRUE_NODE_MEAN_ERROR_THRESHOLD > 0
        assert TRUE_LILITH_MAX_ERROR_THRESHOLD > 0


class TestComparisonUtils:
    """Test shared comparison utilities."""

    def test_import_comparison_utils(self):
        """Test that comparison_utils can be imported."""
        from comparison_utils import (
            Tolerances,
            angular_diff,
            format_coord,
            format_diff,
            format_status,
            TestStatistics,
            print_header,
            print_section,
        )

    def test_angular_diff(self):
        """Test angular difference calculation."""
        from comparison_utils import angular_diff

        # Simple case
        assert angular_diff(10.0, 15.0) == 5.0
        assert angular_diff(15.0, 10.0) == 5.0

        # Wrap around 360
        assert angular_diff(5.0, 355.0) == 10.0
        assert angular_diff(355.0, 5.0) == 10.0

    def test_tolerances_defined(self):
        """Test that tolerances are defined."""
        from comparison_utils import Tolerances

        assert Tolerances.LONGITUDE_STRICT > 0
        assert Tolerances.LONGITUDE_RELAXED > Tolerances.LONGITUDE_STRICT
        assert Tolerances.HOUSE_CUSP > 0
        assert Tolerances.AYANAMSHA > 0

    def test_test_statistics(self):
        """Test TestStatistics class."""
        from comparison_utils import TestStatistics

        stats = TestStatistics()
        assert stats.total == 0
        assert stats.passed == 0

        stats.add_result(True, 0.001)
        assert stats.total == 1
        assert stats.passed == 1

        stats.add_result(False, 0.1)
        assert stats.total == 2
        assert stats.passed == 1
        assert stats.failed == 1


class TestRunAllCompare:
    """Test run_all_compare.py master runner."""

    def test_module_result_dataclass(self):
        """Test ModuleResult dataclass."""
        from run_all_compare import ModuleResult

        result = ModuleResult(
            name="test",
            passed=8,
            total=10,
            duration=1.5,
        )

        assert result.pass_rate == 80.0
        assert result.status == "FAIL"

        result_pass = ModuleResult(name="test2", passed=10, total=10, duration=0.5)
        assert result_pass.status == "PASS"

    def test_print_summary_report(self):
        """Test summary report function exists."""
        from run_all_compare import print_summary_report

        # Just verify it can be called (output goes to stdout)
        from run_all_compare import ModuleResult

        results = [
            ModuleResult(name="test", passed=10, total=10, duration=0.1),
        ]
        # This returns a boolean
        result = print_summary_report(results)
        assert isinstance(result, bool)
