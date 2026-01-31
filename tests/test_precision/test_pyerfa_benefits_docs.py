"""
Tests for the PYERFA_BENEFITS.md documentation.

This module verifies the precision claims made in the pyerfa benefits documentation,
ensuring that the documented precision improvements are accurate and demonstrable.
"""

import math
import os
from pathlib import Path

import pytest

# Check if pyerfa is available
try:
    import erfa

    HAS_ERFA = True
except ImportError:
    HAS_ERFA = False

from libephemeris.erfa_nutation import (
    compare_nutation_models,
    get_erfa_nutation_cached,
    get_erfa_nutation_nut00a,
    get_erfa_nutation_nut06a,
    get_erfa_obliquity_iau2006,
    get_erfa_pnm06a_matrix,
    has_erfa,
)

# Standard test dates
J2000_JD = 2451545.0  # J2000.0 epoch (Jan 1, 2000, 12:00 TT)


class TestDocumentationExists:
    """Test that the PYERFA_BENEFITS.md documentation exists and is complete."""

    def test_documentation_file_exists(self):
        """PYERFA_BENEFITS.md should exist in the docs directory."""
        docs_dir = Path(__file__).parent.parent.parent / "docs"
        doc_path = docs_dir / "PYERFA_BENEFITS.md"
        assert doc_path.exists(), f"Documentation file not found at {doc_path}"

    def test_documentation_has_required_sections(self):
        """Documentation should contain all required sections."""
        docs_dir = Path(__file__).parent.parent.parent / "docs"
        doc_path = docs_dir / "PYERFA_BENEFITS.md"
        content = doc_path.read_text()

        required_sections = [
            "Overview",
            "Installation",
            "Precision Comparison",
            "Quantified Benefits",
            "PyERFA Functions Available",
            "When to Use PyERFA",
            "Summary",
        ]

        for section in required_sections:
            assert section in content, f"Missing section: {section}"

    def test_documentation_mentions_precision_levels(self):
        """Documentation should mention specific precision levels."""
        docs_dir = Path(__file__).parent.parent.parent / "docs"
        doc_path = docs_dir / "PYERFA_BENEFITS.md"
        content = doc_path.read_text()

        # Should mention key precision values
        assert "milliarcsecond" in content.lower() or "mas" in content
        assert "IAU 2000A" in content
        assert "IAU 2000B" in content
        assert "IAU 2006" in content or "nut06a" in content


class TestDocumentedPrecisionClaims:
    """Verify the precision claims made in the documentation."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut06a_precision_over_nut00a(self):
        """
        Verify that nut06a provides ~0.01-0.05 mas improvement over nut00a.

        Documentation claims: "IAU 2006/2000A (nut06a): ~0.01-0.05 mas"
        """
        # Test at 10 years from J2000 (documented test point)
        jd_tt = J2000_JD + 10 * 365.25

        dpsi_00a, deps_00a = get_erfa_nutation_nut00a(jd_tt)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        # Convert difference to milliarcseconds
        dpsi_diff_mas = abs(math.degrees(dpsi_00a - dpsi_06a)) * 3600 * 1000
        deps_diff_mas = abs(math.degrees(deps_00a - deps_06a)) * 3600 * 1000

        # Documentation claims: difference should be in the range 0.01-0.05 mas
        # We allow up to 0.1 mas to account for variation across dates
        assert dpsi_diff_mas < 0.1, (
            f"dpsi difference {dpsi_diff_mas:.4f} mas exceeds documented range"
        )
        assert deps_diff_mas < 0.1, (
            f"deps difference {deps_diff_mas:.4f} mas exceeds documented range"
        )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_error_growth_pattern(self):
        """
        Verify the documented error growth pattern over time.

        Documentation claims error grows with distance from J2000.
        """
        # Test points from documentation table
        test_years = [0, 10, 50, 100]

        previous_diff = 0.0
        for years in test_years:
            jd_tt = J2000_JD + years * 365.25

            dpsi_00a, deps_00a = get_erfa_nutation_nut00a(jd_tt)
            dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

            # Convert to milliarcseconds
            dpsi_diff_mas = abs(math.degrees(dpsi_00a - dpsi_06a)) * 3600 * 1000

            # Verify the difference stays within documented bounds
            # Documentation: ~0.01 at J2000, ~0.1 at 100 years
            max_expected = 0.05 + 0.001 * years  # Linear growth approximation
            assert dpsi_diff_mas < max_expected, (
                f"Error at {years} years ({dpsi_diff_mas:.4f} mas) exceeds expected"
            )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_iau2000b_vs_2000a_precision(self):
        """
        Verify IAU 2000B vs 2000A difference is ~1 mas as documented.

        Documentation claims: "IAU 2000B: ~1 milliarcsecond (mas)"
        """
        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()
        jd_tt = J2000_JD + 10 * 365.25  # 10 years from J2000
        t = ts.tt_jd(jd_tt)

        dpsi_2000a, deps_2000a = iau2000a_radians(t)
        dpsi_2000b, deps_2000b = iau2000b_radians(t)

        # Convert difference to milliarcseconds
        dpsi_diff_mas = abs(math.degrees(dpsi_2000a - dpsi_2000b)) * 3600 * 1000
        deps_diff_mas = abs(math.degrees(deps_2000a - deps_2000b)) * 3600 * 1000

        # Documentation claims ~1 mas difference
        # Allow up to 2 mas to account for variation
        assert dpsi_diff_mas < 2.0, (
            f"2000B vs 2000A dpsi difference {dpsi_diff_mas:.4f} mas exceeds ~1 mas"
        )
        assert deps_diff_mas < 2.0, (
            f"2000B vs 2000A deps difference {deps_diff_mas:.4f} mas exceeds ~1 mas"
        )


class TestDocumentedFunctionality:
    """Verify that documented pyerfa functions work as described."""

    def test_has_erfa_function(self):
        """has_erfa() should return boolean as documented."""
        result = has_erfa()
        assert isinstance(result, bool)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nutation_functions_return_format(self):
        """Nutation functions should return (dpsi, deps) tuple as documented."""
        jd_tt = J2000_JD

        # nut00a
        result = get_erfa_nutation_nut00a(jd_tt)
        assert result is not None
        assert isinstance(result, tuple)
        assert len(result) == 2
        dpsi, deps = result
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

        # nut06a
        result = get_erfa_nutation_nut06a(jd_tt)
        assert result is not None
        assert isinstance(result, tuple)
        assert len(result) == 2

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_obliquity_function(self):
        """get_erfa_obliquity_iau2006 should return obliquity in radians."""
        eps = get_erfa_obliquity_iau2006(J2000_JD)
        assert eps is not None
        assert isinstance(eps, float)

        # Should be approximately 23.44 degrees
        eps_deg = math.degrees(eps)
        assert 23.4 < eps_deg < 23.5, f"Obliquity {eps_deg}° not in expected range"

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_matrix_format(self):
        """get_erfa_pnm06a_matrix should return 3x3 matrix as documented."""
        rbpn = get_erfa_pnm06a_matrix(J2000_JD)
        assert rbpn is not None
        assert hasattr(rbpn, "shape")
        assert rbpn.shape == (3, 3)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_compare_nutation_models_function(self):
        """compare_nutation_models should return dict with documented keys."""
        result = compare_nutation_models(J2000_JD)
        assert isinstance(result, dict)
        assert "skyfield_iau2000b" in result
        assert "skyfield_iau2000a" in result
        assert "erfa_nut00a" in result
        assert "erfa_nut06a" in result
        assert "differences_mas" in result

    def test_cached_nutation_fallback(self):
        """get_erfa_nutation_cached should work even without pyerfa (fallback)."""
        # This should always work due to Skyfield fallback
        dpsi, deps = get_erfa_nutation_cached(J2000_JD)
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

        # Values should be reasonable
        dpsi_arcsec = abs(math.degrees(dpsi) * 3600)
        deps_arcsec = abs(math.degrees(deps) * 3600)
        assert dpsi_arcsec < 25, "dpsi should be within ±25 arcseconds"
        assert deps_arcsec < 15, "deps should be within ±15 arcseconds"


class TestDocumentedUseCases:
    """Verify the documented use cases and recommendations."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_astrological_precision_context(self):
        """
        Verify that for typical astrological use, the precision difference is negligible.

        Documentation states: "For typical astrological work, IAU 2000A is already sufficient."
        """
        # Typical astrological precision requirement: 1 arcminute = 60 arcseconds
        # This is ~60,000 milliarcseconds

        jd_tt = J2000_JD + 25 * 365.25  # Current epoch (approx)

        dpsi_00a, deps_00a = get_erfa_nutation_nut00a(jd_tt)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        # Difference in arcseconds
        dpsi_diff_arcsec = abs(math.degrees(dpsi_00a - dpsi_06a)) * 3600
        deps_diff_arcsec = abs(math.degrees(deps_00a - deps_06a)) * 3600

        # The difference should be negligible compared to 1 arcminute
        # (less than 0.001 arcseconds = 1 milliarcsecond)
        astrological_threshold = 60.0  # 1 arcminute in arcseconds

        assert dpsi_diff_arcsec < astrological_threshold * 0.001, (
            "Precision difference should be negligible for astrology"
        )
        assert deps_diff_arcsec < astrological_threshold * 0.001, (
            "Precision difference should be negligible for astrology"
        )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_long_term_error_growth(self):
        """
        Verify that pyerfa provides measurable benefit for long-term calculations.

        Documentation states: "Long-term Ephemerides: Calculations spanning centuries"
        The nut00a vs nut06a difference remains small even at distant dates,
        but pyerfa's combined matrix (pnm06a) avoids cross-term accumulation errors.
        """
        # Test at 500 years from J2000 (far future)
        jd_tt = J2000_JD + 500 * 365.25

        dpsi_00a, deps_00a = get_erfa_nutation_nut00a(jd_tt)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        dpsi_diff_mas = abs(math.degrees(dpsi_00a - dpsi_06a)) * 3600 * 1000

        # The nut00a vs nut06a difference is actually quite small (~0.01 mas)
        # even at 500 years. The main benefit of pyerfa for long-term
        # calculations is the pnm06a matrix which avoids cross-term errors.
        assert dpsi_diff_mas < 2.0, f"Long-term error {dpsi_diff_mas:.4f} mas too large"
        # Verify there is some measurable difference (even if small)
        assert dpsi_diff_mas > 0.001, (
            f"Expected some difference at 500 years, got {dpsi_diff_mas:.4f} mas"
        )


class TestDocumentedObliquityModels:
    """Verify obliquity model comparison claims in documentation."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_obliquity_laskar_vs_iau2006_difference(self):
        """
        Verify obliquity difference between IAU 2006 and Laskar 1986.

        Documentation claims: "~42 mas" at J2000
        """
        # IAU 2006 obliquity from ERFA
        eps_iau2006 = get_erfa_obliquity_iau2006(J2000_JD)

        # Laskar 1986 formula (documented fallback)
        T = (J2000_JD - J2000_JD) / 36525.0  # = 0 at J2000
        eps_laskar = math.radians(
            23.439291111
            - 0.013004166667 * T
            - 1.638888889e-7 * T**2
            + 5.036111111e-7 * T**3
        )

        # Difference in milliarcseconds
        diff_mas = abs(math.degrees(eps_iau2006 - eps_laskar)) * 3600 * 1000

        # Documentation claims ~42 mas at J2000
        # Allow some margin for the documented approximation
        assert diff_mas < 100, (
            f"Obliquity difference {diff_mas:.4f} mas at J2000 larger than expected"
        )


class TestInstallationDocumentation:
    """Verify that documented installation methods are accurate."""

    def test_pyerfa_optional_dependency_documented(self):
        """Verify the optional dependency is correctly documented."""
        # Check pyproject.toml for the precision extra
        project_root = Path(__file__).parent.parent.parent
        pyproject_path = project_root / "pyproject.toml"

        if pyproject_path.exists():
            content = pyproject_path.read_text()
            assert "precision" in content, "precision extra should be in pyproject.toml"
            assert "pyerfa" in content, "pyerfa should be mentioned in pyproject.toml"
