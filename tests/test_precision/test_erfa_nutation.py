"""
Tests for ERFA Nutation Integration Evaluation.

This module tests the pyerfa integration and evaluates precision improvements
from using erfa.nut00a(), erfa.nut06a(), and erfa.pnm06a().

Test categories:
1. Basic functionality tests - verify functions work correctly
2. Precision comparison tests - compare ERFA vs Skyfield implementations
3. Edge case tests - extreme dates, cached values
4. Integration tests - verify the module works with existing libephemeris code
"""

import math
from datetime import datetime

import pytest

# Check if pyerfa is available
try:
    import erfa

    HAS_ERFA = True
except ImportError:
    HAS_ERFA = False

from libephemeris.erfa_nutation import (
    clear_erfa_nutation_cache,
    compare_nutation_models,
    get_erfa_nutation_cache_info,
    get_erfa_nutation_cached,
    get_erfa_nutation_nut00a,
    get_erfa_nutation_nut06a,
    get_erfa_obliquity_iau2006,
    get_erfa_pnm06a_matrix,
    has_erfa,
)

# Standard test dates
J2000_JD = 2451545.0  # J2000.0 epoch (Jan 1, 2000, 12:00 TT)


class TestERFAAvailability:
    """Test ERFA availability detection."""

    def test_has_erfa_returns_bool(self):
        """has_erfa() should return a boolean."""
        result = has_erfa()
        assert isinstance(result, bool)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_has_erfa_true_when_installed(self):
        """has_erfa() should return True when pyerfa is installed."""
        assert has_erfa() is True


class TestERFANutation:
    """Test ERFA nutation functions."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut00a_returns_tuple(self):
        """get_erfa_nutation_nut00a should return a tuple of two floats."""
        result = get_erfa_nutation_nut00a(J2000_JD)
        assert result is not None
        assert isinstance(result, tuple)
        assert len(result) == 2
        dpsi, deps = result
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut06a_returns_tuple(self):
        """get_erfa_nutation_nut06a should return a tuple of two floats."""
        result = get_erfa_nutation_nut06a(J2000_JD)
        assert result is not None
        assert isinstance(result, tuple)
        assert len(result) == 2
        dpsi, deps = result
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_returns_matrix(self):
        """get_erfa_pnm06a_matrix should return a 3x3 matrix."""
        result = get_erfa_pnm06a_matrix(J2000_JD)
        assert result is not None
        # pyerfa returns numpy array
        assert hasattr(result, "shape")
        assert result.shape == (3, 3)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_obliquity_iau2006_returns_float(self):
        """get_erfa_obliquity_iau2006 should return a float in radians."""
        result = get_erfa_obliquity_iau2006(J2000_JD)
        assert result is not None
        assert isinstance(result, float)
        # Obliquity at J2000 should be approximately 23.44 degrees
        obliquity_deg = math.degrees(result)
        assert 23.4 < obliquity_deg < 23.5

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nutation_values_reasonable(self):
        """Nutation values should be in reasonable ranges."""
        result = get_erfa_nutation_nut06a(J2000_JD)
        dpsi, deps = result

        # Nutation in longitude is typically ±20 arcseconds
        dpsi_arcsec = math.degrees(dpsi) * 3600
        assert -25 < dpsi_arcsec < 25

        # Nutation in obliquity is typically ±10 arcseconds
        deps_arcsec = math.degrees(deps) * 3600
        assert -15 < deps_arcsec < 15


class TestCachedNutation:
    """Test cached nutation functionality."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_erfa_nutation_cache()

    def test_cached_nutation_returns_tuple(self):
        """get_erfa_nutation_cached should always return a tuple."""
        # This should work even without ERFA (uses Skyfield fallback)
        result = get_erfa_nutation_cached(J2000_JD)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_cache_works(self):
        """Caching should store results."""
        clear_erfa_nutation_cache()
        info_before = get_erfa_nutation_cache_info()

        # First call
        get_erfa_nutation_cached(J2000_JD)
        info_after_first = get_erfa_nutation_cache_info()

        # Second call with same value
        get_erfa_nutation_cached(J2000_JD)
        info_after_second = get_erfa_nutation_cache_info()

        # After first call, should have 1 miss
        assert info_after_first["misses"] == info_before["misses"] + 1

        # After second call, should have 1 hit
        assert info_after_second["hits"] == info_after_first["hits"] + 1

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_model_selection_nut00a(self):
        """Should be able to select nut00a model."""
        result = get_erfa_nutation_cached(J2000_JD, model="nut00a")
        assert isinstance(result, tuple)
        assert len(result) == 2

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_model_selection_nut06a(self):
        """Should be able to select nut06a model."""
        result = get_erfa_nutation_cached(J2000_JD, model="nut06a")
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestPrecisionComparison:
    """
    Test precision comparison between Skyfield and ERFA implementations.

    These tests verify the precision differences between:
    - Skyfield's IAU 2000A (iau2000a_radians)
    - ERFA's nut00a (should be identical to Skyfield's IAU 2000A)
    - ERFA's nut06a (IAU 2000A with IAU 2006 adjustments - most accurate)
    """

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_skyfield_vs_erfa_nut00a_equivalence(self):
        """
        Skyfield's iau2000a and ERFA's nut00a should produce nearly identical results.

        Both implement the same IAU 2000A nutation model (MHB2000).
        Any differences should be due to floating point implementation details,
        not algorithmic differences.
        """
        from skyfield.nutationlib import iau2000a_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()

        test_dates = [
            J2000_JD,  # J2000.0
            J2000_JD + 3652.5,  # ~10 years later
            J2000_JD - 3652.5,  # ~10 years earlier
            J2000_JD + 36525,  # 100 years later
            J2000_JD - 36525,  # 100 years earlier
        ]

        for jd_tt in test_dates:
            t = ts.tt_jd(jd_tt)
            dpsi_skyfield, deps_skyfield = iau2000a_radians(t)
            dpsi_erfa, deps_erfa = get_erfa_nutation_nut00a(jd_tt)

            # Difference in microarcseconds
            dpsi_diff_uas = abs(math.degrees(dpsi_skyfield - dpsi_erfa)) * 3600 * 1e6
            deps_diff_uas = abs(math.degrees(deps_skyfield - deps_erfa)) * 3600 * 1e6

            # Should agree to within 1 microarcsecond (numerical precision)
            assert dpsi_diff_uas < 1.0, (
                f"dpsi diff at JD {jd_tt}: {dpsi_diff_uas:.3f} µas"
            )
            assert deps_diff_uas < 1.0, (
                f"deps diff at JD {jd_tt}: {deps_diff_uas:.3f} µas"
            )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut00a_vs_nut06a_differences(self):
        """
        ERFA's nut00a and nut06a should differ by small amounts.

        nut06a includes:
        - J2 secular variation correction
        - IAU 2006 obliquity frame correction

        These differences should be on the order of 0.01-0.1 milliarcseconds
        for dates within a century of J2000.
        """
        test_dates = [
            J2000_JD,  # J2000.0
            J2000_JD + 3652.5,  # ~10 years later
            J2000_JD + 36525,  # 100 years later
        ]

        for jd_tt in test_dates:
            dpsi_00a, deps_00a = get_erfa_nutation_nut00a(jd_tt)
            dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

            # Difference in milliarcseconds
            dpsi_diff_mas = abs(math.degrees(dpsi_00a - dpsi_06a)) * 3600 * 1000
            deps_diff_mas = abs(math.degrees(deps_00a - deps_06a)) * 3600 * 1000

            # Should differ by less than 1 milliarcsecond for recent dates
            years_from_j2000 = (jd_tt - J2000_JD) / 365.25
            max_expected_mas = 0.1 + 0.005 * abs(years_from_j2000)  # Grows with time

            assert dpsi_diff_mas < max_expected_mas, (
                f"dpsi diff at JD {jd_tt} ({years_from_j2000:.0f} years from J2000): "
                f"{dpsi_diff_mas:.4f} mas (expected < {max_expected_mas:.4f} mas)"
            )
            assert deps_diff_mas < max_expected_mas, (
                f"deps diff at JD {jd_tt} ({years_from_j2000:.0f} years from J2000): "
                f"{deps_diff_mas:.4f} mas (expected < {max_expected_mas:.4f} mas)"
            )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_iau2000b_vs_iau2000a_differences(self):
        """
        IAU 2000B and IAU 2000A should differ by ~1 milliarcsecond.

        IAU 2000B is a truncated version of IAU 2000A with only 77 terms
        (vs 1365 terms). The difference indicates the precision improvement
        from using the full model.
        """
        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()

        test_dates = [
            J2000_JD,
            J2000_JD + 365.25 * 10,  # 10 years
            J2000_JD + 365.25 * 50,  # 50 years
        ]

        for jd_tt in test_dates:
            t = ts.tt_jd(jd_tt)
            dpsi_2000a, deps_2000a = iau2000a_radians(t)
            dpsi_2000b, deps_2000b = iau2000b_radians(t)

            # Difference in milliarcseconds
            dpsi_diff_mas = abs(math.degrees(dpsi_2000a - dpsi_2000b)) * 3600 * 1000
            deps_diff_mas = abs(math.degrees(deps_2000a - deps_2000b)) * 3600 * 1000

            # IAU 2000B should be accurate to ~1 mas, so differences can be 0-1 mas
            assert dpsi_diff_mas < 2.0, (
                f"dpsi diff at JD {jd_tt}: {dpsi_diff_mas:.4f} mas"
            )
            assert deps_diff_mas < 2.0, (
                f"deps diff at JD {jd_tt}: {deps_diff_mas:.4f} mas"
            )


class TestCompareNutationModels:
    """Test the comparison function."""

    def test_compare_function_returns_dict(self):
        """compare_nutation_models should return a dictionary."""
        result = compare_nutation_models(J2000_JD)
        assert isinstance(result, dict)
        assert "jd_tt" in result
        assert result["jd_tt"] == J2000_JD

    def test_compare_function_has_skyfield_results(self):
        """Compare function should always include Skyfield results."""
        result = compare_nutation_models(J2000_JD)
        assert "skyfield_iau2000b" in result
        assert "skyfield_iau2000a" in result
        assert "dpsi_arcsec" in result["skyfield_iau2000b"]
        assert "deps_arcsec" in result["skyfield_iau2000b"]

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_compare_function_has_erfa_results(self):
        """Compare function should include ERFA results when available."""
        result = compare_nutation_models(J2000_JD)
        assert "erfa_nut00a" in result
        assert "erfa_nut06a" in result
        assert "differences_mas" in result


class TestPNM06aMatrix:
    """Test the precession-nutation matrix functionality."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_is_rotation_matrix(self):
        """pnm06a should return a proper rotation matrix."""
        import numpy as np

        rbpn = get_erfa_pnm06a_matrix(J2000_JD)

        # A rotation matrix should be orthogonal: R @ R.T = I
        identity = rbpn @ rbpn.T
        assert np.allclose(identity, np.eye(3), atol=1e-15)

        # Determinant should be +1
        det = np.linalg.det(rbpn)
        assert np.isclose(det, 1.0, atol=1e-15)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_at_j2000_near_identity(self):
        """At J2000.0, the matrix should be close to identity (small corrections)."""
        import numpy as np

        rbpn = get_erfa_pnm06a_matrix(J2000_JD)

        # The off-diagonal elements represent small rotations
        # For J2000, they should be very small but non-zero (frame bias)
        off_diag = [
            rbpn[0, 1],
            rbpn[0, 2],
            rbpn[1, 0],
            rbpn[1, 2],
            rbpn[2, 0],
            rbpn[2, 1],
        ]

        # At J2000, the matrix includes:
        # - Frame bias (~17 mas in x, ~7 mas in y, ~80 mas in z)
        # - Nutation (~10-20 arcseconds, i.e., ~1e-4 radians)
        # Off-diagonal elements represent rotations, so they should be < 1e-4
        for val in off_diag:
            assert abs(val) < 1e-3, f"Off-diagonal element {val} too large for J2000"

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_transforms_icrs_to_true_equator(self):
        """
        pnm06a should correctly transform ICRS coordinates to true equator of date.

        The celestial pole (ICRS z-axis) should be transformed to the
        true pole position (with precession and nutation).
        """
        import numpy as np

        # 50 years from J2000 to see significant precession
        jd = J2000_JD + 50 * 365.25
        rbpn = get_erfa_pnm06a_matrix(jd)

        # ICRS north pole
        icrs_pole = np.array([0.0, 0.0, 1.0])

        # Transform to true equator of date
        true_pole = rbpn @ icrs_pole

        # The pole should have moved slightly due to precession
        # Precession is ~50.3" per year in longitude, so after 50 years
        # the pole has moved significantly
        assert not np.allclose(true_pole, icrs_pole, atol=1e-6)

        # But magnitude should still be 1
        assert np.isclose(np.linalg.norm(true_pole), 1.0, atol=1e-15)


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_distant_past_date(self):
        """Functions should work for dates in the distant past."""
        # 2000 years ago
        jd_tt = J2000_JD - 2000 * 365.25
        result = get_erfa_nutation_nut06a(jd_tt)
        assert result is not None
        dpsi, deps = result
        assert math.isfinite(dpsi)
        assert math.isfinite(deps)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_distant_future_date(self):
        """Functions should work for dates in the distant future."""
        # 2000 years in future
        jd_tt = J2000_JD + 2000 * 365.25
        result = get_erfa_nutation_nut06a(jd_tt)
        assert result is not None
        dpsi, deps = result
        assert math.isfinite(dpsi)
        assert math.isfinite(deps)

    def test_invalid_model_raises_error(self):
        """Invalid model name should raise ValueError."""
        clear_erfa_nutation_cache()
        # This depends on pyerfa being available to reach the error
        if HAS_ERFA:
            with pytest.raises(ValueError, match="Unknown nutation model"):
                get_erfa_nutation_cached(J2000_JD, model="invalid_model")


class TestPrecisionImprovementEvaluation:
    """
    Evaluation tests to quantify precision improvements from ERFA integration.

    These tests document and verify the expected precision improvements.
    """

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_precision_improvement_summary(self):
        """
        Document the precision hierarchy of nutation models.

        Expected precision (in milliarcseconds):
        - IAU 2000B: ~1 mas
        - IAU 2000A: ~0.1 mas (10x improvement)
        - IAU 2006/2000A (nut06a): ~0.01-0.05 mas (2-10x further improvement)
        """
        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()
        jd_tt = J2000_JD + 10 * 365.25  # 10 years from J2000

        t = ts.tt_jd(jd_tt)

        # Get all nutation values
        dpsi_2000b, deps_2000b = iau2000b_radians(t)
        dpsi_2000a, deps_2000a = iau2000a_radians(t)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        # Convert to milliarcseconds for easier interpretation
        def to_mas(rad):
            return math.degrees(rad) * 3600 * 1000

        # Compute differences
        diff_2000b_to_2000a_dpsi = abs(to_mas(dpsi_2000b) - to_mas(dpsi_2000a))
        diff_2000b_to_2000a_deps = abs(to_mas(deps_2000b) - to_mas(deps_2000a))

        diff_2000a_to_06a_dpsi = abs(to_mas(dpsi_2000a) - to_mas(dpsi_06a))
        diff_2000a_to_06a_deps = abs(to_mas(deps_2000a) - to_mas(deps_06a))

        # Document the precision hierarchy
        print(f"\n=== Nutation Precision Evaluation at JD {jd_tt} ===")
        print("\nIAU 2000B (77 terms, ~1 mas precision):")
        print(f"  dpsi = {to_mas(dpsi_2000b):.4f} mas")
        print(f"  deps = {to_mas(deps_2000b):.4f} mas")

        print("\nIAU 2000A (1365 terms, ~0.1 mas precision):")
        print(f"  dpsi = {to_mas(dpsi_2000a):.4f} mas")
        print(f"  deps = {to_mas(deps_2000a):.4f} mas")
        print(f"  Improvement over 2000B: {diff_2000b_to_2000a_dpsi:.4f} mas (dpsi)")
        print(f"                         {diff_2000b_to_2000a_deps:.4f} mas (deps)")

        print("\nIAU 2006/2000A (nut06a, ~0.01-0.05 mas precision):")
        print(f"  dpsi = {to_mas(dpsi_06a):.4f} mas")
        print(f"  deps = {to_mas(deps_06a):.4f} mas")
        print(f"  Improvement over 2000A: {diff_2000a_to_06a_dpsi:.4f} mas (dpsi)")
        print(f"                         {diff_2000a_to_06a_deps:.4f} mas (deps)")

        # Verify expected precision levels
        # 2000B vs 2000A should differ by up to ~1 mas
        assert diff_2000b_to_2000a_dpsi < 2.0
        assert diff_2000b_to_2000a_deps < 2.0

        # 2000A vs 06a should differ by up to ~0.1 mas for dates near J2000
        assert diff_2000a_to_06a_dpsi < 0.5
        assert diff_2000a_to_06a_deps < 0.5
