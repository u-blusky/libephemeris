"""
Tests for Meeus polynomial validity range warnings in lunar module.

The mean lunar node and mean Lilith calculations use polynomial approximations
from Meeus "Astronomical Algorithms" that are optimized for dates near J2000.0.
These tests verify that warnings are issued for dates outside the optimal range.
"""

import pytest
import warnings
from libephemeris.lunar import (
    calc_mean_lunar_node,
    calc_mean_lilith,
    MeeusPolynomialWarning,
    MEEUS_OPTIMAL_CENTURIES,
    MEEUS_VALID_CENTURIES,
    MEEUS_MAX_CENTURIES,
)


class TestMeeusPolynomialValidityConstants:
    """Test validity range constants."""

    def test_optimal_centuries_value(self):
        """Optimal range should be ±200 years (~2 centuries)."""
        assert MEEUS_OPTIMAL_CENTURIES == 2.0

    def test_valid_centuries_value(self):
        """Valid range should be ±1000 years (~10 centuries)."""
        assert MEEUS_VALID_CENTURIES == 10.0

    def test_max_centuries_value(self):
        """Max range should be ±2000 years (~20 centuries)."""
        assert MEEUS_MAX_CENTURIES == 20.0


class TestMeanLunarNodeValidityWarnings:
    """Test warnings for mean lunar node outside validity range."""

    def test_no_warning_within_optimal_range(self):
        """No warning for dates within ±200 years of J2000."""
        # J2000 = 2451545.0, 1 century = 36525 days
        jd_1900 = 2451545.0 - 1.0 * 36525.0  # Year 1900
        jd_2100 = 2451545.0 + 1.0 * 36525.0  # Year 2100

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_1900)
            calc_mean_lunar_node(jd_2100)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0

    def test_no_warning_at_j2000(self):
        """No warning for J2000.0 epoch."""
        jd_j2000 = 2451545.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd_j2000)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            # At J2000, result should be approximately 125° (from formula constant)
            assert 124 < result < 126

    def test_warning_outside_valid_range(self):
        """Warning for dates outside ±1000 years (but within max range)."""
        # Year 500 CE (15 centuries before J2000)
        jd_500 = 2451545.0 - 15.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_500)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the optimal range" in str(meeus_warnings[0].message)

    def test_warning_outside_max_range(self):
        """Warning for dates outside ±2000 years (max range)."""
        # Year -500 BCE (25 centuries before J2000)
        jd_500_bce = 2451545.0 - 25.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_500_bce)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the valid range" in str(meeus_warnings[0].message)
            assert "Error may exceed 1" in str(meeus_warnings[0].message)

    def test_warning_for_distant_future(self):
        """Warning for dates in distant future (year 5000 CE)."""
        # Year 5000 CE (30 centuries after J2000)
        jd_5000 = 2451545.0 + 30.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_5000)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the valid range" in str(meeus_warnings[0].message)

    def test_result_still_valid_with_warning(self):
        """Function should still return a valid result even when warning."""
        # Very ancient date
        jd_ancient = 2451545.0 - 30.0 * 36525.0  # ~1000 BCE

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd_ancient)
            # Result should still be normalized to 0-360
            assert 0 <= result < 360


class TestMeanLilithValidityWarnings:
    """Test warnings for mean Lilith outside validity range."""

    def test_no_warning_within_optimal_range(self):
        """No warning for dates within ±200 years of J2000."""
        jd_1900 = 2451545.0 - 1.0 * 36525.0  # Year 1900
        jd_2100 = 2451545.0 + 1.0 * 36525.0  # Year 2100

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lilith(jd_1900)
            calc_mean_lilith(jd_2100)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0

    def test_no_warning_at_j2000(self):
        """No warning for J2000.0 epoch."""
        jd_j2000 = 2451545.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd_j2000)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            # Result should be valid longitude
            assert 0 <= result < 360

    def test_warning_outside_valid_range(self):
        """Warning for dates outside ±1000 years."""
        # Year 500 CE
        jd_500 = 2451545.0 - 15.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lilith(jd_500)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the optimal range" in str(meeus_warnings[0].message)

    def test_warning_outside_max_range(self):
        """Warning for dates outside ±2000 years."""
        jd_ancient = 2451545.0 - 25.0 * 36525.0  # ~500 BCE

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lilith(jd_ancient)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the valid range" in str(meeus_warnings[0].message)

    def test_result_still_valid_with_warning(self):
        """Function should still return a valid result even when warning."""
        jd_ancient = 2451545.0 - 30.0 * 36525.0

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd_ancient)
            assert 0 <= result < 360


class TestMeeusPolynomialWarningClass:
    """Test the warning class itself."""

    def test_warning_is_user_warning(self):
        """MeeusPolynomialWarning should be a UserWarning subclass."""
        assert issubclass(MeeusPolynomialWarning, UserWarning)

    def test_warning_can_be_filtered(self):
        """Warning can be filtered/suppressed if desired."""
        jd_ancient = 2451545.0 - 25.0 * 36525.0

        # Filter the warning
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("ignore", category=MeeusPolynomialWarning)
            calc_mean_lunar_node(jd_ancient)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0


@pytest.mark.unit
class TestMeeusPolynomialAccuracyDocumentation:
    """Test that documented accuracy claims are reasonable."""

    def test_optimal_range_dates(self):
        """Dates in optimal range (1800-2200) should work without warning."""
        test_dates = [
            2378497.0,  # 1800-01-01
            2451545.0,  # 2000-01-01 (J2000)
            2524594.0,  # 2200-01-01
        ]

        for jd in test_dates:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = calc_mean_lunar_node(jd)
                meeus_warnings = [
                    x for x in w if issubclass(x.category, MeeusPolynomialWarning)
                ]
                assert len(meeus_warnings) == 0
                assert 0 <= result < 360

    def test_valid_range_boundary(self):
        """Test at boundary of valid range (~1000 CE and 3000 CE)."""
        # Year 1000 CE (10 centuries before J2000)
        jd_1000 = 2451545.0 - 10.0 * 36525.0
        # Year 3000 CE (10 centuries after J2000)
        jd_3000 = 2451545.0 + 10.0 * 36525.0

        # These are at the boundary, should not trigger warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_1000)
            calc_mean_lunar_node(jd_3000)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            # At exactly 10 centuries, should not warn (boundary is inclusive)
            assert len(meeus_warnings) == 0

    def test_just_outside_valid_range(self):
        """Test just outside valid range should warn."""
        # Year 999 CE (10.01 centuries before J2000)
        jd_just_outside = 2451545.0 - 10.01 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_just_outside)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
