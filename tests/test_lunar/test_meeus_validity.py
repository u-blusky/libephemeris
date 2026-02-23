"""
Tests for Meeus polynomial validity range warnings in lunar module.

NOTE: With precomputed correction tables from JPL ephemeris, the mean lunar
node and mean Lilith functions no longer emit warnings for dates outside the
polynomial validity range. The tables provide high precision across the full
DE440/DE441 range. These tests now verify that:
1. No warnings are emitted (corrections are applied automatically)
2. Valid results are returned across the full date range
3. The MeeusPolynomialWarning and MeeusRangeError classes still exist for
   backward compatibility and for use by other functions
"""

import pytest
import warnings
from libephemeris.lunar import (
    calc_mean_lunar_node,
    calc_mean_lilith,
    MeeusPolynomialWarning,
    MeeusRangeError,
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
    """Test that mean lunar node works without warnings using correction tables."""

    def test_no_warning_within_optimal_range(self):
        """No warning for dates within ±200 years of J2000."""
        jd_1900 = 2451545.0 - 1.0 * 36525.0
        jd_2100 = 2451545.0 + 1.0 * 36525.0

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
            assert 124 < result < 126

    def test_no_warning_outside_optimal_range(self):
        """No warning for dates outside ±200 years - corrections are applied."""
        jd_1500 = 2451545.0 - 5.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd_1500)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            assert 0 <= result < 360

    def test_warning_outside_valid_range(self):
        """Warning emitted for dates outside ±20 centuries from J2000."""
        jd_0ce = 2451545.0 - 20.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd_0ce)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            # At exactly T=20, boundary is > 20, so no warning at exactly 20
            # but this JD is at T=-20.0 which triggers the > 10 warning
            assert len(meeus_warnings) >= 1
            assert 0 <= result < 360

    def test_no_exception_outside_max_range(self):
        """No exception for dates beyond ±2000 years - corrections cover DE441 range."""
        jd_far = 2451545.0 - 25.0 * 36525.0

        result = calc_mean_lunar_node(jd_far)
        assert 0 <= result < 360

    def test_no_exception_for_distant_future(self):
        """No exception for far future dates - corrections cover DE441 range."""
        jd_future = 2451545.0 + 25.0 * 36525.0

        result = calc_mean_lunar_node(jd_future)
        assert 0 <= result < 360


class TestMeanLilithValidityWarnings:
    """Test that mean Lilith works without warnings using correction tables."""

    def test_no_warning_within_optimal_range(self):
        """No warning for dates within ±200 years of J2000."""
        jd_1900 = 2451545.0 - 1.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lilith(jd_1900)
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
            assert 0 <= result < 360

    def test_no_warning_outside_optimal_range(self):
        """No warning for dates outside ±200 years - corrections are applied."""
        jd_1500 = 2451545.0 - 5.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd_1500)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            assert 0 <= result < 360

    def test_warning_outside_valid_range(self):
        """Warning emitted for dates outside ±20 centuries from J2000."""
        jd_0ce = 2451545.0 - 20.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd_0ce)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            # At T=-20.0 which triggers the > 10 warning
            assert len(meeus_warnings) >= 1
            assert 0 <= result < 360

    def test_no_exception_outside_max_range(self):
        """No exception for dates beyond ±2000 years - corrections cover DE441 range."""
        jd_far = 2451545.0 - 25.0 * 36525.0

        result = calc_mean_lilith(jd_far)
        assert 0 <= result < 360


class TestMeeusPolynomialWarningClass:
    """Test the warning class itself."""

    def test_warning_is_user_warning(self):
        """MeeusPolynomialWarning should be a UserWarning subclass."""
        assert issubclass(MeeusPolynomialWarning, UserWarning)


class TestMeeusRangeErrorClass:
    """Test that MeeusRangeError class exists for backward compatibility."""

    def test_error_is_value_error_subclass(self):
        """MeeusRangeError should be a ValueError subclass."""
        assert issubclass(MeeusRangeError, ValueError)

    def test_error_can_be_instantiated(self):
        """MeeusRangeError can be instantiated."""
        err = MeeusRangeError("test message")
        assert str(err) == "test message"


@pytest.mark.unit
class TestMeeusPolynomialAccuracyDocumentation:
    """Test that mean node/lilith produce valid results across ranges."""

    def test_result_in_valid_range(self):
        """Results should always be in [0, 360) range."""
        for centuries in [-15, -10, -5, 0, 5, 10, 15]:
            jd = 2451545.0 + centuries * 36525.0
            result = calc_mean_lunar_node(jd)
            assert 0 <= result < 360

    def test_optimal_range_dates(self):
        """Dates in optimal range (1800-2200) should work without warning."""
        test_dates = [
            2378497.0,
            2451545.0,
            2524594.0,
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
        jd_1000 = 2451545.0 - 10.0 * 36525.0
        jd_3000 = 2451545.0 + 10.0 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd_1000)
            calc_mean_lunar_node(jd_3000)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0

    def test_just_outside_valid_range(self):
        """Warning emitted just outside valid range (>10 centuries)."""
        jd_just_outside = 2451545.0 - 10.01 * 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd_just_outside)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            # T=-10.01 is outside optimal range, warning expected
            assert len(meeus_warnings) >= 1
            assert 0 <= result < 360
