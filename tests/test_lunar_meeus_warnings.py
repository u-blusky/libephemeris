"""
Tests for Meeus polynomial validity range warnings and exceptions in lunar.py.

These tests verify that:
1. Warnings are issued when dates are outside the optimal range (±200 years)
2. Warnings are issued when dates are outside the valid range (±1000 years)
3. Exceptions are raised when dates are beyond the maximum range (±2000 years)
"""

import warnings
import pytest
from libephemeris.lunar import (
    calc_mean_lunar_node,
    calc_mean_lilith,
    MeeusPolynomialWarning,
    MeeusRangeError,
    MEEUS_OPTIMAL_CENTURIES,
    MEEUS_VALID_CENTURIES,
    MEEUS_MAX_CENTURIES,
)


# Julian Date constants for test dates
# J2000.0 = 2451545.0 (Jan 1, 2000 12:00 TT)
JD_J2000 = 2451545.0
# One Julian century = 36525 days
CENTURY_DAYS = 36525.0


def jd_for_year(year: int) -> float:
    """Approximate Julian Date for January 1 of a given year."""
    # T = (year - 2000) / 100 centuries
    # JD = J2000 + T * 36525
    return JD_J2000 + (year - 2000) / 100.0 * CENTURY_DAYS


class TestMeeusConstants:
    """Test that Meeus range constants have expected values."""

    def test_optimal_centuries_is_two(self):
        """Optimal range is ±200 years (±2 centuries)."""
        assert MEEUS_OPTIMAL_CENTURIES == 2.0

    def test_valid_centuries_is_ten(self):
        """Valid range is ±1000 years (±10 centuries)."""
        assert MEEUS_VALID_CENTURIES == 10.0

    def test_max_centuries_is_twenty(self):
        """Maximum range is ±2000 years (±20 centuries)."""
        assert MEEUS_MAX_CENTURIES == 20.0


class TestCalcMeanLunarNodeWarnings:
    """Test warnings for calc_mean_lunar_node at various date ranges."""

    def test_no_warning_within_optimal_range(self):
        """No warning for dates within ±200 years of J2000."""
        # Test year 1900 (100 years before J2000)
        jd = jd_for_year(1900)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)
            # Filter for MeeusPolynomialWarning only
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            assert 0 <= result < 360

    def test_no_warning_at_optimal_boundary(self):
        """No warning for dates exactly at ±200 years."""
        # Test year 1800 (exactly 200 years before J2000)
        jd = jd_for_year(1800)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            assert 0 <= result < 360

    def test_warning_outside_optimal_range(self):
        """Warning for dates between ±200 and ±1000 years."""
        # Test year 1500 (500 years before J2000)
        jd = jd_for_year(1500)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the optimal range" in str(meeus_warnings[0].message)
            assert "1800-2200" in str(meeus_warnings[0].message)
            assert 0 <= result < 360

    def test_warning_outside_valid_range(self):
        """Warning for dates between ±1000 and ±2000 years."""
        # Test year 500 (1500 years before J2000)
        jd = jd_for_year(500)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the recommended range" in str(meeus_warnings[0].message)
            assert "1000-3000" in str(meeus_warnings[0].message)
            assert 0 <= result < 360

    def test_exception_beyond_max_range(self):
        """Exception for dates beyond ±2000 years."""
        # Test year -500 (2500 years before J2000)
        jd = jd_for_year(-500)
        with pytest.raises(MeeusRangeError) as exc_info:
            calc_mean_lunar_node(jd)
        assert "outside the valid range" in str(exc_info.value)
        assert "0-4000 CE" in str(exc_info.value)

    def test_exception_far_future(self):
        """Exception for dates far in the future."""
        # Test year 5000 (3000 years after J2000)
        jd = jd_for_year(5000)
        with pytest.raises(MeeusRangeError) as exc_info:
            calc_mean_lunar_node(jd)
        assert "outside the valid range" in str(exc_info.value)


class TestCalcMeanLilithWarnings:
    """Test warnings for calc_mean_lilith at various date ranges."""

    def test_no_warning_within_optimal_range(self):
        """No warning for dates within ±200 years of J2000."""
        # Test year 2100 (100 years after J2000)
        jd = jd_for_year(2100)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 0
            assert 0 <= result < 360

    def test_warning_outside_optimal_range(self):
        """Warning for dates between ±200 and ±1000 years."""
        # Test year 2500 (500 years after J2000)
        jd = jd_for_year(2500)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the optimal range" in str(meeus_warnings[0].message)
            assert "1800-2200" in str(meeus_warnings[0].message)
            assert 0 <= result < 360

    def test_warning_outside_valid_range(self):
        """Warning for dates between ±1000 and ±2000 years."""
        # Test year 3500 (1500 years after J2000)
        jd = jd_for_year(3500)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lilith(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            assert "outside the recommended range" in str(meeus_warnings[0].message)
            assert "1000-3000" in str(meeus_warnings[0].message)
            assert 0 <= result < 360

    def test_exception_beyond_max_range(self):
        """Exception for dates beyond ±2000 years."""
        # Test year -1000 (3000 years before J2000)
        jd = jd_for_year(-1000)
        with pytest.raises(MeeusRangeError) as exc_info:
            calc_mean_lilith(jd)
        assert "outside the valid range" in str(exc_info.value)
        assert "0-4000 CE" in str(exc_info.value)

    def test_exception_far_future(self):
        """Exception for dates far in the future."""
        # Test year 6000 (4000 years after J2000)
        jd = jd_for_year(6000)
        with pytest.raises(MeeusRangeError) as exc_info:
            calc_mean_lilith(jd)
        assert "outside the valid range" in str(exc_info.value)


class TestMeeusWarningMessages:
    """Test that warning messages contain useful information."""

    def test_optimal_warning_contains_precision_info(self):
        """Optimal range warning mentions precision degradation."""
        jd = jd_for_year(1600)  # 400 years before J2000
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            # Should mention precision degradation
            assert "Precision degrades" in str(meeus_warnings[0].message)

    def test_valid_warning_contains_error_estimate(self):
        """Valid range warning mentions error estimate."""
        jd = jd_for_year(200)  # 1800 years before J2000
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd)
            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]
            assert len(meeus_warnings) == 1
            # Should mention error estimate
            assert "0.1-1" in str(meeus_warnings[0].message)

    def test_exception_suggests_alternatives(self):
        """Max range exception suggests alternative approaches."""
        jd = jd_for_year(-500)
        with pytest.raises(MeeusRangeError) as exc_info:
            calc_mean_lunar_node(jd)
        # Should suggest numerical integration or alternatives
        assert "numerical integration" in str(exc_info.value).lower()


class TestMeeusRangeErrorIsValueError:
    """Test that MeeusRangeError can be caught as ValueError."""

    def test_meeus_range_error_is_value_error(self):
        """MeeusRangeError inherits from ValueError."""
        jd = jd_for_year(-500)
        with pytest.raises(ValueError):
            calc_mean_lunar_node(jd)

    def test_meeus_range_error_type(self):
        """MeeusRangeError is the specific type raised."""
        assert issubclass(MeeusRangeError, ValueError)


class TestMeeusPolynomialWarningType:
    """Test that MeeusPolynomialWarning is a proper warning type."""

    def test_is_user_warning(self):
        """MeeusPolynomialWarning inherits from UserWarning."""
        assert issubclass(MeeusPolynomialWarning, UserWarning)

    def test_can_be_filtered(self):
        """MeeusPolynomialWarning can be filtered with warnings module."""
        jd = jd_for_year(1500)
        # Test that we can ignore the warning
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=MeeusPolynomialWarning)
            # Should not raise
            result = calc_mean_lunar_node(jd)
            assert 0 <= result < 360
