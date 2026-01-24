"""
Comprehensive tests for longitude crossing functions.

Tests solcross_ut, mooncross_ut, and cross_ut.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestSolcrossBasic:
    """Basic tests for Sun crossing function."""

    @pytest.mark.unit
    def test_solcross_vernal_equinox(self):
        """Find when Sun crosses 0° (vernal equinox)."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        jd_cross = ephem.swe_solcross_ut(0.0, jd_start, 0)

        # Should be around March 20, 2024
        year, month, day, hour = ephem.swe_revjul(jd_cross)
        assert year == 2024
        assert month == 3
        assert 19 <= day <= 21

    @pytest.mark.unit
    def test_solcross_summer_solstice(self):
        """Find when Sun crosses 90° (summer solstice)."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        jd_cross = ephem.swe_solcross_ut(90.0, jd_start, 0)

        year, month, day, hour = ephem.swe_revjul(jd_cross)
        assert year == 2024
        assert month == 6
        assert 20 <= day <= 22

    @pytest.mark.unit
    def test_solcross_precision(self):
        """Sun should be very close to target at crossing time."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 45.0
        jd_cross = ephem.swe_solcross_ut(target, jd_start, 0)

        # Check Sun position at crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sun at {pos[0]}, target {target}, diff {diff}"


class TestSolcrossVsPyswisseph:
    """Compare solcross with pyswisseph."""

    @pytest.mark.comparison
    def test_solcross_equinox_timing(self):
        """Crossing time should match pyswisseph."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_lib = ephem.swe_solcross_ut(0.0, jd_start, 0)
        jd_swe = swe.solcross_ut(0.0, jd_start, 0)

        # Difference should be less than 1 minute
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 60, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_all_signs(self, target):
        """All zodiac sign ingresses should match."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_lib = ephem.swe_solcross_ut(float(target), jd_start, 0)
        jd_swe = swe.solcross_ut(float(target), jd_start, 0)

        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 120, f"Target {target}° diff {diff_seconds} seconds"


class TestMooncrossBasic:
    """Basic tests for Moon crossing function."""

    @pytest.mark.unit
    def test_mooncross_basic(self):
        """Moon crossing should return valid JD."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 90.0
        jd_cross = ephem.swe_mooncross_ut(target, jd_start, 0)

        assert jd_cross > jd_start
        # Should be within one lunar orbit (~27 days)
        assert jd_cross < jd_start + 28

    @pytest.mark.unit
    def test_mooncross_precision(self):
        """Moon should be close to target at crossing."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 123.456
        jd_cross = ephem.swe_mooncross_ut(target, jd_start, 0)

        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Moon at {pos[0]}, target {target}, diff {diff}"


class TestMooncrossVsPyswisseph:
    """Compare mooncross with pyswisseph."""

    @pytest.mark.comparison
    def test_mooncross_timing(self):
        """Moon crossing should match pyswisseph."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 180.0

        jd_lib = ephem.swe_mooncross_ut(target, jd_start, 0)
        jd_swe = swe.mooncross_ut(target, jd_start, 0)

        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon crossing diff {diff_seconds} seconds"


class TestCrossingConsecutive:
    """Test finding consecutive crossings."""

    @pytest.mark.unit
    def test_consecutive_sun_crossings(self):
        """Should find 12 consecutive Sun crossings in a year."""
        jd = ephem.swe_julday(2024, 1, 1, 0.0)
        crossings = []

        for target in range(0, 360, 30):
            jd_cross = ephem.swe_solcross_ut(float(target), jd, 0)
            crossings.append(jd_cross)

        # All should be in chronological order (after accounting for year wrap)
        # First 0° crossing should be in March
        assert crossings[0] > ephem.swe_julday(2024, 3, 1, 0.0)

    @pytest.mark.unit
    def test_12_moon_crossings_in_month(self):
        """Should find ~12 crossings in 28 days."""
        jd = ephem.swe_julday(2024, 1, 1, 0.0)
        crossings = []

        for target in range(0, 360, 30):
            jd_cross = ephem.swe_mooncross_ut(float(target), jd, 0)
            if jd_cross < jd + 28:
                crossings.append(jd_cross)

        # Should find all 12 in about 27 days
        assert len(crossings) == 12


class TestCrossingEdgeCases:
    """Test edge cases for crossing functions."""

    @pytest.mark.edge_case
    def test_target_360_equals_0(self):
        """360° should be same as 0°."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_360 = ephem.swe_solcross_ut(360.0, jd_start, 0)
        jd_0 = ephem.swe_solcross_ut(0.0, jd_start, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001

    @pytest.mark.edge_case
    def test_target_negative_normalized(self):
        """Negative target should be normalized."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_neg = ephem.swe_solcross_ut(-30.0, jd_start, 0)  # = 330°
        jd_pos = ephem.swe_solcross_ut(330.0, jd_start, 0)

        assert abs(jd_neg - jd_pos) < 0.001


class TestCrossUtGeneric:
    """Test generic planet crossing function."""

    @pytest.mark.unit
    def test_cross_ut_mercury(self):
        """Should work for Mercury."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 0.0

        jd_cross = ephem.swe_cross_ut(SE_MERCURY, target, jd_start, 0)

        assert jd_cross > jd_start

        # Check Mercury position at crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MERCURY, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1  # Mercury moves fast, allow more tolerance
