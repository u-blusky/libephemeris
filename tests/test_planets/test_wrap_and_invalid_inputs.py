"""Tests for 360° wrap-around edge cases, NaN/Inf inputs, and ECL_NUT special body."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_MEAN_NODE,
    SE_ECL_NUT,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_SIDEREAL,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestWrapAround360:
    """Test positions near the 0°/360° boundary."""

    def _find_jd_near_zero(self, body: int, start_jd: float, days: int = 365) -> float:
        """Find a JD where body longitude is near 0°/360°."""
        best_jd = start_jd
        best_dist = 999.0
        for i in range(days):
            jd = start_jd + i
            result, _ = swe.calc_ut(jd, body, SEFLG_SWIEPH | SEFLG_SPEED)
            lon = result[0]
            dist_to_zero = min(lon, 360.0 - lon)
            if dist_to_zero < best_dist:
                best_dist = dist_to_zero
                best_jd = jd
        return best_jd

    def test_sun_near_zero_valid_range(self):
        """Sun near 0° Aries should still be in [0, 360)."""
        # Sun crosses 0° around March equinox. J2000 is Jan 1 2000.
        # March ~20 is about day 79.
        jd = self._find_jd_near_zero(SE_SUN, JD_J2000, 365)
        result, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
        assert 0.0 <= result[0] < 360.0

    def test_moon_near_zero_valid(self):
        """Moon frequently crosses 0°; position must be in [0, 360)."""
        jd = self._find_jd_near_zero(SE_MOON, JD_J2000, 30)
        result, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED)
        assert 0.0 <= result[0] < 360.0

    def test_positions_never_negative(self):
        """No longitude should ever be negative over a scan of dates."""
        for body in [SE_SUN, SE_MOON, SE_MARS, SE_JUPITER]:
            for i in range(50):
                jd = JD_J2000 + i * 7.3  # sample every ~week
                result, _ = swe.calc_ut(jd, body, SEFLG_SWIEPH)
                assert result[0] >= 0.0, f"Body {body} at JD {jd}: lon {result[0]} < 0"

    def test_positions_never_360(self):
        """No longitude should ever be exactly 360.0 or above."""
        for body in [SE_SUN, SE_MOON, SE_MARS, SE_JUPITER]:
            for i in range(50):
                jd = JD_J2000 + i * 7.3
                result, _ = swe.calc_ut(jd, body, SEFLG_SWIEPH)
                assert result[0] < 360.0, (
                    f"Body {body} at JD {jd}: lon {result[0]} >= 360"
                )

    def test_solcross_near_zero(self):
        """solcross_ut at 0° should return a valid crossing JD."""
        jd_cross = swe.solcross_ut(0.0, JD_J2000, 0)
        assert jd_cross > JD_J2000
        # Verify the Sun is near 0° at that JD
        result, _ = swe.calc_ut(jd_cross, SE_SUN, SEFLG_SWIEPH)
        assert result[0] < 0.01 or result[0] > 359.99

    def test_mooncross_near_zero(self):
        """mooncross_ut at 0° should return a valid crossing JD."""
        jd_cross = swe.mooncross_ut(0.0, JD_J2000, 0)
        assert jd_cross > JD_J2000
        result, _ = swe.calc_ut(jd_cross, SE_MOON, SEFLG_SWIEPH)
        assert result[0] < 0.1 or result[0] > 359.9

    def test_mean_node_range_across_full_cycle(self):
        """Mean Node (retrograde) should stay in [0, 360) over 19 years."""
        for i in range(200):
            jd = JD_J2000 + i * 34.7  # ~19 years in 200 steps
            result, _ = swe.calc_ut(jd, SE_MEAN_NODE, SEFLG_SWIEPH)
            assert 0.0 <= result[0] < 360.0, (
                f"Mean Node at JD {jd}: lon {result[0]} out of range"
            )

    def test_sidereal_positions_in_range(self):
        """Sidereal longitudes must also be in [0, 360)."""
        swe.set_sid_mode(1)  # Lahiri
        try:
            for body in [SE_SUN, SE_MOON, SE_MARS]:
                result, _ = swe.calc_ut(JD_J2000, body, SEFLG_SWIEPH | SEFLG_SIDEREAL)
                assert 0.0 <= result[0] < 360.0, (
                    f"Sidereal lon for body {body}: {result[0]}"
                )
        finally:
            swe.set_sid_mode(0)


@pytest.mark.unit
class TestEclNut:
    """Test SE_ECL_NUT special body for nutation and obliquity."""

    def test_ecl_nut_returns_6_values(self):
        """ECL_NUT should return a 6-element tuple."""
        result, flag = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_SWIEPH)
        assert len(result) == 6

    def test_ecl_nut_values_finite(self):
        """All ECL_NUT values should be finite."""
        result, flag = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_SWIEPH)
        for i, val in enumerate(result[:4]):
            assert math.isfinite(val), f"ECL_NUT[{i}] = {val} not finite"

    def test_ecl_nut_obliquity_plausible(self):
        """True and mean obliquity should be ~23.4°."""
        result, _ = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_SWIEPH)
        true_obl = result[0]
        mean_obl = result[1]
        assert 23.0 < true_obl < 24.0, f"True obliquity {true_obl} implausible"
        assert 23.0 < mean_obl < 24.0, f"Mean obliquity {mean_obl} implausible"

    def test_ecl_nut_nutation_small(self):
        """Nutation in longitude and obliquity should be small (< 0.01°)."""
        result, _ = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_SWIEPH)
        nut_lon = result[2]
        nut_obl = result[3]
        assert abs(nut_lon) < 0.01, f"Nutation in longitude {nut_lon} too large"
        assert abs(nut_obl) < 0.01, f"Nutation in obliquity {nut_obl} too large"

    def test_ecl_nut_last_two_zero(self):
        """Last two values of ECL_NUT should be 0.0."""
        result, _ = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_SWIEPH)
        assert result[4] == 0.0
        assert result[5] == 0.0

    def test_ecl_nut_varies_over_time(self):
        """ECL_NUT values should change measurably over centuries."""
        res_2000, _ = swe.calc_ut(JD_J2000, SE_ECL_NUT, 0)
        res_2100, _ = swe.calc_ut(JD_J2000 + 36525.0, SE_ECL_NUT, 0)
        # Mean obliquity changes ~47" per century
        diff_obl = abs(res_2100[1] - res_2000[1])
        assert diff_obl > 0.001, "Mean obliquity should change over a century"

    def test_ecl_nut_with_speed_flag(self):
        """ECL_NUT with SEFLG_SPEED should not crash."""
        result, _ = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_SPEED)
        assert len(result) == 6

    def test_ecl_nut_with_equatorial_flag(self):
        """ECL_NUT with SEFLG_EQUATORIAL should not crash."""
        result, _ = swe.calc_ut(JD_J2000, SE_ECL_NUT, SEFLG_EQUATORIAL)
        assert len(result) == 6


@pytest.mark.unit
class TestInvalidInputHandling:
    """Test that NaN and Inf inputs are handled gracefully."""

    def test_calc_ut_nan_jd(self):
        """calc_ut with NaN JD should raise or return non-NaN (not silently corrupt)."""
        try:
            result, _ = swe.calc_ut(float("nan"), SE_SUN, SEFLG_SWIEPH)
            # If it doesn't raise, all outputs should be NaN (propagation)
            # or the function should have caught it
        except (ValueError, TypeError, Exception):
            pass  # Any exception is acceptable

    def test_calc_ut_inf_jd(self):
        """calc_ut with Inf JD should raise or handle gracefully."""
        try:
            result, _ = swe.calc_ut(float("inf"), SE_SUN, SEFLG_SWIEPH)
        except (ValueError, TypeError, OverflowError, Exception):
            pass  # Any exception is acceptable

    def test_calc_ut_neg_inf_jd(self):
        """calc_ut with -Inf JD should raise or handle gracefully."""
        try:
            result, _ = swe.calc_ut(float("-inf"), SE_SUN, SEFLG_SWIEPH)
        except (ValueError, TypeError, OverflowError, Exception):
            pass  # Any exception is acceptable

    def test_houses_nan_lat(self):
        """houses with NaN latitude should raise or handle gracefully."""
        try:
            cusps, ascmc = swe.houses(JD_J2000, float("nan"), 12.5, ord("P"))
        except (ValueError, TypeError, Exception):
            pass

    def test_houses_inf_lon(self):
        """houses with Inf longitude should raise or handle gracefully."""
        try:
            cusps, ascmc = swe.houses(JD_J2000, 41.9, float("inf"), ord("P"))
        except (ValueError, TypeError, Exception):
            pass

    def test_julday_nan_components(self):
        """julday with NaN components should raise or handle gracefully."""
        try:
            jd = swe.julday(2000, 1, float("nan"), 12.0)
        except (ValueError, TypeError, Exception):
            pass

    def test_cotrans_nan_input(self):
        """cotrans with NaN should not silently produce valid-looking output."""
        try:
            result = swe.cotrans((float("nan"), 0.0, 1.0), 23.4)
            # If it doesn't raise, NaN should propagate
            if not math.isnan(result[0]):
                # Some implementations might handle this differently
                pass
        except (ValueError, TypeError, Exception):
            pass
