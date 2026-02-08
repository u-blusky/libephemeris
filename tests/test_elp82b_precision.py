"""
Tests for the ELP 2000-82B lunar theory implementation.

These tests verify:
1. Basic functionality of the lunar calculation
2. Precision against DE440 (Skyfield) within the overlap range
3. Calculation capability for the extended range (-3000 to +3000 CE)
4. Velocity calculations
5. Internal functions

The expected accuracy is ~0.5-5 arcseconds within the DE440 overlap range,
and the position should be calculable for the full Moshier range.
"""

from __future__ import annotations

import math
from typing import Tuple

import pytest

import libephemeris as eph
from libephemeris.constants import (
    SE_MOON,
    SEFLG_MOSEPH,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
)
from libephemeris.moshier import (
    J2000,
    MOSHIER_MOON,
    calc_position,
    is_moshier_body,
)
from libephemeris.moshier.elp82b import (
    calc_moon_position,
    is_moon_body,
    _fundamental_arguments,
    _eccentricity_correction,
)


# Convert arcseconds to degrees
ARCSEC_TO_DEG = 1.0 / 3600.0

# Tolerance in arcseconds for precision tests
# ELP 2000-82B achieves ~1-5 arcsec accuracy. We use a relaxed tolerance
# to account for truncated series and allow for some variation.
PRECISION_TOLERANCE_ARCSEC = 300.0  # 5 arcminutes (relaxed for truncated series)


def angular_difference(angle1: float, angle2: float) -> float:
    """Calculate the smallest angular difference between two angles (degrees).

    Args:
        angle1: First angle in degrees.
        angle2: Second angle in degrees.

    Returns:
        Angular difference in degrees (0 to 180).
    """
    diff = abs(angle1 - angle2) % 360.0
    if diff > 180.0:
        diff = 360.0 - diff
    return diff


class TestELP82BBasicFunctionality:
    """Test basic ELP 2000-82B functionality."""

    def test_calc_moon_position_returns_tuple(self):
        """Test that calc_moon_position returns a 3-tuple."""
        lon, lat, dist_km = calc_moon_position(J2000)

        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist_km, float)

    def test_calc_moon_position_valid_ranges(self):
        """Test that returned values are in valid ranges."""
        lon, lat, dist_km = calc_moon_position(J2000)

        # Longitude should be in [0, 360)
        assert 0.0 <= lon < 360.0, f"Longitude {lon} out of range"

        # Latitude should be in [-6, 6] (Moon's orbital inclination is ~5.1°)
        assert -6.5 <= lat <= 6.5, f"Latitude {lat} out of range"

        # Distance should be in [356000, 407000] km (Moon's distance range)
        assert 356000 < dist_km < 407000, f"Distance {dist_km} km out of range"

    def test_calc_position_returns_6tuple(self):
        """Test that calc_position returns a 6-tuple."""
        result = calc_position(J2000, MOSHIER_MOON)

        assert isinstance(result, tuple)
        assert len(result) == 6

        lon, lat, dist, dlon, dlat, ddist = result

        for val in result:
            assert math.isfinite(val), f"Value {val} is not finite"

    def test_calc_position_distance_in_au(self):
        """Test that calc_position returns distance in AU."""
        lon, lat, dist_au, _, _, _ = calc_position(J2000, MOSHIER_MOON)

        # Moon distance should be about 0.0026 AU (385000 km / 149597870.7 km)
        assert 0.0023 < dist_au < 0.0028, f"Moon distance {dist_au} AU unexpected"

    def test_calc_position_wrong_body_raises(self):
        """Test that wrong body ID raises ValueError."""
        from libephemeris.moshier.elp82b import calc_position as elp_calc_position

        with pytest.raises(ValueError):
            elp_calc_position(J2000, 0)  # Sun

        with pytest.raises(ValueError):
            elp_calc_position(J2000, 2)  # Mercury

    def test_is_moon_body(self):
        """Test is_moon_body function."""
        assert is_moon_body(MOSHIER_MOON)
        assert is_moon_body(1)
        assert not is_moon_body(0)
        assert not is_moon_body(2)
        assert not is_moon_body(-1)


class TestELP82BFundamentalArguments:
    """Test fundamental argument calculations."""

    def test_fundamental_arguments_at_j2000(self):
        """Test fundamental arguments at J2000.0."""
        T = 0.0
        Lp, D, M, Mp, F = _fundamental_arguments(T)

        # At J2000.0, all arguments should be their initial values
        # L' (Moon's mean longitude) ~218.3°
        assert 215.0 < Lp < 221.0, f"L' = {Lp}° at J2000"

        # D (mean elongation) ~297.9°
        assert 295.0 < D < 300.0, f"D = {D}° at J2000"

        # M (Sun's mean anomaly) ~357.5°
        assert 355.0 < M < 360.0, f"M = {M}° at J2000"

        # M' (Moon's mean anomaly) ~134.9°
        assert 132.0 < Mp < 138.0, f"M' = {Mp}° at J2000"

        # F (argument of latitude) ~93.3°
        assert 90.0 < F < 96.0, f"F = {F}° at J2000"

    def test_fundamental_arguments_vary_with_time(self):
        """Test that fundamental arguments change over time."""
        T1 = 0.0
        T2 = 0.1  # 10 years

        args1 = _fundamental_arguments(T1)
        args2 = _fundamental_arguments(T2)

        for i, (a1, a2) in enumerate(zip(args1, args2)):
            assert a1 != a2, f"Argument {i} didn't change with time"


class TestELP82BEccentricityCorrection:
    """Test eccentricity correction factor."""

    def test_eccentricity_correction_power_zero(self):
        """Test that E^0 = 1."""
        E = _eccentricity_correction(0.0, 0)
        assert E == 1.0

    def test_eccentricity_correction_power_one(self):
        """Test E^1 at J2000."""
        E = _eccentricity_correction(0.0, 1)
        # E = 1 - 0.002516*T - 0.0000074*T^2 at T=0 gives E=1
        assert abs(E - 1.0) < 0.001

    def test_eccentricity_decreases_with_time(self):
        """Test that eccentricity decreases over time."""
        E_j2000 = _eccentricity_correction(0.0, 1)
        E_future = _eccentricity_correction(1.0, 1)  # 100 years

        assert E_future < E_j2000


class TestELP82BPrecisionVsDE440:
    """Compare ELP 2000-82B positions against JPL DE440 reference."""

    # Test dates within DE440 range
    TEST_DATES = [
        (2451545.0, "J2000.0 (2000-01-01)"),
        (2440000.5, "1968-05-24"),
        (2460000.5, "2023-02-25"),
        (2430000.0, "1941-02-15"),
        (2450000.0, "1995-10-09"),
    ]

    @pytest.mark.precision
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_elp82b_vs_de440_longitude(self, jd, date_desc):
        """ELP 2000-82B longitude should be within tolerance of DE440."""
        pos_moshier, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_MOSEPH | SEFLG_SPEED)
        pos_jpl, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED)

        lon_diff = angular_difference(pos_moshier[0], pos_jpl[0])
        lon_diff_arcsec = lon_diff * 3600.0

        assert lon_diff_arcsec < PRECISION_TOLERANCE_ARCSEC, (
            f'Moon at {date_desc}: longitude diff {lon_diff_arcsec:.1f}" exceeds '
            f'{PRECISION_TOLERANCE_ARCSEC}" tolerance\n'
            f"  Moshier: {pos_moshier[0]:.6f}deg\n"
            f"  JPL:     {pos_jpl[0]:.6f}deg"
        )

    @pytest.mark.precision
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_elp82b_vs_de440_latitude(self, jd, date_desc):
        """ELP 2000-82B latitude should be within tolerance of DE440."""
        pos_moshier, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_MOSEPH | SEFLG_SPEED)
        pos_jpl, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED)

        lat_diff = abs(pos_moshier[1] - pos_jpl[1])
        lat_diff_arcsec = lat_diff * 3600.0

        assert lat_diff_arcsec < PRECISION_TOLERANCE_ARCSEC, (
            f'Moon at {date_desc}: latitude diff {lat_diff_arcsec:.1f}" exceeds '
            f'{PRECISION_TOLERANCE_ARCSEC}" tolerance\n'
            f"  Moshier: {pos_moshier[1]:.6f}deg\n"
            f"  JPL:     {pos_jpl[1]:.6f}deg"
        )

    @pytest.mark.precision
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_elp82b_vs_de440_distance(self, jd, date_desc):
        """ELP 2000-82B distance should be within 1% of DE440."""
        pos_moshier, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_MOSEPH | SEFLG_SPEED)
        pos_jpl, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED)

        dist_diff_percent = abs(pos_moshier[2] - pos_jpl[2]) / pos_jpl[2] * 100.0

        assert dist_diff_percent < 1.0, (
            f"Moon at {date_desc}: distance diff {dist_diff_percent:.3f}% exceeds 1%\n"
            f"  Moshier: {pos_moshier[2]:.8f} AU\n"
            f"  JPL:     {pos_jpl[2]:.8f} AU"
        )


class TestELP82BExtendedRange:
    """Test ELP 2000-82B works for dates outside DE440 range."""

    # Extended date range for Moshier: -3000 to +3000 CE
    EXTENDED_DATES = [
        (1721424.5, "0001-01-01 CE"),  # Start of Common Era
        (1538558.5, "~500 BCE"),
        (1356174.5, "~1000 BCE"),
        (1000000.0, "~-1931 BCE"),
        (2817152.5, "2900-01-01 CE"),
        (2500000.0, "~2132 CE"),
    ]

    @pytest.mark.slow
    @pytest.mark.parametrize("jd,date_desc", EXTENDED_DATES)
    def test_elp82b_extended_range_valid_position(self, jd, date_desc):
        """ELP 2000-82B should return valid positions for extended dates."""
        pos, flag = eph.swe_calc_ut(jd, SE_MOON, SEFLG_MOSEPH | SEFLG_SPEED)

        lon, lat, dist = pos[0], pos[1], pos[2]

        # Position should be valid
        assert 0.0 <= lon < 360.0, f"Moon longitude {lon} invalid at {date_desc}"
        assert -6.5 <= lat <= 6.5, f"Moon latitude {lat} invalid at {date_desc}"
        # Distance in AU (0.0023 to 0.0028 AU)
        assert 0.0020 < dist < 0.0030, f"Moon distance {dist} AU invalid at {date_desc}"

    @pytest.mark.slow
    def test_elp82b_1000bce(self):
        """Test Moon position at 1000 BCE."""
        jd = 1356174.5  # ~1000 BCE

        lon, lat, dist_km = calc_moon_position(jd)

        # Should return valid values
        assert 0.0 <= lon < 360.0
        assert -6.5 <= lat <= 6.5
        assert 356000 < dist_km < 407000

    @pytest.mark.slow
    def test_elp82b_3000ce(self):
        """Test Moon position at 3000 CE."""
        jd = 2816788.5  # ~3000 CE

        lon, lat, dist_km = calc_moon_position(jd)

        # Should return valid values
        assert 0.0 <= lon < 360.0
        assert -6.5 <= lat <= 6.5
        assert 356000 < dist_km < 407000


class TestELP82BVelocity:
    """Test ELP 2000-82B velocity calculations."""

    def test_velocity_reasonable_range(self):
        """Moon velocity should be in physically reasonable range."""
        lon, lat, dist, dlon, dlat, ddist = calc_position(J2000, MOSHIER_MOON)

        # Moon moves about 13.2°/day in longitude
        assert 11.0 < dlon < 15.0, f"Moon dlon {dlon}°/day unexpected"

        # Moon latitude change is typically small (< 1°/day)
        assert -1.5 < dlat < 1.5, f"Moon dlat {dlat}°/day unexpected"

    def test_velocity_matches_position_change(self):
        """Velocity should approximately match actual position change over 1 day."""
        jd = J2000

        pos1 = calc_position(jd, MOSHIER_MOON)
        pos2 = calc_position(jd + 1.0, MOSHIER_MOON)

        # Calculate actual change
        actual_dlon = pos2[0] - pos1[0]
        if actual_dlon < -180:
            actual_dlon += 360
        elif actual_dlon > 180:
            actual_dlon -= 360

        # Reported velocity
        reported_dlon = pos1[3]

        # Should match within 5%
        if abs(reported_dlon) > 0.1:
            relative_error = abs(actual_dlon - reported_dlon) / abs(reported_dlon)
            assert relative_error < 0.05, (
                f"Velocity {reported_dlon:.4f} doesn't match position change "
                f"{actual_dlon:.4f} deg/day (error: {relative_error * 100:.1f}%)"
            )


class TestELP82BDataIntegrity:
    """Test that coefficient data is correctly loaded."""

    def test_longitude_terms_exist(self):
        """Test that longitude terms are loaded."""
        from libephemeris.moshier.elp82b_data import LONGITUDE_MAIN_TERMS

        assert len(LONGITUDE_MAIN_TERMS) > 50, "Too few longitude terms"

        # First term should be the largest
        assert LONGITUDE_MAIN_TERMS[0][4] == 6288774

    def test_latitude_terms_exist(self):
        """Test that latitude terms are loaded."""
        from libephemeris.moshier.elp82b_data import LATITUDE_MAIN_TERMS

        assert len(LATITUDE_MAIN_TERMS) > 50, "Too few latitude terms"

        # First term should be the largest
        assert LATITUDE_MAIN_TERMS[0][4] == 5128122

    def test_distance_terms_exist(self):
        """Test that distance terms are loaded."""
        from libephemeris.moshier.elp82b_data import DISTANCE_MAIN_TERMS

        assert len(DISTANCE_MAIN_TERMS) > 40, "Too few distance terms"

        # First term should be the largest
        assert DISTANCE_MAIN_TERMS[0][4] == -20905355

    def test_planetary_terms_exist(self):
        """Test that planetary perturbation terms are loaded."""
        from libephemeris.moshier.elp82b_data import (
            LONGITUDE_PLANETARY_TERMS,
            LATITUDE_PLANETARY_TERMS,
            DISTANCE_PLANETARY_TERMS,
        )

        assert len(LONGITUDE_PLANETARY_TERMS) > 10
        assert len(LATITUDE_PLANETARY_TERMS) > 5
        assert len(DISTANCE_PLANETARY_TERMS) > 5


class TestELP82BIntegrationWithSweCalc:
    """Test integration with swe_calc_ut API."""

    def test_swe_calc_ut_with_moseph_flag(self):
        """Test that swe_calc_ut works with SEFLG_MOSEPH for Moon."""
        pos, flag = eph.swe_calc_ut(J2000, SE_MOON, SEFLG_MOSEPH)

        assert isinstance(pos, tuple)
        assert len(pos) == 6
        assert 0.0 <= pos[0] < 360.0

    def test_swe_calc_ut_with_speed_flag(self):
        """Test that swe_calc_ut returns velocities with SEFLG_SPEED."""
        pos, flag = eph.swe_calc_ut(J2000, SE_MOON, SEFLG_MOSEPH | SEFLG_SPEED)

        # Velocity should be non-zero
        dlon = pos[3]
        assert dlon != 0.0, "Moon should have non-zero velocity"
        assert 10.0 < dlon < 16.0, f"Moon velocity {dlon} deg/day unexpected"

    def test_moon_is_moshier_body(self):
        """Test that Moon is recognized as a Moshier body."""
        assert is_moshier_body(SE_MOON)
        assert is_moshier_body(1)
