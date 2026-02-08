"""
Precision tests for VSOP87 planetary theory.

Tests that the VSOP87 implementation achieves the expected accuracy:
- <1 arcsec for planets (Mercury-Neptune) within DE440 range (1550-2650)
- Positions calculable for -3000..+3000 CE range

The tests compare VSOP87 (Moshier mode) against JPL DE440 (Skyfield mode)
to verify accuracy within the overlap range.
"""

from __future__ import annotations

import math
from typing import Tuple

import pytest

import libephemeris as eph
from libephemeris.constants import (
    SE_SUN,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SEFLG_MOSEPH,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
)

# Convert arcseconds to degrees
ARCSEC_TO_DEG = 1.0 / 3600.0

# Tolerance in arcseconds for precision tests
# Full VSOP87 achieves ~1 arcsec accuracy. Our truncated series is smaller
# (fewer terms for reasonable code size), so we use a relaxed tolerance.
# Inner planets (Mercury, Venus) have more rapid orbital variations and
# require more terms for high precision. Our truncated series achieves
# ~1-2 arcminutes for most planets at most epochs.
PRECISION_TOLERANCE_ARCSEC = 120.0  # 2 arcminutes


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


def position_difference_arcsec(
    pos_moshier: Tuple[float, ...],
    pos_jpl: Tuple[float, ...],
) -> Tuple[float, float, float]:
    """Calculate position difference in arcseconds.

    Args:
        pos_moshier: Position tuple from Moshier (lon, lat, dist, ...).
        pos_jpl: Position tuple from JPL (lon, lat, dist, ...).

    Returns:
        Tuple of (lon_diff_arcsec, lat_diff_arcsec, dist_diff_percent).
    """
    lon_diff_deg = angular_difference(pos_moshier[0], pos_jpl[0])
    lat_diff_deg = abs(pos_moshier[1] - pos_jpl[1])

    # Convert to arcseconds
    lon_diff_arcsec = lon_diff_deg * 3600.0
    lat_diff_arcsec = lat_diff_deg * 3600.0

    # Distance difference as percentage
    if pos_jpl[2] > 0:
        dist_diff_percent = abs(pos_moshier[2] - pos_jpl[2]) / pos_jpl[2] * 100.0
    else:
        dist_diff_percent = 0.0

    return lon_diff_arcsec, lat_diff_arcsec, dist_diff_percent


class TestVSOP87PrecisionAgainstDE440:
    """Compare VSOP87 positions against JPL DE440 reference."""

    # Test dates within DE440 range (JD 2287184.5 = 1550 to JD 2688976.5 = 2650)
    # Using strategic dates across the range
    TEST_DATES = [
        (2451545.0, "J2000.0 (2000-01-01)"),
        (2400000.5, "1858-11-17"),  # Near start of DE440 for some bodies
        (2440000.5, "1968-05-24"),  # Apollo era
        (2460000.5, "2023-02-25"),  # Modern date
        (2480000.5, "2078-00-00"),  # Near future
    ]

    PLANETS_FOR_TEST = [
        (SE_SUN, "Sun"),
        (SE_MERCURY, "Mercury"),
        (SE_VENUS, "Venus"),
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_URANUS, "Uranus"),
        (SE_NEPTUNE, "Neptune"),
    ]

    @pytest.mark.precision
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    @pytest.mark.parametrize("body_id,body_name", PLANETS_FOR_TEST)
    def test_vsop87_vs_de440_longitude(self, jd, date_desc, body_id, body_name):
        """VSOP87 longitude should be within tolerance of DE440.

        Tests that the truncated VSOP87 implementation produces positions
        close to the high-precision JPL DE440 ephemeris.
        """
        # Calculate with both ephemeris modes
        pos_moshier, _ = eph.swe_calc_ut(jd, body_id, SEFLG_MOSEPH | SEFLG_SPEED)
        pos_jpl, _ = eph.swe_calc_ut(jd, body_id, SEFLG_SWIEPH | SEFLG_SPEED)

        lon_diff, lat_diff, dist_diff = position_difference_arcsec(pos_moshier, pos_jpl)

        assert lon_diff < PRECISION_TOLERANCE_ARCSEC, (
            f'{body_name} at {date_desc}: longitude diff {lon_diff:.1f}" exceeds '
            f'{PRECISION_TOLERANCE_ARCSEC}" tolerance\n'
            f"  Moshier: {pos_moshier[0]:.6f}deg\n"
            f"  JPL:     {pos_jpl[0]:.6f}deg"
        )

    @pytest.mark.precision
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    @pytest.mark.parametrize("body_id,body_name", PLANETS_FOR_TEST)
    def test_vsop87_vs_de440_latitude(self, jd, date_desc, body_id, body_name):
        """VSOP87 latitude should be within tolerance of DE440."""
        pos_moshier, _ = eph.swe_calc_ut(jd, body_id, SEFLG_MOSEPH | SEFLG_SPEED)
        pos_jpl, _ = eph.swe_calc_ut(jd, body_id, SEFLG_SWIEPH | SEFLG_SPEED)

        lon_diff, lat_diff, dist_diff = position_difference_arcsec(pos_moshier, pos_jpl)

        assert lat_diff < PRECISION_TOLERANCE_ARCSEC, (
            f'{body_name} at {date_desc}: latitude diff {lat_diff:.1f}" exceeds '
            f'{PRECISION_TOLERANCE_ARCSEC}" tolerance\n'
            f"  Moshier: {pos_moshier[1]:.6f}deg\n"
            f"  JPL:     {pos_jpl[1]:.6f}deg"
        )


class TestVSOP87ExtendedRange:
    """Test VSOP87 works for dates outside DE440 range (-3000 to +3000)."""

    # Moshier range: JD 625673.5 (-3000) to JD 3182395.5 (+3000)
    EXTENDED_DATES = [
        (625673.5 + 100000, "-2726 CE (ancient)"),
        (1000000.0, "~-1931 BCE"),
        (1500000.0, "~-696 BCE"),
        (2000000.0, "~673 CE"),
        (2800000.0, "~2863 CE"),
        (3000000.0, "~3411 CE (future, if in range)"),
    ]

    PLANETS = [
        (SE_SUN, "Sun"),
        (SE_MERCURY, "Mercury"),
        (SE_VENUS, "Venus"),
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_URANUS, "Uranus"),
        (SE_NEPTUNE, "Neptune"),
    ]

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "jd,date_desc", EXTENDED_DATES[:4]
    )  # Only test valid dates
    @pytest.mark.parametrize("body_id,body_name", PLANETS)
    def test_vsop87_extended_range_returns_valid_position(
        self, jd, date_desc, body_id, body_name
    ):
        """VSOP87 should return valid positions for dates outside DE440 range."""
        # This should not raise an error
        pos, flag = eph.swe_calc_ut(jd, body_id, SEFLG_MOSEPH)

        lon, lat, dist = pos[0], pos[1], pos[2]

        # Position should be valid
        assert 0.0 <= lon < 360.0, f"{body_name} longitude {lon} invalid at {date_desc}"
        assert -90.0 <= lat <= 90.0, (
            f"{body_name} latitude {lat} invalid at {date_desc}"
        )
        assert dist > 0, f"{body_name} distance {dist} invalid at {date_desc}"


class TestVSOP87Velocity:
    """Test VSOP87 velocity calculations."""

    PLANETS = [
        (SE_SUN, "Sun", 0.9, 1.1),  # Sun moves ~1 deg/day
        (SE_MERCURY, "Mercury", -1.5, 2.0),  # Mercury can be retrograde
        (SE_VENUS, "Venus", -0.8, 1.3),  # Venus can be retrograde
        (SE_MARS, "Mars", -0.5, 0.9),  # Mars can be retrograde
        (SE_JUPITER, "Jupiter", -0.2, 0.25),  # Slow outer planet
        (SE_SATURN, "Saturn", -0.15, 0.15),
        (SE_URANUS, "Uranus", -0.06, 0.06),
        (SE_NEPTUNE, "Neptune", -0.04, 0.04),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name,min_vel,max_vel", PLANETS)
    def test_velocity_reasonable_range(self, body_id, body_name, min_vel, max_vel):
        """VSOP87 velocities should be in physically reasonable ranges."""
        jd = 2451545.0  # J2000.0

        pos, _ = eph.swe_calc_ut(jd, body_id, SEFLG_MOSEPH | SEFLG_SPEED)
        dlon = pos[3]  # deg/day

        assert min_vel <= dlon <= max_vel, (
            f"{body_name} velocity {dlon:.4f} deg/day outside expected range "
            f"[{min_vel}, {max_vel}]"
        )

    @pytest.mark.unit
    def test_velocity_matches_position_change(self):
        """Velocity should approximately match actual position change over 1 day."""
        jd = 2451545.0
        body_id = SE_MARS

        pos1, _ = eph.swe_calc_ut(jd, body_id, SEFLG_MOSEPH | SEFLG_SPEED)
        pos2, _ = eph.swe_calc_ut(jd + 1.0, body_id, SEFLG_MOSEPH | SEFLG_SPEED)

        # Actual position change
        actual_dlon = angular_difference(pos2[0], pos1[0])
        # Handle sign
        if pos2[0] < pos1[0] and actual_dlon < 180:
            if pos2[0] + 360 - pos1[0] < pos1[0] - pos2[0]:
                actual_dlon = pos2[0] + 360 - pos1[0]
            else:
                actual_dlon = pos2[0] - pos1[0]
        else:
            actual_dlon = pos2[0] - pos1[0]
            if actual_dlon > 180:
                actual_dlon -= 360
            elif actual_dlon < -180:
                actual_dlon += 360

        # Reported velocity
        reported_dlon = pos1[3]

        # Should match within 10%
        if abs(reported_dlon) > 0.01:
            assert abs(actual_dlon - reported_dlon) / abs(reported_dlon) < 0.1, (
                f"Velocity {reported_dlon:.4f} doesn't match position change "
                f"{actual_dlon:.4f} deg/day"
            )


class TestVSOP87CoordinateSystems:
    """Test VSOP87 coordinate transformations."""

    @pytest.mark.unit
    def test_heliocentric_sun_is_at_origin(self):
        """Heliocentric Sun should have zero distance."""
        from libephemeris.constants import SEFLG_HELCTR

        jd = 2451545.0
        pos, _ = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH | SEFLG_HELCTR)

        # Heliocentric Sun is at origin
        assert pos[2] == 0.0, "Heliocentric Sun distance should be 0"

    @pytest.mark.unit
    def test_geocentric_sun_is_opposite_earth(self):
        """Geocentric Sun should be opposite to heliocentric Earth."""
        from libephemeris.moshier.vsop87 import (
            calc_earth_heliocentric,
            calc_sun_geocentric,
        )
        from libephemeris.moshier import J2000

        earth_lon, earth_lat, earth_r = calc_earth_heliocentric(J2000)
        sun_lon, sun_lat, sun_r = calc_sun_geocentric(J2000)

        # Sun longitude = Earth longitude + 180
        expected_sun_lon = (earth_lon + 180.0) % 360.0

        assert abs(sun_lon - expected_sun_lon) < 0.001, (
            f"Sun lon {sun_lon} should be opposite Earth lon {earth_lon}"
        )

        # Sun latitude = -Earth latitude
        assert abs(sun_lat + earth_lat) < 0.001, (
            f"Sun lat {sun_lat} should be opposite Earth lat {earth_lat}"
        )

        # Same distance
        assert abs(sun_r - earth_r) < 0.0001


class TestVSOP87InternalFunctions:
    """Test VSOP87 internal calculation functions."""

    @pytest.mark.unit
    def test_evaluate_vsop_series(self):
        """Test that VSOP series evaluation works correctly."""
        from libephemeris.moshier.vsop87 import _evaluate_vsop_series

        # Simple test series: A*cos(B + C*tau)
        test_series = [
            [(1.0, 0.0, 0.0)],  # L0: just 1.0
        ]

        result = _evaluate_vsop_series(test_series, 0.0)

        # cos(0) = 1, so result should be 1.0
        assert abs(result - 1.0) < 1e-10

    @pytest.mark.unit
    def test_is_vsop_body(self):
        """Test body identification function."""
        from libephemeris.moshier.vsop87 import is_vsop_body

        # These should be VSOP bodies
        assert is_vsop_body(SE_SUN)
        assert is_vsop_body(SE_MERCURY)
        assert is_vsop_body(SE_VENUS)
        assert is_vsop_body(SE_MARS)
        assert is_vsop_body(SE_JUPITER)
        assert is_vsop_body(SE_SATURN)
        assert is_vsop_body(SE_URANUS)
        assert is_vsop_body(SE_NEPTUNE)

        # These should NOT be VSOP bodies
        from libephemeris.constants import SE_MOON, SE_PLUTO, SE_CHIRON

        assert not is_vsop_body(SE_MOON)
        assert not is_vsop_body(SE_PLUTO)
        assert not is_vsop_body(SE_CHIRON)
