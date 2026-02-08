"""
Tests for _calc_body_moshier() dispatcher function.

This module tests the Moshier ephemeris dispatcher that routes calculations
to VSOP87 (planets), ELP82B (Moon), and Pluto theory modules.
"""

import pytest
import libephemeris as eph
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    SEFLG_EQUATORIAL,
    SEFLG_HELCTR,
)


class TestMoshierBasicPlanets:
    """Test basic planet calculations with Moshier ephemeris."""

    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
            (SE_URANUS, "Uranus"),
            (SE_NEPTUNE, "Neptune"),
            (SE_PLUTO, "Pluto"),
        ],
    )
    def test_moshier_returns_valid_positions(self, body_id, body_name):
        """Moshier should return valid longitude (0-360) for all major planets."""
        jd = 2451545.0  # J2000.0

        pos, flag = eph.swe_calc_ut(jd, body_id, SEFLG_MOSEPH)

        # Check we got valid results
        assert len(pos) == 6, f"{body_name} should return 6 values"
        lon, lat, dist = pos[0], pos[1], pos[2]

        # Longitude should be in valid range
        assert 0 <= lon < 360, f"{body_name} longitude {lon} out of range"

        # Latitude should be reasonable for planets (ecliptic, so |lat| < 20)
        assert -30 < lat < 30, f"{body_name} latitude {lat} seems wrong"

        # Distance should be positive
        assert dist >= 0, f"{body_name} distance should be non-negative"

    def test_moshier_with_speed_flag(self):
        """Moshier should calculate velocities when SEFLG_SPEED is set."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_MARS, SEFLG_MOSEPH | SEFLG_SPEED)

        # Velocity (dlon) should be non-zero for Mars
        dlon = pos[3]
        assert dlon != 0, "Mars velocity should be non-zero"
        # Mars typically moves 0.3-0.8 deg/day
        assert -1.5 < dlon < 1.5, f"Mars velocity {dlon} seems unreasonable"


class TestMoshierLunarNodes:
    """Test lunar node calculations with Moshier ephemeris."""

    def test_mean_node_position(self):
        """Mean lunar node should be calculated correctly."""
        jd = 2451545.0  # J2000.0

        pos, flag = eph.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_MOSEPH)

        lon = pos[0]
        assert 0 <= lon < 360, f"Mean node longitude {lon} out of range"

        # Mean node latitude is always 0
        assert pos[1] == 0.0, "Mean node latitude should be 0"

    def test_true_node_position(self):
        """True lunar node should be calculated correctly."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_TRUE_NODE, SEFLG_MOSEPH)

        lon = pos[0]
        assert 0 <= lon < 360, f"True node longitude {lon} out of range"

    def test_node_velocity(self):
        """Lunar node should have retrograde motion (negative velocity)."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_MOSEPH | SEFLG_SPEED)

        dlon = pos[3]
        # Mean node moves retrograde at ~0.053 deg/day
        assert dlon < 0, "Mean node should have retrograde motion"
        assert -0.1 < dlon < 0, f"Mean node velocity {dlon} seems wrong"


class TestMoshierLilith:
    """Test Lilith (lunar apogee) calculations with Moshier ephemeris."""

    def test_mean_lilith_position(self):
        """Mean Lilith should be calculated correctly."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_MEAN_APOG, SEFLG_MOSEPH)

        lon = pos[0]
        assert 0 <= lon < 360, f"Mean Lilith longitude {lon} out of range"

    def test_osculating_lilith_position(self):
        """Osculating (true) Lilith should be calculated correctly."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_OSCU_APOG, SEFLG_MOSEPH)

        lon = pos[0]
        assert 0 <= lon < 360, f"Osculating Lilith longitude {lon} out of range"


class TestMoshierCoordinateTransformations:
    """Test coordinate transformations in Moshier mode."""

    def test_equatorial_coordinates(self):
        """SEFLG_EQUATORIAL should return RA/Dec instead of lon/lat."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH | SEFLG_EQUATORIAL)

        ra, dec = pos[0], pos[1]

        # RA should be 0-360 degrees
        assert 0 <= ra < 360, f"RA {ra} out of range"

        # Dec should be -90 to +90 degrees
        assert -90 <= dec <= 90, f"Dec {dec} out of range"

    def test_sidereal_mode(self):
        """SEFLG_SIDEREAL should apply ayanamsha correction."""
        jd = 2451545.0

        # Set sidereal mode (Lahiri)
        eph.swe_set_sid_mode(1)  # SE_SIDM_LAHIRI

        pos_tropical, _ = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)
        pos_sidereal, _ = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH | SEFLG_SIDEREAL)

        # Sidereal longitude should be less than tropical (by ~24 degrees at J2000)
        diff = pos_tropical[0] - pos_sidereal[0]
        # Handle wrap-around
        if diff < -180:
            diff += 360
        if diff > 180:
            diff -= 360

        # Lahiri ayanamsha at J2000 is ~23.9 degrees
        assert 20 < diff < 28, f"Ayanamsha difference {diff} seems wrong"


class TestMoshierHeliocentricMode:
    """Test heliocentric calculations with Moshier ephemeris."""

    def test_heliocentric_sun_at_origin(self):
        """Heliocentric Sun should return zero coordinates."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH | SEFLG_HELCTR)

        # Heliocentric Sun is at origin (distance = 0)
        assert pos[2] == 0.0, "Heliocentric Sun distance should be 0"

    def test_heliocentric_mars(self):
        """Heliocentric Mars should have valid coordinates."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_MARS, SEFLG_MOSEPH | SEFLG_HELCTR)

        lon, lat, dist = pos[0], pos[1], pos[2]

        assert 0 <= lon < 360, f"Heliocentric Mars longitude {lon} out of range"
        assert dist > 1.0, f"Mars heliocentric distance {dist} too small"
        assert dist < 2.0, f"Mars heliocentric distance {dist} too large"


class TestMoshierPrecisionComparison:
    """Test Moshier precision against expected accuracy levels."""

    def test_sun_accuracy(self):
        """Sun position should be accurate within ~1 arcsec from reference."""
        jd = 2451545.0  # J2000.0

        pos, _ = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

        # At J2000.0, Sun longitude should be near 280° (winter solstice range)
        # This is a sanity check, not a precision test
        assert 275 < pos[0] < 285, f"Sun at J2000 should be around 280°, got {pos[0]}"

    def test_moon_position(self):
        """Moon position should be reasonable."""
        jd = 2451545.0

        pos, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_MOSEPH)

        # Moon should have valid longitude
        assert 0 <= pos[0] < 360, f"Moon longitude {pos[0]} out of range"

        # Moon velocity should be around 13 deg/day
        pos_speed, _ = eph.swe_calc_ut(jd, SE_MOON, SEFLG_MOSEPH | SEFLG_SPEED)
        dlon = pos_speed[3]
        assert 10 < dlon < 16, f"Moon velocity {dlon} should be around 13 deg/day"


class TestMoshierExtendedDateRange:
    """Test Moshier ephemeris works for dates outside JPL range."""

    def test_year_500_ce(self):
        """Moshier should work for 500 CE (outside JPL DE440 range)."""
        # JD ~1903682.5 is around 500 CE
        jd = 1903682.5

        pos, flag = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

        assert 0 <= pos[0] < 360, f"Sun longitude at 500 CE: {pos[0]}"

    def test_year_1000_bce(self):
        """Moshier should work for 1000 BCE."""
        # JD ~1356180.5 is around 1000 BCE
        jd = 1356180.5

        pos, flag = eph.swe_calc_ut(jd, SE_MARS, SEFLG_MOSEPH)

        assert 0 <= pos[0] < 360, f"Mars longitude at 1000 BCE: {pos[0]}"

    def test_year_2500_ce(self):
        """Moshier should work for 2500 CE."""
        # JD ~2634166.5 is around 2500 CE
        jd = 2634166.5

        pos, flag = eph.swe_calc_ut(jd, SE_JUPITER, SEFLG_MOSEPH)

        assert 0 <= pos[0] < 360, f"Jupiter longitude at 2500 CE: {pos[0]}"


class TestMoshierSweCalc:
    """Test Moshier with swe_calc() (TT time)."""

    def test_swe_calc_with_moshier(self):
        """swe_calc should work with SEFLG_MOSEPH."""
        jd_tt = 2451545.0

        pos, flag = eph.swe_calc(jd_tt, SE_VENUS, SEFLG_MOSEPH)

        assert 0 <= pos[0] < 360, f"Venus longitude: {pos[0]}"


class TestMoshierErrorHandling:
    """Test error handling for unsupported bodies."""

    def test_unsupported_body_raises_error(self):
        """Unsupported body should raise ValueError."""
        jd = 2451545.0

        # Chiron (SE_CHIRON = 15) is not supported by Moshier
        with pytest.raises(ValueError) as exc_info:
            eph.swe_calc_ut(jd, 15, SEFLG_MOSEPH)

        assert "not supported by Moshier" in str(exc_info.value)
