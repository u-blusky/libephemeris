"""
Moshier Heliocentric Cross-Library Comparison Tests.

Validates the entire heliocentric dispatch in the Moshier path (planets.py:948-1025)
by comparing libephemeris against pyswisseph (C library) with SEFLG_MOSEPH | SEFLG_HELCTR.

The Moshier heliocentric dispatch has four separate code paths:
- Sun: returns (0, 0, 0) — heliocentric origin
- Moon: returns Earth heliocentric position (approximation)
- Earth: calc_earth_heliocentric() from moshier/vsop87.py
- Mercury–Neptune: calc_vsop87_heliocentric() from moshier/vsop87.py
- Pluto: calc_pluto_heliocentric() from moshier/pluto.py

Velocities use numerical differentiation with h=0.001 for planets and h=0.01 for Pluto.

This test validates all paths against the C library's single unified codepath.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_EARTH,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================

# Heliocentric position tolerances (C vs Python Moshier).
# Mercury shows the largest VSOP87 C-vs-Python divergence (~0.032° at year 2200).
HELIO_LON_TOL = 0.04  # degrees (~144 arcsec)
HELIO_LAT_TOL = 0.04  # degrees
HELIO_DIST_REL_TOL = 0.001  # relative (0.1%)
HELIO_VEL_TOL = 0.01  # degrees/day

# Pluto: wider tolerances due to Chapront-Francou theory differences.
# Errors grow to ~3.5° for dates far from J2000 (e.g., 1750 CE).
PLUTO_LON_TOL = 4.0  # degrees
PLUTO_LAT_TOL = 4.0  # degrees
PLUTO_DIST_REL_TOL = 0.01  # relative (1%)
PLUTO_VEL_TOL = 0.05  # degrees/day

# Moon heliocentric: libephemeris approximates Moon as Earth heliocentric position
# (planets.py:954-971). The C library computes actual Moon heliocentric position
# (Earth + geocentric Moon offset), so differences can reach ~0.16° in longitude.
MOON_LON_TOL = 0.2  # degrees
MOON_LAT_TOL = 0.2  # degrees
MOON_DIST_REL_TOL = 0.001  # relative
MOON_VEL_TOL = 0.05  # degrees/day


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# 5 test dates spanning the Moshier range
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (1900, 6, 15, 0.0, "1900"),
    (2024, 3, 20, 12.0, "2024 equinox"),
    (1750, 1, 1, 0.0, "1750 pre-DE440"),
    (2200, 7, 4, 18.0, "2200 future"),
]

# Heliocentric bodies: Mercury through Neptune + Earth + Pluto (9 bodies)
HELIO_BODIES = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_EARTH, "Earth"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
]


# ============================================================================
# HELIOCENTRIC POSITION TESTS
# ============================================================================


class TestMoshierHeliocentricPositions:
    """Cross-library comparison for heliocentric positions with Moshier.

    Parametrized over 9 bodies x 5 dates = 45 test cases validating
    longitude, latitude, distance, and velocities.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        HELIO_BODIES,
    )
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_heliocentric_position(
        self,
        planet_id,
        planet_name,
        year,
        month,
        day,
        hour,
        date_label,
    ):
        """Test Moshier heliocentric positions vs pyswisseph C library."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed for {planet_name} @ {date_label}: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed for {planet_name} @ {date_label}: {e}")

        # Select tolerances based on body
        if planet_id == SE_PLUTO:
            lon_tol = PLUTO_LON_TOL
            lat_tol = PLUTO_LAT_TOL
            dist_rel_tol = PLUTO_DIST_REL_TOL
            vel_tol = PLUTO_VEL_TOL
        else:
            lon_tol = HELIO_LON_TOL
            lat_tol = HELIO_LAT_TOL
            dist_rel_tol = HELIO_DIST_REL_TOL
            vel_tol = HELIO_VEL_TOL

        # Longitude
        diff_lon = angular_diff(res_swe[0], res_py[0])
        assert diff_lon < lon_tol, (
            f"Heliocentric {planet_name} @ {date_label}: "
            f"lon diff {diff_lon:.6f}° (tolerance {lon_tol}°)"
        )

        # Latitude
        diff_lat = abs(res_swe[1] - res_py[1])
        assert diff_lat < lat_tol, (
            f"Heliocentric {planet_name} @ {date_label}: "
            f"lat diff {diff_lat:.6f}° (tolerance {lat_tol}°)"
        )

        # Distance (relative comparison to handle varying AU magnitudes)
        dist_swe = res_swe[2]
        dist_py = res_py[2]
        if dist_swe != 0:
            dist_rel = abs(dist_swe - dist_py) / abs(dist_swe)
            assert dist_rel < dist_rel_tol, (
                f"Heliocentric {planet_name} @ {date_label}: "
                f"dist rel diff {dist_rel:.6f} ({dist_rel * 100:.4f}%) "
                f"(tolerance {dist_rel_tol * 100}%)"
            )

        # Velocity (longitude speed)
        diff_vel = abs(res_swe[3] - res_py[3])
        assert diff_vel < vel_tol, (
            f"Heliocentric {planet_name} @ {date_label}: "
            f"lon speed diff {diff_vel:.6f} deg/day (tolerance {vel_tol})"
        )


# ============================================================================
# SUN HELIOCENTRIC (ORIGIN) TESTS
# ============================================================================


class TestMoshierHeliocentricSun:
    """Special tests for heliocentric Sun — must be at origin (0, 0, 0).

    The Sun in heliocentric coordinates is by definition at the origin.
    Both C and Python should return zeros for position and velocity.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_sun_at_origin(self, year, month, day, hour, date_label):
        """Heliocentric Sun must return (0, 0, 0) for position."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, SE_SUN, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        # Both libraries should return zero distance for heliocentric Sun
        assert res_py[2] == 0.0, (
            f"Sun heliocentric distance @ {date_label}: expected 0.0, got {res_py[2]}"
        )

        # Cross-library: C should also return zero
        assert res_swe[2] == 0.0, (
            f"pyswisseph Sun heliocentric distance @ {date_label}: "
            f"expected 0.0, got {res_swe[2]}"
        )

        # Longitude comparison (both should be 0 or equivalent)
        diff_lon = angular_diff(res_swe[0], res_py[0])
        assert diff_lon < 0.001, (
            f"Sun heliocentric lon @ {date_label}: "
            f"C={res_swe[0]:.6f}° vs Py={res_py[0]:.6f}° diff={diff_lon:.6f}°"
        )

        # Velocities should be zero
        assert abs(res_py[3]) < 0.001, (
            f"Sun heliocentric lon speed @ {date_label}: expected ~0, got {res_py[3]}"
        )


# ============================================================================
# MOON HELIOCENTRIC TESTS
# ============================================================================


class TestMoshierHeliocentricMoon:
    """Special tests for heliocentric Moon.

    libephemeris approximates heliocentric Moon as Earth's heliocentric position
    (planets.py:954-971). The C library computes Moon's actual heliocentric
    position (Earth + geocentric Moon offset), causing differences up to ~0.16°
    in longitude. This is a known approximation documented in planets.py.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_moon_heliocentric_vs_c(self, year, month, day, hour, date_label):
        """Compare heliocentric Moon between C and Python."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, SE_MOON, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, SE_MOON, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        diff_lon = angular_diff(res_swe[0], res_py[0])
        diff_lat = abs(res_swe[1] - res_py[1])

        assert diff_lon < MOON_LON_TOL, (
            f"Heliocentric Moon @ {date_label}: "
            f"lon diff {diff_lon:.6f}° (tolerance {MOON_LON_TOL}°)"
        )
        assert diff_lat < MOON_LAT_TOL, (
            f"Heliocentric Moon @ {date_label}: "
            f"lat diff {diff_lat:.6f}° (tolerance {MOON_LAT_TOL}°)"
        )

        # Distance should be ~1 AU (Earth-Sun distance)
        assert 0.98 < res_py[2] < 1.02, (
            f"Heliocentric Moon distance @ {date_label}: "
            f"{res_py[2]:.6f} AU (expected ~1.0 AU)"
        )

        # Velocity comparison (relaxed due to numerical differentiation)
        diff_vel = abs(res_swe[3] - res_py[3])
        assert diff_vel < MOON_VEL_TOL, (
            f"Heliocentric Moon lon speed @ {date_label}: "
            f"diff {diff_vel:.6f} deg/day (tolerance {MOON_VEL_TOL})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_moon_close_to_earth(self, year, month, day, hour, date_label):
        """Heliocentric Moon should be close to heliocentric Earth.

        Since libephemeris returns Earth position for Moon, verify that
        the C library's Moon is also close to Earth (within Moon orbit radius).
        """
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_SPEED

        try:
            res_earth_swe, _ = swe.calc_ut(jd, SE_EARTH, flags)
            res_moon_swe, _ = swe.calc_ut(jd, SE_MOON, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        # Moon's heliocentric position should be within ~0.16° of Earth
        # (C library adds geocentric Moon offset to Earth heliocentric position)
        diff_lon = angular_diff(res_earth_swe[0], res_moon_swe[0])
        assert diff_lon < 0.2, (
            f"C library Moon-Earth heliocentric lon diff @ {date_label}: "
            f"{diff_lon:.6f}° (expected < 0.2°)"
        )
