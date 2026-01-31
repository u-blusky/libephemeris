"""
Planetary Calculations Comparison Tests.

Validates all planetary calculations between pyswisseph and libephemeris:
- All major planets (Sun through Pluto)
- Multiple calculation modes (geocentric, heliocentric, barycentric)
- Position (longitude, latitude, distance) and velocity
- Various time periods and edge cases
"""

import pytest
import swisseph as swe
import libephemeris as ephem
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
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_BARYCTR,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

PLANETS = [
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
]

CALC_MODES = [
    (SEFLG_SWIEPH, "Geocentric"),
    (SEFLG_SWIEPH | SEFLG_SPEED, "Geocentric+Speed"),
    (SEFLG_SWIEPH | SEFLG_HELCTR, "Heliocentric"),
    (SEFLG_SWIEPH | SEFLG_BARYCTR, "Barycentric"),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
    (2100, 12, 31, 23.999, "Late 2100s"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]

# Tolerances
LONGITUDE_STRICT = 0.001  # degrees
LONGITUDE_RELAXED = 0.03  # degrees (heliocentric/barycentric)
LATITUDE_STRICT = 0.001
LATITUDE_RELAXED = 0.03
DISTANCE_STRICT = 0.0001  # AU
DISTANCE_RELAXED = 0.01
VELOCITY_ANGULAR = 0.01  # degrees/day
VELOCITY_RADIAL = 0.001  # AU/day


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestPlanetaryPositions:
    """Compare planetary position calculations between swisseph and libephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_geocentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test geocentric planetary positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = SEFLG_SWIEPH

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_STRICT, (
            f"{planet_name} at {date_desc}: distance diff {diff_dist:.8f} AU exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_geocentric_with_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test geocentric positions with velocity calculations."""
        jd = swe.julday(year, month, day, hour)
        flag = SEFLG_SWIEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag)

        # Position checks
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_STRICT, (
            f"{planet_name} at {date_desc}: distance diff {diff_dist:.8f} AU exceeds tolerance"
        )

        # Velocity checks
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])
        diff_dist_speed = abs(pos_swe[5] - pos_py[5])

        assert diff_lon_speed < VELOCITY_ANGULAR, (
            f"{planet_name} at {date_desc}: lon velocity diff {diff_lon_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_lat_speed < VELOCITY_ANGULAR, (
            f"{planet_name} at {date_desc}: lat velocity diff {diff_lat_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_dist_speed < VELOCITY_RADIAL, (
            f"{planet_name} at {date_desc}: dist velocity diff {diff_dist_speed:.8f} AU/day exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", TEST_DATES[:3]
    )  # Fewer dates for these
    def test_heliocentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test heliocentric planetary positions match pyswisseph."""
        # Skip Sun for heliocentric (meaningless)
        if planet_id == SE_SUN:
            pytest.skip("Sun position is not meaningful in heliocentric coordinates")

        jd = swe.julday(year, month, day, hour)
        flag = SEFLG_SWIEPH | SEFLG_HELCTR

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_RELAXED, (
            f"{planet_name} helio at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_RELAXED, (
            f"{planet_name} helio at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_RELAXED, (
            f"{planet_name} helio at {date_desc}: distance diff {diff_dist:.6f} AU exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES[:3])
    def test_barycentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test barycentric planetary positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = SEFLG_SWIEPH | SEFLG_BARYCTR

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_RELAXED, (
            f"{planet_name} bary at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_RELAXED, (
            f"{planet_name} bary at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_RELAXED, (
            f"{planet_name} bary at {date_desc}: distance diff {diff_dist:.6f} AU exceeds tolerance"
        )


class TestInnerPlanets:
    """Specific tests for inner planets (faster moving, more precision needed)."""

    INNER_PLANETS = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MERCURY, "Mercury"),
        (SE_VENUS, "Venus"),
        (SE_MARS, "Mars"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", INNER_PLANETS)
    def test_inner_planet_daily_motion(self, planet_id, planet_name):
        """Test that daily motion (velocity) matches for inner planets."""
        jd = swe.julday(2024, 6, 15, 12.0)
        flag = SEFLG_SWIEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag)

        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_speed < 0.001, (
            f"{planet_name} daily motion diff {diff_speed:.6f}°/day exceeds tight tolerance"
        )


class TestOuterPlanets:
    """Specific tests for outer planets (slower moving)."""

    OUTER_PLANETS = [
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_URANUS, "Uranus"),
        (SE_NEPTUNE, "Neptune"),
        (SE_PLUTO, "Pluto"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OUTER_PLANETS)
    def test_outer_planet_positions_extended_range(self, planet_id, planet_name):
        """Test outer planets across extended date range."""
        dates = [
            swe.julday(1600, 1, 1, 12.0),
            swe.julday(1800, 1, 1, 12.0),
            swe.julday(2200, 1, 1, 12.0),
            swe.julday(2500, 1, 1, 12.0),
        ]

        for jd in dates:
            pos_swe, _ = swe.calc_ut(jd, planet_id, SEFLG_SWIEPH)
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, SEFLG_SWIEPH)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            assert diff_lon < LONGITUDE_STRICT, (
                f"{planet_name} at JD {jd}: longitude diff {diff_lon:.6f}° exceeds tolerance"
            )


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_j2000_epoch(self, planet_id, planet_name):
        """Test exact J2000.0 epoch positions."""
        jd = 2451545.0  # J2000.0

        pos_swe, _ = swe.calc_ut(jd, planet_id, SEFLG_SWIEPH | SEFLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, SEFLG_SWIEPH | SEFLG_SPEED)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < 0.0001, (
            f"{planet_name} at J2000: longitude diff {diff_lon:.8f}° exceeds tight tolerance"
        )

    @pytest.mark.comparison
    def test_moon_at_midnight(self):
        """Test Moon position at midnight (common use case)."""
        jd = swe.julday(2024, 1, 1, 0.0)

        pos_swe, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < 0.001, (
            f"Moon longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_speed < 0.01, (
            f"Moon speed diff {diff_speed:.6f}°/day exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_sun_at_vernal_equinox(self):
        """Test Sun position at vernal equinox 2024."""
        # 2024 vernal equinox approximately March 20, 03:06 UT
        jd = swe.julday(2024, 3, 20, 3.1)

        pos_swe, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH)

        # Sun should be near 0° (Aries ingress)
        assert pos_swe[0] < 1.0 or pos_swe[0] > 359.0, (
            "Sun should be near 0° at vernal equinox"
        )

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < 0.001, (
            f"Sun at equinox diff {diff_lon:.6f}° exceeds tolerance"
        )
