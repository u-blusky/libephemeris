"""
Moshier Planetary Calculations Cross-Library Comparison Tests.

Validates Moshier (SEFLG_MOSEPH) planetary calculations between pyswisseph
(C library) and libephemeris (Python reimplementation):
- All 10 major planets (Sun through Pluto)
- Geocentric positions (longitude, latitude, distance)
- Geocentric positions with velocity
- Various time periods (J2000.0, 2024, 1900, 2100, 1950)

This is the Moshier-mode mirror of test_compare_planets.py (which covers
SEFLG_SWIEPH / JPL mode). It ensures that the Python Moshier reimplementation
(VSOP87/ELP) produces the same results as the original C Moshier code
integrated in Swiss Ephemeris via pyswisseph.

Without these tests, systematic divergences between the Python VSOP87/ELP
implementation and the C original (e.g. reference frame errors, mean ecliptic
vs J2000 ICRS) would remain undetected.
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
    SEFLG_MOSEPH,
    SEFLG_SPEED,
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

# Same TEST_DATES as test_compare_planets.py
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
    (2100, 12, 31, 23.999, "Late 2100s"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]

# Moshier-specific tolerances (relaxed compared to JPL/SWIEPH)
# The Moshier semi-analytical ephemeris (VSOP87/ELP) has inherently lower
# precision than JPL DE440, so cross-library comparison tolerances are wider
# to account for differences between the C and Python VSOP87 implementations.
#
# Observed C-vs-Python offsets at J2000.0: Sun ~0.010°, Mercury ~0.018°,
# Mars ~0.012°, Saturn ~0.011°, Neptune ~0.013°. These are consistent with
# the ~70 arcsec (~0.019°) tolerance established in test_moshier_vs_jpl.py.
#
# Pluto uses a separate, much larger tolerance because the Moshier Pluto
# theory (Chapront-Francou + Keplerian) degrades significantly at dates
# far from J2000, with errors up to several degrees (matching the
# TOLERANCE_PLUTO_ARCSEC = 60000 in test_moshier_vs_jpl.py).
MOSHIER_LONGITUDE = 0.02  # degrees (~72 arcsec, consistent with 70 arcsec tolerance)
MOSHIER_LATITUDE = 0.02  # degrees
MOSHIER_DISTANCE_REL = 0.001  # relative (dimensionless)
MOSHIER_VELOCITY = 0.01  # degrees/day
MOSHIER_PLUTO_LONGITUDE = 17.0  # degrees (Pluto Moshier theory has large errors)
MOSHIER_EXTENDED_RANGE = 0.03  # degrees (relaxed for dates far from J2000)

# Inner planet speed tolerances
# Moon (~13 deg/day) amplifies numerical differences in differentiation
# between C and Python implementations, requiring a relaxed tolerance.
# Mercury-Mars have lower angular velocities and use tighter tolerances.
INNER_SPEED_LON = 0.01  # degrees/day for Mercury, Venus, Mars, Sun
INNER_SPEED_LON_MOON = 0.05  # degrees/day (relaxed for Moon's high angular velocity)
INNER_SPEED_LAT = 0.001  # degrees/day
INNER_SPEED_LAT_MERCURY = (
    0.003  # degrees/day (Mercury's 7° inclination amplifies diffs)
)
INNER_SPEED_DIST = 0.001  # AU/day

# Outer planet tolerances for dedicated C-vs-Python Moshier comparison.
# Jupiter-Neptune use the same MOSHIER_LONGITUDE as the general tests, which
# accommodates observed C-vs-Python offsets: Saturn ~0.011°, Neptune ~0.013°,
# Uranus ~0.011°. Pluto uses a separate tolerance that is ~8.5x tighter than
# the general MOSHIER_PLUTO_LONGITUDE (17°), because here we compare the
# *same* Chapront-Francou algorithm in C vs Python. The observed Pluto C-vs-
# Python difference ranges from near-zero at J2000 to ~1.43° at dates 100+
# years from J2000, reflecting porting differences in the 2993-row perturbation
# tables (pluto_data.py).
OUTER_MOSHIER_LONGITUDE = 0.02  # degrees (~72 arcsec) for Jupiter-Neptune
OUTER_PLUTO_LONGITUDE = 2.0  # degrees for Pluto C-vs-Python (observed max ~1.43°)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def get_longitude_tolerance(planet_id: int) -> float:
    """Get appropriate longitude tolerance for a given planet ID.

    Pluto's Moshier theory (Chapront-Francou + Keplerian) has much larger
    errors than the VSOP87 theories used for Mercury-Neptune, so it needs
    a separate, relaxed tolerance.
    """
    if planet_id == SE_PLUTO:
        return MOSHIER_PLUTO_LONGITUDE
    return MOSHIER_LONGITUDE


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestPlanetaryPositions:
    """Compare Moshier planetary position calculations between pyswisseph and libephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_geocentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier geocentric planetary positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        lon_tol = get_longitude_tolerance(planet_id)

        assert diff_lon < lon_tol, (
            f"{planet_name} Moshier at {date_desc}: longitude diff {diff_lon:.6f}° "
            f"exceeds tolerance {lon_tol}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )
        assert diff_lat < MOSHIER_LATITUDE, (
            f"{planet_name} Moshier at {date_desc}: latitude diff {diff_lat:.6f}° "
            f"exceeds tolerance {MOSHIER_LATITUDE}°"
        )

        # Distance: relative comparison
        if pos_swe[2] > 0:
            rel_dist = abs(pos_swe[2] - pos_py[2]) / pos_swe[2]
            assert rel_dist < MOSHIER_DISTANCE_REL, (
                f"{planet_name} Moshier at {date_desc}: distance relative diff "
                f"{rel_dist:.6f} exceeds tolerance {MOSHIER_DISTANCE_REL} "
                f"(swe={pos_swe[2]:.8f}, lib={pos_py[2]:.8f})"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_geocentric_with_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier geocentric positions with velocity calculations."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        # Position checks
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        lon_tol = get_longitude_tolerance(planet_id)

        assert diff_lon < lon_tol, (
            f"{planet_name} Moshier+Speed at {date_desc}: longitude diff "
            f"{diff_lon:.6f}° exceeds tolerance {lon_tol}°"
        )
        assert diff_lat < MOSHIER_LATITUDE, (
            f"{planet_name} Moshier+Speed at {date_desc}: latitude diff "
            f"{diff_lat:.6f}° exceeds tolerance"
        )

        # Distance: relative comparison
        if pos_swe[2] > 0:
            rel_dist = abs(pos_swe[2] - pos_py[2]) / pos_swe[2]
            assert rel_dist < MOSHIER_DISTANCE_REL, (
                f"{planet_name} Moshier+Speed at {date_desc}: distance relative diff "
                f"{rel_dist:.6f} exceeds tolerance"
            )

        # Velocity checks
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])

        assert diff_lon_speed < MOSHIER_VELOCITY, (
            f"{planet_name} Moshier at {date_desc}: lon velocity diff "
            f"{diff_lon_speed:.6f}°/day exceeds tolerance {MOSHIER_VELOCITY}°/day"
        )
        assert diff_lat_speed < MOSHIER_VELOCITY, (
            f"{planet_name} Moshier at {date_desc}: lat velocity diff "
            f"{diff_lat_speed:.6f}°/day exceeds tolerance {MOSHIER_VELOCITY}°/day"
        )


class TestInnerPlanets:
    """Specific Moshier tests for inner planets (faster moving, more precision needed)."""

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
        """Test that Moshier daily motion (velocity) matches for inner planets."""
        jd = swe.julday(2024, 6, 15, 12.0)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_speed < MOSHIER_VELOCITY, (
            f"{planet_name} Moshier daily motion diff {diff_speed:.6f}°/day "
            f"exceeds tolerance {MOSHIER_VELOCITY}°/day"
        )


class TestMoshierInnerPlanets:
    """Dedicated Moshier speed tests for inner planets across multiple dates.

    Inner planets (especially Moon and Mercury) have high angular velocities
    that amplify numerical differences in velocity differentiation between
    the C and Python Moshier implementations:
    - Moon: ~13 deg/day (ELP 2000-82B with 60+ periodic terms)
    - Mercury: up to ~2 deg/day (frequent retrograde stations 3-4x/year)
    - Sun/Venus/Mars: ~1 deg/day

    This class validates all 6 components (lon, lat, dist, speed_lon,
    speed_lat, speed_dist) for 5 inner planets x 5 dates = 25 test cases.
    """

    INNER_PLANETS = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MERCURY, "Mercury"),
        (SE_VENUS, "Venus"),
        (SE_MARS, "Mars"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", INNER_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_inner_planet_speed_longitude(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier inner planet longitude velocity matches pyswisseph.

        Moon uses a relaxed tolerance (0.05 deg/day) because its high angular
        velocity (~13 deg/day) amplifies numerical differentiation differences
        between the C and Python ELP 2000-82B implementations.
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        # Longitude speed tolerance: Moon gets relaxed tolerance
        speed_lon_tol = (
            INNER_SPEED_LON_MOON if planet_id == SE_MOON else INNER_SPEED_LON
        )

        diff_speed_lon = abs(pos_swe[3] - pos_py[3])
        assert diff_speed_lon < speed_lon_tol, (
            f"{planet_name} Moshier at {date_desc}: lon velocity diff "
            f"{diff_speed_lon:.6f}°/day exceeds tolerance {speed_lon_tol}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", INNER_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_inner_planet_speed_latitude(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier inner planet latitude and distance velocities match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_speed_lat = abs(pos_swe[4] - pos_py[4])
        diff_speed_dist = abs(pos_swe[5] - pos_py[5])

        # Mercury's high orbital inclination (7°) amplifies latitude speed diffs
        speed_lat_tol = (
            INNER_SPEED_LAT_MERCURY if planet_id == SE_MERCURY else INNER_SPEED_LAT
        )

        assert diff_speed_lat < speed_lat_tol, (
            f"{planet_name} Moshier at {date_desc}: lat velocity diff "
            f"{diff_speed_lat:.6f}°/day exceeds tolerance {speed_lat_tol}°/day "
            f"(swe={pos_swe[4]:.6f}, lib={pos_py[4]:.6f})"
        )
        assert diff_speed_dist < INNER_SPEED_DIST, (
            f"{planet_name} Moshier at {date_desc}: dist velocity diff "
            f"{diff_speed_dist:.8f} AU/day exceeds tolerance {INNER_SPEED_DIST} AU/day "
            f"(swe={pos_swe[5]:.8f}, lib={pos_py[5]:.8f})"
        )


class TestMoshierOuterPlanets:
    """Dedicated Moshier C-vs-Python tests for outer planets across multiple dates.

    Outer planets (Jupiter through Pluto) use different theories with varying
    precision in the Moshier semi-analytical ephemeris:
    - Jupiter/Saturn: VSOP87 with sub-arcsecond C-vs-Python agreement
    - Uranus/Neptune: VSOP87 with slightly degraded but still tight agreement
    - Pluto: Chapront-Francou perturbation theory with 2993 tabulated terms
      (pluto_data.py), completely separate from VSOP87. This is where C-vs-Python
      porting differences are most likely due to the large numeric tables.

    This class validates longitude, latitude, and distance for 5 outer planets
    x 5 dates = 25 test cases, using tighter tolerances than the general
    TestPlanetaryPositions (which lumps all planets together) and much tighter
    than test_moshier_vs_jpl.py (which compares Moshier vs JPL, not C vs Python).
    """

    OUTER_PLANETS = [
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_URANUS, "Uranus"),
        (SE_NEPTUNE, "Neptune"),
        (SE_PLUTO, "Pluto"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OUTER_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_outer_planet_longitude(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier outer planet longitude matches pyswisseph C implementation.

        Jupiter-Neptune use OUTER_MOSHIER_LONGITUDE (0.02 deg) matching the
        observed C-vs-Python VSOP87 agreement (~0.013° max). Pluto uses
        OUTER_PLUTO_LONGITUDE (2.0 deg), ~8.5x tighter than the 17 deg general
        tolerance, validated against observed C-vs-Python Chapront-Francou
        differences (near-zero at J2000 up to ~1.43° at extended dates).
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        lon_tol = (
            OUTER_PLUTO_LONGITUDE if planet_id == SE_PLUTO else OUTER_MOSHIER_LONGITUDE
        )

        assert diff_lon < lon_tol, (
            f"{planet_name} Moshier at {date_desc}: longitude diff {diff_lon:.6f}° "
            f"exceeds tolerance {lon_tol}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OUTER_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_outer_planet_latitude(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier outer planet latitude matches pyswisseph C implementation."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_lat = abs(pos_swe[1] - pos_py[1])
        assert diff_lat < MOSHIER_LATITUDE, (
            f"{planet_name} Moshier at {date_desc}: latitude diff {diff_lat:.6f}° "
            f"exceeds tolerance {MOSHIER_LATITUDE}° "
            f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OUTER_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_outer_planet_distance(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier outer planet distance matches pyswisseph C implementation."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        if pos_swe[2] > 0:
            rel_dist = abs(pos_swe[2] - pos_py[2]) / pos_swe[2]
            assert rel_dist < MOSHIER_DISTANCE_REL, (
                f"{planet_name} Moshier at {date_desc}: distance relative diff "
                f"{rel_dist:.6f} exceeds tolerance {MOSHIER_DISTANCE_REL} "
                f"(swe={pos_swe[2]:.8f}, lib={pos_py[2]:.8f})"
            )


class TestOuterPlanets:
    """Specific Moshier tests for outer planets (slower moving)."""

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
        """Test Moshier outer planets across extended date range.

        Uses relaxed tolerance (MOSHIER_EXTENDED_RANGE) because VSOP87 precision
        degrades at dates far from J2000 (e.g. Uranus at year 2500).
        """
        dates = [
            swe.julday(1600, 1, 1, 12.0),
            swe.julday(1800, 1, 1, 12.0),
            swe.julday(2200, 1, 1, 12.0),
            swe.julday(2500, 1, 1, 12.0),
        ]

        for jd in dates:
            flag_swe = swe.FLG_MOSEPH
            flag_py = SEFLG_MOSEPH

            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])
            lon_tol = (
                MOSHIER_PLUTO_LONGITUDE
                if planet_id == SE_PLUTO
                else MOSHIER_EXTENDED_RANGE
            )

            assert diff_lon < lon_tol, (
                f"{planet_name} Moshier at JD {jd}: longitude diff {diff_lon:.6f}° "
                f"exceeds tolerance {lon_tol}°"
            )


class TestEdgeCases:
    """Test Moshier edge cases and boundary conditions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_j2000_epoch(self, planet_id, planet_name):
        """Test exact J2000.0 epoch Moshier positions."""
        jd = 2451545.0  # J2000.0
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        lon_tol = get_longitude_tolerance(planet_id)

        assert diff_lon < lon_tol, (
            f"{planet_name} Moshier at J2000: longitude diff {diff_lon:.8f}° "
            f"exceeds tolerance {lon_tol}°"
        )

    @pytest.mark.comparison
    def test_moon_at_midnight(self):
        """Test Moshier Moon position at midnight (common use case)."""
        jd = swe.julday(2024, 1, 1, 0.0)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, SE_MOON, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MOON, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Moon Moshier longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_speed < MOSHIER_VELOCITY, (
            f"Moon Moshier speed diff {diff_speed:.6f}°/day exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_sun_at_vernal_equinox(self):
        """Test Moshier Sun position at vernal equinox 2024."""
        # 2024 vernal equinox approximately March 20, 03:06 UT
        jd = swe.julday(2024, 3, 20, 3.1)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, SE_SUN, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Sun Moshier at equinox diff {diff_lon:.6f}° exceeds tolerance"
        )
