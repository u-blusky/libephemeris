"""
Moshier TT Path (swe_calc) Cross-Library Comparison Tests.

Validates Moshier (SEFLG_MOSEPH) planetary calculations via the TT (Terrestrial
Time) path between pyswisseph (C library) and libephemeris (Python reimplementation).

This is complementary to test_moshier_compare_planets.py which uses swe_calc_ut
(UT1 path). The TT path exercises the inverse conversion jd_ut = jd_tt - delta_t
in libephemeris's _calc_body_moshier(is_ut=False) at planets.py:849-852, which
can diverge from the C library for historical dates where Delta T is large
(>7000 seconds before 1600).

Why this matters:
- Astronomical/scientific applications prefer TT over UT1
- swe_calc() is the native TT interface
- The inverse TT->UT estimate (jd_ut = jd_tt - swe_deltat(jd_tt)) is an
  approximation that amplifies for dates with large Delta T
- test_moshier_compare_planets.py only tests swe_calc_ut (UT1 path)
- test_moshier_routing.py only verifies swe_calc works, not cross-library match

Test coverage:
- 5 planets (Sun, Moon, Mars, Jupiter, Saturn) x 5 TT dates = 25 position tests
  comparing all 6 components (lon, lat, dist, speed_lon, speed_lat, speed_dist)
- 5 planets x 2 dates = 10 consistency tests validating
  swe_calc_ut(jd_ut) ~= swe_calc(jd_ut + deltat) for both implementations
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Dates in TT, including historical dates with large Delta T.
# swe.julday() is a pure calendar->JD conversion (no time scale semantics).
# When passed to swe.calc() / swe_calc(), the JD is interpreted as TT.
TT_DATES = [
    (2000, 1, 1, 12.0, "J2000.0 (DeltaT ~64s)"),
    (2024, 6, 15, 12.0, "Modern (DeltaT ~69s)"),
    (1000, 1, 1, 12.0, "Medieval (DeltaT ~1574s)"),
    (0, 1, 1, 12.0, "0 CE (DeltaT ~16800s)"),
    (-1000, 1, 1, 12.0, "-1000 CE (DeltaT ~46000s)"),
]

# Cross-library tolerances for the TT path.
# These match test_moshier_compare_planets.py since both C and Python use the
# same inverse Delta T approximation (jd_ut ~= jd_tt - deltat(jd_tt)).
MOSHIER_LONGITUDE = 0.02  # degrees (~72 arcsec)
MOSHIER_LATITUDE = 0.02  # degrees
MOSHIER_DISTANCE_REL = 0.001  # relative (dimensionless)
MOSHIER_SPEED_LON = 0.01  # degrees/day
MOSHIER_SPEED_LON_MOON = 0.05  # degrees/day (Moon ~13 deg/day)
MOSHIER_SPEED_LAT = 0.01  # degrees/day
MOSHIER_SPEED_DIST = 0.001  # AU/day

# For pre-modern dates (year < 1500), VSOP87/ELP polynomial evaluations
# span larger time intervals from J2000, potentially amplifying C-vs-Python
# floating-point divergence in the series summations and the inverse Delta T
# conversion where Delta T exceeds hundreds of seconds.
# Outer planets (Saturn) at -1000 CE can diverge by ~0.23° in longitude.
HISTORICAL_LONGITUDE = 0.3  # degrees (~1080 arcsec)
HISTORICAL_LATITUDE = 0.05  # degrees
HISTORICAL_DISTANCE_REL = 0.005  # relative
HISTORICAL_SPEED_LON = 0.1  # degrees/day
HISTORICAL_SPEED_LAT = 0.1  # degrees/day
HISTORICAL_SPEED_DIST = 0.01  # AU/day

# Consistency test tolerances: swe_calc_ut(jd_ut) ~= swe_calc(jd_ut + deltat)
# Within each library, both paths should yield near-identical results.
CONSISTENCY_LONGITUDE = 0.001  # degrees (~3.6 arcsec)
CONSISTENCY_LATITUDE = 0.001  # degrees
CONSISTENCY_DISTANCE_REL = 0.0001  # relative

# Dates for consistency tests (moderate Delta T, within Moshier range)
CONSISTENCY_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (1900, 1, 1, 12.0, "Year 1900"),
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def is_deep_historical(year: int) -> bool:
    """Check if a date is pre-modern (year < 1500, significant Delta T).

    Before 1500 CE, Delta T exceeds ~100s and grows rapidly for earlier dates.
    VSOP87/ELP polynomial extrapolation also degrades beyond ~1000 years from
    J2000, and the inverse TT->UT conversion amplifies C-vs-Python divergence.
    Medieval (1000 CE, DeltaT ~1574s) and ancient dates need relaxed tolerances.
    """
    return year < 1500


def get_lon_tolerance(year: int) -> float:
    """Get longitude tolerance based on date era."""
    return HISTORICAL_LONGITUDE if is_deep_historical(year) else MOSHIER_LONGITUDE


def get_lat_tolerance(year: int) -> float:
    """Get latitude tolerance based on date era."""
    return HISTORICAL_LATITUDE if is_deep_historical(year) else MOSHIER_LATITUDE


def get_dist_tolerance(year: int) -> float:
    """Get distance relative tolerance based on date era."""
    return HISTORICAL_DISTANCE_REL if is_deep_historical(year) else MOSHIER_DISTANCE_REL


def get_speed_lon_tolerance(year: int, planet_id: int) -> float:
    """Get longitude speed tolerance based on date era and planet."""
    if is_deep_historical(year):
        return HISTORICAL_SPEED_LON
    return MOSHIER_SPEED_LON_MOON if planet_id == SE_MOON else MOSHIER_SPEED_LON


def get_speed_lat_tolerance(year: int) -> float:
    """Get latitude speed tolerance based on date era."""
    return HISTORICAL_SPEED_LAT if is_deep_historical(year) else MOSHIER_SPEED_LAT


def get_speed_dist_tolerance(year: int) -> float:
    """Get distance speed tolerance based on date era."""
    return HISTORICAL_SPEED_DIST if is_deep_historical(year) else MOSHIER_SPEED_DIST


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMoshierCalcTT:
    """Compare Moshier swe_calc (TT path) between pyswisseph and libephemeris.

    Validates that the TT input path in libephemeris produces the same
    positions as pyswisseph's C implementation for 5 planets x 5 dates = 25
    test cases, comparing all 6 components (lon, lat, dist, speed_lon,
    speed_lat, speed_dist).

    This exercises the inverse Delta T conversion (TT -> UT estimate) in
    _calc_body_moshier(is_ut=False) at planets.py:849-852.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TT_DATES)
    def test_calc_tt_all_components(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier swe_calc (TT) all 6 components match pyswisseph.

        Compares longitude, latitude, distance, and all three velocity
        components between pyswisseph swe.calc() and libephemeris swe_calc()
        when given the same JD in Terrestrial Time.
        """
        jd_tt = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc(jd_tt, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc(jd_tt, planet_id, flag_py)

        lon_tol = get_lon_tolerance(year)
        lat_tol = get_lat_tolerance(year)
        dist_tol = get_dist_tolerance(year)
        speed_lon_tol = get_speed_lon_tolerance(year, planet_id)
        speed_lat_tol = get_speed_lat_tolerance(year)
        speed_dist_tol = get_speed_dist_tolerance(year)

        # Component 0: Longitude
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < lon_tol, (
            f"{planet_name} Moshier TT at {date_desc}: longitude diff {diff_lon:.6f}° "
            f"exceeds tolerance {lon_tol}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )

        # Component 1: Latitude
        diff_lat = abs(pos_swe[1] - pos_py[1])
        assert diff_lat < lat_tol, (
            f"{planet_name} Moshier TT at {date_desc}: latitude diff {diff_lat:.6f}° "
            f"exceeds tolerance {lat_tol}° "
            f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
        )

        # Component 2: Distance (relative comparison)
        if pos_swe[2] > 0:
            rel_dist = abs(pos_swe[2] - pos_py[2]) / pos_swe[2]
            assert rel_dist < dist_tol, (
                f"{planet_name} Moshier TT at {date_desc}: distance relative diff "
                f"{rel_dist:.6f} exceeds tolerance {dist_tol} "
                f"(swe={pos_swe[2]:.8f}, lib={pos_py[2]:.8f})"
            )

        # Component 3: Longitude velocity
        diff_speed_lon = abs(pos_swe[3] - pos_py[3])
        assert diff_speed_lon < speed_lon_tol, (
            f"{planet_name} Moshier TT at {date_desc}: lon velocity diff "
            f"{diff_speed_lon:.6f}°/day exceeds tolerance {speed_lon_tol}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

        # Component 4: Latitude velocity
        diff_speed_lat = abs(pos_swe[4] - pos_py[4])
        assert diff_speed_lat < speed_lat_tol, (
            f"{planet_name} Moshier TT at {date_desc}: lat velocity diff "
            f"{diff_speed_lat:.6f}°/day exceeds tolerance {speed_lat_tol}°/day "
            f"(swe={pos_swe[4]:.6f}, lib={pos_py[4]:.6f})"
        )

        # Component 5: Distance velocity
        diff_speed_dist = abs(pos_swe[5] - pos_py[5])
        assert diff_speed_dist < speed_dist_tol, (
            f"{planet_name} Moshier TT at {date_desc}: dist velocity diff "
            f"{diff_speed_dist:.8f} AU/day exceeds tolerance {speed_dist_tol} "
            f"(swe={pos_swe[5]:.8f}, lib={pos_py[5]:.8f})"
        )


class TestMoshierTTConsistency:
    """Test swe_calc_ut(jd_ut) ~= swe_calc(jd_ut + deltat) for both libraries.

    Validates that the TT and UT1 paths produce consistent results when
    properly related by Delta T. For 5 planets x 2 dates = 10 test cases.

    For each library independently:
      1. Compute jd_ut from calendar date
      2. Compute delta_t = swe_deltat(jd_ut)
      3. Compute jd_tt = jd_ut + delta_t
      4. Verify: swe_calc_ut(jd_ut) ~= swe_calc(jd_tt)

    If the TT -> UT inverse conversion introduces systematic errors,
    these tests will detect the discrepancy.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", CONSISTENCY_DATES)
    def test_tt_ut_consistency(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Both libraries: calc_ut(jd_ut) should match calc(jd_ut + deltat).

        Verifies internal TT/UT consistency for pyswisseph and libephemeris
        independently. Both implementations should produce near-identical
        results from their UT and TT entry points for the same physical instant.
        """
        jd_ut = swe.julday(year, month, day, hour)

        # === pyswisseph consistency ===
        dt_swe = swe.deltat(jd_ut)
        jd_tt_swe = jd_ut + dt_swe

        pos_ut_swe, _ = swe.calc_ut(jd_ut, planet_id, swe.FLG_MOSEPH | swe.FLG_SPEED)
        pos_tt_swe, _ = swe.calc(jd_tt_swe, planet_id, swe.FLG_MOSEPH | swe.FLG_SPEED)

        diff_lon_swe = angular_diff(pos_ut_swe[0], pos_tt_swe[0])
        diff_lat_swe = abs(pos_ut_swe[1] - pos_tt_swe[1])

        assert diff_lon_swe < CONSISTENCY_LONGITUDE, (
            f"pyswisseph {planet_name} at {date_desc}: calc_ut(UT) vs calc(TT) "
            f"longitude diff {diff_lon_swe:.8f}° exceeds {CONSISTENCY_LONGITUDE}° "
            f"(UT={pos_ut_swe[0]:.6f}°, TT={pos_tt_swe[0]:.6f}°, "
            f"dt={dt_swe * 86400:.1f}s)"
        )
        assert diff_lat_swe < CONSISTENCY_LATITUDE, (
            f"pyswisseph {planet_name} at {date_desc}: calc_ut(UT) vs calc(TT) "
            f"latitude diff {diff_lat_swe:.8f}° exceeds {CONSISTENCY_LATITUDE}°"
        )

        if pos_ut_swe[2] > 0:
            rel_dist_swe = abs(pos_ut_swe[2] - pos_tt_swe[2]) / pos_ut_swe[2]
            assert rel_dist_swe < CONSISTENCY_DISTANCE_REL, (
                f"pyswisseph {planet_name} at {date_desc}: calc_ut(UT) vs calc(TT) "
                f"distance rel diff {rel_dist_swe:.8f} exceeds "
                f"{CONSISTENCY_DISTANCE_REL}"
            )

        # === libephemeris consistency ===
        dt_py = ephem.swe_deltat(jd_ut)
        jd_tt_py = jd_ut + dt_py

        pos_ut_py, _ = ephem.swe_calc_ut(jd_ut, planet_id, SEFLG_MOSEPH | SEFLG_SPEED)
        pos_tt_py, _ = ephem.swe_calc(jd_tt_py, planet_id, SEFLG_MOSEPH | SEFLG_SPEED)

        diff_lon_py = angular_diff(pos_ut_py[0], pos_tt_py[0])
        diff_lat_py = abs(pos_ut_py[1] - pos_tt_py[1])

        assert diff_lon_py < CONSISTENCY_LONGITUDE, (
            f"libephemeris {planet_name} at {date_desc}: swe_calc_ut(UT) vs "
            f"swe_calc(TT) longitude diff {diff_lon_py:.8f}° exceeds "
            f"{CONSISTENCY_LONGITUDE}° "
            f"(UT={pos_ut_py[0]:.6f}°, TT={pos_tt_py[0]:.6f}°, "
            f"dt={dt_py * 86400:.1f}s)"
        )
        assert diff_lat_py < CONSISTENCY_LATITUDE, (
            f"libephemeris {planet_name} at {date_desc}: swe_calc_ut(UT) vs "
            f"swe_calc(TT) latitude diff {diff_lat_py:.8f}° exceeds "
            f"{CONSISTENCY_LATITUDE}°"
        )

        if pos_ut_py[2] > 0:
            rel_dist_py = abs(pos_ut_py[2] - pos_tt_py[2]) / pos_ut_py[2]
            assert rel_dist_py < CONSISTENCY_DISTANCE_REL, (
                f"libephemeris {planet_name} at {date_desc}: swe_calc_ut(UT) vs "
                f"swe_calc(TT) distance rel diff {rel_dist_py:.8f} exceeds "
                f"{CONSISTENCY_DISTANCE_REL}"
            )
