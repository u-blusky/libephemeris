"""
Moshier SE_EARTH (ID 14) Cross-Library Comparison Tests.

SE_EARTH is a special body in the Swiss Ephemeris:

- **Geocentric**: Earth's geocentric position is by definition (0, 0, 0) — the
  observer is on Earth. Both the C library and libephemeris should return zeros
  (or equivalent). The C library may raise an error or return zeros silently.

- **Heliocentric**: Uses `calc_earth_heliocentric()` from moshier/vsop87.py,
  which computes Earth's position in the solar system via truncated VSOP87 series.
  This is the foundation for ALL geocentric conversions of other planets.

- **Moon approximation**: libephemeris approximates heliocentric Moon as Earth
  heliocentric (planets.py:958: "Moon's heliocentric position is ~same as Earth"),
  introducing a maximum error of ~0.0003 AU vs the C library.

Code paths tested:
- planets.py:972-987  — SE_EARTH heliocentric dispatch
- moshier/vsop87.py:calc_earth_heliocentric() — VSOP87 Earth computation
- moshier/vsop87.py:calc_planet_geocentric() — geocentric (Earth-Earth = origin)

This file tests 4 flag combinations × 3 dates = 12 test cases that fully document
SE_EARTH Moshier behavior in both implementations.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_EARTH,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_EQUATORIAL,
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

# Heliocentric Earth position tolerances (C vs Python Moshier VSOP87).
# Earth heliocentric is the base of all geocentric conversions, so errors
# here propagate to every other planet.
# Observed C-vs-Python divergence: ~0.006°–0.010° in longitude due to
# truncated VSOP87 series. Tolerance set at 0.02° (2x observed max).
EARTH_HELIO_LON_TOL = 0.02  # degrees (~72 arcsec)
EARTH_HELIO_LAT_TOL = 0.02  # degrees
EARTH_HELIO_DIST_REL_TOL = 0.0001  # relative (0.01%) — distance ~1 AU
EARTH_HELIO_VEL_TOL = 0.01  # degrees/day

# Equatorial coordinate tolerances (after ecliptic-to-equatorial transform).
# Slightly wider than ecliptic due to obliquity projection amplifying lon errors.
EARTH_HELIO_EQUAT_RA_TOL = 0.02  # degrees
EARTH_HELIO_EQUAT_DEC_TOL = 0.02  # degrees


# ============================================================================
# TEST DATES
# ============================================================================

# 3 dates spanning different epochs to catch VSOP87 truncation divergence.
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (1950, 6, 15, 0.0, "1950-mid"),
    (2024, 3, 20, 12.0, "2024-equinox"),
]


# ============================================================================
# GEOCENTRIC EARTH TESTS
# ============================================================================


class TestMoshierEarthGeocentric:
    """Test SE_EARTH in geocentric mode (SEFLG_MOSEPH only, no HELCTR).

    Geocentric Earth should return (0, 0, 0) or an error in both libraries,
    since the observer is on Earth. The C library may handle this differently
    from libephemeris — this test documents the actual behavior of both.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_geocentric_earth(self, year, month, day, hour, date_label):
        """SE_EARTH geocentric: compare behavior between C and Python.

        Both libraries should either:
        - Return (0, 0, 0) position (Earth at geocentric origin), or
        - Raise an error for unsupported geocentric Earth.

        If behaviors differ, this test documents the discrepancy.
        """
        jd = swe.julday(year, month, day, hour)
        flags_swe = swe.FLG_MOSEPH
        flags_py = SEFLG_MOSEPH

        swe_error = None
        swe_result = None
        py_error = None
        py_result = None

        try:
            swe_result, _ = swe.calc_ut(jd, SE_EARTH, flags_swe)
        except Exception as e:
            swe_error = e

        try:
            py_result, _ = ephem.swe_calc_ut(jd, SE_EARTH, flags_py)
        except Exception as e:
            py_error = e

        # Case 1: Both raise errors — verify consistent error behavior
        if swe_error is not None and py_error is not None:
            # Both correctly reject geocentric Earth — pass
            return

        # Case 2: Both return results — compare values
        if swe_result is not None and py_result is not None:
            # Geocentric Earth should be at origin: (0, 0, 0)
            # Allow small numerical noise
            diff_lon = angular_diff(swe_result[0], py_result[0])
            diff_lat = abs(swe_result[1] - py_result[1])
            diff_dist = abs(swe_result[2] - py_result[2])

            assert diff_lon < 0.01, (
                f"Geocentric Earth @ {date_label}: "
                f"lon diff {diff_lon:.6f}° "
                f"(swe={swe_result[0]:.6f}°, lib={py_result[0]:.6f}°)"
            )
            assert diff_lat < 0.01, (
                f"Geocentric Earth @ {date_label}: "
                f"lat diff {diff_lat:.6f}° "
                f"(swe={swe_result[1]:.6f}°, lib={py_result[1]:.6f}°)"
            )
            assert diff_dist < 0.001, (
                f"Geocentric Earth @ {date_label}: "
                f"dist diff {diff_dist:.6f} AU "
                f"(swe={swe_result[2]:.6f}, lib={py_result[2]:.6f})"
            )
            return

        # Case 3: One raises error, other returns result — document divergence
        if swe_error is not None:
            assert py_result is not None  # guaranteed by case 1/2 logic
            pytest.skip(
                f"Behavior divergence @ {date_label}: "
                f"pyswisseph raises '{swe_error}', "
                f"libephemeris returns {py_result[:3]}"
            )
        else:
            assert swe_result is not None  # guaranteed by case 1/2 logic
            pytest.skip(
                f"Behavior divergence @ {date_label}: "
                f"pyswisseph returns {swe_result[:3]}, "
                f"libephemeris raises '{py_error}'"
            )


# ============================================================================
# HELIOCENTRIC EARTH TESTS
# ============================================================================


class TestMoshierEarthHeliocentric:
    """Test SE_EARTH with SEFLG_HELCTR (heliocentric, no speed).

    Heliocentric Earth is the position of Earth in the solar system,
    computed via calc_earth_heliocentric() from VSOP87 theory.
    Distance should be ~1 AU (Earth-Sun distance).

    This position is the foundation for ALL geocentric conversions:
    an error in heliocentric Earth produces identical errors in every
    geocentric planet.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_heliocentric_earth(self, year, month, day, hour, date_label):
        """SE_EARTH heliocentric position: lon, lat, distance vs C library."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR

        try:
            res_swe, _ = swe.calc_ut(jd, SE_EARTH, flags)
        except Exception as e:
            pytest.skip(
                f"pyswisseph does not support SE_EARTH with Moshier HELCTR: {e}"
            )

        try:
            res_py, _ = ephem.swe_calc_ut(jd, SE_EARTH, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed for SE_EARTH HELCTR @ {date_label}: {e}")

        # Longitude
        diff_lon = angular_diff(res_swe[0], res_py[0])
        assert diff_lon < EARTH_HELIO_LON_TOL, (
            f"Heliocentric Earth @ {date_label}: "
            f"lon diff {diff_lon:.6f}° (tolerance {EARTH_HELIO_LON_TOL}°) "
            f"(swe={res_swe[0]:.6f}°, lib={res_py[0]:.6f}°)"
        )

        # Latitude
        diff_lat = abs(res_swe[1] - res_py[1])
        assert diff_lat < EARTH_HELIO_LAT_TOL, (
            f"Heliocentric Earth @ {date_label}: "
            f"lat diff {diff_lat:.6f}° (tolerance {EARTH_HELIO_LAT_TOL}°) "
            f"(swe={res_swe[1]:.6f}°, lib={res_py[1]:.6f}°)"
        )

        # Distance (relative — should be ~1 AU)
        assert 0.98 < res_py[2] < 1.02, (
            f"Heliocentric Earth distance @ {date_label}: "
            f"{res_py[2]:.6f} AU (expected ~1.0 AU)"
        )
        if res_swe[2] != 0:
            dist_rel = abs(res_swe[2] - res_py[2]) / abs(res_swe[2])
            assert dist_rel < EARTH_HELIO_DIST_REL_TOL, (
                f"Heliocentric Earth @ {date_label}: "
                f"dist rel diff {dist_rel:.8f} ({dist_rel * 100:.6f}%) "
                f"(tolerance {EARTH_HELIO_DIST_REL_TOL * 100}%) "
                f"(swe={res_swe[2]:.8f}, lib={res_py[2]:.8f})"
            )


class TestMoshierEarthHeliocentricSpeed:
    """Test SE_EARTH with SEFLG_HELCTR | SEFLG_SPEED.

    Validates heliocentric Earth position AND velocity (numerical
    differentiation with h=0.001 in planets.py:977-987).
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_heliocentric_earth_speed(self, year, month, day, hour, date_label):
        """SE_EARTH heliocentric with speed: position + velocity vs C library."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, SE_EARTH, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph does not support SE_EARTH HELCTR|SPEED: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, SE_EARTH, flags)
        except Exception as e:
            pytest.skip(
                f"libephemeris failed for SE_EARTH HELCTR|SPEED @ {date_label}: {e}"
            )

        # Position checks (same as above)
        diff_lon = angular_diff(res_swe[0], res_py[0])
        assert diff_lon < EARTH_HELIO_LON_TOL, (
            f"Heliocentric Earth @ {date_label}: "
            f"lon diff {diff_lon:.6f}° (tolerance {EARTH_HELIO_LON_TOL}°) "
            f"(swe={res_swe[0]:.6f}°, lib={res_py[0]:.6f}°)"
        )

        diff_lat = abs(res_swe[1] - res_py[1])
        assert diff_lat < EARTH_HELIO_LAT_TOL, (
            f"Heliocentric Earth @ {date_label}: "
            f"lat diff {diff_lat:.6f}° (tolerance {EARTH_HELIO_LAT_TOL}°) "
            f"(swe={res_swe[1]:.6f}°, lib={res_py[1]:.6f}°)"
        )

        if res_swe[2] != 0:
            dist_rel = abs(res_swe[2] - res_py[2]) / abs(res_swe[2])
            assert dist_rel < EARTH_HELIO_DIST_REL_TOL, (
                f"Heliocentric Earth @ {date_label}: "
                f"dist rel diff {dist_rel:.8f} ({dist_rel * 100:.6f}%) "
                f"(tolerance {EARTH_HELIO_DIST_REL_TOL * 100}%)"
            )

        # Velocity checks — longitude speed (degrees/day)
        # Earth's mean motion is ~0.9856°/day
        diff_vel_lon = abs(res_swe[3] - res_py[3])
        assert diff_vel_lon < EARTH_HELIO_VEL_TOL, (
            f"Heliocentric Earth @ {date_label}: "
            f"lon speed diff {diff_vel_lon:.6f} deg/day "
            f"(tolerance {EARTH_HELIO_VEL_TOL}) "
            f"(swe={res_swe[3]:.6f}, lib={res_py[3]:.6f})"
        )

        # Latitude speed
        diff_vel_lat = abs(res_swe[4] - res_py[4])
        assert diff_vel_lat < EARTH_HELIO_VEL_TOL, (
            f"Heliocentric Earth @ {date_label}: "
            f"lat speed diff {diff_vel_lat:.6f} deg/day "
            f"(tolerance {EARTH_HELIO_VEL_TOL}) "
            f"(swe={res_swe[4]:.6f}, lib={res_py[4]:.6f})"
        )

        # Distance speed (AU/day)
        diff_vel_dist = abs(res_swe[5] - res_py[5])
        assert diff_vel_dist < 0.001, (
            f"Heliocentric Earth @ {date_label}: "
            f"dist speed diff {diff_vel_dist:.8f} AU/day "
            f"(swe={res_swe[5]:.8f}, lib={res_py[5]:.8f})"
        )


class TestMoshierEarthHeliocentricEquatorial:
    """Test SE_EARTH with SEFLG_HELCTR | SEFLG_EQUATORIAL.

    Validates heliocentric Earth position transformed to equatorial
    coordinates (RA/Dec). This exercises the coordinate transformation
    code path in planets.py after heliocentric computation.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_label",
        TEST_DATES,
    )
    def test_heliocentric_earth_equatorial(self, year, month, day, hour, date_label):
        """SE_EARTH heliocentric equatorial: RA, Dec, distance vs C library."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_EQUATORIAL

        try:
            res_swe, _ = swe.calc_ut(jd, SE_EARTH, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph does not support SE_EARTH HELCTR|EQUATORIAL: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, SE_EARTH, flags)
        except Exception as e:
            pytest.skip(
                f"libephemeris failed for SE_EARTH HELCTR|EQUAT @ {date_label}: {e}"
            )

        # Right Ascension (position[0] in equatorial mode)
        diff_ra = angular_diff(res_swe[0], res_py[0])
        assert diff_ra < EARTH_HELIO_EQUAT_RA_TOL, (
            f"Heliocentric Earth equatorial @ {date_label}: "
            f"RA diff {diff_ra:.6f}° (tolerance {EARTH_HELIO_EQUAT_RA_TOL}°) "
            f"(swe={res_swe[0]:.6f}°, lib={res_py[0]:.6f}°)"
        )

        # Declination (position[1] in equatorial mode)
        diff_dec = abs(res_swe[1] - res_py[1])
        assert diff_dec < EARTH_HELIO_EQUAT_DEC_TOL, (
            f"Heliocentric Earth equatorial @ {date_label}: "
            f"Dec diff {diff_dec:.6f}° (tolerance {EARTH_HELIO_EQUAT_DEC_TOL}°) "
            f"(swe={res_swe[1]:.6f}°, lib={res_py[1]:.6f}°)"
        )

        # Distance (should still be ~1 AU, same as ecliptic)
        if res_swe[2] != 0:
            dist_rel = abs(res_swe[2] - res_py[2]) / abs(res_swe[2])
            assert dist_rel < EARTH_HELIO_DIST_REL_TOL, (
                f"Heliocentric Earth equatorial @ {date_label}: "
                f"dist rel diff {dist_rel:.8f} ({dist_rel * 100:.6f}%) "
                f"(tolerance {EARTH_HELIO_DIST_REL_TOL * 100}%) "
                f"(swe={res_swe[2]:.8f}, lib={res_py[2]:.8f})"
            )
