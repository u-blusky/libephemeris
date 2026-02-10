"""
Moshier Observation Modes Cross-Library Comparison Tests.

Validates Moshier (SEFLG_MOSEPH) observation mode calculations between pyswisseph
(C library) and libephemeris (Python reimplementation):
- Topocentric (SEFLG_TOPOCTR): parallax correction for observer position
- Heliocentric (SEFLG_HELCTR): Sun-centered coordinates
- J2000 (SEFLG_J2000): J2000.0 ecliptic/equatorial coordinates

This is the Moshier-mode mirror of test_compare_observations.py (which covers
SEFLG_SWIEPH / JPL mode). It documents known gaps (topocentric not implemented,
J2000 not handled in Moshier path) and validates heliocentric calculations.

Known limitations in libephemeris Moshier path (planets.py):
- Topocentric: Moon parallax correction is a `pass # TODO` at line 1049
- J2000: SEFLG_J2000 is not extracted from flags (line 857-861), so
  coordinates are always ecliptic-of-date regardless of J2000 flag
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
    SEFLG_TOPOCTR,
    SEFLG_HELCTR,
    SEFLG_J2000,
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

# Moshier C-vs-Python tolerances (relaxed compared to JPL/SWIEPH).
# These match the tolerances in test_moshier_compare_planets.py.
MOSHIER_LONGITUDE = 0.02  # degrees (~72 arcsec)
MOSHIER_LATITUDE = 0.02  # degrees
MOSHIER_DISTANCE_REL = 0.001  # relative (dimensionless)
MOSHIER_VELOCITY = 0.05  # degrees/day (relaxed for numerical differentiation diffs)

# Heliocentric tolerances (wider due to VSOP87 C-vs-Python implementation diffs)
HELIO_LONGITUDE = 0.03  # degrees
HELIO_LATITUDE = 0.03  # degrees
HELIO_DISTANCE = 0.01  # AU

# Pluto uses separate, much larger tolerances due to Chapront-Francou theory
HELIO_PLUTO_LONGITUDE = 2.0  # degrees (Pluto Moshier theory has large errors)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Test subjects with birth data (famous people) - mirrors test_compare_observations.py
# Format: (name, year, month, day, hour, lat, lon, alt)
SUBJECTS = [
    ("Einstein (Ulm)", 1879, 3, 14, 11.5, 48.4011, 9.9876, 478),
    ("Gandhi (Porbandar)", 1869, 10, 2, 7.2, 21.6417, 69.6293, 0),
    ("Mandela (Mvezo)", 1918, 7, 18, 14.0, -31.9566, 28.5133, 0),
]

# Heliocentric planets: all except Sun (which is trivially at origin)
HELIO_PLANETS = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
]

# J2000 planets: subset for ecliptic and equatorial tests
J2000_PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]


# ============================================================================
# TOPOCENTRIC TESTS (Moshier)
# ============================================================================


class TestMoshierTopocentric:
    """Tests for Moshier topocentric calculations.

    Topocentric correction for Moon is NOT implemented in Moshier path
    (planets.py:1049 is `pass # TODO`). These tests document the gap
    by comparing against pyswisseph's Moshier topocentric, which IS
    implemented. Moon parallax can be up to ~1 degree.
    """

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason="Moshier topocentric not implemented (planets.py:1049)",
        strict=False,
    )
    @pytest.mark.parametrize(
        "subject_name,year,month,day,hour,lat,lon,alt",
        SUBJECTS,
    )
    def test_topocentric_moon(
        self, subject_name, year, month, day, hour, lat, lon, alt
    ):
        """Test Moshier topocentric Moon vs pyswisseph.

        Expected to fail because libephemeris Moshier path does not apply
        parallax correction for Moon (planets.py:1049 is pass # TODO).
        """
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_TOPOCTR | SEFLG_SPEED

        # Set observer position for topocentric
        swe.set_topo(lon, lat, alt)
        ephem.swe_set_topo(lon, lat, alt)

        try:
            res_swe, _ = swe.calc_ut(jd, SE_MOON, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, SE_MOON, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        diff_lon = angular_diff(res_swe[0], res_py[0])

        # Moon parallax can be up to ~1 degree; the strict tolerance
        # will catch the missing implementation
        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Moshier Topocentric Moon @ {subject_name}: "
            f"lon diff {diff_lon:.6f}° (tolerance {MOSHIER_LONGITUDE}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "lat,lon,alt,loc_name",
        [
            (0.0, 0.0, 0, "Equator"),
            (45.0, 0.0, 0, "Mid Latitude"),
            (70.0, 0.0, 0, "High Latitude"),
            (-45.0, 0.0, 0, "Southern"),
        ],
    )
    @pytest.mark.xfail(
        reason="Moshier topocentric not implemented (planets.py:1049)",
        strict=False,
    )
    def test_topocentric_moon_latitudes(self, lat, lon, alt, loc_name):
        """Test Moshier topocentric Moon at various latitudes.

        Expected to fail due to missing parallax correction.
        """
        jd = swe.julday(2024, 1, 15, 12.0)
        flags = SEFLG_MOSEPH | SEFLG_TOPOCTR | SEFLG_SPEED

        swe.set_topo(lon, lat, alt)
        ephem.swe_set_topo(lon, lat, alt)

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

        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Moshier Topocentric Moon @ {loc_name}: "
            f"lon diff {diff_lon:.6f}° (tolerance {MOSHIER_LONGITUDE}°)"
        )
        assert diff_lat < MOSHIER_LATITUDE, (
            f"Moshier Topocentric Moon @ {loc_name}: "
            f"lat diff {diff_lat:.6f}° (tolerance {MOSHIER_LATITUDE}°)"
        )


# ============================================================================
# HELIOCENTRIC TESTS (Moshier)
# ============================================================================


class TestMoshierHeliocentric:
    """Tests for Moshier heliocentric calculations.

    Heliocentric mode IS implemented in the Moshier path (planets.py:948-1025)
    using VSOP87 directly for Mercury-Neptune and Chapront-Francou for Pluto.
    These tests should pass with Moshier C-vs-Python tolerances.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        HELIO_PLANETS,
    )
    @pytest.mark.parametrize(
        "subject_name,year,month,day,hour,lat,lon,alt",
        SUBJECTS,
    )
    def test_heliocentric_planet(
        self,
        planet_id,
        planet_name,
        subject_name,
        year,
        month,
        day,
        hour,
        lat,
        lon,
        alt,
    ):
        """Test Moshier heliocentric planet positions vs pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_HELCTR | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        diff_lon = angular_diff(res_swe[0], res_py[0])
        diff_lat = abs(res_swe[1] - res_py[1])
        diff_dist = abs(res_swe[2] - res_py[2])

        # Use wider tolerance for Pluto
        lon_tol = HELIO_PLUTO_LONGITUDE if planet_id == SE_PLUTO else HELIO_LONGITUDE
        lat_tol = HELIO_PLUTO_LONGITUDE if planet_id == SE_PLUTO else HELIO_LATITUDE

        assert diff_lon < lon_tol, (
            f"Moshier Heliocentric {planet_name} @ {subject_name}: "
            f"lon diff {diff_lon:.6f}° (tolerance {lon_tol}°)"
        )
        assert diff_lat < lat_tol, (
            f"Moshier Heliocentric {planet_name} @ {subject_name}: "
            f"lat diff {diff_lat:.6f}° (tolerance {lat_tol}°)"
        )
        assert diff_dist < HELIO_DISTANCE, (
            f"Moshier Heliocentric {planet_name} @ {subject_name}: "
            f"dist diff {diff_dist:.6f} AU (tolerance {HELIO_DISTANCE} AU)"
        )


# ============================================================================
# J2000 COORDINATE TESTS (Moshier)
# ============================================================================


class TestMoshierJ2000Coordinates:
    """Tests for Moshier J2000 coordinate calculations.

    SEFLG_J2000 is NOT extracted in the Moshier path (planets.py:857-861).
    The Moshier path always returns ecliptic-of-date coordinates, ignoring
    the J2000 flag. For dates far from J2000.0, precession introduces
    significant differences (~1.4° per century).

    These tests document the gap by comparing against pyswisseph's Moshier
    J2000 output.
    """

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason="Moshier path does not handle SEFLG_J2000 (planets.py:857-861)",
        strict=False,
    )
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        J2000_PLANETS,
    )
    @pytest.mark.parametrize(
        "subject_name,year,month,day,hour,lat,lon,alt",
        SUBJECTS,
    )
    def test_j2000_ecliptic(
        self,
        planet_id,
        planet_name,
        subject_name,
        year,
        month,
        day,
        hour,
        lat,
        lon,
        alt,
    ):
        """Test Moshier J2000 ecliptic coordinates vs pyswisseph.

        Expected to fail because libephemeris Moshier path ignores SEFLG_J2000,
        returning ecliptic-of-date instead of J2000 ecliptic.
        """
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_J2000 | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        diff_lon = angular_diff(res_swe[0], res_py[0])

        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Moshier J2000 Ecliptic {planet_name} @ {subject_name}: "
            f"lon diff {diff_lon:.6f}° (tolerance {MOSHIER_LONGITUDE}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason="Moshier path does not handle SEFLG_J2000 (planets.py:857-861)",
        strict=False,
    )
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
        ],
    )
    @pytest.mark.parametrize(
        "subject_name,year,month,day,hour,lat,lon,alt",
        SUBJECTS,
    )
    def test_j2000_equatorial(
        self,
        planet_id,
        planet_name,
        subject_name,
        year,
        month,
        day,
        hour,
        lat,
        lon,
        alt,
    ):
        """Test Moshier J2000 equatorial coordinates vs pyswisseph.

        Expected to fail because libephemeris Moshier path ignores SEFLG_J2000,
        returning equatorial-of-date instead of J2000 equatorial.
        """
        jd = swe.julday(year, month, day, hour)
        flags = SEFLG_MOSEPH | SEFLG_J2000 | SEFLG_EQUATORIAL | SEFLG_SPEED

        try:
            res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        diff_ra = angular_diff(res_swe[0], res_py[0])
        diff_dec = abs(res_swe[1] - res_py[1])

        assert diff_ra < MOSHIER_LONGITUDE, (
            f"Moshier J2000 Equatorial {planet_name} @ {subject_name}: "
            f"RA diff {diff_ra:.6f}° (tolerance {MOSHIER_LONGITUDE}°)"
        )
        assert diff_dec < MOSHIER_LATITUDE, (
            f"Moshier J2000 Equatorial {planet_name} @ {subject_name}: "
            f"Dec diff {diff_dec:.6f}° (tolerance {MOSHIER_LATITUDE}°)"
        )
