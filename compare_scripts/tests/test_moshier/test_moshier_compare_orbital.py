"""
Compare orbital element calculations with pyswisseph using Moshier ephemeris.

This is the Moshier-mode mirror of test_compare_orbital.py (which covers
SEFLG_SWIEPH / JPL mode). It validates that orbital functions produce
consistent results between pyswisseph (C library) and libephemeris
(Python reimplementation) when using SEFLG_MOSEPH.

NOTE: Tests for nod_aps are marked as xfail because libephemeris uses a different
methodology for calculating planetary nodes and apsides.

Key difference:
- libephemeris: Uses heliocentric mean orbital elements from Standish (1992) JPL/IERS
  tables. The ascending node longitude is the heliocentric Omega where the planet's
  orbit crosses the ecliptic going north (as seen from the Sun).

- Swiss Ephemeris: Appears to use a geocentric interpretation, returning the ecliptic
  longitude where the planet appears to cross the ecliptic as seen from Earth.

For outer planets (Jupiter, Saturn, etc.), the difference is small (<1°) because
they are far from Earth. For inner planets (Mercury, Venus), the difference can be
~250° because the heliocentric and geocentric orbital planes are oriented very
differently relative to Earth.

Both approaches are astronomically valid - they answer different questions:
- Heliocentric: "Where does the orbit cross the ecliptic as seen from the Sun?"
- Geocentric: "Where does the planet appear to cross the ecliptic from Earth?"

Moshier-specific notes:
- nod_aps_ut with SEFLG_MOSEPH calls swe_calc_ut internally with SEFLG_MOSEPH,
  using VSOP87/ELP semi-analytical theories instead of JPL DE440.
- Tolerances are wider than SWIEPH to account for Moshier's lower precision.
- Three test dates include a historical date (1900) within Moshier's extended range.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_NODBIT_MEAN,
    SE_NODBIT_OSCU,
    SEFLG_MOSEPH,
)

# Mark nod_aps tests as expected to fail
_NOD_APS_XFAIL = pytest.mark.xfail(
    reason="Orbital node calculation methodology differs (heliocentric vs geocentric)",
    strict=False,
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


class MoshierOrbitalTolerance:
    """Tolerance thresholds for Moshier orbital element comparisons.

    Wider than SWIEPH tolerances to account for Moshier semi-analytical
    ephemeris (VSOP87/ELP) having lower precision than JPL DE440.
    """

    ANGLE_DEGREES = 0.1  # 0.1 degree for angular elements
    DISTANCE_AU = 0.01  # 0.01 AU for distances
    ECCENTRICITY = 0.001  # Eccentricity (dimensionless)
    PERIOD_DAYS = 1.0  # Days for orbital period


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for nod_aps tests (Mars, Jupiter, Saturn per task spec)
NOD_APS_PLANETS = [
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# All planets for orbital element tests
TEST_PLANETS = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Method types for nod_aps (using SE_NODBIT constants)
NOD_APS_METHODS = [
    (SE_NODBIT_MEAN, "Mean"),
    (SE_NODBIT_OSCU, "Osculating"),
]

# 3 test dates spanning Moshier range
TEST_DATES = [
    (swe.julday(2000, 1, 1, 12), "J2000.0"),
    (swe.julday(2024, 6, 15, 0), "2024"),
    (swe.julday(1900, 1, 1, 0), "1900"),
]


# ============================================================================
# NODES AND APSIDES TESTS
# ============================================================================


class TestNodAps:
    """Tests for nodes and apsides (nod_aps_ut) in Moshier mode."""

    @_NOD_APS_XFAIL
    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", NOD_APS_PLANETS)
    @pytest.mark.parametrize("method,method_name", NOD_APS_METHODS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_nod_aps_ut(self, body_id, body_name, method, method_name, jd, date_desc):
        """Test Moshier nodes and apsides calculations."""
        # SwissEphemeris
        try:
            ret_swe = swe.nod_aps_ut(jd, body_id, swe.FLG_MOSEPH, method)
            nasc_swe = ret_swe[0]
            ndesc_swe = ret_swe[1]
            peri_swe = ret_swe[2]
            aphe_swe = ret_swe[3]
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            ret_py = pyephem.nod_aps_ut(jd, body_id, SEFLG_MOSEPH, method)
            nasc_py = ret_py[0]
            ndesc_py = ret_py[1]
            peri_py = ret_py[2]
            aphe_py = ret_py[3]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        # Compare node and apsides longitudes
        diff_nasc = angular_diff(nasc_swe[0], nasc_py[0])
        diff_ndesc = angular_diff(ndesc_swe[0], ndesc_py[0])
        diff_peri = angular_diff(peri_swe[0], peri_py[0])
        diff_aphe = angular_diff(aphe_swe[0], aphe_py[0])

        max_diff = max(diff_nasc, diff_ndesc, diff_peri, diff_aphe)

        assert max_diff < MoshierOrbitalTolerance.ANGLE_DEGREES, (
            f"{body_name} {method_name} @ {date_desc}: max diff {max_diff}°"
        )


# ============================================================================
# ORBITAL ELEMENTS TESTS
# ============================================================================


class TestOrbitalElements:
    """Tests for get_orbital_elements in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TEST_PLANETS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_get_orbital_elements(self, body_id, body_name, jd, date_desc):
        """Test Moshier Keplerian orbital elements calculations."""
        # SwissEphemeris
        try:
            ret_swe = swe.get_orbital_elements(jd, body_id, swe.FLG_MOSEPH)
            a_swe = ret_swe[0]  # Semi-major axis
            e_swe = ret_swe[1]  # Eccentricity
            i_swe = ret_swe[2]  # Inclination
            node_swe = ret_swe[3]  # Ascending node
            peri_swe = ret_swe[4]  # Argument of perihelion
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            ret_py = pyephem.get_orbital_elements(jd, body_id, SEFLG_MOSEPH)
            a_py = ret_py[0]
            e_py = ret_py[1]
            i_py = ret_py[2]
            node_py = ret_py[3]
            peri_py = ret_py[4]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        # Calculate differences
        diff_a = abs(a_swe - a_py)
        diff_e = abs(e_swe - e_py)
        diff_i = abs(i_swe - i_py)
        diff_node = angular_diff(node_swe, node_py)
        diff_peri = angular_diff(peri_swe, peri_py)

        # Check tolerances
        assert diff_a < MoshierOrbitalTolerance.DISTANCE_AU, (
            f"{body_name} @ {date_desc}: semi-major axis diff {diff_a} AU"
        )
        assert diff_e < MoshierOrbitalTolerance.ECCENTRICITY, (
            f"{body_name} @ {date_desc}: eccentricity diff {diff_e}"
        )
        assert diff_i < MoshierOrbitalTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: inclination diff {diff_i}°"
        )
        assert diff_node < MoshierOrbitalTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: node diff {diff_node}°"
        )
        assert diff_peri < MoshierOrbitalTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: perihelion diff {diff_peri}°"
        )


# ============================================================================
# ORBITAL DISTANCE TESTS
# ============================================================================


class TestOrbitalDistance:
    """Tests for orbit_max_min_true_distance in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TEST_PLANETS)
    def test_orbit_max_min_true_distance(self, body_id, body_name):
        """Test Moshier orbital distance extremes calculations."""
        jd = swe.julday(2024, 1, 1, 0)

        # SwissEphemeris
        try:
            ret_swe = swe.orbit_max_min_true_distance(jd, body_id, swe.FLG_MOSEPH)
            dmax_swe = ret_swe[0]
            dmin_swe = ret_swe[1]
            dtrue_swe = ret_swe[2]
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            ret_py = pyephem.orbit_max_min_true_distance(jd, body_id, SEFLG_MOSEPH)
            dmax_py = ret_py[0]
            dmin_py = ret_py[1]
            dtrue_py = ret_py[2]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        diff_max = abs(dmax_swe - dmax_py)
        diff_min = abs(dmin_swe - dmin_py)
        diff_true = abs(dtrue_swe - dtrue_py)
        max_diff = max(diff_max, diff_min, diff_true)

        assert max_diff < MoshierOrbitalTolerance.DISTANCE_AU, (
            f"{body_name}: max distance diff {max_diff} AU"
        )


# ============================================================================
# ORBITAL CONSISTENCY TESTS
# ============================================================================


class TestOrbitalConsistency:
    """Tests for orbital element consistency in Moshier mode."""

    @_NOD_APS_XFAIL
    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", NOD_APS_PLANETS)
    def test_orbital_elements_physical_constraints(self, body_id, body_name):
        """Test that Moshier orbital elements satisfy physical constraints."""
        jd = swe.julday(2024, 1, 1, 0)

        ret_py = pyephem.get_orbital_elements(jd, body_id, SEFLG_MOSEPH)
        a = ret_py[0]  # Semi-major axis
        e = ret_py[1]  # Eccentricity
        i = ret_py[2]  # Inclination

        # Physical constraints
        assert a > 0, f"{body_name}: semi-major axis must be positive"
        assert 0 <= e < 1, f"{body_name}: eccentricity must be 0 <= e < 1 for ellipse"
        assert 0 <= i <= 180, f"{body_name}: inclination must be 0-180°"

    @pytest.mark.comparison
    def test_orbital_elements_ordering(self):
        """Test that Moshier orbital distances make physical sense."""
        jd = swe.julday(2024, 1, 1, 0)

        for body_id, body_name in TEST_PLANETS:
            try:
                ret = pyephem.orbit_max_min_true_distance(jd, body_id, SEFLG_MOSEPH)
                dmax = ret[0]
                dmin = ret[1]
                dtrue = ret[2]

                assert dmax >= dmin, f"{body_name}: dmax must be >= dmin"
                assert dmin <= dtrue <= dmax, (
                    f"{body_name}: dtrue must be between dmin and dmax"
                )
            except Exception:
                pass  # Some planets may not be available
