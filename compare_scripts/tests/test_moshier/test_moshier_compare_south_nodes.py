"""
Moshier South Lunar Node (Ketu) Cross-Library Comparison Tests.

Validates the behavior of negative body IDs (-SE_MEAN_NODE = -10,
-SE_TRUE_NODE = -11) with SEFLG_MOSEPH between pyswisseph (C library)
and libephemeris (Python reimplementation).

libephemeris supports negative body IDs in Moshier mode: -SE_MEAN_NODE
and -SE_TRUE_NODE return north_node_lon + 180° with negated latitude
(planets.py:901-912). pyswisseph may not support negative body IDs in
Moshier mode (raising SwissephError) or may implement them differently.

This is critical for Jyotish (Vedic astrology) where Ketu (south node)
is one of the most important points. Applications using -10 for Ketu
(a common pattern in pyswisseph) need to know whether the behavior
differs between implementations.

API compatibility findings are documented via xfail markers and
assertion messages.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TOLERANCES
# ============================================================================

# Reuse tolerances from test_moshier_compare_lunar.py
MEAN_NODE_MOSHIER = 0.02  # degrees (~72 arcsec)
TRUE_NODE_MOSHIER = 0.2  # degrees (~720 arcsec)


# ============================================================================
# TEST DATA - 3 representative dates spanning different eras
# ============================================================================

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (1980, 5, 20, 0.0, "Past"),
    (2024, 11, 5, 18.0, "Recent"),
]

# South node body IDs (negative of north node IDs)
SOUTH_NODE_IDS = [
    (-SE_MEAN_NODE, "Mean South Node (Ketu)", SE_MEAN_NODE, MEAN_NODE_MOSHIER),
    (-SE_TRUE_NODE, "True South Node (Ketu)", SE_TRUE_NODE, TRUE_NODE_MOSHIER),
]


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def _pyswisseph_supports_negative_ids() -> bool:
    """Probe whether pyswisseph supports negative body IDs with SEFLG_MOSEPH.

    Returns True if swe.calc_ut succeeds with body ID -10, False if it raises.
    """
    jd = swe.julday(2000, 1, 1, 12.0)
    try:
        swe.calc_ut(jd, -10, swe.FLG_MOSEPH)
        return True
    except Exception:
        return False


# Evaluate once at module load
_SWE_SUPPORTS_NEGATIVE = _pyswisseph_supports_negative_ids()


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestSouthNodeLibephemeris:
    """Verify libephemeris handles negative body IDs in Moshier mode.

    libephemeris explicitly supports -SE_MEAN_NODE and -SE_TRUE_NODE
    (planets.py:901-912): it computes the north node, then returns
    (lon + 180) % 360 with negated latitude.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    @pytest.mark.parametrize("south_id,south_name,north_id,tol", SOUTH_NODE_IDS)
    def test_south_node_is_180_from_north(
        self, year, month, day, hour, desc, south_id, south_name, north_id, tol
    ):
        """Test that libephemeris south node lon = north node lon + 180°."""
        jd = swe.julday(year, month, day, hour)

        pos_north, _ = ephem.swe_calc_ut(jd, north_id, SEFLG_MOSEPH)
        pos_south, _ = ephem.swe_calc_ut(jd, south_id, SEFLG_MOSEPH)

        expected_south_lon = (pos_north[0] + 180.0) % 360.0
        diff = angular_diff(pos_south[0], expected_south_lon)

        assert diff < 1e-10, (
            f"{south_name} at {desc}: lon {pos_south[0]:.6f}° != "
            f"north+180 {expected_south_lon:.6f}° (diff={diff:.2e}°)"
        )

        # Latitude should be negated
        assert abs(pos_south[1] + pos_north[1]) < 1e-10, (
            f"{south_name} at {desc}: lat {pos_south[1]:.6f}° != "
            f"-north_lat {-pos_north[1]:.6f}°"
        )

        # Distance should be identical
        assert abs(pos_south[2] - pos_north[2]) < 1e-10, (
            f"{south_name} at {desc}: dist differs from north node"
        )


class TestSouthNodeCrossLibrary:
    """Compare south node behavior between pyswisseph and libephemeris.

    pyswisseph may not support negative body IDs in Moshier mode. If it
    raises an exception, tests are marked xfail to document the API
    incompatibility. If both support negative IDs, longitudes are compared
    with the same tolerances as north node comparisons.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    @pytest.mark.parametrize("south_id,south_name,north_id,tol", SOUTH_NODE_IDS)
    def test_south_node_longitude_cross_library(
        self, year, month, day, hour, desc, south_id, south_name, north_id, tol
    ):
        """Compare south node longitude between pyswisseph and libephemeris.

        If pyswisseph does not support negative body IDs with SEFLG_MOSEPH,
        this test is marked xfail to document the API difference. If both
        support it, longitudes are compared with angular_diff.
        """
        if not _SWE_SUPPORTS_NEGATIVE:
            pytest.xfail(
                f"pyswisseph does not support negative body ID {south_id} "
                f"with SEFLG_MOSEPH — libephemeris extension for Ketu"
            )

        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, south_id, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, south_id, SEFLG_MOSEPH)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{south_name} Moshier at {desc}: diff {diff:.6f}° exceeds "
            f"tolerance {tol}° (swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestSouthNodeVsNorthPlusPi:
    """Cross-library: verify south node = north node + 180° in both libs.

    Even if pyswisseph doesn't support negative IDs directly, we can
    compute south node from north node + 180° in both libraries and
    compare the results. This validates the geometric relationship
    independent of negative-ID support.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    @pytest.mark.parametrize("south_id,south_name,north_id,tol", SOUTH_NODE_IDS)
    def test_south_derived_from_north_matches(
        self, year, month, day, hour, desc, south_id, south_name, north_id, tol
    ):
        """Compare libephemeris south node with pyswisseph north+180°.

        Computes south node longitude as (north_lon + 180) % 360 using
        pyswisseph's north node, then compares with libephemeris's native
        south node calculation via negative body ID.
        """
        jd = swe.julday(year, month, day, hour)

        # pyswisseph north node
        pos_swe_north, _ = swe.calc_ut(jd, north_id, swe.FLG_MOSEPH)
        swe_south_lon = (pos_swe_north[0] + 180.0) % 360.0

        # libephemeris native south node via negative ID
        pos_py_south, _ = ephem.swe_calc_ut(jd, south_id, SEFLG_MOSEPH)

        diff = angular_diff(swe_south_lon, pos_py_south[0])

        assert diff < tol, (
            f"{south_name} at {desc}: libephemeris south lon "
            f"{pos_py_south[0]:.6f}° vs pyswisseph north+180° "
            f"{swe_south_lon:.6f}° diff={diff:.6f}° (tol={tol}°)"
        )


class TestPyswissephNegativeIdBehavior:
    """Document pyswisseph behavior with negative body IDs in Moshier mode.

    This class explicitly tests what happens when pyswisseph receives
    negative body IDs with SEFLG_MOSEPH, documenting the result for
    API migration guidance.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("south_id,south_name,north_id,tol", SOUTH_NODE_IDS)
    def test_pyswisseph_negative_id_behavior(self, south_id, south_name, north_id, tol):
        """Document whether pyswisseph supports or rejects negative body IDs."""
        jd = swe.julday(2000, 1, 1, 12.0)

        swe_error = None
        swe_result = None
        try:
            swe_result, _ = swe.calc_ut(jd, south_id, swe.FLG_MOSEPH)
        except Exception as e:
            swe_error = e

        # libephemeris always supports negative IDs
        py_result, _ = ephem.swe_calc_ut(jd, south_id, SEFLG_MOSEPH)
        assert 0 <= py_result[0] < 360, (
            f"libephemeris {south_name} lon out of range: {py_result[0]}"
        )

        if swe_error is not None:
            # pyswisseph rejects negative IDs — document the incompatibility
            pytest.xfail(
                f"pyswisseph raises {type(swe_error).__name__} for body ID "
                f"{south_id} with SEFLG_MOSEPH: {swe_error}. "
                f"libephemeris supports this as an extension (planets.py:901-912). "
                f"Applications using -10/-11 for Ketu must handle this difference."
            )
        else:
            # Both support it — verify results match
            diff = angular_diff(swe_result[0], py_result[0])
            assert diff < tol, (
                f"Both support {south_name}: diff {diff:.6f}° exceeds {tol}°"
            )
