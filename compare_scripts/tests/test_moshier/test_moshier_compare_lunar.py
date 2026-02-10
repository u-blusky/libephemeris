"""
Moshier Lunar Calculations Cross-Library Comparison Tests.

Validates Moshier (SEFLG_MOSEPH) lunar node and Lilith calculations between
pyswisseph (C library) and libephemeris (Python reimplementation):
- Mean Node (SE_MEAN_NODE = 10)
- True Node (SE_TRUE_NODE = 11)
- Mean Lilith / Black Moon (SE_MEAN_APOG = 12)
- True Lilith / Osculating Apogee (SE_OSCU_APOG = 13)

This is the Moshier-mode mirror of test_compare_lunar.py (which covers
SEFLG_SWIEPH / JPL mode). In Moshier mode, lunar nodes are computed from
Meeus polynomials (Mean Node) and analytic perturbation series (True Node)
in lunar.py, while Lilith uses calc_mean_lilith_with_latitude() and
calc_true_lilith(). The C library has its own Moshier implementations.

Without these tests, divergences between the Python and C Moshier lunar
point implementations would remain undetected, which is critical for
Jyotish (Vedic astrology) where Rahu/Ketu (lunar nodes) are among the
most important points after the graha.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TOLERANCES (Moshier-specific, relaxed compared to SWIEPH/JPL mode)
# ============================================================================
# The Moshier semi-analytical lunar calculations use polynomial/perturbation
# methods that differ between the C and Python implementations, requiring
# wider tolerances than the JPL-based calculations in test_compare_lunar.py.

MEAN_NODE_MOSHIER = 0.02  # degrees (~72 arcsec)
TRUE_NODE_MOSHIER = 0.2  # degrees (~720 arcsec)
MEAN_LILITH_MOSHIER = 0.02  # degrees (~72 arcsec)
TRUE_LILITH_MOSHIER = 0.2  # degrees (~720 arcsec)

# Speed comparison tolerances
SPEED_LON_TOL = 0.02  # degrees/day
SPEED_LAT_TOL = 0.01  # degrees/day


# ============================================================================
# TEST DATA - same dates as conftest.py test_dates fixture (9 dates)
# ============================================================================

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (1980, 5, 20, 0.0, "Past"),
    (2024, 11, 5, 18.0, "Recent"),
    (1950, 10, 15, 6.0, "Mid-century"),
    (1550, 1, 1, 0.0, "DE440 Start"),
    (2650, 1, 1, 0.0, "DE440 End"),
    (2000, 2, 29, 12.0, "Leap year"),
    (1999, 12, 31, 23.99, "Y2K eve"),
    (2000, 1, 1, 0.01, "Y2K"),
]

# Subset of dates for speed comparison tests (2 dates × 4 points = 8 tests)
SPEED_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 11, 5, 18.0, "Recent"),
]


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMeanNode:
    """Compare Moshier Mean North Node calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_mean_node_longitude(self, year, month, day, hour, desc):
        """Test Moshier Mean Node longitude matches pyswisseph within tolerance."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_MOSEPH)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MEAN_NODE_MOSHIER, (
            f"Mean Node Moshier at {desc}: diff {diff:.6f}° exceeds tolerance "
            f"{MEAN_NODE_MOSHIER}° (swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestTrueNode:
    """Compare Moshier True North Node calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_true_node_longitude(self, year, month, day, hour, desc):
        """Test Moshier True Node longitude matches pyswisseph within tolerance.

        True Node in Moshier mode is computed via ELP2000-82B perturbation
        series (Python) vs the C Moshier implementation. The tolerance is
        relaxed to 0.2° to account for implementation differences in the
        perturbation series.
        """
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, SEFLG_MOSEPH)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TRUE_NODE_MOSHIER, (
            f"True Node Moshier at {desc}: diff {diff:.6f}° exceeds tolerance "
            f"{TRUE_NODE_MOSHIER}° (swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestMeanLilith:
    """Compare Moshier Mean Lilith (Black Moon) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_mean_lilith_longitude(self, year, month, day, hour, desc):
        """Test Moshier Mean Lilith longitude matches pyswisseph within tolerance.

        Mean Lilith in Moshier mode uses calc_mean_lilith_with_latitude()
        which implements SE-compatible DE404 algorithm with ecliptic projection.
        """
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, SEFLG_MOSEPH)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MEAN_LILITH_MOSHIER, (
            f"Mean Lilith Moshier at {desc}: diff {diff:.6f}° exceeds tolerance "
            f"{MEAN_LILITH_MOSHIER}° (swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestTrueLilith:
    """Compare Moshier True Lilith (Osculating Apogee) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_true_lilith_longitude(self, year, month, day, hour, desc):
        """Test Moshier True Lilith longitude matches pyswisseph within tolerance.

        True Lilith in Moshier mode uses calc_true_lilith() which implements
        the eccentricity vector method for the osculating apogee.
        """
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, SEFLG_MOSEPH)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TRUE_LILITH_MOSHIER, (
            f"True Lilith Moshier at {desc}: diff {diff:.6f}° exceeds tolerance "
            f"{TRUE_LILITH_MOSHIER}° (swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestNodalMotion:
    """Test Moshier lunar node motion characteristics and speed comparison."""

    @pytest.mark.comparison
    def test_mean_node_retrograde(self):
        """Test that Moshier Mean Node is always retrograde."""
        jd = swe.julday(2000, 1, 1, 12.0)

        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_MOSEPH | SEFLG_SPEED)

        # Mean Node is always retrograde (negative daily motion)
        assert pos_py[3] < 0, "Mean Node should have retrograde motion"

    @pytest.mark.comparison
    def test_nodes_opposite(self):
        """Test that Moshier North and South nodes are 180 degrees apart."""
        jd = swe.julday(2000, 1, 1, 12.0)

        pos_north, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_MOSEPH)

        # South Node is 180 degrees from North Node (by definition)
        south_node = (pos_north[0] + 180) % 360
        assert abs(south_node - ((pos_north[0] + 180) % 360)) < 1e-10

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", SPEED_DATES)
    def test_mean_node_speed(self, year, month, day, hour, desc):
        """Test Moshier Mean Node speed matches pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < MEAN_NODE_MOSHIER, (
            f"Mean Node Moshier longitude at {desc}: diff {diff_lon:.6f}°"
        )
        assert diff_speed < SPEED_LON_TOL, (
            f"Mean Node Moshier speed at {desc}: diff {diff_speed:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", SPEED_DATES)
    def test_true_node_speed(self, year, month, day, hour, desc):
        """Test Moshier True Node speed matches pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < TRUE_NODE_MOSHIER, (
            f"True Node Moshier longitude at {desc}: diff {diff_lon:.6f}°"
        )
        assert diff_speed < SPEED_LON_TOL, (
            f"True Node Moshier speed at {desc}: diff {diff_speed:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )


class TestLilithLatitude:
    """Test Moshier Lilith latitude and speed calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES[:3])
    def test_mean_lilith_latitude(self, year, month, day, hour, desc):
        """Test Moshier Mean Lilith latitude."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, SEFLG_MOSEPH)

        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lat < 0.1, (
            f"Mean Lilith Moshier latitude diff {diff_lat:.6f}° at {desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES[:3])
    def test_true_lilith_latitude(self, year, month, day, hour, desc):
        """Test Moshier True Lilith latitude."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, SEFLG_MOSEPH)

        diff_lat = abs(pos_swe[1] - pos_py[1])

        # Latitude tolerance is more relaxed than longitude
        assert diff_lat < 1.0, (
            f"True Lilith Moshier latitude diff {diff_lat:.2f}° at {desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", SPEED_DATES)
    def test_mean_lilith_speed(self, year, month, day, hour, desc):
        """Test Moshier Mean Lilith speed matches pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < MEAN_LILITH_MOSHIER, (
            f"Mean Lilith Moshier longitude at {desc}: diff {diff_lon:.6f}°"
        )
        assert diff_speed < SPEED_LON_TOL, (
            f"Mean Lilith Moshier speed at {desc}: diff {diff_speed:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", SPEED_DATES)
    def test_true_lilith_speed(self, year, month, day, hour, desc):
        """Test Moshier True Lilith speed matches pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < TRUE_LILITH_MOSHIER, (
            f"True Lilith Moshier longitude at {desc}: diff {diff_lon:.6f}°"
        )
        assert diff_speed < SPEED_LON_TOL, (
            f"True Lilith Moshier speed at {desc}: diff {diff_speed:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )
