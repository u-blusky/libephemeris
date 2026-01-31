"""
Utility Functions Comparison Tests.

Compares utility functions between pyswisseph and libephemeris:
- Angle normalization functions
- Degree/centiseconds conversion
- String formatting functions
"""

import pytest
import swisseph as swe
import libephemeris as ephem


# ============================================================================
# TOLERANCES
# ============================================================================

ANGLE_TOL = 1e-10


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestDegnorm:
    """Compare degnorm function."""

    TEST_VALUES = [
        0.0,
        90.0,
        180.0,
        270.0,
        360.0,
        -90.0,
        -180.0,
        -360.0,
        450.0,
        720.0,
        -720.0,
        359.9999,
        0.0001,
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("angle", TEST_VALUES)
    def test_degnorm(self, angle):
        """Test degnorm normalizes to 0-360."""
        result_swe = swe.degnorm(angle)
        result_py = ephem.degnorm(angle)

        diff = abs(result_swe - result_py)

        assert diff < ANGLE_TOL, f"degnorm({angle}): diff {diff:.15f} exceeds tolerance"

        # Result should be in [0, 360)
        assert 0 <= result_py < 360, (
            f"degnorm({angle}) = {result_py}, expected [0, 360)"
        )


class TestRadnorm:
    """Compare radnorm function."""

    import math

    PI = math.pi

    TEST_VALUES = [
        0.0,
        PI / 2,
        PI,
        3 * PI / 2,
        2 * PI,
        -PI / 2,
        -PI,
        -2 * PI,
        3 * PI,
        4 * PI,
        -4 * PI,
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("angle", TEST_VALUES)
    def test_radnorm(self, angle):
        """Test radnorm normalizes to 0-2pi."""
        result_swe = swe.radnorm(angle)
        result_py = ephem.radnorm(angle)

        diff = abs(result_swe - result_py)

        assert diff < ANGLE_TOL, f"radnorm({angle}): diff {diff:.15f} exceeds tolerance"


class TestCsnorm:
    """Compare csnorm function (centiseconds normalization)."""

    TEST_VALUES = [
        0,
        360 * 3600 * 100,
        -360 * 3600 * 100,
        180 * 3600 * 100,
        -180 * 3600 * 100,
        90 * 3600 * 100,
        450 * 3600 * 100,
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("cs", TEST_VALUES)
    def test_csnorm(self, cs):
        """Test csnorm normalizes centiseconds."""
        result_swe = swe.csnorm(cs)
        result_py = ephem.csnorm(cs)

        diff = abs(result_swe - result_py)

        assert diff < 1, f"csnorm({cs}): diff {diff} exceeds tolerance"


class TestD2l:
    """Compare d2l function (degrees to centiseconds)."""

    TEST_VALUES = [0.0, 1.0, 45.0, 90.0, 180.0, 270.0, 359.9999]

    @pytest.mark.comparison
    @pytest.mark.parametrize("deg", TEST_VALUES)
    def test_d2l(self, deg):
        """Test d2l conversion."""
        result_swe = swe.d2l(deg)
        result_py = ephem.d2l(deg)

        diff = abs(result_swe - result_py)

        assert diff < 1, f"d2l({deg}): diff {diff} exceeds tolerance"


class TestDifcsn:
    """Compare difcsn function (centisecond difference)."""

    TEST_PAIRS = [
        (0, 360 * 3600 * 100),
        (90 * 3600 * 100, 270 * 3600 * 100),
        (180 * 3600 * 100, 0),
        (350 * 3600 * 100, 10 * 3600 * 100),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("cs1,cs2", TEST_PAIRS)
    def test_difcsn(self, cs1, cs2):
        """Test difcsn calculates shortest arc difference."""
        result_swe = swe.difcsn(cs1, cs2)
        result_py = ephem.difcsn(cs1, cs2)

        diff = abs(result_swe - result_py)

        assert diff < 1, f"difcsn({cs1}, {cs2}): diff {diff} exceeds tolerance"


class TestDifdegn:
    """Compare difdegn function (degree difference)."""

    TEST_PAIRS = [
        (0.0, 360.0),
        (90.0, 270.0),
        (180.0, 0.0),
        (350.0, 10.0),
        (1.0, 359.0),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("deg1,deg2", TEST_PAIRS)
    def test_difdegn(self, deg1, deg2):
        """Test difdegn calculates shortest arc difference in degrees."""
        result_swe = swe.difdegn(deg1, deg2)
        result_py = ephem.difdegn(deg1, deg2)

        diff = abs(result_swe - result_py)

        assert diff < ANGLE_TOL, (
            f"difdegn({deg1}, {deg2}): diff {diff:.15f} exceeds tolerance"
        )


class TestSplitDeg:
    """Compare split_deg function."""

    TEST_VALUES = [
        (0.0, 0),
        (45.5, 0),
        (123.456789, 0),
        (359.999999, 0),
        (-45.5, 0),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("deg,roundflag", TEST_VALUES)
    def test_split_deg(self, deg, roundflag):
        """Test split_deg breaks degrees into components."""
        result_swe = swe.split_deg(deg, roundflag)
        result_py = ephem.split_deg(deg, roundflag)

        # Compare degree component
        assert result_swe[0] == result_py[0], (
            f"split_deg({deg}): degrees differ {result_swe[0]} vs {result_py[0]}"
        )
        # Compare minutes
        assert result_swe[1] == result_py[1], (
            f"split_deg({deg}): minutes differ {result_swe[1]} vs {result_py[1]}"
        )
        # Compare seconds (with small tolerance)
        diff_sec = abs(result_swe[2] - result_py[2])
        assert diff_sec < 0.001, f"split_deg({deg}): seconds diff {diff_sec:.6f}"
