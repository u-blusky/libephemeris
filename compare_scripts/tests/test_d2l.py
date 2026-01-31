"""
Tests for the d2l (double to long) conversion function.

Tests verify that d2l correctly rounds floating-point numbers to integers
using "round half away from zero" semantics, matching pyswisseph's swe.d2l() behavior.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


def swe_d2l_to_signed(value: int) -> int:
    """
    Convert pyswisseph's d2l result (unsigned 32-bit) to a signed integer.

    pyswisseph's swe.d2l() returns unsigned 32-bit integers due to C binding behavior.
    This function converts them back to signed integers for comparison.
    """
    if value >= 2**31:
        return value - 2**32
    return value


class TestD2lBasic:
    """Basic functionality tests for d2l."""

    def test_d2l_exported(self):
        """Test that d2l is exported from the package."""
        assert hasattr(ephem, "d2l")
        assert callable(ephem.d2l)

    def test_d2l_zero(self):
        """Test conversion of zero."""
        result = ephem.d2l(0.0)
        assert result == 0
        assert isinstance(result, int)

    def test_d2l_positive_integer(self):
        """Test conversion of positive integers (no rounding needed)."""
        assert ephem.d2l(1.0) == 1
        assert ephem.d2l(10.0) == 10
        assert ephem.d2l(100.0) == 100

    def test_d2l_negative_integer(self):
        """Test conversion of negative integers (no rounding needed)."""
        assert ephem.d2l(-1.0) == -1
        assert ephem.d2l(-10.0) == -10
        assert ephem.d2l(-100.0) == -100

    def test_d2l_returns_int(self):
        """Test that d2l returns an integer."""
        result = ephem.d2l(1.5)
        assert isinstance(result, int)


class TestD2lRoundingPositive:
    """Tests for rounding positive values."""

    def test_d2l_round_down(self):
        """Test rounding down for values < 0.5 fractional part."""
        assert ephem.d2l(1.4) == 1
        assert ephem.d2l(1.49) == 1
        assert ephem.d2l(1.499999) == 1
        assert ephem.d2l(2.3) == 2

    def test_d2l_round_up(self):
        """Test rounding up for values > 0.5 fractional part."""
        assert ephem.d2l(1.6) == 2
        assert ephem.d2l(1.51) == 2
        assert ephem.d2l(1.9) == 2
        assert ephem.d2l(2.7) == 3

    def test_d2l_half_rounds_up(self):
        """Test that 0.5 rounds away from zero (up for positive)."""
        assert ephem.d2l(0.5) == 1
        assert ephem.d2l(1.5) == 2
        assert ephem.d2l(2.5) == 3
        assert ephem.d2l(10.5) == 11


class TestD2lRoundingNegative:
    """Tests for rounding negative values."""

    def test_d2l_negative_round_toward_zero(self):
        """Test rounding toward zero for values < 0.5 fractional part."""
        assert ephem.d2l(-1.4) == -1
        assert ephem.d2l(-1.49) == -1
        assert ephem.d2l(-2.3) == -2

    def test_d2l_negative_round_away_from_zero(self):
        """Test rounding away from zero for values > 0.5 fractional part."""
        assert ephem.d2l(-1.6) == -2
        assert ephem.d2l(-1.51) == -2
        assert ephem.d2l(-1.9) == -2
        assert ephem.d2l(-2.7) == -3

    def test_d2l_negative_half_rounds_away_from_zero(self):
        """Test that -0.5 rounds away from zero (down for negative)."""
        assert ephem.d2l(-0.5) == -1
        assert ephem.d2l(-1.5) == -2
        assert ephem.d2l(-2.5) == -3
        assert ephem.d2l(-10.5) == -11


class TestD2lVsSwisseph:
    """Comparison tests with pyswisseph's swe.d2l()."""

    @pytest.mark.parametrize(
        "value",
        [
            0,
            0.0,
            1,
            -1,
            0.4,
            0.5,
            0.6,
            -0.4,
            -0.5,
            -0.6,
            1.4,
            1.5,
            1.6,
            -1.4,
            -1.5,
            -1.6,
            2.5,
            -2.5,
            3.5,
            -3.5,
            10.5,
            -10.5,
            100.49,
            100.5,
            100.51,
            -100.49,
            -100.5,
            -100.51,
            1234.567,
            -1234.567,
            0.001,
            -0.001,
            0.999,
            -0.999,
        ],
    )
    def test_d2l_matches_swisseph(self, value):
        """Test that d2l matches pyswisseph's swe.d2l()."""
        result_lib = ephem.d2l(value)
        result_swe = swe_d2l_to_signed(swe.d2l(value))

        assert result_lib == result_swe, (
            f"Mismatch for value {value}: lib={result_lib}, swe={result_swe}"
        )


class TestD2lEdgeCases:
    """Edge case tests for d2l."""

    def test_d2l_very_small_positive(self):
        """Test very small positive values."""
        result = ephem.d2l(0.0001)
        expected = swe_d2l_to_signed(swe.d2l(0.0001))
        assert result == expected

    def test_d2l_very_small_negative(self):
        """Test very small negative values."""
        result = ephem.d2l(-0.0001)
        expected = swe_d2l_to_signed(swe.d2l(-0.0001))
        assert result == expected

    def test_d2l_large_values(self):
        """Test large values."""
        result = ephem.d2l(1000000.5)
        expected = swe_d2l_to_signed(swe.d2l(1000000.5))
        assert result == expected

        result_neg = ephem.d2l(-1000000.5)
        expected_neg = swe_d2l_to_signed(swe.d2l(-1000000.5))
        assert result_neg == expected_neg

    def test_d2l_near_integer_boundaries(self):
        """Test values very close to integers."""
        # Just below integer
        assert ephem.d2l(0.9999999) == swe_d2l_to_signed(swe.d2l(0.9999999))
        assert ephem.d2l(1.9999999) == swe_d2l_to_signed(swe.d2l(1.9999999))

        # Just above integer
        assert ephem.d2l(1.0000001) == swe_d2l_to_signed(swe.d2l(1.0000001))
        assert ephem.d2l(2.0000001) == swe_d2l_to_signed(swe.d2l(2.0000001))


class TestD2lRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_values(self, random_longitudes):
        """Test with random values."""
        lons = random_longitudes(100)

        for lon in lons:
            # Test positive value
            result_lib = ephem.d2l(lon)
            result_swe = swe_d2l_to_signed(swe.d2l(lon))
            assert result_lib == result_swe, f"Mismatch for value {lon}"

            # Test negative value
            result_lib_neg = ephem.d2l(-lon)
            result_swe_neg = swe_d2l_to_signed(swe.d2l(-lon))
            assert result_lib_neg == result_swe_neg, f"Mismatch for value {-lon}"

            # Test with fractional offset
            for offset in [0.3, 0.5, 0.7]:
                val_pos = lon + offset
                val_neg = -lon - offset
                assert ephem.d2l(val_pos) == swe_d2l_to_signed(swe.d2l(val_pos)), (
                    f"Mismatch for value {val_pos}"
                )
                assert ephem.d2l(val_neg) == swe_d2l_to_signed(swe.d2l(val_neg)), (
                    f"Mismatch for value {val_neg}"
                )
