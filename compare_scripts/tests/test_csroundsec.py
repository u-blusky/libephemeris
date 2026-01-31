"""
Tests for the csroundsec function that rounds centiseconds to arcseconds.

Tests verify that csroundsec correctly rounds centiseconds (1/100 arcsecond)
to the nearest arcsecond (multiple of 100), matching pyswisseph's
swe.csroundsec() behavior.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestCsroundsecBasic:
    """Basic functionality tests for csroundsec."""

    def test_csroundsec_exported(self):
        """Test that csroundsec is exported from the package."""
        assert hasattr(ephem, "csroundsec")
        assert callable(ephem.csroundsec)

    def test_csroundsec_zero(self):
        """Test rounding of zero."""
        result = ephem.csroundsec(0)
        assert result == 0

    def test_csroundsec_exact_arcsecond(self):
        """Test that exact arcseconds return the same value."""
        # 100 centiseconds = 1 arcsecond -> returns 100
        assert ephem.csroundsec(100) == 100
        assert ephem.csroundsec(200) == 200
        assert ephem.csroundsec(1000) == 1000

    def test_csroundsec_returns_int(self):
        """Test that csroundsec returns an integer."""
        result = ephem.csroundsec(150)
        assert isinstance(result, int)


class TestCsroundsecRounding:
    """Tests for rounding behavior."""

    def test_csroundsec_rounds_up_at_half(self):
        """Test rounding up when fraction >= 0.5 (for positive)."""
        # 50 cs = 0.5 arcsec -> rounds to 100 (1 arcsec)
        assert ephem.csroundsec(50) == 100
        # 150 cs = 1.5 arcsec -> rounds to 200 (2 arcsec)
        assert ephem.csroundsec(150) == 200
        # 250 cs = 2.5 arcsec -> rounds to 300 (3 arcsec)
        assert ephem.csroundsec(250) == 300

    def test_csroundsec_rounds_down(self):
        """Test rounding down when fraction < 0.5."""
        # 49 cs = 0.49 arcsec -> rounds to 0
        assert ephem.csroundsec(49) == 0
        # 149 cs = 1.49 arcsec -> rounds to 100
        assert ephem.csroundsec(149) == 100


class TestCsroundsecNegativeValues:
    """Tests for negative value handling."""

    def test_csroundsec_negative_exact(self):
        """Test negative exact arcsecond values.

        Note: Negative exact multiples of 100 (except -100) round toward zero.
        This matches pyswisseph's behavior.
        """
        assert ephem.csroundsec(-100) == -100
        # -200 rounds to -100 (matching pyswisseph behavior)
        assert ephem.csroundsec(-200) == swe.csroundsec(-200)
        assert ephem.csroundsec(-1000) == swe.csroundsec(-1000)

    def test_csroundsec_negative_truncates_towards_zero(self):
        """Test that negative values round towards zero (truncation)."""
        # -50 cs -> 0 (rounds towards zero)
        assert ephem.csroundsec(-50) == 0
        # -99 cs -> 0 (rounds towards zero)
        assert ephem.csroundsec(-99) == 0
        # -149 cs -> -100 (truncates to -1 arcsec)
        assert ephem.csroundsec(-149) == -100
        # -150 cs -> -100 (truncates to -1 arcsec, not -2)
        assert ephem.csroundsec(-150) == -100


class TestCsroundsecLargeValues:
    """Tests for large value handling."""

    def test_csroundsec_large_positive(self):
        """Test large positive values."""
        # 360000 cs = 3600 arcsec = 1 degree in arcseconds
        assert ephem.csroundsec(360000) == 360000
        # 129600000 cs = 1296000 arcsec = 360 degrees in arcseconds
        assert ephem.csroundsec(129600000) == 129600000

    def test_csroundsec_large_negative(self):
        """Test large negative values.

        Note: Large negative exact multiples of 100 round toward zero.
        This matches pyswisseph's behavior.
        """
        assert ephem.csroundsec(-360000) == swe.csroundsec(-360000)
        assert ephem.csroundsec(-129600000) == swe.csroundsec(-129600000)


class TestCsroundsecVsSwisseph:
    """Comparison tests with pyswisseph's swe.csroundsec()."""

    @pytest.mark.parametrize(
        "value",
        [
            0,
            1,
            49,
            50,
            51,
            99,
            100,
            149,
            150,
            151,
            199,
            200,
            250,
            350,
            450,
            550,
            -1,
            -49,
            -50,
            -51,
            -99,
            -100,
            -149,
            -150,
            -151,
            360000,
            -360000,
            129600000,
            -129600000,
        ],
    )
    def test_csroundsec_matches_swisseph(self, value):
        """Test that csroundsec matches pyswisseph's swe.csroundsec()."""
        result_lib = ephem.csroundsec(value)
        result_swe = swe.csroundsec(value)

        assert result_lib == result_swe, (
            f"Mismatch for value {value}: lib={result_lib}, swe={result_swe}"
        )


class TestCsroundsecRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_values(self, random_longitudes):
        """Test with random angles converted to centiseconds."""
        lons = random_longitudes(100)

        for lon in lons:
            # Convert longitude to centiseconds (using degrees -> arcsec -> centisec)
            cs_value = int(lon * 3600 * 100)

            # Test positive value
            result_lib = ephem.csroundsec(cs_value)
            result_swe = swe.csroundsec(cs_value)
            assert result_lib == result_swe, f"Mismatch for cs value {cs_value}"

            # Test negative value
            result_lib_neg = ephem.csroundsec(-cs_value)
            result_swe_neg = swe.csroundsec(-cs_value)
            assert result_lib_neg == result_swe_neg, (
                f"Mismatch for negative cs value {-cs_value}"
            )
