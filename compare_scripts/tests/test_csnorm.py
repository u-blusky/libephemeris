"""
Tests for the csnorm centiseconds normalization function.

Tests verify that csnorm correctly normalizes angles in centiseconds to the
range [0, 129600000) (equivalent to [0°, 360°)), matching pyswisseph's
swe.csnorm() behavior.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


# Constants for centiseconds calculations
CS360 = 360 * 3600 * 100  # 129600000 centiseconds in a full circle
CS180 = 180 * 3600 * 100  # 64800000 centiseconds in a half circle
CS1 = 3600 * 100  # 360000 centiseconds per degree


class TestCsnormBasic:
    """Basic functionality tests for csnorm."""

    def test_csnorm_exported(self):
        """Test that csnorm is exported from the package."""
        assert hasattr(ephem, "csnorm")
        assert callable(ephem.csnorm)

    def test_csnorm_zero(self):
        """Test normalization of zero."""
        result = ephem.csnorm(0)
        assert result == 0

    def test_csnorm_positive_in_range(self):
        """Test that values already in range [0, CS360) are unchanged."""
        test_values = [0, CS1, CS1 * 45, CS1 * 90, CS1 * 180, CS1 * 270, CS360 - 1]
        for value in test_values:
            result = ephem.csnorm(value)
            assert result == value, f"Failed for value {value}"

    def test_csnorm_360_degrees(self):
        """Test that 360° (CS360) is normalized to 0."""
        result = ephem.csnorm(CS360)
        assert result == 0

    def test_csnorm_returns_int(self):
        """Test that csnorm returns an integer."""
        result = ephem.csnorm(CS1 * 45)
        assert isinstance(result, int)


class TestCsnormNegativeValues:
    """Tests for negative value handling."""

    def test_csnorm_negative_1_degree(self):
        """Test normalization of -1° in centiseconds."""
        result = ephem.csnorm(-CS1)
        assert result == CS360 - CS1  # 359° = 129240000 cs

    def test_csnorm_negative_45_degrees(self):
        """Test normalization of -45° in centiseconds."""
        result = ephem.csnorm(-CS1 * 45)
        assert result == CS360 - (CS1 * 45)  # 315°

    def test_csnorm_negative_180_degrees(self):
        """Test normalization of -180° in centiseconds."""
        result = ephem.csnorm(-CS180)
        assert result == CS180  # 180°

    def test_csnorm_negative_large(self):
        """Test normalization of large negative values."""
        # -370° = 360° - 10° = 350°
        result = ephem.csnorm(-CS1 * 370)
        assert result == CS1 * 350

    def test_csnorm_negative_multiple_rotations(self):
        """Test normalization of multiple negative rotations."""
        result = ephem.csnorm(-CS360 * 2)  # -720°
        assert result == 0


class TestCsnormLargeValues:
    """Tests for large value handling."""

    def test_csnorm_greater_than_360(self):
        """Test normalization of values greater than 360°."""
        # 370° = 10°
        result = ephem.csnorm(CS1 * 370)
        assert result == CS1 * 10

    def test_csnorm_multiple_rotations(self):
        """Test normalization of multiple rotations."""
        result = ephem.csnorm(CS360 * 2)  # 720°
        assert result == 0

    def test_csnorm_many_rotations(self):
        """Test normalization of many rotations."""
        result = ephem.csnorm(CS360 * 10)  # 3600°
        assert result == 0

    def test_csnorm_large_with_remainder(self):
        """Test normalization of large values with remainder."""
        # 725° = 2 rotations + 5°
        result = ephem.csnorm(CS1 * 725)
        assert result == CS1 * 5


class TestCsnormVsSwisseph:
    """Comparison tests with pyswisseph's swe.csnorm()."""

    @pytest.mark.parametrize(
        "value",
        [
            0,
            CS1,  # 1°
            CS1 * 45,  # 45°
            CS1 * 90,  # 90°
            CS180,  # 180°
            CS1 * 270,  # 270°
            CS360,  # 360°
            -CS1,  # -1°
            -CS1 * 45,  # -45°
            -CS180,  # -180°
            -CS1 * 270,  # -270°
            -CS360,  # -360°
            CS1 * 370,  # 370°
            -CS1 * 370,  # -370°
            CS360 * 2,  # 720°
            -CS360 * 2,  # -720°
            1,  # 1 centisecond
            -1,  # -1 centisecond
            CS360 - 1,  # Just under 360°
            -CS360 + 1,  # Just over -360°
            CS1 * 1000,  # 1000°
            -CS1 * 1000,  # -1000°
        ],
    )
    def test_csnorm_matches_swisseph(self, value):
        """Test that csnorm matches pyswisseph's swe.csnorm()."""
        result_lib = ephem.csnorm(value)
        result_swe = swe.csnorm(value)

        assert result_lib == result_swe, (
            f"Mismatch for value {value}: lib={result_lib}, swe={result_swe}"
        )


class TestCsnormEdgeCases:
    """Edge case tests for csnorm."""

    def test_csnorm_at_boundary(self):
        """Test at exactly 360°."""
        result = ephem.csnorm(CS360)
        expected = swe.csnorm(CS360)
        assert result == expected

    def test_csnorm_just_below_boundary(self):
        """Test just below 360°."""
        result = ephem.csnorm(CS360 - 1)
        expected = swe.csnorm(CS360 - 1)
        assert result == expected

    def test_csnorm_just_above_boundary(self):
        """Test just above 360°."""
        result = ephem.csnorm(CS360 + 1)
        expected = swe.csnorm(CS360 + 1)
        assert result == expected

    def test_csnorm_negative_boundary(self):
        """Test at exactly -360°."""
        result = ephem.csnorm(-CS360)
        expected = swe.csnorm(-CS360)
        assert result == expected

    def test_csnorm_negative_just_above(self):
        """Test just above -360°."""
        result = ephem.csnorm(-CS360 + 1)
        expected = swe.csnorm(-CS360 + 1)
        assert result == expected


class TestCsnormRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles_in_centiseconds(self, random_longitudes):
        """Test with random angles converted to centiseconds."""
        lons = random_longitudes(100)

        for lon in lons:
            # Convert longitude to centiseconds
            cs_value = int(lon * CS1)

            # Test positive value
            result_lib = ephem.csnorm(cs_value)
            result_swe = swe.csnorm(cs_value)
            assert result_lib == result_swe, f"Mismatch for cs value {cs_value}"

            # Test negative value
            result_lib_neg = ephem.csnorm(-cs_value)
            result_swe_neg = swe.csnorm(-cs_value)
            assert result_lib_neg == result_swe_neg, (
                f"Mismatch for negative cs value {-cs_value}"
            )

    def test_random_large_angles(self, random_longitudes):
        """Test with random large angles."""
        lons = random_longitudes(50)

        for lon in lons:
            # Test large positive angles
            large_cs = int(lon * 10 * CS1)  # Up to 3600 degrees
            result_lib = ephem.csnorm(large_cs)
            result_swe = swe.csnorm(large_cs)
            assert result_lib == result_swe, f"Mismatch for large cs value {large_cs}"

            # Test large negative angles
            result_lib_neg = ephem.csnorm(-large_cs)
            result_swe_neg = swe.csnorm(-large_cs)
            assert result_lib_neg == result_swe_neg, (
                f"Mismatch for large negative cs value {-large_cs}"
            )
