"""
Tests for the degnorm angle normalization function.

Tests verify that degnorm correctly normalizes angles to the range [0, 360),
matching pyswisseph's swe.degnorm() behavior.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestDegnormBasic:
    """Basic functionality tests for degnorm."""

    def test_degnorm_exported(self):
        """Test that degnorm is exported from the package."""
        assert hasattr(ephem, "degnorm")
        assert callable(ephem.degnorm)

    def test_degnorm_zero(self):
        """Test normalization of zero."""
        result = ephem.degnorm(0)
        assert result == 0.0

    def test_degnorm_positive_in_range(self):
        """Test that angles already in range [0, 360) are unchanged."""
        test_values = [0.0, 45.0, 90.0, 180.0, 270.0, 359.9]
        for angle in test_values:
            result = ephem.degnorm(angle)
            assert result == pytest.approx(angle), f"Failed for angle {angle}"

    def test_degnorm_360(self):
        """Test that 360 is normalized to 0."""
        result = ephem.degnorm(360)
        assert result == 0.0

    def test_degnorm_returns_float(self):
        """Test that degnorm returns a float."""
        result = ephem.degnorm(45)
        assert isinstance(result, float)


class TestDegnormNegativeAngles:
    """Tests for negative angle handling."""

    def test_degnorm_negative_small(self):
        """Test normalization of small negative angles."""
        result = ephem.degnorm(-45)
        assert result == pytest.approx(315.0)

    def test_degnorm_negative_180(self):
        """Test normalization of -180 degrees."""
        result = ephem.degnorm(-180)
        assert result == pytest.approx(180.0)

    def test_degnorm_negative_large(self):
        """Test normalization of large negative angles."""
        result = ephem.degnorm(-370)
        assert result == pytest.approx(350.0)

    def test_degnorm_negative_multiple_rotations(self):
        """Test normalization of multiple negative rotations."""
        result = ephem.degnorm(-720)
        assert result == pytest.approx(0.0)

    def test_degnorm_negative_fractional(self):
        """Test normalization of negative fractional angles."""
        result = ephem.degnorm(-359.5)
        assert result == pytest.approx(0.5)


class TestDegnormLargeAngles:
    """Tests for large angle handling."""

    def test_degnorm_greater_than_360(self):
        """Test normalization of angles greater than 360."""
        result = ephem.degnorm(370)
        assert result == pytest.approx(10.0)

    def test_degnorm_multiple_rotations(self):
        """Test normalization of multiple rotations."""
        result = ephem.degnorm(720)
        assert result == pytest.approx(0.0)

    def test_degnorm_many_rotations(self):
        """Test normalization of many rotations."""
        result = ephem.degnorm(3600)  # 10 full rotations
        assert result == pytest.approx(0.0)

    def test_degnorm_large_with_remainder(self):
        """Test normalization of large angles with remainder."""
        result = ephem.degnorm(725)  # 2 rotations + 5 degrees
        assert result == pytest.approx(5.0)


class TestDegnormVsSwisseph:
    """Comparison tests with pyswisseph's swe.degnorm()."""

    @pytest.mark.parametrize(
        "angle",
        [
            0,
            45,
            90,
            180,
            270,
            360,
            -45,
            -90,
            -180,
            -270,
            -360,
            370,
            -370,
            720,
            -720,
            359.5,
            -359.5,
            0.001,
            -0.001,
            359.999,
            -359.999,
            1000,
            -1000,
        ],
    )
    def test_degnorm_matches_swisseph(self, angle):
        """Test that degnorm matches pyswisseph's swe.degnorm()."""
        result_lib = ephem.degnorm(angle)
        result_swe = swe.degnorm(angle)

        assert result_lib == pytest.approx(result_swe, abs=1e-10), (
            f"Mismatch for angle {angle}: lib={result_lib}, swe={result_swe}"
        )


class TestDegnormEdgeCases:
    """Edge case tests for degnorm."""

    def test_degnorm_very_small_positive(self):
        """Test very small positive angles."""
        result = ephem.degnorm(1e-10)
        expected = swe.degnorm(1e-10)
        assert result == pytest.approx(expected, abs=1e-15)

    def test_degnorm_very_small_negative(self):
        """Test very small negative angles."""
        result = ephem.degnorm(-1e-10)
        expected = swe.degnorm(-1e-10)
        assert result == pytest.approx(expected, abs=1e-10)

    def test_degnorm_float_precision(self):
        """Test float precision edge cases."""
        # Near 360 boundary
        result = ephem.degnorm(359.9999999999)
        expected = swe.degnorm(359.9999999999)
        assert result == pytest.approx(expected, abs=1e-10)


class TestDegnormRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles(self, random_longitudes):
        """Test with random angles including negative values."""
        lons = random_longitudes(100)

        for lon in lons:
            # Test positive angle
            result_lib = ephem.degnorm(lon)
            result_swe = swe.degnorm(lon)
            assert result_lib == pytest.approx(result_swe, abs=1e-10), (
                f"Mismatch for angle {lon}"
            )

            # Test negative angle
            result_lib_neg = ephem.degnorm(-lon)
            result_swe_neg = swe.degnorm(-lon)
            assert result_lib_neg == pytest.approx(result_swe_neg, abs=1e-10), (
                f"Mismatch for angle {-lon}"
            )

    def test_random_large_angles(self, random_longitudes):
        """Test with random large angles."""
        lons = random_longitudes(50)

        for lon in lons:
            # Test large positive angles
            large_angle = lon * 10  # Up to 3600 degrees
            result_lib = ephem.degnorm(large_angle)
            result_swe = swe.degnorm(large_angle)
            assert result_lib == pytest.approx(result_swe, abs=1e-10), (
                f"Mismatch for large angle {large_angle}"
            )

            # Test large negative angles
            result_lib_neg = ephem.degnorm(-large_angle)
            result_swe_neg = swe.degnorm(-large_angle)
            assert result_lib_neg == pytest.approx(result_swe_neg, abs=1e-10), (
                f"Mismatch for large negative angle {-large_angle}"
            )
