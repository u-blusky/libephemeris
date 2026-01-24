"""
Tests for the radnorm angle normalization function.

Tests verify that radnorm correctly normalizes angles to the range [0, 2*pi),
matching pyswisseph's swe.radnorm() behavior.
"""

import math

import pytest
import swisseph as swe
import libephemeris as ephem


TWO_PI = 2.0 * math.pi


class TestRadnormBasic:
    """Basic functionality tests for radnorm."""

    def test_radnorm_exported(self):
        """Test that radnorm is exported from the package."""
        assert hasattr(ephem, "radnorm")
        assert callable(ephem.radnorm)

    def test_radnorm_zero(self):
        """Test normalization of zero."""
        result = ephem.radnorm(0)
        assert result == 0.0

    def test_radnorm_positive_in_range(self):
        """Test that angles already in range [0, 2*pi) are unchanged."""
        test_values = [
            0.0,
            math.pi / 4,
            math.pi / 2,
            math.pi,
            3 * math.pi / 2,
            TWO_PI - 0.001,
        ]
        for angle in test_values:
            result = ephem.radnorm(angle)
            assert result == pytest.approx(angle), f"Failed for angle {angle}"

    def test_radnorm_two_pi(self):
        """Test that 2*pi is normalized to 0."""
        result = ephem.radnorm(TWO_PI)
        assert result == pytest.approx(0.0, abs=1e-15)

    def test_radnorm_returns_float(self):
        """Test that radnorm returns a float."""
        result = ephem.radnorm(math.pi / 4)
        assert isinstance(result, float)


class TestRadnormNegativeAngles:
    """Tests for negative angle handling."""

    def test_radnorm_negative_small(self):
        """Test normalization of small negative angles."""
        result = ephem.radnorm(-math.pi / 4)  # -45 degrees
        expected = TWO_PI - math.pi / 4  # 315 degrees
        assert result == pytest.approx(expected)

    def test_radnorm_negative_pi(self):
        """Test normalization of -pi radians."""
        result = ephem.radnorm(-math.pi)
        assert result == pytest.approx(math.pi)

    def test_radnorm_negative_large(self):
        """Test normalization of large negative angles."""
        result = ephem.radnorm(-3 * math.pi)  # -540 degrees
        expected = math.pi  # 180 degrees
        assert result == pytest.approx(expected)

    def test_radnorm_negative_multiple_rotations(self):
        """Test normalization of multiple negative rotations."""
        result = ephem.radnorm(-4 * math.pi)
        assert result == pytest.approx(0.0, abs=1e-15)

    def test_radnorm_negative_fractional(self):
        """Test normalization of negative fractional angles."""
        result = ephem.radnorm(-TWO_PI + 0.1)  # Slightly more than -360 degrees
        assert result == pytest.approx(0.1)


class TestRadnormLargeAngles:
    """Tests for large angle handling."""

    def test_radnorm_greater_than_two_pi(self):
        """Test normalization of angles greater than 2*pi."""
        result = ephem.radnorm(3 * math.pi)  # 540 degrees
        expected = math.pi  # 180 degrees
        assert result == pytest.approx(expected)

    def test_radnorm_multiple_rotations(self):
        """Test normalization of multiple rotations."""
        result = ephem.radnorm(4 * math.pi)  # 720 degrees
        assert result == pytest.approx(0.0, abs=1e-15)

    def test_radnorm_many_rotations(self):
        """Test normalization of many rotations."""
        result = ephem.radnorm(20 * math.pi)  # 10 full rotations
        assert result == pytest.approx(0.0, abs=1e-14)

    def test_radnorm_large_with_remainder(self):
        """Test normalization of large angles with remainder."""
        result = ephem.radnorm(4 * math.pi + 0.5)  # 2 rotations + 0.5 radians
        assert result == pytest.approx(0.5)


class TestRadnormVsSwisseph:
    """Comparison tests with pyswisseph's swe.radnorm()."""

    @pytest.mark.parametrize(
        "angle",
        [
            0,
            math.pi / 4,
            math.pi / 2,
            math.pi,
            3 * math.pi / 2,
            TWO_PI,
            -math.pi / 4,
            -math.pi / 2,
            -math.pi,
            -3 * math.pi / 2,
            -TWO_PI,
            3 * math.pi,
            -3 * math.pi,
            4 * math.pi,
            -4 * math.pi,
            TWO_PI - 0.001,
            -TWO_PI + 0.001,
            1e-10,
            -1e-10,
            TWO_PI - 1e-10,
            -TWO_PI + 1e-10,
            10 * math.pi,
            -10 * math.pi,
        ],
    )
    def test_radnorm_matches_swisseph(self, angle):
        """Test that radnorm matches pyswisseph's swe.radnorm()."""
        result_lib = ephem.radnorm(angle)
        result_swe = swe.radnorm(angle)

        assert result_lib == pytest.approx(result_swe, abs=1e-10), (
            f"Mismatch for angle {angle}: lib={result_lib}, swe={result_swe}"
        )


class TestRadnormEdgeCases:
    """Edge case tests for radnorm."""

    def test_radnorm_very_small_positive(self):
        """Test very small positive angles."""
        result = ephem.radnorm(1e-10)
        expected = swe.radnorm(1e-10)
        assert result == pytest.approx(expected, abs=1e-15)

    def test_radnorm_very_small_negative(self):
        """Test very small negative angles."""
        result = ephem.radnorm(-1e-10)
        expected = swe.radnorm(-1e-10)
        assert result == pytest.approx(expected, abs=1e-10)

    def test_radnorm_float_precision(self):
        """Test float precision edge cases."""
        # Near 2*pi boundary
        result = ephem.radnorm(TWO_PI - 1e-12)
        expected = swe.radnorm(TWO_PI - 1e-12)
        assert result == pytest.approx(expected, abs=1e-10)


class TestRadnormRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles(self, random_longitudes):
        """Test with random angles converted to radians."""
        lons = random_longitudes(100)

        for lon in lons:
            # Convert degrees to radians
            angle = math.radians(lon)

            # Test positive angle
            result_lib = ephem.radnorm(angle)
            result_swe = swe.radnorm(angle)
            assert result_lib == pytest.approx(result_swe, abs=1e-10), (
                f"Mismatch for angle {angle}"
            )

            # Test negative angle
            result_lib_neg = ephem.radnorm(-angle)
            result_swe_neg = swe.radnorm(-angle)
            assert result_lib_neg == pytest.approx(result_swe_neg, abs=1e-10), (
                f"Mismatch for angle {-angle}"
            )

    def test_random_large_angles(self, random_longitudes):
        """Test with random large angles in radians."""
        lons = random_longitudes(50)

        for lon in lons:
            # Convert to radians and multiply for large angles
            large_angle = math.radians(lon) * 10  # Up to ~60 radians
            result_lib = ephem.radnorm(large_angle)
            result_swe = swe.radnorm(large_angle)
            assert result_lib == pytest.approx(result_swe, abs=1e-10), (
                f"Mismatch for large angle {large_angle}"
            )

            # Test large negative angles
            result_lib_neg = ephem.radnorm(-large_angle)
            result_swe_neg = swe.radnorm(-large_angle)
            assert result_lib_neg == pytest.approx(result_swe_neg, abs=1e-10), (
                f"Mismatch for large negative angle {-large_angle}"
            )
