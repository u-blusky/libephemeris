"""
Tests for the difrad2n angular difference function.

Tests verify that difrad2n correctly calculates the angular difference
in radians (p1 - p2) normalized to the range [-π, π].
"""

import math
import pytest
import swisseph as swe
import libephemeris as ephem


class TestDifrad2nBasic:
    """Basic functionality tests for difrad2n."""

    def test_difrad2n_exported(self):
        """Test that difrad2n is exported from the package."""
        assert hasattr(ephem, "difrad2n")
        assert callable(ephem.difrad2n)

    def test_difrad2n_same_angle(self):
        """Test distance between same angles is zero."""
        assert ephem.difrad2n(0, 0) == 0.0
        assert ephem.difrad2n(math.pi, math.pi) == 0.0
        assert ephem.difrad2n(1.5, 1.5) == 0.0

    def test_difrad2n_basic_difference(self):
        """Test basic angle differences."""
        # 0.2 - 0.1 = 0.1
        assert ephem.difrad2n(0.2, 0.1) == pytest.approx(0.1)
        # 0.1 - 0.2 = -0.1
        assert ephem.difrad2n(0.1, 0.2) == pytest.approx(-0.1)

    def test_difrad2n_returns_float(self):
        """Test that difrad2n returns a float."""
        result = ephem.difrad2n(0.5, 1.0)
        assert isinstance(result, float)


class TestDifrad2nNormalization:
    """Tests to verify difrad2n returns values in [-π, π]."""

    def test_difrad2n_positive_result(self):
        """Test positive results when p1 > p2."""
        result = ephem.difrad2n(1.0, 0.5)
        assert result == pytest.approx(0.5)
        assert -math.pi <= result <= math.pi

    def test_difrad2n_negative_result(self):
        """Test negative results when p1 < p2."""
        result = ephem.difrad2n(0.5, 1.0)
        assert result == pytest.approx(-0.5)
        assert -math.pi <= result <= math.pi

    def test_difrad2n_at_pi(self):
        """Test exactly π difference."""
        # When diff is exactly π, pyswisseph returns -π
        result = ephem.difrad2n(math.pi, 0)
        assert result == pytest.approx(-math.pi)

    def test_difrad2n_at_minus_pi(self):
        """Test that differences beyond π are wrapped."""
        # 0 - π = -π
        result = ephem.difrad2n(0, math.pi)
        assert result == pytest.approx(-math.pi)


class TestDifrad2nWraparound:
    """Tests for wraparound handling."""

    def test_difrad2n_across_zero(self):
        """Test angles that wrap around 0/2π."""
        # Near 2*pi, crossing zero
        # 6.0 is just under 2*pi (6.28...), 0.2 is just above 0
        # 6.0 - 0.2 = 5.8, which is > π, so wrap: 5.8 - 2π ≈ -0.483
        result = ephem.difrad2n(6.0, 0.2)
        expected = (6.0 - 0.2) % (2 * math.pi)
        if expected > math.pi:
            expected -= 2 * math.pi
        assert result == pytest.approx(expected)

        # Opposite direction: 0.2 - 6.0 = -5.8 -> mod -> ~0.483
        result = ephem.difrad2n(0.2, 6.0)
        expected = (0.2 - 6.0) % (2 * math.pi)
        if expected > math.pi:
            expected -= 2 * math.pi
        assert result == pytest.approx(expected)

    def test_difrad2n_near_pi(self):
        """Test angles near π radian separation."""
        # 0 - (π - 0.1) should give a small positive value after wrapping
        result = ephem.difrad2n(0, math.pi - 0.1)
        expected = (0 - (math.pi - 0.1)) % (2 * math.pi)
        if expected > math.pi:
            expected -= 2 * math.pi
        assert result == pytest.approx(expected)

    def test_difrad2n_normalized_angles(self):
        """Test with angles beyond 0-2π range."""
        two_pi = 2 * math.pi

        # 2π + 0.5 should be equivalent to 0.5
        result = ephem.difrad2n(two_pi + 0.5, 0.5)
        assert result == pytest.approx(0.0)

        # -π/2 should be equivalent to 3π/2
        result = ephem.difrad2n(-math.pi / 2, math.pi)
        expected = (-math.pi / 2 - math.pi) % two_pi
        if expected > math.pi:
            expected -= two_pi
        assert result == pytest.approx(expected)


class TestDifrad2nVsSwisseph:
    """Comparison tests with pyswisseph's swe.difrad2n()."""

    @pytest.mark.parametrize(
        "p1,p2",
        [
            (0, 0),
            (0.1, 0.2),
            (0.2, 0.1),
            (6.0, 0.2),
            (0.2, 6.0),
            (0, math.pi),
            (math.pi, 0),
            (math.pi / 2, 3 * math.pi / 2),
            (math.pi / 4, 5 * math.pi / 4),
            (0, 2 * math.pi - 0.1),
            (2 * math.pi - 0.1, 0),
            (-0.5, 0.5),
            (7.0, 0.5),
            (0.5, 7.0),
            (-6.0, 6.0),
            (4 * math.pi, 2 * math.pi),
        ],
    )
    def test_difrad2n_matches_swisseph(self, p1, p2):
        """Test that difrad2n matches pyswisseph's swe.difrad2n()."""
        result_lib = ephem.difrad2n(p1, p2)
        result_swe = swe.difrad2n(p1, p2)

        assert result_lib == pytest.approx(result_swe, abs=1e-10), (
            f"Mismatch for ({p1}, {p2}): lib={result_lib}, swe={result_swe}"
        )


class TestDifrad2nVsDifdeg2n:
    """Tests comparing difrad2n with difdeg2n (degree equivalent)."""

    @pytest.mark.parametrize(
        "deg1,deg2",
        [
            (0, 0),
            (10, 20),
            (20, 10),
            (350, 10),
            (10, 350),
            (90, 0),
            (0, 90),
            (270, 0),
            (0, 270),
        ],
    )
    def test_difrad2n_consistent_with_difdeg2n(self, deg1, deg2):
        """Test that difrad2n gives equivalent results to difdeg2n when converted.

        Note: We exclude 180-degree differences from this test because
        difdeg2n returns +180 while difrad2n returns -π (matching pyswisseph).
        Both are valid representations of the same angular distance.
        """
        # Convert degrees to radians
        rad1 = math.radians(deg1)
        rad2 = math.radians(deg2)

        # Calculate differences
        result_deg = ephem.difdeg2n(deg1, deg2)
        result_rad = ephem.difrad2n(rad1, rad2)

        # Convert radian result to degrees for comparison
        result_rad_as_deg = math.degrees(result_rad)

        assert result_rad_as_deg == pytest.approx(result_deg, abs=1e-10), (
            f"Inconsistency for ({deg1}°, {deg2}°): "
            f"difdeg2n={result_deg}, difrad2n in deg={result_rad_as_deg}"
        )


class TestDifrad2nEdgeCases:
    """Edge case tests for difrad2n."""

    def test_difrad2n_exactly_pi(self):
        """Test exactly π radian separation."""
        # When diff is exactly π, pyswisseph returns -π
        assert ephem.difrad2n(0, math.pi) == pytest.approx(-math.pi)
        assert ephem.difrad2n(math.pi, 2 * math.pi) == pytest.approx(-math.pi)
        assert ephem.difrad2n(math.pi, 0) == pytest.approx(-math.pi)

    def test_difrad2n_very_small_difference(self):
        """Test very small angular differences."""
        result = ephem.difrad2n(0, 1e-10)
        expected = swe.difrad2n(0, 1e-10)
        assert result == pytest.approx(expected, abs=1e-15)

    def test_difrad2n_range_always_minus_pi_to_pi(self):
        """Test that result is always in range [-π, π]."""
        test_pairs = [
            (0, 0),
            (0, math.pi / 2),
            (0, math.pi),
            (0, 3 * math.pi / 2),
            (0, 2 * math.pi - 0.1),
            (math.pi / 4, 3 * math.pi / 4),
            (math.pi / 4, 5 * math.pi / 4),
            (math.pi / 4, 7 * math.pi / 4),
            (math.pi / 2, 3 * math.pi / 2),
            (math.pi, 0),
            (3 * math.pi / 2, math.pi / 2),
            (6.0, 0.2),
            (0.2, 6.0),
        ]
        for p1, p2 in test_pairs:
            result = ephem.difrad2n(p1, p2)
            assert -math.pi <= result <= math.pi, (
                f"difrad2n({p1}, {p2}) = {result} is outside [-π, π]"
            )


class TestDifrad2nRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles(self, random_longitudes):
        """Test with random angles converted to radians."""
        lons = random_longitudes(100)

        for i in range(0, len(lons) - 1, 2):
            # Convert longitude degrees to radians
            p1 = math.radians(lons[i])
            p2 = math.radians(lons[i + 1])

            result_lib = ephem.difrad2n(p1, p2)
            result_swe = swe.difrad2n(p1, p2)

            assert result_lib == pytest.approx(result_swe, abs=1e-10), (
                f"Mismatch for ({p1}, {p2})"
            )
            assert -math.pi <= result_lib <= math.pi, (
                f"Result outside [-π, π] for ({p1}, {p2})"
            )

    def test_consistency_with_formula(self, random_longitudes):
        """Test that difrad2n follows the expected formula."""
        lons = random_longitudes(50)
        two_pi = 2 * math.pi

        for i in range(0, len(lons) - 1, 2):
            p1 = math.radians(lons[i])
            p2 = math.radians(lons[i + 1])

            result = ephem.difrad2n(p1, p2)

            # Expected: (p1 - p2) mod 2π, then subtract 2π if > π
            expected = (p1 - p2) % two_pi
            if expected > math.pi:
                expected -= two_pi

            assert result == pytest.approx(expected, abs=1e-10), (
                f"difrad2n({p1}, {p2}) = {result} != expected {expected}"
            )
