"""
Tests for the difcs2n centiseconds angle difference function.

Tests verify that difcs2n correctly calculates the signed angular difference
between two angles in centiseconds (1/100 arcsecond), normalized to the
equivalent of [-180°, +180°].

Constants:
- 1 centisecond = 1/100 arcsecond = 1/360000 degree
- 360° = 129,600,000 centiseconds (CS360)
- 180° = 64,800,000 centiseconds (CS180)
- 1° = 360,000 centiseconds
"""

import pytest
import swisseph as swe
import libephemeris as ephem

# Constants for centiseconds
CS360 = 360 * 3600 * 100  # 129600000 - full circle
CS180 = 180 * 3600 * 100  # 64800000 - half circle
CS_PER_DEG = 3600 * 100  # 360000 - centiseconds per degree


class TestDifcs2nBasic:
    """Basic functionality tests for difcs2n."""

    def test_difcs2n_exported(self):
        """Test that difcs2n is exported from the package."""
        assert hasattr(ephem, "difcs2n")
        assert callable(ephem.difcs2n)

    def test_difcs2n_zero_difference(self):
        """Test that identical angles return 0."""
        assert ephem.difcs2n(0, 0) == 0
        assert ephem.difcs2n(CS180, CS180) == 0
        assert ephem.difcs2n(CS360, CS360) == 0

    def test_difcs2n_small_positive_difference(self):
        """Test small positive differences."""
        # 2° - 1° = 1° = 360000 cs
        result = ephem.difcs2n(2 * CS_PER_DEG, 1 * CS_PER_DEG)
        assert result == CS_PER_DEG

    def test_difcs2n_small_negative_difference(self):
        """Test small negative differences."""
        # 1° - 2° = -1° = -360000 cs
        result = ephem.difcs2n(1 * CS_PER_DEG, 2 * CS_PER_DEG)
        assert result == -CS_PER_DEG

    def test_difcs2n_returns_int(self):
        """Test that difcs2n returns an integer."""
        result = ephem.difcs2n(1000, 500)
        assert isinstance(result, int)


class TestDifcs2nWraparound:
    """Tests for 360° wraparound handling."""

    def test_difcs2n_wraparound_positive(self):
        """Test wraparound when going from high to low angle."""
        # 359° - 1° should give -2° (shorter path), not 358°
        a = 359 * CS_PER_DEG  # 359°
        b = 1 * CS_PER_DEG  # 1°
        result = ephem.difcs2n(a, b)
        expected = -2 * CS_PER_DEG  # -2°
        assert result == expected

    def test_difcs2n_wraparound_negative(self):
        """Test wraparound when going from low to high angle."""
        # 1° - 359° should give +2° (shorter path), not -358°
        a = 1 * CS_PER_DEG  # 1°
        b = 359 * CS_PER_DEG  # 359°
        result = ephem.difcs2n(a, b)
        expected = 2 * CS_PER_DEG  # +2°
        assert result == expected

    def test_difcs2n_at_180_boundary(self):
        """Test behavior at exactly 180° difference."""
        # 180° - 0° = -180° (pyswisseph convention: prefer negative at 180° boundary)
        result = ephem.difcs2n(CS180, 0)
        assert result == -CS180

    def test_difcs2n_just_over_180(self):
        """Test behavior just over 180° difference."""
        # 181° - 0° should wrap to -179°
        a = 181 * CS_PER_DEG
        result = ephem.difcs2n(a, 0)
        expected = -179 * CS_PER_DEG
        assert result == expected

    def test_difcs2n_just_under_180(self):
        """Test behavior just under 180° difference."""
        # 179° - 0° = 179°
        a = 179 * CS_PER_DEG
        result = ephem.difcs2n(a, 0)
        expected = 179 * CS_PER_DEG
        assert result == expected


class TestDifcs2nNormalization:
    """Tests for input normalization."""

    def test_difcs2n_large_inputs(self):
        """Test with inputs larger than 360°."""
        # (360° + 10°) - 5° = 5°
        a = CS360 + 10 * CS_PER_DEG  # 370°
        b = 5 * CS_PER_DEG  # 5°
        result = ephem.difcs2n(a, b)
        expected = 5 * CS_PER_DEG
        assert result == expected

    def test_difcs2n_negative_inputs(self):
        """Test with negative inputs."""
        # -10° is equivalent to 350°
        # -10° - 5° = 350° - 5° = 345° -> normalized to -15°
        a = -10 * CS_PER_DEG  # -10° = 350°
        b = 5 * CS_PER_DEG  # 5°
        result = ephem.difcs2n(a, b)
        expected = -15 * CS_PER_DEG
        assert result == expected

    def test_difcs2n_both_negative(self):
        """Test with both inputs negative."""
        # -10° - (-20°) = -10° + 20° = 10°
        a = -10 * CS_PER_DEG
        b = -20 * CS_PER_DEG
        result = ephem.difcs2n(a, b)
        expected = 10 * CS_PER_DEG
        assert result == expected


class TestDifcs2nVsSwisseph:
    """Comparison tests with pyswisseph's swe.difcs2n()."""

    @pytest.mark.parametrize(
        "a_deg, b_deg",
        [
            (0, 0),
            (10, 20),
            (20, 10),
            (350, 10),
            (10, 350),
            (180, 0),
            (0, 180),
            (181, 0),
            (0, 181),
            (179, 0),
            (0, 179),
            (90, 270),
            (270, 90),
            (45, 315),
            (315, 45),
            (1, 359),
            (359, 1),
        ],
    )
    def test_difcs2n_matches_swisseph(self, a_deg, b_deg):
        """Test that difcs2n matches pyswisseph's swe.difcs2n()."""
        a_cs = a_deg * CS_PER_DEG
        b_cs = b_deg * CS_PER_DEG

        result_lib = ephem.difcs2n(a_cs, b_cs)
        result_swe = swe.difcs2n(a_cs, b_cs)

        assert result_lib == result_swe, (
            f"Mismatch for ({a_deg}°, {b_deg}°): lib={result_lib}, swe={result_swe}"
        )

    @pytest.mark.parametrize(
        "a_cs, b_cs",
        [
            (0, 0),
            (100, 200),
            (1000, 500),
            (CS180, 0),
            (CS360, CS180),
            (CS360 - 1, 1),
            (1, CS360 - 1),
            (CS180 + 1, 0),
            (CS180 - 1, 0),
        ],
    )
    def test_difcs2n_raw_centiseconds_match_swisseph(self, a_cs, b_cs):
        """Test with raw centisecond values against pyswisseph."""
        result_lib = ephem.difcs2n(a_cs, b_cs)
        result_swe = swe.difcs2n(a_cs, b_cs)

        assert result_lib == result_swe, (
            f"Mismatch for ({a_cs}, {b_cs}): lib={result_lib}, swe={result_swe}"
        )


class TestDifcs2nEdgeCases:
    """Edge case tests for difcs2n."""

    def test_difcs2n_single_centisecond(self):
        """Test with single centisecond difference."""
        result = ephem.difcs2n(1, 0)
        assert result == 1

        result = ephem.difcs2n(0, 1)
        assert result == -1

    def test_difcs2n_at_full_circle(self):
        """Test at exactly one full circle."""
        # 360° - 0° = 0° (full circle wraps to 0)
        result = ephem.difcs2n(CS360, 0)
        assert result == 0

    def test_difcs2n_multiple_circles(self):
        """Test with multiple full circles."""
        # 720° - 0° = 0°
        result = ephem.difcs2n(2 * CS360, 0)
        assert result == 0

        # 721° - 1° = 0°
        result = ephem.difcs2n(2 * CS360 + CS_PER_DEG, CS_PER_DEG)
        assert result == 0

    def test_difcs2n_very_large_values(self):
        """Test with very large centisecond values."""
        # 1000° - 0° normalized (use moderate values to avoid pyswisseph overflow)
        large_angle = 1000 * CS_PER_DEG
        result = ephem.difcs2n(large_angle, 0)
        result_swe = swe.difcs2n(large_angle, 0)
        assert result == result_swe


class TestDifcs2nRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles(self, random_longitudes):
        """Test with random angles converted to centiseconds."""
        lons = random_longitudes(100)

        for lon1 in lons[:10]:
            for lon2 in lons[10:20]:
                a_cs = int(lon1 * CS_PER_DEG)
                b_cs = int(lon2 * CS_PER_DEG)

                result_lib = ephem.difcs2n(a_cs, b_cs)
                result_swe = swe.difcs2n(a_cs, b_cs)

                assert result_lib == result_swe, f"Mismatch for ({a_cs}, {b_cs})"

    def test_random_large_centiseconds(self, random_longitudes):
        """Test with random large centisecond values."""
        lons = random_longitudes(50)

        for lon in lons:
            # Test with large values (up to 10 full circles)
            a_cs = int(lon * CS_PER_DEG * 10)
            b_cs = int((lon / 2) * CS_PER_DEG)

            result_lib = ephem.difcs2n(a_cs, b_cs)
            result_swe = swe.difcs2n(a_cs, b_cs)

            assert result_lib == result_swe, (
                f"Mismatch for large values ({a_cs}, {b_cs})"
            )


class TestDifcs2nSymmetry:
    """Tests for symmetry properties."""

    def test_antisymmetric(self):
        """Test that difcs2n(a, b) = -difcs2n(b, a) except at 180° boundary."""
        # Test cases that don't result in exactly 180° difference
        test_cases = [
            (10 * CS_PER_DEG, 20 * CS_PER_DEG),
            (350 * CS_PER_DEG, 10 * CS_PER_DEG),
            (1 * CS_PER_DEG, 179 * CS_PER_DEG),
            (45 * CS_PER_DEG, 200 * CS_PER_DEG),
        ]

        for a, b in test_cases:
            result1 = ephem.difcs2n(a, b)
            result2 = ephem.difcs2n(b, a)
            assert result1 == -result2, (
                f"Antisymmetry failed for ({a}, {b}): {result1} != -{result2}"
            )

    def test_180_degree_special_case(self):
        """Test that 180° boundary always returns -180° (pyswisseph convention)."""
        # At exactly 180° difference, both directions give -180°
        test_cases = [
            (CS180, 0),
            (0, CS180),
            (90 * CS_PER_DEG, 270 * CS_PER_DEG),
            (270 * CS_PER_DEG, 90 * CS_PER_DEG),
        ]

        for a, b in test_cases:
            result = ephem.difcs2n(a, b)
            assert result == -CS180, (
                f"180° boundary should return -180° for ({a}, {b}), got {result}"
            )

    def test_range_bounds(self):
        """Test that result is always in [-180°, +180°] range."""
        test_cases = [
            0,
            1,
            CS_PER_DEG,
            90 * CS_PER_DEG,
            CS180,
            CS180 + 1,
            270 * CS_PER_DEG,
            CS360 - 1,
            CS360,
        ]

        for a in test_cases:
            for b in test_cases:
                result = ephem.difcs2n(a, b)
                assert -CS180 <= result <= CS180, (
                    f"Result {result} out of range for ({a}, {b})"
                )
