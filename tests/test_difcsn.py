"""
Tests for the difcsn centiseconds angle difference function.

Tests verify that difcsn correctly calculates the angular difference
between two angles in centiseconds (1/100 arcsecond), normalized to the
equivalent of [0°, 360°) - always returning a positive value.

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


class TestDifcsnBasic:
    """Basic functionality tests for difcsn."""

    def test_difcsn_exported(self):
        """Test that difcsn is exported from the package."""
        assert hasattr(ephem, "difcsn")
        assert callable(ephem.difcsn)

    def test_difcsn_zero_difference(self):
        """Test that identical angles return 0."""
        assert ephem.difcsn(0, 0) == 0
        assert ephem.difcsn(CS180, CS180) == 0
        assert ephem.difcsn(CS360, CS360) == 0

    def test_difcsn_small_positive_difference(self):
        """Test small positive differences."""
        # 2° - 1° = 1° = 360000 cs
        result = ephem.difcsn(2 * CS_PER_DEG, 1 * CS_PER_DEG)
        assert result == CS_PER_DEG

    def test_difcsn_small_negative_becomes_positive(self):
        """Test that negative differences become positive (wrap to 360°)."""
        # 1° - 2° = -1° -> 359° = 129240000 cs
        result = ephem.difcsn(1 * CS_PER_DEG, 2 * CS_PER_DEG)
        expected = 359 * CS_PER_DEG
        assert result == expected

    def test_difcsn_returns_int(self):
        """Test that difcsn returns an integer."""
        result = ephem.difcsn(1000, 500)
        assert isinstance(result, int)


class TestDifcsnWraparound:
    """Tests for 360° wraparound handling."""

    def test_difcsn_large_angle_difference(self):
        """Test large angle difference."""
        # 359° - 1° = 358°
        a = 359 * CS_PER_DEG  # 359°
        b = 1 * CS_PER_DEG  # 1°
        result = ephem.difcsn(a, b)
        expected = 358 * CS_PER_DEG  # 358°
        assert result == expected

    def test_difcsn_wraparound_positive(self):
        """Test wraparound when going from low to high angle."""
        # 1° - 359° = 2° (wrapping around 0)
        a = 1 * CS_PER_DEG  # 1°
        b = 359 * CS_PER_DEG  # 359°
        result = ephem.difcsn(a, b)
        expected = 2 * CS_PER_DEG  # 2°
        assert result == expected

    def test_difcsn_at_180_boundary(self):
        """Test behavior at exactly 180° difference."""
        # 180° - 0° = 180°
        result = ephem.difcsn(CS180, 0)
        assert result == CS180

    def test_difcsn_just_over_180(self):
        """Test behavior just over 180° difference."""
        # 181° - 0° = 181°
        a = 181 * CS_PER_DEG
        result = ephem.difcsn(a, 0)
        expected = 181 * CS_PER_DEG
        assert result == expected

    def test_difcsn_just_under_180(self):
        """Test behavior just under 180° difference."""
        # 179° - 0° = 179°
        a = 179 * CS_PER_DEG
        result = ephem.difcsn(a, 0)
        expected = 179 * CS_PER_DEG
        assert result == expected


class TestDifcsnNormalization:
    """Tests for input normalization (angles > 360° or < 0)."""

    def test_difcsn_input_over_360(self):
        """Test that inputs over 360° are handled correctly."""
        # (360° + 10°) - 5° = 5° (370° normalized to 10°)
        a = CS360 + 10 * CS_PER_DEG
        b = 5 * CS_PER_DEG
        result = ephem.difcsn(a, b)
        expected = 5 * CS_PER_DEG
        assert result == expected

    def test_difcsn_negative_input(self):
        """Test that negative inputs are handled correctly."""
        # -10° - 5° should work correctly (normalizes to 350° - 5° = 345°)
        a = -10 * CS_PER_DEG
        b = 5 * CS_PER_DEG
        result = ephem.difcsn(a, b)
        expected = 345 * CS_PER_DEG
        assert result == expected

    def test_difcsn_both_negative(self):
        """Test with both negative inputs."""
        # -10° - (-20°) = 10°
        a = -10 * CS_PER_DEG
        b = -20 * CS_PER_DEG
        result = ephem.difcsn(a, b)
        expected = 10 * CS_PER_DEG
        assert result == expected


class TestDifcsnVsSwisseph:
    """Tests comparing difcsn with pyswisseph's swe.difcsn."""

    def test_difcsn_matches_swisseph_zero(self):
        """Test that difcsn matches swisseph for zero inputs."""
        assert ephem.difcsn(0, 0) == swe.difcsn(0, 0)

    def test_difcsn_matches_swisseph_simple(self):
        """Test that difcsn matches swisseph for simple cases."""
        for a_deg in [0, 1, 10, 45, 90, 180, 270, 359]:
            for b_deg in [0, 1, 10, 45, 90, 180, 270, 359]:
                a_cs = a_deg * CS_PER_DEG
                b_cs = b_deg * CS_PER_DEG
                result_lib = ephem.difcsn(a_cs, b_cs)
                result_swe = swe.difcsn(a_cs, b_cs)
                assert result_lib == result_swe, (
                    f"Mismatch for {a_deg}° - {b_deg}°: "
                    f"lib={result_lib}, swe={result_swe}"
                )

    def test_difcsn_matches_swisseph_wraparound(self):
        """Test that difcsn matches swisseph for wraparound cases."""
        test_cases = [
            (1, 359),
            (359, 1),
            (0, 180),
            (180, 0),
            (90, 270),
            (270, 90),
        ]
        for a_deg, b_deg in test_cases:
            a_cs = a_deg * CS_PER_DEG
            b_cs = b_deg * CS_PER_DEG
            result_lib = ephem.difcsn(a_cs, b_cs)
            result_swe = swe.difcsn(a_cs, b_cs)
            assert result_lib == result_swe, (
                f"Mismatch for {a_deg}° - {b_deg}°: lib={result_lib}, swe={result_swe}"
            )


class TestDifcsnEdgeCases:
    """Tests for edge cases."""

    def test_difcsn_full_circle(self):
        """Test that full circle difference is 0."""
        # 360° - 0° should be 0 (360° normalizes to 0°)
        assert ephem.difcsn(CS360, 0) == 0

    def test_difcsn_very_large_values(self):
        """Test with very large centisecond values."""
        # Multiple full circles
        a = 10 * CS360 + 45 * CS_PER_DEG  # 10 full circles + 45°
        b = 3 * CS360 + 30 * CS_PER_DEG  # 3 full circles + 30°
        result = ephem.difcsn(a, b)
        expected = 15 * CS_PER_DEG  # 45° - 30° = 15°
        assert result == expected

    def test_difcsn_fractional_arcseconds(self):
        """Test with sub-arcsecond precision."""
        # 1 centisecond difference
        a = 100  # 1 arcsecond
        b = 99  # 0.99 arcseconds
        result = ephem.difcsn(a, b)
        assert result == 1


class TestDifcsnRandomValues:
    """Tests with random values to ensure comprehensive coverage."""

    @pytest.mark.parametrize(
        "a_deg,b_deg",
        [
            (0, 0),
            (45, 45),
            (90, 90),
            (180, 180),
            (270, 270),
            (10, 20),
            (20, 10),
            (350, 10),
            (10, 350),
            (179, 181),
            (181, 179),
            (1, 359),
            (359, 1),
            (123, 234),
            (234, 123),
        ],
    )
    def test_difcsn_parametrized_vs_swisseph(self, a_deg, b_deg):
        """Parametrized test comparing difcsn with swisseph."""
        a_cs = a_deg * CS_PER_DEG
        b_cs = b_deg * CS_PER_DEG
        result_lib = ephem.difcsn(a_cs, b_cs)
        result_swe = swe.difcsn(a_cs, b_cs)
        assert result_lib == result_swe


class TestDifcsnSymmetry:
    """Tests for mathematical properties."""

    def test_difcsn_not_symmetric(self):
        """Test that difcsn(a, b) != difcsn(b, a) in general."""
        a = 10 * CS_PER_DEG
        b = 20 * CS_PER_DEG
        # difcsn(10, 20) = 350°, difcsn(20, 10) = 10°
        result_ab = ephem.difcsn(a, b)
        result_ba = ephem.difcsn(b, a)
        assert result_ab != result_ba
        assert result_ab + result_ba == CS360

    def test_difcsn_complements_to_360(self):
        """Test that difcsn(a, b) + difcsn(b, a) = 360° for non-equal angles."""
        test_cases = [(10, 20), (1, 359), (90, 270), (45, 180)]
        for a_deg, b_deg in test_cases:
            a_cs = a_deg * CS_PER_DEG
            b_cs = b_deg * CS_PER_DEG
            result_ab = ephem.difcsn(a_cs, b_cs)
            result_ba = ephem.difcsn(b_cs, a_cs)
            assert result_ab + result_ba == CS360, (
                f"Complement check failed for {a_deg}°, {b_deg}°"
            )

    def test_difcsn_always_non_negative(self):
        """Test that difcsn always returns non-negative values."""
        test_cases = [
            (0, 0),
            (10, 20),
            (20, 10),
            (1, 359),
            (359, 1),
            (180, 0),
            (0, 180),
        ]
        for a_deg, b_deg in test_cases:
            a_cs = a_deg * CS_PER_DEG
            b_cs = b_deg * CS_PER_DEG
            result = ephem.difcsn(a_cs, b_cs)
            assert result >= 0, f"Negative result for {a_deg}° - {b_deg}°: {result}"
            assert result < CS360, f"Result >= 360° for {a_deg}° - {b_deg}°: {result}"
