"""
Tests for the _houses_equal_mc function.

Verifies that the function correctly uses the true Ascendant instead of an approximation.
"""

import pytest


class TestHousesEqualMc:
    """Test the _houses_equal_mc function."""

    @pytest.mark.unit
    def test_uses_true_ascendant_not_approximation(self):
        """Verify the function uses the provided Ascendant, not MC+90."""
        from libephemeris.houses import _houses_equal_mc

        # Test with Asc and MC that are NOT 90 degrees apart
        # At non-equatorial latitudes, Asc != MC + 90
        asc = 45.0  # True Ascendant
        mc = 0.0  # Midheaven

        # MC + 90 would be 90.0, but we pass asc = 45.0
        cusps = _houses_equal_mc(asc, mc)

        # First cusp should be the Ascendant (45.0), not MC+90 (90.0)
        assert abs(cusps[1] - 45.0) < 0.0001, f"Cusp 1 should be 45.0, got {cusps[1]}"

    @pytest.mark.unit
    def test_equal_30_degree_divisions(self):
        """Verify cusps are 30 degrees apart."""
        from libephemeris.houses import _houses_equal_mc

        asc = 0.0
        mc = 270.0  # Typical MC for Aries rising

        cusps = _houses_equal_mc(asc, mc)

        # Check each cusp is 30 degrees from the previous
        for i in range(2, 13):
            expected = (asc + (i - 1) * 30.0) % 360.0
            assert abs(cusps[i] - expected) < 0.0001, (
                f"Cusp {i} expected {expected}, got {cusps[i]}"
            )

    @pytest.mark.unit
    def test_returns_13_element_list(self):
        """Verify the function returns a list with 13 elements (0-indexed with 0 unused)."""
        from libephemeris.houses import _houses_equal_mc

        cusps = _houses_equal_mc(0.0, 270.0)

        assert len(cusps) == 13
        assert cusps[0] == 0.0  # Element 0 is unused

    @pytest.mark.unit
    def test_various_ascendant_values(self):
        """Test with various Ascendant values."""
        from libephemeris.houses import _houses_equal_mc

        test_cases = [
            (0.0, 270.0),  # Aries rising
            (90.0, 0.0),  # Cancer rising
            (180.0, 90.0),  # Libra rising
            (270.0, 180.0),  # Capricorn rising
            (123.456, 33.456),  # Arbitrary values
        ]

        for asc, mc in test_cases:
            cusps = _houses_equal_mc(asc, mc)

            # First cusp should always be the Ascendant
            assert abs(cusps[1] - asc) < 0.0001, (
                f"For asc={asc}, cusp 1 should be {asc}, got {cusps[1]}"
            )

            # House 7 should be opposite of House 1
            expected_dsc = (asc + 180.0) % 360.0
            assert abs(cusps[7] - expected_dsc) < 0.0001, (
                f"For asc={asc}, cusp 7 should be {expected_dsc}, got {cusps[7]}"
            )

    @pytest.mark.unit
    def test_mc_parameter_is_accepted(self):
        """Verify the function accepts mc parameter even though it doesn't use it directly."""
        from libephemeris.houses import _houses_equal_mc

        # The mc parameter is kept for API consistency with other house functions
        # It might be used for validation or future enhancements
        asc = 45.0
        mc = 315.0

        # Should not raise any error
        cusps = _houses_equal_mc(asc, mc)
        assert cusps[1] == asc

    @pytest.mark.unit
    def test_wrapping_around_360(self):
        """Verify cusps wrap around 360 degrees correctly."""
        from libephemeris.houses import _houses_equal_mc

        asc = 350.0  # Near end of circle
        mc = 260.0

        cusps = _houses_equal_mc(asc, mc)

        # Cusp 1 = 350
        assert abs(cusps[1] - 350.0) < 0.0001
        # Cusp 2 = 350 + 30 = 380 % 360 = 20
        assert abs(cusps[2] - 20.0) < 0.0001
        # Cusp 3 = 350 + 60 = 410 % 360 = 50
        assert abs(cusps[3] - 50.0) < 0.0001


class TestEqualMcVsOldApproximation:
    """Test that the new implementation differs from the old approximation."""

    @pytest.mark.unit
    def test_differs_from_old_approximation_at_non_equatorial_lat(self):
        """
        At non-equatorial latitudes, the true Asc differs from MC+90.
        This test verifies we now use the true Asc instead of approximation.
        """
        from libephemeris.houses import _houses_equal_mc

        # At latitude 45N, typical difference between true Asc and MC+90
        # can be several degrees
        true_asc = 45.0  # Example: true calculated Ascendant
        mc = 320.0  # This would give old approx asc = 50.0

        cusps = _houses_equal_mc(true_asc, mc)

        old_approx_asc = (mc + 90.0) % 360.0  # 50.0

        # The cusp should be 45.0 (true_asc), not 50.0 (old approximation)
        assert abs(cusps[1] - true_asc) < 0.0001
        assert (
            abs(cusps[1] - old_approx_asc) > 0.1
        )  # Should differ from old approximation
