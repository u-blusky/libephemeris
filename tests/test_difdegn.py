"""
Tests for the difdegn angular difference function.

Tests verify that difdegn correctly calculates the angular difference
(p1 - p2) normalized to the range [0, 360), always returning a positive value.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestDifdegnBasic:
    """Basic functionality tests for difdegn."""

    def test_difdegn_exported(self):
        """Test that difdegn is exported from the package."""
        assert hasattr(ephem, "difdegn")
        assert callable(ephem.difdegn)

    def test_difdegn_same_angle(self):
        """Test distance between same angles is zero."""
        assert ephem.difdegn(0, 0) == 0.0
        assert ephem.difdegn(180, 180) == 0.0
        assert ephem.difdegn(359.5, 359.5) == 0.0

    def test_difdegn_basic_difference(self):
        """Test basic angle differences."""
        # 20 - 10 = 10
        assert ephem.difdegn(20, 10) == pytest.approx(10.0)
        # 10 - 20 = -10 -> 350
        assert ephem.difdegn(10, 20) == pytest.approx(350.0)

    def test_difdegn_returns_float(self):
        """Test that difdegn returns a float."""
        result = ephem.difdegn(45, 90)
        assert isinstance(result, float)


class TestDifdegnAlwaysPositive:
    """Tests to verify difdegn always returns positive values in [0, 360)."""

    def test_difdegn_always_positive_simple(self):
        """Test that result is always positive for simple cases."""
        # difdegn(10, 20) = (10 - 20) % 360 = 350
        assert ephem.difdegn(10, 20) == pytest.approx(350.0)
        assert ephem.difdegn(10, 20) >= 0

        # difdegn(20, 10) = (20 - 10) % 360 = 10
        assert ephem.difdegn(20, 10) == pytest.approx(10.0)
        assert ephem.difdegn(20, 10) >= 0

    def test_difdegn_always_positive_wraparound(self):
        """Test that result is always positive for wraparound cases."""
        # (350 - 10) % 360 = 340
        assert ephem.difdegn(350, 10) == pytest.approx(340.0)
        # (10 - 350) % 360 = 20
        assert ephem.difdegn(10, 350) == pytest.approx(20.0)

    def test_difdegn_180_degree_difference(self):
        """Test 180 degree differences."""
        assert ephem.difdegn(0, 180) == pytest.approx(180.0)
        assert ephem.difdegn(180, 0) == pytest.approx(180.0)
        assert ephem.difdegn(90, 270) == pytest.approx(180.0)


class TestDifdegnWraparound:
    """Tests for wraparound handling."""

    def test_difdegn_across_zero(self):
        """Test angles that wrap around 0/360."""
        # (350 - 10) % 360 = 340
        assert ephem.difdegn(350, 10) == pytest.approx(340.0)
        # (10 - 350) % 360 = 20
        assert ephem.difdegn(10, 350) == pytest.approx(20.0)
        # (359 - 1) % 360 = 358
        assert ephem.difdegn(359, 1) == pytest.approx(358.0)
        # (1 - 359) % 360 = 2
        assert ephem.difdegn(1, 359) == pytest.approx(2.0)

    def test_difdegn_near_180(self):
        """Test angles near 180 degree separation."""
        # (0 - 179) % 360 = 181
        assert ephem.difdegn(0, 179) == pytest.approx(181.0)
        # (0 - 181) % 360 = 179
        assert ephem.difdegn(0, 181) == pytest.approx(179.0)

    def test_difdegn_normalized_angles(self):
        """Test with angles beyond 0-360 range."""
        # (370 - 10) % 360 = 0
        assert ephem.difdegn(370, 10) == pytest.approx(0.0)
        # (-10 - 10) % 360 = 340
        assert ephem.difdegn(-10, 10) == pytest.approx(340.0)
        # (720 - 0) % 360 = 0
        assert ephem.difdegn(720, 0) == pytest.approx(0.0)


class TestDifdegnVsDifdeg2n:
    """Tests comparing difdegn with difdeg2n."""

    @pytest.mark.parametrize(
        "p1,p2",
        [
            (10, 20),
            (20, 10),
            (350, 10),
            (10, 350),
            (0, 180),
            (180, 0),
            (90, 270),
            (270, 90),
            (0, 0),
            (45, 45),
            (0, 90),
            (90, 0),
            (359, 1),
            (1, 359),
        ],
    )
    def test_difdegn_relationship_to_difdeg2n(self, p1, p2):
        """Test that difdegn is the positive normalization of the difference."""
        result_n = ephem.difdegn(p1, p2)
        result_2n = ephem.difdeg2n(p1, p2)

        # difdegn should be in [0, 360)
        assert 0 <= result_n < 360, (
            f"difdegn({p1}, {p2}) = {result_n} is outside [0, 360)"
        )

        # difdegn and difdeg2n should differ by 360 if difdeg2n is negative
        if result_2n < 0:
            assert result_n == pytest.approx(result_2n + 360), (
                f"difdegn({p1}, {p2}) = {result_n} != difdeg2n({p1}, {p2}) + 360 = {result_2n + 360}"
            )
        else:
            assert result_n == pytest.approx(result_2n), (
                f"difdegn({p1}, {p2}) = {result_n} != difdeg2n({p1}, {p2}) = {result_2n}"
            )


class TestDifdegnVsSwisseph:
    """Comparison tests with pyswisseph's swe.difdegn()."""

    @pytest.mark.parametrize(
        "p1,p2",
        [
            (0, 0),
            (10, 20),
            (20, 10),
            (350, 10),
            (10, 350),
            (0, 180),
            (180, 0),
            (90, 270),
            (45, 315),
            (0, 359),
            (359, 0),
            (179, 181),
            (181, 179),
            (-10, 10),
            (370, 10),
            (10, 370),
            (-350, 350),
            (720, 360),
        ],
    )
    def test_difdegn_matches_swisseph(self, p1, p2):
        """Test that difdegn matches pyswisseph's swe.difdegn()."""
        result_lib = ephem.difdegn(p1, p2)
        result_swe = swe.difdegn(p1, p2)

        assert result_lib == pytest.approx(result_swe, abs=1e-10), (
            f"Mismatch for ({p1}, {p2}): lib={result_lib}, swe={result_swe}"
        )


class TestDifdegnEdgeCases:
    """Edge case tests for difdegn."""

    def test_difdegn_exactly_180(self):
        """Test exactly 180 degree separation."""
        assert ephem.difdegn(0, 180) == pytest.approx(180.0)
        assert ephem.difdegn(180, 360) == pytest.approx(180.0)
        assert ephem.difdegn(90, 270) == pytest.approx(180.0)

    def test_difdegn_very_small_difference(self):
        """Test very small angular differences."""
        result = ephem.difdegn(0, 1e-10)
        expected = swe.difdegn(0, 1e-10)
        assert result == pytest.approx(expected, abs=1e-10)

    def test_difdegn_range_always_0_to_360(self):
        """Test that result is always in range [0, 360)."""
        test_pairs = [
            (0, 0),
            (0, 90),
            (0, 180),
            (0, 270),
            (0, 359),
            (45, 135),
            (45, 225),
            (45, 315),
            (90, 270),
            (180, 0),
            (270, 90),
            (350, 10),
            (10, 350),
        ]
        for p1, p2 in test_pairs:
            result = ephem.difdegn(p1, p2)
            assert 0 <= result < 360, (
                f"difdegn({p1}, {p2}) = {result} is outside [0, 360)"
            )


class TestDifdegnRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles(self, random_longitudes):
        """Test with random angles."""
        lons = random_longitudes(100)

        for i in range(0, len(lons) - 1, 2):
            p1 = lons[i]
            p2 = lons[i + 1]

            result_lib = ephem.difdegn(p1, p2)
            result_swe = swe.difdegn(p1, p2)

            assert result_lib == pytest.approx(result_swe, abs=1e-10), (
                f"Mismatch for ({p1}, {p2})"
            )
            assert 0 <= result_lib < 360, f"Result outside [0, 360) for ({p1}, {p2})"

    def test_consistency_with_modulo(self, random_longitudes):
        """Test that difdegn(a, b) == (a - b) % 360."""
        lons = random_longitudes(50)

        for i in range(0, len(lons) - 1, 2):
            p1 = lons[i]
            p2 = lons[i + 1]

            result = ephem.difdegn(p1, p2)
            expected = (p1 - p2) % 360.0

            assert result == pytest.approx(expected, abs=1e-10), (
                f"difdegn({p1}, {p2}) = {result} != ({p1} - {p2}) % 360 = {expected}"
            )
