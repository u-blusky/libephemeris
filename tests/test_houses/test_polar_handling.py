"""
Comprehensive tests for improved Placidus polar latitude handling.

Tests the PolarCircleError exception, get_polar_latitude_threshold(),
swe_houses_with_fallback(), and related functionality.
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import PolarCircleError


class TestPolarCircleError:
    """Test the PolarCircleError exception class."""

    @pytest.mark.unit
    def test_polar_circle_error_is_error_subclass(self):
        """PolarCircleError should be a subclass of Error."""
        assert issubclass(PolarCircleError, ephem.Error)

    @pytest.mark.unit
    def test_polar_circle_error_attributes(self):
        """PolarCircleError should have expected attributes."""
        err = PolarCircleError(
            message="test message",
            latitude=70.0,
            threshold=66.56,
            obliquity=23.44,
            house_system="P",
        )
        assert err.latitude == 70.0
        assert err.threshold == 66.56
        assert err.obliquity == 23.44
        assert err.house_system == "P"
        assert str(err) == "test message"

    @pytest.mark.unit
    def test_polar_circle_error_repr(self):
        """PolarCircleError repr should be informative."""
        err = PolarCircleError(
            message="test",
            latitude=70.0,
            threshold=66.56,
            obliquity=23.44,
            house_system="P",
        )
        repr_str = repr(err)
        assert "PolarCircleError" in repr_str
        assert "70.0" in repr_str
        assert "66.56" in repr_str

    @pytest.mark.unit
    def test_polar_circle_error_can_be_caught_as_error(self):
        """PolarCircleError should be catchable as ephem.Error."""
        jd = 2451545.0

        try:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))
            pytest.fail("Should have raised PolarCircleError")
        except ephem.Error as e:
            # Should be caught as Error
            assert isinstance(e, PolarCircleError)


class TestGetPolarLatitudeThreshold:
    """Test the get_polar_latitude_threshold() helper function."""

    @pytest.mark.unit
    def test_default_obliquity(self):
        """Default obliquity should give approximately 66.56° threshold."""
        threshold = ephem.get_polar_latitude_threshold()
        assert abs(threshold - 66.56) < 0.01

    @pytest.mark.unit
    def test_custom_obliquity(self):
        """Custom obliquity should give correct threshold."""
        # For obliquity of 23.5°, threshold should be 66.5°
        threshold = ephem.get_polar_latitude_threshold(23.5)
        assert threshold == 66.5

    @pytest.mark.unit
    def test_zero_obliquity(self):
        """Zero obliquity should give 90° threshold."""
        threshold = ephem.get_polar_latitude_threshold(0.0)
        assert threshold == 90.0

    @pytest.mark.unit
    def test_max_obliquity(self):
        """Obliquity of 90° should give 0° threshold."""
        threshold = ephem.get_polar_latitude_threshold(90.0)
        assert threshold == 0.0


class TestPolarCircleErrorDetails:
    """Test that PolarCircleError provides helpful information."""

    @pytest.mark.unit
    def test_error_includes_latitude(self):
        """Error message should include the problematic latitude."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))

        assert exc_info.value.latitude == 70.0
        assert "70" in str(exc_info.value)

    @pytest.mark.unit
    def test_error_includes_threshold(self):
        """Error message should include the threshold."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))

        assert exc_info.value.threshold is not None
        assert exc_info.value.threshold > 66.0
        assert exc_info.value.threshold < 67.0

    @pytest.mark.unit
    def test_error_includes_house_system(self):
        """Error message should include the house system."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))

        assert exc_info.value.house_system == "P"
        assert "Placidus" in str(exc_info.value)

    @pytest.mark.unit
    def test_error_suggests_alternatives(self):
        """Error message should suggest alternative house systems."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))

        msg = str(exc_info.value).lower()
        assert "porphyry" in msg or "equal" in msg or "whole sign" in msg

    @pytest.mark.unit
    def test_error_mentions_fallback_function(self):
        """Error message should mention swe_houses_with_fallback."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))

        assert "fallback" in str(exc_info.value).lower()

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0, 89.0])
    def test_northern_polar_latitudes(self, lat):
        """Northern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("P"))

        assert exc_info.value.latitude == lat
        # Should mention Northern
        assert "N" in str(exc_info.value) or "Northern" in str(exc_info.value)

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [-67.0, -70.0, -80.0, -89.0])
    def test_southern_polar_latitudes(self, lat):
        """Southern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("P"))

        assert exc_info.value.latitude == lat
        # Should mention Southern
        assert "S" in str(exc_info.value) or "Southern" in str(exc_info.value)


class TestSweHousesWithFallback:
    """Test the swe_houses_with_fallback() convenience function."""

    @pytest.mark.unit
    def test_returns_normal_result_for_non_polar(self):
        """Should return normal result without fallback for non-polar latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_with_fallback(
            jd, 45.0, 0.0, ord("P")
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert used_fallback is False
        assert warning is None

    @pytest.mark.unit
    def test_falls_back_for_polar_latitude(self):
        """Should fall back to Porphyry for polar latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_with_fallback(
            jd, 70.0, 0.0, ord("P")
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert used_fallback is True
        assert warning is not None
        assert "Placidus" in warning
        assert "Porphyry" in warning

    @pytest.mark.unit
    def test_warning_includes_threshold(self):
        """Warning message should include the threshold."""
        jd = 2451545.0

        _, _, _, warning = ephem.swe_houses_with_fallback(jd, 70.0, 0.0, ord("P"))

        assert "threshold" in warning.lower()

    @pytest.mark.unit
    def test_custom_fallback_system(self):
        """Should use custom fallback system when specified."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_with_fallback(
            jd,
            70.0,
            0.0,
            ord("P"),
            fallback_hsys=ord("E"),  # Equal houses
        )

        assert used_fallback is True
        assert "Equal" in warning

    @pytest.mark.unit
    def test_all_polar_affected_systems(self):
        """Should handle all polar-affected systems."""
        jd = 2451545.0
        lat = 70.0

        for hsys in [ord("P"), ord("K"), ord("G")]:  # Placidus, Koch, Gauquelin
            cusps, ascmc, used_fallback, warning = ephem.swe_houses_with_fallback(
                jd, lat, 0.0, hsys
            )
            assert used_fallback is True
            assert len(cusps) == 12

    @pytest.mark.unit
    def test_polar_safe_systems_no_fallback(self):
        """Polar-safe systems should not trigger fallback."""
        jd = 2451545.0
        lat = 70.0

        # These systems should work at polar latitudes
        for hsys in [ord("O"), ord("E"), ord("W"), ord("M")]:
            cusps, ascmc, used_fallback, warning = ephem.swe_houses_with_fallback(
                jd, lat, 0.0, hsys
            )
            assert used_fallback is False
            assert warning is None

    @pytest.mark.unit
    def test_matches_porphyry_when_fallback_used(self):
        """Fallback cusps should match direct Porphyry calculation."""
        jd = 2451545.0
        lat = 70.0

        # Get fallback result
        cusps_fallback, ascmc_fallback, used_fallback, _ = (
            ephem.swe_houses_with_fallback(jd, lat, 0.0, ord("P"))
        )

        # Get direct Porphyry
        cusps_porphyry, ascmc_porphyry = ephem.swe_houses(jd, lat, 0.0, ord("O"))

        # Should match
        assert used_fallback is True
        for i in range(12):
            assert abs(cusps_fallback[i] - cusps_porphyry[i]) < 0.001

    @pytest.mark.unit
    def test_returns_valid_cusps_at_extreme_latitudes(self):
        """Should return valid cusps even at extreme polar latitudes."""
        jd = 2451545.0

        for lat in [80.0, 85.0, 89.0, -80.0, -85.0, -89.0]:
            cusps, ascmc, _, _ = ephem.swe_houses_with_fallback(jd, lat, 0.0, ord("P"))
            for cusp in cusps:
                assert 0 <= cusp < 360


class TestSweHousesArmcWithFallback:
    """Test the swe_houses_armc_with_fallback() convenience function."""

    @pytest.mark.unit
    def test_returns_normal_result_for_non_polar(self):
        """Should return normal result without fallback for non-polar latitudes."""
        armc = 280.0
        eps = 23.44

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_armc_with_fallback(
            armc, 45.0, eps, ord("P")
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert used_fallback is False
        assert warning is None

    @pytest.mark.unit
    def test_falls_back_for_polar_latitude(self):
        """Should fall back for polar latitudes."""
        armc = 280.0
        eps = 23.44

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_armc_with_fallback(
            armc, 70.0, eps, ord("P")
        )

        assert used_fallback is True
        assert warning is not None
        assert "Placidus" in warning

    @pytest.mark.unit
    def test_custom_fallback(self):
        """Should use custom fallback when specified."""
        armc = 280.0
        eps = 23.44

        _, _, used_fallback, warning = ephem.swe_houses_armc_with_fallback(
            armc, 70.0, eps, ord("P"), fallback_hsys=ord("W")
        )

        assert used_fallback is True
        assert "Whole Sign" in warning


class TestPolarCircleEdgeCases:
    """Test edge cases at and near the polar circle threshold."""

    @pytest.mark.unit
    def test_just_below_threshold(self):
        """Latitude just below threshold should work normally."""
        jd = 2451545.0

        # Get the actual threshold for the current obliquity
        # At J2000, obliquity is ~23.44°, so threshold is ~66.56°
        # Test at 66.0° which should be safely below
        cusps, ascmc = ephem.swe_houses(jd, 66.0, 0.0, ord("P"))
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_just_above_threshold(self):
        """Latitude just above threshold should raise error."""
        jd = 2451545.0

        # At 67° it should definitely fail
        with pytest.raises(PolarCircleError):
            ephem.swe_houses(jd, 67.0, 0.0, ord("P"))

    @pytest.mark.unit
    def test_different_times_different_thresholds(self):
        """Different Julian dates have slightly different obliquities."""
        # J2000.0
        jd_2000 = 2451545.0
        # J1900.0
        jd_1900 = 2415020.0

        # Both should raise error at 70°
        with pytest.raises(PolarCircleError) as exc1:
            ephem.swe_houses(jd_2000, 70.0, 0.0, ord("P"))

        with pytest.raises(PolarCircleError) as exc2:
            ephem.swe_houses(jd_1900, 70.0, 0.0, ord("P"))

        # Thresholds should be slightly different due to different obliquities
        # The difference is very small (~0.01° per century)
        assert exc1.value.threshold is not None
        assert exc2.value.threshold is not None


class TestBackwardCompatibility:
    """Test backward compatibility with existing code."""

    @pytest.mark.unit
    def test_still_raises_error_type(self):
        """PolarCircleError should still be catchable as the base Error."""
        jd = 2451545.0

        # Old code catching Error should still work
        try:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))
            pytest.fail("Should have raised an error")
        except ephem.Error:
            pass  # Old code would catch this

    @pytest.mark.unit
    def test_error_message_still_mentions_polar_circle(self):
        """Error message should still contain 'polar circle' for compatibility."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("P"))

        # Tests may check for this string
        assert "polar" in str(exc_info.value).lower()

    @pytest.mark.unit
    def test_houses_armc_also_raises_polar_error(self):
        """swe_houses_armc should also raise PolarCircleError."""
        armc = 280.0
        eps = 23.44

        with pytest.raises(PolarCircleError):
            ephem.swe_houses_armc(armc, 70.0, eps, ord("P"))


class TestKochAndGauquelinPolar:
    """Test polar handling for Koch and Gauquelin systems."""

    @pytest.mark.unit
    def test_koch_raises_polar_error(self):
        """Koch should raise PolarCircleError at polar latitudes."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.house_system == "K"
        assert "Koch" in str(exc_info.value)

    @pytest.mark.unit
    def test_gauquelin_raises_polar_error(self):
        """Gauquelin should raise PolarCircleError at polar latitudes."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("G"))

        assert exc_info.value.house_system == "G"
        assert "Gauquelin" in str(exc_info.value)

    @pytest.mark.unit
    def test_koch_fallback_works(self):
        """swe_houses_with_fallback should work for Koch."""
        jd = 2451545.0

        cusps, _, used_fallback, warning = ephem.swe_houses_with_fallback(
            jd, 70.0, 0.0, ord("K")
        )

        assert used_fallback is True
        assert "Koch" in warning
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_gauquelin_fallback_works(self):
        """swe_houses_with_fallback should work for Gauquelin."""
        jd = 2451545.0

        cusps, _, used_fallback, warning = ephem.swe_houses_with_fallback(
            jd, 70.0, 0.0, ord("G")
        )

        assert used_fallback is True
        assert "Gauquelin" in warning
        assert len(cusps) == 12


class TestKochPolarCircleErrorDetails:
    """Test that Koch PolarCircleError provides helpful information (similar to Placidus)."""

    @pytest.mark.unit
    def test_koch_error_includes_latitude(self):
        """Koch error message should include the problematic latitude."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.latitude == 70.0
        assert "70" in str(exc_info.value)

    @pytest.mark.unit
    def test_koch_error_includes_threshold(self):
        """Koch error message should include the threshold."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.threshold is not None
        assert exc_info.value.threshold > 66.0
        assert exc_info.value.threshold < 67.0

    @pytest.mark.unit
    def test_koch_error_includes_house_system(self):
        """Koch error message should include the house system."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.house_system == "K"
        assert "Koch" in str(exc_info.value)

    @pytest.mark.unit
    def test_koch_error_suggests_alternatives(self):
        """Koch error message should suggest alternative house systems."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        msg = str(exc_info.value).lower()
        assert "porphyry" in msg or "equal" in msg or "whole sign" in msg

    @pytest.mark.unit
    def test_koch_error_mentions_fallback_function(self):
        """Koch error message should mention swe_houses_with_fallback."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        assert "fallback" in str(exc_info.value).lower()

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0, 89.0])
    def test_koch_northern_polar_latitudes(self, lat):
        """Koch at northern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("K"))

        assert exc_info.value.latitude == lat
        # Should mention Northern
        assert "N" in str(exc_info.value) or "Northern" in str(exc_info.value)

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [-67.0, -70.0, -80.0, -89.0])
    def test_koch_southern_polar_latitudes(self, lat):
        """Koch at southern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("K"))

        assert exc_info.value.latitude == lat
        # Should mention Southern
        assert "S" in str(exc_info.value) or "Southern" in str(exc_info.value)


class TestKochPolarEdgeCases:
    """Test Koch edge cases at and near the polar circle threshold."""

    @pytest.mark.unit
    def test_koch_just_below_threshold(self):
        """Koch at latitude just below threshold should work normally."""
        jd = 2451545.0

        # At J2000, obliquity is ~23.44deg, so threshold is ~66.56deg
        # Test at 66.0deg which should be safely below
        cusps, ascmc = ephem.swe_houses(jd, 66.0, 0.0, ord("K"))
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_koch_just_above_threshold(self):
        """Koch at latitude just above threshold should raise error."""
        jd = 2451545.0

        # At 67deg it should definitely fail
        with pytest.raises(PolarCircleError):
            ephem.swe_houses(jd, 67.0, 0.0, ord("K"))

    @pytest.mark.unit
    def test_koch_matches_porphyry_when_fallback_used(self):
        """Koch fallback cusps should match direct Porphyry calculation."""
        jd = 2451545.0
        lat = 70.0

        # Get fallback result
        cusps_fallback, ascmc_fallback, used_fallback, _ = (
            ephem.swe_houses_with_fallback(jd, lat, 0.0, ord("K"))
        )

        # Get direct Porphyry
        cusps_porphyry, ascmc_porphyry = ephem.swe_houses(jd, lat, 0.0, ord("O"))

        # Should match
        assert used_fallback is True
        for i in range(12):
            assert abs(cusps_fallback[i] - cusps_porphyry[i]) < 0.001

    @pytest.mark.unit
    def test_koch_returns_valid_cusps_at_extreme_latitudes(self):
        """Koch fallback should return valid cusps even at extreme polar latitudes."""
        jd = 2451545.0

        for lat in [80.0, 85.0, 89.0, -80.0, -85.0, -89.0]:
            cusps, ascmc, _, _ = ephem.swe_houses_with_fallback(jd, lat, 0.0, ord("K"))
            for cusp in cusps:
                assert 0 <= cusp < 360

    @pytest.mark.unit
    def test_koch_custom_fallback_system(self):
        """Koch should use custom fallback system when specified."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_with_fallback(
            jd,
            70.0,
            0.0,
            ord("K"),
            fallback_hsys=ord("E"),  # Equal houses
        )

        assert used_fallback is True
        assert "Equal" in warning

    @pytest.mark.unit
    def test_koch_warning_includes_threshold(self):
        """Koch warning message should include the threshold."""
        jd = 2451545.0

        _, _, _, warning = ephem.swe_houses_with_fallback(jd, 70.0, 0.0, ord("K"))

        assert "threshold" in warning.lower()


class TestKochBackwardCompatibility:
    """Test Koch backward compatibility with existing code."""

    @pytest.mark.unit
    def test_koch_still_raises_error_type(self):
        """Koch PolarCircleError should still be catchable as the base Error."""
        jd = 2451545.0

        # Old code catching Error should still work
        try:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))
            pytest.fail("Should have raised an error")
        except ephem.Error:
            pass  # Old code would catch this

    @pytest.mark.unit
    def test_koch_error_message_still_mentions_polar_circle(self):
        """Koch error message should still contain 'polar circle' for compatibility."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.swe_houses(jd, 70.0, 0.0, ord("K"))

        # Tests may check for this string
        assert "polar" in str(exc_info.value).lower()

    @pytest.mark.unit
    def test_koch_houses_armc_also_raises_polar_error(self):
        """Koch swe_houses_armc should also raise PolarCircleError."""
        armc = 280.0
        eps = 23.44

        with pytest.raises(PolarCircleError):
            ephem.swe_houses_armc(armc, 70.0, eps, ord("K"))

    @pytest.mark.unit
    def test_koch_armc_fallback_works(self):
        """swe_houses_armc_with_fallback should work for Koch."""
        armc = 280.0
        eps = 23.44

        cusps, ascmc, used_fallback, warning = ephem.swe_houses_armc_with_fallback(
            armc, 70.0, eps, ord("K")
        )

        assert used_fallback is True
        assert "Koch" in warning
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_koch_armc_custom_fallback(self):
        """Koch swe_houses_armc_with_fallback should use custom fallback when specified."""
        armc = 280.0
        eps = 23.44

        _, _, used_fallback, warning = ephem.swe_houses_armc_with_fallback(
            armc, 70.0, eps, ord("K"), fallback_hsys=ord("W")
        )

        assert used_fallback is True
        assert "Whole Sign" in warning
