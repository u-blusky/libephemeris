"""
Tests for rise_trans_true_hor function in libephemeris.

Tests the calculation of rise, set, and transit times with custom horizon altitudes.
This is useful for locations with mountains or buildings that occlude the horizon.
"""

import pytest

from libephemeris import (
    julday,
    revjul,
    rise_trans,
    rise_trans_true_hor,
    swe_rise_trans_true_hor,
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_CALC_RISE,
    SE_CALC_SET,
    SE_CALC_MTRANSIT,
    SE_BIT_DISC_CENTER,
    SE_BIT_NO_REFRACTION,
)


class TestRiseTransTrueHorBasic:
    """Basic tests for rise_trans_true_hor function."""

    def test_zero_horizon_matches_rise_trans(self):
        """Test that horizon_altitude=0 matches regular rise_trans."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Regular rise_trans
        jd_rise_std, flag_std = rise_trans(
            jd_start, SE_SUN, lat, lon, rsmi=SE_CALC_RISE
        )

        # rise_trans_true_hor with zero horizon
        jd_rise_hor, flag_hor = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )

        # Should match to within a few seconds (0.0001 JD = ~8.6 seconds)
        assert abs(jd_rise_std - jd_rise_hor) < 0.0001
        assert flag_std == flag_hor

    def test_swe_alias_works(self):
        """Test that swe_rise_trans_true_hor is an alias for rise_trans_true_hor."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964

        result1 = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=5.0, rsmi=SE_CALC_RISE
        )
        result2 = swe_rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=5.0, rsmi=SE_CALC_RISE
        )

        assert result1 == result2

    def test_positive_horizon_delays_sunrise(self):
        """Test that positive horizon_altitude (mountains) delays sunrise."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Standard sunrise (horizon at 0)
        jd_rise_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )

        # Sunrise with mountains at 10 degrees
        jd_rise_mountain, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=10.0, rsmi=SE_CALC_RISE
        )

        # Sun needs to rise higher before it appears, so sunrise is later
        assert jd_rise_mountain > jd_rise_std

        # Difference should be significant (roughly 30-60 minutes for 10 degrees)
        diff_minutes = (jd_rise_mountain - jd_rise_std) * 24 * 60
        assert 20 < diff_minutes < 120

    def test_positive_horizon_advances_sunset(self):
        """Test that positive horizon_altitude makes sunset earlier."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Standard sunset (horizon at 0)
        jd_set_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_SET
        )

        # Sunset with mountains at 10 degrees
        jd_set_mountain, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=10.0, rsmi=SE_CALC_SET
        )

        # Sun sets behind mountain earlier
        assert jd_set_mountain < jd_set_std

        # Difference should be significant
        diff_minutes = (jd_set_std - jd_set_mountain) * 24 * 60
        assert 20 < diff_minutes < 120

    def test_negative_horizon_advances_sunrise(self):
        """Test that negative horizon_altitude (observer on mountain) advances sunrise."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Standard sunrise (horizon at 0)
        jd_rise_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )

        # Sunrise with depressed horizon (observer on high mountain)
        jd_rise_elevated, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=-5.0, rsmi=SE_CALC_RISE
        )

        # Observer sees sunrise earlier when horizon is depressed
        assert jd_rise_elevated < jd_rise_std

    def test_transit_unaffected_by_horizon(self):
        """Test that transit times are not affected by horizon altitude."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Transit with standard horizon
        jd_transit_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_MTRANSIT
        )

        # Transit with elevated horizon
        jd_transit_mountain, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=10.0, rsmi=SE_CALC_MTRANSIT
        )

        # Transit time should be identical (within computational precision)
        assert abs(jd_transit_std - jd_transit_mountain) < 0.0001


class TestRiseTransTrueHorMoon:
    """Tests for Moon rise/set with custom horizon."""

    def test_moonrise_with_horizon(self):
        """Test moonrise with elevated horizon."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 51.5074, -0.1278  # London

        jd_rise_std, _ = rise_trans_true_hor(
            jd_start, SE_MOON, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )

        jd_rise_hor, flag = rise_trans_true_hor(
            jd_start, SE_MOON, lat, lon, horizon_altitude=5.0, rsmi=SE_CALC_RISE
        )

        assert jd_rise_hor > jd_start
        assert jd_rise_hor > jd_rise_std  # Delayed by horizon
        assert flag == SE_CALC_RISE


class TestRiseTransTrueHorPlanets:
    """Tests for planet rise/set with custom horizon."""

    def test_mars_rise_with_horizon(self):
        """Test Mars rise with elevated horizon."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 40.7128, -74.0060  # New York

        jd_rise_std, _ = rise_trans_true_hor(
            jd_start, SE_MARS, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )

        jd_rise_hor, flag = rise_trans_true_hor(
            jd_start, SE_MARS, lat, lon, horizon_altitude=8.0, rsmi=SE_CALC_RISE
        )

        assert jd_rise_hor > jd_start
        assert jd_rise_hor > jd_rise_std
        assert flag == SE_CALC_RISE


class TestRiseTransTrueHorCircumpolar:
    """Tests for circumpolar conditions with custom horizon."""

    def test_circumpolar_with_high_horizon(self):
        """Test that very high horizon can make object circumpolar."""
        # At a moderate latitude, with a very high horizon, the sun might never rise
        jd_start = julday(2024, 12, 21, 0)  # Winter solstice
        lat, lon = 55.0, 10.0  # Denmark

        # With standard horizon, sun rises and sets
        jd_rise_std, flag_std = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )
        assert flag_std == SE_CALC_RISE  # Sun does rise

        # With very high horizon (unrealistic, but tests the logic)
        jd_rise_high, flag_high = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=15.0, rsmi=SE_CALC_RISE
        )

        # Sun might not rise above 15 degrees at this latitude in winter
        # Either it's circumpolar below horizon (-2) or rises very late
        # This depends on the sun's maximum altitude that day
        if flag_high == -2:
            assert jd_rise_high == 0.0
        else:
            # If it does rise, it should be much later than standard
            assert jd_rise_high > jd_rise_std


class TestRiseTransTrueHorFlags:
    """Tests for flag combinations with custom horizon."""

    def test_disc_center_with_horizon(self):
        """Test SE_BIT_DISC_CENTER flag with custom horizon."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome
        horizon = 5.0

        # With upper limb (default)
        jd_rise_limb, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=horizon, rsmi=SE_CALC_RISE
        )

        # With disc center
        jd_rise_center, _ = rise_trans_true_hor(
            jd_start,
            SE_SUN,
            lat,
            lon,
            horizon_altitude=horizon,
            rsmi=SE_CALC_RISE | SE_BIT_DISC_CENTER,
        )

        # Center rise should be later than upper limb rise
        assert jd_rise_center > jd_rise_limb

    def test_no_refraction_with_horizon(self):
        """Test SE_BIT_NO_REFRACTION flag with custom horizon."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome
        horizon = 5.0

        # With refraction (default)
        jd_rise_refr, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=horizon, rsmi=SE_CALC_RISE
        )

        # Without refraction
        jd_rise_no_refr, _ = rise_trans_true_hor(
            jd_start,
            SE_SUN,
            lat,
            lon,
            horizon_altitude=horizon,
            rsmi=SE_CALC_RISE | SE_BIT_NO_REFRACTION,
        )

        # Without refraction, rise should be later
        assert jd_rise_no_refr > jd_rise_refr


class TestRiseTransTrueHorLocations:
    """Tests for different geographic locations with custom horizon."""

    def test_equator_with_horizon(self):
        """Test rise/set at equator with custom horizon."""
        jd_start = julday(2024, 3, 20, 0)  # Near equinox
        lat, lon = 0.0, 0.0  # Equator
        horizon = 5.0

        jd_rise, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=horizon, rsmi=SE_CALC_RISE
        )
        jd_set, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=horizon, rsmi=SE_CALC_SET
        )

        # Both events should occur
        assert jd_rise > jd_start
        assert jd_set > jd_rise

        # Day length with elevated horizon should be shorter
        jd_rise_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )
        jd_set_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_SET
        )

        day_length_std = (jd_set_std - jd_rise_std) * 24
        day_length_hor = (jd_set - jd_rise) * 24

        assert day_length_hor < day_length_std

    def test_southern_hemisphere_with_horizon(self):
        """Test rise/set in southern hemisphere with custom horizon."""
        jd_start = julday(2024, 12, 21, 0)  # Southern summer
        lat, lon = -33.8688, 151.2093  # Sydney
        horizon = 3.0

        jd_rise, flag_rise = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=horizon, rsmi=SE_CALC_RISE
        )

        assert jd_rise > jd_start
        assert flag_rise == SE_CALC_RISE


class TestRiseTransTrueHorErrors:
    """Tests for error handling."""

    def test_invalid_planet_raises_error(self):
        """Test that invalid planet ID raises ValueError."""
        jd_start = julday(2024, 6, 21, 0)

        with pytest.raises(ValueError, match="illegal planet number"):
            rise_trans_true_hor(
                jd_start, 9999, 41.9, 12.5, horizon_altitude=5.0, rsmi=SE_CALC_RISE
            )

    def test_invalid_rsmi_raises_error(self):
        """Test that invalid rsmi raises ValueError."""
        jd_start = julday(2024, 6, 21, 0)

        with pytest.raises(ValueError, match="Invalid event type"):
            rise_trans_true_hor(
                jd_start, SE_SUN, 41.9, 12.5, horizon_altitude=5.0, rsmi=0
            )

        with pytest.raises(ValueError, match="Invalid event type"):
            rise_trans_true_hor(
                jd_start, SE_SUN, 41.9, 12.5, horizon_altitude=5.0, rsmi=16
            )


class TestRiseTransTrueHorVariousAngles:
    """Tests for various horizon angles."""

    def test_small_horizon_adjustment(self):
        """Test small horizon adjustment (typical building)."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Typical building might add 2-3 degrees to horizon
        jd_rise_std, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )
        jd_rise_bldg, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=3.0, rsmi=SE_CALC_RISE
        )

        # Should be slightly later
        assert jd_rise_bldg > jd_rise_std

        # But not by too much (few minutes)
        diff_minutes = (jd_rise_bldg - jd_rise_std) * 24 * 60
        assert 5 < diff_minutes < 30

    def test_moderate_horizon(self):
        """Test moderate horizon (hills)."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        jd_rise, flag = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=5.0, rsmi=SE_CALC_RISE
        )

        assert flag == SE_CALC_RISE
        assert jd_rise > jd_start

        # Verify the date
        year, month, day, hour = revjul(jd_rise)
        assert year == 2024
        assert month == 6
        assert day == 21

    def test_significant_horizon(self):
        """Test significant horizon (mountains)."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 46.0, 7.0  # Swiss Alps region

        # With flat horizon
        jd_rise_flat, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=0.0, rsmi=SE_CALC_RISE
        )

        # With mountain horizon (10 degrees)
        jd_rise_mountain, _ = rise_trans_true_hor(
            jd_start, SE_SUN, lat, lon, horizon_altitude=10.0, rsmi=SE_CALC_RISE
        )

        # Significant delay expected
        diff_minutes = (jd_rise_mountain - jd_rise_flat) * 24 * 60
        assert diff_minutes > 30
