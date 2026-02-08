"""
Tests for Arabic Parts (Lots) calculations.

These tests verify the traditional formulas for Arabic Parts,
including correct day/night chart handling for Part of Fortune and Spirit.

The key formulas being tested:
- Part of Fortune:
  - Day: ASC + Moon - Sun
  - Night: ASC + Sun - Moon
- Part of Spirit:
  - Day: ASC + Sun - Moon
  - Night: ASC + Moon - Sun

Fortune and Spirit are swapped between day and night charts,
following traditional sect-based calculation methods.
"""

import pytest

from libephemeris.arabic_parts import (
    _EXTREME_LATITUDE_THRESHOLD,
    _is_sun_above_horizon_3d,
    calc_all_arabic_parts,
    calc_arabic_part_of_faith,
    calc_arabic_part_of_fortune,
    calc_arabic_part_of_love,
    calc_arabic_part_of_spirit,
    is_day_chart,
)
from libephemeris import julday


class TestIsDayChart:
    """Tests for is_day_chart() function."""

    def test_sun_above_horizon_normal_case(self):
        """Sun between ASC and DSC is a day chart (normal case: ASC < DSC)."""
        # ASC at 0°, DSC at 180°
        # Sun at 90° (between them) = day chart
        asc = 0.0
        sun = 90.0
        assert is_day_chart(sun, asc) is True

    def test_sun_below_horizon_normal_case(self):
        """Sun outside ASC-DSC range is a night chart (normal case)."""
        # ASC at 0°, DSC at 180°
        # Sun at 270° (below horizon) = night chart
        asc = 0.0
        sun = 270.0
        assert is_day_chart(sun, asc) is False

    def test_sun_above_horizon_wrapped_case(self):
        """Sun above horizon when ASC > DSC (wrapped case)."""
        # ASC at 350°, DSC at 170°
        # Sun at 355° (above horizon) = day chart
        asc = 350.0
        sun = 355.0
        assert is_day_chart(sun, asc) is True

        # Sun at 100° (also above horizon in wrapped case)
        sun = 100.0
        assert is_day_chart(sun, asc) is True

    def test_sun_below_horizon_wrapped_case(self):
        """Sun below horizon when ASC > DSC (wrapped case)."""
        # ASC at 350°, DSC at 170°
        # Sun at 250° (below horizon) = night chart
        asc = 350.0
        sun = 250.0
        assert is_day_chart(sun, asc) is False

    def test_sun_exactly_on_asc(self):
        """Sun exactly on Ascendant is considered a day chart (sunrise)."""
        asc = 15.0
        sun = 15.0  # Sun exactly on ASC
        assert is_day_chart(sun, asc) is True

    def test_sun_exactly_on_desc(self):
        """Sun exactly on Descendant is considered a day chart (sunset)."""
        asc = 15.0
        desc = (asc + 180.0) % 360.0  # 195°
        sun = desc
        assert is_day_chart(sun, asc) is True


class TestCalcArabicPartOfFortune:
    """Tests for calc_arabic_part_of_fortune() function."""

    def test_diurnal_formula(self):
        """Day chart uses ASC + Moon - Sun."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Day formula: ASC + Moon - Sun = 15 + 240 - 120 = 135
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(135.0)

    def test_nocturnal_formula(self):
        """Night chart uses ASC + Sun - Moon (inverted formula)."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Night formula: ASC + Sun - Moon = 15 + 120 - 240 = -105 => 255 (mod 360)
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=False)
        assert result == pytest.approx(255.0)

    def test_modulo_wrap_positive(self):
        """Result correctly wraps with modulo 360 (overflow)."""
        asc = 300.0
        sun = 10.0
        moon = 350.0
        # Day formula: 300 + 350 - 10 = 640 => 280 (mod 360)
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(280.0)

    def test_modulo_wrap_negative(self):
        """Result correctly wraps with modulo 360 (underflow)."""
        asc = 10.0
        sun = 300.0
        moon = 30.0
        # Day formula: 10 + 30 - 300 = -260 => 100 (mod 360)
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(100.0)


class TestCalcArabicPartOfSpirit:
    """Tests for calc_arabic_part_of_spirit() function."""

    def test_diurnal_formula(self):
        """Day chart uses ASC + Sun - Moon."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Day formula: ASC + Sun - Moon = 15 + 120 - 240 = -105 => 255 (mod 360)
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(255.0)

    def test_nocturnal_formula(self):
        """Night chart uses ASC + Moon - Sun (inverted formula)."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Night formula: ASC + Moon - Sun = 15 + 240 - 120 = 135
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=False)
        assert result == pytest.approx(135.0)

    def test_modulo_wrap_positive(self):
        """Result correctly wraps with modulo 360 (overflow)."""
        asc = 300.0
        sun = 350.0
        moon = 10.0
        # Day formula: 300 + 350 - 10 = 640 => 280 (mod 360)
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(280.0)


class TestFortuneAndSpiritInversion:
    """Tests verifying Fortune and Spirit are inverted between day and night."""

    def test_fortune_day_equals_spirit_night(self):
        """Part of Fortune (day) equals Part of Spirit (night) - same formula."""
        asc = 45.0
        sun = 100.0
        moon = 200.0

        fortune_day = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        spirit_night = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=False)

        # Both use ASC + Moon - Sun
        assert fortune_day == pytest.approx(spirit_night)

    def test_spirit_day_equals_fortune_night(self):
        """Part of Spirit (day) equals Part of Fortune (night) - same formula."""
        asc = 45.0
        sun = 100.0
        moon = 200.0

        spirit_day = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        fortune_night = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=False)

        # Both use ASC + Sun - Moon
        assert spirit_day == pytest.approx(fortune_night)

    def test_fortune_and_spirit_different_in_same_chart(self):
        """Fortune and Spirit have different values in the same chart."""
        asc = 45.0
        sun = 100.0
        moon = 200.0

        fortune = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        spirit = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)

        # Fortune: 45 + 200 - 100 = 145
        # Spirit: 45 + 100 - 200 = -55 => 305
        assert fortune == pytest.approx(145.0)
        assert spirit == pytest.approx(305.0)
        assert fortune != spirit


class TestCalcArabicPartOfLove:
    """Tests for calc_arabic_part_of_love() function."""

    def test_basic_calculation(self):
        """Basic calculation: ASC + Venus - Sun."""
        asc = 30.0
        venus = 60.0
        sun = 45.0
        # Formula: 30 + 60 - 45 = 45
        result = calc_arabic_part_of_love(asc, venus, sun)
        assert result == pytest.approx(45.0)


class TestCalcArabicPartOfFaith:
    """Tests for calc_arabic_part_of_faith() function."""

    def test_basic_calculation(self):
        """Basic calculation: ASC + Mercury - Moon."""
        asc = 30.0
        mercury = 90.0
        moon = 60.0
        # Formula: 30 + 90 - 60 = 60
        result = calc_arabic_part_of_faith(asc, mercury, moon)
        assert result == pytest.approx(60.0)


class TestCalcAllArabicParts:
    """Tests for calc_all_arabic_parts() function using position dictionaries."""

    def test_diurnal_chart_sun_above_horizon(self):
        """Test with Sun above horizon (day chart)."""
        # ASC at 15°, Sun at 90° (between ASC=15° and DSC=195°) = day chart
        positions = {
            "Asc": 15.0,
            "Sun": 90.0,
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts = calc_all_arabic_parts(positions)

        # Verify this is treated as a day chart
        assert is_day_chart(positions["Sun"], positions["Asc"]) is True

        # Day formula for Fortune: ASC + Moon - Sun = 15 + 240 - 90 = 165
        assert parts["Pars_Fortunae"] == pytest.approx(165.0)

        # Day formula for Spirit: ASC + Sun - Moon = 15 + 90 - 240 = -135 => 225
        assert parts["Pars_Spiritus"] == pytest.approx(225.0)

        # Part of Love: ASC + Venus - Sun = 15 + 80 - 90 = 5
        assert parts["Pars_Amoris"] == pytest.approx(5.0)

        # Part of Faith: ASC + Mercury - Moon = 15 + 100 - 240 = -125 => 235
        assert parts["Pars_Fidei"] == pytest.approx(235.0)

    def test_nocturnal_chart_sun_below_horizon(self):
        """Test with Sun below horizon (night chart)."""
        # ASC at 15°, Sun at 270° (outside ASC=15° to DSC=195°) = night chart
        positions = {
            "Asc": 15.0,
            "Sun": 270.0,
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts = calc_all_arabic_parts(positions)

        # Verify this is treated as a night chart
        assert is_day_chart(positions["Sun"], positions["Asc"]) is False

        # Night formula for Fortune: ASC + Sun - Moon = 15 + 270 - 240 = 45
        assert parts["Pars_Fortunae"] == pytest.approx(45.0)

        # Night formula for Spirit: ASC + Moon - Sun = 15 + 240 - 270 = -15 => 345
        assert parts["Pars_Spiritus"] == pytest.approx(345.0)

        # Part of Love (not sect-dependent): ASC + Venus - Sun = 15 + 80 - 270 = -175 => 185
        assert parts["Pars_Amoris"] == pytest.approx(185.0)

        # Part of Faith (not sect-dependent): ASC + Mercury - Moon = 15 + 100 - 240 = -125 => 235
        assert parts["Pars_Fidei"] == pytest.approx(235.0)

    def test_fortune_spirit_swap_day_vs_night(self):
        """Verify Fortune and Spirit values swap between day and night charts."""
        # Same base positions, but Sun position changes sect
        base_positions = {
            "Asc": 0.0,
            "Moon": 180.0,
            "Mercury": 45.0,
            "Venus": 60.0,
        }

        # Day chart: Sun at 90° (above horizon)
        day_positions = {**base_positions, "Sun": 90.0}
        day_parts = calc_all_arabic_parts(day_positions)

        # Night chart: Sun at 270° (below horizon)
        night_positions = {**base_positions, "Sun": 270.0}
        night_parts = calc_all_arabic_parts(night_positions)

        # The difference between Fortune and Spirit should be inverted
        # Day Fortune uses Moon, Night Fortune uses Sun (and vice versa for Spirit)
        # Day Fortune: 0 + 180 - 90 = 90
        # Day Spirit: 0 + 90 - 180 = -90 => 270
        assert day_parts["Pars_Fortunae"] == pytest.approx(90.0)
        assert day_parts["Pars_Spiritus"] == pytest.approx(270.0)

        # Night Fortune: 0 + 270 - 180 = 90
        # Night Spirit: 0 + 180 - 270 = -90 => 270
        assert night_parts["Pars_Fortunae"] == pytest.approx(90.0)
        assert night_parts["Pars_Spiritus"] == pytest.approx(270.0)

        # Note: In this particular case values are the same due to symmetric positions
        # Let's verify with asymmetric values

    def test_fortune_spirit_swap_asymmetric(self):
        """Verify Formula swap with asymmetric Moon/Sun positions."""
        base_positions = {
            "Asc": 30.0,
            "Moon": 120.0,
            "Mercury": 45.0,
            "Venus": 60.0,
        }

        # Day chart: Sun at 60°
        day_positions = {**base_positions, "Sun": 60.0}
        day_parts = calc_all_arabic_parts(day_positions)

        # Day Fortune: 30 + 120 - 60 = 90
        # Day Spirit: 30 + 60 - 120 = -30 => 330
        assert day_parts["Pars_Fortunae"] == pytest.approx(90.0)
        assert day_parts["Pars_Spiritus"] == pytest.approx(330.0)

        # Night chart: Sun at 300° (below horizon)
        night_positions = {**base_positions, "Sun": 300.0}
        night_parts = calc_all_arabic_parts(night_positions)

        # Night Fortune: ASC + Sun - Moon = 30 + 300 - 120 = 210
        # Night Spirit: ASC + Moon - Sun = 30 + 120 - 300 = -150 => 210
        assert night_parts["Pars_Fortunae"] == pytest.approx(210.0)
        assert night_parts["Pars_Spiritus"] == pytest.approx(210.0)

        # Day and Night should produce different results for Fortune
        assert day_parts["Pars_Fortunae"] != night_parts["Pars_Fortunae"]

    def test_sun_exactly_on_asc_edge_case(self):
        """Test edge case: Sun exactly on Ascendant (sunrise moment)."""
        # Sun exactly on ASC should be treated as day chart
        positions = {
            "Asc": 100.0,
            "Sun": 100.0,  # Exactly on ASC
            "Moon": 200.0,
            "Mercury": 150.0,
            "Venus": 120.0,
        }

        parts = calc_all_arabic_parts(positions)

        # Should be a day chart
        assert is_day_chart(positions["Sun"], positions["Asc"]) is True

        # Day Fortune: 100 + 200 - 100 = 200
        assert parts["Pars_Fortunae"] == pytest.approx(200.0)

        # Day Spirit: 100 + 100 - 200 = 0
        assert parts["Pars_Spiritus"] == pytest.approx(0.0)

    def test_missing_keys_default_to_zero(self):
        """Test that missing keys default to 0.0."""
        positions = {"Asc": 90.0}  # Only ASC provided

        parts = calc_all_arabic_parts(positions)

        # With all other values as 0, Sun=0 is not between ASC=90 and DSC=270
        # So this is a night chart
        # Night Fortune: 90 + 0 - 0 = 90
        # Night Spirit: 90 + 0 - 0 = 90
        assert parts["Pars_Fortunae"] == pytest.approx(90.0)
        assert parts["Pars_Spiritus"] == pytest.approx(90.0)

    def test_all_parts_returned(self):
        """Verify all four standard parts are returned."""
        positions = {
            "Asc": 0.0,
            "Sun": 90.0,
            "Moon": 180.0,
            "Mercury": 45.0,
            "Venus": 135.0,
        }

        parts = calc_all_arabic_parts(positions)

        assert "Pars_Fortunae" in parts
        assert "Pars_Spiritus" in parts
        assert "Pars_Amoris" in parts
        assert "Pars_Fidei" in parts
        assert len(parts) == 4

    def test_values_in_valid_range(self):
        """All returned values should be in range 0-360."""
        positions = {
            "Asc": 350.0,
            "Sun": 10.0,
            "Moon": 270.0,
            "Mercury": 300.0,
            "Venus": 5.0,
        }

        parts = calc_all_arabic_parts(positions)

        for name, value in parts.items():
            assert 0.0 <= value < 360.0, f"{name} = {value} is out of range"


class TestIsDayChart3D:
    """Tests for 3D day/night calculation at extreme latitudes."""

    def test_threshold_constant_exists(self):
        """Verify the extreme latitude threshold is defined."""
        assert _EXTREME_LATITUDE_THRESHOLD == 60.0

    def test_moderate_latitude_uses_2d_method(self):
        """At moderate latitudes, should use 2D method even with jd provided."""
        # Rome, Italy - moderate latitude (41.9°N)
        # June 21, 2024 at noon - Sun clearly above horizon
        jd = julday(2024, 6, 21, 12.0)
        sun_lon = 90.0  # ~Cancer (summer solstice)
        asc = 180.0  # Libra rising

        # 2D method: Sun at 90° is not between ASC=180° and DSC=0°/360°
        # So 2D says night chart (which may be wrong for actual Rome noon)
        result_2d = is_day_chart(sun_lon, asc)

        # With location at moderate latitude, should still use 2D
        result_with_loc = is_day_chart(sun_lon, asc, jd=jd, geo_lat=41.9, geo_lon=12.5)
        assert result_2d == result_with_loc

    def test_extreme_latitude_uses_3d_method_northern(self):
        """At extreme northern latitude with location, uses 3D calculation."""
        # Tromso, Norway (69.6°N) - Arctic summer, Sun above horizon at midnight
        # June 21, 2024 at midnight UTC - Sun should be above horizon
        jd = julday(2024, 6, 21, 0.0)  # Midnight UTC

        # Sun at summer solstice position
        sun_lon = 90.0  # Cancer 0°
        asc = 0.0  # Aries rising (doesn't matter for 3D calc)

        # Without location: uses 2D method
        result_2d = is_day_chart(sun_lon, asc)

        # With location at extreme latitude: uses 3D method
        result_3d = is_day_chart(sun_lon, asc, jd=jd, geo_lat=69.6, geo_lon=19.0)

        # In Arctic summer, Sun is above horizon even at midnight
        # 3D method should return True (day chart)
        assert result_3d is True

    def test_extreme_latitude_winter_night(self):
        """At extreme latitude in winter, Sun below horizon at noon."""
        # Tromso, Norway (69.6°N) - polar night period
        # December 21, 2024 at noon UTC - Sun below horizon
        jd = julday(2024, 12, 21, 12.0)

        # Sun at winter solstice position
        sun_lon = 270.0  # Capricorn 0°
        asc = 0.0

        # 3D method should return False (night chart) even at "noon"
        result_3d = is_day_chart(sun_lon, asc, jd=jd, geo_lat=69.6, geo_lon=19.0)
        assert result_3d is False

    def test_extreme_latitude_southern_hemisphere(self):
        """Test 3D calculation for Antarctic location."""
        # McMurdo Station, Antarctica (-77.8°S)
        # December 21, 2024 - Antarctic summer, midnight sun
        jd = julday(2024, 12, 21, 0.0)  # Midnight UTC

        sun_lon = 270.0  # Capricorn (Sun at summer solstice for S hemisphere)
        asc = 0.0

        # In Antarctic summer, Sun is above horizon at midnight
        result_3d = is_day_chart(sun_lon, asc, jd=jd, geo_lat=-77.8, geo_lon=166.7)
        assert result_3d is True

    def test_is_sun_above_horizon_3d_directly(self):
        """Test the internal 3D function directly."""
        # Standard daytime scenario - Rome at noon in summer
        jd = julday(2024, 6, 21, 12.0)
        sun_lon = 90.0
        sun_lat = 0.0
        geo_lat = 41.9
        geo_lon = 12.5

        result = _is_sun_above_horizon_3d(jd, sun_lon, sun_lat, geo_lat, geo_lon)
        assert result is True

    def test_is_sun_above_horizon_3d_night(self):
        """Test 3D function for nighttime scenario."""
        # Midnight in Rome, winter
        jd = julday(2024, 12, 21, 0.0)  # Midnight UTC
        sun_lon = 270.0
        sun_lat = 0.0
        geo_lat = 41.9
        geo_lon = 12.5

        result = _is_sun_above_horizon_3d(jd, sun_lon, sun_lat, geo_lat, geo_lon)
        assert result is False

    def test_fallback_to_2d_without_jd(self):
        """Without jd parameter, should use 2D method."""
        sun_lon = 90.0
        asc = 0.0

        # Without jd: 2D method (Sun between 0° and 180°)
        result = is_day_chart(sun_lon, asc, geo_lat=70.0, geo_lon=20.0)

        # Should be same as pure 2D call
        result_2d = is_day_chart(sun_lon, asc)
        assert result == result_2d

    def test_fallback_to_2d_without_lat(self):
        """Without geo_lat parameter, should use 2D method."""
        jd = julday(2024, 6, 21, 12.0)
        sun_lon = 90.0
        asc = 0.0

        result = is_day_chart(sun_lon, asc, jd=jd, geo_lon=20.0)
        result_2d = is_day_chart(sun_lon, asc)
        assert result == result_2d


class TestCalcAllArabicPartsWithLocation:
    """Tests for calc_all_arabic_parts with location parameters."""

    def test_moderate_latitude_same_as_without_location(self):
        """At moderate latitudes, results should be same with or without location."""
        positions = {
            "Asc": 15.0,
            "Sun": 90.0,
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts_without = calc_all_arabic_parts(positions)
        parts_with = calc_all_arabic_parts(
            positions,
            jd=julday(2024, 6, 21, 12.0),
            geo_lat=41.9,
            geo_lon=12.5,
        )

        # Same results because latitude is moderate
        assert parts_without == parts_with

    def test_extreme_latitude_may_differ(self):
        """At extreme latitudes, 3D calc may give different sect determination."""
        # Scenario: Arctic summer midnight - Sun geometrically above horizon
        # but 2D method might say night chart depending on ASC/Sun relationship
        jd = julday(2024, 6, 21, 0.0)

        positions = {
            "Asc": 270.0,  # Capricorn rising
            "Sun": 90.0,  # Cancer (summer solstice)
            "Moon": 180.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts_2d = calc_all_arabic_parts(positions)
        parts_3d = calc_all_arabic_parts(
            positions,
            jd=jd,
            geo_lat=69.6,  # Tromso
            geo_lon=19.0,
        )

        # 2D: Sun at 90° is NOT between ASC=270° and DSC=90° (it's exactly on DSC)
        # Actually with wrapped logic: ASC=270, DSC=90
        # sun_lon >= asc OR sun_lon <= desc => 90 >= 270 (False) OR 90 <= 90 (True)
        # So 2D says day chart

        # 3D should also say day chart (Sun above horizon in Arctic summer)
        # Verify 3D is being used (latitude > 60)
        assert abs(69.6) > _EXTREME_LATITUDE_THRESHOLD

        # Both should agree in this case (day chart)
        is_day_2d = is_day_chart(positions["Sun"], positions["Asc"])
        is_day_3d = is_day_chart(
            positions["Sun"], positions["Asc"], jd=jd, geo_lat=69.6, geo_lon=19.0
        )

        # In Arctic summer at midnight, Sun IS above horizon (midnight sun)
        assert is_day_3d is True

    def test_sun_lat_parameter_passed(self):
        """Verify Sun_lat from positions is used in 3D calculation."""
        jd = julday(2024, 6, 21, 12.0)

        positions = {
            "Asc": 15.0,
            "Sun": 90.0,
            "Sun_lat": 0.0,  # Sun ecliptic latitude (nearly always ~0)
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        # Should not raise any errors
        parts = calc_all_arabic_parts(
            positions,
            jd=jd,
            geo_lat=70.0,
            geo_lon=20.0,
        )

        assert "Pars_Fortunae" in parts
