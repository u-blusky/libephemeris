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
    calc_all_arabic_parts,
    calc_arabic_part_of_faith,
    calc_arabic_part_of_fortune,
    calc_arabic_part_of_love,
    calc_arabic_part_of_spirit,
    is_day_chart,
)


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
