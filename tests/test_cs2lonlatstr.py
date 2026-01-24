"""
Tests for the cs2lonlatstr centiseconds to longitude/latitude string conversion function.

Tests verify that cs2lonlatstr correctly converts angles in centiseconds to
formatted degree strings with directional characters (e.g., "45°30'00\" N").
"""

import libephemeris as ephem


# Constants for centiseconds calculations
CS1 = 3600 * 100  # 360000 centiseconds per degree


class TestCs2lonlatstrBasic:
    """Basic functionality tests for cs2lonlatstr."""

    def test_cs2lonlatstr_exported(self):
        """Test that cs2lonlatstr is exported from the package."""
        assert hasattr(ephem, "cs2lonlatstr")
        assert callable(ephem.cs2lonlatstr)

    def test_cs2lonlatstr_zero(self):
        """Test conversion of zero."""
        result = ephem.cs2lonlatstr(0, "N", "S")
        assert "0" in result
        assert "°" in result
        assert "'" in result
        assert '"' in result
        assert "N" in result

    def test_cs2lonlatstr_returns_string(self):
        """Test that cs2lonlatstr returns a string."""
        result = ephem.cs2lonlatstr(CS1, "E", "W")
        assert isinstance(result, str)

    def test_cs2lonlatstr_one_degree(self):
        """Test conversion of exactly 1 degree."""
        result = ephem.cs2lonlatstr(CS1, "N", "S")
        # Should contain 1 degree and N
        assert "1°" in result
        assert "N" in result

    def test_cs2lonlatstr_contains_degree_symbol(self):
        """Test that result contains the degree symbol."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert "°" in result

    def test_cs2lonlatstr_contains_minute_symbol(self):
        """Test that result contains the minute symbol (')."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert "'" in result

    def test_cs2lonlatstr_contains_second_symbol(self):
        """Test that result contains the second symbol (\")."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert '"' in result


class TestCs2lonlatstrDirectionalCharacters:
    """Tests for directional character handling."""

    def test_positive_uses_plus_char(self):
        """Test that positive values use plus_char."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert "N" in result
        assert "S" not in result

    def test_negative_uses_minus_char(self):
        """Test that negative values use minus_char."""
        result = ephem.cs2lonlatstr(-CS1 * 45, "N", "S")
        assert "S" in result
        assert "N" not in result

    def test_zero_uses_plus_char(self):
        """Test that zero uses plus_char."""
        result = ephem.cs2lonlatstr(0, "N", "S")
        assert "N" in result
        assert "S" not in result

    def test_longitude_east_west(self):
        """Test longitude with E/W characters."""
        result_east = ephem.cs2lonlatstr(CS1 * 122, "E", "W")
        assert "E" in result_east
        assert "W" not in result_east

        result_west = ephem.cs2lonlatstr(-CS1 * 122, "E", "W")
        assert "W" in result_west
        assert "E" not in result_west

    def test_custom_direction_chars(self):
        """Test with custom directional characters."""
        result = ephem.cs2lonlatstr(CS1 * 30, "+", "-")
        assert "+" in result

        result_neg = ephem.cs2lonlatstr(-CS1 * 30, "+", "-")
        assert "-" in result_neg


class TestCs2lonlatstrValues:
    """Tests for various value conversions."""

    def test_45_degrees_30_minutes(self):
        """Test conversion of 45°30'."""
        # 45 degrees 30 minutes = 45*360000 + 30*6000 = 16380000 cs
        cs_value = 45 * CS1 + 30 * 6000
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "45°" in result
        assert "30'" in result
        assert "N" in result

    def test_122_degrees_15_minutes_30_seconds(self):
        """Test conversion of 122°15'30\"."""
        # 122 degrees 15 minutes 30 seconds
        cs_value = 122 * CS1 + 15 * 6000 + 30 * 100
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert "122°" in result
        assert "15'" in result
        assert '30"' in result
        assert "E" in result

    def test_negative_longitude(self):
        """Test negative longitude (West)."""
        cs_value = -(122 * CS1 + 15 * 6000 + 30 * 100)
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert "122°" in result
        assert "15'" in result
        assert '30"' in result
        assert "W" in result

    def test_90_degrees_latitude(self):
        """Test 90° latitude (North Pole)."""
        result = ephem.cs2lonlatstr(90 * CS1, "N", "S")
        assert "90°" in result
        assert "N" in result

    def test_negative_90_degrees_latitude(self):
        """Test -90° latitude (South Pole)."""
        result = ephem.cs2lonlatstr(-90 * CS1, "N", "S")
        assert "90°" in result
        assert "S" in result


class TestCs2lonlatstrRounding:
    """Tests for seconds rounding behavior."""

    def test_round_up_centiseconds(self):
        """Test that centiseconds >= 50 round up to next second."""
        # 1 degree 0 minutes 0.50 seconds should round to 1 second
        cs_value = CS1 + 50
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "1°" in result
        assert ' 1"' in result or "'1\"" in result or "' 1\"" in result

    def test_round_down_centiseconds(self):
        """Test that centiseconds < 50 round down."""
        # 1 degree 0 minutes 0.49 seconds should round to 0 seconds
        cs_value = CS1 + 49
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "1°" in result
        assert ' 0"' in result

    def test_seconds_carry_to_minutes(self):
        """Test that 60 seconds rounds to next minute."""
        # 1 degree 59 minutes 59.50 seconds should round to 2 degrees 0 minutes
        cs_value = CS1 + 59 * 6000 + 59 * 100 + 50
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "2°" in result
        assert " 0'" in result

    def test_no_centiseconds_in_output(self):
        """Test that output does not contain centisecond decimals."""
        # Unlike cs2degstr, cs2lonlatstr should show whole seconds only
        cs_value = CS1 + 45 * 100 + 78  # 1°0'45.78"
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        # Should show "46" (rounded from 45.78) not "45.78"
        assert ".78" not in result


class TestCs2lonlatstrFormat:
    """Tests for output format consistency."""

    def test_format_ends_with_direction(self):
        """Test that the format ends with directional character."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert result.rstrip().endswith("N")

        result_neg = ephem.cs2lonlatstr(-CS1 * 45, "N", "S")
        assert result_neg.rstrip().endswith("S")

    def test_format_contains_all_parts(self):
        """Test that the format contains degrees, minutes, and seconds."""
        result = ephem.cs2lonlatstr(45 * CS1 + 30 * 6000 + 15 * 100, "N", "S")
        # Should have degree symbol, minute symbol, double quote, and direction
        assert "°" in result
        assert "'" in result
        assert '"' in result
        assert "N" in result


class TestCs2lonlatstrEdgeCases:
    """Edge case tests for cs2lonlatstr."""

    def test_small_positive_value(self):
        """Test very small positive value."""
        result = ephem.cs2lonlatstr(1, "N", "S")
        assert "0°" in result
        assert "N" in result

    def test_small_negative_value(self):
        """Test very small negative value."""
        result = ephem.cs2lonlatstr(-1, "N", "S")
        assert "0°" in result
        assert "S" in result

    def test_180_degrees_longitude(self):
        """Test 180° longitude (dateline)."""
        result = ephem.cs2lonlatstr(180 * CS1, "E", "W")
        assert "180°" in result
        assert "E" in result

    def test_large_value(self):
        """Test large value beyond normal coordinate range."""
        result = ephem.cs2lonlatstr(360 * CS1, "E", "W")
        assert "360°" in result


class TestCs2lonlatstrTypicalUseCase:
    """Tests for typical geographic coordinate use cases."""

    def test_rome_latitude(self):
        """Test latitude of Rome (~41°54' N)."""
        # 41 degrees 54 minutes
        cs_value = 41 * CS1 + 54 * 6000
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "41°" in result
        assert "54'" in result
        assert "N" in result

    def test_rome_longitude(self):
        """Test longitude of Rome (~12°30' E)."""
        # 12 degrees 30 minutes
        cs_value = 12 * CS1 + 30 * 6000
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert "12°" in result
        assert "30'" in result
        assert "E" in result

    def test_new_york_latitude(self):
        """Test latitude of New York (~40°43' N)."""
        cs_value = 40 * CS1 + 43 * 6000
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "40°" in result
        assert "43'" in result
        assert "N" in result

    def test_new_york_longitude(self):
        """Test longitude of New York (~74°00' W)."""
        cs_value = -(74 * CS1)
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert "74°" in result
        assert "W" in result

    def test_sydney_latitude(self):
        """Test latitude of Sydney (~33°52' S)."""
        cs_value = -(33 * CS1 + 52 * 6000)
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert "33°" in result
        assert "52'" in result
        assert "S" in result

    def test_sydney_longitude(self):
        """Test longitude of Sydney (~151°12' E)."""
        cs_value = 151 * CS1 + 12 * 6000
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert "151°" in result
        assert "12'" in result
        assert "E" in result
