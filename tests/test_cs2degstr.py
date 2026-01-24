"""
Tests for the cs2degstr centiseconds to degrees string conversion function.

Tests verify that cs2degstr correctly converts angles in centiseconds to
formatted degree strings (e.g., "123°45'06.78\"").
"""

import libephemeris as ephem


# Constants for centiseconds calculations
CS360 = 360 * 3600 * 100  # 129600000 centiseconds in a full circle
CS180 = 180 * 3600 * 100  # 64800000 centiseconds in a half circle
CS1 = 3600 * 100  # 360000 centiseconds per degree


class TestCs2degstrBasic:
    """Basic functionality tests for cs2degstr."""

    def test_cs2degstr_exported(self):
        """Test that cs2degstr is exported from the package."""
        assert hasattr(ephem, "cs2degstr")
        assert callable(ephem.cs2degstr)

    def test_cs2degstr_zero(self):
        """Test conversion of zero."""
        result = ephem.cs2degstr(0)
        assert "0" in result
        assert "°" in result
        assert "'" in result
        assert '"' in result

    def test_cs2degstr_returns_string(self):
        """Test that cs2degstr returns a string."""
        result = ephem.cs2degstr(CS1)
        assert isinstance(result, str)

    def test_cs2degstr_one_degree(self):
        """Test conversion of exactly 1 degree."""
        result = ephem.cs2degstr(CS1)
        # Should be "  1° 0' 0.00\""
        assert "1°" in result
        assert "0'" in result
        assert '0.00"' in result

    def test_cs2degstr_contains_degree_symbol(self):
        """Test that result contains the degree symbol."""
        result = ephem.cs2degstr(CS1 * 45)
        assert "°" in result

    def test_cs2degstr_contains_minute_symbol(self):
        """Test that result contains the minute symbol (')."""
        result = ephem.cs2degstr(CS1 * 45)
        assert "'" in result

    def test_cs2degstr_contains_second_symbol(self):
        """Test that result contains the second symbol (\")."""
        result = ephem.cs2degstr(CS1 * 45)
        assert '"' in result


class TestCs2degstrValues:
    """Tests for various value conversions."""

    def test_cs2degstr_10_degrees(self):
        """Test conversion of 10 degrees."""
        result = ephem.cs2degstr(CS1 * 10)
        assert "10°" in result

    def test_cs2degstr_90_degrees(self):
        """Test conversion of 90 degrees."""
        result = ephem.cs2degstr(CS1 * 90)
        assert "90°" in result

    def test_cs2degstr_180_degrees(self):
        """Test conversion of 180 degrees."""
        result = ephem.cs2degstr(CS180)
        assert "180°" in result

    def test_cs2degstr_359_degrees(self):
        """Test conversion of 359 degrees."""
        result = ephem.cs2degstr(CS1 * 359)
        assert "359°" in result

    def test_cs2degstr_with_minutes(self):
        """Test conversion with minutes component."""
        # 1 degree 30 minutes = 360000 + 30*6000 = 360000 + 180000 = 540000 cs
        cs_value = CS1 + (30 * 6000)
        result = ephem.cs2degstr(cs_value)
        assert "1°" in result
        assert "30'" in result

    def test_cs2degstr_with_seconds(self):
        """Test conversion with seconds component."""
        # 1 degree 0 minutes 45 seconds = 360000 + 0 + 4500 = 364500 cs
        cs_value = CS1 + (45 * 100)
        result = ephem.cs2degstr(cs_value)
        assert "1°" in result
        assert "45." in result  # seconds part

    def test_cs2degstr_with_centiseconds(self):
        """Test conversion with centiseconds component."""
        # 1 degree 0 minutes 0.50 seconds = 360000 + 50 = 360050 cs
        cs_value = CS1 + 50
        result = ephem.cs2degstr(cs_value)
        assert "1°" in result
        assert ".50" in result  # centiseconds part

    def test_cs2degstr_complex_value(self):
        """Test conversion of a complex value with all components."""
        # 123 degrees 45 minutes 06.78 seconds
        # = 123*360000 + 45*6000 + 6*100 + 78
        # = 44280000 + 270000 + 600 + 78
        # = 44550678
        cs_value = 123 * CS1 + 45 * 6000 + 6 * 100 + 78
        result = ephem.cs2degstr(cs_value)
        assert "123°" in result
        assert "45'" in result
        assert "6.78" in result or " 6.78" in result


class TestCs2degstrNegativeValues:
    """Tests for negative value handling."""

    def test_cs2degstr_negative_1_degree(self):
        """Test conversion of -1 degree."""
        result = ephem.cs2degstr(-CS1)
        assert "-1°" in result
        assert "0'" in result
        assert '0.00"' in result

    def test_cs2degstr_negative_45_degrees(self):
        """Test conversion of -45 degrees."""
        result = ephem.cs2degstr(-CS1 * 45)
        assert "-45°" in result

    def test_cs2degstr_negative_with_minutes(self):
        """Test conversion of negative value with minutes."""
        cs_value = -(CS1 + 30 * 6000)  # -1 degree 30 minutes
        result = ephem.cs2degstr(cs_value)
        assert "-1°" in result
        assert "30'" in result


class TestCs2degstrLargeValues:
    """Tests for large value handling."""

    def test_cs2degstr_360_degrees(self):
        """Test conversion of 360 degrees (full circle)."""
        result = ephem.cs2degstr(CS360)
        assert "360°" in result

    def test_cs2degstr_greater_than_360(self):
        """Test conversion of value greater than 360 degrees."""
        result = ephem.cs2degstr(CS1 * 400)  # 400 degrees
        assert "400°" in result

    def test_cs2degstr_1000_degrees(self):
        """Test conversion of 1000 degrees."""
        result = ephem.cs2degstr(CS1 * 1000)
        assert "1000°" in result


class TestCs2degstrEdgeCases:
    """Edge case tests for cs2degstr."""

    def test_cs2degstr_1_centisecond(self):
        """Test conversion of 1 centisecond."""
        result = ephem.cs2degstr(1)
        assert "0°" in result
        assert "0'" in result
        assert ".01" in result

    def test_cs2degstr_99_centiseconds(self):
        """Test conversion of 99 centiseconds (just under 1 arcsecond)."""
        result = ephem.cs2degstr(99)
        assert "0°" in result
        assert ".99" in result

    def test_cs2degstr_100_centiseconds(self):
        """Test conversion of 100 centiseconds (exactly 1 arcsecond)."""
        result = ephem.cs2degstr(100)
        assert "0°" in result
        assert "1.00" in result

    def test_cs2degstr_5999_centiseconds(self):
        """Test conversion of 5999 centiseconds (just under 1 arcminute)."""
        result = ephem.cs2degstr(5999)
        # 5999 cs = 0 degrees, 0 minutes, 59.99 seconds
        assert "0°" in result
        assert "0'" in result
        assert "59.99" in result

    def test_cs2degstr_6000_centiseconds(self):
        """Test conversion of 6000 centiseconds (exactly 1 arcminute)."""
        result = ephem.cs2degstr(6000)
        # 6000 cs = 0 degrees, 1 minute, 0.00 seconds
        assert "0°" in result
        assert "1'" in result
        assert "0.00" in result


class TestCs2degstrFormat:
    """Tests for output format consistency."""

    def test_cs2degstr_format_contains_all_parts(self):
        """Test that the format contains degrees, minutes, and seconds."""
        result = ephem.cs2degstr(45 * CS1 + 30 * 6000 + 15 * 100 + 25)
        # Should have degree symbol, minute symbol, and double quote
        assert "°" in result
        assert "'" in result
        assert '"' in result

    def test_cs2degstr_seconds_has_two_decimals(self):
        """Test that seconds always has 2 decimal places."""
        result = ephem.cs2degstr(CS1)  # 1 degree exactly
        # Should end with ".00\""
        assert '.00"' in result

    def test_cs2degstr_leading_zeros_in_centiseconds(self):
        """Test that centiseconds have leading zeros when needed."""
        result = ephem.cs2degstr(CS1 + 5)  # 1 degree + 0.05 seconds
        # Should have ".05" not ".5"
        assert ".05" in result


class TestCs2degstrRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_angles_in_centiseconds(self, random_longitudes):
        """Test with random angles converted to centiseconds."""
        lons = random_longitudes(50)

        for lon in lons:
            # Convert longitude to centiseconds
            cs_value = int(lon * CS1)

            # Test positive value
            result = ephem.cs2degstr(cs_value)
            assert isinstance(result, str), f"Failed for cs value {cs_value}"
            assert "°" in result, f"Missing degree symbol for {cs_value}"
            assert "'" in result, f"Missing minute symbol for {cs_value}"
            assert '"' in result, f"Missing second symbol for {cs_value}"

            # Test negative value
            result_neg = ephem.cs2degstr(-cs_value)
            assert isinstance(result_neg, str), f"Failed for cs value {-cs_value}"
            if cs_value > 0:
                assert "-" in result_neg, f"Missing minus sign for {-cs_value}"

    def test_round_trip_consistency(self, random_longitudes):
        """Test that values can be extracted from the string correctly."""
        lons = random_longitudes(20)

        for lon in lons:
            # Convert longitude to centiseconds (integer)
            cs_value = int(lon * CS1)

            # Get the string representation
            result = ephem.cs2degstr(cs_value)

            # Parse back the components
            # Expected format: "XXX°YY'ZZ.cc\""
            import re

            match = re.match(r"\s*(-?\d+)°\s*(\d+)'\s*(\d+)\.(\d+)\"", result)
            assert match is not None, f"Could not parse result: {result}"

            deg = int(match.group(1))
            mins = int(match.group(2))
            secs = int(match.group(3))
            cs = int(match.group(4))

            # Reconstruct the value
            if deg >= 0:
                reconstructed = deg * CS1 + mins * 6000 + secs * 100 + cs
            else:
                reconstructed = deg * CS1 - mins * 6000 - secs * 100 - cs

            assert reconstructed == cs_value, (
                f"Round trip failed: original={cs_value}, "
                f"reconstructed={reconstructed}, string={result}"
            )
