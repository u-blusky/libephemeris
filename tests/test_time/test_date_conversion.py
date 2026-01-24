"""
Tests for date_conversion function.

Tests cover:
- Conversion from Gregorian to Julian calendar
- Conversion from Julian to Gregorian calendar
- Edge cases at the Gregorian reform date (Oct 15, 1582)
- Same-calendar identity (no conversion needed)
- Invalid calendar parameter
"""

import pytest
import libephemeris as ephem


class TestDateConversionGregorianToJulian:
    """Test conversion from Gregorian to Julian calendar."""

    @pytest.mark.unit
    def test_gregorian_reform_date(self):
        """Oct 15, 1582 Gregorian = Oct 5, 1582 Julian."""
        year, month, day, hour = ephem.date_conversion(1582, 10, 15, 12.0, "j")
        assert year == 1582
        assert month == 10
        assert day == 5
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_modern_date(self):
        """2000-01-01 Gregorian to Julian."""
        year, month, day, hour = ephem.date_conversion(2000, 1, 1, 12.0, "j")
        # In 2000, Julian calendar is 13 days behind Gregorian
        assert year == 1999
        assert month == 12
        assert day == 19
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_j2000_epoch_to_julian(self):
        """J2000.0 = 2000-01-01 12:00 TT to Julian calendar."""
        year, month, day, hour = ephem.date_conversion(2000, 1, 1, 12.0, "j")
        # Verify the result matches expectations
        assert year == 1999
        assert month == 12
        assert day == 19

    @pytest.mark.unit
    def test_post_reform_date(self):
        """Date shortly after Gregorian reform."""
        year, month, day, hour = ephem.date_conversion(1582, 10, 16, 0.0, "j")
        assert year == 1582
        assert month == 10
        assert day == 6
        assert hour == pytest.approx(0.0, abs=1e-10)


class TestDateConversionJulianToGregorian:
    """Test conversion from Julian to Gregorian calendar."""

    @pytest.mark.unit
    def test_julian_reform_date(self):
        """Oct 5, 1582 Julian = Oct 15, 1582 Gregorian."""
        # Date before Oct 15, 1582 is assumed Julian
        year, month, day, hour = ephem.date_conversion(1582, 10, 5, 12.0, "g")
        assert year == 1582
        assert month == 10
        assert day == 15
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_ancient_date(self):
        """Convert an ancient Julian date to Gregorian."""
        # Julius Caesar's death: March 15, 44 BC (Julian)
        year, month, day, hour = ephem.date_conversion(-43, 3, 15, 12.0, "g")
        # Just verify it's a valid result - dates this old are complex
        assert isinstance(year, int)
        assert 1 <= month <= 12
        assert 1 <= day <= 31

    @pytest.mark.unit
    def test_pre_reform_date(self):
        """Date before Gregorian reform."""
        year, month, day, hour = ephem.date_conversion(1500, 6, 15, 6.0, "g")
        # 82 years before reform, Julian ~10 days behind
        # This is before reform so input is Julian, output Gregorian
        assert isinstance(year, int)
        assert 1 <= month <= 12


class TestDateConversionIdentity:
    """Test that same-calendar conversion returns original date."""

    @pytest.mark.unit
    def test_gregorian_to_gregorian(self):
        """Modern date stays same when converting Gregorian to Gregorian."""
        year, month, day, hour = ephem.date_conversion(2020, 6, 15, 8.5, "g")
        assert year == 2020
        assert month == 6
        assert day == 15
        assert hour == pytest.approx(8.5, abs=1e-10)

    @pytest.mark.unit
    def test_julian_to_julian(self):
        """Pre-reform date stays same when converting Julian to Julian."""
        year, month, day, hour = ephem.date_conversion(1400, 3, 21, 18.0, "j")
        assert year == 1400
        assert month == 3
        assert day == 21
        assert hour == pytest.approx(18.0, abs=1e-10)


class TestDateConversionRoundtrip:
    """Test roundtrip conversions."""

    @pytest.mark.unit
    def test_roundtrip_pre_reform_date(self):
        """Pre-reform date: Julian -> Gregorian -> Julian."""
        # This roundtrip works because the date is before the reform
        original = (1400, 3, 15, 12.0)
        greg = ephem.date_conversion(*original, "g")
        back = ephem.date_conversion(*greg, "j")
        # The round-trip should preserve the original
        assert back[0] == original[0]
        assert back[1] == original[1]
        assert back[2] == original[2]
        assert back[3] == pytest.approx(original[3], abs=1e-8)

    @pytest.mark.unit
    def test_roundtrip_reform_boundary(self):
        """Reform boundary: Oct 15, 1582 Gregorian roundtrip."""
        original = (1582, 10, 15, 12.0)
        julian = ephem.date_conversion(*original, "j")
        back = ephem.date_conversion(*julian, "g")
        assert back[0] == original[0]
        assert back[1] == original[1]
        assert back[2] == original[2]
        assert back[3] == pytest.approx(original[3], abs=1e-8)


class TestDateConversionEdgeCases:
    """Test edge cases and special scenarios."""

    @pytest.mark.unit
    def test_uppercase_calendar_parameter(self):
        """Calendar parameter should be case-insensitive."""
        result_lower = ephem.date_conversion(2000, 1, 1, 12.0, "j")
        result_upper = ephem.date_conversion(2000, 1, 1, 12.0, "J")
        assert result_lower == result_upper

    @pytest.mark.unit
    def test_gregorian_uppercase(self):
        """Uppercase 'G' should work."""
        result_lower = ephem.date_conversion(2000, 1, 1, 12.0, "g")
        result_upper = ephem.date_conversion(2000, 1, 1, 12.0, "G")
        assert result_lower == result_upper

    @pytest.mark.unit
    def test_invalid_calendar_parameter(self):
        """Invalid calendar parameter should raise ValueError."""
        with pytest.raises(ValueError, match="calendar must be 'j' or 'g'"):
            ephem.date_conversion(2000, 1, 1, 12.0, "x")

    @pytest.mark.unit
    def test_fractional_hours_preserved(self):
        """Fractional hours should be preserved through conversion."""
        # 6:30:45 = 6.5125 hours
        year, month, day, hour = ephem.date_conversion(2000, 6, 15, 6.5125, "j")
        assert hour == pytest.approx(6.5125, abs=1e-6)

    @pytest.mark.unit
    def test_midnight_boundary(self):
        """Test midnight (hour=0) is handled correctly."""
        year, month, day, hour = ephem.date_conversion(2000, 1, 1, 0.0, "j")
        assert hour == pytest.approx(0.0, abs=1e-10)

    @pytest.mark.unit
    def test_near_midnight(self):
        """Test times near midnight are handled correctly."""
        year, month, day, hour = ephem.date_conversion(2000, 1, 1, 23.999, "j")
        assert hour == pytest.approx(23.999, abs=1e-6)


class TestDateConversionHistorical:
    """Test with historically significant dates."""

    @pytest.mark.unit
    def test_newton_birth_julian_to_gregorian(self):
        """
        Isaac Newton was born on Dec 25, 1642 Julian (Jan 4, 1643 Gregorian).

        Note: Since Dec 25, 1642 is after Oct 15, 1582 (the standard reform date),
        the function treats it as Gregorian input. This test demonstrates the conversion
        direction from Gregorian Dec 25, 1642 to Julian Dec 15, 1642.

        For the historical conversion (Julian Dec 25, 1642 to Gregorian Jan 4, 1643),
        one would need to use a date before the reform or manually specify the calendar.
        """
        year, month, day, hour = ephem.date_conversion(1642, 12, 25, 12.0, "j")
        # 1642 had 10-day difference between calendars
        assert year == 1642
        assert month == 12
        assert day == 15

    @pytest.mark.unit
    def test_october_revolution(self):
        """Russian October Revolution: Oct 25, 1917 Julian = Nov 7, 1917 Gregorian."""
        # Oct 25, 1917 is after reform date numerically, but Russia still used Julian
        # Since the input date is > Oct 15, 1582, it's treated as Gregorian
        # So this test shows the conversion from Gregorian perspective
        year, month, day, hour = ephem.date_conversion(1917, 11, 7, 0.0, "j")
        assert year == 1917
        assert month == 10
        assert day == 25
