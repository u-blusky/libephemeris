"""
Comprehensive tests for utc_time_zone() function.

Tests cover:
- Basic timezone offset application
- Positive offsets (east of UTC: CET, JST, etc.)
- Negative offsets (west of UTC: EST, PST, etc.)
- Day boundary crossings
- Month boundary crossings
- Year boundary crossings
- Fractional hour offsets
- Edge cases (leap years, UTC+0)
"""

import pytest
import libephemeris as ephem


class TestUtcTimeZoneBasic:
    """Test basic timezone offset application."""

    @pytest.mark.unit
    def test_returns_tuple(self):
        """Function should return a tuple of 6 elements."""
        result = ephem.utc_time_zone(2024, 1, 15, 10, 30, 0.0, 1)
        assert isinstance(result, tuple)
        assert len(result) == 6
        # year, month, day, hour, minute should be int
        assert isinstance(result[0], int)
        assert isinstance(result[1], int)
        assert isinstance(result[2], int)
        assert isinstance(result[3], int)
        assert isinstance(result[4], int)
        # second should be float
        assert isinstance(result[5], float)

    @pytest.mark.unit
    def test_zero_offset(self):
        """UTC+0 should return the same time."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 30, 45.5, 0
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 30
        assert second == pytest.approx(45.5, abs=0.001)

    @pytest.mark.unit
    def test_positive_offset_no_boundary(self):
        """Positive offset without crossing day boundary."""
        # 2024-01-15 10:30:00 UTC -> CET (UTC+1) = 11:30:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 30, 0.0, 1
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 11
        assert minute == 30
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_no_boundary(self):
        """Negative offset without crossing day boundary."""
        # 2024-01-15 15:30:00 UTC -> EST (UTC-5) = 10:30:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 15, 30, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 30
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneDayBoundary:
    """Test day boundary crossings."""

    @pytest.mark.unit
    def test_positive_offset_next_day(self):
        """Positive offset crossing into next day."""
        # 2024-01-15 23:00:00 UTC -> UTC+2 = 2024-01-16 01:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 23, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 1
        assert day == 16
        assert hour == 1
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_previous_day(self):
        """Negative offset crossing into previous day."""
        # 2024-01-15 02:00:00 UTC -> EST (UTC-5) = 2024-01-14 21:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 2, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 14
        assert hour == 21
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneMonthBoundary:
    """Test month boundary crossings."""

    @pytest.mark.unit
    def test_positive_offset_next_month(self):
        """Positive offset crossing into next month."""
        # 2024-01-31 23:00:00 UTC -> UTC+2 = 2024-02-01 01:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 31, 23, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 2
        assert day == 1
        assert hour == 1
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_previous_month(self):
        """Negative offset crossing into previous month."""
        # 2024-02-01 02:00:00 UTC -> EST (UTC-5) = 2024-01-31 21:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 2, 1, 2, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 31
        assert hour == 21
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneYearBoundary:
    """Test year boundary crossings."""

    @pytest.mark.unit
    def test_positive_offset_next_year(self):
        """Positive offset crossing into next year."""
        # 2023-12-31 23:00:00 UTC -> UTC+2 = 2024-01-01 01:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2023, 12, 31, 23, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 1
        assert day == 1
        assert hour == 1
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_previous_year(self):
        """Negative offset crossing into previous year."""
        # 2024-01-01 02:00:00 UTC -> EST (UTC-5) = 2023-12-31 21:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 1, 2, 0, 0.0, -5
        )
        assert year == 2023
        assert month == 12
        assert day == 31
        assert hour == 21
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneFractionalOffsets:
    """Test fractional hour timezone offsets."""

    @pytest.mark.unit
    def test_half_hour_offset(self):
        """Test 30-minute offset (e.g., India UTC+5:30)."""
        # 2024-01-15 10:00:00 UTC -> IST (UTC+5.5) = 15:30:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 0, 0.0, 5.5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 15
        assert minute == 30
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_quarter_hour_offset(self):
        """Test 45-minute offset (e.g., Nepal UTC+5:45)."""
        # 2024-01-15 10:00:00 UTC -> Nepal Time (UTC+5.75) = 15:45:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 0, 0.0, 5.75
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 15
        assert minute == 45
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_half_hour_offset(self):
        """Test negative 30-minute offset (e.g., Newfoundland UTC-3:30)."""
        # 2024-01-15 15:00:00 UTC -> NST (UTC-3.5) = 11:30:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 15, 0, 0.0, -3.5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 11
        assert minute == 30
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneCommonTimezones:
    """Test common real-world timezone conversions."""

    @pytest.mark.unit
    def test_utc_to_cet(self):
        """Test UTC to Central European Time (UTC+1)."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 6, 15, 12, 0, 0.0, 1
        )
        assert hour == 13
        assert minute == 0

    @pytest.mark.unit
    def test_utc_to_est(self):
        """Test UTC to Eastern Standard Time (UTC-5)."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 18, 0, 0.0, -5
        )
        assert hour == 13
        assert minute == 0

    @pytest.mark.unit
    def test_utc_to_pst(self):
        """Test UTC to Pacific Standard Time (UTC-8)."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 20, 0, 0.0, -8
        )
        assert hour == 12
        assert minute == 0

    @pytest.mark.unit
    def test_utc_to_jst(self):
        """Test UTC to Japan Standard Time (UTC+9)."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 3, 0, 0.0, 9
        )
        assert hour == 12
        assert minute == 0


class TestUtcTimeZoneEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.unit
    def test_leap_year_feb29(self):
        """Test Feb 29 in leap year."""
        # 2024-02-29 12:00:00 UTC -> UTC+5 = 17:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 2, 29, 12, 0, 0.0, 5
        )
        assert year == 2024
        assert month == 2
        assert day == 29
        assert hour == 17
        assert minute == 0

    @pytest.mark.unit
    def test_leap_year_crossing_to_mar1(self):
        """Test crossing from Feb 29 to Mar 1 in leap year."""
        # 2024-02-29 23:00:00 UTC -> UTC+2 = 2024-03-01 01:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 2, 29, 23, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 3
        assert day == 1
        assert hour == 1
        assert minute == 0

    @pytest.mark.unit
    def test_fractional_seconds_preserved(self):
        """Test that fractional seconds are preserved."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 30, 45.123, 1
        )
        assert second == pytest.approx(45.123, abs=0.01)

    @pytest.mark.unit
    def test_midnight_utc(self):
        """Test midnight UTC conversion."""
        # 2024-01-15 00:00:00 UTC -> EST (UTC-5) = 2024-01-14 19:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 0, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 14
        assert hour == 19
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_large_positive_offset(self):
        """Test UTC+14 (Line Islands, Kiribati)."""
        # 2024-01-15 10:00:00 UTC -> UTC+14 = 2024-01-16 00:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 0, 0.0, 14
        )
        assert year == 2024
        assert month == 1
        assert day == 16
        assert hour == 0
        assert minute == 0

    @pytest.mark.unit
    def test_large_negative_offset(self):
        """Test UTC-12 (Baker Island)."""
        # 2024-01-15 10:00:00 UTC -> UTC-12 = 2024-01-14 22:00:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 0, 0.0, -12
        )
        assert year == 2024
        assert month == 1
        assert day == 14
        assert hour == 22
        assert minute == 0


class TestUtcTimeZoneConsistency:
    """Test internal consistency."""

    @pytest.mark.unit
    def test_roundtrip_offset(self):
        """Applying +N then -N offset should return original time."""
        orig_year, orig_month, orig_day = 2024, 6, 15
        orig_hour, orig_minute, orig_second = 12, 30, 45.5

        # Apply +5 offset
        y1, m1, d1, h1, min1, s1 = ephem.utc_time_zone(
            orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second, 5
        )

        # Apply -5 offset to get back
        y2, m2, d2, h2, min2, s2 = ephem.utc_time_zone(y1, m1, d1, h1, min1, s1, -5)

        assert y2 == orig_year
        assert m2 == orig_month
        assert d2 == orig_day
        assert h2 == orig_hour
        assert min2 == orig_minute
        assert s2 == pytest.approx(orig_second, abs=0.01)

    @pytest.mark.unit
    def test_consecutive_offsets(self):
        """Test that consecutive integer offsets give consecutive hours."""
        base_hour = 12
        for offset in range(-5, 6):
            year, month, day, hour, minute, second = ephem.utc_time_zone(
                2024, 1, 15, base_hour, 0, 0.0, offset
            )
            expected_hour = (base_hour + offset) % 24
            assert hour == expected_hour, f"Failed for offset {offset}"
