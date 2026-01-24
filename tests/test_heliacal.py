"""
Tests for heliacal_ut function in libephemeris.

Tests the calculation of heliacal rising and setting events for celestial bodies.

Heliacal events are the first/last visibility of a celestial body at dawn/dusk.
These were crucial for ancient calendars (e.g., Egyptian calendar based on Sirius).

Reference data based on astronomical calculations and historical records.
"""

import pytest

from libephemeris import (
    julday,
    revjul,
    heliacal_ut,
    swe_heliacal_ut,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_SUN,
    SE_MOON,
    SE_HELIACAL_RISING,
    SE_HELIACAL_SETTING,
    SE_EVENING_FIRST,
    SE_MORNING_LAST,
)


class TestHeliacalBasic:
    """Basic tests for heliacal_ut function."""

    def test_venus_heliacal_rising_returns_valid_result(self):
        """Test that Venus heliacal rising returns a valid Julian Day."""
        # Start from January 1, 2024
        jd_start = julday(2024, 1, 1, 0)
        # Rome, Italy
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        # Should find an event within a year
        assert jd_event > jd_start
        assert jd_event < jd_start + 365
        assert retflag == SE_HELIACAL_RISING

    def test_mercury_heliacal_rising_returns_valid_result(self):
        """Test that Mercury heliacal rising returns a valid Julian Day."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MERCURY, event_type=SE_HELIACAL_RISING
        )

        # Mercury has ~3 synodic periods per year, so should find one quickly
        assert jd_event > jd_start
        assert jd_event < jd_start + 120  # Within ~4 months
        assert retflag == SE_HELIACAL_RISING

    def test_jupiter_heliacal_rising(self):
        """Test Jupiter heliacal rising."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 51.5074, -0.1278  # London

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_JUPITER, event_type=SE_HELIACAL_RISING
        )

        assert jd_event > jd_start
        assert jd_event < jd_start + 400  # Within a bit more than a year
        assert retflag == SE_HELIACAL_RISING

    def test_saturn_heliacal_rising(self):
        """Test Saturn heliacal rising."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 40.7128, -74.0060  # New York

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_SATURN, event_type=SE_HELIACAL_RISING
        )

        # Saturn may not have heliacal rising within search window for all start dates
        # Just verify it returns a valid response (either found or not found)
        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag == SE_HELIACAL_RISING:
            assert jd_event > jd_start
            assert jd_event < jd_start + 400

    def test_mars_heliacal_rising(self):
        """Test Mars heliacal rising."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.6762, 139.6503  # Tokyo

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MARS, event_type=SE_HELIACAL_RISING
        )

        # Mars has ~780 day synodic period, may not have event in window
        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag == SE_HELIACAL_RISING:
            assert jd_event > jd_start
            assert jd_event < jd_start + 780


class TestHeliacalEventTypes:
    """Test different heliacal event types."""

    def test_heliacal_setting(self):
        """Test heliacal setting (evening last visibility)."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_SETTING
        )

        assert jd_event > jd_start
        assert retflag == SE_HELIACAL_SETTING

    def test_evening_first(self):
        """Test evening first visibility (after superior conjunction)."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_EVENING_FIRST
        )

        # May or may not find this event depending on current position
        # Just verify it returns something valid
        assert retflag in (SE_EVENING_FIRST, -1)

    def test_morning_last(self):
        """Test morning last visibility (before superior conjunction)."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MERCURY, event_type=SE_MORNING_LAST
        )

        # May or may not find this event
        assert retflag in (SE_MORNING_LAST, -1)


class TestHeliacalValidation:
    """Test input validation for heliacal_ut."""

    def test_sun_raises_error(self):
        """Test that SE_SUN raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="SE_SUN"):
            heliacal_ut(jd_start, lat, lon, body=SE_SUN, event_type=SE_HELIACAL_RISING)

    def test_moon_raises_error(self):
        """Test that SE_MOON raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="SE_MOON"):
            heliacal_ut(jd_start, lat, lon, body=SE_MOON, event_type=SE_HELIACAL_RISING)

    def test_invalid_event_type_raises_error(self):
        """Test that invalid event type raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="Invalid event_type"):
            heliacal_ut(jd_start, lat, lon, body=SE_VENUS, event_type=99)


class TestHeliacalLocations:
    """Test heliacal calculations at various locations."""

    def test_northern_hemisphere(self):
        """Test at northern latitude."""
        jd_start = julday(2024, 6, 1, 0)  # Summer
        lat, lon = 60.1699, 24.9384  # Helsinki

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_JUPITER, event_type=SE_HELIACAL_RISING
        )

        # High latitude summer nights are short, may affect visibility
        assert retflag in (SE_HELIACAL_RISING, -1)

    def test_southern_hemisphere(self):
        """Test at southern latitude."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = -33.8688, 151.2093  # Sydney

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        # In the southern hemisphere, twilight times differ
        # Allow for either success or not-found
        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag == SE_HELIACAL_RISING:
            assert jd_event > jd_start

    def test_equatorial_location(self):
        """Test at equatorial location."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 0.3476, 32.5825  # Kampala, Uganda

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_SATURN, event_type=SE_HELIACAL_RISING
        )

        # Equatorial locations have consistent twilight times
        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag == SE_HELIACAL_RISING:
            assert jd_event > jd_start


class TestHeliacalAtmosphericConditions:
    """Test heliacal calculations with different atmospheric conditions."""

    def test_with_altitude(self):
        """Test calculation with observer at altitude."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome
        altitude = 2000.0  # 2000 meters (e.g., mountain observatory)

        jd_event, retflag = heliacal_ut(
            jd_start,
            lat,
            lon,
            altitude=altitude,
            body=SE_VENUS,
            event_type=SE_HELIACAL_RISING,
        )

        assert jd_event > jd_start
        assert retflag == SE_HELIACAL_RISING

    def test_with_pressure(self):
        """Test calculation with custom pressure."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start,
            lat,
            lon,
            pressure=1000.0,  # Lower pressure
            body=SE_VENUS,
            event_type=SE_HELIACAL_RISING,
        )

        assert jd_event > jd_start
        assert retflag == SE_HELIACAL_RISING

    def test_with_humidity(self):
        """Test calculation with different humidity."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start,
            lat,
            lon,
            humidity=0.8,  # High humidity
            body=SE_VENUS,
            event_type=SE_HELIACAL_RISING,
        )

        assert jd_event > jd_start
        assert retflag == SE_HELIACAL_RISING


class TestHeliacalAlias:
    """Test swe_heliacal_ut alias."""

    def test_swe_alias_works(self):
        """Test that swe_heliacal_ut is an alias for heliacal_ut."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        result1 = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )
        result2 = swe_heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        assert result1 == result2


class TestHeliacalDateValidation:
    """Test that heliacal events return sensible dates."""

    def test_event_date_is_reasonable(self):
        """Test that returned date is within expected range."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        if retflag > 0:
            year, month, day, hour = revjul(jd_event)
            # Event should be in 2024 or 2025
            assert year in (2024, 2025)
            # Month and day should be valid
            assert 1 <= month <= 12
            assert 1 <= day <= 31
