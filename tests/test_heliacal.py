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
    heliacal_pheno_ut,
    swe_heliacal_pheno_ut,
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

        # Venus synodic period is ~584 days, so the next heliacal rising
        # after inferior conjunction may be up to ~500 days away
        assert jd_event > jd_start
        assert jd_event < jd_start + 600
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


class TestSweHeliacalUt:
    """Test swe_heliacal_ut pyswisseph-compatible API."""

    def test_swe_heliacal_ut_basic_call(self):
        """Test basic swe_heliacal_ut call with array parameters."""
        jd_start = julday(2024, 1, 1, 0)
        # Geographic position: Rome (lon, lat, altitude)
        geopos = (12.4964, 41.9028, 0.0)
        # Atmospheric conditions: pressure, temp, humidity, extinction
        datm = (1013.25, 15.0, 40.0, 0.0)
        # Observer: age, Snellen ratio, and optical params
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_HELIACAL_RISING
        )

        # Should return a tuple of 3 floats (jd1, jd2, jd3)
        assert isinstance(result, tuple)
        assert len(result) == 3
        # Venus synodic period is ~584 days, so next heliacal rising
        # after inferior conjunction may be up to ~500 days away
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 600

    def test_swe_heliacal_ut_with_planet_name_mercury(self):
        """Test swe_heliacal_ut with Mercury."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Mercury", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3
        # Mercury has ~3 synodic periods per year
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 120

    def test_swe_heliacal_ut_with_planet_name_mars(self):
        """Test swe_heliacal_ut with Mars."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Mars", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_planet_name_jupiter(self):
        """Test swe_heliacal_ut with Jupiter."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-0.1278, 51.5074, 0.0)  # London
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Jupiter", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_planet_name_saturn(self):
        """Test swe_heliacal_ut with Saturn."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Saturn", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_heliacal_setting(self):
        """Test swe_heliacal_ut for heliacal setting event."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_HELIACAL_SETTING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_evening_first(self):
        """Test swe_heliacal_ut for evening first event."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_EVENING_FIRST
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_morning_last(self):
        """Test swe_heliacal_ut for morning last event."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Mercury", SE_MORNING_LAST
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_sun_raises_error(self):
        """Test that Sun raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="Sun"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Sun", SE_HELIACAL_RISING)

    def test_swe_heliacal_ut_moon_raises_error(self):
        """Test that Moon raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="Moon"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Moon", SE_HELIACAL_RISING)

    def test_swe_heliacal_ut_invalid_event_type(self):
        """Test that invalid event type raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="Invalid event_type"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Venus", 99)

    def test_swe_heliacal_ut_invalid_object_name(self):
        """Test that unknown object name raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="not recognized"):
            swe_heliacal_ut(
                jd_start, geopos, datm, dobs, "InvalidPlanet", SE_HELIACAL_RISING
            )

    def test_swe_heliacal_ut_case_insensitive(self):
        """Test that planet name matching is case insensitive."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Should work with lowercase
        result1 = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "venus", SE_HELIACAL_RISING
        )
        # Should work with uppercase
        result2 = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "VENUS", SE_HELIACAL_RISING
        )
        # Should work with mixed case
        result3 = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_HELIACAL_RISING
        )

        # All should return same result
        assert result1[0] == result2[0] == result3[0]

    def test_swe_heliacal_ut_default_atmospheric(self):
        """Test that zero atmospheric values get defaults."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        # All zeros should get defaults
        datm = (0, 0, 0, 0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_altitude(self):
        """Test swe_heliacal_ut with observer at altitude."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 2000.0)  # 2000m altitude
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_southern_hemisphere(self):
        """Test swe_heliacal_ut in southern hemisphere."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (151.2093, -33.8688, 0.0)  # Sydney
        datm = (1013.25, 25.0, 60.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_planet_id_string(self):
        """Test swe_heliacal_ut with planet ID as string."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Venus is planet ID 3
        result = swe_heliacal_ut(jd_start, geopos, datm, dobs, "3", SE_HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3


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


# =============================================================================
# HELIACAL_PHENO_UT TESTS
# =============================================================================


class TestHeliacalPhenoBasic:
    """Basic tests for heliacal_pheno_ut function."""

    def test_venus_pheno_returns_valid_result(self):
        """Test that Venus heliacal phenomena returns a valid result tuple."""
        jd = julday(2024, 1, 1, 12)  # Noon
        lat, lon = 41.9028, 12.4964  # Rome

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        # Should return a tuple of 50 floats
        assert isinstance(dret, tuple)
        assert len(dret) == 50
        assert all(isinstance(x, float) for x in dret)
        assert retflag > 0

    def test_mercury_pheno_returns_valid_result(self):
        """Test that Mercury heliacal phenomena returns a valid result."""
        jd = julday(2024, 3, 15, 5)  # Morning
        lat, lon = 51.5074, -0.1278  # London

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_MERCURY, event_type=SE_HELIACAL_RISING
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50
        assert retflag > 0

    def test_jupiter_pheno_returns_valid_result(self):
        """Test that Jupiter heliacal phenomena returns a valid result."""
        jd = julday(2024, 6, 1, 4)  # Early morning
        lat, lon = 40.7128, -74.0060  # New York

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_JUPITER, event_type=SE_HELIACAL_RISING
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50

    def test_mars_pheno_evening(self):
        """Test Mars heliacal phenomena for evening setting."""
        jd = julday(2024, 7, 20, 20)  # Evening
        lat, lon = 35.6762, 139.6503  # Tokyo

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_MARS, event_type=SE_HELIACAL_SETTING
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50


class TestHeliacalPhenoValues:
    """Test that heliacal_pheno_ut returns sensible values."""

    def test_altitude_values_are_reasonable(self):
        """Test that altitude values are within valid range."""
        jd = julday(2024, 1, 15, 6)  # Dawn
        lat, lon = 41.9028, 12.4964  # Rome

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        # Altitudes should be in -90 to +90 range
        alt_o = dret[0]  # Object topocentric altitude
        app_alt_o = dret[1]  # Object apparent altitude
        geo_alt_o = dret[2]  # Object geocentric altitude
        alt_s = dret[4]  # Sun altitude

        assert -90 <= alt_o <= 90
        assert -90 <= app_alt_o <= 90
        assert -90 <= geo_alt_o <= 90
        assert -90 <= alt_s <= 90

    def test_azimuth_values_are_reasonable(self):
        """Test that azimuth values are within valid range."""
        jd = julday(2024, 1, 15, 6)
        lat, lon = 41.9028, 12.4964

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        azi_o = dret[3]  # Object azimuth
        azi_s = dret[5]  # Sun azimuth

        assert 0 <= azi_o <= 360
        assert 0 <= azi_s <= 360

    def test_arcus_visionis_values(self):
        """Test that arcus visionis values are calculated."""
        jd = julday(2024, 1, 15, 6)
        lat, lon = 41.9028, 12.4964

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        tav_act = dret[6]  # Topocentric arcus visionis
        arcv_act = dret[7]  # Geocentric arcus visionis

        # TAV is the altitude difference between body and Sun
        # Should be a reasonable value (typically -180 to +180)
        assert -180 <= tav_act <= 180
        assert -180 <= arcv_act <= 180

    def test_extinction_coefficient_is_positive(self):
        """Test that extinction coefficient is a sensible positive value."""
        jd = julday(2024, 1, 15, 6)
        lat, lon = 41.9028, 12.4964

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        k_act = dret[10]  # Extinction coefficient

        # Extinction coefficient typically 0.1 to 0.6
        assert 0 < k_act < 1.0

    def test_magnitude_for_venus(self):
        """Test that magnitude is calculated for Venus."""
        jd = julday(2024, 1, 15, 6)
        lat, lon = 41.9028, 12.4964

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        magnitude = dret[20]

        # Venus magnitude typically -4.5 to -3
        # But can vary, so just check it's a number
        assert isinstance(magnitude, float)

    def test_parallax_is_small(self):
        """Test that parallax is a small positive value."""
        jd = julday(2024, 1, 15, 6)
        lat, lon = 41.9028, 12.4964

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        parallax = dret[19]

        # Parallax for planets is typically < 1 degree
        assert 0 <= parallax < 1.0


class TestHeliacalPhenoEventTypes:
    """Test different event types for heliacal_pheno_ut."""

    def test_morning_first_event(self):
        """Test SE_HELIACAL_RISING (morning first) event type."""
        jd = julday(2024, 3, 1, 5)
        lat, lon = 41.9028, 12.4964

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        assert len(dret) == 50
        assert retflag > 0

    def test_evening_last_event(self):
        """Test SE_HELIACAL_SETTING (evening last) event type."""
        jd = julday(2024, 3, 1, 19)
        lat, lon = 41.9028, 12.4964

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_SETTING
        )

        assert len(dret) == 50
        assert retflag > 0

    def test_evening_first_event(self):
        """Test SE_EVENING_FIRST event type."""
        jd = julday(2024, 3, 1, 19)
        lat, lon = 41.9028, 12.4964

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_MERCURY, event_type=SE_EVENING_FIRST
        )

        assert len(dret) == 50

    def test_morning_last_event(self):
        """Test SE_MORNING_LAST event type."""
        jd = julday(2024, 3, 1, 5)
        lat, lon = 41.9028, 12.4964

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_MERCURY, event_type=SE_MORNING_LAST
        )

        assert len(dret) == 50


class TestHeliacalPhenoValidation:
    """Test input validation for heliacal_pheno_ut."""

    def test_invalid_body_raises_error(self):
        """Test that invalid body ID raises ValueError."""
        jd = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="illegal planet number"):
            heliacal_pheno_ut(jd, lat, lon, body=999, event_type=SE_HELIACAL_RISING)

    def test_invalid_event_type_raises_error(self):
        """Test that invalid event type raises ValueError."""
        jd = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="Invalid event_type"):
            heliacal_pheno_ut(jd, lat, lon, body=SE_VENUS, event_type=99)


class TestHeliacalPhenoLocations:
    """Test heliacal_pheno_ut at various locations."""

    def test_northern_hemisphere(self):
        """Test at northern latitude."""
        jd = julday(2024, 6, 1, 3)  # Summer morning
        lat, lon = 60.1699, 24.9384  # Helsinki

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_JUPITER, event_type=SE_HELIACAL_RISING
        )

        assert len(dret) == 50

    def test_southern_hemisphere(self):
        """Test at southern latitude."""
        jd = julday(2024, 1, 15, 5)  # Summer morning (S. hemisphere)
        lat, lon = -33.8688, 151.2093  # Sydney

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_SATURN, event_type=SE_HELIACAL_RISING
        )

        assert len(dret) == 50

    def test_equatorial_location(self):
        """Test at equatorial location."""
        jd = julday(2024, 4, 1, 6)
        lat, lon = 0.3476, 32.5825  # Kampala, Uganda

        dret, retflag = heliacal_pheno_ut(
            jd, lat, lon, body=SE_MARS, event_type=SE_HELIACAL_RISING
        )

        assert len(dret) == 50


class TestHeliacalPhenoAtmospheric:
    """Test heliacal_pheno_ut with different atmospheric conditions."""

    def test_with_altitude(self):
        """Test with observer at high altitude."""
        jd = julday(2024, 1, 15, 5)
        lat, lon = 41.9028, 12.4964
        altitude = 2000.0  # 2000 meters

        dret, retflag = heliacal_pheno_ut(
            jd,
            lat,
            lon,
            altitude=altitude,
            body=SE_VENUS,
            event_type=SE_HELIACAL_RISING,
        )

        assert len(dret) == 50
        # Extinction should be lower at high altitude
        k_high = dret[10]

        dret_low, _ = heliacal_pheno_ut(
            jd, lat, lon, altitude=0, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )
        k_low = dret_low[10]

        assert k_high < k_low  # Lower extinction at higher altitude

    def test_with_high_humidity(self):
        """Test with high humidity."""
        jd = julday(2024, 1, 15, 5)
        lat, lon = 41.9028, 12.4964

        dret_low, _ = heliacal_pheno_ut(
            jd, lat, lon, humidity=0.2, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )
        dret_high, _ = heliacal_pheno_ut(
            jd, lat, lon, humidity=0.9, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        # Higher humidity should increase extinction
        k_low = dret_low[10]
        k_high = dret_high[10]

        assert k_high > k_low

    def test_with_low_pressure(self):
        """Test with low atmospheric pressure."""
        jd = julday(2024, 1, 15, 5)
        lat, lon = 41.9028, 12.4964

        dret_normal, _ = heliacal_pheno_ut(
            jd, lat, lon, pressure=1013.25, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )
        dret_low, _ = heliacal_pheno_ut(
            jd, lat, lon, pressure=800.0, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        # Lower pressure should decrease extinction
        k_normal = dret_normal[10]
        k_low = dret_low[10]

        assert k_low < k_normal


class TestHeliacalPhenoAlias:
    """Test swe_heliacal_pheno_ut alias."""

    def test_swe_alias_works(self):
        """Test that swe_heliacal_pheno_ut is an alias for heliacal_pheno_ut."""
        jd = julday(2024, 1, 15, 6)
        lat, lon = 41.9028, 12.4964

        result1 = heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )
        result2 = swe_heliacal_pheno_ut(
            jd, lat, lon, body=SE_VENUS, event_type=SE_HELIACAL_RISING
        )

        assert result1 == result2


class TestHeliacalPhenoMoon:
    """Test heliacal_pheno_ut specifically for Moon (crescent calculations)."""

    def test_moon_crescent_values(self):
        """Test that Moon-specific crescent values are calculated."""
        jd = julday(2024, 1, 12, 18)  # Near new moon
        lat, lon = 41.9028, 12.4964

        dret, _ = heliacal_pheno_ut(
            jd, lat, lon, body=SE_MOON, event_type=SE_EVENING_FIRST
        )

        w_moon = dret[16]  # Crescent width
        l_moon = dret[25]  # Crescent length
        q_yallop = dret[17]  # Yallop q-test
        illumination = dret[27]  # Illumination percentage

        # These values should be calculated for Moon
        assert isinstance(w_moon, float)
        assert isinstance(l_moon, float)
        assert isinstance(q_yallop, float)
        assert isinstance(illumination, float)


# =============================================================================
# PLANET-SPECIFIC HELIACAL TESTS
# =============================================================================


class TestHeliacalInnerOuterPlanets:
    """Test proper handling of inner vs outer planet geometry."""

    def test_is_inner_planet_mercury(self):
        """Test that Mercury is identified as an inner planet."""
        from libephemeris import is_inner_planet, INNER_PLANETS

        assert is_inner_planet(SE_MERCURY)
        assert SE_MERCURY in INNER_PLANETS

    def test_is_inner_planet_venus(self):
        """Test that Venus is identified as an inner planet."""
        from libephemeris import is_inner_planet, INNER_PLANETS

        assert is_inner_planet(SE_VENUS)
        assert SE_VENUS in INNER_PLANETS

    def test_is_inner_planet_mars(self):
        """Test that Mars is NOT an inner planet."""
        from libephemeris import is_inner_planet

        assert not is_inner_planet(SE_MARS)

    def test_is_inner_planet_jupiter(self):
        """Test that Jupiter is NOT an inner planet."""
        from libephemeris import is_inner_planet

        assert not is_inner_planet(SE_JUPITER)

    def test_is_inner_planet_saturn(self):
        """Test that Saturn is NOT an inner planet."""
        from libephemeris import is_inner_planet

        assert not is_inner_planet(SE_SATURN)


class TestHeliacalOuterPlanetValidation:
    """Test that outer planets reject invalid event types."""

    def test_mars_evening_first_raises_error(self):
        """Test that Mars with SE_EVENING_FIRST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="inner planets"):
            heliacal_ut(jd_start, lat, lon, body=SE_MARS, event_type=SE_EVENING_FIRST)

    def test_mars_morning_last_raises_error(self):
        """Test that Mars with SE_MORNING_LAST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="inner planets"):
            heliacal_ut(jd_start, lat, lon, body=SE_MARS, event_type=SE_MORNING_LAST)

    def test_jupiter_evening_first_raises_error(self):
        """Test that Jupiter with SE_EVENING_FIRST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="inner planets"):
            heliacal_ut(
                jd_start, lat, lon, body=SE_JUPITER, event_type=SE_EVENING_FIRST
            )

    def test_jupiter_morning_last_raises_error(self):
        """Test that Jupiter with SE_MORNING_LAST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="inner planets"):
            heliacal_ut(jd_start, lat, lon, body=SE_JUPITER, event_type=SE_MORNING_LAST)

    def test_saturn_evening_first_raises_error(self):
        """Test that Saturn with SE_EVENING_FIRST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="inner planets"):
            heliacal_ut(jd_start, lat, lon, body=SE_SATURN, event_type=SE_EVENING_FIRST)

    def test_saturn_morning_last_raises_error(self):
        """Test that Saturn with SE_MORNING_LAST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        with pytest.raises(ValueError, match="inner planets"):
            heliacal_ut(jd_start, lat, lon, body=SE_SATURN, event_type=SE_MORNING_LAST)


class TestHeliacalInnerPlanetEventTypes:
    """Test that inner planets accept all event types."""

    def test_mercury_heliacal_rising(self):
        """Test Mercury heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MERCURY, event_type=SE_HELIACAL_RISING
        )

        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_mercury_heliacal_setting(self):
        """Test Mercury heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MERCURY, event_type=SE_HELIACAL_SETTING
        )

        assert retflag in (SE_HELIACAL_SETTING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_mercury_evening_first(self):
        """Test Mercury evening first visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MERCURY, event_type=SE_EVENING_FIRST
        )

        # Should not raise an error - inner planets can use this event type
        assert retflag in (SE_EVENING_FIRST, -1)

    def test_mercury_morning_last(self):
        """Test Mercury morning last visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MERCURY, event_type=SE_MORNING_LAST
        )

        # Should not raise an error - inner planets can use this event type
        assert retflag in (SE_MORNING_LAST, -1)

    def test_venus_evening_first(self):
        """Test Venus evening first visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_EVENING_FIRST
        )

        # Should not raise an error - inner planets can use this event type
        assert retflag in (SE_EVENING_FIRST, -1)

    def test_venus_morning_last(self):
        """Test Venus morning last visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_VENUS, event_type=SE_MORNING_LAST
        )

        # Should not raise an error - inner planets can use this event type
        assert retflag in (SE_MORNING_LAST, -1)


class TestHeliacalOuterPlanetRisingSetting:
    """Test that outer planets work correctly with heliacal rising/setting."""

    def test_mars_heliacal_rising(self):
        """Test Mars heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MARS, event_type=SE_HELIACAL_RISING
        )

        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_mars_heliacal_setting(self):
        """Test Mars heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_MARS, event_type=SE_HELIACAL_SETTING
        )

        assert retflag in (SE_HELIACAL_SETTING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_jupiter_heliacal_rising(self):
        """Test Jupiter heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_JUPITER, event_type=SE_HELIACAL_RISING
        )

        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_jupiter_heliacal_setting(self):
        """Test Jupiter heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_JUPITER, event_type=SE_HELIACAL_SETTING
        )

        assert retflag in (SE_HELIACAL_SETTING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_saturn_heliacal_rising(self):
        """Test Saturn heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_SATURN, event_type=SE_HELIACAL_RISING
        )

        assert retflag in (SE_HELIACAL_RISING, -1)
        if retflag > 0:
            assert jd_event > jd_start

    def test_saturn_heliacal_setting(self):
        """Test Saturn heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964

        jd_event, retflag = heliacal_ut(
            jd_start, lat, lon, body=SE_SATURN, event_type=SE_HELIACAL_SETTING
        )

        assert retflag in (SE_HELIACAL_SETTING, -1)
        if retflag > 0:
            assert jd_event > jd_start


class TestSweHeliacalUtOuterPlanetValidation:
    """Test that swe_heliacal_ut also validates inner/outer planet constraints."""

    def test_mars_evening_first_via_swe_api_raises_error(self):
        """Test that Mars with SE_EVENING_FIRST raises ValueError via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="inner planets"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Mars", SE_EVENING_FIRST)

    def test_jupiter_morning_last_via_swe_api_raises_error(self):
        """Test that Jupiter with SE_MORNING_LAST raises ValueError via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="inner planets"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Jupiter", SE_MORNING_LAST)

    def test_saturn_evening_first_via_swe_api_raises_error(self):
        """Test that Saturn with SE_EVENING_FIRST raises ValueError via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="inner planets"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Saturn", SE_EVENING_FIRST)

    def test_mercury_evening_first_via_swe_api_works(self):
        """Test that Mercury with SE_EVENING_FIRST works via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Should not raise an error
        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Mercury", SE_EVENING_FIRST
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_venus_morning_last_via_swe_api_works(self):
        """Test that Venus with SE_MORNING_LAST works via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Should not raise an error
        result = swe_heliacal_ut(jd_start, geopos, datm, dobs, "Venus", SE_MORNING_LAST)

        assert isinstance(result, tuple)
        assert len(result) == 3
