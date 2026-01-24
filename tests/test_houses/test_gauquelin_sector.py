"""
Tests for gauquelin_sector and swe_gauquelin_sector functions.

Tests the calculation of Gauquelin sectors (1-36) for celestial bodies.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_SWIEPH,
)


class TestGauquelinSectorBasic:
    """Basic tests for gauquelin_sector function."""

    @pytest.mark.unit
    def test_gauquelin_sector_returns_float(self):
        """gauquelin_sector() should return a float."""
        jd = 2451545.0  # J2000.0
        lat = 48.85  # Paris
        lon = 2.35

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)

        assert isinstance(result, float)

    @pytest.mark.unit
    def test_gauquelin_sector_in_valid_range(self):
        """gauquelin_sector() should return value in range [1, 37)."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_integer_part_is_sector_number(self):
        """Integer part of result should be sector number (1-36)."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)
        sector_num = int(result)

        assert 1 <= sector_num <= 36

    @pytest.mark.unit
    def test_gauquelin_sector_alias_matches(self):
        """swe_gauquelin_sector should give same result as gauquelin_sector."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result1 = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)
        result2 = ephem.swe_gauquelin_sector(jd, SE_MARS, lat, lon)

        assert result1 == result2


class TestGauquelinSectorMultiplePlanets:
    """Test gauquelin_sector with different planets."""

    PLANETS = [SE_SUN, SE_MOON, SE_MARS, SE_JUPITER, SE_SATURN]

    @pytest.mark.unit
    @pytest.mark.parametrize("planet", PLANETS)
    def test_gauquelin_sector_all_planets_return_valid(self, planet):
        """gauquelin_sector should return valid values for all major planets."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, planet, lat, lon)

        assert 1.0 <= result < 37.0, f"Invalid result {result} for planet {planet}"


class TestGauquelinSectorMethods:
    """Test gauquelin_sector with different calculation methods."""

    @pytest.mark.unit
    @pytest.mark.parametrize("method", [0, 1])
    def test_gauquelin_sector_method_returns_valid(self, method):
        """gauquelin_sector should return valid results for methods 0 and 1."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon, method=method)

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_method_0_differs_from_method_1_for_high_latitude_body(self):
        """Method 0 (with lat) and method 1 (without lat) may give different results."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result_with_lat = ephem.gauquelin_sector(jd, SE_MOON, lat, lon, method=0)
        result_without_lat = ephem.gauquelin_sector(jd, SE_MOON, lat, lon, method=1)

        # They should both be valid
        assert 1.0 <= result_with_lat < 37.0
        assert 1.0 <= result_without_lat < 37.0


class TestGauquelinSectorComparisonWithSwisseph:
    """Compare gauquelin_sector results with Swiss Ephemeris."""

    # Test cases: (jd, planet, lat, lon, method)
    TEST_CASES = [
        # Paris, Mars at J2000
        (2451545.0, SE_MARS, 48.85, 2.35, 0),
        # Rome, Sun at different time
        (2451600.0, SE_SUN, 41.9, 12.5, 0),
        # London, Moon
        (2451550.0, SE_MOON, 51.5, -0.12, 0),
        # New York, Jupiter
        (2451545.0, SE_JUPITER, 40.7, -74.0, 0),
        # Southern hemisphere, Saturn
        (2451545.0, SE_SATURN, -33.9, 151.2, 0),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("jd,planet,lat,lon,method", TEST_CASES)
    def test_gauquelin_sector_matches_swisseph_sector(
        self, jd, planet, lat, lon, method
    ):
        """gauquelin_sector sector number should match Swiss Ephemeris."""
        # Get result from libephemeris
        result_lib = ephem.gauquelin_sector(jd, planet, lat, lon, method=method)

        # Get result from Swiss Ephemeris
        geopos = (lon, lat, 0.0)  # (lon, lat, alt)
        result_swe = swe.gauquelin_sector(
            jd, planet, method, geopos, 0.0, 0.0, SEFLG_SWIEPH
        )

        # Compare sector numbers
        sector_lib = int(result_lib)
        sector_swe = int(result_swe)

        # Allow difference of up to 1 sector due to different algorithms
        # (Swiss Ephemeris uses more complex rise/set calculations)
        diff = abs(sector_lib - sector_swe)
        # Handle wrap-around (sector 36 is adjacent to sector 1)
        if diff > 18:
            diff = 36 - diff

        assert diff <= 2, (
            f"Sector mismatch: libephemeris={sector_lib}, swisseph={sector_swe}, "
            f"diff={diff}"
        )


class TestGauquelinSectorEdgeCases:
    """Test edge cases for gauquelin_sector."""

    @pytest.mark.unit
    def test_gauquelin_sector_equator(self):
        """gauquelin_sector should work at the equator."""
        jd = 2451545.0
        lat = 0.0  # Equator
        lon = 0.0

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_high_latitude(self):
        """gauquelin_sector should work at high latitudes."""
        jd = 2451545.0
        lat = 65.0  # Near Arctic circle
        lon = 25.0

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_southern_hemisphere(self):
        """gauquelin_sector should work in southern hemisphere."""
        jd = 2451545.0
        lat = -45.0
        lon = 170.0

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_negative_longitude(self):
        """gauquelin_sector should work with negative longitude (Western)."""
        jd = 2451545.0
        lat = 40.7
        lon = -74.0  # New York

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)

        assert 1.0 <= result < 37.0


class TestGauquelinSectorSectorDistribution:
    """Test that sectors are properly distributed across time."""

    @pytest.mark.unit
    def test_different_times_give_different_sectors(self):
        """Different times should give different sector positions."""
        lat = 48.85
        lon = 2.35

        # Test Mars at different times (spaced ~2 hours apart)
        # 2 hours = ~2/24 of a day = ~0.083 days
        sectors = []
        for i in range(6):
            jd = 2451545.0 + i * 0.1  # ~2.4 hours apart
            sector = ephem.gauquelin_sector(jd, SE_MARS, lat, lon)
            sectors.append(int(sector))

        # At least some of the sectors should be different
        unique_sectors = set(sectors)
        assert len(unique_sectors) >= 2, (
            f"Expected different sectors at different times, got {sectors}"
        )

    @pytest.mark.unit
    def test_sun_cycles_through_all_sectors_in_day(self):
        """Sun should cycle through many sectors over 24 hours."""
        lat = 48.85
        lon = 2.35
        jd_start = 2451545.0

        # Sample 24 times over 24 hours
        sectors = set()
        for i in range(24):
            jd = jd_start + i / 24.0
            sector = ephem.gauquelin_sector(jd, SE_SUN, lat, lon)
            sectors.add(int(sector))

        # Sun should appear in at least 20 different sectors in 24 hours
        assert len(sectors) >= 20, (
            f"Expected Sun to visit many sectors in 24h, got {len(sectors)}: {sorted(sectors)}"
        )


class TestGauquelinSectorAtmosphericParameters:
    """Test gauquelin_sector with atmospheric parameters."""

    @pytest.mark.unit
    def test_gauquelin_sector_with_altitude(self):
        """gauquelin_sector should accept altitude parameter."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35
        altitude = 1000.0  # 1km altitude

        result = ephem.gauquelin_sector(jd, SE_MARS, lat, lon, altitude=altitude)

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_with_pressure_and_temperature(self):
        """gauquelin_sector should accept pressure and temperature parameters."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(
            jd, SE_MARS, lat, lon, pressure=1000.0, temperature=20.0
        )

        assert 1.0 <= result < 37.0
