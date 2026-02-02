"""
Tests for lunar eclipse location-based calculations in libephemeris.

Tests the lun_eclipse_when_loc function which finds lunar eclipses
visible from a specific geographic location.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

from libephemeris import (
    julday,
    revjul,
    lun_eclipse_when_loc,
    swe_lun_eclipse_when_loc,
    SE_ECL_TOTAL,
    SE_ECL_VISIBLE,
    SE_ECL_MAX_VISIBLE,
    SE_ECL_1ST_VISIBLE,
    SE_ECL_2ND_VISIBLE,
    SE_ECL_3RD_VISIBLE,
    SE_ECL_4TH_VISIBLE,
)


class TestLunEclipseWhenLoc:
    """Test suite for lun_eclipse_when_loc function."""

    def test_finds_lunar_eclipse_from_rome(self):
        """Test that function finds a lunar eclipse visible from Rome."""
        # Start from Jan 1, 2024
        jd_start = julday(2024, 1, 1, 0)
        rome_lat, rome_lon = 41.9028, 12.4964

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, rome_lat, rome_lon)

        # Should find an eclipse
        assert ecl_type != 0
        assert times[0] > jd_start  # Maximum should be after start
        assert times[0] > 0  # Should have valid maximum time
        # Should be marked as visible
        assert ecl_type & SE_ECL_VISIBLE

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2024, 1, 1, 0)
        london_lat, london_lon = 51.5074, -0.1278

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, london_lat, london_lon)

        # Should return 10-element times tuple
        assert len(times) == 10
        # All elements should be floats
        assert all(isinstance(t, float) for t in times)

        # Should return 11-element attr tuple
        assert len(attr) == 11
        # All elements should be floats
        assert all(isinstance(a, float) for a in attr)

        # Eclipse type should be int
        assert isinstance(ecl_type, int)

    def test_attributes_have_valid_values(self):
        """Test that eclipse attributes have reasonable values."""
        jd_start = julday(2024, 1, 1, 0)
        tokyo_lat, tokyo_lon = 35.6762, 139.6503

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, tokyo_lat, tokyo_lon)

        # Umbral magnitude should be between 0 and ~2 for visible eclipses
        assert 0 <= attr[0] <= 2.5

        # Penumbral magnitude should be positive
        assert attr[1] >= 0

        # Azimuth should be 0-360
        assert 0 <= attr[3] <= 360

        # Altitude could be any value but should be reasonable
        assert -90 <= attr[4] <= 90

        # Moon diameter should be approximately 0.5 degrees
        assert 0.4 < attr[5] < 0.6

    def test_phase_times_ordering(self):
        """Test that phase times are in correct chronological order."""
        jd_start = julday(2024, 1, 1, 0)
        sydney_lat, sydney_lon = -33.8688, 151.2093

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, sydney_lat, sydney_lon)

        # Penumbral begin should be first (if present)
        if times[5] > 0:
            assert times[5] < times[0]  # Penumbral begin < maximum

        # Penumbral end should be last (if present)
        if times[6] > 0:
            assert times[6] > times[0]  # Penumbral end > maximum

        # Partial phases should be between penumbral phases
        if times[1] > 0 and times[5] > 0:
            assert times[1] >= times[5]  # Partial begin >= penumbral begin
        if times[4] > 0 and times[6] > 0:
            assert times[4] <= times[6]  # Partial end <= penumbral end

    def test_different_locations_may_see_different_eclipses(self):
        """Test that different locations might see different eclipses."""
        jd_start = julday(2024, 1, 1, 0)

        # Two locations on opposite sides of the Earth
        new_york_lat, new_york_lon = 40.7128, -74.0060
        beijing_lat, beijing_lon = 39.9042, 116.4074

        ecl_type_ny, times_ny, attr_ny = lun_eclipse_when_loc(
            jd_start, new_york_lat, new_york_lon
        )
        ecl_type_bj, times_bj, attr_bj = lun_eclipse_when_loc(
            jd_start, beijing_lat, beijing_lon
        )

        # Both should find eclipses (may be same or different)
        assert ecl_type_ny & SE_ECL_VISIBLE
        assert ecl_type_bj & SE_ECL_VISIBLE

    def test_visibility_flags_are_set(self):
        """Test that visibility flags are properly set."""
        jd_start = julday(2024, 1, 1, 0)
        berlin_lat, berlin_lon = 52.5200, 13.4050

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, berlin_lat, berlin_lon)

        # SE_ECL_VISIBLE should always be set for returned eclipses
        assert ecl_type & SE_ECL_VISIBLE

        # At least some visibility flags should be set
        visibility_flags = (
            SE_ECL_MAX_VISIBLE
            | SE_ECL_1ST_VISIBLE
            | SE_ECL_2ND_VISIBLE
            | SE_ECL_3RD_VISIBLE
            | SE_ECL_4TH_VISIBLE
        )
        # At least one phase should be visible
        assert ecl_type & visibility_flags

    def test_moon_altitude_is_valid(self):
        """Test that Moon altitude is reasonable for visible eclipse."""
        jd_start = julday(2024, 1, 1, 0)
        paris_lat, paris_lon = 48.8566, 2.3522

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, paris_lat, paris_lon)

        moon_alt = attr[4]

        # For a visible eclipse, Moon should be above horizon at some point
        # The altitude at maximum might be below horizon if eclipse is partially visible
        # during moonrise/moonset, but it should still be reasonable
        assert -90 <= moon_alt <= 90

    def test_swe_alias(self):
        """Test that swe_lun_eclipse_when_loc is an alias for lun_eclipse_when_loc."""
        jd_start = julday(2024, 1, 1, 0)
        cape_town_lat, cape_town_lon = -33.9249, 18.4241

        ecl_type1, times1, attr1 = lun_eclipse_when_loc(
            jd_start, cape_town_lat, cape_town_lon
        )
        ecl_type2, times2, attr2 = swe_lun_eclipse_when_loc(
            jd_start, cape_town_lat, cape_town_lon
        )

        assert times1 == times2
        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_sequential_visible_eclipses(self):
        """Test finding multiple sequential visible lunar eclipses."""
        jd = julday(2024, 1, 1, 0)
        mumbai_lat, mumbai_lon = 19.0760, 72.8777
        eclipses = []

        # Find 3 sequential visible eclipses
        for _ in range(3):
            ecl_type, times, attr = lun_eclipse_when_loc(jd, mumbai_lat, mumbai_lon)
            eclipses.append((times[0], ecl_type))
            jd = times[0] + 1  # Start after this eclipse

        # Each eclipse should be later than the previous
        for i in range(1, len(eclipses)):
            assert eclipses[i][0] > eclipses[i - 1][0]

    def test_known_eclipse_may_2022_visible_from_americas(self):
        """Test the total lunar eclipse of May 16, 2022 visible from Americas.

        Reference: NASA - Total Lunar Eclipse of May 16, 2022
        Maximum eclipse: approximately 04:12 UTC
        Visible from North and South America, Europe, Africa
        """
        # Start searching from May 1, 2022
        jd_start = julday(2022, 5, 1, 0)
        # Rio de Janeiro - should see this eclipse
        rio_lat, rio_lon = -22.9068, -43.1729

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, rio_lat, rio_lon)

        # Should be a total eclipse
        assert ecl_type & SE_ECL_TOTAL
        assert ecl_type & SE_ECL_VISIBLE

        # Check maximum is on May 16, 2022
        year, month, day, hour = revjul(times[0])
        assert year == 2022
        assert month == 5
        assert day == 16

    def test_altitude_parameter(self):
        """Test that altitude parameter is accepted and doesn't cause errors."""
        jd_start = julday(2024, 1, 1, 0)
        denver_lat, denver_lon = 39.7392, -104.9903
        denver_altitude = 1609  # 1 mile high

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, denver_lat, denver_lon, altitude=denver_altitude
        )

        # Should find an eclipse
        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_southern_hemisphere(self):
        """Test finding eclipses from southern hemisphere."""
        jd_start = julday(2024, 1, 1, 0)
        # Melbourne, Australia
        melbourne_lat, melbourne_lon = -37.8136, 144.9631

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, melbourne_lat, melbourne_lon
        )

        # Should find an eclipse
        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_high_latitude_location(self):
        """Test finding eclipses from high latitude location."""
        jd_start = julday(2024, 1, 1, 0)
        # Reykjavik, Iceland
        reykjavik_lat, reykjavik_lon = 64.1466, -21.9426

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, reykjavik_lat, reykjavik_lon
        )

        # Should find an eclipse (might take longer due to polar location)
        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start


class TestLunEclipseWhenLocEdgeCases:
    """Test edge cases for location-based lunar eclipse calculations."""

    def test_equatorial_location(self):
        """Test finding eclipses from equatorial location."""
        jd_start = julday(2024, 1, 1, 0)
        # Quito, Ecuador (nearly on the equator)
        quito_lat, quito_lon = -0.1807, -78.4678

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, quito_lat, quito_lon)

        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_prime_meridian(self):
        """Test finding eclipses from location on prime meridian."""
        jd_start = julday(2024, 1, 1, 0)
        # Greenwich, UK
        greenwich_lat, greenwich_lon = 51.4772, 0.0

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, greenwich_lat, greenwich_lon
        )

        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_international_date_line(self):
        """Test finding eclipses from location near date line."""
        jd_start = julday(2024, 1, 1, 0)
        # Fiji
        fiji_lat, fiji_lon = -18.1416, 178.4419

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, fiji_lat, fiji_lon)

        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_early_20th_century(self):
        """Test finding eclipse in early 20th century."""
        # 1920 CE (within ephemeris range)
        jd_start = julday(1920, 1, 1, 0)
        london_lat, london_lon = 51.5074, -0.1278

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, london_lat, london_lon)

        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_mid_21st_century(self):
        """Test finding eclipse in mid 21st century."""
        # 2050 CE (within ephemeris range)
        jd_start = julday(2050, 1, 1, 0)
        new_york_lat, new_york_lon = 40.7128, -74.0060

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, new_york_lat, new_york_lon
        )

        assert ecl_type & SE_ECL_VISIBLE
        assert times[0] > jd_start

    def test_moonrise_moonset_during_eclipse(self):
        """Test that moonrise/moonset times are calculated when applicable."""
        # This test verifies the function doesn't crash when calculating
        # moonrise/moonset times. The actual times depend on specific eclipse
        # geometry and location.
        jd_start = julday(2024, 1, 1, 0)
        moscow_lat, moscow_lon = 55.7558, 37.6173

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, moscow_lat, moscow_lon)

        # times[7] is moonrise, times[8] is moonset
        # These may be 0 if Moon doesn't rise/set during eclipse
        assert isinstance(times[7], float)
        assert isinstance(times[8], float)
        # If set, they should be within the eclipse window
        if times[7] > 0:
            assert times[5] <= times[7] <= times[6]  # Within penumbral phase
        if times[8] > 0:
            assert times[5] <= times[8] <= times[6]  # Within penumbral phase
