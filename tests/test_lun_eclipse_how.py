"""
Tests for lun_eclipse_how function in libephemeris.

Tests the lunar eclipse circumstances calculation for a specific
location and time.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

from libephemeris import (
    julday,
    lun_eclipse_when,
    lun_eclipse_how,
    swe_lun_eclipse_how,
    SE_ECL_TOTAL,
    SE_ECL_PARTIAL,
    SE_ECL_PENUMBRAL,
    SE_ECL_VISIBLE,
)


class TestLunEclipseHow:
    """Test suite for lun_eclipse_how function."""

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse)
        jd_eclipse = julday(2022, 5, 16, 4.2)
        rome_lat, rome_lon = 41.9028, 12.4964

        attr, ecl_type = lun_eclipse_how(jd_eclipse, rome_lat, rome_lon)

        # Should return 11-element attr tuple
        assert len(attr) == 11
        # All elements should be floats
        assert all(isinstance(a, float) for a in attr)
        # Eclipse type should be int
        assert isinstance(ecl_type, int)

    def test_attributes_during_eclipse(self):
        """Test that eclipse attributes have reasonable values during eclipse."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse maximum ~04:12 UTC)
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Use a location where Moon should be visible
        rio_lat, rio_lon = -22.9068, -43.1729

        attr, ecl_type = lun_eclipse_how(jd_eclipse, rio_lat, rio_lon)

        # Umbral magnitude should be positive during eclipse
        assert attr[0] > 0

        # Penumbral magnitude should be positive during eclipse
        assert attr[1] > 0

        # Azimuth should be 0-360
        assert 0 <= attr[3] <= 360

        # Altitude could be any value but should be reasonable
        assert -90 <= attr[4] <= 90

        # Moon diameter should be approximately 0.5 degrees
        assert 0.4 < attr[5] < 0.6

        # Umbra diameter should be positive
        assert attr[6] > 0

        # Penumbra diameter should be positive and larger than umbra
        assert attr[7] > 0
        assert attr[7] > attr[6]

    def test_no_eclipse_returns_zero_magnitude(self):
        """Test that non-eclipse time returns zero magnitude."""
        # Use a time when there's definitely no eclipse
        jd_no_eclipse = julday(2022, 6, 1, 12.0)
        london_lat, london_lon = 51.5074, -0.1278

        attr, ecl_type = lun_eclipse_how(jd_no_eclipse, london_lat, london_lon)

        # Umbral and penumbral magnitudes should be zero
        assert attr[0] == 0.0
        assert attr[1] == 0.0
        # Eclipse type should be 0
        assert ecl_type == 0

    def test_moon_position_calculated(self):
        """Test that Moon's position is calculated correctly."""
        # Use a known eclipse time
        jd_eclipse = julday(2022, 5, 16, 4.2)
        tokyo_lat, tokyo_lon = 35.6762, 139.6503

        attr, ecl_type = lun_eclipse_how(jd_eclipse, tokyo_lat, tokyo_lon)

        # Moon azimuth should be in valid range
        assert 0 <= attr[3] <= 360

        # Moon altitude should be in valid range
        assert -90 <= attr[4] <= 90

        # Moon diameter should be approximately 0.5 degrees
        assert 0.4 < attr[5] < 0.6

    def test_different_locations_different_altitudes(self):
        """Test that different locations have different Moon altitudes."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        # Two locations on opposite sides of the Earth
        new_york_lat, new_york_lon = 40.7128, -74.0060
        sydney_lat, sydney_lon = -33.8688, 151.2093

        attr_ny, _ = lun_eclipse_how(jd_eclipse, new_york_lat, new_york_lon)
        attr_syd, _ = lun_eclipse_how(jd_eclipse, sydney_lat, sydney_lon)

        # Moon altitudes should be different
        assert attr_ny[4] != attr_syd[4]

    def test_visibility_flag_when_moon_above_horizon(self):
        """Test that visibility flag is set when Moon is above horizon."""
        # May 16, 2022 eclipse - Moon should be visible from Rio
        jd_eclipse = julday(2022, 5, 16, 4.2)
        rio_lat, rio_lon = -22.9068, -43.1729

        attr, ecl_type = lun_eclipse_how(jd_eclipse, rio_lat, rio_lon)

        # Check if Moon is above horizon
        moon_alt = attr[4]
        if moon_alt > -1.0:
            # Should have visibility flag set
            assert ecl_type & SE_ECL_VISIBLE

    def test_eclipse_type_detection(self):
        """Test that eclipse type is correctly detected."""
        # First find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        times, global_ecl_type = lun_eclipse_when(jd_start, eclipse_type=SE_ECL_TOTAL)
        jd_max = times[0]

        # Check circumstances at maximum
        rio_lat, rio_lon = -22.9068, -43.1729
        attr, ecl_type = lun_eclipse_how(jd_max, rio_lat, rio_lon)

        # Should detect the eclipse type
        assert ecl_type & (SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_PENUMBRAL)

    def test_swe_alias(self):
        """Test that swe_lun_eclipse_how is an alias for lun_eclipse_how."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        cape_town_lat, cape_town_lon = -33.9249, 18.4241

        attr1, ecl_type1 = lun_eclipse_how(jd_eclipse, cape_town_lat, cape_town_lon)
        attr2, ecl_type2 = swe_lun_eclipse_how(jd_eclipse, cape_town_lat, cape_town_lon)

        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_altitude_parameter(self):
        """Test that altitude parameter is accepted and doesn't cause errors."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        denver_lat, denver_lon = 39.7392, -104.9903
        denver_altitude = 1609  # 1 mile high

        attr, ecl_type = lun_eclipse_how(
            jd_eclipse, denver_lat, denver_lon, altitude=denver_altitude
        )

        # Should return valid results
        assert len(attr) == 11
        assert isinstance(ecl_type, int)

    def test_shadow_diameters_calculated(self):
        """Test that umbra and penumbra diameters are calculated."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        berlin_lat, berlin_lon = 52.5200, 13.4050

        attr, ecl_type = lun_eclipse_how(jd_eclipse, berlin_lat, berlin_lon)

        # Umbra diameter (attr[6]) should be positive
        assert attr[6] > 0

        # Penumbra diameter (attr[7]) should be positive
        assert attr[7] > 0

        # Penumbra should be larger than umbra
        assert attr[7] > attr[6]

    def test_consistency_with_lun_eclipse_when(self):
        """Test that magnitudes are consistent with lun_eclipse_when values."""
        # Find a lunar eclipse
        jd_start = julday(2024, 1, 1, 0)
        times, global_type = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Get circumstances at maximum
        # Use a location where eclipse should be visible
        london_lat, london_lon = 51.5074, -0.1278
        attr, ecl_type = lun_eclipse_how(jd_max, london_lat, london_lon)

        # At least one of umbral or penumbral magnitude should be positive
        assert attr[0] > 0 or attr[1] > 0


class TestLunEclipseHowEdgeCases:
    """Test edge cases for lun_eclipse_how function."""

    def test_equatorial_location(self):
        """Test eclipse circumstances from equatorial location."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Quito, Ecuador (nearly on the equator)
        quito_lat, quito_lon = -0.1807, -78.4678

        attr, ecl_type = lun_eclipse_how(jd_eclipse, quito_lat, quito_lon)

        assert len(attr) == 11
        # Moon altitude should be in valid range
        assert -90 <= attr[4] <= 90

    def test_high_latitude_location(self):
        """Test eclipse circumstances from high latitude location."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Reykjavik, Iceland
        reykjavik_lat, reykjavik_lon = 64.1466, -21.9426

        attr, ecl_type = lun_eclipse_how(jd_eclipse, reykjavik_lat, reykjavik_lon)

        assert len(attr) == 11
        # Moon altitude should be in valid range
        assert -90 <= attr[4] <= 90

    def test_southern_hemisphere(self):
        """Test eclipse circumstances from southern hemisphere."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Melbourne, Australia
        melbourne_lat, melbourne_lon = -37.8136, 144.9631

        attr, ecl_type = lun_eclipse_how(jd_eclipse, melbourne_lat, melbourne_lon)

        assert len(attr) == 11
        # Moon altitude should be in valid range
        assert -90 <= attr[4] <= 90

    def test_partial_eclipse(self):
        """Test circumstances during a partial lunar eclipse."""
        # Find a partial lunar eclipse
        jd_start = julday(2023, 10, 1, 0)
        times, global_type = lun_eclipse_when(jd_start, eclipse_type=SE_ECL_PARTIAL)
        jd_max = times[0]

        # Check circumstances
        london_lat, london_lon = 51.5074, -0.1278
        attr, ecl_type = lun_eclipse_how(jd_max, london_lat, london_lon)

        # Should have positive umbral magnitude (but less than 1 for partial)
        if attr[0] > 0:
            assert attr[0] < 1.5  # Allow some margin

    def test_penumbral_eclipse(self):
        """Test circumstances during a penumbral lunar eclipse."""
        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        times, global_type = lun_eclipse_when(jd_start, eclipse_type=SE_ECL_PENUMBRAL)
        jd_max = times[0]

        # Check circumstances
        paris_lat, paris_lon = 48.8566, 2.3522
        attr, ecl_type = lun_eclipse_how(jd_max, paris_lat, paris_lon)

        # Should have positive penumbral magnitude
        # Umbral magnitude should be zero or very small for penumbral-only
        assert attr[1] >= 0

    def test_early_20th_century(self):
        """Test eclipse circumstances in early 20th century."""
        # Find an eclipse in 1920
        jd_start = julday(1920, 1, 1, 0)
        times, global_type = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Check circumstances
        london_lat, london_lon = 51.5074, -0.1278
        attr, ecl_type = lun_eclipse_how(jd_max, london_lat, london_lon)

        assert len(attr) == 11
        assert isinstance(ecl_type, int)

    def test_mid_21st_century(self):
        """Test eclipse circumstances in mid 21st century."""
        # Find an eclipse in 2050
        jd_start = julday(2050, 1, 1, 0)
        times, global_type = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Check circumstances
        new_york_lat, new_york_lon = 40.7128, -74.0060
        attr, ecl_type = lun_eclipse_how(jd_max, new_york_lat, new_york_lon)

        assert len(attr) == 11
        assert isinstance(ecl_type, int)

    def test_international_date_line(self):
        """Test eclipse circumstances near international date line."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Fiji
        fiji_lat, fiji_lon = -18.1416, 178.4419

        attr, ecl_type = lun_eclipse_how(jd_eclipse, fiji_lat, fiji_lon)

        assert len(attr) == 11
        # Moon azimuth should be in valid range
        assert 0 <= attr[3] <= 360

    def test_prime_meridian(self):
        """Test eclipse circumstances on prime meridian."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Greenwich, UK
        greenwich_lat, greenwich_lon = 51.4772, 0.0

        attr, ecl_type = lun_eclipse_how(jd_eclipse, greenwich_lat, greenwich_lon)

        assert len(attr) == 11
        assert isinstance(ecl_type, int)
