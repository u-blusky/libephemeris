"""
Tests for lunar occultation location calculations in libephemeris.

Tests the lun_occult_where function which calculates where on Earth
a lunar occultation is visible at a given time.

Lunar occultations occur when the Moon passes in front of a planet or star.
This function determines the geographic coordinates where the occultation
is visible at a specific moment.
"""

import pytest
from libephemeris import (
    julday,
    lun_occult_where,
    swe_lun_occult_where,
    lun_occult_when_glob,
    SE_ECL_TOTAL,
    SE_ECL_PARTIAL,
)


class TestLunOccultWhere:
    """Test suite for lun_occult_where function."""

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        # First find a known occultation time
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Now get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        # geopos should be 10-element tuple
        assert len(geopos) == 10
        assert all(isinstance(g, float) for g in geopos)

        # attr should be 20-element tuple (pyswisseph compatible)
        assert len(attr) == 20
        assert all(isinstance(a, float) for a in attr)

        # Occultation type should be int
        assert isinstance(ocl_type, int)

    def test_finds_valid_location_during_occultation(self):
        """Test that function finds valid location during a known occultation."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        # Should find an occultation
        assert ocl_type != 0

        # Longitude should be valid (-180 to 180)
        assert -180.0 <= geopos[0] <= 180.0

        # Latitude should be valid (-90 to 90)
        assert -90.0 <= geopos[1] <= 90.0

    def test_no_occultation_when_moon_far_from_star(self):
        """Test that function returns 0 when no occultation is happening."""
        # Use a random time when Moon is not near Regulus
        jd = julday(2024, 1, 1, 12)  # Random date

        ocl_type, geopos, attr = lun_occult_where(jd, 0, "Regulus")

        # Should return 0 for no occultation
        assert ocl_type == 0

        # geopos should be zeros
        assert all(g == 0.0 for g in geopos)

    def test_occultation_type_flags(self):
        """Test that occultation type flags are set correctly."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        # Should be either total or partial
        assert (ocl_type & SE_ECL_TOTAL) or (ocl_type & SE_ECL_PARTIAL)

    def test_geographic_limits_reasonable(self):
        """Test that geographic limits are reasonable."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        if ocl_type != 0:
            # Central line
            central_lon, central_lat = geopos[0], geopos[1]
            # Northern limit
            north_lon, north_lat = geopos[2], geopos[3]
            # Southern limit
            south_lon, south_lat = geopos[4], geopos[5]

            # All latitudes should be valid
            assert -90.0 <= central_lat <= 90.0
            assert -90.0 <= north_lat <= 90.0
            assert -90.0 <= south_lat <= 90.0

            # All longitudes should be valid
            assert -180.0 <= central_lon <= 180.0
            assert -180.0 <= north_lon <= 180.0
            assert -180.0 <= south_lon <= 180.0

            # Northern limit should be north of southern limit
            assert north_lat >= south_lat

    def test_attributes_reasonable(self):
        """Test that occultation attributes are reasonable."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        if ocl_type != 0:
            # Fraction covered should be between 0 and 1
            assert 0.0 <= attr[0] <= 1.0

            # Diameter ratio should be positive (and small for stars)
            assert attr[1] >= 0.0

            # Obscuration (attr[2]) should be between 0 and 1
            assert 0.0 <= attr[2] <= 1.0

            # Path width should be reasonable (0 to 1000 km)
            assert 0.0 <= attr[3] <= 1000.0

            # Moon azimuth should be 0-360
            assert 0.0 <= attr[4] < 360.0 or -180.0 <= attr[4] <= 180.0

            # Moon altitude should be -90 to 90
            assert -90.0 <= attr[5] <= 90.0

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_where(jd, 0, "")

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd = julday(2017, 6, 28, 10)

        with pytest.raises(ValueError):
            lun_occult_where(jd, 0, "UnknownStar123")

    def test_swe_alias(self):
        """Test that swe_lun_occult_where is an alias."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        ocl_type1, geopos1, attr1 = lun_occult_where(jd_max, 0, "Regulus")
        ocl_type2, geopos2, attr2 = swe_lun_occult_where(jd_max, 0, "Regulus")

        assert geopos1 == geopos2
        assert attr1 == attr2
        assert ocl_type1 == ocl_type2

    def test_central_latitude_near_moon_declination(self):
        """Test that central latitude is near Moon's declination.

        The sub-lunar point should have latitude approximately equal
        to the Moon's declination at the time of occultation.
        """
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        if ocl_type != 0:
            central_lat = geopos[1]
            # Moon's declination ranges from about -28 to +28 degrees
            # Central latitude should be within this range
            assert -30.0 <= central_lat <= 30.0


class TestLunOccultWhereEdgeCases:
    """Test edge cases for lun_occult_where function."""

    def test_occultation_at_different_times_same_event(self):
        """Test occultation location changes during an event."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get location at maximum
        ocl_type_max, geopos_max, _ = lun_occult_where(jd_max, 0, "Regulus")

        if ocl_type_max != 0:
            # The longitude should change as the Earth rotates
            # At maximum, we should find a valid location
            assert geopos_max[0] != 0.0 or geopos_max[1] != 0.0

    def test_star_occultation_total(self):
        """Test that star occultations are typically total.

        Since stars have negligible angular size compared to the Moon,
        they should produce total occultations when they happen.
        """
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        if ocl_type != 0:
            # For a star, it should be total when visible
            # (since star angular size is negligible)
            assert ocl_type & SE_ECL_TOTAL

    def test_fraction_covered_during_total_occultation(self):
        """Test that fraction covered is 1.0 for total occultation."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        if ocl_type & SE_ECL_TOTAL:
            # Fraction covered should be 1.0 for total occultation
            assert attr[0] >= 0.99  # Allow small numerical error


class TestLunOccultWhereIntegration:
    """Integration tests for lun_occult_where with other functions."""

    def test_consistency_with_lun_occult_when_glob(self):
        """Test that lun_occult_where is consistent with lun_occult_when_glob."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        glob_type, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Get where it's visible
        ocl_type, geopos, attr = lun_occult_where(jd_max, 0, "Regulus")

        # Both should indicate an occultation
        assert glob_type != 0
        assert ocl_type != 0

        # The type should be consistent (both total or both partial)
        if glob_type & SE_ECL_TOTAL:
            assert ocl_type & SE_ECL_TOTAL
        if glob_type & SE_ECL_PARTIAL:
            assert ocl_type & SE_ECL_PARTIAL

    def test_multiple_calls_same_result(self):
        """Test that calling the function multiple times gives same result."""
        # Find a Regulus occultation
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(jd_start, "Regulus")
        jd_max = times[0]

        # Call multiple times
        result1 = lun_occult_where(jd_max, 0, "Regulus")
        result2 = lun_occult_where(jd_max, 0, "Regulus")
        result3 = lun_occult_where(jd_max, 0, "Regulus")

        # All results should be identical
        assert result1 == result2
        assert result2 == result3
