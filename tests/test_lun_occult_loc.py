"""
Tests for lunar occultation location-specific calculations in libephemeris.

Tests the lun_occult_when_loc function which finds lunar occultations
of planets and fixed stars visible from a specific geographic location.

Lunar occultations occur when the Moon passes in front of a planet or star.
This function adds location-specific visibility checking to ensure the
occultation is observable from the given coordinates.
"""

import pytest
from libephemeris import (
    julday,
    lun_occult_when_loc,
    swe_lun_occult_when_loc,
    SE_ECL_TOTAL,
    SE_ECL_PARTIAL,
    SE_ECL_VISIBLE,
    SE_ECL_MAX_VISIBLE,
    SE_ECL_1ST_VISIBLE,
    SE_ECL_4TH_VISIBLE,
)


class TestLunOccultWhenLoc:
    """Test suite for lun_occult_when_loc function."""

    def test_finds_star_occultation_regulus_from_location(self):
        """Test finding a lunar occultation of Regulus visible from a location.

        Regulus (Alpha Leonis) is frequently occulted by the Moon
        because it lies very close to the ecliptic.
        """
        # Start from Jan 1, 2017 - known occultation series in 2017
        jd_start = julday(2017, 1, 1, 0)

        # Search from a location in Europe (central location for good visibility)
        lat, lon = 45.0, 10.0  # Northern Italy

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # Should find an occultation
        assert ocl_type != 0
        assert times[0] > jd_start  # Maximum should be after start
        assert times[0] > 0  # Should have valid maximum time

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # Times should be 10-element tuple
        assert len(times) == 10
        assert all(isinstance(t, float) for t in times)

        # Attributes should be 20-element tuple (pyswisseph compatible)
        assert len(attr) == 20
        assert all(isinstance(a, float) for a in attr)

        # Occultation type should be int
        assert isinstance(ocl_type, int)

    def test_contact_times_ordering(self):
        """Test that contact times are in correct chronological order."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        jd_max = times[0]
        jd_first = times[1]
        jd_second = times[2]
        jd_third = times[3]
        jd_fourth = times[4]

        # First contact should be before maximum
        if jd_first > 0:
            assert jd_first < jd_max

        # Fourth contact should be after maximum
        if jd_fourth > 0:
            assert jd_fourth > jd_max

        # If we have second/third contacts (total occultation)
        if jd_second > 0 and jd_third > 0:
            assert jd_first <= jd_second <= jd_max
            assert jd_max <= jd_third <= jd_fourth

    def test_occultation_type_flags(self):
        """Test that occultation type flags are set correctly."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # Should be either total or partial
        assert (ocl_type & SE_ECL_TOTAL) or (ocl_type & SE_ECL_PARTIAL)

        # Should have visibility flags
        assert ocl_type & SE_ECL_VISIBLE

    def test_visibility_flags_set(self):
        """Test that visibility flags are properly set."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # At least one visibility flag should be set
        visibility_flags = (
            SE_ECL_VISIBLE
            | SE_ECL_MAX_VISIBLE
            | SE_ECL_1ST_VISIBLE
            | SE_ECL_4TH_VISIBLE
        )
        assert ocl_type & visibility_flags

    def test_moon_attributes_returned(self):
        """Test that occultation attributes are returned.

        Uses pyswisseph compatible format.
        """
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # Fraction covered should be between 0 and 1
        assert 0.0 <= attr[0] <= 1.0

        # Target azimuth should be between 0 and 360
        assert 0.0 <= attr[4] < 360.0

        # Target altitude should be between -90 and 90
        assert -90.0 <= attr[5] <= 90.0

        # Angular separation (elongation) should be small for an occultation
        assert 0.0 <= attr[7] < 1.0

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 45.0, 10.0

        with pytest.raises(ValueError):
            lun_occult_when_loc(jd_start, 0, "", lat, lon)

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        with pytest.raises(ValueError):
            lun_occult_when_loc(jd_start, 0, "UnknownStar123", lat, lon)

    def test_swe_alias(self):
        """Test that swe_lun_occult_when_loc is an alias."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times1, attr1, ocl_type1 = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)
        times2, attr2, ocl_type2 = swe_lun_occult_when_loc(
            jd_start, 0, "Regulus", lat, lon
        )

        assert times1 == times2
        assert attr1 == attr2
        assert ocl_type1 == ocl_type2

    def test_occultation_duration_reasonable(self):
        """Test that occultation duration is reasonable.

        Lunar occultations typically last a few seconds to about an hour,
        depending on the geometry.
        """
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        jd_first = times[1]
        jd_fourth = times[4]

        if jd_first > 0 and jd_fourth > 0:
            duration_hours = (jd_fourth - jd_first) * 24
            # Duration should be between 0.01 hours (~30 seconds) and 2 hours
            assert 0.001 < duration_hours < 2.0


class TestLunOccultLocDifferentLocations:
    """Test occultation visibility from different locations."""

    def test_same_occultation_different_visibility(self):
        """Test that the same occultation may be found but with different attributes."""
        jd_start = julday(2017, 1, 1, 0)

        # Two different locations
        lat1, lon1 = 51.5, -0.1  # London
        lat2, lon2 = 35.7, 139.7  # Tokyo

        # Both should find occultations
        times1, attr1, ocl_type1 = lun_occult_when_loc(
            jd_start, 0, "Regulus", lat1, lon1
        )
        times2, attr2, ocl_type2 = lun_occult_when_loc(
            jd_start, 0, "Regulus", lat2, lon2
        )

        # Both should find occultations
        assert times1[0] > 0
        assert times2[0] > 0

        # Moon altitude will be different at different locations
        # (even for the same event, if visible from both)
        # Just verify they're valid altitudes
        assert -90.0 <= attr1[4] <= 90.0
        assert -90.0 <= attr2[4] <= 90.0


class TestLunOccultLocEdgeCases:
    """Test edge cases for lunar occultation location calculations."""

    def test_altitude_parameter(self):
        """Test that altitude parameter is accepted."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0
        altitude = 1000.0  # 1000m above sea level

        # Should not raise an error
        times, attr, ocl_type = lun_occult_when_loc(
            jd_start, 0, "Regulus", lat, lon, altitude
        )

        assert times[0] > 0

    def test_polar_region(self):
        """Test occultation search from polar region."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 70.0, 25.0  # Northern Norway

        # Should find an occultation or raise RuntimeError if none visible
        try:
            times, attr, ocl_type = lun_occult_when_loc(
                jd_start, 0, "Regulus", lat, lon
            )
            assert times[0] > jd_start
        except RuntimeError as e:
            # It's valid if no occultation is visible from this location
            assert "No lunar occultation" in str(e)

    def test_southern_hemisphere(self):
        """Test occultation search from southern hemisphere."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = -33.9, 18.4  # Cape Town, South Africa

        # Should find an occultation or raise RuntimeError if none visible
        try:
            times, attr, ocl_type = lun_occult_when_loc(
                jd_start, 0, "Regulus", lat, lon
            )
            assert times[0] > jd_start
        except RuntimeError as e:
            assert "No lunar occultation" in str(e)

    def test_moonrise_moonset_times(self):
        """Test that moonrise/moonset times are provided when applicable."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # Moonrise time is at index 7
        moonrise = times[7]
        # Moonset time is at index 8
        moonset = times[8]

        # These should be 0 or valid JD times
        assert moonrise >= 0
        assert moonset >= 0

        # If moonrise is set, it should be within the occultation window
        if moonrise > 0 and times[1] > 0 and times[4] > 0:
            assert times[1] <= moonrise <= times[4]

        # If moonset is set, it should be within the occultation window
        if moonset > 0 and times[1] > 0 and times[4] > 0:
            assert times[1] <= moonset <= times[4]

    def test_angular_separation_at_maximum(self):
        """Test that angular separation at maximum is provided."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        times, attr, ocl_type = lun_occult_when_loc(jd_start, 0, "Regulus", lat, lon)

        # Angular separation at maximum is at index 7
        min_separation = attr[7]

        # Should be a reasonable value during an occultation
        # Moon's angular radius is about 0.25 degrees, but topocentric
        # effects can cause some variation in the separation value
        assert 0 <= min_separation < 1.0
