"""
Tests for lunar occultation calculations in libephemeris.

Tests the lun_occult_when_glob function which finds lunar occultations
of planets and fixed stars.

Lunar occultations occur when the Moon passes in front of a planet or star.
These events are relatively rare and happen in specific cycles depending on
the Moon's nodal position relative to the target body's ecliptic latitude.

Known occultation series:
- Regulus: 2017 series (June 28, 2017 was an occultation)
"""

import pytest
from libephemeris import (
    julday,
    revjul,
    lun_occult_when_glob,
    swe_lun_occult_when_glob,
    SE_ECL_TOTAL,
    SE_ECL_PARTIAL,
)


class TestLunOccultWhenGlob:
    """Test suite for lun_occult_when_glob function."""

    def test_finds_star_occultation_regulus(self):
        """Test that function finds a lunar occultation of Regulus.

        Regulus (Alpha Leonis) is frequently occulted by the Moon
        because it lies very close to the ecliptic.
        Known occultation: June 28, 2017.
        """
        # Start from Jan 1, 2017 - there's a known occultation on June 28, 2017
        jd_start = julday(2017, 1, 1, 0)

        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")

        # Should find an occultation
        assert ocl_type != 0
        assert times[0] > jd_start  # Maximum should be after start
        assert times[0] > 0  # Should have valid maximum time

        # Verify it's around June 2017
        year, month, day, hour = revjul(times[0])
        assert year == 2017
        assert month == 6

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2017, 1, 1, 0)

        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")

        # Should return 8-element tuple
        assert len(times) == 8
        # All elements should be floats
        assert all(isinstance(t, float) for t in times)
        # Occultation type should be int
        assert isinstance(ocl_type, int)

    def test_contact_times_ordering(self):
        """Test that contact times are in correct chronological order."""
        jd_start = julday(2017, 1, 1, 0)

        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")

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
        """Test that occultation type flags are set."""
        jd_start = julday(2017, 1, 1, 0)

        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")

        # Should be either total or partial
        assert (ocl_type & SE_ECL_TOTAL) or (ocl_type & SE_ECL_PARTIAL)

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_when_glob(jd_start, 0, "")

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd_start = julday(2017, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_when_glob(jd_start, 0, "UnknownStar123")

    def test_swe_alias(self):
        """Test that swe_lun_occult_when_glob is an alias."""
        jd_start = julday(2017, 1, 1, 0)

        times1, ocl_type1 = lun_occult_when_glob(jd_start, 0, "Regulus")
        times2, ocl_type2 = swe_lun_occult_when_glob(jd_start, 0, "Regulus")

        assert times1 == times2
        assert ocl_type1 == ocl_type2

    def test_occultation_duration_reasonable(self):
        """Test that occultation duration is reasonable.

        Lunar occultations typically last a few seconds to about an hour,
        depending on the geometry.
        """
        jd_start = julday(2017, 1, 1, 0)

        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")

        jd_first = times[1]
        jd_fourth = times[4]

        if jd_first > 0 and jd_fourth > 0:
            duration_hours = (jd_fourth - jd_first) * 24
            # Duration should be between 0.01 hours (~30 seconds) and 2 hours
            assert 0.001 < duration_hours < 2.0


class TestLunOccultEdgeCases:
    """Test edge cases for lunar occultation calculations."""

    def test_known_regulus_occultation_2017(self):
        """Test finding the known Regulus occultation of June 28, 2017.

        This is a well-documented lunar occultation event.
        """
        # Search from June 1, 2017
        jd_start = julday(2017, 6, 1, 0)

        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")

        # Verify it's on June 28
        year, month, day, hour = revjul(times[0])
        assert year == 2017
        assert month == 6
        assert 27 <= day <= 28  # Allow for timezone differences

    def test_search_continues_after_no_occultation(self):
        """Test that search properly advances when no immediate occultation."""
        # Start from a date with no immediate occultation
        jd_start = julday(2020, 1, 1, 0)

        # The function should either find an occultation or raise RuntimeError
        # after exhausting the search limit
        try:
            times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")
            # If found, verify it's after start date
            assert times[0] > jd_start
        except RuntimeError as e:
            # Expected if no occultation in 20 years
            assert "No lunar occultation" in str(e)

    def test_raises_error_for_invalid_planet(self):
        """Test that invalid planet ID raises appropriate error."""
        jd_start = julday(2020, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_when_glob(jd_start, 999, "")  # Invalid planet ID


class TestLunOccultResolveStarId:
    """Regression tests for _resolve_star_id tuple unpacking fix.

    These tests verify that _resolve_star_id returning 3 values
    (star_id, error_message, canonical_name) is correctly handled
    throughout the eclipse module.
    """

    def test_valid_star_does_not_raise_unpack_error(self):
        """Test that valid star names don't cause unpacking errors.

        Regression test for ValueError: too many values to unpack (expected 2).
        The _resolve_star_id function returns 3 values, not 2.
        """
        jd_start = julday(2017, 1, 1, 0)

        # This should not raise "ValueError: too many values to unpack"
        times, ocl_type = lun_occult_when_glob(jd_start, 0, "Regulus")
        assert times[0] > jd_start

    def test_invalid_star_raises_value_error_with_message(self):
        """Test that invalid stars raise ValueError with descriptive message."""
        jd_start = julday(2017, 1, 1, 0)

        with pytest.raises(ValueError) as exc_info:
            lun_occult_when_glob(jd_start, 0, "InvalidStarXYZ123")

        # Error message should mention the star name
        assert "invalidstarxyz123" in str(exc_info.value).lower()
