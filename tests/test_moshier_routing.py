"""
Tests for SEFLG_MOSEPH routing gate in swe_calc_ut() and swe_calc().

This module tests the routing behavior when SEFLG_MOSEPH flag is passed:
- SEFLG_MOSEPH should route to the Moshier ephemeris system
- Moshier range validation should work independently from JPL range validation
- Moshier provides extended date range (-3000 to +3000 CE) via semi-analytical algorithms
"""

import pytest
import libephemeris as eph
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SEFLG_MOSEPH,
    SEFLG_SWIEPH,
    SEFLG_JPLEPH,
    SEFLG_SPEED,
)
from libephemeris.exceptions import (
    EphemerisRangeError,
    validate_jd_range_moshier,
    MOSHIER_JD_START,
    MOSHIER_JD_END,
)


class TestMoshierRangeConstants:
    """Test the Moshier ephemeris range constants."""

    def test_moshier_range_constants_exist(self):
        """Moshier range constants should be defined."""
        assert MOSHIER_JD_START == 625673.5  # -3000 Jan 1
        assert MOSHIER_JD_END == 3182395.5  # +3000 Dec 31

    def test_moshier_range_wider_than_jpl(self):
        """Moshier range should be much wider than JPL DE440 range."""
        # DE440 covers ~1550-2650, Moshier covers -3000 to +3000
        # This ensures the test is meaningful
        jpl_start = 2287184.5  # ~1550
        jpl_end = 2688976.5  # ~2650

        assert MOSHIER_JD_START < jpl_start
        assert MOSHIER_JD_END > jpl_end


class TestValidateJdRangeMoshier:
    """Test the validate_jd_range_moshier() function."""

    def test_valid_jd_within_moshier_range(self):
        """JD within Moshier range should not raise."""
        # J2000.0 - well within range
        validate_jd_range_moshier(2451545.0)  # Should not raise

    def test_valid_jd_at_range_edges(self):
        """JD at edges of Moshier range should not raise."""
        validate_jd_range_moshier(MOSHIER_JD_START)  # Should not raise
        validate_jd_range_moshier(MOSHIER_JD_END)  # Should not raise

    def test_jd_before_moshier_range_raises_error(self):
        """JD before Moshier start should raise EphemerisRangeError."""
        jd_before = MOSHIER_JD_START - 1.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range_moshier(jd_before)

        err = exc_info.value
        assert err.requested_jd == jd_before
        assert err.start_jd == MOSHIER_JD_START
        assert err.end_jd == MOSHIER_JD_END
        assert "Moshier" in err.message

    def test_jd_after_moshier_range_raises_error(self):
        """JD after Moshier end should raise EphemerisRangeError."""
        jd_after = MOSHIER_JD_END + 1.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range_moshier(jd_after)

        err = exc_info.value
        assert err.requested_jd == jd_after
        assert "Moshier" in err.message

    def test_error_includes_body_info_when_provided(self):
        """Error should include body information when body_id is provided."""
        jd_before = MOSHIER_JD_START - 1.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range_moshier(jd_before, body_id=SE_SUN)

        err = exc_info.value
        assert err.body_id == SE_SUN
        assert err.body_name == "Sun"

    def test_error_includes_func_name_when_provided(self):
        """Error message should include function name when provided."""
        jd_before = MOSHIER_JD_START - 1.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range_moshier(jd_before, func_name="test_func")

        err = exc_info.value
        assert "test_func:" in err.message


class TestSweCalcUtMoshierRouting:
    """Test SEFLG_MOSEPH routing in swe_calc_ut()."""

    def test_moseph_flag_routes_to_moshier(self):
        """swe_calc_ut with SEFLG_MOSEPH should route to Moshier and return valid results."""
        jd = 2451545.0  # J2000.0

        # Moshier is now implemented - should return valid position
        pos, flag = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

        # Check we got valid results
        assert len(pos) == 6  # lon, lat, dist, speed_lon, speed_lat, speed_dist
        assert 0 <= pos[0] < 360  # Valid longitude
        assert -90 <= pos[1] <= 90  # Valid latitude
        assert pos[2] > 0  # Positive distance

    def test_moseph_with_speed_flag(self):
        """SEFLG_MOSEPH combined with SEFLG_SPEED should return velocities."""
        jd = 2451545.0

        pos, flag = eph.swe_calc_ut(jd, SE_MARS, SEFLG_MOSEPH | SEFLG_SPEED)

        # Should return valid position with velocities
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Valid longitude
        # Mars typically moves ~0.3-0.8 deg/day in longitude
        assert abs(pos[3]) < 2.0  # Reasonable speed

    def test_moseph_range_validation_different_from_jpl(self):
        """SEFLG_MOSEPH should use Moshier range validation, not JPL range."""
        # This date is outside JPL DE440 range but inside Moshier range
        # JD 2000000.0 is ~764 CE - outside DE440 (1550-2650) but inside Moshier (-3000 to +3000)
        jd_outside_jpl = 2000000.0

        # With SEFLG_MOSEPH, should work (Moshier covers this date)
        pos, flag = eph.swe_calc_ut(jd_outside_jpl, SE_SUN, SEFLG_MOSEPH)
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Valid longitude

    def test_moseph_outside_moshier_range_raises_ephemeris_range_error(self):
        """SEFLG_MOSEPH with date outside Moshier range should raise EphemerisRangeError."""
        # JD 100000.0 is ~4660 BCE - outside Moshier range (-3000 to +3000)
        jd_outside_moshier = 100000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.swe_calc_ut(jd_outside_moshier, SE_SUN, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in err.message

    def test_without_moseph_flag_uses_jpl(self):
        """Without SEFLG_MOSEPH, calculation should use JPL ephemeris."""
        jd = 2451545.0  # J2000.0 - within JPL range

        # Should succeed with SEFLG_SWIEPH (uses JPL/Skyfield)
        pos, flag = eph.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH)
        assert len(pos) == 6  # lon, lat, dist, speed_lon, speed_lat, speed_dist

    def test_default_flags_use_jpl(self):
        """Default flags (0) should use JPL ephemeris, not Moshier."""
        jd = 2451545.0

        # Should succeed without SEFLG_MOSEPH
        pos, flag = eph.swe_calc_ut(jd, SE_SUN, 0)
        assert len(pos) == 6


class TestSweCalcMoshierRouting:
    """Test SEFLG_MOSEPH routing in swe_calc() (TT time)."""

    def test_moseph_flag_routes_to_moshier(self):
        """swe_calc with SEFLG_MOSEPH should route to Moshier and return valid results."""
        jd = 2451545.0  # J2000.0

        # Moshier is now implemented - should return valid position
        pos, flag = eph.swe_calc(jd, SE_MOON, SEFLG_MOSEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Valid longitude

    def test_moseph_range_validation_in_swe_calc(self):
        """swe_calc should also validate Moshier range when SEFLG_MOSEPH is set."""
        # Outside Moshier range
        jd_outside = 100000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.swe_calc(jd_outside, SE_JUPITER, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in err.message

    def test_without_moseph_uses_jpl(self):
        """swe_calc without SEFLG_MOSEPH should use JPL ephemeris."""
        jd = 2451545.0

        pos, flag = eph.swe_calc(jd, SE_MOON, SEFLG_JPLEPH)
        assert len(pos) == 6


class TestMoshierVsJplRanges:
    """Test the difference between Moshier and JPL range behavior."""

    def test_date_in_jpl_range_works_without_moseph(self):
        """Dates within JPL range should work without SEFLG_MOSEPH."""
        jd = 2451545.0  # J2000.0

        # Should work with default flags
        pos, _ = eph.swe_calc_ut(jd, SE_SUN, 0)
        assert abs(pos[0]) < 360  # Valid longitude

    def test_date_outside_jpl_but_in_moshier_routes_to_moshier(self):
        """Date outside JPL but inside Moshier should work with SEFLG_MOSEPH."""
        # ~500 CE - outside JPL (1550-2650) but inside Moshier (-3000 to +3000)
        jd = 1903682.5

        # Without SEFLG_MOSEPH, should raise EphemerisRangeError for JPL
        with pytest.raises(EphemerisRangeError):
            eph.swe_calc_ut(jd, SE_SUN, 0)

        # With SEFLG_MOSEPH, should work and return valid results
        pos, flag = eph.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Valid longitude


class TestMoshierFlagConstants:
    """Test that SEFLG_MOSEPH constant is properly defined."""

    def test_seflg_moseph_value(self):
        """SEFLG_MOSEPH should have correct value (4)."""
        assert SEFLG_MOSEPH == 4

    def test_seflg_moseph_distinct_from_other_ephemeris_flags(self):
        """SEFLG_MOSEPH should be distinct from SEFLG_JPLEPH and SEFLG_SWIEPH."""
        assert SEFLG_MOSEPH != SEFLG_JPLEPH
        assert SEFLG_MOSEPH != SEFLG_SWIEPH

        # And they should be different bits
        assert SEFLG_JPLEPH == 1
        assert SEFLG_SWIEPH == 2
        assert SEFLG_MOSEPH == 4
