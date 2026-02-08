"""
Moshier Ephemeris Range Validation Tests.

Validates that:
1. Moshier ephemeris (SEFLG_MOSEPH) correctly handles dates outside JPL/DE440
   range (~1549.5-2650.8), working for dates within its own range (~-2999 to +4999 CE)
2. JPL ephemeris (SEFLG_SWIEPH) correctly raises EphemerisRangeError for dates
   outside DE440 coverage

Test dates:
- Historical: -2999, -1000, 0, 500, 1000 CE (before DE440)
- Edge: 1549 mid-year, 2651 CE (just outside DE440)
- Extended: 2700, 4000 CE (far future, Moshier only)
"""

from __future__ import annotations

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_PLUTO,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SEFLG_MOSEPH,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
)
from libephemeris.exceptions import (
    EphemerisRangeError,
    MOSHIER_JD_START,
    MOSHIER_JD_END,
)


# =============================================================================
# TEST CONFIGURATIONS
# =============================================================================

# Dates outside DE440 range (~1549.5-2650.8) but inside Moshier range (~-2999 to +4999)
# Note: DE440 actual boundaries are approximately JD 2287184.5 to 2688976.5
# Moshier boundaries: JD 625673.5 to 3182395.5 (~-2999 to ~+4999)
DATES_OUTSIDE_DE440_INSIDE_MOSHIER = [
    (-2999, 1, 1, 12.0, "-2999 (near Moshier start)"),
    (-1000, 6, 15, 12.0, "-1000 BCE"),
    (0, 1, 1, 12.0, "0 CE (1 BCE)"),
    (500, 6, 21, 12.0, "500 CE (early medieval)"),
    (1000, 12, 25, 12.0, "1000 CE (medieval)"),
    (1549, 6, 1, 12.0, "1549 mid-year (just before DE440)"),
    (2651, 1, 1, 12.0, "2651 (just after DE440)"),
    (2700, 6, 15, 12.0, "2700 (far future)"),
    (4000, 12, 31, 12.0, "4000 (near Moshier end)"),
]

# Dates inside DE440 range (should work with both SEFLG_SWIEPH and SEFLG_MOSEPH)
DATES_INSIDE_DE440 = [
    (1550, 1, 1, 12.0, "1550 (DE440 start)"),
    (1800, 6, 15, 12.0, "1800"),
    (2000, 1, 1, 12.0, "2000 (J2000.0)"),
    (2400, 6, 15, 12.0, "2400"),
    (2650, 1, 1, 12.0, "2650 (DE440 end)"),
]

# Dates outside both DE440 and Moshier ranges
# Moshier range: JD 625673.5 (~-2999) to 3182395.5 (~+4999)
DATES_OUTSIDE_MOSHIER = [
    (-4000, 1, 1, 12.0, "-4000 (before Moshier)"),
    (-3001, 6, 15, 12.0, "-3001 (just before Moshier)"),
    (5000, 1, 1, 12.0, "5000 (just after Moshier)"),
    (6000, 6, 15, 12.0, "6000 (far future)"),
]

# Planets to test with
TEST_PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_PLUTO, "Pluto"),
]


# =============================================================================
# TEST CLASSES
# =============================================================================


class TestMoshierWorksOutsideDE440Range:
    """Verify SEFLG_MOSEPH works for dates outside DE440 but inside Moshier range."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_OUTSIDE_DE440_INSIDE_MOSHIER
    )
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_moshier_returns_valid_position(
        self, year, month, day, hour, date_desc, planet_id, planet_name
    ):
        """SEFLG_MOSEPH should return valid positions for dates outside DE440."""
        jd = ephem.swe_julday(year, month, day, hour)

        # With SEFLG_MOSEPH, should work
        pos, flag = ephem.swe_calc_ut(jd, planet_id, SEFLG_MOSEPH)

        # Basic validation
        assert len(pos) == 6, f"{planet_name} at {date_desc}: expected 6 elements"
        assert 0 <= pos[0] < 360, (
            f"{planet_name} at {date_desc}: invalid longitude {pos[0]}"
        )
        assert -90 <= pos[1] <= 90, (
            f"{planet_name} at {date_desc}: invalid latitude {pos[1]}"
        )
        assert pos[2] >= 0, f"{planet_name} at {date_desc}: invalid distance {pos[2]}"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_OUTSIDE_DE440_INSIDE_MOSHIER
    )
    def test_moshier_with_speed_flag(self, year, month, day, hour, date_desc):
        """SEFLG_MOSEPH | SEFLG_SPEED should return valid velocities."""
        jd = ephem.swe_julday(year, month, day, hour)

        pos, flag = ephem.swe_calc_ut(jd, SE_MARS, SEFLG_MOSEPH | SEFLG_SPEED)

        # Position checks
        assert len(pos) == 6
        assert 0 <= pos[0] < 360

        # Velocity should be reasonable (Mars: ~0.3-0.8 deg/day typical)
        assert abs(pos[3]) < 2.0, f"Mars velocity at {date_desc} unreasonable: {pos[3]}"


class TestJPLRaisesErrorOutsideDE440Range:
    """Verify SEFLG_SWIEPH raises EphemerisRangeError for dates outside DE440."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_OUTSIDE_DE440_INSIDE_MOSHIER
    )
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_swieph_raises_ephemeris_range_error(
        self, year, month, day, hour, date_desc, planet_id, planet_name
    ):
        """SEFLG_SWIEPH should raise EphemerisRangeError for dates outside DE440."""
        jd = ephem.swe_julday(year, month, day, hour)

        with pytest.raises(EphemerisRangeError) as exc_info:
            ephem.swe_calc_ut(jd, planet_id, SEFLG_SWIEPH)

        err = exc_info.value
        assert err.requested_jd == jd
        # Should mention JPL or DE440, not Moshier
        assert (
            "Moshier" not in err.message
            or "JPL" in err.message
            or "DE440" in err.message
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_OUTSIDE_DE440_INSIDE_MOSHIER
    )
    def test_default_flags_raise_error(self, year, month, day, hour, date_desc):
        """Default flags (0) should use JPL and raise error for dates outside DE440."""
        jd = ephem.swe_julday(year, month, day, hour)

        with pytest.raises(EphemerisRangeError):
            ephem.swe_calc_ut(jd, SE_SUN, 0)


class TestBothEphemerisWorkInsideDE440Range:
    """Verify both JPL and Moshier work for dates inside DE440 range."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_INSIDE_DE440)
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_swieph_works(
        self, year, month, day, hour, date_desc, planet_id, planet_name
    ):
        """SEFLG_SWIEPH should work for dates inside DE440 range."""
        jd = ephem.swe_julday(year, month, day, hour)

        pos, flag = ephem.swe_calc_ut(jd, planet_id, SEFLG_SWIEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_INSIDE_DE440)
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_moseph_works(
        self, year, month, day, hour, date_desc, planet_id, planet_name
    ):
        """SEFLG_MOSEPH should also work for dates inside DE440 range."""
        jd = ephem.swe_julday(year, month, day, hour)

        pos, flag = ephem.swe_calc_ut(jd, planet_id, SEFLG_MOSEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360


class TestMoshierRaisesErrorOutsideMoshierRange:
    """Verify SEFLG_MOSEPH raises EphemerisRangeError for dates outside Moshier range."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_OUTSIDE_MOSHIER)
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_moseph_raises_ephemeris_range_error(
        self, year, month, day, hour, date_desc, planet_id, planet_name
    ):
        """SEFLG_MOSEPH should raise EphemerisRangeError for dates outside Moshier range."""
        jd = ephem.swe_julday(year, month, day, hour)

        with pytest.raises(EphemerisRangeError) as exc_info:
            ephem.swe_calc_ut(jd, planet_id, SEFLG_MOSEPH)

        err = exc_info.value
        assert err.requested_jd == jd
        # Should mention Moshier
        assert "Moshier" in err.message


class TestMoshierRangeConstants:
    """Verify Moshier range constants are correct."""

    def test_moshier_start_jd_value(self):
        """Moshier start JD should be 625673.5 (approximately -2999 CE)."""
        # MOSHIER_JD_START = 625673.5 corresponds to approximately -2999 CE
        assert MOSHIER_JD_START == 625673.5

        # Verify -2999 Jan 1 is inside range
        jd_neg2999 = ephem.swe_julday(-2999, 1, 1, 0.0)
        assert jd_neg2999 >= MOSHIER_JD_START

    def test_moshier_end_jd_value(self):
        """Moshier end JD should be 3182395.5 (approximately year 4000 CE)."""
        # MOSHIER_JD_END = 3182395.5 corresponds to approximately year 4000 CE
        assert MOSHIER_JD_END == 3182395.5

        # Verify year 4000 Jan 1 is inside range
        jd_4000 = ephem.swe_julday(4000, 1, 1, 12.0)
        assert jd_4000 <= MOSHIER_JD_END

        # Verify year 4001 is outside range
        jd_4001 = ephem.swe_julday(4001, 1, 1, 12.0)
        assert jd_4001 > MOSHIER_JD_END


# Dates valid for lunar nodes in Moshier mode
# Mean Node uses Meeus polynomial: valid for years ~0-4000 CE
# True Node requires Moon position - but in Moshier mode uses Moshier's ELP
# which should work for the full Moshier range, but True Node still needs
# the Sun position from Skyfield for the full calculation, so it's limited
# to DE440 range in practice.
DATES_FOR_MEAN_NODE_MOSHIER = [
    (500, 6, 21, 12.0, "500 CE (early medieval)"),
    (1000, 12, 25, 12.0, "1000 CE (medieval)"),
    (1549, 6, 1, 12.0, "1549 mid-year"),
    (2651, 1, 1, 12.0, "2651"),
    (2700, 6, 15, 12.0, "2700 (far future)"),
    (3000, 6, 15, 12.0, "3000"),
]

# True Node in Moshier mode: still limited to DE440 range because calc_true_lunar_node
# uses Skyfield for the Moon position internally
DATES_FOR_TRUE_NODE_MOSHIER = [
    (1550, 6, 1, 12.0, "1550 (in DE440)"),
    (2000, 1, 1, 12.0, "2000 (J2000)"),
    (2650, 1, 1, 12.0, "2650 (near DE440 end)"),
]


class TestLunarNodesInMoshierRange:
    """Verify lunar nodes work correctly in Moshier mode.

    Note: Mean Node uses Meeus polynomial (valid ~0-4000 CE).
    True Node requires Moon position from Skyfield (limited to DE440 range).
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_FOR_MEAN_NODE_MOSHIER
    )
    def test_mean_node_moshier(self, year, month, day, hour, date_desc):
        """Mean lunar node should work in Moshier mode for supported dates."""
        jd = ephem.swe_julday(year, month, day, hour)

        pos, flag = ephem.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_MOSEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360
        # Latitude should be 0 for mean node
        assert abs(pos[1]) < 0.001

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_FOR_TRUE_NODE_MOSHIER
    )
    def test_true_node_moshier(self, year, month, day, hour, date_desc):
        """True lunar node should work in Moshier mode for dates in DE440 range.

        Note: True Node calculation in Moshier mode still uses Skyfield for Moon
        position, so it's limited to DE440 date range.
        """
        jd = ephem.swe_julday(year, month, day, hour)

        pos, flag = ephem.swe_calc_ut(jd, SE_TRUE_NODE, SEFLG_MOSEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360


class TestSweCalcMoshierRange:
    """Test range handling for swe_calc (TT time) in addition to swe_calc_ut."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", DATES_OUTSIDE_DE440_INSIDE_MOSHIER[:3]
    )
    def test_swe_calc_moshier_works(self, year, month, day, hour, date_desc):
        """swe_calc with SEFLG_MOSEPH should work for extended dates."""
        jd = ephem.swe_julday(year, month, day, hour)

        pos, flag = ephem.swe_calc(jd, SE_SUN, SEFLG_MOSEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_OUTSIDE_MOSHIER[:2])
    def test_swe_calc_moshier_raises_error(self, year, month, day, hour, date_desc):
        """swe_calc with SEFLG_MOSEPH should raise error outside Moshier range."""
        jd = ephem.swe_julday(year, month, day, hour)

        with pytest.raises(EphemerisRangeError) as exc_info:
            ephem.swe_calc(jd, SE_SUN, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in err.message


class TestRangeEdgeCases:
    """Test edge cases at range boundaries."""

    @pytest.mark.comparison
    def test_exactly_at_moshier_start(self):
        """Date exactly at Moshier start should work."""
        pos, flag = ephem.swe_calc_ut(MOSHIER_JD_START, SE_SUN, SEFLG_MOSEPH)
        assert len(pos) == 6

    @pytest.mark.comparison
    def test_exactly_at_moshier_end(self):
        """Date exactly at Moshier end should work."""
        pos, flag = ephem.swe_calc_ut(MOSHIER_JD_END, SE_SUN, SEFLG_MOSEPH)
        assert len(pos) == 6

    @pytest.mark.comparison
    def test_just_before_moshier_start_fails(self):
        """Date just before Moshier start should fail."""
        jd = MOSHIER_JD_START - 1.0

        with pytest.raises(EphemerisRangeError):
            ephem.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

    @pytest.mark.comparison
    def test_just_after_moshier_end_fails(self):
        """Date just after Moshier end should fail."""
        jd = MOSHIER_JD_END + 1.0

        with pytest.raises(EphemerisRangeError):
            ephem.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)
