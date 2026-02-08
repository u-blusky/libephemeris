"""
Tests for Moshier mode support in house calculations.

This module tests the propagation of SEFLG_MOSEPH flag through
swe_houses() and swe_houses_ex() for extended date range support.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_MOSEPH,
    SEFLG_SIDEREAL,
    SEFLG_SPEED,
    SE_SUN,
)
from libephemeris.exceptions import (
    EphemerisRangeError,
    MOSHIER_JD_START,
    MOSHIER_JD_END,
)


# Test dates outside DE440 range but within Moshier range
# DE440 covers ~1550-2650, Moshier covers -3000 to +3000
JD_YEAR_1000 = 2086667.5  # ~1000 CE (outside DE440 range)
JD_YEAR_1200 = 2159712.5  # ~1200 CE (outside DE440 range)
JD_YEAR_3000 = 2816787.5  # ~3000 CE (outside DE440 range)

# Standard test date within DE440 range for comparison
JD_J2000 = 2451545.0  # 2000-01-01

# Test location
ROME_LAT = 41.9
ROME_LON = 12.5


class TestSweHousesMoshierFlag:
    """Test swe_houses() with SEFLG_MOSEPH flag."""

    @pytest.mark.unit
    def test_swe_houses_accepts_iflag_parameter(self):
        """swe_houses() should accept optional iflag parameter."""
        # Should work without iflag
        cusps1, ascmc1 = ephem.swe_houses(JD_J2000, ROME_LAT, ROME_LON, ord("P"))

        # Should work with iflag=0
        cusps2, ascmc2 = ephem.swe_houses(JD_J2000, ROME_LAT, ROME_LON, ord("P"), 0)

        # Results should be essentially identical
        for i in range(12):
            assert abs(cusps1[i] - cusps2[i]) < 0.0001

    @pytest.mark.unit
    def test_swe_houses_moshier_within_de440_range(self):
        """swe_houses() with SEFLG_MOSEPH should work within DE440 range."""
        # Calculate with JPL (default)
        cusps_jpl, ascmc_jpl = ephem.swe_houses(JD_J2000, ROME_LAT, ROME_LON, ord("P"))

        # Calculate with Moshier
        cusps_mosh, ascmc_mosh = ephem.swe_houses(
            JD_J2000, ROME_LAT, ROME_LON, ord("P"), SEFLG_MOSEPH
        )

        # Both should return valid cusps (12 houses)
        assert len(cusps_jpl) == 12
        assert len(cusps_mosh) == 12

        # Results should be similar (within a few arcminutes)
        # Moshier is less precise than JPL but should be close
        for i in range(12):
            diff = abs(cusps_jpl[i] - cusps_mosh[i])
            # Handle 360 wraparound
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.1, f"House {i + 1} differs by {diff:.4f}°"

    @pytest.mark.unit
    def test_swe_houses_moshier_ancient_date(self):
        """swe_houses() with SEFLG_MOSEPH should work for ancient dates."""
        # This date is outside DE440 range (~1000 CE)
        cusps, ascmc = ephem.swe_houses(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), SEFLG_MOSEPH
        )

        # Should return valid cusps
        assert len(cusps) == 12
        assert len(ascmc) == 8

        # All cusps should be in valid range
        for i, cusp in enumerate(cusps):
            assert 0.0 <= cusp < 360.0, f"House {i + 1} cusp {cusp} out of range"

    @pytest.mark.unit
    def test_swe_houses_moshier_far_future_date(self):
        """swe_houses() with SEFLG_MOSEPH should work for far future dates."""
        # This date is outside DE440 range (~3000 CE)
        cusps, ascmc = ephem.swe_houses(
            JD_YEAR_3000, ROME_LAT, ROME_LON, ord("E"), SEFLG_MOSEPH
        )

        # Should return valid cusps (Equal houses)
        assert len(cusps) == 12
        assert len(ascmc) == 8

        # Equal houses: each cusp should be 30° apart
        for i in range(11):
            diff = (cusps[i + 1] - cusps[i]) % 360.0
            assert abs(diff - 30.0) < 0.001, f"House spacing incorrect: {diff}°"

    @pytest.mark.unit
    def test_swe_houses_jpl_sunshine_fallback_outside_de440_range(self):
        """swe_houses() Sunshine houses fallback to sun_dec=0 outside DE440 range."""
        # Sunshine houses ('I') requires Sun position calculation
        # Without SEFLG_MOSEPH, this falls back to sun_dec=0.0 (equinox behavior)
        # This is by design - the function has exception handling to avoid failures
        cusps, ascmc = ephem.swe_houses(JD_YEAR_1000, ROME_LAT, ROME_LON, ord("I"), 0)

        # Should still return valid cusps (with fallback sun_dec=0)
        assert len(cusps) == 12
        assert len(ascmc) == 8

    @pytest.mark.unit
    def test_swe_houses_moshier_sunshine_houses(self):
        """swe_houses() with SEFLG_MOSEPH should work for Sunshine houses outside DE440."""
        # Sunshine houses require Sun position - this tests the internal swe_calc_ut propagation
        cusps, ascmc = ephem.swe_houses(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("I"), SEFLG_MOSEPH
        )

        # Should return valid cusps
        assert len(cusps) == 12
        assert len(ascmc) == 8

    @pytest.mark.unit
    def test_swe_houses_multiple_systems_with_moshier(self):
        """swe_houses() with SEFLG_MOSEPH should work for various house systems."""
        house_systems = ["P", "K", "R", "C", "E", "W", "O", "B", "M", "T"]

        for hsys in house_systems:
            cusps, ascmc = ephem.swe_houses(
                JD_YEAR_1200, ROME_LAT, ROME_LON, ord(hsys), SEFLG_MOSEPH
            )
            assert len(cusps) == 12, f"House system {hsys} returned wrong cusp count"
            assert len(ascmc) == 8, f"House system {hsys} returned wrong ascmc count"


class TestSweHousesExMoshierFlag:
    """Test swe_houses_ex() with SEFLG_MOSEPH flag."""

    @pytest.mark.unit
    def test_swe_houses_ex_moshier_propagation(self):
        """swe_houses_ex() should propagate SEFLG_MOSEPH to swe_houses()."""
        # This should work for ancient dates with SEFLG_MOSEPH
        cusps, ascmc = ephem.swe_houses_ex(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), SEFLG_MOSEPH
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8

    @pytest.mark.unit
    def test_swe_houses_ex_sidereal_with_moshier(self):
        """swe_houses_ex() should work with both SEFLG_SIDEREAL and SEFLG_MOSEPH."""
        # Set up sidereal mode (Lahiri)
        ephem.swe_set_sid_mode(1)  # SE_SIDM_LAHIRI

        flags = SEFLG_SIDEREAL | SEFLG_MOSEPH
        cusps, ascmc = ephem.swe_houses_ex(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), flags
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8

    @pytest.mark.unit
    def test_swe_houses_ex_moshier_sunshine_sidereal(self):
        """swe_houses_ex() with SEFLG_MOSEPH should work for Sunshine houses in sidereal mode."""
        ephem.swe_set_sid_mode(1)  # SE_SIDM_LAHIRI

        flags = SEFLG_SIDEREAL | SEFLG_MOSEPH
        cusps, ascmc = ephem.swe_houses_ex(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("I"), flags
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8


class TestSweHousesEx2MoshierFlag:
    """Test swe_houses_ex2() with SEFLG_MOSEPH flag."""

    @pytest.mark.unit
    def test_swe_houses_ex2_moshier_with_speed(self):
        """swe_houses_ex2() should propagate SEFLG_MOSEPH and calculate velocities."""
        flags = SEFLG_MOSEPH | SEFLG_SPEED
        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), flags
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert len(cusps_speed) == 12
        assert len(ascmc_speed) == 8

        # Velocities should be non-zero when SEFLG_SPEED is set
        # Ascendant velocity is typically around 360°/day (diurnal motion)
        assert cusps_speed[0] != 0.0

    @pytest.mark.unit
    def test_swe_houses_ex2_moshier_without_speed(self):
        """swe_houses_ex2() with SEFLG_MOSEPH but without SEFLG_SPEED returns zero velocities."""
        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("E"), SEFLG_MOSEPH
        )

        assert len(cusps) == 12
        assert len(cusps_speed) == 12

        # Without SEFLG_SPEED, velocities should be zero
        for v in cusps_speed:
            assert v == 0.0


class TestMoshierHousesConsistency:
    """Test consistency between different function variants with Moshier."""

    @pytest.mark.unit
    def test_swe_houses_and_swe_houses_ex_consistency(self):
        """swe_houses() and swe_houses_ex() should return same results with SEFLG_MOSEPH."""
        # Without sidereal flag, results should be identical
        cusps1, ascmc1 = ephem.swe_houses(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), SEFLG_MOSEPH
        )
        cusps2, ascmc2 = ephem.swe_houses_ex(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), SEFLG_MOSEPH
        )

        for i in range(12):
            assert abs(cusps1[i] - cusps2[i]) < 0.0001

    @pytest.mark.unit
    def test_moshier_houses_reasonable_values(self):
        """Moshier houses should produce astronomically reasonable values."""
        cusps, ascmc = ephem.swe_houses(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("P"), SEFLG_MOSEPH
        )

        asc = ascmc[0]
        mc = ascmc[1]
        armc = ascmc[2]
        vertex = ascmc[3]

        # All angles should be in valid range
        assert 0.0 <= asc < 360.0
        assert 0.0 <= mc < 360.0
        assert 0.0 <= armc < 360.0
        assert 0.0 <= vertex < 360.0

        # MC and IC should be opposite (roughly)
        ic_expected = (mc + 180.0) % 360.0
        # Cusp 4 is IC
        assert (
            abs(cusps[3] - ic_expected) < 1.0 or abs(cusps[3] - ic_expected - 360) < 1.0
        )


class TestMoshierEdgeCases:
    """Test edge cases for Moshier mode in houses."""

    @pytest.mark.unit
    def test_moshier_at_range_boundary(self):
        """Houses should work at Moshier range boundaries."""
        # Near start of Moshier range
        jd_near_start = MOSHIER_JD_START + 1.0
        cusps, ascmc = ephem.swe_houses(
            jd_near_start, ROME_LAT, ROME_LON, ord("E"), SEFLG_MOSEPH
        )
        assert len(cusps) == 12

        # Near end of Moshier range
        jd_near_end = MOSHIER_JD_END - 1.0
        cusps, ascmc = ephem.swe_houses(
            jd_near_end, ROME_LAT, ROME_LON, ord("E"), SEFLG_MOSEPH
        )
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_moshier_gauquelin_sectors(self):
        """Gauquelin sectors (36) should work with Moshier."""
        cusps, ascmc = ephem.swe_houses(
            JD_YEAR_1000, ROME_LAT, ROME_LON, ord("G"), SEFLG_MOSEPH
        )

        # Gauquelin returns 36 sectors
        assert len(cusps) == 36
        assert len(ascmc) == 8
