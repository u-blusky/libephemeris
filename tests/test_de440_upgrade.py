"""
Tests for DE440 ephemeris upgrade.

These tests verify that the upgrade from DE421 to DE440 as the default
ephemeris file works correctly and provides the expected benefits.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_TIDAL_DE440,
    SE_TIDAL_DE441,
    SE_TIDAL_DEFAULT,
    SEFLG_SPEED,
)


class TestDE440Default:
    """Test that DE440 is correctly set as the default ephemeris."""

    def test_default_ephemeris_is_de440(self):
        """Verify that the default ephemeris file is de440.bsp."""
        from libephemeris.state import _EPHEMERIS_FILE

        assert _EPHEMERIS_FILE == "de440.bsp"

    def test_tidal_default_is_de440(self):
        """Verify that SE_TIDAL_DEFAULT matches DE440."""
        assert SE_TIDAL_DEFAULT == SE_TIDAL_DE440
        assert SE_TIDAL_DE440 == -25.936
        assert SE_TIDAL_DE441 == -25.936  # DE441 uses same value

    def test_get_tid_acc_returns_de440_default(self):
        """Verify that get_tid_acc() returns DE440 value by default."""
        # Reset to ensure default state
        ephem.close()

        tid_acc = ephem.get_tid_acc()
        assert tid_acc == SE_TIDAL_DE440


class TestDE440DateRange:
    """Test the extended date range of DE440."""

    def test_date_1550_works(self):
        """DE440 should support dates from 1550."""
        jd = ephem.swe_julday(1550, 6, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0  # Distance should be positive

    def test_date_2600_works(self):
        """DE440 should support dates until 2650."""
        jd = ephem.swe_julday(2600, 6, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0

    def test_historical_date_1600(self):
        """Test calculation at a historical date (year 1600)."""
        jd = ephem.swe_julday(1600, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_JUPITER, SEFLG_SPEED)
        assert 0 <= pos[0] < 360

    def test_future_date_2500(self):
        """Test calculation at a future date (year 2500)."""
        jd = ephem.swe_julday(2500, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_SATURN, SEFLG_SPEED)
        assert 0 <= pos[0] < 360


class TestDE440OuterPlanets:
    """Test outer planet calculations with DE440."""

    def test_jupiter_position(self):
        """Jupiter position calculation should work."""
        jd = ephem.swe_julday(2024, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_JUPITER, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0  # Distance should be positive

    def test_saturn_position(self):
        """Saturn position calculation should work."""
        jd = ephem.swe_julday(2024, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_SATURN, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0

    def test_uranus_position(self):
        """Uranus position calculation should work."""
        jd = ephem.swe_julday(2024, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_URANUS, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0

    def test_neptune_position(self):
        """Neptune position calculation should work."""
        jd = ephem.swe_julday(2024, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_NEPTUNE, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0

    def test_pluto_position(self):
        """Pluto position calculation should work."""
        jd = ephem.swe_julday(2024, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_PLUTO, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        assert pos[2] > 0


class TestDE440ModernDates:
    """Test that DE440 works correctly for modern dates."""

    def test_j2000_sun(self):
        """Test Sun position at J2000.0 epoch."""
        jd = 2451545.0  # J2000.0
        pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
        # Sun should be around 280° longitude at J2000
        assert 270 <= pos[0] <= 290

    def test_j2000_moon(self):
        """Test Moon position at J2000.0 epoch."""
        jd = 2451545.0  # J2000.0
        pos, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
        # Moon distance should be around 0.0025 AU (384,400 km)
        assert 0.002 < pos[2] < 0.003

    def test_modern_chart_calculation(self):
        """Test a full modern chart calculation works."""
        jd = ephem.swe_julday(2024, 6, 15, 12.0)

        # Calculate all major planets
        bodies = [
            SE_SUN,
            SE_MOON,
            SE_JUPITER,
            SE_SATURN,
            SE_URANUS,
            SE_NEPTUNE,
            SE_PLUTO,
        ]
        for body in bodies:
            pos, _ = ephem.swe_calc_ut(jd, body, SEFLG_SPEED)
            assert 0 <= pos[0] < 360, f"Body {body} has invalid longitude"


class TestDE440TidalConstants:
    """Test the tidal acceleration constants for DE440."""

    def test_se_tidal_de440_exists(self):
        """SE_TIDAL_DE440 constant should exist."""
        assert hasattr(ephem, "SE_TIDAL_DE440")
        assert ephem.SE_TIDAL_DE440 == -25.936

    def test_tidal_de440_alias_exists(self):
        """TIDAL_DE440 alias should exist."""
        assert hasattr(ephem, "TIDAL_DE440")
        assert ephem.TIDAL_DE440 == -25.936

    def test_set_tid_acc_de440(self):
        """set_tid_acc should accept DE440 value."""
        ephem.set_tid_acc(ephem.SE_TIDAL_DE440)
        assert ephem.get_tid_acc() == ephem.SE_TIDAL_DE440


class TestDE440VsDE421Compatibility:
    """Test that calculations remain compatible for dates in DE421 range."""

    def test_modern_date_still_works(self):
        """Dates within the old DE421 range should still work."""
        # 2000 was in DE421 range (1900-2050)
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
        assert 0 <= pos[0] < 360

    def test_early_1900s_works(self):
        """Early 1900s dates should work (were in DE421)."""
        jd = ephem.swe_julday(1920, 6, 21, 0.0)
        pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
        assert 0 <= pos[0] < 360

    def test_2040s_works(self):
        """2040s dates should work (were in DE421)."""
        jd = ephem.swe_julday(2045, 12, 31, 23.99)
        pos, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)
        assert 0 <= pos[0] < 360
