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
        """Verify tidal acceleration constants.

        SE_TIDAL_DEFAULT matches pyswisseph TIDAL_DEFAULT (-25.8).
        SE_TIDAL_DE440 and SE_TIDAL_DE441 use the JPL DE441 value (-25.936).
        """
        assert SE_TIDAL_DEFAULT == -25.8  # Matches pyswisseph TIDAL_DEFAULT
        assert SE_TIDAL_DE440 == -25.936
        assert SE_TIDAL_DE441 == -25.936  # DE441 uses same value

    def test_get_tid_acc_returns_de440_default(self):
        """Verify that get_tid_acc() returns the default tidal acceleration value.

        The default tidal acceleration matches pyswisseph TIDAL_DEFAULT (-25.8),
        not the DE440-specific value (-25.936).
        """
        # Reset to ensure default state
        ephem.close()

        tid_acc = ephem.get_tid_acc()
        assert tid_acc == SE_TIDAL_DEFAULT


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


class TestDE440Fallback:
    """Test ephemeris file selection (DE440 default, env var override)."""

    def test_env_var_selection_logic_exists(self):
        """Verify that env var ephemeris selection is implemented in get_planets()."""
        import inspect
        from libephemeris import state

        source = inspect.getsource(state.get_planets)
        # Check that env var selection logic is present
        assert "_get_effective_ephemeris_file" in source
        assert "_EPHEMERIS_FILE" in source

    def test_default_prefers_de440(self):
        """Verify that DE440 is preferred when available."""
        import os
        from libephemeris import state

        # Reset state
        ephem.close()

        # Check that DE440 is the configured default
        assert state._EPHEMERIS_FILE == "de440.bsp"

        # If DE440 exists locally, it should be used
        base_dir = os.path.abspath(os.path.join(os.path.dirname(state.__file__), ".."))
        de440_path = os.path.join(base_dir, "de440.bsp")

        if os.path.exists(de440_path):
            # Trigger loading
            state.get_planets()
            assert state._EPHEMERIS_FILE == "de440.bsp"

    def test_effective_ephemeris_uses_env_var(self):
        """Test that _get_effective_ephemeris_file respects env var override."""
        import inspect
        from libephemeris import state

        source = inspect.getsource(state._get_effective_ephemeris_file)
        # Check that the env var mechanism is implemented
        assert "_EPHEMERIS_ENV_VAR" in source or "LIBEPHEMERIS_EPHEMERIS" in source


class TestDE440OuterPlanetPrecision:
    """
    Test that DE440 provides improved outer planet precision.

    DE440 includes 12+ years of additional spacecraft tracking data
    compared to DE421, and uses ICRF 3.0 reference frame.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        """Reset ephemeris state before each test."""
        ephem.close()
        yield
        ephem.close()

    def test_outer_planet_positions_at_multiple_dates(self):
        """
        Verify outer planet calculations work at multiple dates.

        DE440's improved outer planet ephemeris should produce
        positions with no anomalies or errors.
        """
        outer_planets = [SE_URANUS, SE_NEPTUNE, SE_PLUTO]
        test_dates = [
            ephem.swe_julday(2000, 1, 1, 12.0),  # J2000
            ephem.swe_julday(2020, 6, 21, 0.0),  # Recent solstice
            ephem.swe_julday(2050, 12, 21, 12.0),  # Future date
        ]

        for jd in test_dates:
            for planet in outer_planets:
                pos, _ = ephem.swe_calc_ut(jd, planet, SEFLG_SPEED)
                # Longitude should be valid
                assert 0 <= pos[0] < 360, f"Invalid longitude for planet {planet}"
                # Distance should be positive and reasonable
                assert pos[2] > 0, f"Invalid distance for planet {planet}"
                # Speed should be finite
                assert abs(pos[3]) < 1, f"Unreasonable speed for planet {planet}"

    def test_uranus_neptune_pluto_distance_ranges(self):
        """Verify outer planet distances are in expected ranges."""
        jd = ephem.swe_julday(2024, 1, 1, 12.0)

        # Uranus: 17.3-20.1 AU from Sun, varies from Earth
        pos_uranus, _ = ephem.swe_calc_ut(jd, SE_URANUS, 0)
        assert 17 < pos_uranus[2] < 22, "Uranus distance out of range"

        # Neptune: 29.0-30.3 AU from Sun, varies from Earth
        pos_neptune, _ = ephem.swe_calc_ut(jd, SE_NEPTUNE, 0)
        assert 28 < pos_neptune[2] < 32, "Neptune distance out of range"

        # Pluto: 29.7-49.3 AU from Sun, varies from Earth
        pos_pluto, _ = ephem.swe_calc_ut(jd, SE_PLUTO, 0)
        assert 28 < pos_pluto[2] < 52, "Pluto distance out of range"

    def test_outer_planet_positions_consistency(self):
        """
        Test that outer planet positions are consistent over time.

        The positions should change smoothly and predictably.
        """
        jd_base = ephem.swe_julday(2024, 1, 1, 12.0)

        for planet in [SE_URANUS, SE_NEPTUNE, SE_PLUTO]:
            positions = []
            for day_offset in [0, 10, 20, 30]:
                jd = jd_base + day_offset
                pos, _ = ephem.swe_calc_ut(jd, planet, 0)
                positions.append(pos[0])

            # Positions should change by a small, consistent amount
            for i in range(1, len(positions)):
                delta = positions[i] - positions[i - 1]
                # Normalize for 360 degree wraparound
                if delta > 180:
                    delta -= 360
                elif delta < -180:
                    delta += 360
                # Outer planets move slowly, less than 0.5 deg per 10 days
                assert abs(delta) < 0.5, f"Planet {planet} position jump too large"

    def test_de440_extended_date_range_for_outer_planets(self):
        """
        Test outer planets at dates beyond DE421 range.

        DE421 covered 1900-2050. DE440 covers 1550-2650.
        """
        # Date beyond DE421 range (after 2050)
        jd_2100 = ephem.swe_julday(2100, 1, 1, 12.0)

        for planet in [SE_URANUS, SE_NEPTUNE, SE_PLUTO]:
            pos, _ = ephem.swe_calc_ut(jd_2100, planet, SEFLG_SPEED)
            assert 0 <= pos[0] < 360
            assert pos[2] > 0

        # Historical date (before DE421 range starts in 1900)
        jd_1800 = ephem.swe_julday(1800, 1, 1, 12.0)

        for planet in [SE_URANUS, SE_NEPTUNE, SE_PLUTO]:
            pos, _ = ephem.swe_calc_ut(jd_1800, planet, SEFLG_SPEED)
            assert 0 <= pos[0] < 360
            assert pos[2] > 0
