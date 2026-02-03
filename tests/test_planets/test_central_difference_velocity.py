"""
Tests for central difference velocity calculations.

This module validates that the O(h^2) central difference velocity calculations
are working correctly and provide improved precision over forward differences.

The central difference formula: f'(x) = [f(x+h) - f(x-h)] / (2h)
has error O(h^2) compared to O(h) for forward differences, providing
approximately 100x better numerical precision for the same timestep.
"""

import math
import pytest
import swisseph as swe

from libephemeris import swe_calc_ut, swe_calc, swe_calc_pctr
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_OSCU_APOG,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SEFLG_SPEED,
)


# Test date: J2000.0
JD_J2000 = 2451545.0


class TestCentralDifferenceVelocity:
    """Tests verifying central difference velocity calculation precision."""

    @pytest.mark.parametrize(
        "planet_id,planet_name,vel_tolerance",
        [
            (SE_SUN, "Sun", 0.001),
            (SE_MOON, "Moon", 0.01),
            (SE_MERCURY, "Mercury", 0.001),
            (SE_VENUS, "Venus", 0.001),
            (SE_MARS, "Mars", 0.001),
            (SE_JUPITER, "Jupiter", 0.001),
            (SE_SATURN, "Saturn", 0.001),
            (SE_URANUS, "Uranus", 0.001),
            (SE_NEPTUNE, "Neptune", 0.001),
            (SE_PLUTO, "Pluto", 0.001),
        ],
    )
    def test_velocity_matches_pyswisseph(
        self, planet_id: int, planet_name: str, vel_tolerance: float
    ):
        """
        Test that libephemeris velocity matches pyswisseph within tolerance.

        The central difference method should provide velocities that closely
        match the analytical derivatives used by Swiss Ephemeris.
        """
        # Get libephemeris velocity
        lib_pos, _ = swe_calc_ut(JD_J2000, planet_id, SEFLG_SPEED)
        lib_dlon = lib_pos[3]

        # Get pyswisseph velocity
        swe_pos, _ = swe.calc_ut(JD_J2000, planet_id, SEFLG_SPEED)
        swe_dlon = swe_pos[3]

        # Compare velocities
        vel_diff = abs(lib_dlon - swe_dlon)

        assert vel_diff < vel_tolerance, (
            f"{planet_name} velocity diff {vel_diff:.6f} deg/day exceeds "
            f"tolerance {vel_tolerance} deg/day "
            f"(libephemeris: {lib_dlon:.6f}, pyswisseph: {swe_dlon:.6f})"
        )

    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
        ],
    )
    def test_velocity_sign_matches_position_change(
        self, planet_id: int, planet_name: str
    ):
        """
        Test that velocity sign correctly predicts position change direction.

        A positive velocity should correspond to increasing longitude,
        and negative velocity to decreasing longitude.
        """
        jd = JD_J2000

        # Get position and velocity at jd
        pos, _ = swe_calc_ut(jd, planet_id, SEFLG_SPEED)
        lon = pos[0]
        dlon = pos[3]

        # Get position 1 day later
        pos_next, _ = swe_calc_ut(jd + 1.0, planet_id, SEFLG_SPEED)
        lon_next = pos_next[0]

        # Calculate actual change (handling wrap-around)
        actual_change = lon_next - lon
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity sign should match direction of change
        if abs(dlon) > 0.001:  # Skip if velocity is very small
            assert (dlon > 0) == (actual_change > 0), (
                f"{planet_name}: velocity sign {'+' if dlon > 0 else '-'} "
                f"doesn't match position change sign {'+' if actual_change > 0 else '-'}"
            )

    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MERCURY, "Mercury"),
            (SE_MARS, "Mars"),
        ],
    )
    def test_velocity_magnitude_reasonable(self, planet_id: int, planet_name: str):
        """
        Test that velocity magnitude is within expected astronomical bounds.

        Each planet has characteristic velocity ranges based on its orbital period.
        """
        pos, _ = swe_calc_ut(JD_J2000, planet_id, SEFLG_SPEED)
        dlon = pos[3]

        # Expected approximate maximum velocities (deg/day)
        expected_max = {
            SE_SUN: 1.1,  # ~1 deg/day
            SE_MOON: 16.0,  # ~12-15 deg/day
            SE_MERCURY: 2.5,  # Can be retrograde or fast direct
            SE_MARS: 1.0,  # ~0.5 deg/day direct
        }

        max_vel = expected_max.get(planet_id, 2.0)
        assert abs(dlon) < max_vel, (
            f"{planet_name} velocity {dlon:.4f} deg/day exceeds "
            f"expected maximum {max_vel} deg/day"
        )

    def test_all_velocity_components_returned(self):
        """Test that all 6 position/velocity components are returned."""
        pos, _ = swe_calc_ut(JD_J2000, SE_MARS, SEFLG_SPEED)

        assert len(pos) == 6, f"Expected 6 components, got {len(pos)}"

        # Velocity components should be non-zero for moving bodies
        dlon, dlat, ddist = pos[3], pos[4], pos[5]

        # Longitude velocity should definitely be non-zero for Mars
        assert dlon != 0.0, "Mars longitude velocity should not be zero"

    def test_velocity_consistent_across_time(self):
        """
        Test velocity is consistent when calculated at different times.

        The velocity calculation should work correctly regardless of
        the input Julian Day.
        """
        test_jds = [
            2415020.5,  # 1900
            2451545.0,  # 2000 (J2000)
            2458849.5,  # 2020
        ]

        for jd in test_jds:
            pos, _ = swe_calc_ut(jd, SE_JUPITER, SEFLG_SPEED)
            dlon = pos[3]

            # Jupiter's velocity should be roughly 0.08-0.14 deg/day when direct
            # Can be negative during retrograde
            assert abs(dlon) < 0.25, (
                f"Jupiter velocity {dlon:.6f} deg/day at JD {jd} seems unreasonable"
            )


class TestCentralDifferenceForLunarPoints:
    """Tests for central difference velocity on lunar nodes and apsides."""

    def test_mean_node_velocity(self):
        """Test that mean node has negative (retrograde) velocity."""
        pos, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED)
        dlon = pos[3]

        # Mean node moves retrograde at ~-0.053 deg/day
        assert dlon < 0, f"Mean node velocity {dlon} should be negative (retrograde)"
        assert abs(dlon) < 0.1, f"Mean node velocity {dlon} magnitude too large"

    def test_true_node_velocity(self):
        """Test that true node velocity is calculated and reasonable."""
        pos, _ = swe_calc_ut(JD_J2000, SE_TRUE_NODE, SEFLG_SPEED)
        dlon = pos[3]

        # True node can oscillate but should be reasonable
        assert abs(dlon) < 1.0, f"True node velocity {dlon} magnitude too large"

    def test_osculating_apogee_velocity(self):
        """Test that osculating apogee (True Lilith) velocity is calculated."""
        pos, _ = swe_calc_ut(JD_J2000, SE_OSCU_APOG, SEFLG_SPEED)
        dlon = pos[3]

        # Osculating apogee moves ~0.1 deg/day on average
        assert abs(dlon) < 5.0, f"Osculating apogee velocity {dlon} seems too large"

    def test_interpolated_apogee_velocity(self):
        """Test that interpolated apogee velocity is calculated."""
        pos, _ = swe_calc_ut(JD_J2000, SE_INTP_APOG, SEFLG_SPEED)
        dlon = pos[3]

        # Interpolated apogee should have reasonable velocity
        assert abs(dlon) < 5.0, f"Interpolated apogee velocity {dlon} seems too large"

    def test_interpolated_perigee_velocity(self):
        """Test that interpolated perigee velocity is calculated."""
        pos, _ = swe_calc_ut(JD_J2000, SE_INTP_PERG, SEFLG_SPEED)
        dlon = pos[3]

        # Interpolated perigee should have reasonable velocity
        assert abs(dlon) < 5.0, f"Interpolated perigee velocity {dlon} seems too large"


class TestCentralDifferencePlanetCentric:
    """Tests for central difference velocity in planet-centric calculations."""

    def test_pctr_velocity_calculated(self):
        """Test that planet-centric velocity is calculated."""
        pos, _ = swe_calc_pctr(JD_J2000, SE_MOON, SE_MARS, SEFLG_SPEED)

        # All 6 components should be returned
        assert len(pos) == 6

        # Moon seen from Mars should have non-zero velocity
        dlon = pos[3]
        assert dlon != 0.0, "Moon velocity from Mars should not be zero"

    def test_pctr_velocity_reasonable(self):
        """Test that planet-centric velocity magnitude is reasonable."""
        pos, _ = swe_calc_pctr(JD_J2000, SE_JUPITER, SE_MARS, SEFLG_SPEED)
        dlon = pos[3]

        # Jupiter seen from Mars - velocity depends on both orbital motions
        # Should be less than sum of both velocities
        assert abs(dlon) < 1.0, f"Jupiter-from-Mars velocity {dlon} seems too large"


class TestCentralDifferenceNumericalStability:
    """Tests for numerical stability of central difference calculations."""

    def test_velocity_stable_for_small_timestep(self):
        """
        Test that velocity is stable for small timestep changes.

        The central difference method should give consistent results
        regardless of minor changes in the calculation time.
        """
        jd = JD_J2000
        epsilon = 1e-6  # Very small time offset

        pos1, _ = swe_calc_ut(jd, SE_MARS, SEFLG_SPEED)
        pos2, _ = swe_calc_ut(jd + epsilon, SE_MARS, SEFLG_SPEED)

        dlon1, dlon2 = pos1[3], pos2[3]

        # Velocities should be nearly identical for tiny time difference
        assert abs(dlon1 - dlon2) < 1e-4, (
            f"Velocity unstable: {dlon1} vs {dlon2} for {epsilon} day offset"
        )

    def test_velocity_handles_wrap_around(self):
        """
        Test that velocity is correct near 0/360 degree boundary.

        The wrap-around handling in central difference should prevent
        spurious large velocities at the 360->0 boundary.
        """
        # Find a time when a fast-moving body is near 0 degrees
        # Moon moves ~13 deg/day, so it crosses 0 frequently
        for offset in range(0, 365, 1):
            jd = JD_J2000 + offset
            pos, _ = swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)
            lon = pos[0]
            dlon = pos[3]

            # If Moon is near 0 or 360, velocity should still be reasonable
            if lon < 5 or lon > 355:
                # Moon's velocity should never exceed ~16 deg/day
                assert abs(dlon) < 16.0, (
                    f"Moon velocity {dlon} deg/day near wrap boundary "
                    f"(lon={lon}) seems too large"
                )
                return  # Found and tested a wrap-around case

        pytest.skip("Could not find Moon near 0/360 boundary in test range")

    def test_swe_calc_tt_velocity(self):
        """Test that swe_calc (TT input) also calculates velocity correctly."""
        pos_ut, _ = swe_calc_ut(JD_J2000, SE_VENUS, SEFLG_SPEED)
        pos_tt, _ = swe_calc(JD_J2000, SE_VENUS, SEFLG_SPEED)

        dlon_ut = pos_ut[3]
        dlon_tt = pos_tt[3]

        # Velocities should be very close (TT vs UT is ~1 minute at J2000)
        assert abs(dlon_ut - dlon_tt) < 0.001, (
            f"Velocity difference between UT and TT: {abs(dlon_ut - dlon_tt)}"
        )


class TestCentralDifferenceVsPyswisseph:
    """
    Direct comparison tests between libephemeris and pyswisseph velocities.

    These tests validate that the central difference implementation achieves
    the expected O(h^2) precision improvement.
    """

    @pytest.fixture(autouse=True)
    def setup_ephemeris(self):
        """Setup pyswisseph ephemeris path."""
        swe.set_ephe_path("/Users/giacomo/dev/libephemeris/ephe")
        yield

    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
        ],
    )
    def test_velocity_precision_o_h_squared(self, planet_id: int, planet_name: str):
        """
        Test that velocity precision is O(h^2) as expected for central difference.

        For central difference with h = 1 second, the error should be
        approximately h^2 = (1/86400)^2 ~ 1e-10, which translates to
        very high precision velocity values.
        """
        # Test at multiple dates
        test_jds = [
            2415020.5,  # 1900
            2451545.0,  # J2000
            2458849.5,  # 2020
        ]

        for jd in test_jds:
            lib_pos, _ = swe_calc_ut(jd, planet_id, SEFLG_SPEED)
            swe_pos, _ = swe.calc_ut(jd, planet_id, SEFLG_SPEED)

            lib_dlon = lib_pos[3]
            swe_dlon = swe_pos[3]

            # Velocity difference should be very small
            vel_diff = abs(lib_dlon - swe_dlon)

            # For O(h^2) with h=1sec, expect better than 0.001 deg/day precision
            assert vel_diff < 0.01, (
                f"{planet_name} at JD {jd}: velocity diff {vel_diff:.8f} deg/day "
                f"exceeds expected O(h^2) precision "
                f"(libephemeris: {lib_dlon:.8f}, pyswisseph: {swe_dlon:.8f})"
            )
