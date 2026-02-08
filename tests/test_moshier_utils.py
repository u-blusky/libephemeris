"""
Tests for moshier/utils.py - Aberration and light-time correction functions.

This module tests the standalone aberration and light-time correction functions
used by the Moshier ephemeris for computing apparent planetary positions.
"""

from __future__ import annotations

import math

import pytest

from libephemeris.moshier.utils import (
    # Constants
    C_LIGHT_AU_DAY,
    DEG_TO_RAD,
    RAD_TO_DEG,
    ARCSEC_TO_RAD,
    J2000,
    # Light-time functions
    light_time_correction,
    iterative_light_time_correction,
    # Aberration functions
    annual_aberration,
    annual_aberration_cartesian,
    apply_aberration_to_position,
    aberration_in_longitude_latitude,
    apply_light_time_and_aberration,
    compute_earth_velocity,
    # Coordinate conversions
    spherical_to_cartesian,
    cartesian_to_spherical,
    spherical_to_cartesian_velocity,
    normalize_angle,
)


# =============================================================================
# LIGHT-TIME CORRECTION TESTS
# =============================================================================


class TestLightTimeCorrection:
    """Tests for light_time_correction function."""

    @pytest.mark.unit
    def test_light_time_one_au(self):
        """Light-time for 1 AU should be ~499 seconds (~0.00578 days)."""
        lt = light_time_correction(1.0)
        # 1 AU / c = 499.004784 seconds = 0.00577551 days
        expected_days = 499.004784 / 86400.0
        assert abs(lt - expected_days) < 0.0001

    @pytest.mark.unit
    def test_light_time_zero_distance(self):
        """Light-time for zero distance should be zero."""
        lt = light_time_correction(0.0)
        assert lt == 0.0

    @pytest.mark.unit
    def test_light_time_saturn(self):
        """Light-time for Saturn (~10 AU) should be ~83 minutes."""
        # Saturn average distance ~9.5 AU, light-time ~79 minutes
        lt = light_time_correction(10.0)
        lt_minutes = lt * 24 * 60
        assert 80 < lt_minutes < 90  # ~83 minutes for 10 AU


class TestIterativeLightTimeCorrection:
    """Tests for iterative_light_time_correction function."""

    @pytest.mark.unit
    def test_iterative_convergence(self):
        """Light-time correction should converge with iterations."""
        jd_tt = J2000

        # Mock Earth position at origin
        earth_pos = (0.0, 0.0, 0.0)

        # Mock target function returning fixed position at 5 AU
        def mock_target(jd):
            return (5.0, 0.0, 0.0)

        corrected_pos, lt = iterative_light_time_correction(
            jd_tt, earth_pos, mock_target, iterations=3
        )

        # Distance should be 5 AU
        dist = math.sqrt(sum(x * x for x in corrected_pos))
        assert abs(dist - 5.0) < 1e-10

        # Light-time should be 5 AU / c
        expected_lt = 5.0 / C_LIGHT_AU_DAY
        assert abs(lt - expected_lt) < 1e-10

    @pytest.mark.unit
    def test_iterative_with_moving_target(self):
        """Light-time iteration should work with a moving target."""
        jd_tt = J2000
        earth_pos = (0.0, 0.0, 0.0)

        # Target moves 0.01 AU/day in x direction
        initial_x = 2.0
        velocity = 0.01  # AU/day

        def moving_target(jd):
            dt = jd - jd_tt
            return (initial_x + velocity * dt, 0.0, 0.0)

        corrected_pos, lt = iterative_light_time_correction(
            jd_tt, earth_pos, moving_target, iterations=3
        )

        # Light-time is approximately 2 AU / c
        expected_lt = 2.0 / C_LIGHT_AU_DAY

        # Position should be slightly less than 2 AU because we look back in time
        # When we look back by light_time, target was at x = 2 - velocity * lt
        expected_x = initial_x - velocity * lt

        assert abs(corrected_pos[0] - expected_x) < 0.001
        assert abs(lt - expected_lt) < 0.001

    @pytest.mark.unit
    def test_single_iteration_vs_multiple(self):
        """More iterations should improve accuracy."""
        jd_tt = J2000
        earth_pos = (0.0, 0.0, 0.0)

        # Fast-moving target to exaggerate light-time effect
        initial_x = 5.0
        velocity = 0.1  # AU/day

        def moving_target(jd):
            dt = jd - jd_tt
            return (initial_x + velocity * dt, 0.0, 0.0)

        pos_1iter, lt_1iter = iterative_light_time_correction(
            jd_tt, earth_pos, moving_target, iterations=1
        )
        pos_3iter, lt_3iter = iterative_light_time_correction(
            jd_tt, earth_pos, moving_target, iterations=3
        )

        # Results should differ slightly (3 iterations is more accurate)
        # The difference is small but measurable
        diff = abs(pos_1iter[0] - pos_3iter[0])
        assert diff > 0  # There is a difference
        assert diff < 0.01  # But not huge


# =============================================================================
# ABERRATION TESTS
# =============================================================================


class TestAnnualAberration:
    """Tests for annual_aberration function (spherical coordinates)."""

    @pytest.mark.unit
    def test_aberration_magnitude(self):
        """Maximum aberration should be ~20.5 arcsec."""
        # Maximum aberration occurs when object is perpendicular to Earth's motion
        # Earth orbital velocity ~30 km/s = ~0.0172 AU/day

        # Object at ecliptic longitude 0, latitude 0
        # Earth at longitude 270 (Sun at 90), velocity toward longitude 0
        lon, lat = 0.0, 0.0

        # Earth position and velocity (simplified: circular orbit)
        # When Sun is at lon=90, Earth is at lon=270, moving toward lon=0
        earth_x, earth_y, earth_z = spherical_to_cartesian(270.0, 0.0, 1.0)
        earth_vx = 0.0172  # AU/day toward +x
        earth_vy = 0.0
        earth_vz = 0.0

        d_lon, d_lat = annual_aberration(
            lon, lat, earth_x, earth_y, earth_z, earth_vx, earth_vy, earth_vz
        )

        # Maximum aberration ~20.5 arcsec = 0.0057 degrees
        assert abs(d_lon) < 0.01  # degrees
        assert abs(d_lat) < 0.01  # degrees

        # Convert to arcsec for comparison
        d_lon_arcsec = d_lon * 3600
        assert abs(d_lon_arcsec) < 25  # Should be ~20 arcsec or less


class TestAnnualAberrationCartesian:
    """Tests for annual_aberration_cartesian function."""

    @pytest.mark.unit
    def test_aberration_perpendicular_velocity(self):
        """Aberration is maximum when velocity is perpendicular to line of sight."""
        # Object in +x direction
        target_dir = (1.0, 0.0, 0.0)
        # Earth velocity in +y direction (perpendicular)
        # Typical velocity ~30 km/s = ~0.0172 AU/day
        v_earth = 0.0172
        earth_vel = (0.0, v_earth, 0.0)

        dx, dy, dz = annual_aberration_cartesian(target_dir, earth_vel)

        # Aberration should be primarily in y direction
        # Magnitude: v/c = 0.0172 / 173.14 = ~10^-4 radians = ~20 arcsec
        expected_dy = v_earth / C_LIGHT_AU_DAY
        assert abs(dx) < 1e-10  # No x component (r·v = 0)
        assert abs(dy - expected_dy) < 1e-10
        assert abs(dz) < 1e-10

    @pytest.mark.unit
    def test_aberration_parallel_velocity(self):
        """Aberration is zero when velocity is parallel to line of sight."""
        # Object in +x direction
        target_dir = (1.0, 0.0, 0.0)
        # Earth velocity also in +x direction (parallel)
        v_earth = 0.0172
        earth_vel = (v_earth, 0.0, 0.0)

        dx, dy, dz = annual_aberration_cartesian(target_dir, earth_vel)

        # Aberration should be zero (r·v/c = v/c, so Δr = v/c - v/c*r = 0)
        assert abs(dx) < 1e-10
        assert abs(dy) < 1e-10
        assert abs(dz) < 1e-10

    @pytest.mark.unit
    def test_aberration_constant(self):
        """Verify aberration constant κ ≈ 20.49552 arcsec."""
        # When velocity is perpendicular, aberration = v/c in radians
        # For Earth's orbital velocity, this is the aberration constant

        # Earth's mean orbital velocity: 29.78 km/s
        # In AU/day: 29.78 * 86400 / 149597870.7 = 0.01720 AU/day
        v_earth = 0.01720  # AU/day (approximate mean)

        target_dir = (1.0, 0.0, 0.0)
        earth_vel = (0.0, v_earth, 0.0)

        dx, dy, dz = annual_aberration_cartesian(target_dir, earth_vel)

        # Convert to arcsec
        aberration_rad = math.sqrt(dx * dx + dy * dy + dz * dz)
        aberration_arcsec = aberration_rad / ARCSEC_TO_RAD

        # Should be close to 20.49552 arcsec
        assert 19 < aberration_arcsec < 21


class TestApplyAberrationToPosition:
    """Tests for apply_aberration_to_position function."""

    @pytest.mark.unit
    def test_preserves_distance(self):
        """Aberration should preserve distance (only change direction)."""
        position = (5.0, 3.0, 1.0)
        earth_vel = (0.0, 0.0172, 0.0)

        apparent = apply_aberration_to_position(position, earth_vel)

        original_dist = math.sqrt(sum(x * x for x in position))
        apparent_dist = math.sqrt(sum(x * x for x in apparent))

        assert abs(original_dist - apparent_dist) < 1e-10

    @pytest.mark.unit
    def test_zero_velocity(self):
        """With zero velocity, position should be unchanged."""
        position = (5.0, 3.0, 1.0)
        earth_vel = (0.0, 0.0, 0.0)

        apparent = apply_aberration_to_position(position, earth_vel)

        assert abs(apparent[0] - position[0]) < 1e-10
        assert abs(apparent[1] - position[1]) < 1e-10
        assert abs(apparent[2] - position[2]) < 1e-10


class TestAberrationInLongitudeLatitude:
    """Tests for aberration_in_longitude_latitude function."""

    @pytest.mark.unit
    def test_simplified_aberration(self):
        """Test the simplified Bradley aberration formula."""
        # Object at lon=0, lat=0
        # Sun at lon=90 (so Earth velocity is toward lon=0)
        lon, lat = 0.0, 0.0
        sun_lon = 90.0
        obliquity = 23.4

        d_lon, d_lat = aberration_in_longitude_latitude(lon, lat, sun_lon, obliquity)

        # Maximum aberration ~20 arcsec = 0.0057 degrees
        assert abs(d_lon) < 0.01  # degrees


# =============================================================================
# COMBINED LIGHT-TIME AND ABERRATION TESTS
# =============================================================================


class TestApplyLightTimeAndAberration:
    """Tests for apply_light_time_and_aberration function."""

    @pytest.mark.unit
    def test_combined_correction_without_target_func(self):
        """Test with simple light-time (no iteration)."""
        jd_tt = J2000
        geometric_pos = (5.0, 0.0, 0.0)
        earth_pos = (0.0, 1.0, 0.0)  # Not at origin
        earth_vel = (0.0, 0.0172, 0.0)

        apparent, lt = apply_light_time_and_aberration(
            jd_tt, geometric_pos, earth_pos, earth_vel, target_func=None
        )

        # Distance should be preserved (approximately)
        orig_dist = math.sqrt(sum(x * x for x in geometric_pos))
        app_dist = math.sqrt(sum(x * x for x in apparent))
        assert abs(orig_dist - app_dist) < 0.001

        # Light-time should be calculated from geometric position
        expected_lt = orig_dist / C_LIGHT_AU_DAY
        assert abs(lt - expected_lt) < 1e-10


class TestComputeEarthVelocity:
    """Tests for compute_earth_velocity function."""

    @pytest.mark.unit
    def test_earth_velocity_magnitude(self):
        """Earth's velocity should be approximately 0.0172 AU/day."""
        jd_tt = J2000
        vx, vy, vz = compute_earth_velocity(jd_tt)

        v_mag = math.sqrt(vx * vx + vy * vy + vz * vz)

        # Earth's mean orbital velocity: ~0.0172 AU/day
        assert 0.015 < v_mag < 0.020

    @pytest.mark.unit
    def test_earth_velocity_primarily_ecliptic(self):
        """Earth's velocity should be primarily in the ecliptic plane."""
        jd_tt = J2000
        vx, vy, vz = compute_earth_velocity(jd_tt)

        v_ecliptic = math.sqrt(vx * vx + vy * vy)
        v_total = math.sqrt(vx * vx + vy * vy + vz * vz)

        # z component should be very small compared to total
        assert abs(vz) < 0.001 * v_ecliptic  # < 0.1% of ecliptic velocity


# =============================================================================
# COORDINATE CONVERSION TESTS
# =============================================================================


class TestSphericalToCartesianVelocity:
    """Tests for spherical_to_cartesian_velocity function."""

    @pytest.mark.unit
    def test_zero_velocity(self):
        """Zero spherical velocity should give zero Cartesian velocity."""
        x, y, z, vx, vy, vz = spherical_to_cartesian_velocity(
            lon=45.0, lat=30.0, dist=2.0, dlon=0.0, dlat=0.0, ddist=0.0
        )

        # Position should be correct
        assert abs(z - 2.0 * math.sin(30 * DEG_TO_RAD)) < 1e-10

        # Velocity should be zero
        assert abs(vx) < 1e-10
        assert abs(vy) < 1e-10
        assert abs(vz) < 1e-10

    @pytest.mark.unit
    def test_radial_velocity(self):
        """Pure radial velocity (ddist only)."""
        lon, lat, dist = 0.0, 0.0, 1.0
        ddist = 1.0  # 1 AU/day radial velocity

        x, y, z, vx, vy, vz = spherical_to_cartesian_velocity(
            lon=lon, lat=lat, dist=dist, dlon=0.0, dlat=0.0, ddist=ddist
        )

        # At lon=0, lat=0, radial velocity is purely in +x
        assert abs(vx - ddist) < 1e-10
        assert abs(vy) < 1e-10
        assert abs(vz) < 1e-10


# =============================================================================
# INTEGRATION TESTS
# =============================================================================


class TestAberrationAccuracy:
    """Integration tests comparing aberration corrections with known values."""

    @pytest.mark.unit
    def test_aberration_order_of_magnitude(self):
        """Verify aberration correction is on the order of 20 arcsec."""
        # Typical Earth velocity
        earth_vel = (0.0, 0.0172, 0.0)

        # Object at various positions
        for lon in [0, 90, 180, 270]:
            x, y, z = spherical_to_cartesian(lon, 0.0, 5.0)
            position = (x, y, z)

            apparent = apply_aberration_to_position(position, earth_vel)

            # Convert back to spherical
            app_lon, app_lat, app_dist = cartesian_to_spherical(*apparent)
            orig_lon, orig_lat, orig_dist = cartesian_to_spherical(*position)

            # Angular difference should be ~20 arcsec or less
            d_lon = abs(app_lon - orig_lon)
            if d_lon > 180:
                d_lon = 360 - d_lon
            d_lon_arcsec = d_lon * 3600

            assert d_lon_arcsec < 25  # Max ~20.5 arcsec


class TestLightTimeAccuracy:
    """Integration tests for light-time correction accuracy."""

    @pytest.mark.unit
    def test_light_time_for_outer_planets(self):
        """Light-time for outer planets should be significant."""
        # Jupiter ~5 AU: ~42 minutes
        lt_jupiter = light_time_correction(5.0)
        lt_jupiter_min = lt_jupiter * 24 * 60
        assert 35 < lt_jupiter_min < 50

        # Saturn ~10 AU: ~83 minutes
        lt_saturn = light_time_correction(10.0)
        lt_saturn_min = lt_saturn * 24 * 60
        assert 75 < lt_saturn_min < 95

        # Neptune ~30 AU: ~4 hours
        lt_neptune = light_time_correction(30.0)
        lt_neptune_hours = lt_neptune * 24
        assert 3 < lt_neptune_hours < 5
