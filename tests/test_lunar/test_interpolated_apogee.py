"""
Tests for the interpolated lunar apogee (SE_INTP_APOG) and perigee (SE_INTP_PERG).

The interpolated apogee is a smoothed version of the osculating apogee that
removes spurious short-period oscillations by polynomial interpolation through
successive osculating apogee positions.

These tests verify:
1. Basic functionality and API compatibility
2. Return value format and ranges
3. Relationship between interpolated and osculating apogee
4. Velocity calculations
5. Perigee as opposite of apogee
"""

import math
import pytest
import libephemeris as swe


class TestInterpolatedApogeeBasic:
    """Basic functionality tests for interpolated apogee."""

    def test_calc_ut_returns_valid_position(self):
        """Test that swe_calc_ut returns a valid position for SE_INTP_APOG."""
        jd_ut = 2451545.0  # J2000.0
        result, retflag = swe.swe_calc_ut(jd_ut, swe.SE_INTP_APOG, swe.SEFLG_SPEED)

        # Result should be a tuple of 6 values
        assert len(result) == 6
        lon, lat, dist, speed_lon, speed_lat, speed_dist = result

        # Longitude should be in [0, 360)
        assert 0 <= lon < 360, f"Longitude {lon} out of range"

        # Latitude should be reasonable (typically small for lunar apogee)
        assert -10 < lat < 10, f"Latitude {lat} seems unreasonable"

        # Eccentricity (returned as dist) should be reasonable (~0.055)
        assert 0.04 < dist < 0.07, f"Eccentricity {dist} seems unreasonable"

    def test_calc_returns_valid_position(self):
        """Test that swe_calc (TT) returns a valid position for SE_INTP_APOG."""
        jd_tt = 2451545.0  # J2000.0
        result, retflag = swe.swe_calc(jd_tt, swe.SE_INTP_APOG, 0)

        assert len(result) == 6
        lon = result[0]
        assert 0 <= lon < 360

    def test_interpolated_perigee_basic(self):
        """Test that SE_INTP_PERG returns valid position."""
        jd_ut = 2451545.0
        result, retflag = swe.swe_calc_ut(jd_ut, swe.SE_INTP_PERG, 0)

        assert len(result) == 6
        lon = result[0]
        assert 0 <= lon < 360

    def test_perigee_opposite_apogee(self):
        """Test that interpolated perigee is 180 degrees from apogee."""
        jd_ut = 2451545.0

        apogee_result, _ = swe.swe_calc_ut(jd_ut, swe.SE_INTP_APOG, 0)
        perigee_result, _ = swe.swe_calc_ut(jd_ut, swe.SE_INTP_PERG, 0)

        apogee_lon = apogee_result[0]
        perigee_lon = perigee_result[0]

        # Calculate angular difference
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Should be very close to 180 degrees
        assert abs(diff - 180) < 0.1, f"Difference from 180: {abs(diff - 180)}"


class TestInterpolatedApogeeVelocity:
    """Tests for velocity/speed calculations."""

    def test_velocity_is_calculated(self):
        """Test that velocity is calculated when SEFLG_SPEED is set."""
        jd_ut = 2451545.0
        result, _ = swe.swe_calc_ut(jd_ut, swe.SE_INTP_APOG, swe.SEFLG_SPEED)

        speed_lon = result[3]
        # Velocity should be non-zero
        assert speed_lon != 0.0, "Velocity should be non-zero"

    def test_velocity_reasonable_magnitude(self):
        """Test that velocity has reasonable magnitude."""
        jd_ut = 2451545.0
        result, _ = swe.swe_calc_ut(jd_ut, swe.SE_INTP_APOG, swe.SEFLG_SPEED)

        speed_lon = result[3]
        # Mean apogee moves at about 40.7 degrees/year = 0.111 deg/day
        # Interpolated can show larger oscillations due to the polynomial fit
        # The osculating apogee varies by several degrees per day
        # Allow wide tolerance
        assert abs(speed_lon) < 10, f"Velocity {speed_lon} deg/day seems unreasonable"

    def test_no_velocity_without_flag(self):
        """Test that velocity is zero when SEFLG_SPEED is not set."""
        jd_ut = 2451545.0
        result, _ = swe.swe_calc_ut(jd_ut, swe.SE_INTP_APOG, 0)

        speed_lon = result[3]
        assert speed_lon == 0.0


class TestInterpolatedVsOsculating:
    """Tests comparing interpolated and osculating apogee."""

    def test_interpolated_smoother_than_osculating(self):
        """Test that interpolated apogee varies more smoothly than osculating."""
        # Sample both over several days and check that interpolated is smoother
        jd_start = 2451545.0
        dt = 1.0  # 1 day step

        oscu_positions = []
        intp_positions = []

        for i in range(10):
            jd = jd_start + i * dt

            oscu_result, _ = swe.swe_calc_ut(jd, swe.SE_OSCU_APOG, 0)
            intp_result, _ = swe.swe_calc_ut(jd, swe.SE_INTP_APOG, 0)

            oscu_positions.append(oscu_result[0])
            intp_positions.append(intp_result[0])

        # Calculate day-to-day changes
        oscu_changes = []
        intp_changes = []

        for i in range(1, len(oscu_positions)):
            oscu_diff = oscu_positions[i] - oscu_positions[i - 1]
            intp_diff = intp_positions[i] - intp_positions[i - 1]

            # Handle wrap-around
            if oscu_diff > 180:
                oscu_diff -= 360
            elif oscu_diff < -180:
                oscu_diff += 360
            if intp_diff > 180:
                intp_diff -= 360
            elif intp_diff < -180:
                intp_diff += 360

            oscu_changes.append(oscu_diff)
            intp_changes.append(intp_diff)

        # Calculate variance of changes (smoother = lower variance)
        oscu_mean = sum(oscu_changes) / len(oscu_changes)
        intp_mean = sum(intp_changes) / len(intp_changes)

        oscu_variance = sum((c - oscu_mean) ** 2 for c in oscu_changes) / len(
            oscu_changes
        )
        intp_variance = sum((c - intp_mean) ** 2 for c in intp_changes) / len(
            intp_changes
        )

        # Interpolated should have lower variance (smoother)
        assert intp_variance < oscu_variance, (
            f"Interpolated variance {intp_variance} should be less than osculating {oscu_variance}"
        )

    def test_interpolated_differs_from_osculating(self):
        """Test that interpolated and osculating apogee give different results."""
        jd_ut = 2451545.0

        oscu_result, _ = swe.swe_calc_ut(jd_ut, swe.SE_OSCU_APOG, 0)
        intp_result, _ = swe.swe_calc_ut(jd_ut, swe.SE_INTP_APOG, 0)

        oscu_lon = oscu_result[0]
        intp_lon = intp_result[0]

        # They should differ but not too much
        diff = abs(oscu_lon - intp_lon)
        if diff > 180:
            diff = 360 - diff

        # Difference is typically a few degrees due to short-period oscillations
        # being smoothed out
        assert diff > 0.01, "Interpolated and osculating should differ"
        assert diff < 30, f"Difference {diff} degrees seems too large"


class TestInterpolatedApogeeAtMultipleDates:
    """Tests at various dates to ensure stability."""

    @pytest.mark.parametrize(
        "jd,description",
        [
            (2451545.0, "J2000.0"),
            (2440587.5, "1970-01-01 (Unix epoch)"),
            (2458849.5, "2020-01-01"),
            (2460676.5, "2025-01-01"),
        ],
    )
    def test_valid_at_various_dates(self, jd, description):
        """Test that interpolated apogee works at various dates."""
        result, _ = swe.swe_calc_ut(jd, swe.SE_INTP_APOG, swe.SEFLG_SPEED)

        lon = result[0]
        speed = result[3]

        assert 0 <= lon < 360, f"Invalid longitude {lon} at {description}"
        # Speed can be positive or negative due to oscillations in the interpolation
        assert abs(speed) < 10, f"Speed {speed} seems unreasonable at {description}"


class TestInterpolatedApogeeDirectFunctions:
    """Test the direct lunar.py functions."""

    def test_calc_interpolated_apogee_directly(self):
        """Test calling calc_interpolated_apogee directly."""
        from libephemeris import lunar

        jd_tt = 2451545.0
        lon, lat, ecc = lunar.calc_interpolated_apogee(jd_tt)

        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.04 < ecc < 0.07

    def test_calc_interpolated_perigee_directly(self):
        """Test calling calc_interpolated_perigee directly."""
        from libephemeris import lunar

        jd_tt = 2451545.0
        lon, lat, ecc = lunar.calc_interpolated_perigee(jd_tt)

        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.04 < ecc < 0.07


class TestInterpolatedApogeeConsistency:
    """Test internal consistency of the implementation."""

    def test_continuous_motion(self):
        """Test that interpolated apogee moves continuously (no jumps)."""
        jd_start = 2451545.0
        dt = 0.1  # 2.4 hour steps

        prev_lon = None
        for i in range(50):
            jd = jd_start + i * dt
            result, _ = swe.swe_calc_ut(jd, swe.SE_INTP_APOG, 0)
            lon = result[0]

            if prev_lon is not None:
                diff = lon - prev_lon
                # Handle wrap-around
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

                # Maximum daily motion ~0.2 deg, so in 2.4 hours should be < 0.02 deg
                # Allow more tolerance since interpolation can cause small shifts
                assert abs(diff) < 0.5, (
                    f"Jump of {diff} degrees at step {i} seems too large"
                )

            prev_lon = lon

    def test_linear_solver_helper(self):
        """Test the internal linear system solver."""
        from libephemeris.lunar import _solve_linear_system

        # Simple 2x2 system: [[2, 1], [1, 3]] x = [4, 5]
        # From 2x + y = 4 and x + 3y = 5:
        # Solving: x = 1.4, y = 1.2
        # Verification: 2*1.4 + 1.2 = 4, 1.4 + 3*1.2 = 5
        A = [[2.0, 1.0], [1.0, 3.0]]
        b = [4.0, 5.0]

        x = _solve_linear_system(A, b)

        assert abs(x[0] - 1.4) < 1e-10
        assert abs(x[1] - 1.2) < 1e-10
