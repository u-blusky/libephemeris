"""
Tests for house cusp velocity calculations with SEFLG_SPEED flag.

Tests that swe_houses_ex2 and swe_houses_armc_ex2 respect the SEFLG_SPEED flag:
- When SEFLG_SPEED is not set, velocities should be zero
- When SEFLG_SPEED is set, velocities should be calculated
- Positions should be identical regardless of flag
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED, SEFLG_SIDEREAL, SE_SIDM_LAHIRI


class TestHousesEx2SpeedFlag:
    """Tests for SEFLG_SPEED flag in swe_houses_ex2."""

    @pytest.mark.unit
    def test_houses_ex2_without_speed_flag_returns_zero_velocities(self):
        """When SEFLG_SPEED is not set, velocities should be zero."""
        jd = 2451545.0  # J2000.0
        lat, lon = 41.9, 12.5

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd,
            lat,
            lon,
            ord("P"),
            0,  # No SEFLG_SPEED
        )

        # All cusp velocities should be zero
        for i, speed in enumerate(cusps_speed):
            assert speed == 0.0, f"Cusp {i + 1} speed should be 0, got {speed}"

        # All ascmc velocities should be zero
        for i, speed in enumerate(ascmc_speed):
            assert speed == 0.0, f"ASCMC {i} speed should be 0, got {speed}"

    @pytest.mark.unit
    def test_houses_ex2_with_speed_flag_returns_nonzero_velocities(self):
        """When SEFLG_SPEED is set, velocities should be calculated."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SPEED
        )

        # All cusp velocities should be non-zero
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"
            # Reasonable range: ~100-600 deg/day at mid-latitudes
            assert 100 < abs(speed) < 600, f"Cusp {i + 1} speed {speed} out of range"

        # ASC and MC velocities should be non-zero
        assert ascmc_speed[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed[1] != 0.0, "MC speed should be non-zero"
        # MC should be close to 360 deg/day
        assert 300 < abs(ascmc_speed[1]) < 400, (
            f"MC speed {ascmc_speed[1]} out of range"
        )

    @pytest.mark.unit
    def test_houses_ex2_positions_same_with_or_without_speed(self):
        """Cusp positions should be identical regardless of SEFLG_SPEED flag."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_no_speed, ascmc_no_speed, _, _ = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), 0
        )
        cusps_with_speed, ascmc_with_speed, _, _ = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SPEED
        )

        # Positions should be identical
        for i in range(12):
            assert cusps_no_speed[i] == cusps_with_speed[i], (
                f"Cusp {i + 1} position differs"
            )

        for i in range(8):
            assert ascmc_no_speed[i] == ascmc_with_speed[i], (
                f"ASCMC {i} position differs"
            )

    @pytest.mark.unit
    def test_houses_ex2_with_sidereal_and_speed_flag(self):
        """SEFLG_SPEED should work together with SEFLG_SIDEREAL."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        # Set Lahiri ayanamsa
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        # With both flags
        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SIDEREAL | SEFLG_SPEED
        )

        # Velocities should still be calculated
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"

        # With SEFLG_SIDEREAL but not SEFLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SIDEREAL
        )

        # Velocities should be zero
        for i, speed in enumerate(cusps_speed2):
            assert speed == 0.0, f"Cusp {i + 1} speed should be 0 without SEFLG_SPEED"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("R"), "Regiomontanus"),
            (ord("C"), "Campanus"),
            (ord("E"), "Equal"),
            (ord("B"), "Alcabitius"),
            (ord("M"), "Morinus"),
        ],
    )
    def test_houses_ex2_speed_flag_various_systems(self, hsys, name):
        """SEFLG_SPEED flag should work for all house systems."""
        jd = 2451545.0
        lat, lon = 45.0, 12.0

        # Without SEFLG_SPEED
        cusps1, ascmc1, cusps_speed1, ascmc_speed1 = ephem.swe_houses_ex2(
            jd, lat, lon, hsys, 0
        )

        # With SEFLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.swe_houses_ex2(
            jd, lat, lon, hsys, SEFLG_SPEED
        )

        # Positions should match
        assert cusps1 == cusps2, f"{name}: positions differ"
        assert ascmc1 == ascmc2, f"{name}: ascmc differ"

        # Without flag: zero velocities
        assert all(s == 0.0 for s in cusps_speed1), (
            f"{name}: speed should be 0 without flag"
        )

        # With flag: non-zero velocities
        assert any(s != 0.0 for s in cusps_speed2), (
            f"{name}: speed should be non-zero with flag"
        )

    @pytest.mark.unit
    def test_houses_ex2_speed_flag_whole_sign(self):
        """Whole Sign cusps are fixed at sign boundaries, so velocity is zero.

        Whole Sign houses place cusps at exact sign boundaries (0°, 30°, 60°, etc.).
        Within a short time interval, these don't change, so cusp velocities are zero.
        However, ASC and MC velocities should still be calculated.
        """
        jd = 2451545.0
        lat, lon = 45.0, 12.0

        # Without SEFLG_SPEED
        cusps1, ascmc1, cusps_speed1, ascmc_speed1 = ephem.swe_houses_ex2(
            jd, lat, lon, ord("W"), 0
        )

        # With SEFLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.swe_houses_ex2(
            jd, lat, lon, ord("W"), SEFLG_SPEED
        )

        # Positions should match
        assert cusps1 == cusps2, "Whole Sign: positions differ"
        assert ascmc1 == ascmc2, "Whole Sign: ascmc differ"

        # Without flag: zero velocities
        assert all(s == 0.0 for s in cusps_speed1)
        assert all(s == 0.0 for s in ascmc_speed1)

        # With flag: cusp velocities are zero (fixed at sign boundaries)
        # but ASC/MC velocities should be calculated
        assert all(s == 0.0 for s in cusps_speed2), (
            "Whole Sign cusps are fixed, velocity should be 0"
        )
        # ASC and MC still have velocities
        assert ascmc_speed2[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed2[1] != 0.0, "MC speed should be non-zero"


class TestHousesArmcEx2SpeedFlag:
    """Tests for SEFLG_SPEED flag in swe_houses_armc_ex2."""

    @pytest.mark.unit
    def test_houses_armc_ex2_without_speed_flag_returns_zero_velocities(self):
        """When SEFLG_SPEED is not set, velocities should be zero."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_armc_ex2(
            armc,
            lat,
            eps,
            ord("P"),
            0,  # No SEFLG_SPEED
        )

        # All cusp velocities should be zero
        for i, speed in enumerate(cusps_speed):
            assert speed == 0.0, f"Cusp {i + 1} speed should be 0, got {speed}"

        # All ascmc velocities should be zero
        for i, speed in enumerate(ascmc_speed):
            assert speed == 0.0, f"ASCMC {i} speed should be 0, got {speed}"

    @pytest.mark.unit
    def test_houses_armc_ex2_with_speed_flag_returns_nonzero_velocities(self):
        """When SEFLG_SPEED is set, velocities should be calculated."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_armc_ex2(
            armc, lat, eps, ord("P"), SEFLG_SPEED
        )

        # All cusp velocities should be non-zero
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"
            # Reasonable range: ~100-600 deg/day at mid-latitudes
            assert 100 < abs(speed) < 600, f"Cusp {i + 1} speed {speed} out of range"

        # ASC and MC velocities should be non-zero
        assert ascmc_speed[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed[1] != 0.0, "MC speed should be non-zero"
        # ARMC velocity should be exactly 360 deg/day
        assert abs(ascmc_speed[2] - 360.0) < 0.01, (
            f"ARMC speed {ascmc_speed[2]} should be ~360"
        )

    @pytest.mark.unit
    def test_houses_armc_ex2_positions_same_with_or_without_speed(self):
        """Cusp positions should be identical regardless of SEFLG_SPEED flag."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps_no_speed, ascmc_no_speed, _, _ = ephem.swe_houses_armc_ex2(
            armc, lat, eps, ord("P"), 0
        )
        cusps_with_speed, ascmc_with_speed, _, _ = ephem.swe_houses_armc_ex2(
            armc, lat, eps, ord("P"), SEFLG_SPEED
        )

        # Positions should be identical
        for i in range(12):
            assert cusps_no_speed[i] == cusps_with_speed[i], (
                f"Cusp {i + 1} position differs"
            )

        for i in range(8):
            assert ascmc_no_speed[i] == ascmc_with_speed[i], (
                f"ASCMC {i} position differs"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("R"), "Regiomontanus"),
            (ord("C"), "Campanus"),
            (ord("E"), "Equal"),
            (ord("B"), "Alcabitius"),
            (ord("M"), "Morinus"),
        ],
    )
    def test_houses_armc_ex2_speed_flag_various_systems(self, hsys, name):
        """SEFLG_SPEED flag should work for all house systems."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        # Without SEFLG_SPEED
        cusps1, ascmc1, cusps_speed1, ascmc_speed1 = ephem.swe_houses_armc_ex2(
            armc, lat, eps, hsys, 0
        )

        # With SEFLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.swe_houses_armc_ex2(
            armc, lat, eps, hsys, SEFLG_SPEED
        )

        # Positions should match
        assert cusps1 == cusps2, f"{name}: positions differ"
        assert ascmc1 == ascmc2, f"{name}: ascmc differ"

        # Without flag: zero velocities
        assert all(s == 0.0 for s in cusps_speed1), (
            f"{name}: speed should be 0 without flag"
        )

        # With flag: non-zero velocities
        assert any(s != 0.0 for s in cusps_speed2), (
            f"{name}: speed should be non-zero with flag"
        )

    @pytest.mark.unit
    def test_houses_armc_ex2_speed_flag_whole_sign(self):
        """Whole Sign cusps are fixed at sign boundaries, so velocity is zero.

        Whole Sign houses place cusps at exact sign boundaries (0°, 30°, 60°, etc.).
        Within a short time interval, these don't change, so cusp velocities are zero.
        However, ASC and MC velocities should still be calculated.
        """
        armc = 150.0
        lat = 45.0
        eps = 23.44

        # Without SEFLG_SPEED
        cusps1, ascmc1, cusps_speed1, ascmc_speed1 = ephem.swe_houses_armc_ex2(
            armc, lat, eps, ord("W"), 0
        )

        # With SEFLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.swe_houses_armc_ex2(
            armc, lat, eps, ord("W"), SEFLG_SPEED
        )

        # Positions should match
        assert cusps1 == cusps2, "Whole Sign: positions differ"
        assert ascmc1 == ascmc2, "Whole Sign: ascmc differ"

        # Without flag: zero velocities
        assert all(s == 0.0 for s in cusps_speed1)
        assert all(s == 0.0 for s in ascmc_speed1)

        # With flag: cusp velocities are zero (fixed at sign boundaries)
        # but ASC/MC velocities should be calculated
        assert all(s == 0.0 for s in cusps_speed2), (
            "Whole Sign cusps are fixed, velocity should be 0"
        )
        # ASC and MC still have velocities
        assert ascmc_speed2[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed2[1] != 0.0, "MC speed should be non-zero"


class TestHouseCuspVelocityEdgeCases:
    """Edge case tests for house cusp velocity calculations."""

    @pytest.mark.edge_case
    def test_houses_ex2_equator_with_speed(self):
        """SEFLG_SPEED should work at the equator."""
        jd = 2451545.0
        lat, lon = 0.0, 0.0

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)
        assert ascmc_speed[1] != 0.0  # MC velocity

    @pytest.mark.edge_case
    def test_houses_ex2_high_latitude_with_speed(self):
        """SEFLG_SPEED should work at high latitudes."""
        jd = 2451545.0
        lat, lon = 60.0, 0.0

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)

    @pytest.mark.edge_case
    def test_houses_ex2_southern_hemisphere_with_speed(self):
        """SEFLG_SPEED should work in the southern hemisphere."""
        jd = 2451545.0
        lat, lon = -33.9, 151.2  # Sydney

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)

    @pytest.mark.edge_case
    def test_houses_armc_ex2_equator_with_speed(self):
        """SEFLG_SPEED should work at the equator for ARMC-based calculation."""
        armc = 180.0
        lat = 0.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_armc_ex2(
            armc, lat, eps, ord("P"), SEFLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)


class TestHouseCuspVelocityProgression:
    """Tests for progressed chart applications."""

    @pytest.mark.unit
    def test_velocity_for_progressed_mc_crossing(self):
        """Velocities should allow calculating when MC crosses a point."""
        jd = 2451545.0
        lat, lon = 45.0, 0.0

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SPEED
        )

        # MC velocity should be positive (moving forward)
        mc_speed = ascmc_speed[1]
        assert mc_speed > 0, f"MC speed {mc_speed} should be positive"

        # Calculate time for MC to move 1 degree
        # time = distance / velocity (in days)
        time_for_1_deg = 1.0 / mc_speed

        # Should be approximately 4 minutes (1/360 of a day)
        expected_time = 1.0 / 360.0  # ~0.00278 days
        assert abs(time_for_1_deg - expected_time) < 0.001, (
            f"Time for 1 deg: {time_for_1_deg}, expected ~{expected_time}"
        )

    @pytest.mark.unit
    def test_velocities_consistent_across_time(self):
        """Velocities should be consistent when calculated at different times."""
        lat, lon = 45.0, 0.0
        jd1 = 2451545.0
        jd2 = 2451545.5  # 12 hours later

        _, _, _, ascmc_speed1 = ephem.swe_houses_ex2(
            jd1, lat, lon, ord("P"), SEFLG_SPEED
        )
        _, _, _, ascmc_speed2 = ephem.swe_houses_ex2(
            jd2, lat, lon, ord("P"), SEFLG_SPEED
        )

        # MC velocity should be similar at different times
        # (variation possible due to Earth's orbital motion)
        mc_speed1 = ascmc_speed1[1]
        mc_speed2 = ascmc_speed2[1]

        # Should be within 1% of each other
        rel_diff = abs(mc_speed1 - mc_speed2) / mc_speed1
        assert rel_diff < 0.01, f"MC speed varies too much: {mc_speed1} vs {mc_speed2}"
