"""
Minor Bodies (Asteroids, TNOs) Comparison Tests.

Compares minor body calculations between pyswisseph and libephemeris.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_CHIRON,
    SE_PHOLUS,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

MAIN_ASTEROID_TOL = 0.01  # degrees
CENTAUR_TOL = 0.01  # degrees
TNO_TOL = 0.1  # degrees (relaxed for TNOs)


# ============================================================================
# TEST DATA
# ============================================================================

MAIN_ASTEROIDS = [
    (SE_CERES, "Ceres"),
    (SE_PALLAS, "Pallas"),
    (SE_JUNO, "Juno"),
    (SE_VESTA, "Vesta"),
]

CENTAURS = [
    (SE_CHIRON, "Chiron"),
    (SE_PHOLUS, "Pholus"),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Current"),
    (1980, 5, 20, 14.5, "Past"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMainAsteroids:
    """Compare main belt asteroid calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", MAIN_ASTEROIDS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_main_asteroid_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test main asteroid positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MAIN_ASTEROID_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestCentaurs:
    """Compare centaur calculations (Chiron, Pholus)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", CENTAURS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_centaur_position(self, body_id, body_name, year, month, day, hour, desc):
        """Test centaur positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < CENTAUR_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestChironSpecific:
    """Specific tests for Chiron (most commonly used asteroid)."""

    @pytest.mark.comparison
    def test_chiron_at_j2000(self):
        """Test Chiron position at J2000 epoch."""
        jd = 2451545.0

        pos_swe, _ = swe.calc_ut(jd, swe.CHIRON, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_CHIRON, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < 0.001, f"Chiron at J2000 diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_chiron_with_speed(self):
        """Test Chiron position with velocity."""
        jd = swe.julday(2024, 1, 1, 12.0)

        pos_swe, _ = swe.calc_ut(jd, swe.CHIRON, swe.FLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_CHIRON, 256)  # SEFLG_SPEED

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < CENTAUR_TOL, f"Chiron longitude diff {diff_lon:.6f}°"
        assert diff_speed < 0.01, f"Chiron speed diff {diff_speed:.6f}°/day"


class TestCeresSpecific:
    """Specific tests for Ceres (dwarf planet)."""

    @pytest.mark.comparison
    def test_ceres_at_j2000(self):
        """Test Ceres position at J2000 epoch."""
        jd = 2451545.0

        pos_swe, _ = swe.calc_ut(jd, swe.CERES, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_CERES, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < 0.001, f"Ceres at J2000 diff {diff:.6f}° exceeds tight tolerance"


class TestAsteroidVelocity:
    """Test asteroid velocity calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", MAIN_ASTEROIDS + CENTAURS)
    def test_asteroid_velocity(self, body_id, body_name):
        """Test asteroid velocity calculations match."""
        jd = swe.julday(2024, 6, 15, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 256)  # SEFLG_SPEED
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_speed < 0.01, (
            f"{body_name} velocity diff {diff_speed:.6f}°/day exceeds tolerance"
        )
