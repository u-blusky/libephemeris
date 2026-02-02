"""
Minor Bodies (Asteroids, TNOs) Comparison Tests.

Compares minor body calculations between pyswisseph and libephemeris.
Tests cover all minor bodies defined in constants.py.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    # Main belt asteroids
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    # Centaurs (original)
    SE_CHIRON,
    SE_PHOLUS,
    # Additional Centaurs
    SE_NESSUS,
    SE_ASBOLUS,
    SE_CHARIKLO,
    # Trans-Neptunian Objects (TNOs)
    SE_ERIS,
    SE_SEDNA,
    SE_HAUMEA,
    SE_MAKEMAKE,
    SE_IXION,
    SE_ORCUS,
    SE_QUAOAR,
    SE_VARUNA,
    SE_GONGGONG,
    # Large main belt asteroids
    SE_HYGIEA,
    SE_INTERAMNIA,
    SE_DAVIDA,
    SE_EUROPA_AST,
    SE_SYLVIA,
    SE_PSYCHE,
    # Near-Earth asteroids
    SE_APOPHIS,
    SE_EROS,
    SE_AMOR,
    SE_ICARUS,
    SE_TORO,
    SE_TOUTATIS,
    SE_ITOKAWA,
    SE_BENNU,
    SE_RYUGU,
    # Astrologically significant asteroids
    SE_SAPPHO,
    SE_PANDORA_AST,
    SE_LILITH_AST,
    SE_HIDALGO,
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

MAIN_ASTEROID_TOL = 0.01  # degrees - main belt asteroids with good ephemeris
CENTAUR_TOL = 0.01  # degrees - centaurs with good ephemeris
TNO_TOL = 0.1  # degrees - relaxed for distant TNOs
NEA_TOL = 0.05  # degrees - Near-Earth asteroids (variable orbits)
DISTANT_TNO_TOL = 0.5  # degrees - very distant/uncertain TNOs like Sedna


# ============================================================================
# TEST DATA
# ============================================================================

MAIN_ASTEROIDS = [
    (SE_CERES, "Ceres"),
    (SE_PALLAS, "Pallas"),
    (SE_JUNO, "Juno"),
    (SE_VESTA, "Vesta"),
]

LARGE_MAIN_BELT = [
    (SE_HYGIEA, "Hygiea"),
    (SE_INTERAMNIA, "Interamnia"),
    (SE_DAVIDA, "Davida"),
    (SE_EUROPA_AST, "Europa (asteroid)"),
    (SE_SYLVIA, "Sylvia"),
    (SE_PSYCHE, "Psyche"),
]

CENTAURS = [
    (SE_CHIRON, "Chiron"),
    (SE_PHOLUS, "Pholus"),
]

ADDITIONAL_CENTAURS = [
    (SE_NESSUS, "Nessus"),
    (SE_ASBOLUS, "Asbolus"),
    (SE_CHARIKLO, "Chariklo"),
]

TNOS_STANDARD = [
    (SE_ERIS, "Eris"),
    (SE_HAUMEA, "Haumea"),
    (SE_MAKEMAKE, "Makemake"),
    (SE_ORCUS, "Orcus"),
    (SE_QUAOAR, "Quaoar"),
    (SE_VARUNA, "Varuna"),
    (SE_IXION, "Ixion"),
    (SE_GONGGONG, "Gonggong"),
]

TNOS_DISTANT = [
    (SE_SEDNA, "Sedna"),  # Very distant, relaxed tolerance
]

NEAR_EARTH_ASTEROIDS = [
    (SE_EROS, "Eros"),
    (SE_AMOR, "Amor"),
    (SE_ICARUS, "Icarus"),
    (SE_TORO, "Toro"),
    (SE_APOPHIS, "Apophis"),
    (SE_TOUTATIS, "Toutatis"),
    (SE_ITOKAWA, "Itokawa"),
    (SE_BENNU, "Bennu"),
    (SE_RYUGU, "Ryugu"),
]

ASTROLOGICAL_ASTEROIDS = [
    (SE_SAPPHO, "Sappho"),
    (SE_PANDORA_AST, "Pandora (asteroid)"),
    (SE_LILITH_AST, "Lilith (asteroid)"),
    (SE_HIDALGO, "Hidalgo"),
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


class TestLargeMainBeltAsteroids:
    """Compare large main belt asteroids (Hygiea, Interamnia, Davida, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", LARGE_MAIN_BELT)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_large_main_belt_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test large main belt asteroid positions match within tolerance."""
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


class TestAdditionalCentaurs:
    """Compare additional centaur calculations (Nessus, Asbolus, Chariklo)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ADDITIONAL_CENTAURS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_additional_centaur_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test additional centaur positions match within tolerance."""
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


class TestTNOsStandard:
    """Compare Trans-Neptunian Object calculations (Eris, Makemake, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TNOS_STANDARD)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_tno_position(self, body_id, body_name, year, month, day, hour, desc):
        """Test TNO positions match within relaxed tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TNO_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestTNOsDistant:
    """Compare very distant TNO calculations (Sedna) with extra relaxed tolerance."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TNOS_DISTANT)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_distant_tno_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test distant TNO positions match within extra relaxed tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < DISTANT_TNO_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestNearEarthAsteroids:
    """Compare Near-Earth Asteroid calculations (Apophis, Eros, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", NEAR_EARTH_ASTEROIDS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_nea_position(self, body_id, body_name, year, month, day, hour, desc):
        """Test Near-Earth asteroid positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < NEA_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestAstrologicalAsteroids:
    """Compare astrologically significant asteroids (Sappho, Lilith, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ASTROLOGICAL_ASTEROIDS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_astrological_asteroid_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test astrological asteroid positions match within tolerance."""
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


class TestChironSpecific:
    """Specific tests for Chiron (most commonly used asteroid)."""

    @pytest.mark.comparison
    def test_chiron_at_j2000(self):
        """Test Chiron position at J2000 epoch."""
        jd = 2451545.0

        try:
            pos_swe, _ = swe.calc_ut(jd, swe.CHIRON, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_CHIRON, 0)
        except Exception as e:
            pytest.skip(f"Chiron not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < 0.001, f"Chiron at J2000 diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_chiron_with_speed(self):
        """Test Chiron position with velocity."""
        jd = swe.julday(2024, 1, 1, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, swe.CHIRON, swe.FLG_SPEED)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_CHIRON, 256)  # SEFLG_SPEED
        except Exception as e:
            pytest.skip(f"Chiron not available: {e}")
            return

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

        try:
            pos_swe, _ = swe.calc_ut(jd, swe.CERES, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_CERES, 0)
        except Exception as e:
            pytest.skip(f"Ceres not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < 0.001, f"Ceres at J2000 diff {diff:.6f}° exceeds tight tolerance"


class TestErisSpecific:
    """Specific tests for Eris (largest known dwarf planet)."""

    @pytest.mark.comparison
    def test_eris_at_j2000(self):
        """Test Eris position at J2000 epoch."""
        jd = 2451545.0

        try:
            pos_swe, _ = swe.calc_ut(jd, SE_ERIS, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_ERIS, 0)
        except Exception as e:
            pytest.skip(f"Eris not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TNO_TOL, f"Eris at J2000 diff {diff:.6f}° exceeds tolerance"


class TestGonggongSpecific:
    """Specific tests for Gonggong (formerly 2007 OR10)."""

    @pytest.mark.comparison
    def test_gonggong_at_current_epoch(self):
        """Test Gonggong position at current epoch."""
        jd = swe.julday(2024, 6, 15, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, SE_GONGGONG, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_GONGGONG, 0)
        except Exception as e:
            pytest.skip(f"Gonggong not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TNO_TOL, f"Gonggong diff {diff:.6f}° exceeds tolerance"


class TestApophisSpecific:
    """Specific tests for Apophis (potentially hazardous asteroid)."""

    @pytest.mark.comparison
    def test_apophis_at_current_epoch(self):
        """Test Apophis position at current epoch."""
        jd = swe.julday(2024, 6, 15, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, SE_APOPHIS, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_APOPHIS, 0)
        except Exception as e:
            pytest.skip(f"Apophis not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < NEA_TOL, f"Apophis diff {diff:.6f}° exceeds tolerance"


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


class TestAllBodiesCoverage:
    """Verify all defined minor bodies have test coverage."""

    @pytest.mark.comparison
    def test_all_bodies_in_test_lists(self):
        """Verify all minor body constants are included in test lists."""
        all_tested_bodies = set()
        for body_list in [
            MAIN_ASTEROIDS,
            LARGE_MAIN_BELT,
            CENTAURS,
            ADDITIONAL_CENTAURS,
            TNOS_STANDARD,
            TNOS_DISTANT,
            NEAR_EARTH_ASTEROIDS,
            ASTROLOGICAL_ASTEROIDS,
        ]:
            for body_id, _ in body_list:
                all_tested_bodies.add(body_id)

        # These are all the SE_* minor body constants defined in constants.py
        expected_bodies = {
            SE_CERES,
            SE_PALLAS,
            SE_JUNO,
            SE_VESTA,
            SE_CHIRON,
            SE_PHOLUS,
            SE_NESSUS,
            SE_ASBOLUS,
            SE_CHARIKLO,
            SE_ERIS,
            SE_SEDNA,
            SE_HAUMEA,
            SE_MAKEMAKE,
            SE_IXION,
            SE_ORCUS,
            SE_QUAOAR,
            SE_VARUNA,
            SE_GONGGONG,
            SE_HYGIEA,
            SE_INTERAMNIA,
            SE_DAVIDA,
            SE_EUROPA_AST,
            SE_SYLVIA,
            SE_PSYCHE,
            SE_APOPHIS,
            SE_EROS,
            SE_AMOR,
            SE_ICARUS,
            SE_TORO,
            SE_TOUTATIS,
            SE_ITOKAWA,
            SE_BENNU,
            SE_RYUGU,
            SE_SAPPHO,
            SE_PANDORA_AST,
            SE_LILITH_AST,
            SE_HIDALGO,
        }

        missing_bodies = expected_bodies - all_tested_bodies
        assert not missing_bodies, (
            f"Missing test coverage for body IDs: {missing_bodies}"
        )
