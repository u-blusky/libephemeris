"""
Sidereal/Ayanamsha Comparison Tests.

Compares all ayanamsha (sidereal) modes between pyswisseph and libephemeris.
Tests all 43 ayanamsha systems across different dates and planets.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_SIDM_USER,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
    SE_SIDM_DELUCE,
    SE_SIDM_RAMAN,
    SE_SIDM_USHASHASHI,
    SE_SIDM_KRISHNAMURTI,
    SE_SIDM_DJWHAL_KHUL,
    SE_SIDM_YUKTESHWAR,
    SE_SIDM_JN_BHASIN,
    SE_SIDM_BABYL_KUGLER1,
    SE_SIDM_BABYL_KUGLER2,
    SE_SIDM_BABYL_KUGLER3,
    SE_SIDM_BABYL_HUBER,
    SE_SIDM_BABYL_ETPSC,
    SE_SIDM_ALDEBARAN_15TAU,
    SE_SIDM_HIPPARCHOS,
    SE_SIDM_SASSANIAN,
    SE_SIDM_GALCENT_0SAG,
    SE_SIDM_J2000,
    SE_SIDM_J1900,
    SE_SIDM_B1950,
    SE_SIDM_SURYASIDDHANTA,
    SE_SIDM_SURYASIDDHANTA_MSUN,
    SE_SIDM_ARYABHATA,
    SE_SIDM_ARYABHATA_MSUN,
    SE_SIDM_SS_REVATI,
    SE_SIDM_SS_CITRA,
    SE_SIDM_TRUE_CITRA,
    SE_SIDM_TRUE_REVATI,
    SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_GALCENT_RGILBRAND,
    SE_SIDM_GALEQU_IAU1958,
    SE_SIDM_GALEQU_TRUE,
    SE_SIDM_GALEQU_MULA,
    SE_SIDM_GALALIGN_MARDYKS,
    SE_SIDM_TRUE_MULA,
    SE_SIDM_GALCENT_MULA_WILHELM,
    SE_SIDM_ARYABHATA_522,
    SE_SIDM_BABYL_BRITTON,
    SE_SIDM_TRUE_SHEORAN,
    SE_SIDM_GALCENT_COCHRANE,
    SE_SIDM_GALEQU_FIORENZA,
    SE_SIDM_VALENS_MOON,
    SEFLG_SIDEREAL,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# AYANAMSHA MODES
# ============================================================================

FORMULA_BASED_AYANAMSHAS = [
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_DELUCE, "De Luce"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_USHASHASHI, "Ushashashi"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
    (SE_SIDM_JN_BHASIN, "JN Bhasin"),
    (SE_SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
    (SE_SIDM_BABYL_KUGLER2, "Babylonian Kugler 2"),
    (SE_SIDM_BABYL_KUGLER3, "Babylonian Kugler 3"),
    (SE_SIDM_BABYL_HUBER, "Babylonian Huber"),
    (SE_SIDM_BABYL_ETPSC, "Babylonian ETPSC"),
    (SE_SIDM_ALDEBARAN_15TAU, "Aldebaran 15 Tau"),
    (SE_SIDM_HIPPARCHOS, "Hipparchos"),
    (SE_SIDM_SASSANIAN, "Sassanian"),
    (SE_SIDM_J2000, "J2000"),
    (SE_SIDM_J1900, "J1900"),
    (SE_SIDM_B1950, "B1950"),
    (SE_SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (SE_SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta Mean Sun"),
    (SE_SIDM_ARYABHATA, "Aryabhata"),
    (SE_SIDM_ARYABHATA_MSUN, "Aryabhata Mean Sun"),
    (SE_SIDM_SS_REVATI, "SS Revati"),
    (SE_SIDM_SS_CITRA, "SS Citra"),
    (SE_SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SE_SIDM_BABYL_BRITTON, "Babylonian Britton"),
]

STAR_BASED_AYANAMSHAS = [
    (SE_SIDM_TRUE_CITRA, "True Citra"),
    (SE_SIDM_TRUE_REVATI, "True Revati"),
    (SE_SIDM_TRUE_PUSHYA, "True Pushya"),
    (SE_SIDM_TRUE_MULA, "True Mula"),
    (SE_SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SE_SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    (SE_SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
    (SE_SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
    (SE_SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
    (SE_SIDM_GALEQU_IAU1958, "Galactic Equator IAU 1958"),
    (SE_SIDM_GALEQU_TRUE, "Galactic Equator True"),
    (SE_SIDM_GALEQU_MULA, "Galactic Equator Mula"),
    (SE_SIDM_GALEQU_FIORENZA, "Galactic Equator Fiorenza"),
    (SE_SIDM_GALALIGN_MARDYKS, "Galactic Alignment Mardyks"),
    (SE_SIDM_VALENS_MOON, "Valens Moon"),
]

ALL_AYANAMSHAS = FORMULA_BASED_AYANAMSHAS + STAR_BASED_AYANAMSHAS

TEST_DATES = [
    (2451545.0, "J2000.0"),
    (2415020.5, "1900-01-01"),
    (2433282.5, "1950-01-01"),
    (2469807.5, "2050-01-01"),
]

TEST_PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]

# Tolerances
STRICT_TOLERANCE = 0.001  # 3.6 arcseconds - for formula-based
RELAXED_TOLERANCE = 0.1  # 6 arcminutes - for star-based


def get_tolerance(sid_mode: int) -> float:
    """Get appropriate tolerance for ayanamsha mode."""
    star_modes = {m[0] for m in STAR_BASED_AYANAMSHAS}
    if sid_mode in star_modes:
        return RELAXED_TOLERANCE
    return STRICT_TOLERANCE


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestAyanamshaValues:
    """Compare ayanamsha values between implementations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_formula_based_ayanamsha(self, sid_mode, sid_name, jd, date_desc):
        """Test formula-based ayanamsha values match within strict tolerance."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = angular_diff(ayan_swe, ayan_py)

        assert diff < STRICT_TOLERANCE, (
            f"{sid_name} at {date_desc}: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", STAR_BASED_AYANAMSHAS)
    @pytest.mark.parametrize(
        "jd,date_desc", TEST_DATES[:2]
    )  # Fewer dates for star-based
    def test_star_based_ayanamsha(self, sid_mode, sid_name, jd, date_desc):
        """Test star-based ayanamsha values match within relaxed tolerance."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = angular_diff(ayan_swe, ayan_py)

        assert diff < RELAXED_TOLERANCE, (
            f"{sid_name} at {date_desc}: diff {diff:.6f}° exceeds tolerance"
        )


class TestSiderealPlanets:
    """Compare sidereal planetary positions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS[:10])
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_sidereal_planet_positions(
        self, sid_mode, sid_name, planet_id, planet_name
    ):
        """Test sidereal planetary positions match."""
        jd = 2451545.0  # J2000

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        pos_swe, _ = swe.calc_ut(jd, planet_id, SEFLG_SIDEREAL)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, SEFLG_SIDEREAL)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        tolerance = get_tolerance(sid_mode)

        assert diff_lon < tolerance, (
            f"{planet_name} sidereal ({sid_name}): diff {diff_lon:.6f}° exceeds tolerance"
        )


class TestMajorAyanamshas:
    """Detailed tests for most commonly used ayanamshas."""

    @pytest.mark.comparison
    def test_lahiri_at_j2000(self):
        """Test Lahiri ayanamsha at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SE_SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # Lahiri should be approximately 23.86 degrees at J2000
        assert 23.5 < ayan_py < 24.5, f"Lahiri value {ayan_py}° out of expected range"
        assert diff < 0.001, f"Lahiri diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_fagan_bradley_at_j2000(self):
        """Test Fagan-Bradley ayanamsha at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SE_SIDM_FAGAN_BRADLEY)
        ephem.swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # Fagan-Bradley should be approximately 24.74 degrees at J2000
        assert 24.0 < ayan_py < 25.5, (
            f"Fagan-Bradley value {ayan_py}° out of expected range"
        )
        assert diff < 0.001, f"Fagan-Bradley diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_raman_at_j2000(self):
        """Test Raman ayanamsha at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SE_SIDM_RAMAN)
        ephem.swe_set_sid_mode(SE_SIDM_RAMAN)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # Raman should be approximately 22.4 degrees at J2000
        assert 22.0 < ayan_py < 23.0, f"Raman value {ayan_py}° out of expected range"
        assert diff < 0.001, f"Raman diff {diff:.6f}° exceeds tight tolerance"


class TestAyanamshaProgression:
    """Test that ayanamsha changes correctly over time."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS[:5])
    def test_ayanamsha_increases_over_century(self, sid_mode, sid_name):
        """Test that ayanamsha increases over a century (precession)."""
        jd_1900 = 2415020.5
        jd_2000 = 2451545.0

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        # Get ayanamsha at both epochs from both libraries
        ayan_1900_swe = swe.get_ayanamsa_ut(jd_1900)
        ayan_2000_swe = swe.get_ayanamsa_ut(jd_2000)
        ayan_1900_py = ephem.swe_get_ayanamsa_ut(jd_1900)
        ayan_2000_py = ephem.swe_get_ayanamsa_ut(jd_2000)

        # Change should be positive (precession increases ayanamsha)
        change_swe = ayan_2000_swe - ayan_1900_swe
        change_py = ayan_2000_py - ayan_1900_py

        # Skip J2000 mode which is 0 at J2000 by definition
        if sid_mode != SE_SIDM_J2000:
            assert change_swe > 0.5, f"{sid_name}: expected precession over century"
            assert change_py > 0.5, f"{sid_name}: expected precession over century"

        # Changes should be similar
        change_diff = abs(change_swe - change_py)
        assert change_diff < 0.01, (
            f"{sid_name}: century change differs by {change_diff:.6f}°"
        )


class TestJ2000Mode:
    """Test J2000 ayanamsha mode specifically."""

    @pytest.mark.comparison
    def test_j2000_is_zero_at_epoch(self):
        """J2000 ayanamsha should be exactly 0 at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SE_SIDM_J2000)
        ephem.swe_set_sid_mode(SE_SIDM_J2000)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        assert abs(ayan_swe) < 0.0001, f"SWE J2000 at J2000 = {ayan_swe}, expected 0"
        assert abs(ayan_py) < 0.0001, f"PY J2000 at J2000 = {ayan_py}, expected 0"

    @pytest.mark.comparison
    def test_j2000_changes_over_time(self):
        """J2000 ayanamsha should change from the J2000 epoch."""
        swe.set_sid_mode(SE_SIDM_J2000)
        ephem.swe_set_sid_mode(SE_SIDM_J2000)

        ayan_1900 = ephem.swe_get_ayanamsa_ut(2415020.5)
        ayan_2000 = ephem.swe_get_ayanamsa_ut(2451545.0)
        ayan_2100 = ephem.swe_get_ayanamsa_ut(2488069.5)

        ayan_1900_swe = swe.get_ayanamsa_ut(2415020.5)
        ayan_2100_swe = swe.get_ayanamsa_ut(2488069.5)

        assert abs(ayan_1900 - ayan_1900_swe) < 0.001, (
            f"J2000 at 1900: lib={ayan_1900}, swe={ayan_1900_swe}"
        )
        assert abs(ayan_2000) < 0.001, f"J2000 at 2000 should be ~0, got {ayan_2000}"
        assert abs(ayan_2100 - ayan_2100_swe) < 0.001, (
            f"J2000 at 2100: lib={ayan_2100}, swe={ayan_2100_swe}"
        )


class TestUserDefinedAyanamsha:
    """
    Test SE_SIDM_USER (255) custom ayanamsha mode.

    This mode allows users to define custom ayanamsha with t0 (reference time)
    and ayan_t0 (ayanamsha value at reference time) parameters.
    Many Vedic astrologers use custom ayanamsha values.
    """

    @pytest.mark.comparison
    def test_user_ayanamsha_value_at_j2000(self):
        """Test SE_SIDM_USER ayanamsha value matches pyswisseph at J2000."""
        jd = 2451545.0  # J2000
        t0 = jd
        ayan_t0 = 23.5  # Custom ayanamsha value

        # Set user-defined mode with same parameters
        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # At t0, both should return exactly ayan_t0
        assert abs(ayan_py - ayan_t0) < 0.001, f"Expected {ayan_t0}, got {ayan_py}"
        assert diff < STRICT_TOLERANCE, (
            f"SE_SIDM_USER at J2000: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_user_ayanamsha_at_multiple_dates(self):
        """Test SE_SIDM_USER ayanamsha matches pyswisseph across multiple dates."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 23.5

        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        for jd, date_desc in TEST_DATES:
            ayan_swe = swe.get_ayanamsa_ut(jd)
            ayan_py = ephem.swe_get_ayanamsa_ut(jd)

            diff = abs(ayan_swe - ayan_py)

            assert diff < STRICT_TOLERANCE, (
                f"SE_SIDM_USER at {date_desc}: diff {diff:.6f}° exceeds tolerance "
                f"(swe={ayan_swe:.6f}, py={ayan_py:.6f})"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_user_sidereal_planet_positions(self, planet_id, planet_name):
        """Test sidereal planetary positions with SE_SIDM_USER match pyswisseph."""
        jd = 2451545.0  # J2000
        t0 = jd
        ayan_t0 = 23.5

        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        pos_swe, _ = swe.calc_ut(jd, planet_id, SEFLG_SIDEREAL)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, SEFLG_SIDEREAL)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < STRICT_TOLERANCE, (
            f"{planet_name} sidereal (SE_SIDM_USER): diff {diff_lon:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_changing_parameters_changes_results(self):
        """Verify that changing t0/ayan_t0 parameters changes results."""
        jd = 2451545.0  # J2000

        # First set of parameters
        t0_1 = jd
        ayan_t0_1 = 23.5

        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0_1, ayan_t0=ayan_t0_1)
        ayan_1 = ephem.swe_get_ayanamsa_ut(jd)

        # Second set with different ayan_t0
        t0_2 = jd
        ayan_t0_2 = 25.0

        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0_2, ayan_t0=ayan_t0_2)
        ayan_2 = ephem.swe_get_ayanamsa_ut(jd)

        # Results should be different
        diff = abs(ayan_1 - ayan_2)
        assert diff > 1.0, (
            f"Changing ayan_t0 should change result: got {ayan_1:.6f} vs {ayan_2:.6f}"
        )

        # Third set with different t0
        t0_3 = 2415020.5  # J1900
        ayan_t0_3 = 23.5

        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0_3, ayan_t0=ayan_t0_3)
        ayan_3 = ephem.swe_get_ayanamsa_ut(jd)

        # Should differ from first (same ayan_t0 but different t0)
        diff_t0 = abs(ayan_1 - ayan_3)
        assert diff_t0 > 0.1, (
            f"Changing t0 should change result: got {ayan_1:.6f} vs {ayan_3:.6f}"
        )

    @pytest.mark.comparison
    def test_user_ayanamsha_different_reference_epochs(self):
        """Test SE_SIDM_USER with different reference epochs matches pyswisseph."""
        # Test with t0 at J1900
        t0 = 2415020.5  # J1900
        ayan_t0 = 22.5
        jd = 2451545.0  # Test at J2000

        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        assert diff < STRICT_TOLERANCE, (
            f"SE_SIDM_USER with J1900 reference: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_user_sidereal_positions_at_future_date(self):
        """Test SE_SIDM_USER sidereal positions at future date match pyswisseph."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 23.5
        jd = 2469807.5  # 2050-01-01

        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        for planet_id, planet_name in TEST_PLANETS[:2]:  # Sun and Moon
            pos_swe, _ = swe.calc_ut(jd, planet_id, SEFLG_SIDEREAL)
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, SEFLG_SIDEREAL)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            assert diff_lon < STRICT_TOLERANCE, (
                f"{planet_name} sidereal (SE_SIDM_USER) at 2050: "
                f"diff {diff_lon:.6f}° exceeds tolerance"
            )
