"""
Moshier Sidereal/Ayanamsha Comparison Tests.

Validates that sidereal calculations using the Moshier ephemeris (SEFLG_MOSEPH)
produce consistent results between pyswisseph (C) and libephemeris (Python).

This is the Moshier-mode mirror of test_compare_sidereal.py (which covers
SEFLG_SWIEPH / JPL mode). It ensures that:
  - All 28 formula-based ayanamsha modes produce matching ayanamsha values
  - All 15 star-based ayanamsha modes produce matching ayanamsha values
  - Sidereal planetary positions (SEFLG_MOSEPH | SEFLG_SIDEREAL) match
  - Major ayanamshas (Lahiri, Fagan-Bradley, Raman) are within expected ranges
  - Ayanamsha progresses correctly over time (precession)
  - J2000 mode behaves correctly at epoch
  - User-defined ayanamsha (SE_SIDM_USER) works with Moshier

Note on star-based ayanamshas:
    In pyswisseph's C implementation, SEFLG_MOSEPH with star-based ayanamshas
    (True Citra, True Revati, etc.) uses the Moshier fixed-star routines to
    compute stellar positions. In libephemeris, star-based ayanamshas use
    Skyfield's star catalog regardless of ephemeris flag. Both approaches are
    valid but may produce slightly different results, hence the relaxed
    tolerance for star-based modes.

Test coverage: 140+ ayanamsha tests (28 formula x 5 dates + star-based) plus
125 sidereal planetary position tests (5 planets x 5 ayanamshas x 5 dates).
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
    SEFLG_MOSEPH,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# AYANAMSHA MODES (same as test_compare_sidereal.py)
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

# Test dates within Moshier range (-3000 to +3000 CE)
# Includes dates that exercise Moshier's extended range advantage
TEST_DATES = [
    (2451545.0, "J2000.0"),
    (2415020.5, "1900-01-01"),
    (2433282.5, "1950-01-01"),
    (2469807.5, "2050-01-01"),
    (2488069.5, "2100-01-01"),
]

TEST_PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Representative ayanamshas for sidereal planet tests (5 modes)
REPRESENTATIVE_AYANAMSHAS = [
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_J2000, "J2000"),
]

# Tolerances
# Formula-based ayanamsha: tight tolerance (both use same formulas)
STRICT_TOLERANCE = 0.001  # 3.6 arcseconds

# Star-based ayanamsha: relaxed because C Moshier fixed-star routines differ
# from Skyfield's star catalog approach
RELAXED_TOLERANCE = 0.1  # 6 arcminutes

# Moshier sidereal position tolerance: combines Moshier planet tolerance
# (~0.02 deg) with ayanamsha tolerance (~0.001 deg for formula-based)
MOSHIER_SIDEREAL_TOLERANCE = 0.025  # degrees


def get_tolerance(sid_mode: int) -> float:
    """Get appropriate tolerance for ayanamsha mode."""
    star_modes = {m[0] for m in STAR_BASED_AYANAMSHAS}
    if sid_mode in star_modes:
        return RELAXED_TOLERANCE
    return STRICT_TOLERANCE


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMoshierAyanamshaValues:
    """Compare ayanamsha values between implementations with Moshier context.

    Since swe_get_ayanamsa_ut() does not take an ephemeris flag, the ayanamsha
    computation itself is the same regardless of Moshier/SWIEPH. However, this
    test validates that the ayanamsha values are correct in the context of
    Moshier calculations (same dates, same sid_modes), ensuring consistency
    for the full sidereal pipeline.
    """

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
            f"{sid_name} at {date_desc}: diff {diff:.6f}° exceeds tolerance "
            f"(swe={ayan_swe:.6f}°, lib={ayan_py:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", STAR_BASED_AYANAMSHAS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:2])
    def test_star_based_ayanamsha(self, sid_mode, sid_name, jd, date_desc):
        """Test star-based ayanamsha values match within relaxed tolerance.

        Star-based modes (True Citra, True Revati, etc.) may differ between
        pyswisseph's C Moshier fixed-star routines and libephemeris's Skyfield
        star catalog, hence the relaxed tolerance.
        """
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.swe_get_ayanamsa_ut(jd)

        diff = angular_diff(ayan_swe, ayan_py)

        assert diff < RELAXED_TOLERANCE, (
            f"{sid_name} at {date_desc}: diff {diff:.6f}° exceeds tolerance "
            f"(swe={ayan_swe:.6f}°, lib={ayan_py:.6f}°)"
        )


class TestMoshierSiderealPlanets:
    """Compare sidereal planetary positions with SEFLG_MOSEPH.

    Tests 5 planets x 5 ayanamshas x 5 dates = 125 test cases.
    This validates the full Moshier sidereal pipeline:
      1. Moshier planet calculation (VSOP87/ELP)
      2. Ayanamsha computation (formula-based)
      3. Sidereal correction (tropical_lon - ayanamsha)
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", REPRESENTATIVE_AYANAMSHAS)
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_moshier_sidereal_positions(
        self, sid_mode, sid_name, planet_id, planet_name, jd, date_desc
    ):
        """Test Moshier sidereal planetary positions match pyswisseph."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < MOSHIER_SIDEREAL_TOLERANCE, (
            f"{planet_name} sidereal ({sid_name}) at {date_desc}: "
            f"diff {diff_lon:.6f}° exceeds tolerance {MOSHIER_SIDEREAL_TOLERANCE}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestMoshierMajorAyanamshas:
    """Detailed tests for the most commonly used ayanamshas with Moshier."""

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

    @pytest.mark.comparison
    def test_lahiri_sidereal_sun_moshier(self):
        """Test Lahiri sidereal Sun position using Moshier at J2000."""
        jd = 2451545.0

        swe.set_sid_mode(SE_SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        pos_swe, _ = swe.calc_ut(jd, SE_SUN, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < MOSHIER_SIDEREAL_TOLERANCE, (
            f"Lahiri sidereal Sun (Moshier): diff {diff_lon:.6f}° exceeds tolerance "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )


class TestMoshierAyanamshaProgression:
    """Test that ayanamsha changes correctly over time with Moshier dates."""

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

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", REPRESENTATIVE_AYANAMSHAS[:3])
    def test_moshier_sidereal_progression(self, sid_mode, sid_name):
        """Test that Moshier sidereal positions change consistently over time."""
        jd_1950 = 2433282.5
        jd_2050 = 2469807.5

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        # Sun position at two dates
        pos_1950_swe, _ = swe.calc_ut(jd_1950, SE_SUN, flag_swe)
        pos_2050_swe, _ = swe.calc_ut(jd_2050, SE_SUN, flag_swe)
        pos_1950_py, _ = ephem.swe_calc_ut(jd_1950, SE_SUN, flag_py)
        pos_2050_py, _ = ephem.swe_calc_ut(jd_2050, SE_SUN, flag_py)

        # Both libraries should show similar motion over 100 years
        motion_swe = angular_diff(pos_2050_swe[0], pos_1950_swe[0])
        motion_py = angular_diff(pos_2050_py[0], pos_1950_py[0])

        motion_diff = abs(motion_swe - motion_py)
        assert motion_diff < 0.05, (
            f"{sid_name}: Sun sidereal motion over 100yr differs by {motion_diff:.6f}° "
            f"(swe={motion_swe:.6f}°, lib={motion_py:.6f}°)"
        )


class TestMoshierJ2000Mode:
    """Test J2000 ayanamsha mode with Moshier."""

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

        # Should be negative before J2000, positive after
        assert ayan_1900 < 0, f"J2000 at 1900 should be negative, got {ayan_1900}"
        assert abs(ayan_2000) < 0.001, f"J2000 at 2000 should be ~0, got {ayan_2000}"
        assert ayan_2100 > 0, f"J2000 at 2100 should be positive, got {ayan_2100}"

    @pytest.mark.comparison
    def test_j2000_sidereal_planet_moshier(self):
        """Test J2000 sidereal planet position with Moshier at J2000 epoch.

        At J2000 epoch, the J2000 ayanamsha is 0, so sidereal and tropical
        positions should be nearly identical.
        """
        jd = 2451545.0

        swe.set_sid_mode(SE_SIDM_J2000)
        ephem.swe_set_sid_mode(SE_SIDM_J2000)

        # Moshier sidereal
        flag_sid = SEFLG_MOSEPH | SEFLG_SIDEREAL
        pos_sid, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_sid)

        # Moshier tropical
        flag_trop = SEFLG_MOSEPH
        pos_trop, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_trop)

        # At J2000, J2000 ayanamsha ~0, so positions should be very close
        # (small nutation correction still applies via _get_true_ayanamsa)
        diff = angular_diff(pos_sid[0], pos_trop[0])
        assert diff < 0.01, (
            f"J2000 sidereal vs tropical at epoch: diff {diff:.6f}° "
            f"(sid={pos_sid[0]:.6f}°, trop={pos_trop[0]:.6f}°)"
        )


class TestMoshierUserDefinedAyanamsha:
    """Test SE_SIDM_USER custom ayanamsha mode with Moshier."""

    @pytest.mark.comparison
    def test_user_ayanamsha_value_at_j2000(self):
        """Test SE_SIDM_USER ayanamsha value matches pyswisseph at J2000."""
        jd = 2451545.0  # J2000
        t0 = jd
        ayan_t0 = 23.5  # Custom ayanamsha value

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
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS[:2])
    def test_user_sidereal_planet_moshier(self, planet_id, planet_name):
        """Test sidereal planetary positions with SE_SIDM_USER and SEFLG_MOSEPH."""
        jd = 2451545.0  # J2000
        t0 = jd
        ayan_t0 = 23.5

        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < MOSHIER_SIDEREAL_TOLERANCE, (
            f"{planet_name} sidereal (SE_SIDM_USER, Moshier): "
            f"diff {diff_lon:.6f}° exceeds tolerance "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )

    @pytest.mark.comparison
    def test_user_sidereal_positions_at_future_date(self):
        """Test SE_SIDM_USER Moshier sidereal positions at future date."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 23.5
        jd = 2469807.5  # 2050-01-01

        swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)
        ephem.swe_set_sid_mode(SE_SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        for planet_id, planet_name in TEST_PLANETS[:2]:  # Sun and Moon
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            assert diff_lon < MOSHIER_SIDEREAL_TOLERANCE, (
                f"{planet_name} sidereal (SE_SIDM_USER, Moshier) at 2050: "
                f"diff {diff_lon:.6f}° exceeds tolerance"
            )
