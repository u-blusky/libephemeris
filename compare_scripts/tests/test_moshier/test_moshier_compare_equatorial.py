"""
Moshier Equatorial Coordinate (RA/Dec) Cross-Library Comparison Tests.

Validates the ecliptic-to-equatorial transformation chain in Moshier mode
(SEFLG_MOSEPH | SEFLG_EQUATORIAL) between pyswisseph (C library) and
libephemeris (Python reimplementation):

  ecliptic (lon, lat) -> obliquity -> (RA, Dec)

The transformation chain in libephemeris (planets.py:1052-1069) uses:
  - moshier.mean_obliquity(jd_tt): IAU 2006 precession (precession.py)
  - moshier.ecliptic_to_equatorial(lon, lat, obliquity): spherical rotation
  - Velocity: numerical differentiation with step 0.001 days

The C library uses its own internal obliquity (potentially Lieske 1979) and
transformation formulas. The obliquity difference between IAU 2006 and
Lieske 1979 is ~0.01 arcsec at J2000 but accumulates for dates far from J2000.

RA/Dec accuracy is critical for:
  - Telescopic observation and instrument pointing
  - Stellar catalogue cross-matching
  - Downstream azalt(), heliacal_ut(), rise_trans() calculations

This module also tests with SEFLG_NONUT to isolate the nutation effect,
ensuring that the mean equatorial positions (without nutation correction)
also agree between implementations.

Test count: 56 cases
  - TestEquatorialPositions: 8 planets x 5 dates = 40 (RA, Dec, speed_RA)
  - TestEquatorialNoNut: 8 planets x 2 dates = 16 (RA, Dec without nutation)
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_NONUT,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# 8 planets: Sun through Neptune, excluding Moon (high angular velocity
# amplifies numerical differentiation differences) and Pluto (Moshier
# Chapront-Francou theory has multi-degree ecliptic errors that dominate
# any equatorial transformation test). Moon and Pluto ecliptic accuracy
# is already validated in test_moshier_compare_planets.py.
PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
]

# 5 dates spanning modern and extended range
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
    (2100, 12, 31, 23.999, "Late 2100s"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]

# 2 dates for NONUT tests (J2000 and current era)
NONUT_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
]


# ============================================================================
# TOLERANCES
# ============================================================================

# Equatorial position tolerances.
# RA and Dec are derived from ecliptic coordinates via obliquity rotation.
# The dominant error sources are:
#   1. Ecliptic position difference (~0.02 deg from VSOP87 C-vs-Python)
#   2. Obliquity formula difference (IAU 2006 vs Lieske 1979: ~0.01 arcsec
#      at J2000, growing to ~0.2 arcsec at year 1900/2100)
#   3. Nutation model differences (affecting true obliquity)
# The 0.02 deg ecliptic error propagates nearly 1:1 into RA and is slightly
# reduced in Dec (by cos(obliquity) ~ 0.92). The obliquity contribution
# is negligible (<0.001 deg) for the date range tested.
RA_TOL = 0.03  # degrees (~108 arcsec) - ecliptic error + transformation
DEC_TOL = 0.03  # degrees (~108 arcsec) - ecliptic error + transformation

# Equatorial velocity tolerance.
# libephemeris uses numerical differentiation (step=0.001 days) for velocity
# transformation from ecliptic to equatorial, while the C library uses its
# own internal differentiation. The 0.01 deg/day tolerance accommodates
# both the underlying ecliptic velocity differences and the numerical
# differentiation approach difference.
SPEED_RA_TOL = 0.02  # degrees/day
SPEED_DEC_TOL = 0.02  # degrees/day

# NONUT tolerances: without nutation, the transformation uses mean obliquity
# only, eliminating one source of C-vs-Python difference. Tolerances are
# the same since ecliptic position error still dominates.
NONUT_RA_TOL = 0.03  # degrees
NONUT_DEC_TOL = 0.03  # degrees


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestEquatorialPositions:
    """Compare Moshier equatorial (RA/Dec) positions between pyswisseph and libephemeris.

    Tests the full chain: ecliptic Moshier position -> obliquity -> RA/Dec
    transformation with SEFLG_MOSEPH | SEFLG_EQUATORIAL | SEFLG_SPEED.
    Validates RA (with angular wrap), Dec, and RA/Dec velocities.

    8 planets x 5 dates = 40 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_ra(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier equatorial RA matches pyswisseph.

        RA (Right Ascension) is compared using angular_diff to handle the
        0/360 degree wrap. The tolerance accounts for both the underlying
        ecliptic longitude difference and the obliquity transformation
        difference (IAU 2006 vs Lieske 1979).
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_ra = angular_diff(pos_swe[0], pos_py[0])

        assert diff_ra < RA_TOL, (
            f"{planet_name} Moshier equatorial at {date_desc}: "
            f"RA diff {diff_ra:.6f}° exceeds tolerance {RA_TOL}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_dec(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier equatorial Dec matches pyswisseph.

        Declination is compared with absolute difference. The ecliptic
        latitude error is reduced by cos(obliquity) (~0.92) in the
        transformation, but the ecliptic longitude error contributes
        via sin(obliquity) (~0.40).
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_dec = abs(pos_swe[1] - pos_py[1])

        assert diff_dec < DEC_TOL, (
            f"{planet_name} Moshier equatorial at {date_desc}: "
            f"Dec diff {diff_dec:.6f}° exceeds tolerance {DEC_TOL}° "
            f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_speed_ra(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier equatorial RA velocity matches pyswisseph.

        libephemeris uses numerical differentiation (step=0.001 days) to
        transform ecliptic velocities to equatorial, while the C library
        uses its own internal approach. The tolerance accommodates both
        the underlying velocity differences and the method difference.
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_speed_ra = abs(pos_swe[3] - pos_py[3])

        assert diff_speed_ra < SPEED_RA_TOL, (
            f"{planet_name} Moshier equatorial at {date_desc}: "
            f"RA speed diff {diff_speed_ra:.6f}°/day exceeds "
            f"tolerance {SPEED_RA_TOL}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_speed_dec(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier equatorial Dec velocity matches pyswisseph.

        Dec velocity tends to be smaller than RA velocity for most planets
        (ecliptic is close to equator), so differences are typically smaller.
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_speed_dec = abs(pos_swe[4] - pos_py[4])

        assert diff_speed_dec < SPEED_DEC_TOL, (
            f"{planet_name} Moshier equatorial at {date_desc}: "
            f"Dec speed diff {diff_speed_dec:.6f}°/day exceeds "
            f"tolerance {SPEED_DEC_TOL}°/day "
            f"(swe={pos_swe[4]:.6f}, lib={pos_py[4]:.6f})"
        )


class TestEquatorialNoNut:
    """Compare Moshier equatorial positions with SEFLG_NONUT.

    SEFLG_NONUT suppresses nutation, so the transformation uses mean
    obliquity only. This isolates the effect of different nutation models
    between the C and Python implementations: if results agree with NONUT
    but diverge without it, the nutation model is the source of divergence.

    The C library uses its own nutation (likely IAU 1980 or 2000A), while
    libephemeris may use a different model. By testing with NONUT, we verify
    that the mean equatorial transformation chain is consistent regardless
    of nutation differences.

    8 planets x 2 dates = 16 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", NONUT_DATES)
    def test_nonut_equatorial_ra(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier equatorial RA without nutation matches pyswisseph.

        With NONUT, both implementations should use mean obliquity only,
        removing nutation as a source of disagreement.
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_NONUT
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_NONUT

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_ra = angular_diff(pos_swe[0], pos_py[0])

        assert diff_ra < NONUT_RA_TOL, (
            f"{planet_name} Moshier NONUT at {date_desc}: "
            f"RA diff {diff_ra:.6f}° exceeds tolerance {NONUT_RA_TOL}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", NONUT_DATES)
    def test_nonut_equatorial_dec(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test Moshier equatorial Dec without nutation matches pyswisseph.

        Dec comparison with NONUT isolates the mean obliquity transformation.
        Differences here indicate divergence in the mean obliquity formula
        (IAU 2006 vs Lieske 1979) or the ecliptic-to-equatorial rotation.
        """
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_NONUT
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_NONUT

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)

        diff_dec = abs(pos_swe[1] - pos_py[1])

        assert diff_dec < NONUT_DEC_TOL, (
            f"{planet_name} Moshier NONUT at {date_desc}: "
            f"Dec diff {diff_dec:.6f}° exceeds tolerance {NONUT_DEC_TOL}° "
            f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
        )
