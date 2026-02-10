"""
Moshier Time Functions Comparison Tests.

Compares Delta T and sidereal time between pyswisseph (C library) and
libephemeris (Python) for dates across the extended Moshier range (-3000
to +3000 CE).

This is the Moshier-mode mirror of test_compare_time.py (TestDeltat /
TestSiderealTime sections) and addresses the fact that libephemeris uses
Skyfield's Delta T model (Stephenson-Morrison-Hohenkerk 2016 + IERS)
while the C library uses Morrison-Stephenson with its own polynomial
fit.  For dates before 1620 CE the two models can diverge significantly
(up to hundreds of seconds), so graduated tolerances are applied per era.

Key design notes:
  - libephemeris swe_deltat_ex with SEFLG_MOSEPH returns the same value
    as SEFLG_SWIEPH (Skyfield Delta T). The flag only affects position
    calculations.  See time_utils.py:269.
  - pyswisseph deltat_ex returns a plain float (Delta T in days).
  - libephemeris swe_deltat_ex returns (delta_t, serr) tuple.
  - Both sidtime functions return hours as a float.

Test matrix:
  - 20 Delta T test cases parametrised across the Moshier date range.
  - 10 sidereal time test cases parametrised similarly.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SEFLG_MOSEPH


# ============================================================================
# TOLERANCES
# ============================================================================

# Graduated Delta T tolerances in seconds, keyed by era.
# These reflect known systematic differences between the Skyfield model
# (Stephenson-Morrison-Hohenkerk 2016) and the Swiss Ephemeris C library
# model (Morrison-Stephenson polynomial + Espenak-Meeus).
#
# Empirical observations (C vs Python, Jan 1 noon UT):
#   -2999:  ~780s   |  -2000: ~400s  |  -1000: ~100s  |  -500: ~25s
#   0:      ~40s    |  200:   ~40s   |  500:   ~80s   |  800:  ~183s
#   1000:   ~205s   |  1200:  ~124s  |  1400:  ~26s   |  1582: ~3s
#   1700:   ~0.9s   |  1800:  ~0.7s  |  1900:  ~0.01s |  2000: ~0.02s
#   2024:   ~0.01s  |  2100:  ~2.7s  |  3000:  ~875s
#
# The divergence grows quadratically with distance from the IERS-covered
# era (~1900-2024) because both models extrapolate Earth rotation with
# different polynomial coefficients.
DELTA_T_TOL_BY_ERA = {
    "modern": 0.5,  # 1900-2050 CE  – IERS / observed data overlap
    "classical": 5.0,  # 1620-1900 CE  – telescopic era
    "medieval": 250.0,  # 1000-1620 CE  – sparse historical records
    "ancient": 200.0,  # 0-1000 CE     – very sparse / eclipse records
    "deep": 800.0,  # < 0 CE        – pure extrapolation
    "future": 1000.0,  # > 2050 CE     – forward extrapolation
}

# Sidereal time tolerance in hours.
# Differences arise because GST depends on Delta T through the UT->TT
# conversion. A 1000 s Delta T error maps to ~0.07 h rotation. For ancient
# dates we allow up to 0.15 h to absorb the cascading effect.
SIDTIME_TOL_MODERN = 0.001  # hours – 1900-2050
SIDTIME_TOL_HISTORICAL = 0.05  # hours – 1000-1900
SIDTIME_TOL_ANCIENT = 0.15  # hours – < 1000 CE
SIDTIME_TOL_FUTURE = 0.003  # hours – > 2050 CE


# ============================================================================
# HELPER
# ============================================================================


def _year_to_jd(year: int) -> float:
    """Convert a proleptic year to Julian Day (Jan 1, 12:00 UT).

    Uses Gregorian calendar for year >= 1582, Julian otherwise.
    Astronomical year numbering: 1 BCE = year 0, 2 BCE = year -1, etc.
    """
    cal = 1 if year >= 1582 else 0  # 1 = Gregorian, 0 = Julian
    return swe.julday(year, 1, 1, 12.0, cal)


def _delta_t_tolerance(year: int) -> float:
    """Return the appropriate Delta T tolerance in seconds for *year*."""
    if year < 0:
        return DELTA_T_TOL_BY_ERA["deep"]
    if year < 1000:
        return DELTA_T_TOL_BY_ERA["ancient"]
    if year < 1620:
        return DELTA_T_TOL_BY_ERA["medieval"]
    if year < 1900:
        return DELTA_T_TOL_BY_ERA["classical"]
    if year <= 2050:
        return DELTA_T_TOL_BY_ERA["modern"]
    return DELTA_T_TOL_BY_ERA["future"]


def _sidtime_tolerance(year: int) -> float:
    """Return the appropriate sidereal time tolerance in hours for *year*."""
    if year < 1000:
        return SIDTIME_TOL_ANCIENT
    if year < 1900:
        return SIDTIME_TOL_HISTORICAL
    if year <= 2050:
        return SIDTIME_TOL_MODERN
    return SIDTIME_TOL_FUTURE


# ============================================================================
# TEST DATA
# ============================================================================

# 20 Delta T test dates spanning the Moshier range -2999 .. +3000 CE.
# Each entry: (year, description)
DELTAT_TEST_DATES = [
    # Deep past (< 0 CE)
    (-2999, "3000 BCE – Moshier start"),
    (-2000, "2001 BCE – deep antiquity"),
    (-1000, "1001 BCE – Bronze Age"),
    (-500, "501 BCE – Classical Greece"),
    # Ancient (0-1000 CE)
    (0, "1 BCE (year 0) – turn of era"),
    (200, "200 CE – Roman Empire"),
    (500, "500 CE – fall of Rome"),
    (800, "800 CE – Carolingian"),
    # Medieval (1000-1620 CE)
    (1000, "1000 CE – High Middle Ages"),
    (1200, "1200 CE – Late Middle Ages"),
    (1400, "1400 CE – Renaissance"),
    (1582, "1582 CE – Gregorian reform"),
    # Classical (1620-1900 CE)
    (1700, "1700 CE – early telescopic"),
    (1800, "1800 CE – 19th century"),
    # Modern (1900-2100 CE)
    (1900, "1900 CE – 20th century start"),
    (1950, "1950 CE – mid-century"),
    (2000, "2000 CE – J2000 epoch"),
    (2024, "2024 CE – present"),
    # Future
    (2100, "2100 CE – near future"),
    (3000, "3000 CE – Moshier end"),
]

# 10 sidereal time test dates spanning the Moshier range.
SIDTIME_TEST_DATES = [
    (-2000, "2001 BCE – deep antiquity"),
    (-500, "501 BCE – Classical Greece"),
    (0, "1 BCE (year 0) – turn of era"),
    (500, "500 CE – fall of Rome"),
    (1000, "1000 CE – High Middle Ages"),
    (1582, "1582 CE – Gregorian reform"),
    (1900, "1900 CE – 20th century start"),
    (2000, "2000 CE – J2000 epoch"),
    (2024, "2024 CE – present"),
    (3000, "3000 CE – Moshier end"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestDeltat:
    """Compare Delta T (swe_deltat_ex with SEFLG_MOSEPH) across the Moshier range.

    libephemeris always uses the Skyfield model regardless of the ephemeris
    flag; the C library may use a Moshier-specific polynomial for Delta T.
    Graduated tolerances accommodate the growing model divergence for
    ancient dates.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,desc", DELTAT_TEST_DATES)
    def test_deltat_ex_moseph(self, year, desc):
        """Compare deltat_ex(jd, SEFLG_MOSEPH) between C and Python."""
        jd = _year_to_jd(year)
        tol = _delta_t_tolerance(year)

        # pyswisseph returns float (Delta T in days)
        dt_swe = swe.deltat_ex(jd, swe.FLG_MOSEPH)

        # libephemeris returns (delta_t, serr) tuple
        dt_py, serr = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

        diff_days = abs(dt_swe - dt_py)
        diff_seconds = diff_days * 86400.0

        # Report values for diagnostic purposes
        dt_swe_sec = dt_swe * 86400.0
        dt_py_sec = dt_py * 86400.0

        assert diff_seconds < tol, (
            f"{desc}: Delta T diff {diff_seconds:.3f}s exceeds "
            f"tolerance {tol}s "
            f"(swe={dt_swe_sec:.3f}s, lib={dt_py_sec:.3f}s)"
        )


class TestSiderealTime:
    """Compare sidereal time across the Moshier range.

    Both implementations compute GAST (Greenwich Apparent Sidereal Time)
    from GMST + equation of the equinoxes.  Differences in Delta T
    propagate through the UT1->TT conversion and affect the result.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,desc", SIDTIME_TEST_DATES)
    def test_sidtime_moseph(self, year, desc):
        """Compare swe.sidtime(jd) vs ephem.sidtime(jd) for Moshier dates."""
        jd = _year_to_jd(year)
        tol = _sidtime_tolerance(year)

        st_swe = swe.sidtime(jd)
        st_py = ephem.sidtime(jd)

        diff = abs(st_swe - st_py)
        # Handle wrap-around at 24 hours
        if diff > 12:
            diff = 24 - diff

        assert diff < tol, (
            f"{desc}: sidtime diff {diff:.6f}h exceeds "
            f"tolerance {tol}h "
            f"(swe={st_swe:.6f}h, lib={st_py:.6f}h)"
        )
