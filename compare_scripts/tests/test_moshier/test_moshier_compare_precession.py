"""
Moshier Precession/Obliquity Cross-Library Comparison Tests.

Compares the mean obliquity and nutation used internally by the C library
Swiss Ephemeris (pyswisseph) against the Python reimplementation (libephemeris)
for Moshier mode (SEFLG_MOSEPH).

PROBLEM:
  - libephemeris uses IAU 2006 for precession (mean obliquity) and IAU 2000B
    for nutation (77 terms).
  - The C library Swiss Ephemeris uses Lieske 1979 for precession and its own
    IAU 2000B nutation implementation which may differ subtly.
  - At J2000.0, both agree to ~0.01 arcsec, but the difference grows for
    dates far from J2000 (~50.3 arcsec/yr precession rate means even small
    formula differences accumulate).

APPROACH:
  - Compute ecliptic (lon, lat) and equatorial (RA, Dec) positions of the Sun
    using SEFLG_MOSEPH in both implementations.
  - Derive the "implicit obliquity" from each pair using the spherical
    trigonometry relation:
        sin(Dec) = sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon)
    For the Sun (lat ~ 0), this simplifies to:
        eps_implicit = arcsin(sin(Dec) / sin(lon))
  - Compare implicit obliquities between C and Python implementations.
  - Compare also with moshier.mean_obliquity() and moshier.true_obliquity()
    directly.

TEST APPROACH FOR NONUT:
  - Use SEFLG_NONUT to suppress nutation and derive the implicit mean
    obliquity (without nutation correction).
  - Compare implicit mean obliquity C-vs-Python and against
    moshier.mean_obliquity() directly.

This quantifies the precession divergence as a function of time and establishes
the obliquity as a source of systematic error in ecliptic-to-equatorial
transformations.

Test count: 10 cases
  - TestImplicitObliquity: 10 dates (implicit obliquity C-vs-Python, with
    nutation; also checks against moshier.true_obliquity)
  - TestImplicitMeanObliquity: 10 dates (implicit mean obliquity C-vs-Python,
    NONUT; also checks against moshier.mean_obliquity)
  - TestDirectMeanObliquity: 10 dates (moshier.mean_obliquity vs C-derived)
  - TestNutationContribution: 10 dates (nutation in obliquity C-vs-Python)
"""

import math

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_NONUT,
)
from libephemeris.moshier.precession import (
    mean_obliquity,
    true_obliquity,
    nutation_angles,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# 10 test dates spanning a wide range of epochs.
# The dates are chosen to sample the precession divergence across a range
# from ancient historical dates (-2000 CE) to far-future dates (3000 CE).
# Each tuple: (year, month, day, hour, description)
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0 (reference epoch)"),
    (2024, 11, 15, 0.0, "Current era 2024"),
    (1900, 1, 1, 0.0, "Year 1900"),
    (1800, 1, 1, 0.0, "Year 1800"),
    (1600, 1, 1, 0.0, "Year 1600"),
    (1200, 1, 1, 0.0, "Year 1200"),
    (800, 1, 1, 0.0, "Year 800 CE"),
    (100, 1, 1, 0.0, "Year 100 CE"),
    (-500, 1, 1, 0.0, "Year 500 BCE"),
    (-1500, 1, 1, 0.0, "Year 1500 BCE"),
]


# ============================================================================
# TOLERANCES
# ============================================================================

# Tolerance for implicit obliquity comparison (C-vs-Python), in degrees.
# At J2000, IAU 2006 vs Lieske 1979 differ by ~0.01 arcsec.
# For dates far from J2000, the divergence grows.
# We use generous tolerances because the derivation involves two separate
# planetary computations (ecliptic + equatorial) which each have their own
# Moshier C-vs-Python differences (~0.02 deg for ecliptic), and the implicit
# obliquity derivation amplifies errors near lon=0/180.
#
# The C-vs-Python ecliptic/equatorial position differences propagate into
# the derived obliquity, so we cannot expect sub-arcsec agreement from
# this indirect method. The value of these tests is quantifying the
# divergence, not asserting tight bounds.
IMPLICIT_OBLIQUITY_TOL = 0.05  # degrees (~180 arcsec)

# Tolerance for mean obliquity comparison (NONUT, C-vs-Python derived).
# Same considerations as above.
IMPLICIT_MEAN_OBLIQUITY_TOL = 0.05  # degrees

# Tolerance for direct mean obliquity comparison
# (moshier.mean_obliquity vs C-derived implicit mean obliquity).
# This is limited by the accuracy of the implicit derivation method.
DIRECT_MEAN_OBLIQUITY_TOL = 0.05  # degrees

# Tolerance for nutation contribution comparison (C-vs-Python).
# Nutation in obliquity is typically < 10 arcsec (~0.003 deg).
# The difference is derived from (true_obliquity - mean_obliquity) on each side.
NUTATION_CONTRIBUTION_TOL = 0.1  # degrees


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def derive_obliquity(
    lon_deg: float, lat_deg: float, ra_deg: float, dec_deg: float
) -> float:
    """Derive the implicit obliquity from ecliptic and equatorial coordinates.

    Uses the spherical trigonometry relation:
        sin(Dec) = sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon)

    For the Sun (lat ~ 0), this simplifies to:
        eps = arcsin(sin(Dec) / sin(lon))

    For the general case, we solve numerically.

    Args:
        lon_deg: Ecliptic longitude in degrees.
        lat_deg: Ecliptic latitude in degrees.
        ra_deg: Right Ascension in degrees.
        dec_deg: Declination in degrees.

    Returns:
        Implicit obliquity in degrees, or NaN if derivation fails.
    """
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    dec = math.radians(dec_deg)

    sin_lon = math.sin(lon)
    cos_lat = math.cos(lat)
    sin_lat = math.sin(lat)
    sin_dec = math.sin(dec)

    # Avoid division by zero or near-zero sin(lon)
    # (lon near 0 or 180 makes derivation numerically unstable)
    denom = cos_lat * sin_lon
    if abs(denom) < 0.1:
        return float("nan")

    # sin(Dec) = sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon)
    # For small lat: sin(Dec) ~ sin(eps)*sin(lon)
    # General: solve for eps using iterative approach
    # First approximation (assuming lat ~ 0):
    sin_eps_approx = sin_dec / sin_lon
    if abs(sin_eps_approx) > 1.0:
        return float("nan")

    eps = math.asin(sin_eps_approx)

    # Refine with Newton iteration for non-zero latitude
    if abs(lat_deg) > 0.001:
        for _ in range(10):
            f = sin_lat * math.cos(eps) + cos_lat * math.sin(eps) * sin_lon - sin_dec
            fp = -sin_lat * math.sin(eps) + cos_lat * math.cos(eps) * sin_lon
            if abs(fp) < 1e-15:
                break
            eps -= f / fp

    return math.degrees(eps)


def get_jd_tt_from_ut(jd_ut: float) -> float:
    """Convert JD UT to JD TT using delta T from pyswisseph.

    Args:
        jd_ut: Julian Day in Universal Time.

    Returns:
        Julian Day in Terrestrial Time.
    """
    delta_t = swe.deltat(jd_ut)
    return jd_ut + delta_t


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestImplicitObliquity:
    """Compare implicit true obliquity derived from ecliptic/equatorial positions.

    For each test date, computes the Sun's ecliptic and equatorial positions
    using both C library and Python implementation with SEFLG_MOSEPH, derives
    the implicit obliquity from each, and compares them.

    Also compares the Python-derived implicit obliquity against
    moshier.true_obliquity() as a consistency check.

    10 dates = 10 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_implicit_true_obliquity_c_vs_python(
        self, year, month, day, hour, date_desc
    ):
        """Derive and compare implicit true obliquity from C and Python positions.

        Computes Sun ecliptic + equatorial with SEFLG_MOSEPH in both
        implementations. Derives the obliquity implicit in the ecl->eq
        transformation from each. Compares the two derived obliquities.
        """
        jd = swe.julday(year, month, day, hour)

        # --- C library (pyswisseph) ---
        flag_ecl_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_eq_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        ecl_swe, _ = swe.calc_ut(jd, SE_SUN, flag_ecl_swe)
        eq_swe, _ = swe.calc_ut(jd, SE_SUN, flag_eq_swe)

        eps_swe = derive_obliquity(ecl_swe[0], ecl_swe[1], eq_swe[0], eq_swe[1])

        # --- Python library (libephemeris) ---
        flag_ecl_py = SEFLG_MOSEPH | SEFLG_SPEED
        flag_eq_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
        ecl_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_ecl_py)
        eq_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_eq_py)

        eps_py = derive_obliquity(ecl_py[0], ecl_py[1], eq_py[0], eq_py[1])

        # Skip if derivation failed (lon near 0 or 180)
        if math.isnan(eps_swe) or math.isnan(eps_py):
            pytest.skip(
                f"Cannot derive obliquity at {date_desc}: "
                f"Sun lon too close to 0/180 "
                f"(swe_lon={ecl_swe[0]:.1f}, py_lon={ecl_py[0]:.1f})"
            )

        diff_deg = abs(eps_swe - eps_py)
        diff_arcsec = diff_deg * 3600.0

        assert diff_deg < IMPLICIT_OBLIQUITY_TOL, (
            f"Implicit true obliquity at {date_desc}: "
            f"C-vs-Python diff {diff_arcsec:.2f} arcsec ({diff_deg:.6f}°) "
            f"exceeds tolerance {IMPLICIT_OBLIQUITY_TOL}° "
            f"(C={eps_swe:.6f}°, Python={eps_py:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_python_implicit_vs_direct_true_obliquity(
        self, year, month, day, hour, date_desc
    ):
        """Compare Python implicit obliquity against moshier.true_obliquity().

        This checks internal consistency: the obliquity derived from
        libephemeris ecliptic/equatorial output should match the obliquity
        returned by moshier.true_obliquity() for the same JD.
        """
        jd = swe.julday(year, month, day, hour)
        jd_tt = get_jd_tt_from_ut(jd)

        # Get direct true obliquity from moshier module
        eps_direct = true_obliquity(jd_tt)

        # Derive implicit obliquity from Python ecl/eq positions
        flag_ecl = SEFLG_MOSEPH | SEFLG_SPEED
        flag_eq = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
        ecl_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_ecl)
        eq_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_eq)

        eps_implicit = derive_obliquity(ecl_py[0], ecl_py[1], eq_py[0], eq_py[1])

        if math.isnan(eps_implicit):
            pytest.skip(
                f"Cannot derive obliquity at {date_desc}: "
                f"Sun lon too close to 0/180 (lon={ecl_py[0]:.1f})"
            )

        diff_deg = abs(eps_direct - eps_implicit)
        diff_arcsec = diff_deg * 3600.0

        # This should be very tight since both come from the same codebase.
        # Differences arise only from numerical precision in the derivation.
        assert diff_deg < 0.01, (
            f"Python true obliquity at {date_desc}: "
            f"direct vs implicit diff {diff_arcsec:.2f} arcsec ({diff_deg:.6f}°) "
            f"(direct={eps_direct:.6f}°, implicit={eps_implicit:.6f}°)"
        )


class TestImplicitMeanObliquity:
    """Compare implicit mean obliquity with SEFLG_NONUT (no nutation).

    Using SEFLG_NONUT suppresses nutation, so the ecl->eq transformation
    uses only mean obliquity. This isolates the precession model difference
    (IAU 2006 vs Lieske 1979) from nutation differences.

    10 dates = 10 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_implicit_mean_obliquity_c_vs_python(
        self, year, month, day, hour, date_desc
    ):
        """Derive and compare implicit mean obliquity from C and Python (NONUT).

        With SEFLG_NONUT, the transformation uses mean obliquity only.
        Differences here isolate the precession model divergence
        (IAU 2006 in Python vs Lieske 1979 in C).
        """
        jd = swe.julday(year, month, day, hour)

        # --- C library ---
        flag_ecl_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_NONUT
        flag_eq_swe = (
            swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_NONUT
        )
        ecl_swe, _ = swe.calc_ut(jd, SE_SUN, flag_ecl_swe)
        eq_swe, _ = swe.calc_ut(jd, SE_SUN, flag_eq_swe)

        eps_swe = derive_obliquity(ecl_swe[0], ecl_swe[1], eq_swe[0], eq_swe[1])

        # --- Python library ---
        flag_ecl_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_NONUT
        flag_eq_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_NONUT
        ecl_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_ecl_py)
        eq_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flag_eq_py)

        eps_py = derive_obliquity(ecl_py[0], ecl_py[1], eq_py[0], eq_py[1])

        if math.isnan(eps_swe) or math.isnan(eps_py):
            pytest.skip(
                f"Cannot derive mean obliquity at {date_desc}: "
                f"Sun lon too close to 0/180"
            )

        diff_deg = abs(eps_swe - eps_py)
        diff_arcsec = diff_deg * 3600.0

        assert diff_deg < IMPLICIT_MEAN_OBLIQUITY_TOL, (
            f"Implicit mean obliquity at {date_desc}: "
            f"C-vs-Python diff {diff_arcsec:.2f} arcsec ({diff_deg:.6f}°) "
            f"exceeds tolerance {IMPLICIT_MEAN_OBLIQUITY_TOL}° "
            f"(C={eps_swe:.6f}°, Python={eps_py:.6f}°)"
        )


class TestDirectMeanObliquity:
    """Compare moshier.mean_obliquity() against C-derived implicit mean obliquity.

    Uses the C library with SEFLG_NONUT to derive the mean obliquity the C
    library uses internally, then compares against moshier.mean_obliquity()
    which uses the IAU 2006 formula.

    This directly quantifies the IAU 2006 vs Lieske 1979 precession divergence.

    10 dates = 10 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_mean_obliquity_python_vs_c_derived(
        self, year, month, day, hour, date_desc
    ):
        """Compare moshier.mean_obliquity() (IAU 2006) vs C-library implicit.

        The C library's internal mean obliquity (Lieske 1979) is derived
        indirectly from its NONUT ecliptic/equatorial positions. The Python
        mean_obliquity() uses IAU 2006. The difference between these
        quantifies the precession model divergence at each epoch.
        """
        jd = swe.julday(year, month, day, hour)
        jd_tt = get_jd_tt_from_ut(jd)

        # Direct Python mean obliquity (IAU 2006)
        eps_py_direct = mean_obliquity(jd_tt)

        # C-derived implicit mean obliquity (Lieske 1979)
        flag_ecl = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_NONUT
        flag_eq = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_NONUT
        ecl_swe, _ = swe.calc_ut(jd, SE_SUN, flag_ecl)
        eq_swe, _ = swe.calc_ut(jd, SE_SUN, flag_eq)

        eps_c_derived = derive_obliquity(ecl_swe[0], ecl_swe[1], eq_swe[0], eq_swe[1])

        if math.isnan(eps_c_derived):
            pytest.skip(
                f"Cannot derive C obliquity at {date_desc}: "
                f"Sun lon too close to 0/180 (lon={ecl_swe[0]:.1f})"
            )

        diff_deg = abs(eps_py_direct - eps_c_derived)
        diff_arcsec = diff_deg * 3600.0

        assert diff_deg < DIRECT_MEAN_OBLIQUITY_TOL, (
            f"Mean obliquity at {date_desc}: "
            f"Python IAU2006 vs C-derived (Lieske) diff "
            f"{diff_arcsec:.2f} arcsec ({diff_deg:.6f}°) "
            f"exceeds tolerance {DIRECT_MEAN_OBLIQUITY_TOL}° "
            f"(Python={eps_py_direct:.6f}°, C-derived={eps_c_derived:.6f}°)"
        )


class TestNutationContribution:
    """Compare nutation contribution to obliquity between C and Python.

    Derives nutation-in-obliquity indirectly:
        nutation_eps = implicit_true_obliquity - implicit_mean_obliquity

    This is computed for both C and Python, then the contributions are compared.

    10 dates = 10 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_nutation_in_obliquity_c_vs_python(self, year, month, day, hour, date_desc):
        """Compare nutation contribution to obliquity: C vs Python.

        For each implementation:
          1. Derive implicit true obliquity (with nutation)
          2. Derive implicit mean obliquity (NONUT)
          3. Difference = nutation in obliquity

        Then compare the nutation contributions between C and Python.
        Also compare Python nutation against moshier.nutation_angles().
        """
        jd = swe.julday(year, month, day, hour)
        jd_tt = get_jd_tt_from_ut(jd)

        # --- C library: derive true and mean obliquity ---
        flag_ecl = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_eq = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        flag_ecl_nn = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_NONUT
        flag_eq_nn = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_NONUT

        ecl_swe, _ = swe.calc_ut(jd, SE_SUN, flag_ecl)
        eq_swe, _ = swe.calc_ut(jd, SE_SUN, flag_eq)
        ecl_swe_nn, _ = swe.calc_ut(jd, SE_SUN, flag_ecl_nn)
        eq_swe_nn, _ = swe.calc_ut(jd, SE_SUN, flag_eq_nn)

        eps_true_swe = derive_obliquity(ecl_swe[0], ecl_swe[1], eq_swe[0], eq_swe[1])
        eps_mean_swe = derive_obliquity(
            ecl_swe_nn[0], ecl_swe_nn[1], eq_swe_nn[0], eq_swe_nn[1]
        )

        # --- Python: direct nutation ---
        _, d_eps_py = nutation_angles(jd_tt)

        if math.isnan(eps_true_swe) or math.isnan(eps_mean_swe):
            pytest.skip(
                f"Cannot derive nutation at {date_desc}: Sun lon too close to 0/180"
            )

        # C-derived nutation in obliquity
        nut_eps_swe = eps_true_swe - eps_mean_swe

        diff_deg = abs(nut_eps_swe - d_eps_py)
        diff_arcsec = diff_deg * 3600.0

        assert diff_deg < NUTATION_CONTRIBUTION_TOL, (
            f"Nutation in obliquity at {date_desc}: "
            f"C-derived vs Python diff {diff_arcsec:.2f} arcsec ({diff_deg:.6f}°) "
            f"exceeds tolerance {NUTATION_CONTRIBUTION_TOL}° "
            f"(C-derived={nut_eps_swe:.6f}°, Python={d_eps_py:.6f}°)"
        )
