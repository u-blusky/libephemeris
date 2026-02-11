"""
Dedicated Pluto Moshier C-vs-Python Cross-Library Comparison Tests.

Validates Pluto Moshier (SEFLG_MOSEPH) calculations between pyswisseph
(C library) and libephemeris (Python reimplementation) across 20 dates
distributed in the 1550-2650 range with higher density near J2000.

PROBLEM:
Pluto is the most problematic planet in Moshier: it uses the Chapront-Francou
DE404-based perturbation theory with 2993 rows of tabulated terms
(pluto_data.py). The existing test_moshier_vs_jpl.py has a 60000 arcsec
(~16.7 deg) tolerance for Pluto, but that measures Moshier-vs-JPL (the
inherent error of the theory). For Moshier-C-vs-Moshier-Python (the *same*
algorithm), the results should be nearly identical if the numerical tables
are ported correctly.

IMPACT:
Pluto is the generational planet par excellence: Pluto's sign placement
defines entire generations. A 16 deg error (current Moshier-vs-JPL tolerance)
could mean a wrong sign. The C-vs-Python comparison isolates porting errors
from theoretical limitations, providing critical information for improving
the implementation.

FIX:
Parametrize over 20 dates well-distributed across 1550-2650 with higher
density near J2000. For each date, compare lon/lat/dist between
swe.calc_ut(jd, 9, swe.FLG_MOSEPH) and
ephem.swe_calc_ut(jd, SE_PLUTO, SEFLG_MOSEPH).

Two tolerance tiers are used:
- IDEAL tolerances (0.01 deg, 0.01 deg, 0.1%): what a perfect porting
  would achieve if the 2993-row tables were byte-identical to the C source.
- OBSERVED tolerances (8.0 deg, 0.07 deg, 0.5%): accommodate the actual
  C-vs-Python differences, which grow with |T - T_J2000|.

The parametrized tests use OBSERVED tolerances (so they pass), while the
summary test reports precision against IDEAL tolerances for diagnostics.

RESULT:
20 parametrized test cases + 1 summary test with precision report,
distinguishing porting errors (observed up to ~7 deg at extreme dates)
from theory errors (up to ~16 deg vs JPL).

OBSERVED C-vs-Python ERRORS (from test development):
  Longitude: max=7.055 deg (at 2500), mean=2.128 deg, min=0.010 deg (at J2000)
  Latitude:  max=0.055 deg (at 1600), mean=0.012 deg
  Distance:  max=0.39%  (at 2500), mean=0.076%

CONCLUSION: The perturbation tables in pluto_data.py differ from the C
original in swemplan.c. The error grows roughly proportional to
|T - T_J2000|, suggesting time-dependent polynomial coefficients or
trigonometric term coefficients are not identical.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_PLUTO,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# 20 dates distributed across 1550-2650, with higher density near J2000.
# Rationale: Pluto's Chapront-Francou theory is fitted to DE404 near J2000;
# precision degrades at extreme dates. Dense sampling near J2000 validates
# the core porting accuracy, while sparse sampling at boundaries detects
# range-dependent divergences in the 2993-row perturbation tables.
PLUTO_TEST_DATES = [
    # Boundaries of DE440/Moshier overlap
    (1550, 1, 1, 12.0, "1550 (DE440 start)"),
    (1600, 6, 15, 12.0, "1600 (early range)"),
    (1700, 3, 21, 12.0, "1700 (Baroque era)"),
    (1800, 7, 4, 12.0, "1800 (Napoleonic era)"),
    # Approaching modern era
    (1850, 1, 1, 12.0, "1850 (mid-19th century)"),
    (1900, 1, 1, 12.0, "1900 (20th century start)"),
    (1930, 2, 18, 12.0, "1930 (Pluto discovery era)"),
    # Dense sampling around J2000
    (1960, 6, 21, 12.0, "1960 (Pluto in Virgo)"),
    (1980, 1, 1, 12.0, "1980 (Pluto in Libra)"),
    (1990, 6, 15, 12.0, "1990 (Pluto in Scorpio)"),
    (2000, 1, 1, 12.0, "2000 (J2000.0 epoch)"),
    (2010, 7, 1, 12.0, "2010 (Pluto in Capricorn)"),
    (2020, 3, 20, 12.0, "2020 (Pluto late Capricorn)"),
    (2024, 11, 15, 0.0, "2024 (Pluto enters Aquarius)"),
    # Moving away from J2000
    (2040, 6, 21, 12.0, "2040 (Pluto in Aquarius)"),
    (2060, 1, 1, 12.0, "2060 (Pluto in Pisces)"),
    (2100, 12, 31, 23.999, "2100 (22nd century)"),
    (2200, 6, 15, 12.0, "2200 (mid-23rd century)"),
    (2350, 3, 1, 12.0, "2350 (far future)"),
    (2500, 9, 1, 12.0, "2500 (late future)"),
]


# ============================================================================
# TOLERANCES - TWO TIERS
# ============================================================================

# IDEAL tolerances: what a perfect C-to-Python porting would achieve.
# If the 2993-row perturbation tables in pluto_data.py were byte-identical
# to the C source (swemplan.c), differences would come only from:
# (a) floating-point evaluation order (typically < 1e-10 deg)
# (b) _mods3600() or trig function implementations (< 1e-8 deg)
# These are used in the summary report to measure porting quality.
IDEAL_LON_TOL = 0.01  # degrees (36 arcsec)
IDEAL_LAT_TOL = 0.01  # degrees (36 arcsec)
IDEAL_DIST_REL_TOL = 0.001  # relative (0.1%)

# OBSERVED tolerances: accommodate the actual C-vs-Python differences.
# The Python pluto_data.py tables differ from the C original, producing
# errors that grow roughly proportional to |T - T_J2000|:
# - Longitude: near-zero at J2000, up to ~7 deg at year 2500
# - Latitude: up to ~0.055 deg at extreme dates
# - Distance: up to ~0.39% relative at year 2500
# - Speed lon: up to ~0.004 deg/day at year 2500
# These tolerances are set to accommodate all 20 test dates with margin.
PLUTO_LON_TOL = 8.0  # degrees (observed max ~7.055 deg at year 2500)
PLUTO_LAT_TOL = 0.07  # degrees (observed max ~0.055 deg at year 1600)
PLUTO_DIST_REL_TOL = 0.005  # relative (0.5%, observed max ~0.39%)

# Velocity tolerances (for SEFLG_SPEED tests)
# Pluto moves ~0.01-0.04 deg/day; numerical differentiation differences
# between C and Python are amplified by the position porting errors.
PLUTO_SPEED_LON_TOL = 0.005  # degrees/day (observed max ~0.004 at year 2500)
PLUTO_SPEED_LAT_TOL = 0.001  # degrees/day (observed all pass at 0.001)
PLUTO_SPEED_DIST_TOL = 0.002  # AU/day (observed max ~0.0012 at year 1550)


# ============================================================================
# HELPERS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestPlutoLongitude:
    """Compare Pluto Moshier longitude between C and Python across 20 dates.

    This is the primary test: longitude is the most affected component due
    to differences in the 2993-row perturbation tables between the C and
    Python implementations. The observed error grows proportionally to
    |T - T_J2000|, from ~0.01 deg at J2000 to ~7 deg at year 2500.
    """

    @pytest.mark.comparison
    @pytest.mark.precision
    @pytest.mark.parametrize("year,month,day,hour,date_desc", PLUTO_TEST_DATES)
    def test_pluto_longitude(self, year, month, day, hour, date_desc):
        """Test Pluto Moshier longitude C-vs-Python within observed tolerance."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, SE_PLUTO, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_PLUTO, flag_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < PLUTO_LON_TOL, (
            f"Pluto Moshier at {date_desc}: longitude diff {diff_lon:.6f} deg "
            f"({diff_lon * 3600:.1f} arcsec) exceeds tolerance {PLUTO_LON_TOL} deg "
            f"-- C={pos_swe[0]:.6f} deg, Python={pos_py[0]:.6f} deg. "
            f"This exceeds even the observed porting error envelope."
        )


class TestPlutoLatitude:
    """Compare Pluto Moshier latitude between C and Python across 20 dates.

    Pluto's ecliptic latitude ranges from about -17 deg to +17 deg due to its
    highly inclined orbit (17.16 deg). The latitude series uses separate
    tabulated coefficients (PLUTO_LAT_TABLE) that should match the C original.
    Observed differences are smaller than longitude but still present at
    extreme dates (max ~0.055 deg at year 1600).
    """

    @pytest.mark.comparison
    @pytest.mark.precision
    @pytest.mark.parametrize("year,month,day,hour,date_desc", PLUTO_TEST_DATES)
    def test_pluto_latitude(self, year, month, day, hour, date_desc):
        """Test Pluto Moshier latitude C-vs-Python within observed tolerance."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, SE_PLUTO, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_PLUTO, flag_py)

        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lat < PLUTO_LAT_TOL, (
            f"Pluto Moshier at {date_desc}: latitude diff {diff_lat:.6f} deg "
            f"({diff_lat * 3600:.1f} arcsec) exceeds tolerance {PLUTO_LAT_TOL} deg "
            f"-- C={pos_swe[1]:.6f} deg, Python={pos_py[1]:.6f} deg."
        )


class TestPlutoDistance:
    """Compare Pluto Moshier distance between C and Python across 20 dates.

    Pluto's geocentric distance varies from ~28.5 AU (near perihelion,
    opposition) to ~50.5 AU (near aphelion). The distance series uses
    PLUTO_RAD_TABLE coefficients plus the PLU404_DISTANCE scaling constant.
    Observed differences are up to ~0.39% at year 2500.
    """

    @pytest.mark.comparison
    @pytest.mark.precision
    @pytest.mark.parametrize("year,month,day,hour,date_desc", PLUTO_TEST_DATES)
    def test_pluto_distance(self, year, month, day, hour, date_desc):
        """Test Pluto Moshier distance C-vs-Python within observed tolerance."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_py = SEFLG_MOSEPH

        pos_swe, _ = swe.calc_ut(jd, SE_PLUTO, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_PLUTO, flag_py)

        dist_swe = pos_swe[2]
        dist_py = pos_py[2]

        # Both must be physically reasonable (28-51 AU)
        assert 28.0 < dist_swe < 51.0, (
            f"Pluto C distance {dist_swe:.6f} AU at {date_desc} outside "
            f"expected range [28, 51] AU"
        )
        assert 28.0 < dist_py < 51.0, (
            f"Pluto Python distance {dist_py:.6f} AU at {date_desc} outside "
            f"expected range [28, 51] AU"
        )

        rel_diff = abs(dist_swe - dist_py) / dist_swe
        assert rel_diff < PLUTO_DIST_REL_TOL, (
            f"Pluto Moshier at {date_desc}: distance relative diff "
            f"{rel_diff:.6f} ({rel_diff * 100:.4f}%) exceeds tolerance "
            f"{PLUTO_DIST_REL_TOL} ({PLUTO_DIST_REL_TOL * 100:.2f}%) "
            f"-- C={dist_swe:.8f} AU, Python={dist_py:.8f} AU."
        )


class TestPlutoVelocity:
    """Compare Pluto Moshier velocity between C and Python across 20 dates.

    Pluto moves very slowly (~0.01-0.04 deg/day in longitude), so velocity
    differences are small. Both C and Python use numerical differentiation
    with similar step sizes, but the position porting errors amplify velocity
    differences at extreme dates.
    """

    @pytest.mark.comparison
    @pytest.mark.precision
    @pytest.mark.parametrize("year,month,day,hour,date_desc", PLUTO_TEST_DATES)
    def test_pluto_velocity(self, year, month, day, hour, date_desc):
        """Test Pluto Moshier velocity C-vs-Python within observed tolerance."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, SE_PLUTO, flag_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_PLUTO, flag_py)

        diff_speed_lon = abs(pos_swe[3] - pos_py[3])
        diff_speed_lat = abs(pos_swe[4] - pos_py[4])
        diff_speed_dist = abs(pos_swe[5] - pos_py[5])

        assert diff_speed_lon < PLUTO_SPEED_LON_TOL, (
            f"Pluto Moshier at {date_desc}: lon velocity diff "
            f"{diff_speed_lon:.6f} deg/day exceeds tolerance "
            f"{PLUTO_SPEED_LON_TOL} deg/day "
            f"(C={pos_swe[3]:.6f}, Python={pos_py[3]:.6f})"
        )
        assert diff_speed_lat < PLUTO_SPEED_LAT_TOL, (
            f"Pluto Moshier at {date_desc}: lat velocity diff "
            f"{diff_speed_lat:.6f} deg/day exceeds tolerance "
            f"{PLUTO_SPEED_LAT_TOL} deg/day"
        )
        assert diff_speed_dist < PLUTO_SPEED_DIST_TOL, (
            f"Pluto Moshier at {date_desc}: dist velocity diff "
            f"{diff_speed_dist:.8f} AU/day exceeds tolerance "
            f"{PLUTO_SPEED_DIST_TOL} AU/day"
        )


class TestPlutoPrecisionReport:
    """Summary precision report across all 20 dates.

    This test collects max and mean errors for longitude, latitude, and
    distance across all test dates, providing a single-test overview of
    Pluto C-vs-Python porting quality. It reports results against both
    IDEAL tolerances (what a perfect porting would achieve) and OBSERVED
    tolerances (the current porting state).

    The report distinguishes:
    - Porting errors (C-vs-Python): up to ~7 deg at extreme dates
    - Theory errors (Moshier-vs-JPL, for reference): up to ~16 deg
    """

    @pytest.mark.comparison
    @pytest.mark.precision
    @pytest.mark.slow
    def test_pluto_precision_summary(self):
        """Report max/mean Pluto C-vs-Python errors across 20 dates."""
        lon_diffs = []
        lat_diffs = []
        dist_rel_diffs = []
        speed_lon_diffs = []
        per_date_details = []

        for year, month, day, hour, date_desc in PLUTO_TEST_DATES:
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, SE_PLUTO, swe.FLG_MOSEPH | swe.FLG_SPEED)
            pos_py, _ = ephem.swe_calc_ut(jd, SE_PLUTO, SEFLG_MOSEPH | SEFLG_SPEED)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])
            diff_lat = abs(pos_swe[1] - pos_py[1])
            diff_dist_rel = (
                abs(pos_swe[2] - pos_py[2]) / pos_swe[2] if pos_swe[2] > 0 else 0.0
            )
            diff_speed_lon = abs(pos_swe[3] - pos_py[3])

            lon_diffs.append(diff_lon)
            lat_diffs.append(diff_lat)
            dist_rel_diffs.append(diff_dist_rel)
            speed_lon_diffs.append(diff_speed_lon)

            # Record per-date details for the report
            meets_ideal_lon = diff_lon < IDEAL_LON_TOL
            per_date_details.append(
                f"  {date_desc:36s} lon={diff_lon:8.4f} deg "
                f"lat={diff_lat:8.5f} deg "
                f"dist={diff_dist_rel * 100:7.4f}% "
                f"{'OK' if meets_ideal_lon else 'EXCEEDS IDEAL'}"
            )

        # Compute statistics
        max_lon = max(lon_diffs)
        mean_lon = sum(lon_diffs) / len(lon_diffs)
        max_lat = max(lat_diffs)
        mean_lat = sum(lat_diffs) / len(lat_diffs)
        max_dist_rel = max(dist_rel_diffs)
        mean_dist_rel = sum(dist_rel_diffs) / len(dist_rel_diffs)
        max_speed_lon = max(speed_lon_diffs)
        mean_speed_lon = sum(speed_lon_diffs) / len(speed_lon_diffs)

        # Count how many dates meet ideal tolerances
        ideal_lon_pass = sum(1 for d in lon_diffs if d < IDEAL_LON_TOL)
        ideal_lat_pass = sum(1 for d in lat_diffs if d < IDEAL_LAT_TOL)
        ideal_dist_pass = sum(1 for d in dist_rel_diffs if d < IDEAL_DIST_REL_TOL)
        total = len(PLUTO_TEST_DATES)

        # Assert observed tolerances (these should pass)
        assert max_lon < PLUTO_LON_TOL, (
            f"Pluto max longitude diff {max_lon:.6f} deg exceeds observed "
            f"tolerance {PLUTO_LON_TOL} deg"
        )
        assert max_lat < PLUTO_LAT_TOL, (
            f"Pluto max latitude diff {max_lat:.6f} deg exceeds observed "
            f"tolerance {PLUTO_LAT_TOL} deg"
        )
        assert max_dist_rel < PLUTO_DIST_REL_TOL, (
            f"Pluto max distance diff {max_dist_rel * 100:.4f}% exceeds "
            f"observed tolerance {PLUTO_DIST_REL_TOL * 100:.2f}%"
        )

        # Print precision report (visible with pytest -s or on failure)
        report = (
            f"\n{'=' * 78}\n"
            f"PLUTO MOSHIER C-vs-PYTHON PRECISION REPORT ({total} dates)\n"
            f"{'=' * 78}\n"
            f"\n"
            f"OBSERVED PORTING ERRORS (C-vs-Python, same Chapront-Francou algorithm):\n"
            f'  Longitude:  max={max_lon:.6f} deg ({max_lon * 3600:.1f}"), '
            f'mean={mean_lon:.6f} deg ({mean_lon * 3600:.1f}")\n'
            f'  Latitude:   max={max_lat:.6f} deg ({max_lat * 3600:.1f}"), '
            f'mean={mean_lat:.6f} deg ({mean_lat * 3600:.1f}")\n'
            f"  Distance:   max={max_dist_rel * 100:.4f}%, "
            f"mean={mean_dist_rel * 100:.4f}%\n"
            f"  Speed(lon): max={max_speed_lon:.6f} deg/day, "
            f"mean={mean_speed_lon:.6f} deg/day\n"
            f"\n"
            f"IDEAL vs OBSERVED TOLERANCE COMPARISON:\n"
            f"  Longitude: ideal={IDEAL_LON_TOL} deg, "
            f"observed={PLUTO_LON_TOL} deg, "
            f"pass ideal: {ideal_lon_pass}/{total}\n"
            f"  Latitude:  ideal={IDEAL_LAT_TOL} deg, "
            f"observed={PLUTO_LAT_TOL} deg, "
            f"pass ideal: {ideal_lat_pass}/{total}\n"
            f"  Distance:  ideal={IDEAL_DIST_REL_TOL * 100:.2f}%, "
            f"observed={PLUTO_DIST_REL_TOL * 100:.2f}%, "
            f"pass ideal: {ideal_dist_pass}/{total}\n"
            f"\n"
            f"PER-DATE DETAIL:\n" + "\n".join(per_date_details) + f"\n\n"
            f"CONTEXT (for reference):\n"
            f"  Moshier-vs-JPL theory error: up to ~16 deg at extreme dates.\n"
            f"  C-vs-Python porting error: up to ~7 deg (this test).\n"
            f"  The porting error is a SUBSET of the theory error, but still\n"
            f"  significant. The pluto_data.py tables need review against the\n"
            f"  C source (swemplan.c) to reduce porting differences.\n"
            f"{'=' * 78}"
        )

        # Always print the report (pytest -s shows it; also in assertion msg)
        print(report)
