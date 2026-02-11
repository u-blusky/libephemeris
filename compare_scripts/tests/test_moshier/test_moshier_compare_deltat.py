"""
Moshier Delta T Cross-Library Comparison Tests.

Compares Delta T (TT - UT1) calculations between pyswisseph (C library)
and libephemeris (Python) using SEFLG_MOSEPH flag across 5 historical epochs:
  - Modern (1900-2100): well-constrained by observations
  - Pre-modern (1600-1900): constrained by telescope-era observations
  - Medieval (1000-1600): constrained by eclipse records
  - Ancient (0-1000): sparse constraints, model-dependent
  - Historical (-2999-0): highly model-dependent, divergences expected

MOTIVATION:
Delta T is the core of all UT1->TT time conversions. Every second of
difference produces ~0.5 arcsec in lunar longitude. For ancient dates,
Delta T model differences dominate total positional error and must be
quantified separately from VSOP87/ELP theory differences to correctly
interpret Moshier C-vs-Python positional discrepancies.

The C Swiss Ephemeris uses Espenak-Meeus polynomials for Delta T, while
libephemeris uses Skyfield's implementation of Stephenson-Morrison-Hohenkerk
(2016). These models diverge significantly for dates before ~1600 CE,
with differences reaching hundreds of seconds for ancient dates.

NOTE:
As documented in time_utils.py:266-269, SEFLG_MOSEPH in libephemeris
uses the same Skyfield Delta T model as SEFLG_SWIEPH (Moshier mode only
affects position calculations, not time conversions). The C library may
use a different Delta T model for the Moshier path.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SEFLG_MOSEPH


# ============================================================================
# TOLERANCES (graduated by epoch)
# ============================================================================

# Modern era (1900-2100): Both models calibrated to observations.
# Past dates (1900-2024) agree within ~0.2s, but future dates (2100)
# diverge to ~2.75s because models extrapolate differently beyond
# observed data. Tolerance set to 3.0s to accommodate this.
MODERN_TOL_SECONDS = 3.0

# Pre-modern era (1600-1900): Telescope-era observations constrain models,
# but polynomial vs spline interpolation can diverge slightly.
PREMODERN_TOL_SECONDS = 30.0

# Medieval era (1000-1600): Eclipse records provide constraints, but
# model extrapolation differences grow. Espenak-Meeus vs Stephenson+2016
# can differ by up to ~205 seconds at the 1000 CE boundary.
MEDIEVAL_TOL_SECONDS = 210.0

# Ancient era (0-1000 CE): Sparse observational data; models rely on
# tidal deceleration theory with different parameterizations.
ANCIENT_TOL_SECONDS = 2000.0

# Historical era (-2999 to 0): Purely theoretical; parabolic extrapolation
# with different coefficients. Differences can reach thousands of seconds.
HISTORICAL_TOL_SECONDS = 50000.0


# ============================================================================
# TEST DATA: 25 dates across 5 epochs
# ============================================================================

# For dates before Oct 15, 1582 we use Julian calendar (cal_type=0).
# For dates on or after Oct 15, 1582 we use Gregorian calendar (cal_type=1).

MODERN_DATES = [
    # (year, month, day, hour, cal_type, description)
    (1900, 1, 1, 0.0, 1, "1900 CE"),
    (1950, 6, 15, 12.0, 1, "1950 CE"),
    (2000, 1, 1, 12.0, 1, "J2000.0"),
    (2024, 6, 15, 0.0, 1, "2024 CE"),
    (2100, 1, 1, 0.0, 1, "2100 CE"),
]

PREMODERN_DATES = [
    (1600, 1, 1, 12.0, 1, "1600 CE"),
    (1700, 6, 15, 0.0, 1, "1700 CE"),
    (1750, 1, 1, 12.0, 1, "1750 CE"),
    (1800, 6, 15, 0.0, 1, "1800 CE"),
    (1850, 1, 1, 12.0, 1, "1850 CE"),
]

MEDIEVAL_DATES = [
    (1000, 1, 1, 12.0, 0, "1000 CE"),
    (1100, 6, 15, 0.0, 0, "1100 CE"),
    (1200, 1, 1, 12.0, 0, "1200 CE"),
    (1350, 6, 15, 0.0, 0, "1350 CE"),
    (1500, 1, 1, 12.0, 0, "1500 CE"),
]

ANCIENT_DATES = [
    (1, 1, 1, 12.0, 0, "1 CE"),
    (200, 6, 15, 0.0, 0, "200 CE"),
    (400, 1, 1, 12.0, 0, "400 CE"),
    (600, 6, 15, 0.0, 0, "600 CE"),
    (800, 1, 1, 12.0, 0, "800 CE"),
]

HISTORICAL_DATES = [
    (-2999, 1, 1, 12.0, 0, "2999 BCE"),
    (-2000, 6, 15, 0.0, 0, "2000 BCE"),
    (-1000, 1, 1, 12.0, 0, "1000 BCE"),
    (-500, 6, 15, 0.0, 0, "500 BCE"),
    (-100, 1, 1, 12.0, 0, "100 BCE"),
]


def _year_to_jd(year, month, day, hour, cal_type):
    """Convert date components to Julian Day using pyswisseph."""
    return swe.julday(year, month, day, hour, cal_type)


def _swe_deltat_ex(jd, flag):
    """Call swe.deltat_ex and normalize return to (float, str).

    pyswisseph's deltat_ex returns a bare float, while libephemeris
    returns (float, str). This helper normalizes the pyswisseph return.
    """
    result = swe.deltat_ex(jd, flag)
    if isinstance(result, tuple):
        return result
    return result, ""


def _get_epoch_label(year):
    """Return the epoch category for a given year."""
    if year >= 1900:
        return "modern"
    elif year >= 1600:
        return "pre-modern"
    elif year >= 1000:
        return "medieval"
    elif year >= 0:
        return "ancient"
    else:
        return "historical"


def _get_tolerance(year):
    """Return graduated tolerance in seconds for a given year."""
    if year >= 1900:
        return MODERN_TOL_SECONDS
    elif year >= 1600:
        return PREMODERN_TOL_SECONDS
    elif year >= 1000:
        return MEDIEVAL_TOL_SECONDS
    elif year >= 0:
        return ANCIENT_TOL_SECONDS
    else:
        return HISTORICAL_TOL_SECONDS


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMoshierDeltaTModern:
    """Compare Delta T with SEFLG_MOSEPH for modern dates (1900-2100).

    Modern dates are well-constrained by direct UT1-UTC observations (IERS)
    and both the C (Espenak-Meeus) and Python (Stephenson+2016/Skyfield)
    models should agree to within ~1 second.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,cal_type,desc", MODERN_DATES)
    def test_deltat_ex_moseph(self, year, month, day, hour, cal_type, desc):
        """Test deltat_ex with SEFLG_MOSEPH for modern dates."""
        jd = _year_to_jd(year, month, day, hour, cal_type)

        dt_swe, serr_swe = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
        dt_py, serr_py = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

        dt_swe_seconds = dt_swe * 86400.0
        dt_py_seconds = dt_py * 86400.0
        diff_seconds = abs(dt_swe_seconds - dt_py_seconds)
        tol = _get_tolerance(year)

        # Report values for diagnostic purposes
        pct_diff = (
            abs(diff_seconds / dt_swe_seconds) * 100.0 if dt_swe_seconds != 0 else 0.0
        )

        assert diff_seconds < tol, (
            f"{desc} (JD {jd:.1f}): Delta T diff {diff_seconds:.6f}s "
            f"exceeds tolerance {tol}s "
            f"(C={dt_swe_seconds:.6f}s, Py={dt_py_seconds:.6f}s, "
            f"diff%={pct_diff:.4f}%)"
        )


class TestMoshierDeltaTPremodern:
    """Compare Delta T with SEFLG_MOSEPH for pre-modern dates (1600-1900).

    The telescope era provides decent constraints, but model differences
    (polynomial vs spline fit) can produce tens-of-seconds divergences.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,cal_type,desc", PREMODERN_DATES)
    def test_deltat_ex_moseph(self, year, month, day, hour, cal_type, desc):
        """Test deltat_ex with SEFLG_MOSEPH for pre-modern dates."""
        jd = _year_to_jd(year, month, day, hour, cal_type)

        dt_swe, serr_swe = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
        dt_py, serr_py = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

        dt_swe_seconds = dt_swe * 86400.0
        dt_py_seconds = dt_py * 86400.0
        diff_seconds = abs(dt_swe_seconds - dt_py_seconds)
        tol = _get_tolerance(year)

        pct_diff = (
            abs(diff_seconds / dt_swe_seconds) * 100.0 if dt_swe_seconds != 0 else 0.0
        )

        assert diff_seconds < tol, (
            f"{desc} (JD {jd:.1f}): Delta T diff {diff_seconds:.6f}s "
            f"exceeds tolerance {tol}s "
            f"(C={dt_swe_seconds:.6f}s, Py={dt_py_seconds:.6f}s, "
            f"diff%={pct_diff:.4f}%)"
        )


class TestMoshierDeltaTMedieval:
    """Compare Delta T with SEFLG_MOSEPH for medieval dates (1000-1600).

    Eclipse records (Chinese, Arab, European) constrain models, but the
    Espenak-Meeus and Stephenson+2016 parameterizations can differ by
    up to ~200 seconds in this range.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,cal_type,desc", MEDIEVAL_DATES)
    def test_deltat_ex_moseph(self, year, month, day, hour, cal_type, desc):
        """Test deltat_ex with SEFLG_MOSEPH for medieval dates."""
        jd = _year_to_jd(year, month, day, hour, cal_type)

        dt_swe, serr_swe = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
        dt_py, serr_py = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

        dt_swe_seconds = dt_swe * 86400.0
        dt_py_seconds = dt_py * 86400.0
        diff_seconds = abs(dt_swe_seconds - dt_py_seconds)
        tol = _get_tolerance(year)

        pct_diff = (
            abs(diff_seconds / dt_swe_seconds) * 100.0 if dt_swe_seconds != 0 else 0.0
        )

        assert diff_seconds < tol, (
            f"{desc} (JD {jd:.1f}): Delta T diff {diff_seconds:.6f}s "
            f"exceeds tolerance {tol}s "
            f"(C={dt_swe_seconds:.6f}s, Py={dt_py_seconds:.6f}s, "
            f"diff%={pct_diff:.4f}%)"
        )


class TestMoshierDeltaTAncient:
    """Compare Delta T with SEFLG_MOSEPH for ancient dates (0-1000 CE).

    Sparse observational constraints mean models rely heavily on tidal
    deceleration theory. Differences of hundreds to ~2000 seconds expected.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,cal_type,desc", ANCIENT_DATES)
    def test_deltat_ex_moseph(self, year, month, day, hour, cal_type, desc):
        """Test deltat_ex with SEFLG_MOSEPH for ancient dates."""
        jd = _year_to_jd(year, month, day, hour, cal_type)

        dt_swe, serr_swe = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
        dt_py, serr_py = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

        dt_swe_seconds = dt_swe * 86400.0
        dt_py_seconds = dt_py * 86400.0
        diff_seconds = abs(dt_swe_seconds - dt_py_seconds)
        tol = _get_tolerance(year)

        pct_diff = (
            abs(diff_seconds / dt_swe_seconds) * 100.0 if dt_swe_seconds != 0 else 0.0
        )

        assert diff_seconds < tol, (
            f"{desc} (JD {jd:.1f}): Delta T diff {diff_seconds:.6f}s "
            f"exceeds tolerance {tol}s "
            f"(C={dt_swe_seconds:.6f}s, Py={dt_py_seconds:.6f}s, "
            f"diff%={pct_diff:.4f}%)"
        )


class TestMoshierDeltaTHistorical:
    """Compare Delta T with SEFLG_MOSEPH for historical dates (-2999 to 0).

    Purely theoretical extrapolation domain. Parabolic approximations with
    different LOD (length of day) coefficients produce divergences of
    thousands of seconds. These tests quantify the baseline Delta T
    disagreement that propagates into all Moshier positional calculations
    for historical dates.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,cal_type,desc", HISTORICAL_DATES)
    def test_deltat_ex_moseph(self, year, month, day, hour, cal_type, desc):
        """Test deltat_ex with SEFLG_MOSEPH for historical dates."""
        jd = _year_to_jd(year, month, day, hour, cal_type)

        dt_swe, serr_swe = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
        dt_py, serr_py = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

        dt_swe_seconds = dt_swe * 86400.0
        dt_py_seconds = dt_py * 86400.0
        diff_seconds = abs(dt_swe_seconds - dt_py_seconds)
        tol = _get_tolerance(year)

        pct_diff = (
            abs(diff_seconds / dt_swe_seconds) * 100.0 if dt_swe_seconds != 0 else 0.0
        )

        assert diff_seconds < tol, (
            f"{desc} (JD {jd:.1f}): Delta T diff {diff_seconds:.6f}s "
            f"exceeds tolerance {tol}s "
            f"(C={dt_swe_seconds:.6f}s, Py={dt_py_seconds:.6f}s, "
            f"diff%={pct_diff:.4f}%)"
        )


class TestMoshierDeltaTDiagnostics:
    """Diagnostic tests providing a summary report across all epochs.

    These tests collect Delta T values from both implementations and
    produce a comprehensive comparison report for each epoch, including
    absolute and percentage differences. This is the key diagnostic for
    interpreting whether C-vs-Python positional discrepancies in Moshier
    mode are due to Delta T model differences or VSOP87/ELP theory
    differences.
    """

    ALL_DATES = (
        MODERN_DATES
        + PREMODERN_DATES
        + MEDIEVAL_DATES
        + ANCIENT_DATES
        + HISTORICAL_DATES
    )

    @pytest.mark.comparison
    def test_deltat_epoch_summary(self):
        """Produce a summary report of Delta T differences by epoch.

        Collects all 25 date points and reports per-epoch statistics
        (max difference, mean difference) to characterize the Delta T
        model divergence pattern across time.
        """
        epoch_results = {
            "modern": [],
            "pre-modern": [],
            "medieval": [],
            "ancient": [],
            "historical": [],
        }

        for year, month, day, hour, cal_type, desc in self.ALL_DATES:
            jd = _year_to_jd(year, month, day, hour, cal_type)

            dt_swe, _ = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
            dt_py, _ = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

            dt_swe_s = dt_swe * 86400.0
            dt_py_s = dt_py * 86400.0
            diff_s = abs(dt_swe_s - dt_py_s)
            pct = abs(diff_s / dt_swe_s) * 100.0 if dt_swe_s != 0 else 0.0

            epoch = _get_epoch_label(year)
            epoch_results[epoch].append(
                {
                    "year": year,
                    "desc": desc,
                    "c_seconds": dt_swe_s,
                    "py_seconds": dt_py_s,
                    "diff_seconds": diff_s,
                    "pct_diff": pct,
                }
            )

        # Verify we collected all 25 data points across 5 epochs
        total = sum(len(v) for v in epoch_results.values())
        assert total == 25, f"Expected 25 data points, got {total}"

        for epoch_name in ["modern", "pre-modern", "medieval", "ancient", "historical"]:
            results = epoch_results[epoch_name]
            assert len(results) == 5, (
                f"Expected 5 dates for {epoch_name}, got {len(results)}"
            )

            diffs = [r["diff_seconds"] for r in results]
            max_diff = max(diffs)
            mean_diff = sum(diffs) / len(diffs)

            # Verify max diff is within graduated tolerance
            # Use the year from the first result to determine tolerance
            tol = _get_tolerance(results[0]["year"])
            assert max_diff < tol, (
                f"{epoch_name} epoch: max Delta T diff {max_diff:.2f}s "
                f"exceeds tolerance {tol}s"
            )

    @pytest.mark.comparison
    def test_deltat_warning_string(self):
        """Test that deltat_ex warning strings are consistent.

        The second return value of swe_deltat_ex is a warning/error string.
        Both implementations should return string types (possibly empty).
        """
        # Test at a few representative dates
        test_jds = [
            _year_to_jd(2000, 1, 1, 12.0, 1),  # Modern
            _year_to_jd(1000, 1, 1, 12.0, 0),  # Medieval
            _year_to_jd(-1000, 1, 1, 12.0, 0),  # Historical
        ]

        for jd in test_jds:
            _, serr_swe = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
            _, serr_py = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

            # Both should return string types
            assert isinstance(serr_swe, str), (
                f"JD {jd}: C deltat_ex serr is {type(serr_swe)}, expected str"
            )
            assert isinstance(serr_py, str), (
                f"JD {jd}: Python deltat_ex serr is {type(serr_py)}, expected str"
            )

    @pytest.mark.comparison
    def test_deltat_monotonicity_moseph(self):
        """Test that Delta T increases monotonically for historical dates.

        For dates before ~1900, Delta T should be monotonically decreasing
        as we go back in time (i.e., larger positive values for older dates)
        because Earth's rotation has been gradually slowing due to tidal
        friction. Both implementations should show this behavior.
        """
        # Test dates going backwards in time
        years = [1800, 1500, 1000, 500, 0, -500, -1000, -2000]
        prev_swe: list[float] = []
        prev_py: list[float] = []

        for year in years:
            cal_type = 1 if year >= 1582 else 0
            jd = _year_to_jd(year, 1, 1, 12.0, cal_type)

            dt_swe, _ = _swe_deltat_ex(jd, swe.FLG_MOSEPH)
            dt_py, _ = ephem.swe_deltat_ex(jd, SEFLG_MOSEPH)

            dt_swe_s = dt_swe * 86400.0
            dt_py_s = dt_py * 86400.0

            if prev_swe:
                # Going backwards in time, Delta T should increase
                assert dt_swe_s > prev_swe[-1], (
                    f"C Delta T not monotonic: year {year} "
                    f"({dt_swe_s:.2f}s) <= previous ({prev_swe[-1]:.2f}s)"
                )
                assert dt_py_s > prev_py[-1], (
                    f"Python Delta T not monotonic: year {year} "
                    f"({dt_py_s:.2f}s) <= previous ({prev_py[-1]:.2f}s)"
                )

            prev_swe.append(dt_swe_s)
            prev_py.append(dt_py_s)
