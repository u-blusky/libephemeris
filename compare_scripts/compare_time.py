"""
Time Functions Comparison Script.

Compares all time-related functions between pyswisseph and libephemeris:
- Julian Day conversions (julday, revjul)
- Delta T calculations
- UTC conversions
- Sidereal time
- Equation of time
- Local Mean Time conversions
"""

import swisseph as swe
import libephemeris as pyephem
import sys
import math

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import (
    TestStatistics,
    print_header,
    print_section,
    parse_args,
    format_status,
    STANDARD_SUBJECTS,
    HISTORICAL_SUBJECTS,
    EDGE_CASE_DATES,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class TimeTolerance:
    """Tolerance thresholds for time comparisons."""

    JULIAN_DAY = 1e-10  # JD precision
    DELTA_T = 0.01  # seconds
    SIDEREAL_TIME = 0.0001  # hours
    TIME_EQU = 0.001  # degrees


# ============================================================================
# TEST CASES
# ============================================================================

JULDAY_TEST_CASES = [
    # (year, month, day, hour, calendar_type, description)
    (2000, 1, 1, 12.0, 1, "J2000.0 Epoch"),
    (1900, 1, 1, 0.0, 1, "1900 Start"),
    (2024, 2, 29, 12.0, 1, "Leap Year 2024"),
    (1999, 12, 31, 23.999, 1, "End of 1999"),
    (2050, 6, 15, 6.5, 1, "Future Date"),
    (1582, 10, 15, 12.0, 1, "Gregorian Start"),
    (1582, 10, 4, 12.0, 0, "Julian Calendar Last Day"),
    (-4713, 1, 1, 12.0, 0, "JD Zero Point"),
    (1, 1, 1, 0.0, 0, "Year 1 AD Julian"),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_julday(year, month, day, hour, cal_type, desc, verbose=False):
    """Compare julday function."""
    try:
        jd_swe = swe.julday(year, month, day, hour, cal_type)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        jd_py = pyephem.julday(year, month, day, hour, cal_type)
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    diff = abs(jd_swe - jd_py)
    passed = diff < TimeTolerance.JULIAN_DAY

    if verbose:
        print(f"\n[julday] {desc}")
        print(f"  Input: {year}-{month:02d}-{day:02d} {hour}h (cal={cal_type})")
        print(f"  SWE: {jd_swe:.10f}")
        print(f"  PY:  {jd_py:.10f}")
        print(f"  Diff: {diff:.15f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[julday] {desc}: SWE={jd_swe:.6f} PY={jd_py:.6f} Diff={diff:.12f} {status}"
        )

    return passed, diff, None


def compare_revjul(jd, desc, verbose=False):
    """Compare revjul function."""
    try:
        y_swe, m_swe, d_swe, h_swe = swe.revjul(jd)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        y_py, m_py, d_py, h_py = pyephem.revjul(jd)
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    # Compare all components
    year_match = int(y_swe) == int(y_py)
    month_match = int(m_swe) == int(m_py)
    day_match = int(d_swe) == int(d_py)
    hour_diff = abs(h_swe - h_py)

    passed = year_match and month_match and day_match and hour_diff < 1e-8

    if verbose:
        print(f"\n[revjul] {desc}")
        print(f"  Input JD: {jd}")
        print(f"  SWE: {int(y_swe)}-{int(m_swe):02d}-{int(d_swe):02d} {h_swe:.8f}h")
        print(f"  PY:  {int(y_py)}-{int(m_py):02d}-{int(d_py):02d} {h_py:.8f}h")
        print(f"  Hour Diff: {hour_diff:.12f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[revjul] {desc}: "
            f"SWE={int(y_swe)}-{int(m_swe):02d}-{int(d_swe):02d} "
            f"PY={int(y_py)}-{int(m_py):02d}-{int(d_py):02d} "
            f"HourDiff={hour_diff:.10f} {status}"
        )

    return passed, hour_diff, None


def compare_deltat(jd, desc, verbose=False):
    """Compare deltat function."""
    try:
        dt_swe = swe.deltat(jd)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        dt_py = pyephem.deltat(jd)
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    # Delta T is in days, convert diff to seconds for tolerance check
    diff_days = abs(dt_swe - dt_py)
    diff_seconds = diff_days * 86400

    passed = diff_seconds < TimeTolerance.DELTA_T

    if verbose:
        print(f"\n[deltat] {desc}")
        print(f"  Input JD: {jd}")
        print(f"  SWE: {dt_swe:.10f} days ({dt_swe * 86400:.4f} sec)")
        print(f"  PY:  {dt_py:.10f} days ({dt_py * 86400:.4f} sec)")
        print(f"  Diff: {diff_seconds:.6f} sec {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[deltat] {desc}: "
            f"SWE={dt_swe * 86400:.4f}s PY={dt_py * 86400:.4f}s "
            f"Diff={diff_seconds:.6f}s {status}"
        )

    return passed, diff_seconds, None


def compare_sidtime(jd, lon, desc, verbose=False):
    """Compare sidtime (sidereal time) function."""
    try:
        st_swe = swe.sidtime(jd)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        st_py = pyephem.sidtime(jd)
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    diff = abs(st_swe - st_py)
    # Handle wrap-around at 24 hours
    if diff > 12:
        diff = 24 - diff

    passed = diff < TimeTolerance.SIDEREAL_TIME

    if verbose:
        print(f"\n[sidtime] {desc}")
        print(f"  Input JD: {jd}")
        print(f"  SWE: {st_swe:.8f} hours")
        print(f"  PY:  {st_py:.8f} hours")
        print(f"  Diff: {diff:.10f} hours {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[sidtime] {desc}: "
            f"SWE={st_swe:.6f}h PY={st_py:.6f}h "
            f"Diff={diff:.10f}h {status}"
        )

    return passed, diff, None


def compare_sidtime0(jd, eps, nut, desc, verbose=False):
    """Compare sidtime0 function."""
    try:
        st_swe = swe.sidtime0(jd, eps, nut)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        st_py = pyephem.sidtime0(jd, eps, nut)
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    diff = abs(st_swe - st_py)
    if diff > 12:
        diff = 24 - diff

    passed = diff < TimeTolerance.SIDEREAL_TIME

    if verbose:
        print(f"\n[sidtime0] {desc}")
        print(f"  Input: JD={jd}, eps={eps}, nut={nut}")
        print(f"  SWE: {st_swe:.8f} hours")
        print(f"  PY:  {st_py:.8f} hours")
        print(f"  Diff: {diff:.10f} hours {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[sidtime0] {desc}: "
            f"SWE={st_swe:.6f}h PY={st_py:.6f}h "
            f"Diff={diff:.10f}h {status}"
        )

    return passed, diff, None


def compare_utc_to_jd(year, month, day, hour, minute, second, desc, verbose=False):
    """Compare utc_to_jd function."""
    try:
        jd_et_swe, jd_ut_swe = swe.utc_to_jd(year, month, day, hour, minute, second, 1)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        jd_et_py, jd_ut_py = pyephem.utc_to_jd(
            year, month, day, hour, minute, second, 1
        )
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    diff_et = abs(jd_et_swe - jd_et_py)
    diff_ut = abs(jd_ut_swe - jd_ut_py)
    max_diff = max(diff_et, diff_ut)

    passed = max_diff < TimeTolerance.JULIAN_DAY

    if verbose:
        print(f"\n[utc_to_jd] {desc}")
        print(f"  Input: {year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second}")
        print(f"  JD_ET: SWE={jd_et_swe:.10f} PY={jd_et_py:.10f} Diff={diff_et:.12f}")
        print(f"  JD_UT: SWE={jd_ut_swe:.10f} PY={jd_ut_py:.10f} Diff={diff_ut:.12f}")
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[utc_to_jd] {desc}: "
            f"ET_Diff={diff_et:.12f} UT_Diff={diff_ut:.12f} {status}"
        )

    return passed, max_diff, None


def compare_time_equ(jd, desc, verbose=False):
    """Compare time_equ (equation of time) function."""
    try:
        teq_swe = swe.time_equ(jd)
    except Exception as e:
        return False, 0.0, f"SWE ERROR: {e}"

    try:
        teq_py = pyephem.time_equ(jd)
    except Exception as e:
        return False, 0.0, f"PY ERROR: {e}"

    # time_equ returns (E, ) tuple where E is equation of time in days
    e_swe = teq_swe[0] if isinstance(teq_swe, tuple) else teq_swe
    e_py = teq_py[0] if isinstance(teq_py, tuple) else teq_py

    diff = abs(e_swe - e_py)
    passed = diff < 1e-6  # Very small tolerance for equation of time

    if verbose:
        print(f"\n[time_equ] {desc}")
        print(f"  Input JD: {jd}")
        print(f"  SWE: {e_swe:.10f} days ({e_swe * 1440:.4f} minutes)")
        print(f"  PY:  {e_py:.10f} days ({e_py * 1440:.4f} minutes)")
        print(f"  Diff: {diff:.12f} days {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[time_equ] {desc}: "
            f"SWE={e_swe * 1440:.4f}min PY={e_py * 1440:.4f}min "
            f"Diff={diff * 1440:.6f}min {status}"
        )

    return passed, diff, None


# ============================================================================
# ROUNDTRIP TESTS
# ============================================================================


def test_julday_revjul_roundtrip(year, month, day, hour, cal_type, desc, verbose=False):
    """Test that julday -> revjul gives back original values."""
    # Forward
    jd_py = pyephem.julday(year, month, day, hour, cal_type)

    # Reverse
    y2, m2, d2, h2 = pyephem.revjul(jd_py, cal_type)

    # Compare
    year_match = int(y2) == year
    month_match = int(m2) == month
    day_match = int(d2) == day
    hour_diff = abs(h2 - hour)

    passed = year_match and month_match and day_match and hour_diff < 1e-8

    if verbose:
        print(f"\n[roundtrip] {desc}")
        print(f"  Original: {year}-{month:02d}-{day:02d} {hour}h")
        print(f"  JD: {jd_py:.10f}")
        print(f"  Recovered: {int(y2)}-{int(m2):02d}-{int(d2):02d} {h2:.8f}h")
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[roundtrip] {desc}: "
            f"Original={year}-{month:02d}-{day:02d} "
            f"Recovered={int(y2)}-{int(m2):02d}-{int(d2):02d} {status}"
        )

    return passed, hour_diff, None


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all time function comparisons."""
    print_header("TIME FUNCTIONS COMPARISON")
    stats = TestStatistics()

    # 1. Julian Day
    print_section("JULIAN DAY (julday)")
    for year, month, day, hour, cal_type, desc in JULDAY_TEST_CASES:
        passed, diff, error = compare_julday(
            year, month, day, hour, cal_type, desc, verbose
        )
        stats.add_result(passed, diff, error is not None)

    # 2. Reverse Julian Day
    print_section("REVERSE JULIAN DAY (revjul)")
    test_jds = [
        (2451545.0, "J2000.0"),
        (2460000.0, "JD 2460000"),
        (2415020.5, "1900-01-01"),
        (2488070.0, "2100-01-01"),
        (0.0, "JD Zero"),
    ]
    for jd, desc in test_jds:
        passed, diff, error = compare_revjul(jd, desc, verbose)
        stats.add_result(passed, diff, error is not None)

    # 3. Delta T
    print_section("DELTA T (deltat)")
    delta_t_jds = [
        (swe.julday(1900, 1, 1, 0), "1900"),
        (swe.julday(1950, 1, 1, 0), "1950"),
        (swe.julday(2000, 1, 1, 12), "2000"),
        (swe.julday(2020, 1, 1, 0), "2020"),
        (swe.julday(2024, 6, 15, 12), "2024"),
    ]
    for jd, desc in delta_t_jds:
        passed, diff, error = compare_deltat(jd, desc, verbose)
        stats.add_result(passed, diff, error is not None)

    # 4. Sidereal Time
    print_section("SIDEREAL TIME (sidtime)")
    for name, year, month, day, hour, lat, lon, alt in STANDARD_SUBJECTS[:4]:
        jd = swe.julday(year, month, day, hour)
        passed, diff, error = compare_sidtime(jd, lon, f"{name} ({year})", verbose)
        stats.add_result(passed, diff, error is not None)

    # 5. Sidereal Time with specified obliquity/nutation
    print_section("SIDEREAL TIME (sidtime0)")
    # Get standard obliquity and nutation for J2000
    jd_test = 2451545.0
    eps = 23.439291  # Mean obliquity at J2000
    nut = 0.00256  # Typical nutation in longitude
    test_cases = [
        (2451545.0, eps, nut, "J2000 standard"),
        (2451545.0, 23.5, 0.0, "Modified obliquity"),
        (2460000.0, eps, nut, "JD 2460000"),
    ]
    for jd, e, n, desc in test_cases:
        passed, diff, error = compare_sidtime0(jd, e, n, desc, verbose)
        stats.add_result(passed, diff, error is not None)

    # 6. UTC to JD
    print_section("UTC TO JD (utc_to_jd)")
    utc_cases = [
        (2000, 1, 1, 12, 0, 0.0, "J2000.0"),
        (2024, 6, 15, 14, 30, 45.5, "Modern date"),
        (1999, 12, 31, 23, 59, 59.0, "End of 1999"),
        (2024, 2, 29, 0, 0, 0.0, "Leap day 2024"),
    ]
    for year, month, day, hour, minute, second, desc in utc_cases:
        passed, diff, error = compare_utc_to_jd(
            year, month, day, hour, minute, second, desc, verbose
        )
        stats.add_result(passed, diff, error is not None)

    # 7. Equation of Time
    print_section("EQUATION OF TIME (time_equ)")
    equ_jds = [
        (swe.julday(2024, 2, 11, 12), "Feb 11 (max negative)"),
        (swe.julday(2024, 5, 14, 12), "May 14 (zero)"),
        (swe.julday(2024, 7, 26, 12), "Jul 26 (zero)"),
        (swe.julday(2024, 11, 3, 12), "Nov 3 (max positive)"),
    ]
    for jd, desc in equ_jds:
        passed, diff, error = compare_time_equ(jd, desc, verbose)
        stats.add_result(passed, diff, error is not None)

    # 8. Roundtrip tests
    print_section("ROUNDTRIP TESTS (julday -> revjul)")
    for year, month, day, hour, cal_type, desc in JULDAY_TEST_CASES[:5]:
        passed, diff, error = test_julday_revjul_roundtrip(
            year, month, day, hour, cal_type, desc, verbose
        )
        stats.add_result(passed, diff, error is not None)

    # Summary
    stats.print_summary("TIME FUNCTIONS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_time.py [OPTIONS]")
    print()
    print("Options:")
    print("  -v, --verbose    Show detailed output for each test")
    print("  -h, --help       Show this help message")
    print()


def main():
    """Main entry point."""
    args = parse_args(sys.argv)

    if args["help"]:
        print_help()
        sys.exit(0)

    passed, total = run_all_comparisons(verbose=args["verbose"])
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
