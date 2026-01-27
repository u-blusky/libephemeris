"""
Lunar Occultations Comparison Script.

Compares lunar occultation calculations between pyswisseph and libephemeris:
- lun_occult_when_glob - global occultation search
- lun_occult_when_loc - local occultation
- lun_occult_where - occultation location
"""

import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import (
    TestStatistics,
    print_header,
    print_section,
    parse_args,
    format_status,
    jd_to_date_str,
    time_diff_seconds,
    EventComparisonResult,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class OccultTolerance:
    """Tolerance thresholds for occultation comparisons."""

    TIME_SECONDS = 300.0  # 5 minutes for occultations
    POSITION_DEGREES = 0.5  # 0.5 degree for coordinates


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for occultation tests (Moon can occult these)
OCCULT_BODIES = [
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Test locations
OCCULT_LOCATIONS = [
    ("New York", 40.7128, -74.0060, 0),
    ("London", 51.5074, -0.1278, 0),
    ("Sydney", -33.8688, 151.2093, 0),
]

# Search start dates
SEARCH_DATES = [
    (2024, 1, 1, "2024 Start"),
    (2025, 1, 1, "2025 Start"),
    (2026, 1, 1, "2026 Start"),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_lun_occult_when_glob(jd_start, body_id, body_name, desc, verbose=False):
    """Compare lun_occult_when_glob function."""
    result = EventComparisonResult("lun_occult_when_glob", f"{body_name} {desc}")
    result.tolerance_seconds = OccultTolerance.TIME_SECONDS

    # SwissEphemeris
    try:
        ret_swe = swe.lun_occult_when_glob(jd_start, body_id, "", SEFLG_SWIEPH, 0)
        # Returns (retflag, tret) where tret[0] is max occultation time
        if ret_swe[0] == 0:
            result.error_swe = "no occultation found"
            print(result.format_result(verbose))
            return True, 0.0, False  # Not an error
        result.jd_swe = ret_swe[1][0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.lun_occult_when_glob(jd_start, body_id, "", SEFLG_SWIEPH, 0)
        if ret_py[0] == 0:
            result.error_py = "no occultation found"
            print(result.format_result(verbose))
            return True, 0.0, False
        result.jd_py = ret_py[1][0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_lun_occult_when_loc(
    jd_start, body_id, body_name, lat, lon, alt, loc_name, desc, verbose=False
):
    """Compare lun_occult_when_loc function."""
    result = EventComparisonResult(
        "lun_occult_when_loc", f"{body_name} @ {loc_name} {desc}"
    )
    result.tolerance_seconds = OccultTolerance.TIME_SECONDS

    geopos = (lon, lat, alt)

    # SwissEphemeris
    try:
        ret_swe = swe.lun_occult_when_loc(
            jd_start, body_id, "", SEFLG_SWIEPH, geopos, 0
        )
        if ret_swe[0] == 0:
            result.error_swe = "no occultation visible"
            print(result.format_result(verbose))
            return True, 0.0, False
        result.jd_swe = ret_swe[1][0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.lun_occult_when_loc(
            jd_start, body_id, "", SEFLG_SWIEPH, geopos, 0
        )
        if ret_py[0] == 0:
            result.error_py = "no occultation visible"
            print(result.format_result(verbose))
            return True, 0.0, False
        result.jd_py = ret_py[1][0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_lun_occult_where(jd, body_id, body_name, verbose=False):
    """Compare lun_occult_where function."""
    # SwissEphemeris
    try:
        ret_swe = swe.lun_occult_where(jd, body_id, "", SEFLG_SWIEPH)
        if ret_swe[0] == 0:
            print(f"[lun_occult_where] {body_name}: SWE - no central line")
            return True, 0.0, False
        lon_swe, lat_swe = ret_swe[1][0], ret_swe[1][1]
    except Exception as e:
        print(f"[lun_occult_where] {body_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.lun_occult_where(jd, body_id, "", SEFLG_SWIEPH)
        if ret_py[0] == 0:
            print(f"[lun_occult_where] {body_name}: PY - no central line")
            return True, 0.0, False
        lon_py, lat_py = ret_py[1][0], ret_py[1][1]
    except Exception as e:
        print(f"[lun_occult_where] {body_name}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = abs(lon_swe - lon_py)
    diff_lat = abs(lat_swe - lat_py)
    max_diff = max(diff_lon, diff_lat)
    passed = max_diff < OccultTolerance.POSITION_DEGREES

    if verbose:
        print(f"\n[lun_occult_where] {body_name}")
        print(f"  SWE: Lon={lon_swe:.4f} Lat={lat_swe:.4f}")
        print(f"  PY:  Lon={lon_py:.4f} Lat={lat_py:.4f}")
        print(f"  Diff: {max_diff:.4f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[lun_occult_where] {body_name}: "
            f"SWE=({lon_swe:.2f},{lat_swe:.2f}) "
            f"PY=({lon_py:.2f},{lat_py:.2f}) "
            f"MaxDiff={max_diff:.4f} {status}"
        )

    return passed, max_diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all occultation comparisons."""
    print_header("LUNAR OCCULTATIONS COMPARISON")
    stats = TestStatistics()

    # 1. Global occultation search
    print_section("LUNAR OCCULTATION GLOBAL (lun_occult_when_glob)")
    for year, month, day, date_desc in SEARCH_DATES:
        jd_start = swe.julday(year, month, day, 0.0)

        for body_id, body_name in OCCULT_BODIES:
            passed, diff, error = compare_lun_occult_when_glob(
                jd_start, body_id, body_name, date_desc, verbose
            )
            stats.add_result(passed, diff, error)

    # 2. Local occultation (for first search result)
    print_section("LUNAR OCCULTATION LOCAL (lun_occult_when_loc)")
    jd_start = swe.julday(2024, 1, 1, 0.0)

    for body_id, body_name in OCCULT_BODIES[:2]:  # Venus and Mars
        for loc_name, lat, lon, alt in OCCULT_LOCATIONS[:2]:
            passed, diff, error = compare_lun_occult_when_loc(
                jd_start, body_id, body_name, lat, lon, alt, loc_name, "2024", verbose
            )
            stats.add_result(passed, diff, error)

    # 3. Occultation location (where on Earth)
    print_section("LUNAR OCCULTATION LOCATION (lun_occult_where)")
    # Find an occultation first, then check its location
    for body_id, body_name in OCCULT_BODIES[:2]:
        jd_start = swe.julday(2024, 1, 1, 0.0)
        try:
            ret = swe.lun_occult_when_glob(jd_start, body_id, "", SEFLG_SWIEPH, 0)
            if ret[0] != 0:
                jd_occult = ret[1][0]
                passed, diff, error = compare_lun_occult_where(
                    jd_occult, body_id, body_name, verbose
                )
                stats.add_result(passed, diff, error)
        except Exception:
            pass

    # Summary
    stats.print_summary("LUNAR OCCULTATIONS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_occultations.py [OPTIONS]")
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
