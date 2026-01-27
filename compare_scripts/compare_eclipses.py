"""
Eclipse Functions Comparison Script.

Compares solar and lunar eclipse calculations between pyswisseph and libephemeris:
- Solar eclipse global (when, where)
- Solar eclipse local (when, how)
- Lunar eclipse (when, how)
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
    STANDARD_SUBJECTS,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class EclipseTolerance:
    """Tolerance thresholds for eclipse comparisons."""

    TIME_SECONDS = 120.0  # 2 minutes for eclipse times
    POSITION_DEGREES = 0.1  # 0.1 degree for coordinates
    MAGNITUDE = 0.01  # Eclipse magnitude


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Known eclipses for validation
KNOWN_SOLAR_ECLIPSES = [
    # (start_jd, description) - search from this date
    (swe.julday(2017, 8, 1, 0), "Great American Eclipse 2017"),
    (swe.julday(2024, 4, 1, 0), "North American Eclipse 2024"),
    (swe.julday(2020, 12, 1, 0), "South American Eclipse 2020"),
    (swe.julday(2023, 4, 1, 0), "Hybrid Eclipse 2023"),
    (swe.julday(2026, 8, 1, 0), "2026 Total Eclipse"),
]

KNOWN_LUNAR_ECLIPSES = [
    (swe.julday(2022, 5, 1, 0), "Blood Moon May 2022"),
    (swe.julday(2022, 11, 1, 0), "Blood Moon Nov 2022"),
    (swe.julday(2025, 3, 1, 0), "Lunar Eclipse Mar 2025"),
    (swe.julday(2025, 9, 1, 0), "Lunar Eclipse Sep 2025"),
]

# Locations for local eclipse tests
ECLIPSE_LOCATIONS = [
    ("New York", 40.7128, -74.0060, 0),
    ("London", 51.5074, -0.1278, 0),
    ("Sydney", -33.8688, 151.2093, 0),
    ("Tokyo", 35.6762, 139.6503, 0),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_sol_eclipse_when_glob(jd_start, desc, verbose=False):
    """Compare sol_eclipse_when_glob function."""
    result = EventComparisonResult("sol_eclipse_when_glob", desc)
    result.tolerance_seconds = EclipseTolerance.TIME_SECONDS

    # SwissEphemeris
    try:
        ret_swe = swe.sol_eclipse_when_glob(jd_start, SEFLG_SWIEPH, 0)
        # Returns (retflag, tret) where tret[0] is max eclipse time
        result.jd_swe = ret_swe[1][0]
        result.extra_info["type_swe"] = ret_swe[0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.sol_eclipse_when_glob(jd_start, SEFLG_SWIEPH, 0)
        result.jd_py = ret_py[1][0]
        result.extra_info["type_py"] = ret_py[0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_sol_eclipse_where(jd, desc, verbose=False):
    """Compare sol_eclipse_where function."""
    # SwissEphemeris
    try:
        ret_swe = swe.sol_eclipse_where(jd, SEFLG_SWIEPH)
        # Returns (retflag, geopos, attr) where geopos[0]=lon, geopos[1]=lat
        lon_swe, lat_swe = ret_swe[1][0], ret_swe[1][1]
    except Exception as e:
        print(f"[sol_eclipse_where] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.sol_eclipse_where(jd, SEFLG_SWIEPH)
        lon_py, lat_py = ret_py[1][0], ret_py[1][1]
    except Exception as e:
        print(f"[sol_eclipse_where] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = abs(lon_swe - lon_py)
    diff_lat = abs(lat_swe - lat_py)
    max_diff = max(diff_lon, diff_lat)
    passed = max_diff < EclipseTolerance.POSITION_DEGREES

    if verbose:
        print(f"\n[sol_eclipse_where] {desc}")
        print(f"  SWE: Lon={lon_swe:.4f} Lat={lat_swe:.4f}")
        print(f"  PY:  Lon={lon_py:.4f} Lat={lat_py:.4f}")
        print(f"  Diff: Lon={diff_lon:.4f} Lat={diff_lat:.4f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[sol_eclipse_where] {desc}: "
            f"SWE=({lon_swe:.2f},{lat_swe:.2f}) "
            f"PY=({lon_py:.2f},{lat_py:.2f}) "
            f"MaxDiff={max_diff:.4f} {status}"
        )

    return passed, max_diff, False


def compare_sol_eclipse_when_loc(
    jd_start, lat, lon, alt, loc_name, desc, verbose=False
):
    """Compare sol_eclipse_when_loc function."""
    result = EventComparisonResult("sol_eclipse_when_loc", f"{desc} @ {loc_name}")
    result.tolerance_seconds = (
        EclipseTolerance.TIME_SECONDS * 2
    )  # More relaxed for local

    geopos = (lon, lat, alt)

    # SwissEphemeris
    try:
        ret_swe = swe.sol_eclipse_when_loc(jd_start, SEFLG_SWIEPH, geopos)
        result.jd_swe = ret_swe[1][0]  # Maximum eclipse time
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.sol_eclipse_when_loc(jd_start, SEFLG_SWIEPH, geopos)
        result.jd_py = ret_py[1][0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_sol_eclipse_how(jd, lat, lon, alt, loc_name, verbose=False):
    """Compare sol_eclipse_how function."""
    geopos = (lon, lat, alt)

    # SwissEphemeris
    try:
        ret_swe = swe.sol_eclipse_how(jd, SEFLG_SWIEPH, geopos)
        # Returns (retflag, attr) where attr contains eclipse details
        mag_swe = ret_swe[1][0]  # Eclipse magnitude
    except Exception as e:
        print(f"[sol_eclipse_how] @ {loc_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.sol_eclipse_how(jd, SEFLG_SWIEPH, geopos)
        mag_py = ret_py[1][0]
    except Exception as e:
        print(f"[sol_eclipse_how] @ {loc_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(mag_swe - mag_py)
    passed = diff < EclipseTolerance.MAGNITUDE

    if verbose:
        print(f"\n[sol_eclipse_how] @ {loc_name}")
        print(f"  SWE magnitude: {mag_swe:.6f}")
        print(f"  PY magnitude:  {mag_py:.6f}")
        print(f"  Diff: {diff:.6f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[sol_eclipse_how] @ {loc_name}: "
            f"SWE_mag={mag_swe:.4f} PY_mag={mag_py:.4f} "
            f"Diff={diff:.6f} {status}"
        )

    return passed, diff, False


def compare_lun_eclipse_when(jd_start, desc, verbose=False):
    """Compare lun_eclipse_when function."""
    result = EventComparisonResult("lun_eclipse_when", desc)
    result.tolerance_seconds = EclipseTolerance.TIME_SECONDS

    # SwissEphemeris
    try:
        ret_swe = swe.lun_eclipse_when(jd_start, SEFLG_SWIEPH, 0)
        result.jd_swe = ret_swe[1][0]  # Maximum eclipse time
        result.extra_info["type_swe"] = ret_swe[0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.lun_eclipse_when(jd_start, SEFLG_SWIEPH, 0)
        result.jd_py = ret_py[1][0]
        result.extra_info["type_py"] = ret_py[0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_lun_eclipse_how(jd, lat, lon, alt, loc_name, verbose=False):
    """Compare lun_eclipse_how function."""
    geopos = (lon, lat, alt)

    # SwissEphemeris
    try:
        ret_swe = swe.lun_eclipse_how(jd, SEFLG_SWIEPH, geopos)
        mag_swe = ret_swe[1][0]  # Umbral magnitude
    except Exception as e:
        print(f"[lun_eclipse_how] @ {loc_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.lun_eclipse_how(jd, SEFLG_SWIEPH, geopos)
        mag_py = ret_py[1][0]
    except Exception as e:
        print(f"[lun_eclipse_how] @ {loc_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(mag_swe - mag_py)
    passed = diff < EclipseTolerance.MAGNITUDE

    if verbose:
        print(f"\n[lun_eclipse_how] @ {loc_name}")
        print(f"  SWE magnitude: {mag_swe:.6f}")
        print(f"  PY magnitude:  {mag_py:.6f}")
        print(f"  Diff: {diff:.6f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[lun_eclipse_how] @ {loc_name}: "
            f"SWE_mag={mag_swe:.4f} PY_mag={mag_py:.4f} "
            f"Diff={diff:.6f} {status}"
        )

    return passed, diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all eclipse comparisons."""
    print_header("ECLIPSE FUNCTIONS COMPARISON")
    stats = TestStatistics()

    # 1. Solar Eclipse Global - When
    print_section("SOLAR ECLIPSE GLOBAL (sol_eclipse_when_glob)")
    for jd_start, desc in KNOWN_SOLAR_ECLIPSES:
        passed, diff, error = compare_sol_eclipse_when_glob(jd_start, desc, verbose)
        stats.add_result(passed, diff, error)

    # 2. Solar Eclipse Global - Where (for found eclipses)
    print_section("SOLAR ECLIPSE LOCATION (sol_eclipse_where)")
    for jd_start, desc in KNOWN_SOLAR_ECLIPSES[:3]:
        # First find the eclipse
        try:
            ret = swe.sol_eclipse_when_glob(jd_start, SEFLG_SWIEPH, 0)
            jd_eclipse = ret[1][0]
            passed, diff, error = compare_sol_eclipse_where(jd_eclipse, desc, verbose)
            stats.add_result(passed, diff, error)
        except Exception as e:
            print(f"[sol_eclipse_where] {desc}: SKIP (no eclipse found)")

    # 3. Solar Eclipse Local - When
    print_section("SOLAR ECLIPSE LOCAL (sol_eclipse_when_loc)")
    for jd_start, desc in KNOWN_SOLAR_ECLIPSES[:2]:
        for loc_name, lat, lon, alt in ECLIPSE_LOCATIONS[:2]:
            passed, diff, error = compare_sol_eclipse_when_loc(
                jd_start, lat, lon, alt, loc_name, desc, verbose
            )
            stats.add_result(passed, diff, error)

    # 4. Solar Eclipse How (magnitude at location)
    print_section("SOLAR ECLIPSE DETAILS (sol_eclipse_how)")
    for jd_start, desc in KNOWN_SOLAR_ECLIPSES[:1]:
        try:
            ret = swe.sol_eclipse_when_glob(jd_start, SEFLG_SWIEPH, 0)
            jd_eclipse = ret[1][0]
            for loc_name, lat, lon, alt in ECLIPSE_LOCATIONS[:2]:
                passed, diff, error = compare_sol_eclipse_how(
                    jd_eclipse, lat, lon, alt, loc_name, verbose
                )
                stats.add_result(passed, diff, error)
        except Exception:
            pass

    # 5. Lunar Eclipse - When
    print_section("LUNAR ECLIPSE (lun_eclipse_when)")
    for jd_start, desc in KNOWN_LUNAR_ECLIPSES:
        passed, diff, error = compare_lun_eclipse_when(jd_start, desc, verbose)
        stats.add_result(passed, diff, error)

    # 6. Lunar Eclipse How
    print_section("LUNAR ECLIPSE DETAILS (lun_eclipse_how)")
    for jd_start, desc in KNOWN_LUNAR_ECLIPSES[:2]:
        try:
            ret = swe.lun_eclipse_when(jd_start, SEFLG_SWIEPH, 0)
            jd_eclipse = ret[1][0]
            for loc_name, lat, lon, alt in ECLIPSE_LOCATIONS[:2]:
                passed, diff, error = compare_lun_eclipse_how(
                    jd_eclipse, lat, lon, alt, loc_name, verbose
                )
                stats.add_result(passed, diff, error)
        except Exception:
            pass

    # Summary
    stats.print_summary("ECLIPSE FUNCTIONS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_eclipses.py [OPTIONS]")
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
