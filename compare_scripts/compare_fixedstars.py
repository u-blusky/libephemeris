"""
Fixed Stars Comparison Script.

Compares fixed star calculations between pyswisseph and libephemeris:
- fixstar_ut / fixstar2_ut - positions
- fixstar_mag / fixstar2_mag - magnitudes
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
    angular_diff,
    format_coord,
    format_diff,
    FIXED_STARS,
    STANDARD_SUBJECTS,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class StarTolerance:
    """Tolerance thresholds for fixed star comparisons."""

    POSITION_DEGREES = 0.001  # 0.001 degree (3.6 arcsec)
    MAGNITUDE = 0.1  # 0.1 magnitude


# ============================================================================
# ADDITIONAL STARS FOR TESTING
# ============================================================================

# Star names that work with both implementations
ADDITIONAL_STARS = [
    "Polaris",
    "Castor",
    "Achernar",
    "Canopus",
    "Algol",
    "Pleiades",  # Alcyone
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_fixstar_ut(star_name, jd, desc, verbose=False):
    """Compare fixstar_ut function for a single star."""
    # SwissEphemeris
    try:
        ret_swe = swe.fixstar_ut(star_name, jd, SEFLG_SWIEPH)
        # Returns (name, pos, retflags) where pos is 6-element tuple
        lon_swe = ret_swe[1][0]
        lat_swe = ret_swe[1][1]
        dist_swe = ret_swe[1][2]
    except Exception as e:
        print(f"[fixstar_ut] {star_name} {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.fixstar_ut(star_name, jd, SEFLG_SWIEPH)
        lon_py = ret_py[1][0]
        lat_py = ret_py[1][1]
        dist_py = ret_py[1][2]
    except Exception as e:
        print(f"[fixstar_ut] {star_name} {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = angular_diff(lon_swe, lon_py)
    diff_lat = abs(lat_swe - lat_py)
    max_diff = max(diff_lon, diff_lat)
    passed = max_diff < StarTolerance.POSITION_DEGREES

    if verbose:
        print(f"\n[fixstar_ut] {star_name} - {desc}")
        print(f"  SWE: Lon={lon_swe:.6f} Lat={lat_swe:.6f}")
        print(f"  PY:  Lon={lon_py:.6f} Lat={lat_py:.6f}")
        print(f"  Diff: Lon={diff_lon:.8f} Lat={diff_lat:.8f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[fixstar_ut] {star_name:<12} {desc}: "
            f"Lon={lon_swe:.4f}/{lon_py:.4f} "
            f"Lat={lat_swe:.4f}/{lat_py:.4f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_fixstar2_ut(star_name, jd, desc, verbose=False):
    """Compare fixstar2_ut function (uses catalog number or name)."""
    # SwissEphemeris
    try:
        ret_swe = swe.fixstar2_ut(star_name, jd, SEFLG_SWIEPH)
        lon_swe = ret_swe[1][0]
        lat_swe = ret_swe[1][1]
        full_name_swe = ret_swe[0]
    except Exception as e:
        print(f"[fixstar2_ut] {star_name} {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.fixstar2_ut(star_name, jd, SEFLG_SWIEPH)
        lon_py = ret_py[1][0]
        lat_py = ret_py[1][1]
        full_name_py = ret_py[0]
    except Exception as e:
        print(f"[fixstar2_ut] {star_name} {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = angular_diff(lon_swe, lon_py)
    diff_lat = abs(lat_swe - lat_py)
    max_diff = max(diff_lon, diff_lat)
    passed = max_diff < StarTolerance.POSITION_DEGREES

    if verbose:
        print(f"\n[fixstar2_ut] {star_name} - {desc}")
        print(f"  SWE: {full_name_swe} Lon={lon_swe:.6f} Lat={lat_swe:.6f}")
        print(f"  PY:  {full_name_py} Lon={lon_py:.6f} Lat={lat_py:.6f}")
        print(f"  Diff: Lon={diff_lon:.8f} Lat={diff_lat:.8f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[fixstar2_ut] {star_name:<12} {desc}: "
            f"Lon={lon_swe:.4f}/{lon_py:.4f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_fixstar_mag(star_name, verbose=False):
    """Compare fixstar_mag function."""
    # SwissEphemeris
    try:
        ret_swe = swe.fixstar_mag(star_name)
        mag_swe = ret_swe[1]
    except Exception as e:
        print(f"[fixstar_mag] {star_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.fixstar_mag(star_name)
        mag_py = ret_py[1]
    except Exception as e:
        print(f"[fixstar_mag] {star_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(mag_swe - mag_py)
    passed = diff < StarTolerance.MAGNITUDE

    if verbose:
        print(f"\n[fixstar_mag] {star_name}")
        print(f"  SWE: {mag_swe:.3f}")
        print(f"  PY:  {mag_py:.3f}")
        print(f"  Diff: {diff:.4f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[fixstar_mag] {star_name:<12}: "
            f"SWE={mag_swe:.2f} PY={mag_py:.2f} "
            f"Diff={diff:.4f} {status}"
        )

    return passed, diff, False


def compare_fixstar2_mag(star_name, verbose=False):
    """Compare fixstar2_mag function."""
    # SwissEphemeris
    try:
        ret_swe = swe.fixstar2_mag(star_name)
        mag_swe = ret_swe[1]
        name_swe = ret_swe[0]
    except Exception as e:
        print(f"[fixstar2_mag] {star_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.fixstar2_mag(star_name)
        mag_py = ret_py[1]
        name_py = ret_py[0]
    except Exception as e:
        print(f"[fixstar2_mag] {star_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(mag_swe - mag_py)
    passed = diff < StarTolerance.MAGNITUDE

    if verbose:
        print(f"\n[fixstar2_mag] {star_name}")
        print(f"  SWE: {name_swe} mag={mag_swe:.3f}")
        print(f"  PY:  {name_py} mag={mag_py:.3f}")
        print(f"  Diff: {diff:.4f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[fixstar2_mag] {star_name:<12}: "
            f"SWE={mag_swe:.2f} PY={mag_py:.2f} "
            f"Diff={diff:.4f} {status}"
        )

    return passed, diff, False


def compare_star_with_speed(star_name, jd, verbose=False):
    """Compare star position with proper motion (speed flag)."""
    flags = SEFLG_SWIEPH | SEFLG_SPEED

    # SwissEphemeris
    try:
        ret_swe = swe.fixstar_ut(star_name, jd, flags)
        lon_swe = ret_swe[1][0]
        lon_speed_swe = ret_swe[1][3]  # Proper motion in lon
    except Exception as e:
        print(f"[fixstar+speed] {star_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.fixstar_ut(star_name, jd, flags)
        lon_py = ret_py[1][0]
        lon_speed_py = ret_py[1][3]
    except Exception as e:
        print(f"[fixstar+speed] {star_name}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = angular_diff(lon_swe, lon_py)
    diff_speed = abs(lon_speed_swe - lon_speed_py)
    passed = diff_lon < StarTolerance.POSITION_DEGREES

    if verbose:
        print(f"\n[fixstar+speed] {star_name}")
        print(f"  SWE: Lon={lon_swe:.6f} dLon/dt={lon_speed_swe:.10f}")
        print(f"  PY:  Lon={lon_py:.6f} dLon/dt={lon_speed_py:.10f}")
        print(
            f"  Diff: Lon={diff_lon:.8f} Speed={diff_speed:.12f} {format_status(passed)}"
        )
    else:
        status = format_status(passed)
        print(
            f"[fixstar+speed] {star_name:<12}: "
            f"Lon_Diff={diff_lon:.6f} Speed_Diff={diff_speed:.10f} {status}"
        )

    return passed, diff_lon, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all fixed star comparisons."""
    print_header("FIXED STARS COMPARISON")
    stats = TestStatistics()

    # Test dates
    test_jds = [
        (swe.julday(2000, 1, 1, 12), "J2000.0"),
        (swe.julday(2024, 6, 15, 0), "2024"),
        (swe.julday(1950, 1, 1, 0), "1950"),
    ]

    # 1. fixstar_ut - Main stars
    print_section("FIXED STAR POSITIONS (fixstar_ut)")
    for star_name in FIXED_STARS:
        for jd, desc in test_jds[:2]:  # Test first 2 dates
            passed, diff, error = compare_fixstar_ut(star_name, jd, desc, verbose)
            stats.add_result(passed, diff, error)

    # 2. fixstar2_ut - Additional stars
    print_section("FIXED STAR POSITIONS (fixstar2_ut)")
    jd = swe.julday(2024, 1, 1, 0)
    for star_name in FIXED_STARS[:8] + ADDITIONAL_STARS:
        passed, diff, error = compare_fixstar2_ut(star_name, jd, "2024", verbose)
        stats.add_result(passed, diff, error)

    # 3. fixstar_mag - Magnitudes
    print_section("FIXED STAR MAGNITUDES (fixstar_mag)")
    for star_name in FIXED_STARS:
        passed, diff, error = compare_fixstar_mag(star_name, verbose)
        stats.add_result(passed, diff, error)

    # 4. fixstar2_mag - Magnitudes with name search
    print_section("FIXED STAR MAGNITUDES (fixstar2_mag)")
    for star_name in FIXED_STARS[:8]:
        passed, diff, error = compare_fixstar2_mag(star_name, verbose)
        stats.add_result(passed, diff, error)

    # 5. Stars with proper motion (speed flag)
    print_section("FIXED STARS WITH PROPER MOTION")
    jd = swe.julday(2024, 1, 1, 0)
    for star_name in ["Sirius", "Arcturus", "Barnard's Star", "Vega"]:
        passed, diff, error = compare_star_with_speed(star_name, jd, verbose)
        stats.add_result(passed, diff, error)

    # 6. Precession test - same star at different epochs
    print_section("PRECESSION TEST")
    for star_name in ["Regulus", "Spica", "Aldebaran"]:
        for jd, desc in test_jds:
            passed, diff, error = compare_fixstar_ut(star_name, jd, desc, verbose)
            stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("FIXED STARS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_fixedstars.py [OPTIONS]")
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
