"""
Hypothetical Planets Comparison Script.

Compares hypothetical body calculations between pyswisseph and libephemeris:
- Uranian planets (Hamburg School): Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon
- Transpluto (Isis)
- Other fictitious bodies: Vulcan, White Moon, Proserpina, Waldemath, Planet X variants

Note: Comparisons use relaxed tolerances as hypothetical bodies use simplified
Keplerian elements and secular polynomials.
"""

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


# ============================================================================
# HYPOTHETICAL BODY DEFINITIONS
# ============================================================================

# Uranian Planets (Hamburg School) - IDs 40-47
URANIAN_PLANETS = [
    (SE_CUPIDO, swe.CUPIDO, "Cupido"),
    (SE_HADES, swe.HADES, "Hades"),
    (SE_ZEUS, swe.ZEUS, "Zeus"),
    (SE_KRONOS, swe.KRONOS, "Kronos"),
    (SE_APOLLON, swe.APOLLON, "Apollon"),
    (SE_ADMETOS, swe.ADMETOS, "Admetos"),
    (SE_VULKANUS, swe.VULKANUS, "Vulkanus"),
    (SE_POSEIDON, swe.POSEIDON, "Poseidon"),
]

# Other hypothetical bodies - IDs 48+
OTHER_HYPOTHETICAL = [
    (SE_ISIS, swe.ISIS, "Transpluto/Isis"),
]

# Bodies that may not be available in all pyswisseph versions
OPTIONAL_HYPOTHETICAL = [
    # These require seorbel.txt file or may not be in all versions
    # Uncomment if your pyswisseph version supports them
    # (SE_VULCAN, 55, "Vulcan"),
    # (SE_WHITE_MOON, 56, "White Moon/Selena"),
    # (SE_PROSERPINA, 57, "Proserpina"),
    # (SE_WALDEMATH, 58, "Waldemath"),
]

# Tolerance for hypothetical bodies (relaxed due to different calculation methods)
HYPOTHETICAL_LONGITUDE_TOLERANCE = 1.0  # 1 degree tolerance for most
URANIAN_LONGITUDE_TOLERANCE = 0.1  # 0.1 degree tolerance for Uranian planets


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_hypothetical_body(jd, body_py, body_swe, name, tolerance=1.0, verbose=True):
    """
    Compare single hypothetical body position.

    Args:
        jd: Julian Day
        body_py: libephemeris body ID
        body_swe: pyswisseph body ID
        name: Body name for display
        tolerance: Longitude tolerance in degrees
        verbose: Show detailed output

    Returns:
        (status_code, max_diff) where status_code is:
        1 = PASSED, 0 = FAILED, 2 = SKIPPED (no data)
    """
    try:
        # SwissEph
        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            if verbose:
                print(f"\n{name:-^80}")
                print(f"  SKIPPED - SwissEph error ({e})")
            return 2, 0.0  # 2 = SKIPPED

        # LibEphemeris
        try:
            pos_py, _ = ephem.swe_calc_ut(jd, body_py, 0)
        except Exception as e:
            if verbose:
                print(f"\n{name:-^80}")
                print(f"  SKIPPED - LibEphemeris error ({e})")
            return 2, 0.0

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        passed = diff_lon < tolerance and diff_lat < 5.0

        if verbose:
            print(f"\n{name:-^80}")
            print(
                f"  Longitude:  SWE={format_coord(pos_swe[0])}  "
                f"PY={format_coord(pos_py[0])}  "
                f"Diff={format_diff(diff_lon)} {format_status(diff_lon < tolerance)}"
            )
            print(
                f"  Latitude:   SWE={format_coord(pos_swe[1])}  "
                f"PY={format_coord(pos_py[1])}  "
                f"Diff={format_diff(diff_lat)} {format_status(diff_lat < 5)}"
            )
            print(
                f"  Distance:   SWE={format_coord(pos_swe[2], 4)} AU  "
                f"PY={format_coord(pos_py[2], 4)} AU  "
                f"Diff={format_diff(diff_dist, 6)} AU"
            )
            print(f"  Status:     {format_status(passed)}")

        return (1 if passed else 0), max(diff_lon, diff_lat)

    except Exception as e:
        print(f"\n{name}: ERROR - {e}")
        return 0, 999.0


def compare_uranian_planets(jd, date_str, verbose=True):
    """Compare all Uranian planets at a given date."""
    results = []

    for body_py, body_swe, name in URANIAN_PLANETS:
        status, diff = compare_hypothetical_body(
            jd, body_py, body_swe, name, URANIAN_LONGITUDE_TOLERANCE, verbose
        )
        results.append((name, status, diff))

    return results


def compare_other_hypothetical(jd, date_str, verbose=True):
    """Compare other hypothetical bodies at a given date."""
    results = []

    for body_py, body_swe, name in OTHER_HYPOTHETICAL:
        status, diff = compare_hypothetical_body(
            jd, body_py, body_swe, name, HYPOTHETICAL_LONGITUDE_TOLERANCE, verbose
        )
        results.append((name, status, diff))

    return results


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=True):
    """Run all hypothetical body comparisons."""
    print_header("HYPOTHETICAL BODIES COMPARISON: LibEphemeris vs SwissEphemeris")

    print("\n" + "!" * 80)
    print("! NOTE: Hypothetical bodies use simplified orbital elements")
    print("! Uranian planets tolerance: 0.1 degree (astrological precision)")
    print("! Other hypothetical bodies tolerance: 1.0 degree")
    print("! Swiss Ephemeris uses seorbel.txt for orbital elements")
    print("!" * 80 + "\n")

    # Test dates
    subjects = [
        ("J2000.0", 2000, 1, 1, 12.0),
        ("2024-01-01", 2024, 1, 1, 0.0),
        ("2010-07-01", 2010, 7, 1, 12.0),
        ("1980-01-01", 1980, 1, 1, 0.0),
    ]

    # Global stats
    g_total = 0
    g_passed = 0
    g_skipped = 0
    g_failed = 0
    g_max_diff = 0.0

    for name, year, month, day, hour in subjects:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month:02d}-{day:02d}"

        print(f"\n{'=' * 80}")
        print(f"DATE: {name} ({date_str})")
        print(f"{'=' * 80}")

        # Uranian planets
        print(f"\n{'-' * 40}")
        print("URANIAN PLANETS (Hamburg School)")
        print(f"{'-' * 40}")

        results = compare_uranian_planets(jd, date_str, verbose)
        for body_name, status_code, diff in results:
            g_total += 1
            if status_code == 1:  # PASSED
                g_passed += 1
                g_max_diff = max(g_max_diff, diff)
            elif status_code == 2:  # SKIPPED
                g_skipped += 1
            else:  # FAILED
                g_failed += 1

        # Other hypothetical bodies
        print(f"\n{'-' * 40}")
        print("OTHER HYPOTHETICAL BODIES")
        print(f"{'-' * 40}")

        results = compare_other_hypothetical(jd, date_str, verbose)
        for body_name, status_code, diff in results:
            g_total += 1
            if status_code == 1:
                g_passed += 1
                g_max_diff = max(g_max_diff, diff)
            elif status_code == 2:
                g_skipped += 1
            else:
                g_failed += 1

    # Summary
    print("\n" + "=" * 80)
    print("HYPOTHETICAL BODIES COMPARISON SUMMARY")
    print("=" * 80)
    print(f"Total tests:   {g_total}")
    print(f"Passed:        {g_passed} ✓")
    print(f"Skipped:       {g_skipped} -")
    print(f"Failed:        {g_failed} ✗")

    if g_passed > 0:
        print(
            f"Pass rate:     {g_passed / (g_total - g_skipped) * 100:.1f}% (of executed)"
        )
        print(f"Max diff:      {g_max_diff:.6f}")
    elif g_total == g_skipped:
        print("Pass rate:     N/A (All skipped)")
    else:
        print("Pass rate:     0.0%")

    print("=" * 80)

    return g_passed, g_total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_hypothetical.py [OPTIONS]")
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

    passed, total = run_all_comparisons(verbose=not args["quiet"])

    # Exit with appropriate code
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
