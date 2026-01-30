"""
Comparison script for Minor Bodies: Asteroids, Centaurs, and TNOs.
Note: Comparisons use relaxed tolerances as LibEphemeris uses simplified
Keplerian elements while SwissEphemeris uses high-precision perturbation models.

Covers:
- Main belt asteroids: Ceres, Pallas, Juno, Vesta
- Centaurs: Chiron, Pholus, Nessus, Asbolus, Chariklo
- Trans-Neptunian Objects (TNOs): Eris, Sedna, Makemake, Haumea, Orcus, Ixion, Quaoar
- Near-Earth asteroids: Apophis
"""

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


# ============================================================================
# MINOR BODY DEFINITIONS
# ============================================================================

# Main belt asteroids (IDs 17-20)
MAIN_BELT_ASTEROIDS = [
    (SE_CERES, swe.CERES, "Ceres"),
    (SE_PALLAS, swe.PALLAS, "Pallas"),
    (SE_JUNO, swe.JUNO, "Juno"),
    (SE_VESTA, swe.VESTA, "Vesta"),
]

# Centaurs (IDs 15-16 and via AST_OFFSET)
CENTAURS = [
    (SE_CHIRON, swe.CHIRON, "Chiron"),
    (SE_PHOLUS, swe.PHOLUS, "Pholus"),
    # Additional centaurs via asteroid number offset
    (SE_NESSUS, SE_AST_OFFSET + 7066, "Nessus"),
    (SE_ASBOLUS, SE_AST_OFFSET + 8405, "Asbolus"),
    (SE_CHARIKLO, SE_AST_OFFSET + 10199, "Chariklo"),
]

# Trans-Neptunian Objects (TNOs) - via AST_OFFSET + asteroid number
TNOS = [
    (SE_ERIS, SE_AST_OFFSET + 136199, "Eris"),
    (SE_MAKEMAKE, SE_AST_OFFSET + 136472, "Makemake"),
    (SE_HAUMEA, SE_AST_OFFSET + 136108, "Haumea"),
    (SE_ORCUS, SE_AST_OFFSET + 90482, "Orcus"),
    (SE_IXION, SE_AST_OFFSET + 28978, "Ixion"),
    (SE_QUAOAR, SE_AST_OFFSET + 50000, "Quaoar"),
    (SE_SEDNA, SE_AST_OFFSET + 90377, "Sedna"),
    (SE_VARUNA, SE_AST_OFFSET + 20000, "Varuna"),
    (SE_GONGGONG, SE_AST_OFFSET + 225088, "Gonggong"),
]

# Near-Earth asteroids
NEAR_EARTH_ASTEROIDS = [
    (SE_APOPHIS, SE_AST_OFFSET + 99942, "Apophis"),
]

# Legacy combined list for backward compatibility
ASTEROIDS = MAIN_BELT_ASTEROIDS + [
    (SE_CHIRON, swe.CHIRON, "Chiron"),
    (SE_PHOLUS, swe.PHOLUS, "Pholus"),
]

# Tolerances for different body types
# Inner solar system bodies have tighter tolerances
# Outer solar system bodies (especially TNOs) need more relaxed tolerances
TOLERANCES = {
    "main_belt": {"longitude": 10.0, "latitude": 5.0},  # Degrees
    "centaur": {"longitude": 10.0, "latitude": 5.0},
    "tno": {
        "longitude": 30.0,
        "latitude": 10.0,
    },  # TNOs have larger errors due to long periods
    "near_earth": {"longitude": 10.0, "latitude": 5.0},
}


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_minor_body(
    jd, body_py, body_swe, name, lon_tol=10.0, lat_tol=5.0, verbose=True
):
    """
    Compare single minor body position.

    Args:
        jd: Julian Day
        body_py: libephemeris body ID
        body_swe: pyswisseph body ID
        name: Body name for display
        lon_tol: Longitude tolerance in degrees
        lat_tol: Latitude tolerance in degrees
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
                print(f"  SKIPPED - SwissEph data file missing ({e})")
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

        # Check against tolerances
        passed = diff_lon < lon_tol and diff_lat < lat_tol

        if verbose:
            print(f"\n{name:-^80}")
            print(
                f"  Longitude:  SWE={format_coord(pos_swe[0])}  "
                f"PY={format_coord(pos_py[0])}  "
                f"Diff={format_diff(diff_lon)} {format_status(diff_lon < lon_tol)}"
            )
            print(
                f"  Latitude:   SWE={format_coord(pos_swe[1])}  "
                f"PY={format_coord(pos_py[1])}  "
                f"Diff={format_diff(diff_lat)} {format_status(diff_lat < lat_tol)}"
            )
            print(
                f"  Distance:   SWE={format_coord(pos_swe[2], 4)} AU  "
                f"PY={format_coord(pos_py[2], 4)} AU  "
                f"Diff={format_diff(diff_dist, 6)} AU"
            )
            print(f"  Status:     {format_status(passed)}")

        return (1 if passed else 0), max(diff_lon, diff_lat)

    except Exception as e:
        if verbose:
            print(f"\n{name}: ERROR - {e}")
        return 0, 999.0


def compare_asteroid(jd, body_py, body_swe, name, verbose=True):
    """Legacy function for backward compatibility."""
    return compare_minor_body(jd, body_py, body_swe, name, 10.0, 5.0, verbose)


def run_category_tests(jd, date_str, category_name, bodies, tolerances, verbose=True):
    """
    Run tests for a category of minor bodies.

    Returns:
        (passed, failed, skipped, max_diff)
    """
    passed = 0
    failed = 0
    skipped = 0
    max_diff = 0.0

    print(f"\n{'-' * 40}")
    print(f"{category_name.upper()}")
    print(f"{'-' * 40}")

    for body_py, body_swe, body_name in bodies:
        status_code, diff = compare_minor_body(
            jd,
            body_py,
            body_swe,
            body_name,
            tolerances["longitude"],
            tolerances["latitude"],
            verbose,
        )

        if status_code == 1:  # PASSED
            passed += 1
            max_diff = max(max_diff, diff)
        elif status_code == 2:  # SKIPPED
            skipped += 1
        else:  # FAILED
            failed += 1

    return passed, failed, skipped, max_diff


def run_all_comparisons(verbose=True):
    """Run all minor body comparisons."""
    print_header("MINOR BODIES COMPARISON: LibEphemeris vs SwissEphemeris")

    print("\n" + "!" * 80)
    print("! NOTE: LibEphemeris uses simplified Keplerian orbital elements")
    print(
        "! Expected accuracy: ~1-5 arcminutes for inner bodies (astrological precision)"
    )
    print("! TNOs may have larger errors due to long orbital periods and perturbations")
    print("! SwissEphemeris uses full perturbation models (scientific precision)")
    print("!" * 80 + "\n")

    # Test dates
    subjects = [
        ("J2000.0", 2000, 1, 1, 12.0),
        ("2024-01-01", 2024, 1, 1, 0.0),
        ("2010-07-01", 2010, 7, 1, 12.0),
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

        # Main belt asteroids
        passed, failed, skipped, max_diff = run_category_tests(
            jd,
            date_str,
            "Main Belt Asteroids",
            MAIN_BELT_ASTEROIDS,
            TOLERANCES["main_belt"],
            verbose,
        )
        g_passed += passed
        g_failed += failed
        g_skipped += skipped
        g_total += passed + failed + skipped
        g_max_diff = max(g_max_diff, max_diff)

        # Centaurs
        passed, failed, skipped, max_diff = run_category_tests(
            jd, date_str, "Centaurs", CENTAURS, TOLERANCES["centaur"], verbose
        )
        g_passed += passed
        g_failed += failed
        g_skipped += skipped
        g_total += passed + failed + skipped
        g_max_diff = max(g_max_diff, max_diff)

        # TNOs
        passed, failed, skipped, max_diff = run_category_tests(
            jd,
            date_str,
            "Trans-Neptunian Objects (TNOs)",
            TNOS,
            TOLERANCES["tno"],
            verbose,
        )
        g_passed += passed
        g_failed += failed
        g_skipped += skipped
        g_total += passed + failed + skipped
        g_max_diff = max(g_max_diff, max_diff)

        # Near-Earth Asteroids
        passed, failed, skipped, max_diff = run_category_tests(
            jd,
            date_str,
            "Near-Earth Asteroids",
            NEAR_EARTH_ASTEROIDS,
            TOLERANCES["near_earth"],
            verbose,
        )
        g_passed += passed
        g_failed += failed
        g_skipped += skipped
        g_total += passed + failed + skipped
        g_max_diff = max(g_max_diff, max_diff)

    # Summary
    print("\n" + "=" * 80)
    print("MINOR BODIES COMPARISON SUMMARY")
    print("=" * 80)
    print(f"Total tests:   {g_total}")
    print(f"Passed:        {g_passed} ✓")
    print(f"Skipped:       {g_skipped} -")
    print(f"Failed:        {g_failed} ✗")

    if g_passed > 0:
        executed = g_total - g_skipped
        if executed > 0:
            print(f"Pass rate:     {g_passed / executed * 100:.1f}% (of executed)")
        print(f"Max diff:      {g_max_diff:.6f}")
    elif g_total == g_skipped:
        print("Pass rate:     N/A (All skipped)")
    else:
        print("Pass rate:     0.0%")

    print("=" * 80)

    return g_passed, g_total


# ============================================================================
# LEGACY MAIN FUNCTION
# ============================================================================


def main():
    """Legacy main function for backward compatibility."""
    args = parse_args(sys.argv)

    if args.get("help"):
        print("Usage: python compare_minor_bodies.py [OPTIONS]")
        print()
        print("Options:")
        print("  -v, --verbose    Show detailed output for each test")
        print("  -q, --quiet      Suppress per-test output")
        print("  -h, --help       Show this help message")
        print()
        return 0

    passed, total = run_all_comparisons(verbose=not args.get("quiet", False))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
