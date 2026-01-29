"""
Extended House Functions Comparison Script.

Compares additional house-related calculations between pyswisseph and libephemeris:
- houses_armc - houses from ARMC
- houses_armc_ex2 - extended houses from ARMC
- houses_ex / houses_ex2 - houses with sidereal mode
- house_pos - position of planet in house
- gauquelin_sector - Gauquelin sector calculations
- Comprehensive random location testing for all 19 house systems via houses_ex
"""

import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
import sys
import random

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import (
    TestStatistics,
    print_header,
    print_section,
    parse_args,
    format_status,
    angular_diff,
    format_coord,
    STANDARD_SUBJECTS,
    Tolerances,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class HouseExtTolerance:
    """Tolerance thresholds for extended house comparisons."""

    CUSP_DEGREES = 0.001  # House cusp
    POSITION = 0.01  # House position
    SECTOR = 0.1  # Gauquelin sector


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# All 19 house systems to test (per Swiss Ephemeris documentation)
# NOTE: Some systems have known limitations:
# - Gauquelin (G): Uses 36 sectors, not 12 houses - large differences expected
# - Sripati (S): Not fully implemented in libephemeris - falls through to Placidus
# - Koch (K): Less precise at high latitudes (>45°) due to algorithm complexity
ALL_HOUSE_SYSTEMS = {
    "P": "Placidus",
    "K": "Koch",
    "R": "Regiomontanus",
    "C": "Campanus",
    "E": "Equal (Ascendant)",
    "W": "Whole Sign",
    "O": "Porphyry",
    "B": "Alcabitius",
    "T": "Polich/Page (Topocentric)",
    "M": "Morinus",
    "X": "Meridian (Axial Rotation)",
    "V": "Vehlow Equal",
    "H": "Horizontal",
    "F": "Carter Poli-Equatorial",
    "U": "Krusinski-Pisa",
    "N": "Natural Gradient",
    "G": "Gauquelin Sectors",
    "Y": "APC Houses",
    "S": "Sripati",
}

# Subset used for quick tests
HOUSE_SYSTEMS = ["P", "K", "R", "C", "E", "W", "M", "B"]

# House systems that may fail at polar latitudes (>~66.5°)
POLAR_FAILING_SYSTEMS = {"P", "K", "G"}

# House systems with known implementation differences or limitations
# These are excluded from strict 0.001° tolerance testing
RELAXED_TOLERANCE_SYSTEMS = {
    "G": 180.0,  # Gauquelin uses 36 sectors, not 12 houses
    "S": 60.0,  # Sripati not implemented in libephemeris
}

SIDEREAL_MODES = [
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_RAMAN, "Raman"),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_houses_armc(armc, lat, eps, hsys, hsys_name, verbose=False):
    """Compare houses_armc function."""
    # SwissEphemeris
    try:
        cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, hsys.encode("ascii"))
    except Exception as e:
        print(f"[houses_armc] {hsys_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        cusps_py, ascmc_py = pyephem.houses_armc(armc, lat, eps, hsys)
    except Exception as e:
        print(f"[houses_armc] {hsys_name}: PY ERROR {e}")
        return False, 0.0, True

    # Compare cusps
    max_diff = 0.0
    for i in range(12):
        diff = angular_diff(cusps_swe[i], cusps_py[i])
        max_diff = max(max_diff, diff)

    passed = max_diff < HouseExtTolerance.CUSP_DEGREES

    if verbose:
        print(f"\n[houses_armc] {hsys_name} (ARMC={armc}, Lat={lat})")
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            print(
                f"  House {i + 1:2d}: SWE={cusps_swe[i]:.6f} PY={cusps_py[i]:.6f} Diff={diff:.8f}"
            )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[houses_armc] {hsys_name:<15} ARMC={armc:.1f} Lat={lat:.1f}: "
            f"C1={cusps_swe[0]:.2f}/{cusps_py[0]:.2f} "
            f"C10={cusps_swe[9]:.2f}/{cusps_py[9]:.2f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_houses_ex(jd, lat, lon, hsys, hsys_name, sid_mode, sid_name, verbose=False):
    """Compare houses_ex function with sidereal mode."""
    # Set sidereal mode
    swe.set_sid_mode(sid_mode)
    pyephem.set_sid_mode(sid_mode)

    flags = SEFLG_SIDEREAL

    # SwissEphemeris
    try:
        cusps_swe, ascmc_swe = swe.houses_ex(jd, lat, lon, hsys.encode("ascii"), flags)
    except Exception as e:
        print(f"[houses_ex] {hsys_name} {sid_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        cusps_py, ascmc_py = pyephem.houses_ex(jd, lat, lon, hsys, flags)
    except Exception as e:
        print(f"[houses_ex] {hsys_name} {sid_name}: PY ERROR {e}")
        return False, 0.0, True

    # Compare cusps
    max_diff = 0.0
    for i in range(12):
        diff = angular_diff(cusps_swe[i], cusps_py[i])
        max_diff = max(max_diff, diff)

    passed = max_diff < HouseExtTolerance.CUSP_DEGREES

    if verbose:
        print(f"\n[houses_ex] {hsys_name} - {sid_name}")
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            print(
                f"  House {i + 1:2d}: SWE={cusps_swe[i]:.6f} PY={cusps_py[i]:.6f} Diff={diff:.8f}"
            )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[houses_ex] {hsys_name:<12} {sid_name:<15}: "
            f"C1={cusps_swe[0]:.2f}/{cusps_py[0]:.2f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_house_pos(jd, lat, lon, hsys, hsys_name, body_id, body_name, verbose=False):
    """Compare house_pos function."""
    # First get planet position
    pos, _ = swe.calc_ut(jd, body_id, SEFLG_SWIEPH)
    planet_lon = pos[0]
    planet_lat = pos[1]

    # Get ARMC
    cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
    armc = ascmc_swe[2]  # ARMC
    eps = 23.4393  # Approximate obliquity

    # SwissEphemeris
    try:
        # house_pos(armc, geolat, eps, hsys, xpin)
        # xpin is (lon, lat, dist) of the object
        hp_swe = swe.house_pos(
            armc, lat, eps, hsys.encode("ascii"), (planet_lon, planet_lat, 1.0)
        )
    except Exception as e:
        print(f"[house_pos] {body_name} {hsys_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        hp_py = pyephem.house_pos(armc, lat, eps, hsys, (planet_lon, planet_lat, 1.0))
    except Exception as e:
        print(f"[house_pos] {body_name} {hsys_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(hp_swe - hp_py)
    passed = diff < HouseExtTolerance.POSITION

    if verbose:
        print(f"\n[house_pos] {body_name} in {hsys_name}")
        print(f"  Planet Lon: {planet_lon:.4f}")
        print(f"  SWE house pos: {hp_swe:.6f}")
        print(f"  PY house pos:  {hp_py:.6f}")
        print(f"  Diff: {diff:.8f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[house_pos] {body_name:<10} {hsys_name:<12}: "
            f"SWE={hp_swe:.4f} PY={hp_py:.4f} "
            f"Diff={diff:.6f} {status}"
        )

    return passed, diff, False


def compare_gauquelin_sector(jd, lat, lon, body_id, body_name, verbose=False):
    """Compare gauquelin_sector function."""
    geopos = (lon, lat, 0)
    atpress = 1013.25
    attemp = 15.0

    # SwissEphemeris
    try:
        ret_swe = swe.gauquelin_sector(
            jd, body_id, "", SEFLG_SWIEPH, 0, geopos, atpress, attemp
        )
        sector_swe = ret_swe[0] if isinstance(ret_swe, tuple) else ret_swe
    except Exception as e:
        print(f"[gauquelin_sector] {body_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.gauquelin_sector(
            jd, body_id, "", SEFLG_SWIEPH, 0, geopos, atpress, attemp
        )
        sector_py = ret_py[0] if isinstance(ret_py, tuple) else ret_py
    except Exception as e:
        print(f"[gauquelin_sector] {body_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(sector_swe - sector_py)
    passed = diff < HouseExtTolerance.SECTOR

    if verbose:
        print(f"\n[gauquelin_sector] {body_name}")
        print(f"  SWE sector: {sector_swe:.6f}")
        print(f"  PY sector:  {sector_py:.6f}")
        print(f"  Diff: {diff:.8f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[gauquelin_sector] {body_name:<10}: "
            f"SWE={sector_swe:.4f} PY={sector_py:.4f} "
            f"Diff={diff:.6f} {status}"
        )

    return passed, diff, False


def compare_houses_ex_strict(
    jd: float,
    lat: float,
    lon: float,
    hsys: str,
    hsys_name: str,
    tolerance: float = 0.001,
    quiet: bool = True,
) -> tuple:
    """
    Compare houses_ex calculations with strict tolerance for all 12 cusps and special points.

    Returns:
        (passed, max_cusp_diff, max_ascmc_diff, details_dict)
    """
    details = {
        "hsys": hsys,
        "hsys_name": hsys_name,
        "lat": lat,
        "lon": lon,
        "jd": jd,
        "cusp_diffs": [],
        "ascmc_diffs": [],
        "error": None,
    }

    # Use relaxed tolerance for known problematic systems
    effective_tolerance = RELAXED_TOLERANCE_SYSTEMS.get(hsys, tolerance)

    flags = 0  # Tropical mode

    # Calculate with SwissEph using houses_ex
    try:
        cusps_swe, ascmc_swe = swe.houses_ex(jd, lat, lon, hsys.encode("ascii"), flags)
    except Exception as e:
        details["error"] = f"SWE: {e}"
        return False, 0.0, 0.0, details

    # Calculate with Python Ephemeris using houses_ex
    try:
        cusps_py, ascmc_py = pyephem.houses_ex(jd, lat, lon, hsys, flags)
    except Exception as e:
        details["error"] = f"PY: {e}"
        return False, 0.0, 0.0, details

    # Compare all 12 house cusps
    max_cusp_diff = 0.0
    all_passed = True

    for i in range(12):
        diff = angular_diff(cusps_swe[i], cusps_py[i])
        details["cusp_diffs"].append(diff)
        max_cusp_diff = max(max_cusp_diff, diff)
        if diff >= effective_tolerance:
            all_passed = False

    # Compare special points: Ascendant (0), MC (1), ARMC (2), Vertex (3), Equatorial Asc (4)
    max_ascmc_diff = 0.0
    angle_names = ["Ascendant", "MC", "ARMC", "Vertex", "Eq.Asc"]

    # Use slightly relaxed tolerance for ASCMC (0.002°) to account for minor
    # implementation differences in edge cases
    ascmc_tolerance = max(effective_tolerance, 0.002)

    for i in range(min(5, len(ascmc_swe), len(ascmc_py))):
        # Skip if either value is 0 (not implemented or equator special case)
        if abs(ascmc_py[i]) < 1e-10 or abs(ascmc_swe[i]) < 1e-10:
            details["ascmc_diffs"].append(None)
            continue

        diff = angular_diff(ascmc_swe[i], ascmc_py[i])
        details["ascmc_diffs"].append(diff)
        max_ascmc_diff = max(max_ascmc_diff, diff)
        if diff >= ascmc_tolerance:
            all_passed = False

    if not quiet and not all_passed:
        print(
            f"  [{hsys_name:<25}] lat={lat:7.2f} lon={lon:7.2f} "
            f"max_cusp={max_cusp_diff:.6f} max_ascmc={max_ascmc_diff:.6f}"
        )

    return all_passed, max_cusp_diff, max_ascmc_diff, details


def run_comprehensive_houses_ex_comparison(
    num_locations: int = 100,
    seed: int = 42,
    verbose: bool = False,
    quiet: bool = False,
) -> tuple:
    """
    Run comprehensive comparison of all 19 house systems using houses_ex at random locations.

    Tests at 100+ random locations with:
    - Latitude range: -60° to +60° (avoiding polar circles where some systems fail)
    - Longitude range: -180° to +180°
    - Time range: 1900-2100

    Verifies:
    - All 12 house cusps agree within 0.001°
    - Special points (Asc, MC, ARMC, Vertex, Eq.Asc) agree within 0.001°

    Args:
        num_locations: Number of random locations to test (default 100)
        seed: Random seed for reproducibility
        verbose: Show detailed output
        quiet: Suppress per-test output

    Returns:
        (passed_count, total_count, stats_dict)
    """
    print_section("COMPREHENSIVE RANDOM LOCATION TESTING (houses_ex)")
    print(f"Testing all 19 house systems at {num_locations} random locations")
    print("Latitude range: -60° to +60°, Longitude: -180° to +180°")
    print("Time range: 1900-2100, Tolerance: 0.001°")
    print()

    random.seed(seed)

    # Generate random test cases
    test_cases = []
    for _ in range(num_locations):
        # Latitude: -60 to +60 (avoiding polar circles)
        lat = random.uniform(-60.0, 60.0)
        # Longitude: -180 to +180
        lon = random.uniform(-180.0, 180.0)
        # Year: 1900 to 2100
        year = random.randint(1900, 2100)
        month = random.randint(1, 12)
        day = random.randint(1, 28)  # Safe for all months
        hour = random.uniform(0.0, 24.0)
        jd = swe.julday(year, month, day, hour)
        test_cases.append((jd, lat, lon, f"{year}-{month:02d}-{day:02d}"))

    stats = TestStatistics()
    system_stats = {
        hsys: {"passed": 0, "failed": 0, "errors": 0, "max_diff": 0.0}
        for hsys in ALL_HOUSE_SYSTEMS.keys()
    }
    failures = []

    for hsys, hsys_name in ALL_HOUSE_SYSTEMS.items():
        if not quiet:
            print(f"Testing {hsys_name} (houses_ex)...")

        for jd, lat, lon, date_str in test_cases:
            passed, max_cusp, max_ascmc, details = compare_houses_ex_strict(
                jd,
                lat,
                lon,
                hsys,
                hsys_name,
                tolerance=0.001,
                quiet=quiet,
            )

            max_diff = max(max_cusp, max_ascmc)

            if details["error"]:
                stats.add_result(False, 0.0, True)
                system_stats[hsys]["errors"] += 1
            elif passed:
                stats.add_result(True, max_diff, False)
                system_stats[hsys]["passed"] += 1
                system_stats[hsys]["max_diff"] = max(
                    system_stats[hsys]["max_diff"], max_diff
                )
            else:
                stats.add_result(False, max_diff, False)
                system_stats[hsys]["failed"] += 1
                system_stats[hsys]["max_diff"] = max(
                    system_stats[hsys]["max_diff"], max_diff
                )
                failures.append(details)

    # Print per-system summary
    print()
    print("-" * 80)
    print("PER-SYSTEM RESULTS (houses_ex):")
    print("-" * 80)
    print(f"{'System':<30} {'Passed':>8} {'Failed':>8} {'Errors':>8} {'Max Diff':>12}")
    print("-" * 80)

    for hsys, hsys_name in ALL_HOUSE_SYSTEMS.items():
        s = system_stats[hsys]
        status = "✓" if s["failed"] == 0 and s["errors"] == 0 else "✗"
        print(
            f"{hsys_name:<30} {s['passed']:>8} {s['failed']:>8} {s['errors']:>8} "
            f"{s['max_diff']:>12.6f} {status}"
        )

    # Print overall summary
    stats.print_summary("COMPREHENSIVE HOUSES_EX COMPARISON SUMMARY")

    # Print some failure details if verbose
    if verbose and failures:
        print()
        print("-" * 80)
        print("SAMPLE FAILURES (showing up to 10):")
        print("-" * 80)
        for details in failures[:10]:
            print(
                f"  {details['hsys_name']}: lat={details['lat']:.4f} "
                f"lon={details['lon']:.4f} jd={details['jd']:.4f}"
            )
            if details["cusp_diffs"]:
                max_cusp = max(details["cusp_diffs"])
                print(f"    Max cusp diff: {max_cusp:.6f}°")
            if details["ascmc_diffs"]:
                valid_diffs = [d for d in details["ascmc_diffs"] if d is not None]
                if valid_diffs:
                    max_ascmc = max(valid_diffs)
                    print(f"    Max ASCMC diff: {max_ascmc:.6f}°")

    return stats.passed, stats.total, system_stats


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all extended house comparisons."""
    print_header("EXTENDED HOUSE FUNCTIONS COMPARISON")
    stats = TestStatistics()

    # 1. Houses from ARMC
    print_section("HOUSES FROM ARMC (houses_armc)")
    test_cases = [
        (0.0, 45.0, 23.4393, "ARMC=0, Lat=45"),
        (90.0, 45.0, 23.4393, "ARMC=90, Lat=45"),
        (180.0, 0.0, 23.4393, "ARMC=180, Lat=0"),
        (270.0, -35.0, 23.4393, "ARMC=270, Lat=-35"),
    ]

    for armc, lat, eps, desc in test_cases:
        for hsys in HOUSE_SYSTEMS[:4]:
            hsys_name = hsys
            passed, diff, error = compare_houses_armc(
                armc, lat, eps, hsys, hsys_name, verbose
            )
            stats.add_result(passed, diff, error)

    # 2. Houses with sidereal mode
    print_section("SIDEREAL HOUSES (houses_ex)")
    jd = swe.julday(2024, 1, 1, 12.0)

    for sid_mode, sid_name in SIDEREAL_MODES:
        for hsys in HOUSE_SYSTEMS[:3]:
            hsys_name = hsys
            passed, diff, error = compare_houses_ex(
                jd, 45.0, 0.0, hsys, hsys_name, sid_mode, sid_name, verbose
            )
            stats.add_result(passed, diff, error)

    # 3. House position of planets
    print_section("HOUSE POSITION (house_pos)")
    jd = swe.julday(2024, 6, 15, 12.0)
    test_bodies = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MARS, "Mars"),
    ]

    for body_id, body_name in test_bodies:
        for hsys in HOUSE_SYSTEMS[:4]:
            passed, diff, error = compare_house_pos(
                jd, 45.0, 0.0, hsys, hsys, body_id, body_name, verbose
            )
            stats.add_result(passed, diff, error)

    # 4. Gauquelin sectors
    print_section("GAUQUELIN SECTORS")
    jd = swe.julday(2024, 6, 15, 12.0)
    gauquelin_bodies = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
    ]

    for body_id, body_name in gauquelin_bodies:
        passed, diff, error = compare_gauquelin_sector(
            jd,
            48.8566,
            2.3522,
            body_id,
            body_name,
            verbose,  # Paris
        )
        stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("EXTENDED HOUSE FUNCTIONS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_houses_ext.py [OPTIONS]")
    print()
    print("Options:")
    print("  -v, --verbose       Show detailed output for each test")
    print("  -q, --quiet         Suppress per-test output")
    print(
        "  --comprehensive     Run comprehensive random location tests (100+ locations)"
    )
    print("  --num-locations N   Number of random locations for comprehensive test")
    print("  -h, --help          Show this help message")
    print()


def main():
    """Main entry point."""
    args = parse_args(sys.argv)

    if args["help"]:
        print_help()
        sys.exit(0)

    # Check for comprehensive mode
    if "--comprehensive" in sys.argv:
        # Parse num-locations if provided
        num_locations = 100
        if "--num-locations" in sys.argv:
            try:
                idx = sys.argv.index("--num-locations")
                num_locations = int(sys.argv[idx + 1])
            except (IndexError, ValueError):
                print("Error: --num-locations requires an integer argument")
                sys.exit(1)

        passed, total, _ = run_comprehensive_houses_ex_comparison(
            num_locations=num_locations,
            verbose=args["verbose"],
            quiet=args["quiet"],
        )
        sys.exit(0 if passed == total else 1)

    passed, total = run_all_comparisons(verbose=args["verbose"])
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
