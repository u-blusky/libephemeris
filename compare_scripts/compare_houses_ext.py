"""
Extended House Functions Comparison Script.

Compares additional house-related calculations between pyswisseph and libephemeris:
- houses_armc - houses from ARMC
- houses_armc_ex2 - extended houses from ARMC
- houses_ex / houses_ex2 - houses with sidereal mode
- house_pos - position of planet in house
- gauquelin_sector - Gauquelin sector calculations
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

HOUSE_SYSTEMS = ["P", "K", "R", "C", "E", "W", "M", "B"]

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
