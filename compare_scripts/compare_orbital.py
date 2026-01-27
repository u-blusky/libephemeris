"""
Orbital Elements Comparison Script.

Compares orbital element calculations between pyswisseph and libephemeris:
- nod_aps / nod_aps_ut - nodes and apsides
- get_orbital_elements / get_orbital_elements_ut - Keplerian elements
- orbit_max_min_true_distance - orbital distance extremes
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
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class OrbitalTolerance:
    """Tolerance thresholds for orbital element comparisons."""

    ANGLE_DEGREES = 0.01  # 0.01 degree for angular elements
    DISTANCE_AU = 0.001  # 0.001 AU for distances
    ECCENTRICITY = 0.0001  # Eccentricity (dimensionless)
    PERIOD_DAYS = 0.1  # Days for orbital period


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for orbital element tests
TEST_PLANETS = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Method types for nod_aps
NOD_APS_METHODS = [
    (0, "Mean"),
    (1, "Osculating"),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_nod_aps_ut(jd, body_id, body_name, method, method_name, verbose=False):
    """Compare nod_aps_ut function."""
    # SwissEphemeris
    try:
        ret_swe = swe.nod_aps_ut(jd, body_id, SEFLG_SWIEPH, method)
        # Returns (nasc, ndesc, peri, aphe) each with 6 elements
        nasc_swe = ret_swe[0]
        ndesc_swe = ret_swe[1]
        peri_swe = ret_swe[2]
        aphe_swe = ret_swe[3]
    except Exception as e:
        print(f"[nod_aps_ut] {body_name} {method_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.nod_aps_ut(jd, body_id, SEFLG_SWIEPH, method)
        nasc_py = ret_py[0]
        ndesc_py = ret_py[1]
        peri_py = ret_py[2]
        aphe_py = ret_py[3]
    except Exception as e:
        print(f"[nod_aps_ut] {body_name} {method_name}: PY ERROR {e}")
        return False, 0.0, True

    # Compare ascending node longitude
    diff_nasc = angular_diff(nasc_swe[0], nasc_py[0])
    diff_ndesc = angular_diff(ndesc_swe[0], ndesc_py[0])
    diff_peri = angular_diff(peri_swe[0], peri_py[0])
    diff_aphe = angular_diff(aphe_swe[0], aphe_py[0])

    max_diff = max(diff_nasc, diff_ndesc, diff_peri, diff_aphe)
    passed = max_diff < OrbitalTolerance.ANGLE_DEGREES

    if verbose:
        print(f"\n[nod_aps_ut] {body_name} - {method_name}")
        print(
            f"  Asc Node:   SWE={nasc_swe[0]:.6f} PY={nasc_py[0]:.6f} Diff={diff_nasc:.8f}"
        )
        print(
            f"  Desc Node:  SWE={ndesc_swe[0]:.6f} PY={ndesc_py[0]:.6f} Diff={diff_ndesc:.8f}"
        )
        print(
            f"  Perihelion: SWE={peri_swe[0]:.6f} PY={peri_py[0]:.6f} Diff={diff_peri:.8f}"
        )
        print(
            f"  Aphelion:   SWE={aphe_swe[0]:.6f} PY={aphe_py[0]:.6f} Diff={diff_aphe:.8f}"
        )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[nod_aps_ut] {body_name:<10} {method_name:<10}: "
            f"AscN={nasc_swe[0]:.2f}/{nasc_py[0]:.2f} "
            f"Peri={peri_swe[0]:.2f}/{peri_py[0]:.2f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_get_orbital_elements_ut(jd, body_id, body_name, verbose=False):
    """Compare get_orbital_elements_ut function."""
    # SwissEphemeris
    try:
        ret_swe = swe.get_orbital_elements(jd, body_id, SEFLG_SWIEPH)
        # Returns tuple of 50 orbital elements
        # Key elements: [0] semi-major axis, [1] eccentricity, [2] inclination
        # [3] asc node, [4] arg perihelion, [5] perihelion distance
        a_swe = ret_swe[0]  # Semi-major axis
        e_swe = ret_swe[1]  # Eccentricity
        i_swe = ret_swe[2]  # Inclination
        node_swe = ret_swe[3]  # Ascending node
        peri_swe = ret_swe[4]  # Argument of perihelion
    except Exception as e:
        print(f"[get_orbital_elements] {body_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.get_orbital_elements(jd, body_id, SEFLG_SWIEPH)
        a_py = ret_py[0]
        e_py = ret_py[1]
        i_py = ret_py[2]
        node_py = ret_py[3]
        peri_py = ret_py[4]
    except Exception as e:
        print(f"[get_orbital_elements] {body_name}: PY ERROR {e}")
        return False, 0.0, True

    # Calculate differences
    diff_a = abs(a_swe - a_py)
    diff_e = abs(e_swe - e_py)
    diff_i = abs(i_swe - i_py)
    diff_node = angular_diff(node_swe, node_py)
    diff_peri = angular_diff(peri_swe, peri_py)

    # Check tolerances
    passed_a = diff_a < OrbitalTolerance.DISTANCE_AU
    passed_e = diff_e < OrbitalTolerance.ECCENTRICITY
    passed_i = diff_i < OrbitalTolerance.ANGLE_DEGREES
    passed_node = diff_node < OrbitalTolerance.ANGLE_DEGREES
    passed_peri = diff_peri < OrbitalTolerance.ANGLE_DEGREES

    passed = passed_a and passed_e and passed_i and passed_node and passed_peri
    max_diff = max(diff_a, diff_e * 100, diff_i, diff_node, diff_peri)

    if verbose:
        print(f"\n[get_orbital_elements] {body_name}")
        print(
            f"  Semi-major axis: SWE={a_swe:.8f} PY={a_py:.8f} Diff={diff_a:.10f} AU {format_status(passed_a)}"
        )
        print(
            f"  Eccentricity:    SWE={e_swe:.8f} PY={e_py:.8f} Diff={diff_e:.10f} {format_status(passed_e)}"
        )
        print(
            f"  Inclination:     SWE={i_swe:.8f} PY={i_py:.8f} Diff={diff_i:.10f} deg {format_status(passed_i)}"
        )
        print(
            f"  Asc Node:        SWE={node_swe:.6f} PY={node_py:.6f} Diff={diff_node:.8f} deg {format_status(passed_node)}"
        )
        print(
            f"  Arg Peri:        SWE={peri_swe:.6f} PY={peri_py:.6f} Diff={diff_peri:.8f} deg {format_status(passed_peri)}"
        )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[get_orbital_elements] {body_name:<10}: "
            f"a={a_swe:.4f}/{a_py:.4f} "
            f"e={e_swe:.6f}/{e_py:.6f} "
            f"i={i_swe:.4f}/{i_py:.4f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_orbit_max_min_true_distance(jd, body_id, body_name, verbose=False):
    """Compare orbit_max_min_true_distance function."""
    # SwissEphemeris
    try:
        ret_swe = swe.orbit_max_min_true_distance(jd, body_id, SEFLG_SWIEPH)
        # Returns (dmax, dmin, dtrue)
        dmax_swe = ret_swe[0]
        dmin_swe = ret_swe[1]
        dtrue_swe = ret_swe[2]
    except Exception as e:
        print(f"[orbit_max_min_true_distance] {body_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.orbit_max_min_true_distance(jd, body_id, SEFLG_SWIEPH)
        dmax_py = ret_py[0]
        dmin_py = ret_py[1]
        dtrue_py = ret_py[2]
    except Exception as e:
        print(f"[orbit_max_min_true_distance] {body_name}: PY ERROR {e}")
        return False, 0.0, True

    diff_max = abs(dmax_swe - dmax_py)
    diff_min = abs(dmin_swe - dmin_py)
    diff_true = abs(dtrue_swe - dtrue_py)
    max_diff = max(diff_max, diff_min, diff_true)

    passed = max_diff < OrbitalTolerance.DISTANCE_AU

    if verbose:
        print(f"\n[orbit_max_min_true_distance] {body_name}")
        print(
            f"  Max Distance:  SWE={dmax_swe:.8f} PY={dmax_py:.8f} Diff={diff_max:.10f} AU"
        )
        print(
            f"  Min Distance:  SWE={dmin_swe:.8f} PY={dmin_py:.8f} Diff={diff_min:.10f} AU"
        )
        print(
            f"  True Distance: SWE={dtrue_swe:.8f} PY={dtrue_py:.8f} Diff={diff_true:.10f} AU"
        )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[orbit_max_min_true_distance] {body_name:<10}: "
            f"Max={dmax_swe:.4f}/{dmax_py:.4f} "
            f"Min={dmin_swe:.4f}/{dmin_py:.4f} "
            f"True={dtrue_swe:.4f}/{dtrue_py:.4f} "
            f"MaxDiff={max_diff:.8f} {status}"
        )

    return passed, max_diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all orbital element comparisons."""
    print_header("ORBITAL ELEMENTS COMPARISON")
    stats = TestStatistics()

    # Test dates
    test_jds = [
        (swe.julday(2000, 1, 1, 12), "J2000.0"),
        (swe.julday(2024, 6, 15, 0), "2024"),
    ]

    # 1. Nodes and Apsides
    print_section("NODES AND APSIDES (nod_aps_ut)")
    for jd, date_desc in test_jds:
        for body_id, body_name in TEST_PLANETS:
            for method, method_name in NOD_APS_METHODS:
                passed, diff, error = compare_nod_aps_ut(
                    jd, body_id, body_name, method, method_name, verbose
                )
                stats.add_result(passed, diff, error)

    # 2. Orbital Elements
    print_section("ORBITAL ELEMENTS (get_orbital_elements)")
    for jd, date_desc in test_jds:
        for body_id, body_name in TEST_PLANETS:
            passed, diff, error = compare_get_orbital_elements_ut(
                jd, body_id, body_name, verbose
            )
            stats.add_result(passed, diff, error)

    # 3. Orbital Distance Extremes
    print_section("ORBITAL DISTANCE EXTREMES (orbit_max_min_true_distance)")
    jd = swe.julday(2024, 1, 1, 0)
    for body_id, body_name in TEST_PLANETS:
        passed, diff, error = compare_orbit_max_min_true_distance(
            jd, body_id, body_name, verbose
        )
        stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("ORBITAL ELEMENTS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_orbital.py [OPTIONS]")
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
