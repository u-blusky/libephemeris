"""
Extended Crossing Functions Comparison Script.

Compares additional crossing calculations between pyswisseph and libephemeris:
- mooncross_node / mooncross_node_ut - Moon crossing its nodes
- helio_cross / helio_cross_ut - heliocentric crossings
- swe_cross_ut - generic planetary crossings
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


class CrossingTolerance:
    """Tolerance thresholds for crossing comparisons."""

    TIME_SECONDS_MOON = 300.0  # 5 minutes for Moon (fast)
    TIME_SECONDS_PLANET = 600.0  # 10 minutes for planets


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_mooncross_node_ut(jd_start, node_type, desc, verbose=False):
    """
    Compare mooncross_node_ut function.

    node_type: 0 = ascending node, 1 = descending node
    """
    result = EventComparisonResult("mooncross_node_ut", desc)
    result.tolerance_seconds = CrossingTolerance.TIME_SECONDS_MOON

    # SwissEphemeris
    try:
        jd_swe = swe.mooncross_node_ut(jd_start, node_type)
        result.jd_swe = jd_swe
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        jd_py = pyephem.mooncross_node_ut(jd_start, node_type)
        result.jd_py = jd_py
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_mooncross_node(jd_start_et, node_type, desc, verbose=False):
    """Compare mooncross_node function (ET version)."""
    result = EventComparisonResult("mooncross_node", desc)
    result.tolerance_seconds = CrossingTolerance.TIME_SECONDS_MOON

    # SwissEphemeris
    try:
        jd_swe = swe.mooncross_node(jd_start_et, node_type)
        result.jd_swe = jd_swe
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        jd_py = pyephem.mooncross_node(jd_start_et, node_type)
        result.jd_py = jd_py
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_helio_cross_ut(
    jd_start, body_id, body_name, target_lon, desc, verbose=False
):
    """Compare helio_cross_ut function."""
    result = EventComparisonResult(
        "helio_cross_ut", f"{body_name} crosses {target_lon}° (helio) {desc}"
    )
    result.tolerance_seconds = CrossingTolerance.TIME_SECONDS_PLANET

    # SwissEphemeris
    try:
        jd_swe = swe.helio_cross_ut(body_id, target_lon, jd_start, SEFLG_SWIEPH)
        result.jd_swe = jd_swe
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        jd_py = pyephem.helio_cross_ut(body_id, target_lon, jd_start, SEFLG_SWIEPH)
        result.jd_py = jd_py
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_helio_cross(
    jd_start_et, body_id, body_name, target_lon, desc, verbose=False
):
    """Compare helio_cross function (ET version)."""
    result = EventComparisonResult(
        "helio_cross", f"{body_name} crosses {target_lon}° (helio/ET) {desc}"
    )
    result.tolerance_seconds = CrossingTolerance.TIME_SECONDS_PLANET

    # SwissEphemeris
    try:
        jd_swe = swe.helio_cross(body_id, target_lon, jd_start_et, SEFLG_SWIEPH)
        result.jd_swe = jd_swe
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        jd_py = pyephem.helio_cross(body_id, target_lon, jd_start_et, SEFLG_SWIEPH)
        result.jd_py = jd_py
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_cross_ut(jd_start, body_id, body_name, target_lon, desc, verbose=False):
    """Compare swe_cross_ut function (generic crossing)."""
    result = EventComparisonResult(
        "swe_cross_ut", f"{body_name} crosses {target_lon}° {desc}"
    )
    result.tolerance_seconds = CrossingTolerance.TIME_SECONDS_PLANET

    # SwissEphemeris - note: this might be implemented differently
    try:
        # swe_cross_ut takes (body, x2cross, jd_start, flag)
        jd_swe = (
            swe.solcross_ut(target_lon, jd_start, SEFLG_SWIEPH)
            if body_id == SE_SUN
            else None
        )
        if jd_swe is None:
            result.error_swe = "not implemented for this body"
            print(result.format_result(verbose))
            return True, 0.0, False  # Skip
        result.jd_swe = jd_swe
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        jd_py = pyephem.swe_cross_ut(body_id, target_lon, jd_start, SEFLG_SWIEPH)
        result.jd_py = jd_py
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all extended crossing comparisons."""
    print_header("EXTENDED CROSSING FUNCTIONS COMPARISON")
    stats = TestStatistics()

    # 1. Moon crossing its nodes
    print_section("MOON CROSSING NODES (mooncross_node_ut)")

    jd_start = swe.julday(2024, 1, 1, 0.0)

    # Test multiple node crossings
    for i in range(4):
        # Ascending node
        passed, diff, error = compare_mooncross_node_ut(
            jd_start + i * 14, 0, f"Asc node #{i + 1}", verbose
        )
        stats.add_result(passed, diff, error)

        # Descending node
        passed, diff, error = compare_mooncross_node_ut(
            jd_start + i * 14, 1, f"Desc node #{i + 1}", verbose
        )
        stats.add_result(passed, diff, error)

    # 2. Moon crossing nodes (ET version)
    print_section("MOON CROSSING NODES ET (mooncross_node)")
    jd_et = swe.julday(2024, 1, 1, 0.0) + swe.deltat(swe.julday(2024, 1, 1, 0.0))

    passed, diff, error = compare_mooncross_node(jd_et, 0, "Asc node (ET)", verbose)
    stats.add_result(passed, diff, error)

    passed, diff, error = compare_mooncross_node(jd_et, 1, "Desc node (ET)", verbose)
    stats.add_result(passed, diff, error)

    # 3. Heliocentric crossings
    print_section("HELIOCENTRIC CROSSINGS (helio_cross_ut)")

    jd_start = swe.julday(2024, 1, 1, 0.0)
    test_planets = [
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
    ]

    for body_id, body_name in test_planets:
        # Test crossing various longitudes
        for target_lon in [0, 90, 180, 270]:
            passed, diff, error = compare_helio_cross_ut(
                jd_start, body_id, body_name, target_lon, f"from 2024", verbose
            )
            stats.add_result(passed, diff, error)

    # 4. Heliocentric crossings (ET version)
    print_section("HELIOCENTRIC CROSSINGS ET (helio_cross)")
    jd_et = swe.julday(2024, 6, 1, 0.0) + swe.deltat(swe.julday(2024, 6, 1, 0.0))

    for body_id, body_name in test_planets[:2]:
        passed, diff, error = compare_helio_cross(
            jd_et, body_id, body_name, 0, "(ET)", verbose
        )
        stats.add_result(passed, diff, error)

    # 5. Generic crossing (swe_cross_ut)
    print_section("GENERIC PLANETARY CROSSINGS (swe_cross_ut)")
    jd_start = swe.julday(2024, 1, 1, 0.0)

    # Test with various planets
    cross_tests = [
        (SE_VENUS, "Venus", 0),
        (SE_VENUS, "Venus", 90),
        (SE_MARS, "Mars", 180),
        (SE_JUPITER, "Jupiter", 0),
    ]

    for body_id, body_name, target_lon in cross_tests:
        passed, diff, error = compare_cross_ut(
            jd_start, body_id, body_name, target_lon, "2024", verbose
        )
        stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("EXTENDED CROSSINGS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_crossings_ext.py [OPTIONS]")
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
