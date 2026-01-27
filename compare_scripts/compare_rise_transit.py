"""
Rise/Set/Transit Functions Comparison Script.

Compares rise_trans and rise_trans_true_hor functions between
pyswisseph and libephemeris for various bodies and locations.
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
    RISE_SET_LOCATIONS,
)


# ============================================================================
# CONSTANTS
# ============================================================================

# Rise/set calculation types
SE_CALC_RISE = 1
SE_CALC_SET = 2
SE_CALC_MTRANSIT = 4  # Upper culmination
SE_CALC_ITRANSIT = 8  # Lower culmination

# Bodies to test
TEST_BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]

# Event types to test
EVENT_TYPES = [
    (SE_CALC_RISE, "Rise"),
    (SE_CALC_SET, "Set"),
    (SE_CALC_MTRANSIT, "Transit"),
]


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class RiseSetTolerance:
    """Tolerance thresholds for rise/set comparisons."""

    TIME_SECONDS = 120.0  # 2 minutes
    TIME_SECONDS_MOON = 300.0  # 5 minutes for Moon (fast moving)


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_rise_trans(
    jd_start,
    body_id,
    body_name,
    event_type,
    event_name,
    lat,
    lon,
    alt,
    loc_name,
    verbose=False,
):
    """Compare rise_trans function for a single body/event/location."""
    result = EventComparisonResult(
        "rise_trans", f"{body_name} {event_name} @ {loc_name}"
    )
    result.tolerance_seconds = (
        RiseSetTolerance.TIME_SECONDS_MOON
        if body_id == SE_MOON
        else RiseSetTolerance.TIME_SECONDS
    )

    geopos = (lon, lat, alt)
    atpress = 1013.25  # Standard atmospheric pressure
    attemp = 15.0  # Standard temperature

    # SwissEphemeris
    try:
        ret_swe = swe.rise_trans(
            jd_start, body_id, "", SEFLG_SWIEPH, event_type, geopos, atpress, attemp
        )
        if ret_swe[0] == -2:  # Body never rises/sets
            result.error_swe = "never rises/sets"
            print(result.format_result(verbose))
            return True, 0.0, False  # Not an error, just skip
        result.jd_swe = ret_swe[1][0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.rise_trans(
            jd_start, body_id, "", SEFLG_SWIEPH, event_type, geopos, atpress, attemp
        )
        if ret_py[0] == -2:
            result.error_py = "never rises/sets"
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


def compare_rise_trans_true_hor(
    jd_start,
    body_id,
    body_name,
    event_type,
    event_name,
    lat,
    lon,
    alt,
    loc_name,
    hor_alt,
    verbose=False,
):
    """Compare rise_trans_true_hor function with true horizon altitude."""
    result = EventComparisonResult(
        "rise_trans_true_hor", f"{body_name} {event_name} @ {loc_name} (hor={hor_alt})"
    )
    result.tolerance_seconds = (
        RiseSetTolerance.TIME_SECONDS_MOON
        if body_id == SE_MOON
        else RiseSetTolerance.TIME_SECONDS
    )

    geopos = (lon, lat, alt)
    atpress = 1013.25
    attemp = 15.0

    # SwissEphemeris
    try:
        ret_swe = swe.rise_trans_true_hor(
            jd_start,
            body_id,
            "",
            SEFLG_SWIEPH,
            event_type,
            geopos,
            atpress,
            attemp,
            hor_alt,
        )
        if ret_swe[0] == -2:
            result.error_swe = "never rises/sets"
            print(result.format_result(verbose))
            return True, 0.0, False
        result.jd_swe = ret_swe[1][0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.rise_trans_true_hor(
            jd_start,
            body_id,
            "",
            SEFLG_SWIEPH,
            event_type,
            geopos,
            atpress,
            attemp,
            hor_alt,
        )
        if ret_py[0] == -2:
            result.error_py = "never rises/sets"
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


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all rise/set/transit comparisons."""
    print_header("RISE/SET/TRANSIT FUNCTIONS COMPARISON")
    stats = TestStatistics()

    # Test dates
    test_dates = [
        (2024, 6, 21, "Summer Solstice 2024"),
        (2024, 12, 21, "Winter Solstice 2024"),
        (2024, 3, 20, "Vernal Equinox 2024"),
    ]

    # 1. Basic rise_trans for Sun and Moon at standard locations
    print_section("RISE/SET/TRANSIT - SUN & MOON")
    for year, month, day, date_desc in test_dates[:1]:
        jd_start = swe.julday(year, month, day, 0.0)

        for body_id, body_name in TEST_BODIES[:2]:  # Sun and Moon
            for event_type, event_name in EVENT_TYPES:
                for loc_name, lat, lon, alt in RISE_SET_LOCATIONS[:3]:
                    passed, diff, error = compare_rise_trans(
                        jd_start,
                        body_id,
                        body_name,
                        event_type,
                        event_name,
                        lat,
                        lon,
                        alt,
                        loc_name,
                        verbose,
                    )
                    stats.add_result(passed, diff, error)

    # 2. Planets at mid-latitude
    print_section("RISE/SET/TRANSIT - PLANETS")
    jd_start = swe.julday(2024, 6, 15, 0.0)
    lat, lon, alt = 45.0, 0.0, 0  # Mid-latitude

    for body_id, body_name in TEST_BODIES[2:]:  # Venus, Mars, Jupiter
        for event_type, event_name in EVENT_TYPES:
            passed, diff, error = compare_rise_trans(
                jd_start,
                body_id,
                body_name,
                event_type,
                event_name,
                lat,
                lon,
                alt,
                "Mid-Latitude",
                verbose,
            )
            stats.add_result(passed, diff, error)

    # 3. High latitude tests (polar day/night edge cases)
    print_section("RISE/SET - HIGH LATITUDE")
    for loc_name, lat, lon, alt in RISE_SET_LOCATIONS[4:6]:  # Arctic locations
        for year, month, day, date_desc in test_dates:
            jd_start = swe.julday(year, month, day, 0.0)

            for event_type, event_name in [
                (SE_CALC_RISE, "Rise"),
                (SE_CALC_SET, "Set"),
            ]:
                passed, diff, error = compare_rise_trans(
                    jd_start,
                    SE_SUN,
                    "Sun",
                    event_type,
                    event_name,
                    lat,
                    lon,
                    alt,
                    f"{loc_name} {date_desc}",
                    verbose,
                )
                stats.add_result(passed, diff, error)

    # 4. True horizon tests
    print_section("RISE/SET - TRUE HORIZON")
    jd_start = swe.julday(2024, 6, 15, 0.0)
    lat, lon, alt = 45.0, 0.0, 100  # 100m altitude

    horizon_altitudes = [0.0, -0.5, -1.0, 0.5]  # Various horizon depressions/elevations

    for hor_alt in horizon_altitudes:
        for event_type, event_name in [(SE_CALC_RISE, "Rise"), (SE_CALC_SET, "Set")]:
            passed, diff, error = compare_rise_trans_true_hor(
                jd_start,
                SE_SUN,
                "Sun",
                event_type,
                event_name,
                lat,
                lon,
                alt,
                "Mid-Lat",
                hor_alt,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # 5. Moon at various locations (faster motion = more challenging)
    print_section("RISE/SET - MOON COMPREHENSIVE")
    for loc_name, lat, lon, alt in RISE_SET_LOCATIONS[:4]:
        jd_start = swe.julday(2024, 11, 15, 0.0)

        for event_type, event_name in EVENT_TYPES:
            passed, diff, error = compare_rise_trans(
                jd_start,
                SE_MOON,
                "Moon",
                event_type,
                event_name,
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("RISE/SET/TRANSIT COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_rise_transit.py [OPTIONS]")
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
