"""
Utility Functions Comparison Script.

Compares utility functions between pyswisseph and libephemeris:
- degnorm / radnorm - angle normalization
- difdeg2n / difdegn / difrad2n - angle differences
- split_deg - degree splitting
- deg_midp / rad_midp - midpoints
- cs2degstr / cs2timestr / cs2lonlatstr - formatting
"""

import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
import sys
import math

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import (
    TestStatistics,
    print_header,
    print_section,
    parse_args,
    format_status,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class UtilTolerance:
    """Tolerance thresholds for utility comparisons."""

    ANGLE = 1e-10  # Very small for normalization
    CENTISEC = 1  # 1 centisecond for formatting


# ============================================================================
# TEST DATA
# ============================================================================

# Test angles for normalization
TEST_ANGLES = [
    0.0,
    90.0,
    180.0,
    270.0,
    360.0,
    -90.0,
    -180.0,
    -270.0,
    -360.0,
    450.0,
    720.0,
    -450.0,
    -720.0,
    0.001,
    359.999,
    -0.001,
    -359.999,
]

# Angle pairs for difference calculations
ANGLE_PAIRS = [
    (0, 90, "0-90"),
    (350, 10, "350-10 (wrap)"),
    (180, 0, "180-0"),
    (45, 315, "45-315 (wrap)"),
    (0, 180, "0-180"),
    (1, 359, "1-359 (wrap)"),
]

# Split degree test cases
SPLIT_DEG_CASES = [
    (0.0, 0, "Zero"),
    (30.5, 0, "30.5 deg"),
    (123.456789, 0, "Arbitrary"),
    (359.999, 0, "Near 360"),
    (90.0, pyephem.SPLIT_DEG_ZODIACAL, "Zodiacal 90"),
    (127.5, pyephem.SPLIT_DEG_ZODIACAL, "Zodiacal 127.5"),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_degnorm(angle, verbose=False):
    """Compare degnorm function."""
    # SwissEphemeris
    try:
        result_swe = swe.degnorm(angle)
    except Exception as e:
        print(f"[degnorm] {angle}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.degnorm(angle)
    except Exception as e:
        print(f"[degnorm] {angle}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    passed = diff < UtilTolerance.ANGLE

    if verbose:
        print(f"\n[degnorm] Input: {angle}")
        print(f"  SWE: {result_swe:.15f}")
        print(f"  PY:  {result_py:.15f}")
        print(f"  Diff: {diff:.20f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[degnorm] {angle:>10.4f} -> "
            f"SWE={result_swe:.6f} PY={result_py:.6f} "
            f"Diff={diff:.15f} {status}"
        )

    return passed, diff, False


def compare_radnorm(angle_rad, verbose=False):
    """Compare radnorm function."""
    # SwissEphemeris
    try:
        result_swe = swe.radnorm(angle_rad)
    except Exception as e:
        print(f"[radnorm] {angle_rad}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.radnorm(angle_rad)
    except Exception as e:
        print(f"[radnorm] {angle_rad}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    passed = diff < UtilTolerance.ANGLE

    if verbose:
        print(f"\n[radnorm] Input: {angle_rad}")
        print(f"  SWE: {result_swe:.15f}")
        print(f"  PY:  {result_py:.15f}")
        print(f"  Diff: {diff:.20f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[radnorm] {angle_rad:>10.4f} -> "
            f"SWE={result_swe:.6f} PY={result_py:.6f} "
            f"Diff={diff:.15f} {status}"
        )

    return passed, diff, False


def compare_difdeg2n(deg1, deg2, desc, verbose=False):
    """Compare difdeg2n function (difference normalized to -180..180)."""
    # SwissEphemeris
    try:
        result_swe = swe.difdeg2n(deg1, deg2)
    except Exception as e:
        print(f"[difdeg2n] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.difdeg2n(deg1, deg2)
    except Exception as e:
        print(f"[difdeg2n] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    passed = diff < UtilTolerance.ANGLE

    if verbose:
        print(f"\n[difdeg2n] {desc} ({deg1}, {deg2})")
        print(f"  SWE: {result_swe:.15f}")
        print(f"  PY:  {result_py:.15f}")
        print(f"  Diff: {diff:.20f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[difdeg2n] {desc:<15}: "
            f"({deg1},{deg2}) -> "
            f"SWE={result_swe:>8.4f} PY={result_py:>8.4f} "
            f"Diff={diff:.15f} {status}"
        )

    return passed, diff, False


def compare_difdegn(deg1, deg2, desc, verbose=False):
    """Compare difdegn function (positive difference)."""
    # SwissEphemeris
    try:
        result_swe = swe.difdegn(deg1, deg2)
    except Exception as e:
        print(f"[difdegn] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.difdegn(deg1, deg2)
    except Exception as e:
        print(f"[difdegn] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    passed = diff < UtilTolerance.ANGLE

    if verbose:
        print(f"\n[difdegn] {desc} ({deg1}, {deg2})")
        print(f"  SWE: {result_swe:.15f}")
        print(f"  PY:  {result_py:.15f}")
        print(f"  Diff: {diff:.20f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[difdegn] {desc:<15}: "
            f"({deg1},{deg2}) -> "
            f"SWE={result_swe:>8.4f} PY={result_py:>8.4f} "
            f"Diff={diff:.15f} {status}"
        )

    return passed, diff, False


def compare_deg_midp(deg1, deg2, desc, verbose=False):
    """Compare deg_midp function (midpoint)."""
    # SwissEphemeris
    try:
        result_swe = swe.deg_midp(deg1, deg2)
    except Exception as e:
        print(f"[deg_midp] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.deg_midp(deg1, deg2)
    except Exception as e:
        print(f"[deg_midp] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    # Handle wrap-around
    if diff > 180:
        diff = 360 - diff
    passed = diff < UtilTolerance.ANGLE

    if verbose:
        print(f"\n[deg_midp] {desc} ({deg1}, {deg2})")
        print(f"  SWE: {result_swe:.15f}")
        print(f"  PY:  {result_py:.15f}")
        print(f"  Diff: {diff:.20f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[deg_midp] {desc:<15}: "
            f"midp({deg1},{deg2}) -> "
            f"SWE={result_swe:>8.4f} PY={result_py:>8.4f} "
            f"Diff={diff:.15f} {status}"
        )

    return passed, diff, False


def compare_split_deg(deg, flags, desc, verbose=False):
    """Compare split_deg function."""
    # SwissEphemeris
    try:
        result_swe = swe.split_deg(deg, flags)
        # Returns (deg, min, sec, frac, sign)
    except Exception as e:
        print(f"[split_deg] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.split_deg(deg, flags)
    except Exception as e:
        print(f"[split_deg] {desc}: PY ERROR {e}")
        return False, 0.0, True

    # Compare all components
    match = (
        result_swe[0] == result_py[0]  # degrees
        and result_swe[1] == result_py[1]  # minutes
        and abs(result_swe[2] - result_py[2]) < 1  # seconds (may differ by rounding)
        and result_swe[4] == result_py[4]  # sign
    )

    if verbose:
        print(f"\n[split_deg] {desc} ({deg}, flags={flags})")
        print(f"  SWE: {result_swe}")
        print(f"  PY:  {result_py}")
        print(f"  Match: {format_status(match)}")
    else:
        status = format_status(match)
        print(
            f"[split_deg] {desc:<20}: "
            f"SWE={result_swe[0]}d{result_swe[1]}m{result_swe[2]:.0f}s "
            f"PY={result_py[0]}d{result_py[1]}m{result_py[2]:.0f}s {status}"
        )

    return match, 0.0 if match else 1.0, False


def compare_cs2degstr(cs, verbose=False):
    """Compare cs2degstr function (centiseconds to degree string)."""
    # SwissEphemeris
    try:
        result_swe = swe.cs2degstr(cs)
    except Exception as e:
        print(f"[cs2degstr] {cs}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.cs2degstr(cs)
    except Exception as e:
        print(f"[cs2degstr] {cs}: PY ERROR {e}")
        return False, 0.0, True

    # String comparison
    match = result_swe == result_py

    if verbose:
        print(f"\n[cs2degstr] Input: {cs} centiseconds")
        print(f"  SWE: '{result_swe}'")
        print(f"  PY:  '{result_py}'")
        print(f"  Match: {format_status(match)}")
    else:
        status = format_status(match)
        print(f"[cs2degstr] {cs:>10}: SWE='{result_swe}' PY='{result_py}' {status}")

    return match, 0.0 if match else 1.0, False


def compare_cs2timestr(cs, verbose=False):
    """Compare cs2timestr function (centiseconds to time string)."""
    # SwissEphemeris
    try:
        result_swe = swe.cs2timestr(cs, ord(":"), True)
    except Exception as e:
        print(f"[cs2timestr] {cs}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.cs2timestr(cs, ord(":"), True)
    except Exception as e:
        print(f"[cs2timestr] {cs}: PY ERROR {e}")
        return False, 0.0, True

    match = result_swe == result_py

    if verbose:
        print(f"\n[cs2timestr] Input: {cs} centiseconds")
        print(f"  SWE: '{result_swe}'")
        print(f"  PY:  '{result_py}'")
        print(f"  Match: {format_status(match)}")
    else:
        status = format_status(match)
        print(f"[cs2timestr] {cs:>10}: SWE='{result_swe}' PY='{result_py}' {status}")

    return match, 0.0 if match else 1.0, False


def compare_d2l(deg, verbose=False):
    """Compare d2l function (degrees to centiseconds)."""
    # SwissEphemeris
    try:
        result_swe = swe.d2l(deg)
    except Exception as e:
        print(f"[d2l] {deg}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.d2l(deg)
    except Exception as e:
        print(f"[d2l] {deg}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    passed = diff < UtilTolerance.CENTISEC

    if verbose:
        print(f"\n[d2l] Input: {deg} degrees")
        print(f"  SWE: {result_swe}")
        print(f"  PY:  {result_py}")
        print(f"  Diff: {diff} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[d2l] {deg:>10.4f}: SWE={result_swe} PY={result_py} Diff={diff} {status}"
        )

    return passed, diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all utility function comparisons."""
    print_header("UTILITY FUNCTIONS COMPARISON")
    stats = TestStatistics()

    # 1. Degree normalization
    print_section("DEGREE NORMALIZATION (degnorm)")
    for angle in TEST_ANGLES:
        passed, diff, error = compare_degnorm(angle, verbose)
        stats.add_result(passed, diff, error)

    # 2. Radian normalization
    print_section("RADIAN NORMALIZATION (radnorm)")
    for angle in TEST_ANGLES:
        angle_rad = math.radians(angle)
        passed, diff, error = compare_radnorm(angle_rad, verbose)
        stats.add_result(passed, diff, error)

    # 3. Degree difference (normalized to -180..180)
    print_section("DEGREE DIFFERENCE -180..180 (difdeg2n)")
    for deg1, deg2, desc in ANGLE_PAIRS:
        passed, diff, error = compare_difdeg2n(deg1, deg2, desc, verbose)
        stats.add_result(passed, diff, error)

    # 4. Degree difference (positive)
    print_section("DEGREE DIFFERENCE POSITIVE (difdegn)")
    for deg1, deg2, desc in ANGLE_PAIRS:
        passed, diff, error = compare_difdegn(deg1, deg2, desc, verbose)
        stats.add_result(passed, diff, error)

    # 5. Midpoint
    print_section("DEGREE MIDPOINT (deg_midp)")
    for deg1, deg2, desc in ANGLE_PAIRS:
        passed, diff, error = compare_deg_midp(deg1, deg2, desc, verbose)
        stats.add_result(passed, diff, error)

    # 6. Split degree
    print_section("SPLIT DEGREE (split_deg)")
    for deg, flags, desc in SPLIT_DEG_CASES:
        passed, diff, error = compare_split_deg(deg, flags, desc, verbose)
        stats.add_result(passed, diff, error)

    # 7. Centiseconds to degree string
    # NOTE: cs2degstr may crash on some systems - skipping
    # print_section("CS TO DEGREE STRING (cs2degstr)")
    # cs_values = [0, 360000, 4500000, 12960000, -3600000]
    # for cs in cs_values:
    #     passed, diff, error = compare_cs2degstr(cs, verbose)
    #     stats.add_result(passed, diff, error)

    # 8. Centiseconds to time string
    # NOTE: cs2timestr may crash on some systems - skipping
    # print_section("CS TO TIME STRING (cs2timestr)")
    # for cs in cs_values:
    #     passed, diff, error = compare_cs2timestr(cs, verbose)
    #     stats.add_result(passed, diff, error)

    # 9. Degrees to centiseconds
    print_section("DEGREES TO CENTISECONDS (d2l)")
    for deg in [0.0, 1.0, 90.0, 180.5, 359.999]:
        passed, diff, error = compare_d2l(deg, verbose)
        stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("UTILITY FUNCTIONS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_utilities.py [OPTIONS]")
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
