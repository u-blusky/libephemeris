"""
Planetary Phenomena Comparison Script.

Compares planetary phenomena calculations between pyswisseph and libephemeris:
- pheno / pheno_ut - phase angle, elongation, magnitude, etc.
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
    format_coord,
    format_diff,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class PhenoTolerance:
    """Tolerance thresholds for phenomena comparisons."""

    ANGLE_DEGREES = 0.01  # Phase angle, elongation
    MAGNITUDE = 0.1  # Apparent magnitude
    DIAMETER = 0.01  # Apparent diameter in arcsec


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for phenomena tests
# Note: pheno works for planets visible from Earth
PHENO_PLANETS = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Phenomena array indices
PHENO_PHASE_ANGLE = 0
PHENO_PHASE = 1  # Illuminated fraction
PHENO_ELONGATION = 2
PHENO_DIAMETER = 3
PHENO_MAGNITUDE = 4


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_pheno_ut(jd, body_id, body_name, verbose=False):
    """Compare pheno_ut function."""
    # SwissEphemeris
    try:
        ret_swe = swe.pheno_ut(jd, body_id, SEFLG_SWIEPH)
        # Returns (retflags, attr) where attr contains phenomena data
        attr_swe = (
            ret_swe[1] if isinstance(ret_swe, tuple) and len(ret_swe) > 1 else ret_swe
        )

        phase_angle_swe = attr_swe[PHENO_PHASE_ANGLE]
        phase_swe = attr_swe[PHENO_PHASE]
        elongation_swe = attr_swe[PHENO_ELONGATION]
        diameter_swe = attr_swe[PHENO_DIAMETER]
        magnitude_swe = attr_swe[PHENO_MAGNITUDE]
    except Exception as e:
        print(f"[pheno_ut] {body_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.pheno_ut(jd, body_id, SEFLG_SWIEPH)
        attr_py = ret_py[1] if isinstance(ret_py, tuple) and len(ret_py) > 1 else ret_py

        phase_angle_py = attr_py[PHENO_PHASE_ANGLE]
        phase_py = attr_py[PHENO_PHASE]
        elongation_py = attr_py[PHENO_ELONGATION]
        diameter_py = attr_py[PHENO_DIAMETER]
        magnitude_py = attr_py[PHENO_MAGNITUDE]
    except Exception as e:
        print(f"[pheno_ut] {body_name}: PY ERROR {e}")
        return False, 0.0, True

    # Calculate differences
    diff_phase_angle = abs(phase_angle_swe - phase_angle_py)
    diff_phase = abs(phase_swe - phase_py)
    diff_elongation = abs(elongation_swe - elongation_py)
    diff_diameter = abs(diameter_swe - diameter_py)
    diff_magnitude = abs(magnitude_swe - magnitude_py)

    # Check tolerances
    passed_phase = diff_phase_angle < PhenoTolerance.ANGLE_DEGREES
    passed_elong = diff_elongation < PhenoTolerance.ANGLE_DEGREES
    passed_mag = diff_magnitude < PhenoTolerance.MAGNITUDE

    passed = passed_phase and passed_elong and passed_mag
    max_diff = max(diff_phase_angle, diff_elongation, diff_magnitude)

    if verbose:
        print(f"\n[pheno_ut] {body_name}")
        print(
            f"  Phase Angle: SWE={phase_angle_swe:.6f} PY={phase_angle_py:.6f} Diff={diff_phase_angle:.8f} {format_status(passed_phase)}"
        )
        print(
            f"  Phase:       SWE={phase_swe:.6f} PY={phase_py:.6f} Diff={diff_phase:.8f}"
        )
        print(
            f"  Elongation:  SWE={elongation_swe:.6f} PY={elongation_py:.6f} Diff={diff_elongation:.8f} {format_status(passed_elong)}"
        )
        print(
            f"  Diameter:    SWE={diameter_swe:.6f} PY={diameter_py:.6f} Diff={diff_diameter:.8f} arcsec"
        )
        print(
            f"  Magnitude:   SWE={magnitude_swe:.6f} PY={magnitude_py:.6f} Diff={diff_magnitude:.8f} {format_status(passed_mag)}"
        )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[pheno_ut] {body_name:<10}: "
            f"PhaseAng={phase_angle_swe:.2f}/{phase_angle_py:.2f} "
            f"Elong={elongation_swe:.2f}/{elongation_py:.2f} "
            f"Mag={magnitude_swe:.2f}/{magnitude_py:.2f} "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_pheno(jd_et, body_id, body_name, verbose=False):
    """Compare pheno function (ET version)."""
    # SwissEphemeris
    try:
        ret_swe = swe.pheno(jd_et, body_id, SEFLG_SWIEPH)
        attr_swe = (
            ret_swe[1] if isinstance(ret_swe, tuple) and len(ret_swe) > 1 else ret_swe
        )

        phase_angle_swe = attr_swe[PHENO_PHASE_ANGLE]
        elongation_swe = attr_swe[PHENO_ELONGATION]
        magnitude_swe = attr_swe[PHENO_MAGNITUDE]
    except Exception as e:
        print(f"[pheno] {body_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.pheno(jd_et, body_id, SEFLG_SWIEPH)
        attr_py = ret_py[1] if isinstance(ret_py, tuple) and len(ret_py) > 1 else ret_py

        phase_angle_py = attr_py[PHENO_PHASE_ANGLE]
        elongation_py = attr_py[PHENO_ELONGATION]
        magnitude_py = attr_py[PHENO_MAGNITUDE]
    except Exception as e:
        print(f"[pheno] {body_name}: PY ERROR {e}")
        return False, 0.0, True

    diff_phase_angle = abs(phase_angle_swe - phase_angle_py)
    diff_elongation = abs(elongation_swe - elongation_py)
    diff_magnitude = abs(magnitude_swe - magnitude_py)

    passed = (
        diff_phase_angle < PhenoTolerance.ANGLE_DEGREES
        and diff_elongation < PhenoTolerance.ANGLE_DEGREES
        and diff_magnitude < PhenoTolerance.MAGNITUDE
    )
    max_diff = max(diff_phase_angle, diff_elongation, diff_magnitude)

    if verbose:
        print(f"\n[pheno] {body_name}")
        print(
            f"  Phase Angle: SWE={phase_angle_swe:.6f} PY={phase_angle_py:.6f} Diff={diff_phase_angle:.8f}"
        )
        print(
            f"  Elongation:  SWE={elongation_swe:.6f} PY={elongation_py:.6f} Diff={diff_elongation:.8f}"
        )
        print(
            f"  Magnitude:   SWE={magnitude_swe:.6f} PY={magnitude_py:.6f} Diff={diff_magnitude:.8f}"
        )
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[pheno] {body_name:<10}: "
            f"PhaseAng={phase_angle_swe:.2f}/{phase_angle_py:.2f} "
            f"Elong={elongation_swe:.2f}/{elongation_py:.2f} "
            f"Mag={magnitude_swe:.2f}/{magnitude_py:.2f} {status}"
        )

    return passed, max_diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all phenomena comparisons."""
    print_header("PLANETARY PHENOMENA COMPARISON")
    stats = TestStatistics()

    # Test dates covering different planetary configurations
    test_dates = [
        (2024, 1, 15, "Jan 2024"),
        (2024, 4, 15, "Apr 2024"),
        (2024, 7, 15, "Jul 2024"),
        (2024, 10, 15, "Oct 2024"),
    ]

    # 1. pheno_ut for all planets at various dates
    print_section("PLANETARY PHENOMENA (pheno_ut)")
    for year, month, day, date_desc in test_dates:
        jd = swe.julday(year, month, day, 12.0)

        for body_id, body_name in PHENO_PLANETS:
            passed, diff, error = compare_pheno_ut(jd, body_id, body_name, verbose)
            stats.add_result(passed, diff, error)

    # 2. pheno (ET version) for subset
    print_section("PLANETARY PHENOMENA ET (pheno)")
    jd_et = swe.julday(2024, 6, 15, 12.0) + swe.deltat(swe.julday(2024, 6, 15, 12.0))
    for body_id, body_name in PHENO_PLANETS:
        passed, diff, error = compare_pheno(jd_et, body_id, body_name, verbose)
        stats.add_result(passed, diff, error)

    # 3. Special configurations - inner planets at various elongations
    print_section("INNER PLANETS AT VARIOUS ELONGATIONS")
    # Venus and Mercury have interesting phase variations
    inner_planet_dates = [
        (2024, 3, 22, "Venus western elongation"),
        (2024, 6, 4, "Mercury eastern elongation"),
        (2024, 8, 26, "Mercury western elongation"),
    ]

    for year, month, day, desc in inner_planet_dates:
        jd = swe.julday(year, month, day, 12.0)

        for body_id, body_name in [(SE_VENUS, "Venus"), (SE_MERCURY, "Mercury")]:
            passed, diff, error = compare_pheno_ut(jd, body_id, body_name, verbose)
            stats.add_result(passed, diff, error)

    # 4. Mars at opposition (brightest)
    print_section("OUTER PLANETS AT OPPOSITION")
    opposition_dates = [
        (2025, 1, 16, SE_MARS, "Mars opposition 2025"),
        (2024, 12, 7, SE_JUPITER, "Jupiter opposition 2024"),
        (2024, 9, 8, SE_SATURN, "Saturn opposition 2024"),
    ]

    for year, month, day, body_id, desc in opposition_dates:
        jd = swe.julday(year, month, day, 12.0)
        body_name = desc.split()[0]
        passed, diff, error = compare_pheno_ut(jd, body_id, body_name, verbose)
        stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("PLANETARY PHENOMENA COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_phenomena.py [OPTIONS]")
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
