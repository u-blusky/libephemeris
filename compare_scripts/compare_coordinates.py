"""
Coordinate Transformation Comparison Script.

Compares coordinate transformation functions between pyswisseph and libephemeris:
- cotrans / cotrans_sp - ecliptic/equatorial/horizontal coordinate transforms
- azalt / azalt_rev - azimuth/altitude calculations
- refrac / refrac_extended - atmospheric refraction
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
    angular_diff,
    format_coord,
    format_diff,
    STANDARD_SUBJECTS,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class CoordTolerance:
    """Tolerance thresholds for coordinate comparisons."""

    ANGLE_DEGREES = 0.0001  # 0.0001 degree (0.36 arcsec)
    REFRACTION = 0.001  # 0.001 degree for refraction


# ============================================================================
# TEST DATA
# ============================================================================

# Test coordinates for transformations
TEST_COORDINATES = [
    # (lon, lat, description)
    (0.0, 0.0, "Vernal Equinox"),
    (90.0, 0.0, "Summer Solstice"),
    (180.0, 0.0, "Autumnal Equinox"),
    (270.0, 0.0, "Winter Solstice"),
    (45.0, 23.5, "Mid ecliptic with lat"),
    (120.0, -15.0, "Negative latitude"),
    (315.0, 45.0, "High latitude"),
]

# Obliquity values for testing
OBLIQUITY_VALUES = [
    (23.4393, "Modern obliquity"),
    (23.5, "Round value"),
    (23.0, "Lower obliquity"),
    (24.0, "Higher obliquity"),
]

# Test locations for azalt
AZALT_LOCATIONS = [
    ("Equator", 0.0, 0.0, 0),
    ("Mid-lat North", 45.0, 0.0, 0),
    ("High-lat North", 70.0, 0.0, 0),
    ("Mid-lat South", -35.0, 0.0, 0),
]

# Altitude values for refraction tests
ALTITUDE_VALUES = [0.0, 5.0, 10.0, 20.0, 45.0, 70.0, 85.0, 90.0]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_cotrans(lon, lat, eps, desc, verbose=False):
    """Compare cotrans function (coordinate transformation)."""
    coords = (lon, lat, 1.0)  # 1.0 for distance (unused in transform)

    # SwissEphemeris
    try:
        result_swe = swe.cotrans(coords, -eps)  # Negative eps: ecl->equ
        lon_out_swe, lat_out_swe = result_swe[0], result_swe[1]
    except Exception as e:
        print(f"[cotrans] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.cotrans(coords, -eps)
        lon_out_py, lat_out_py = result_py[0], result_py[1]
    except Exception as e:
        print(f"[cotrans] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = angular_diff(lon_out_swe, lon_out_py)
    diff_lat = abs(lat_out_swe - lat_out_py)
    max_diff = max(diff_lon, diff_lat)
    passed = max_diff < CoordTolerance.ANGLE_DEGREES

    if verbose:
        print(f"\n[cotrans] {desc}")
        print(f"  Input: Lon={lon:.4f} Lat={lat:.4f} Eps={eps:.4f}")
        print(f"  SWE: Lon={lon_out_swe:.8f} Lat={lat_out_swe:.8f}")
        print(f"  PY:  Lon={lon_out_py:.8f} Lat={lat_out_py:.8f}")
        print(
            f"  Diff: Lon={diff_lon:.10f} Lat={diff_lat:.10f} {format_status(passed)}"
        )
    else:
        status = format_status(passed)
        print(
            f"[cotrans] {desc}: "
            f"In=({lon:.1f},{lat:.1f}) "
            f"Out_SWE=({lon_out_swe:.4f},{lat_out_swe:.4f}) "
            f"Out_PY=({lon_out_py:.4f},{lat_out_py:.4f}) "
            f"MaxDiff={max_diff:.8f} {status}"
        )

    return passed, max_diff, False


def compare_cotrans_sp(lon, lat, lon_speed, lat_speed, eps, desc, verbose=False):
    """Compare cotrans_sp function (with speed)."""
    coords = (lon, lat, 1.0, lon_speed, lat_speed, 0.0)

    # SwissEphemeris
    try:
        result_swe = swe.cotrans_sp(coords, -eps)
    except Exception as e:
        print(f"[cotrans_sp] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.cotrans_sp(coords, -eps)
    except Exception as e:
        print(f"[cotrans_sp] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_lon = angular_diff(result_swe[0], result_py[0])
    diff_lat = abs(result_swe[1] - result_py[1])
    diff_lon_speed = abs(result_swe[3] - result_py[3])
    diff_lat_speed = abs(result_swe[4] - result_py[4])

    max_diff = max(diff_lon, diff_lat)
    passed = max_diff < CoordTolerance.ANGLE_DEGREES

    if verbose:
        print(f"\n[cotrans_sp] {desc}")
        print(f"  Input: Lon={lon:.4f} Lat={lat:.4f} Eps={eps:.4f}")
        print(f"  SWE: Lon={result_swe[0]:.6f} Lat={result_swe[1]:.6f}")
        print(f"  PY:  Lon={result_py[0]:.6f} Lat={result_py[1]:.6f}")
        print(f"  SpeedDiff: Lon={diff_lon_speed:.8f} Lat={diff_lat_speed:.8f}")
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[cotrans_sp] {desc}: "
            f"PosDiff={max_diff:.8f} SpeedDiff={max(diff_lon_speed, diff_lat_speed):.8f} {status}"
        )

    return passed, max_diff, False


def compare_azalt(jd, lat, lon, alt, ra, dec, desc, verbose=False):
    """Compare azalt function (equatorial to horizontal)."""
    geopos = (lon, lat, alt)
    # xin is RA/Dec in degrees
    xin = (ra, dec, 1.0)

    # SwissEphemeris
    try:
        result_swe = swe.azalt(jd, pyephem.SE_EQU2HOR, geopos, 1013.25, 15.0, xin)
        az_swe, alt_swe = result_swe[0], result_swe[1]
    except Exception as e:
        print(f"[azalt] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.azalt(jd, pyephem.SE_EQU2HOR, geopos, 1013.25, 15.0, xin)
        az_py, alt_py = result_py[0], result_py[1]
    except Exception as e:
        print(f"[azalt] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_az = angular_diff(az_swe, az_py)
    diff_alt = abs(alt_swe - alt_py)
    max_diff = max(diff_az, diff_alt)
    passed = max_diff < CoordTolerance.ANGLE_DEGREES

    if verbose:
        print(f"\n[azalt] {desc}")
        print(f"  Input: RA={ra:.4f} Dec={dec:.4f} Lat={lat:.1f}")
        print(f"  SWE: Az={az_swe:.6f} Alt={alt_swe:.6f}")
        print(f"  PY:  Az={az_py:.6f} Alt={alt_py:.6f}")
        print(f"  Diff: Az={diff_az:.8f} Alt={diff_alt:.8f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[azalt] {desc}: "
            f"SWE=({az_swe:.2f},{alt_swe:.2f}) "
            f"PY=({az_py:.2f},{alt_py:.2f}) "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_azalt_rev(jd, lat, lon, alt, az, altitude, desc, verbose=False):
    """Compare azalt_rev function (horizontal to equatorial)."""
    geopos = (lon, lat, alt)
    xin = (az, altitude, 1.0)

    # SwissEphemeris
    try:
        result_swe = swe.azalt_rev(jd, pyephem.SE_HOR2EQU, geopos, xin)
        ra_swe, dec_swe = result_swe[0], result_swe[1]
    except Exception as e:
        print(f"[azalt_rev] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.azalt_rev(jd, pyephem.SE_HOR2EQU, geopos, xin)
        ra_py, dec_py = result_py[0], result_py[1]
    except Exception as e:
        print(f"[azalt_rev] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff_ra = angular_diff(ra_swe, ra_py)
    diff_dec = abs(dec_swe - dec_py)
    max_diff = max(diff_ra, diff_dec)
    passed = max_diff < CoordTolerance.ANGLE_DEGREES

    if verbose:
        print(f"\n[azalt_rev] {desc}")
        print(f"  Input: Az={az:.4f} Alt={altitude:.4f} Lat={lat:.1f}")
        print(f"  SWE: RA={ra_swe:.6f} Dec={dec_swe:.6f}")
        print(f"  PY:  RA={ra_py:.6f} Dec={dec_py:.6f}")
        print(f"  Diff: RA={diff_ra:.8f} Dec={diff_dec:.8f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[azalt_rev] {desc}: "
            f"SWE=({ra_swe:.2f},{dec_swe:.2f}) "
            f"PY=({ra_py:.2f},{dec_py:.2f}) "
            f"MaxDiff={max_diff:.6f} {status}"
        )

    return passed, max_diff, False


def compare_refrac(inalt, calc_type, desc, verbose=False):
    """Compare refrac function."""
    # SwissEphemeris
    try:
        result_swe = swe.refrac(inalt, 1013.25, 15.0, calc_type)
    except Exception as e:
        print(f"[refrac] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.refrac(inalt, 1013.25, 15.0, calc_type)
    except Exception as e:
        print(f"[refrac] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(result_swe - result_py)
    passed = diff < CoordTolerance.REFRACTION

    if verbose:
        print(f"\n[refrac] {desc}")
        print(f"  Input Alt: {inalt:.4f}")
        print(f"  SWE: {result_swe:.8f}")
        print(f"  PY:  {result_py:.8f}")
        print(f"  Diff: {diff:.10f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[refrac] {desc}: "
            f"In={inalt:.1f} SWE={result_swe:.6f} PY={result_py:.6f} "
            f"Diff={diff:.8f} {status}"
        )

    return passed, diff, False


def compare_refrac_extended(inalt, geoalt, calc_type, desc, verbose=False):
    """Compare refrac_extended function."""
    # SwissEphemeris
    try:
        result_swe = swe.refrac_extended(
            inalt, geoalt, 1013.25, 15.0, 0.0065, calc_type
        )
        # Returns tuple of values
        refr_swe = result_swe[0] if isinstance(result_swe, tuple) else result_swe
    except Exception as e:
        print(f"[refrac_extended] {desc}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        result_py = pyephem.refrac_extended(
            inalt, geoalt, 1013.25, 15.0, 0.0065, calc_type
        )
        refr_py = result_py[0] if isinstance(result_py, tuple) else result_py
    except Exception as e:
        print(f"[refrac_extended] {desc}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(refr_swe - refr_py)
    passed = diff < CoordTolerance.REFRACTION

    if verbose:
        print(f"\n[refrac_extended] {desc}")
        print(f"  Input Alt: {inalt:.4f}, Geo Alt: {geoalt}")
        print(f"  SWE: {refr_swe:.8f}")
        print(f"  PY:  {refr_py:.8f}")
        print(f"  Diff: {diff:.10f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[refrac_extended] {desc}: "
            f"In={inalt:.1f} SWE={refr_swe:.6f} PY={refr_py:.6f} "
            f"Diff={diff:.8f} {status}"
        )

    return passed, diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose=False):
    """Run all coordinate transformation comparisons."""
    print_header("COORDINATE TRANSFORMATIONS COMPARISON")
    stats = TestStatistics()

    # 1. cotrans - ecliptic to equatorial
    print_section("COORDINATE TRANSFORM (cotrans)")
    for lon, lat, coord_desc in TEST_COORDINATES:
        for eps, eps_desc in OBLIQUITY_VALUES[:2]:
            desc = f"{coord_desc} @ {eps_desc}"
            passed, diff, error = compare_cotrans(lon, lat, eps, desc, verbose)
            stats.add_result(passed, diff, error)

    # 2. cotrans_sp - with speed
    print_section("COORDINATE TRANSFORM WITH SPEED (cotrans_sp)")
    eps = 23.4393
    for lon, lat, coord_desc in TEST_COORDINATES[:4]:
        desc = f"{coord_desc} with speed"
        passed, diff, error = compare_cotrans_sp(lon, lat, 1.0, 0.1, eps, desc, verbose)
        stats.add_result(passed, diff, error)

    # 3. azalt - equatorial to horizontal
    print_section("EQUATORIAL TO HORIZONTAL (azalt)")
    jd = swe.julday(2024, 6, 15, 22.0)  # Evening
    test_objects = [
        (0.0, 0.0, "RA=0, Dec=0"),
        (90.0, 45.0, "RA=90, Dec=45"),
        (180.0, -23.4, "RA=180, Dec=-23.4"),
        (270.0, 60.0, "RA=270, Dec=60"),
    ]
    for loc_name, lat, lon, alt in AZALT_LOCATIONS[:2]:
        for ra, dec, obj_desc in test_objects:
            desc = f"{obj_desc} @ {loc_name}"
            passed, diff, error = compare_azalt(
                jd, lat, lon, alt, ra, dec, desc, verbose
            )
            stats.add_result(passed, diff, error)

    # 4. azalt_rev - horizontal to equatorial
    print_section("HORIZONTAL TO EQUATORIAL (azalt_rev)")
    az_alt_pairs = [
        (0.0, 45.0, "North, 45 alt"),
        (90.0, 30.0, "East, 30 alt"),
        (180.0, 60.0, "South, 60 alt"),
        (270.0, 15.0, "West, 15 alt"),
    ]
    for loc_name, lat, lon, alt in AZALT_LOCATIONS[:2]:
        for az, altitude, obj_desc in az_alt_pairs:
            desc = f"{obj_desc} @ {loc_name}"
            passed, diff, error = compare_azalt_rev(
                jd, lat, lon, alt, az, altitude, desc, verbose
            )
            stats.add_result(passed, diff, error)

    # 5. refrac - true to apparent
    print_section("ATMOSPHERIC REFRACTION (refrac)")
    for alt in ALTITUDE_VALUES:
        desc = f"Alt={alt:.0f} (true->app)"
        passed, diff, error = compare_refrac(alt, pyephem.SE_TRUE_TO_APP, desc, verbose)
        stats.add_result(passed, diff, error)

    # 6. refrac - apparent to true
    print_section("ATMOSPHERIC REFRACTION REVERSE (refrac)")
    for alt in ALTITUDE_VALUES:
        desc = f"Alt={alt:.0f} (app->true)"
        passed, diff, error = compare_refrac(alt, pyephem.SE_APP_TO_TRUE, desc, verbose)
        stats.add_result(passed, diff, error)

    # 7. refrac_extended
    print_section("EXTENDED REFRACTION (refrac_extended)")
    for alt in ALTITUDE_VALUES[1:5]:  # Skip 0 and very high
        for geoalt in [0, 1000, 3000]:
            desc = f"Alt={alt:.0f} @ GeoAlt={geoalt}m"
            passed, diff, error = compare_refrac_extended(
                alt, geoalt, pyephem.SE_TRUE_TO_APP, desc, verbose
            )
            stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("COORDINATE TRANSFORMATIONS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_coordinates.py [OPTIONS]")
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
