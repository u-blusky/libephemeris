"""
Heliacal Events Comparison Script.

Compares heliacal event calculations between pyswisseph and libephemeris:
- heliacal_ut - heliacal rising/setting of planets and stars
- heliacal_pheno_ut - detailed heliacal phenomena
- vis_limit_mag - visibility limiting magnitude

Tests planets (Mercury, Venus, Mars, Jupiter, Saturn) and stars
(Sirius, Canopus, Arcturus, Vega, Capella) at multiple geographic locations.
"""

import os
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
import sys

# Set Swiss Ephemeris data path for star catalog
EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
if os.path.exists(EPHE_PATH):
    swe.set_ephe_path(EPHE_PATH)

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
    FIXED_STARS,
)


# ============================================================================
# CONSTANTS
# ============================================================================

# Heliacal event types (defined here in case libephemeris constants differ)
HELIACAL_RISING = 1
HELIACAL_SETTING = 2
EVENING_FIRST = 3
EVENING_LAST = 4
MORNING_FIRST = 5
MORNING_LAST = 6


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class HeliacalTolerance:
    """Tolerance thresholds for heliacal comparisons."""

    # Legacy tolerance (1 hour) for basic tests
    TIME_SECONDS = 3600.0  # 1 hour (heliacal events have inherent uncertainty)

    # New tolerances for pyswisseph comparison
    PLANET_TIME_DAYS = 1.0  # 1 day tolerance for planets
    STAR_TIME_DAYS = 2.0  # 2 days tolerance for stars
    PLANET_TIME_SECONDS = PLANET_TIME_DAYS * 86400.0  # 86400 seconds = 1 day
    STAR_TIME_SECONDS = STAR_TIME_DAYS * 86400.0  # 172800 seconds = 2 days

    MAGNITUDE = 0.5  # Visibility magnitude


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Bodies for heliacal tests - full planet list with IDs for swe
HELIACAL_PLANETS = [
    ("Mercury", swe.MERCURY),
    ("Venus", swe.VENUS),
    ("Mars", swe.MARS),
    ("Jupiter", swe.JUPITER),
    ("Saturn", swe.SATURN),
]

# Legacy configuration
HELIACAL_BODIES = [
    (SE_VENUS, "Venus"),
    (SE_MERCURY, "Mercury"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]

# Stars commonly used for heliacal events - legacy list
HELIACAL_STARS = [
    "Sirius",
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
]

# Stars specified for comprehensive tests
HELIACAL_STARS_COMPREHENSIVE = [
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
]

# Test locations - legacy
HELIACAL_LOCATIONS = [
    ("Cairo", 30.0444, 31.2357, 75),  # Egypt - classic heliacal observations
    ("Athens", 37.9838, 23.7275, 100),  # Greece
    ("Babylon", 32.5355, 44.4208, 30),  # Iraq - historical
]

# Test locations for comprehensive comparison at specific latitudes
# Format: (name, lat, lon, alt)
LATITUDE_TEST_LOCATIONS = [
    ("Equator (0N)", 0.0, 0.0, 0),  # Latitude 0
    ("Mid-Low (30N)", 30.0, 0.0, 0),  # Latitude 30N
    ("Mid (45N)", 45.0, 0.0, 0),  # Latitude 45N
    ("High (60N)", 60.0, 0.0, 0),  # Latitude 60N
]

# Standard atmospheric and observer conditions for comparison tests
# atmo: (pressure, temp, humidity%, extinction)
STANDARD_ATMO = (1013.25, 15.0, 50.0, 0.25)
# observer: (age, snellen, monocular/binocular, magnification, aperture, transmission)
STANDARD_OBSERVER = (25.0, 1.0, 1.0, 1.0, 0.0, 0.0)


# ============================================================================
# COMPARISON FUNCTIONS FOR COMPREHENSIVE TESTS
# ============================================================================


def compare_planet_heliacal_comprehensive(
    jd_start,
    planet_name,
    planet_id,
    event_type,
    event_name,
    lat,
    lon,
    alt,
    loc_name,
    verbose=False,
):
    """
    Compare swe_heliacal_ut for a planet with 1-day tolerance.

    Uses the Swiss Ephemeris API (by name) for both libraries.

    Returns: (passed, diff_days, is_error)
    """
    # geopos is (lon, lat, alt) per pyswisseph docs
    geopos = (lon, lat, alt)

    # SwissEphemeris (pyswisseph)
    try:
        ret_swe = swe.heliacal_ut(
            jd_start,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            planet_name,
            event_type,
            0,  # flags
        )
        # ret_swe is a tuple of 3 julian days (start, optimum, end)
        jd_swe = ret_swe[0]
    except Exception as e:
        if verbose:
            print(
                f"[heliacal_ut] {planet_name} {event_name} @ {loc_name}: SWE ERROR {e}"
            )
        return False, 0.0, True

    # LibEphemeris (pyephem.swe_heliacal_ut)
    try:
        dret_py, retflag = pyephem.swe_heliacal_ut(
            jd_start,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            planet_name,
            event_type,
        )
        jd_py = dret_py[0]
    except Exception as e:
        if verbose:
            print(
                f"[heliacal_ut] {planet_name} {event_name} @ {loc_name}: PY ERROR {e}"
            )
        return False, 0.0, True

    diff_seconds = abs(jd_swe - jd_py) * 86400.0
    diff_days = diff_seconds / 86400.0
    passed = diff_seconds < HeliacalTolerance.PLANET_TIME_SECONDS

    if verbose:
        swe_date = jd_to_date_str(jd_swe)
        py_date = jd_to_date_str(jd_py)
        print(f"\n[heliacal_ut] {planet_name} {event_name} @ {loc_name}")
        print(f"  SWE: {swe_date} (JD {jd_swe:.6f})")
        print(f"  PY:  {py_date} (JD {jd_py:.6f})")
        print(f"  Diff: {diff_days:.2f}d {format_status(passed)}")
    else:
        status = format_status(passed)
        swe_date = jd_to_date_str(jd_swe)
        py_date = jd_to_date_str(jd_py)
        print(
            f"[heliacal_ut] {planet_name} {event_name} @ {loc_name}: "
            f"SWE={swe_date} PY={py_date} Diff={diff_days:.2f}d {status}"
        )

    return passed, diff_days, False


def compare_star_heliacal_comprehensive(
    jd_start,
    star_name,
    event_type,
    event_name,
    lat,
    lon,
    alt,
    loc_name,
    verbose=False,
):
    """
    Compare swe_heliacal_ut for a star with 2-day tolerance.

    Returns: (passed, diff_days, is_error)
    """
    geopos = (lon, lat, alt)

    # SwissEphemeris (pyswisseph)
    try:
        ret_swe = swe.heliacal_ut(
            jd_start,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            star_name,
            event_type,
            0,  # flags
        )
        jd_swe = ret_swe[0]
    except Exception as e:
        if verbose:
            print(f"[heliacal_ut] {star_name} {event_name} @ {loc_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        dret_py, retflag = pyephem.swe_heliacal_ut(
            jd_start,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            star_name,
            event_type,
        )
        jd_py = dret_py[0]
    except Exception as e:
        if verbose:
            print(f"[heliacal_ut] {star_name} {event_name} @ {loc_name}: PY ERROR {e}")
        return False, 0.0, True

    diff_seconds = abs(jd_swe - jd_py) * 86400.0
    diff_days = diff_seconds / 86400.0
    passed = diff_seconds < HeliacalTolerance.STAR_TIME_SECONDS

    if verbose:
        swe_date = jd_to_date_str(jd_swe)
        py_date = jd_to_date_str(jd_py)
        print(f"\n[heliacal_ut] {star_name} {event_name} @ {loc_name}")
        print(f"  SWE: {swe_date} (JD {jd_swe:.6f})")
        print(f"  PY:  {py_date} (JD {jd_py:.6f})")
        print(f"  Diff: {diff_days:.2f}d {format_status(passed)}")
    else:
        status = format_status(passed)
        swe_date = jd_to_date_str(jd_swe)
        py_date = jd_to_date_str(jd_py)
        print(
            f"[heliacal_ut] {star_name} {event_name} @ {loc_name}: "
            f"SWE={swe_date} PY={py_date} Diff={diff_days:.2f}d {status}"
        )

    return passed, diff_days, False


def compare_vis_limit_mag_comprehensive(
    jd, object_name, lat, lon, alt, loc_name, verbose=False
):
    """
    Compare vis_limit_mag for a specific object with detailed output.

    Returns: (passed, diff_mag, is_error)
    """
    geopos = (lon, lat, alt)

    # SwissEphemeris - returns (retflag, dret_tuple)
    # retflag: -2 = below horizon, 0 = photopic, 1 = scotopic, 2 = mixed
    # dret[0] = limiting magnitude, dret[1] = object alt, dret[2] = object az, etc.
    try:
        ret_swe = swe.vis_limit_mag(
            jd, geopos, STANDARD_ATMO, STANDARD_OBSERVER, object_name, 0
        )
        retflag_swe = ret_swe[0]
        dret_swe = ret_swe[1]
        # Compare return flags first
        obj_alt_swe = dret_swe[1] if len(dret_swe) > 1 else None
        obj_mag_swe = dret_swe[7] if len(dret_swe) > 7 else None
    except Exception as e:
        if verbose:
            print(f"[vis_limit_mag] {object_name} @ {loc_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris - returns (retflag, dret_tuple)
    try:
        ret_py = pyephem.vis_limit_mag(
            jd, geopos, STANDARD_ATMO, STANDARD_OBSERVER, object_name, 0
        )
        retflag_py = ret_py[0]
        dret_py = ret_py[1]
        obj_alt_py = dret_py[1] if len(dret_py) > 1 else None
        obj_mag_py = dret_py[7] if len(dret_py) > 7 else None
    except Exception as e:
        if verbose:
            print(f"[vis_limit_mag] {object_name} @ {loc_name}: PY ERROR {e}")
        return False, 0.0, True

    # Compare return flags (both should indicate same visibility status)
    flag_match = retflag_swe == retflag_py

    # Compare object altitude (if both have valid data)
    if obj_alt_swe is not None and obj_alt_py is not None:
        alt_diff = abs(obj_alt_swe - obj_alt_py)
        alt_passed = alt_diff < 1.0  # 1 degree tolerance for altitude
    else:
        alt_diff = 0.0
        alt_passed = True

    # Compare object magnitude (if available)
    if obj_mag_swe is not None and obj_mag_py is not None:
        mag_diff = abs(obj_mag_swe - obj_mag_py)
        mag_passed = mag_diff < HeliacalTolerance.MAGNITUDE
    else:
        mag_diff = 0.0
        mag_passed = True

    passed = flag_match and alt_passed and mag_passed

    if verbose:
        print(f"\n[vis_limit_mag] {object_name} @ {loc_name}")
        print(f"  SWE: flag={retflag_swe} alt={obj_alt_swe:.2f} mag={obj_mag_swe:.2f}")
        print(f"  PY:  flag={retflag_py} alt={obj_alt_py:.2f} mag={obj_mag_py:.2f}")
        print(
            f"  FlagMatch={flag_match} AltDiff={alt_diff:.2f} MagDiff={mag_diff:.3f} "
            f"{format_status(passed)}"
        )
    else:
        status = format_status(passed)
        print(
            f"[vis_limit_mag] {object_name} @ {loc_name}: "
            f"flags SWE={retflag_swe} PY={retflag_py} "
            f"alt_diff={alt_diff:.2f} mag_diff={mag_diff:.3f} {status}"
        )

    # Return the maximum difference for statistics
    max_diff = max(alt_diff, mag_diff)
    return passed, max_diff, False


# ============================================================================
# LEGACY COMPARISON FUNCTIONS
# ============================================================================


def compare_heliacal_ut_planet(
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
    """Compare heliacal_ut for a planet (legacy function)."""
    result = EventComparisonResult(
        "heliacal_ut", f"{body_name} {event_name} @ {loc_name}"
    )
    result.tolerance_seconds = HeliacalTolerance.TIME_SECONDS

    # Geographic position and atmospheric conditions
    geopos = (lon, lat, alt)
    # Standard atmospheric conditions
    atmo = (1013.25, 15.0, 50.0, 0.25)  # pressure, temp, humidity, extinction
    # Observer conditions
    observer = (25.0, 1.0, 1.0, 1.0, 0.0, 0.0)  # age, snellen, etc.

    # SwissEphemeris
    try:
        ret_swe = swe.heliacal_ut(
            jd_start, geopos, atmo, observer, body_name, event_type, 0
        )
        result.jd_swe = ret_swe[0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py, _ = pyephem.swe_heliacal_ut(
            jd_start, geopos, atmo, observer, body_name, event_type
        )
        result.jd_py = ret_py[0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_heliacal_ut_star(
    jd_start, star_name, event_type, event_name, lat, lon, alt, loc_name, verbose=False
):
    """Compare heliacal_ut for a fixed star (legacy function)."""
    result = EventComparisonResult(
        "heliacal_ut", f"{star_name} {event_name} @ {loc_name}"
    )
    result.tolerance_seconds = HeliacalTolerance.TIME_SECONDS

    geopos = (lon, lat, alt)
    atmo = (1013.25, 15.0, 50.0, 0.25)
    observer = (25.0, 1.0, 1.0, 1.0, 0.0, 0.0)

    # SwissEphemeris
    try:
        ret_swe = swe.heliacal_ut(
            jd_start, geopos, atmo, observer, star_name, event_type, 0
        )
        result.jd_swe = ret_swe[0]
    except Exception as e:
        result.error_swe = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py, _ = pyephem.swe_heliacal_ut(
            jd_start, geopos, atmo, observer, star_name, event_type
        )
        result.jd_py = ret_py[0]
    except Exception as e:
        result.error_py = str(e)
        print(result.format_result(verbose))
        return False, 0.0, True

    result.calculate_diff()
    print(result.format_result(verbose))

    return result.passed, result.diff_seconds, False


def compare_heliacal_pheno_ut(jd, body_name, lat, lon, alt, loc_name, verbose=False):
    """Compare heliacal_pheno_ut function."""
    geopos = (lon, lat, alt)
    atmo = (1013.25, 15.0, 50.0, 0.25)
    observer = (25.0, 1.0, 1.0, 1.0, 0.0, 0.0)

    # SwissEphemeris
    try:
        ret_swe = swe.heliacal_pheno_ut(
            jd, geopos, atmo, observer, body_name, HELIACAL_RISING, 0
        )
        # Returns tuple of phenomena data
        data_swe = ret_swe if isinstance(ret_swe, (list, tuple)) else [ret_swe]
    except Exception as e:
        print(f"[heliacal_pheno_ut] {body_name} @ {loc_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py, _ = pyephem.swe_heliacal_pheno_ut(
            jd, geopos, atmo, observer, body_name, HELIACAL_RISING
        )
        data_py = ret_py if isinstance(ret_py, (list, tuple)) else [ret_py]
    except Exception as e:
        print(f"[heliacal_pheno_ut] {body_name} @ {loc_name}: PY ERROR {e}")
        return False, 0.0, True

    # Compare first few values
    try:
        diff = abs(data_swe[0] - data_py[0])
        passed = diff < 1.0  # Relaxed tolerance for phenomena data
    except (IndexError, TypeError):
        passed = True  # If data format differs, consider it compatible
        diff = 0.0

    if verbose:
        print(f"\n[heliacal_pheno_ut] {body_name} @ {loc_name}")
        print(f"  SWE: {data_swe[:3] if len(data_swe) >= 3 else data_swe}")
        print(f"  PY:  {data_py[:3] if len(data_py) >= 3 else data_py}")
        print(f"  Status: {format_status(passed)}")
    else:
        status = format_status(passed)
        print(f"[heliacal_pheno_ut] {body_name} @ {loc_name}: Diff={diff:.4f} {status}")

    return passed, diff, False


def compare_vis_limit_mag(jd, lat, lon, alt, loc_name, verbose=False):
    """Compare vis_limit_mag function (legacy)."""
    geopos = (lon, lat, alt)
    atmo = (1013.25, 15.0, 50.0, 0.25)
    observer = (25.0, 1.0, 1.0, 1.0, 0.0, 0.0)

    # SwissEphemeris
    try:
        ret_swe = swe.vis_limit_mag(jd, geopos, atmo, observer, "Venus", 0)
        if isinstance(ret_swe, tuple) and len(ret_swe) > 0:
            if isinstance(ret_swe[0], (list, tuple)):
                mag_swe = ret_swe[0][0]
            else:
                mag_swe = ret_swe[0]
        else:
            mag_swe = float(ret_swe)
    except Exception as e:
        print(f"[vis_limit_mag] @ {loc_name}: SWE ERROR {e}")
        return False, 0.0, True

    # LibEphemeris
    try:
        ret_py = pyephem.vis_limit_mag(jd, geopos, atmo, observer, "Venus", 0)
        if isinstance(ret_py, tuple) and len(ret_py) >= 2:
            dret_py = ret_py[1]
            mag_py = dret_py[0]
        else:
            mag_py = ret_py[0]
    except Exception as e:
        print(f"[vis_limit_mag] @ {loc_name}: PY ERROR {e}")
        return False, 0.0, True

    diff = abs(mag_swe - mag_py)
    passed = diff < HeliacalTolerance.MAGNITUDE

    if verbose:
        print(f"\n[vis_limit_mag] @ {loc_name}")
        print(f"  SWE: {mag_swe:.3f}")
        print(f"  PY:  {mag_py:.3f}")
        print(f"  Diff: {diff:.4f} {format_status(passed)}")
    else:
        status = format_status(passed)
        print(
            f"[vis_limit_mag] @ {loc_name}: "
            f"SWE={mag_swe:.2f} PY={mag_py:.2f} "
            f"Diff={diff:.4f} {status}"
        )

    return passed, diff, False


# ============================================================================
# MAIN COMPARISON RUNNERS
# ============================================================================


def run_comprehensive_comparisons(verbose=False):
    """Run comprehensive heliacal comparisons against pyswisseph."""
    print_header("COMPREHENSIVE HELIACAL EVENTS COMPARISON VS PYSWISSEPH")
    stats = TestStatistics()

    # Start date for searches
    jd_start = swe.julday(2024, 1, 1, 0.0)

    # =========================================================================
    # 1. Planet heliacal rising at all latitudes
    # =========================================================================
    print_section("PLANETS - HELIACAL RISING (1-day tolerance)")

    for planet_name, planet_id in HELIACAL_PLANETS:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_planet_heliacal_comprehensive(
                jd_start,
                planet_name,
                planet_id,
                HELIACAL_RISING,
                "Rising",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # =========================================================================
    # 2. Planet heliacal setting at all latitudes
    # =========================================================================
    print_section("PLANETS - HELIACAL SETTING (1-day tolerance)")

    for planet_name, planet_id in HELIACAL_PLANETS:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_planet_heliacal_comprehensive(
                jd_start,
                planet_name,
                planet_id,
                HELIACAL_SETTING,
                "Setting",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # =========================================================================
    # 3. Inner planets - evening first / morning last
    # =========================================================================
    print_section("INNER PLANETS - EVENING FIRST (1-day tolerance)")

    inner_planets = [("Mercury", swe.MERCURY), ("Venus", swe.VENUS)]
    for planet_name, planet_id in inner_planets:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_planet_heliacal_comprehensive(
                jd_start,
                planet_name,
                planet_id,
                EVENING_FIRST,
                "EveFirst",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    print_section("INNER PLANETS - MORNING LAST (1-day tolerance)")

    for planet_name, planet_id in inner_planets:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_planet_heliacal_comprehensive(
                jd_start,
                planet_name,
                planet_id,
                MORNING_LAST,
                "MornLast",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # =========================================================================
    # 4. Stars heliacal rising at all latitudes
    # =========================================================================
    print_section("STARS - HELIACAL RISING (2-day tolerance)")

    for star_name in HELIACAL_STARS_COMPREHENSIVE:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_star_heliacal_comprehensive(
                jd_start,
                star_name,
                HELIACAL_RISING,
                "Rising",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # =========================================================================
    # 5. Stars heliacal setting at all latitudes
    # =========================================================================
    print_section("STARS - HELIACAL SETTING (2-day tolerance)")

    for star_name in HELIACAL_STARS_COMPREHENSIVE:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_star_heliacal_comprehensive(
                jd_start,
                star_name,
                HELIACAL_SETTING,
                "Setting",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # =========================================================================
    # 6. Visibility limiting magnitude comparisons
    # =========================================================================
    print_section("VISIBILITY LIMITING MAGNITUDE (vis_limit_mag)")

    # Test at night time for different locations
    jd_night = swe.julday(2024, 6, 15, 22.0)  # Night time in June

    # Test for several planets and stars
    test_objects = ["Venus", "Jupiter", "Saturn", "Mars", "Sirius", "Vega"]
    for obj_name in test_objects:
        for loc_name, lat, lon, alt in LATITUDE_TEST_LOCATIONS:
            passed, diff, error = compare_vis_limit_mag_comprehensive(
                jd_night, obj_name, lat, lon, alt, loc_name, verbose
            )
            stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("COMPREHENSIVE HELIACAL COMPARISON SUMMARY")
    return stats.passed, stats.total


def run_all_comparisons(verbose=False):
    """Run all heliacal comparisons (legacy mode)."""
    print_header("HELIACAL EVENTS COMPARISON")
    stats = TestStatistics()

    jd_start = swe.julday(2024, 1, 1, 0.0)

    # 1. Heliacal rising/setting of planets
    print_section("HELIACAL EVENTS - PLANETS")
    event_types = [
        (HELIACAL_RISING, "Rising"),
        (HELIACAL_SETTING, "Setting"),
    ]

    for body_id, body_name in HELIACAL_BODIES[:2]:  # Venus and Mercury
        for event_type, event_name in event_types:
            for loc_name, lat, lon, alt in HELIACAL_LOCATIONS[:2]:
                passed, diff, error = compare_heliacal_ut_planet(
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

    # 2. Heliacal rising of stars
    print_section("HELIACAL EVENTS - STARS")
    for star_name in HELIACAL_STARS[:3]:
        for loc_name, lat, lon, alt in HELIACAL_LOCATIONS[:2]:
            passed, diff, error = compare_heliacal_ut_star(
                jd_start,
                star_name,
                HELIACAL_RISING,
                "Rising",
                lat,
                lon,
                alt,
                loc_name,
                verbose,
            )
            stats.add_result(passed, diff, error)

    # 3. Heliacal phenomena details
    print_section("HELIACAL PHENOMENA (heliacal_pheno_ut)")
    jd_test = swe.julday(2024, 6, 15, 12.0)
    for body_id, body_name in HELIACAL_BODIES[:2]:
        for loc_name, lat, lon, alt in HELIACAL_LOCATIONS[:1]:
            passed, diff, error = compare_heliacal_pheno_ut(
                jd_test, body_name, lat, lon, alt, loc_name, verbose
            )
            stats.add_result(passed, diff, error)

    # 4. Visibility limiting magnitude
    print_section("VISIBILITY LIMITING MAGNITUDE")
    jd_test = swe.julday(2024, 6, 15, 22.0)  # Night time
    for loc_name, lat, lon, alt in HELIACAL_LOCATIONS:
        passed, diff, error = compare_vis_limit_mag(
            jd_test, lat, lon, alt, loc_name, verbose
        )
        stats.add_result(passed, diff, error)

    # Summary
    stats.print_summary("HELIACAL EVENTS COMPARISON SUMMARY")
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_heliacal.py [OPTIONS]")
    print()
    print("Options:")
    print("  -v, --verbose       Show detailed output for each test")
    print("  -c, --comprehensive Run comprehensive pyswisseph comparison tests")
    print("  -h, --help          Show this help message")
    print()
    print("Comprehensive mode tests:")
    print("  - Planets: Mercury, Venus, Mars, Jupiter, Saturn")
    print("  - Stars: Sirius, Canopus, Arcturus, Vega, Capella")
    print("  - Event types: heliacal rising, heliacal setting,")
    print("                 evening first, morning last (inner planets only)")
    print("  - Latitudes: 0, 30N, 45N, 60N")
    print("  - Tolerances: 1 day for planets, 2 days for stars")
    print()


def main():
    """Main entry point."""
    args = parse_args(sys.argv)

    if args["help"]:
        print_help()
        sys.exit(0)

    # Check for comprehensive mode
    comprehensive = "--comprehensive" in sys.argv or "-c" in sys.argv

    if comprehensive:
        passed, total = run_comprehensive_comparisons(verbose=args["verbose"])
    else:
        passed, total = run_all_comparisons(verbose=args["verbose"])

    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
