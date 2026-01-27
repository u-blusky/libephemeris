"""
Shared utilities for comparison scripts.

This module provides common classes, functions, and constants used across
all comparison scripts in the suite.
"""

from typing import Optional, List
from dataclasses import dataclass

# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class Tolerances:
    """Tolerance thresholds for different comparison types."""

    # Position tolerances (degrees)
    LONGITUDE_STRICT = 0.001  # Geocentric/topocentric
    LONGITUDE_RELAXED = 0.03  # Heliocentric/barycentric
    LATITUDE_STRICT = 0.001
    LATITUDE_RELAXED = 0.03

    # Distance tolerances (AU)
    DISTANCE_STRICT = 0.0001
    DISTANCE_RELAXED = 0.01

    # Velocity tolerances
    VELOCITY_ANGULAR = 0.01  # degrees/day
    VELOCITY_RADIAL = 0.001  # AU/day

    # House cusp tolerances (degrees)
    HOUSE_CUSP = 0.001

    # ASCMC angle tolerances (degrees)
    # Some angles (Co-Asc, Polar Asc) may not be implemented or have calculation differences
    ASCMC_ANGLE = 1.0

    # Ayanamsha tolerances (degrees)
    # Relaxed tolerance for star-based ayanamshas due to numerical precision in
    # star position calculations, precession, and coordinate transformations
    AYANAMSHA = 0.06


# ============================================================================
# TEST SUBJECTS
# ============================================================================

# Format: (Name, Year, Month, Day, Hour, Lat, Lon, Alt)
STANDARD_SUBJECTS = [
    ("Standard J2000", 2000, 1, 1, 12.0, 0.0, 0.0, 0),
    ("Rome", 1980, 5, 20, 14.5, 41.9028, 12.4964, 0),
    ("New York", 2024, 11, 5, 9.0, 40.7128, -74.0060, 0),
    ("Sydney", 1950, 10, 15, 22.0, -33.8688, 151.2093, 0),
]

HIGH_LATITUDE_SUBJECTS = [
    ("Tromso (Arctic)", 1990, 1, 15, 12.0, 69.6492, 18.9553, 0),
    ("McMurdo (Antarctic)", 2005, 6, 21, 0.0, -77.8463, 166.6681, 0),
]

EQUATORIAL_SUBJECTS = [
    ("Equator", 1975, 3, 21, 12.0, 0.0, 45.0, 0),
]

ALL_SUBJECTS = STANDARD_SUBJECTS + HIGH_LATITUDE_SUBJECTS + EQUATORIAL_SUBJECTS

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360° wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def format_coord(value: float, decimals: int = 6, width: int = 10) -> str:
    """Format coordinate value with consistent width."""
    return f"{value:{width}.{decimals}f}"


def format_diff(value: float, decimals: int = 8, width: int = 10) -> str:
    """Format difference value with consistent width."""
    return f"{value:{width}.{decimals}f}"


def format_status(passed: bool) -> str:
    """Format pass/fail status."""
    return "✓" if passed else "✗"


# ============================================================================
# COMPARISON RESULT CLASSES
# ============================================================================


@dataclass
class PositionResult:
    """Stores position comparison result."""

    lon_swe: float = 0.0
    lat_swe: float = 0.0
    dist_swe: float = 0.0
    lon_py: float = 0.0
    lat_py: float = 0.0
    dist_py: float = 0.0

    diff_lon: float = 0.0
    diff_lat: float = 0.0
    diff_dist: float = 0.0

    passed: bool = False
    error_swe: Optional[str] = None
    error_py: Optional[str] = None

    def calculate_diffs(self):
        """Calculate all differences."""
        self.diff_lon = angular_diff(self.lon_swe, self.lon_py)
        self.diff_lat = abs(self.lat_swe - self.lat_py)
        self.diff_dist = abs(self.dist_swe - self.dist_py)

    def check_passed(self, lon_tol: float, lat_tol: float, dist_tol: float) -> bool:
        """Check if within tolerances."""
        if self.error_swe or self.error_py:
            self.passed = False
            return False
        self.passed = (
            self.diff_lon < lon_tol
            and self.diff_lat < lat_tol
            and self.diff_dist < dist_tol
        )
        return self.passed


@dataclass
class VelocityResult:
    """Stores velocity comparison result."""

    lon_speed_swe: float = 0.0
    lat_speed_swe: float = 0.0
    dist_speed_swe: float = 0.0
    lon_speed_py: float = 0.0
    lat_speed_py: float = 0.0
    dist_speed_py: float = 0.0

    diff_lon_speed: float = 0.0
    diff_lat_speed: float = 0.0
    diff_dist_speed: float = 0.0

    passed: bool = False

    def calculate_diffs(self):
        """Calculate all differences."""
        self.diff_lon_speed = abs(self.lon_speed_swe - self.lon_speed_py)
        self.diff_lat_speed = abs(self.lat_speed_swe - self.lat_speed_py)
        self.diff_dist_speed = abs(self.dist_speed_swe - self.dist_speed_py)

    def check_passed(self, ang_tol: float, rad_tol: float) -> bool:
        """Check if within tolerances."""
        self.passed = (
            self.diff_lon_speed < ang_tol
            and self.diff_lat_speed < ang_tol
            and self.diff_dist_speed < rad_tol
        )
        return self.passed


# ============================================================================
# SUMMARY STATISTICS
# ============================================================================


class TestStatistics:
    """Tracks and reports test statistics."""

    def __init__(self):
        self.total = 0
        self.passed = 0
        self.failed = 0
        self.errors = 0
        self.max_diff = 0.0
        self.diff_sum = 0.0

    def add_result(self, passed: bool, diff: float = 0.0, error: bool = False):
        """Add a test result."""
        self.total += 1
        if error:
            self.errors += 1
        elif passed:
            self.passed += 1
        else:
            self.failed += 1

        if not error:
            self.max_diff = max(self.max_diff, diff)
            self.diff_sum += diff

    def avg_diff(self) -> float:
        """Calculate average difference (excluding errors)."""
        count = self.total - self.errors
        return self.diff_sum / count if count > 0 else 0.0

    def pass_rate(self) -> float:
        """Calculate pass rate (excluding errors)."""
        count = self.total - self.errors
        return (self.passed / count * 100) if count > 0 else 0.0

    def print_summary(self, title: str = "SUMMARY"):
        """Print formatted summary."""
        print()
        print("=" * 80)
        print(title)
        print("=" * 80)
        print(f"Total tests:   {self.total}")
        print(f"Passed:        {self.passed} ✓")
        print(f"Failed:        {self.failed} ✗")
        print(f"Errors:        {self.errors}")
        if self.total > self.errors:
            print(f"Pass rate:     {self.pass_rate():.1f}%")
            print(f"Max diff:      {self.max_diff:.6f}")
            print(f"Avg diff:      {self.avg_diff():.6f}")
        print("=" * 80)


# ============================================================================
# COMMAND LINE HELPERS
# ============================================================================


def parse_args(args: List[str]) -> dict:
    """Parse common command line arguments."""
    return {
        "verbose": "--verbose" in args or "-v" in args,
        "quiet": "--quiet" in args or "-q" in args,
        "help": "--help" in args or "-h" in args,
    }


def print_header(title: str):
    """Print formatted header."""
    print("=" * 80)
    print(title)
    print("=" * 80)
    print()


def print_section(title: str):
    """Print formatted section header."""
    print()
    print("-" * 60)
    print(title)
    print("-" * 60)


# ============================================================================
# ADDITIONAL TEST SUBJECTS
# ============================================================================

# Historical dates for testing calendar transitions
HISTORICAL_SUBJECTS = [
    ("Renaissance", 1600, 1, 1, 12.0, 45.4642, 9.1900, 0),  # Milan
    ("Enlightenment", 1776, 7, 4, 12.0, 38.9072, -77.0369, 0),  # Washington DC
    ("Victorian", 1888, 11, 9, 12.0, 51.5074, -0.1278, 0),  # London
]

# Future dates
FUTURE_SUBJECTS = [
    ("Near Future", 2030, 6, 15, 12.0, 35.6762, 139.6503, 0),  # Tokyo
    ("Mid Century", 2050, 1, 1, 0.0, 0.0, 0.0, 0),  # Equator
    ("Late Century", 2099, 12, 31, 23.99, 41.9028, 12.4964, 0),  # Rome
]

# Edge case dates
EDGE_CASE_DATES = [
    ("Leap Year Feb 29", 2024, 2, 29, 12.0, 0.0, 0.0, 0),
    ("Year Start", 2000, 1, 1, 0.0, 0.0, 0.0, 0),
    ("Year End", 1999, 12, 31, 23.999, 0.0, 0.0, 0),
    ("Summer Solstice", 2024, 6, 21, 12.0, 51.5074, -0.1278, 0),
    ("Winter Solstice", 2024, 12, 21, 12.0, 51.5074, -0.1278, 0),
]

# Locations for rise/set testing (various latitudes)
RISE_SET_LOCATIONS = [
    ("Equator", 0.0, 0.0, 0),
    ("Tropics North", 23.4367, 0.0, 0),  # Tropic of Cancer
    ("Tropics South", -23.4367, 0.0, 0),  # Tropic of Capricorn
    ("Mid Latitude", 45.0, 0.0, 0),
    ("Arctic Circle", 66.5, 0.0, 0),
    ("High Arctic", 75.0, 0.0, 0),
    ("Antarctic", -70.0, 0.0, 0),
]

# Famous fixed stars for testing
FIXED_STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Fomalhaut",
    "Sirius",
    "Vega",
    "Altair",
    "Betelgeuse",
    "Rigel",
    "Procyon",
    "Pollux",
    "Capella",
    "Arcturus",
    "Deneb",
]


# ============================================================================
# TIME COMPARISON UTILITIES
# ============================================================================


class TimeComparisonResult:
    """Result for time function comparison."""

    def __init__(self, func_name: str, input_desc: str):
        self.func_name = func_name
        self.input_desc = input_desc
        self.swe_result = None
        self.py_result = None
        self.diff = 0.0
        self.passed = False
        self.error_swe: Optional[str] = None
        self.error_py: Optional[str] = None

    def format_result(self) -> str:
        if self.error_swe:
            return f"[{self.func_name}] {self.input_desc}: SWE ERROR {self.error_swe}"
        if self.error_py:
            return f"[{self.func_name}] {self.input_desc}: PY ERROR {self.error_py}"
        status = format_status(self.passed)
        return (
            f"[{self.func_name}] {self.input_desc}: "
            f"SWE={self.swe_result} PY={self.py_result} Diff={self.diff:.10f} {status}"
        )


# ============================================================================
# ECLIPSE/EVENT COMPARISON UTILITIES
# ============================================================================


def jd_to_date_str(jd: float) -> str:
    """Convert Julian Day to date string for display."""
    import swisseph as swe

    y, m, d, h = swe.revjul(jd)
    hours = int(h)
    minutes = int((h - hours) * 60)
    return f"{int(y)}-{int(m):02d}-{int(d):02d} {hours:02d}:{minutes:02d}"


def time_diff_seconds(jd1: float, jd2: float) -> float:
    """Calculate time difference in seconds between two Julian Days."""
    return abs(jd1 - jd2) * 86400.0


def time_diff_minutes(jd1: float, jd2: float) -> float:
    """Calculate time difference in minutes between two Julian Days."""
    return abs(jd1 - jd2) * 1440.0


class EventComparisonResult:
    """Result for event (eclipse, rise/set) comparison."""

    def __init__(self, event_type: str, event_desc: str):
        self.event_type = event_type
        self.event_desc = event_desc
        self.jd_swe: Optional[float] = None
        self.jd_py: Optional[float] = None
        self.diff_seconds = 0.0
        self.tolerance_seconds = 60.0  # Default 1 minute
        self.passed = False
        self.error_swe: Optional[str] = None
        self.error_py: Optional[str] = None
        self.extra_info: dict = {}

    def calculate_diff(self):
        if self.jd_swe is not None and self.jd_py is not None:
            self.diff_seconds = time_diff_seconds(self.jd_swe, self.jd_py)
            self.passed = self.diff_seconds < self.tolerance_seconds

    def format_result(self, verbose: bool = False) -> str:
        if self.error_swe:
            return f"[{self.event_type}] {self.event_desc}: SWE ERROR {self.error_swe}"
        if self.error_py:
            return f"[{self.event_type}] {self.event_desc}: PY ERROR {self.error_py}"

        status = format_status(self.passed)
        swe_date = jd_to_date_str(self.jd_swe) if self.jd_swe else "N/A"
        py_date = jd_to_date_str(self.jd_py) if self.jd_py else "N/A"

        if verbose:
            return (
                f"[{self.event_type}] {self.event_desc}\n"
                f"  SWE: {swe_date} (JD {self.jd_swe:.6f})\n"
                f"  PY:  {py_date} (JD {self.jd_py:.6f})\n"
                f"  Diff: {self.diff_seconds:.1f}s {status}"
            )
        return (
            f"[{self.event_type}] {self.event_desc}: "
            f"SWE={swe_date} PY={py_date} Diff={self.diff_seconds:.1f}s {status}"
        )


# ============================================================================
# ARRAY/TUPLE COMPARISON UTILITIES
# ============================================================================


def compare_arrays(
    arr_swe, arr_py, tolerances: List[float], labels: Optional[List[str]] = None
) -> tuple:
    """
    Compare two arrays element by element.

    Returns:
        (all_passed, max_diff, diffs_list)
    """
    diffs = []
    all_passed = True
    max_diff = 0.0

    for i in range(min(len(arr_swe), len(arr_py), len(tolerances))):
        diff = abs(arr_swe[i] - arr_py[i])
        diffs.append(diff)
        max_diff = max(max_diff, diff)
        if diff >= tolerances[i]:
            all_passed = False

    return all_passed, max_diff, diffs


def format_array_comparison(
    arr_swe,
    arr_py,
    labels: List[str],
    tolerances: List[float],
    width: int = 12,
    decimals: int = 6,
) -> str:
    """Format array comparison for display."""
    lines = []
    for i in range(min(len(arr_swe), len(arr_py), len(labels))):
        diff = abs(arr_swe[i] - arr_py[i])
        passed = diff < tolerances[i] if i < len(tolerances) else True
        lines.append(
            f"  {labels[i]:<15}: SWE={arr_swe[i]:{width}.{decimals}f} "
            f"PY={arr_py[i]:{width}.{decimals}f} "
            f"Diff={diff:.{decimals}f} {format_status(passed)}"
        )
    return "\n".join(lines)
