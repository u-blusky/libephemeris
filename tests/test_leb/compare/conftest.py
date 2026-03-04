"""
Shared fixtures and utilities for LEB vs Skyfield comparison tests.

This module provides:
- CompareHelper class for executing functions in both modes
- Tolerances dataclass with configurable thresholds
- Test date generators
- Helper functions for angular calculations
"""

from __future__ import annotations

import dataclasses
import os
from typing import Any, Callable, Generator

import pytest

import libephemeris as ephem
from libephemeris.time_utils import swe_julday


# =============================================================================
# TOLERANCES
# =============================================================================


@dataclasses.dataclass
class Tolerances:
    """Configurable tolerance thresholds for LEB vs Skyfield comparison.

    All values can be overridden via environment variables with the prefix
    LEB_TOL_, e.g. LEB_TOL_POSITION_ARCSEC=0.5

    Note: These tolerances reflect the actual LEB Chebyshev approximation precision.
    LEB precision varies by body type and can be higher for outer planets.
    See tests/test_leb/test_leb_precision.py for reference values.

    Typical precision by body type:
    - Inner planets (Mercury-Venus): ~0.1-0.2 arcsec
    - Outer planets (Mars-Pluto): ~0.1-5 arcsec (varies by date)
    - Moon: ~0.1 arcsec position, ~0.01 deg/day speed
    - Asteroids: ~0.1 arcsec
    - Ecliptic bodies (nodes, Lilith): ~0.5 arcsec
    """

    # Core position (ICRS planets can vary 0.1-5 arcsec depending on body/date)
    POSITION_ARCSEC: float = 5.0
    ECLIPTIC_ARCSEC: float = 0.5
    HYPOTHETICAL_ARCSEC: float = 0.5
    EQUATORIAL_ARCSEC: float = 5.0
    J2000_ARCSEC: float = 5.0
    SIDEREAL_ARCSEC: float = 20.0  # Some modes (41, 255) have larger errors
    DISTANCE_AU: float = 1e-4

    # Velocity (ecliptic bodies like OscuApogee can have higher speed errors)
    SPEED_LON_DEG_DAY: float = 0.05
    SPEED_LAT_DEG_DAY: float = 0.05
    SPEED_DIST_AU_DAY: float = 1e-4

    # Ayanamsha
    AYANAMSHA_ARCSEC: float = 0.001

    # Timing (indirect callers)
    CROSSING_SUN_SEC: float = 1.0
    CROSSING_MOON_SEC: float = 5.0
    CROSSING_MOON_NODE_SEC: float = 10.0
    CROSSING_PLANET_SEC: float = 30.0
    ECLIPSE_TIMING_SEC: float = 1.0
    ECLIPSE_MAGNITUDE: float = 1e-4
    ECLIPSE_POSITION_DEG: float = 0.001
    ELONGATION_ARCSEC: float = 0.5  # Elongation uses two planet positions
    STATION_TIMING_SEC: float = 1.0
    RISE_TRANSIT_SEC: float = 1.0
    HOUSE_SUNSHINE_ARCSEC: float = 0.001
    GAUQUELIN_SECTOR: float = 0.001

    @classmethod
    def from_env(cls) -> "Tolerances":
        """Load tolerance overrides from environment variables."""
        kwargs = {}
        for field in dataclasses.fields(cls):
            env_name = f"LEB_TOL_{field.name.upper()}"
            env_val = os.environ.get(env_name)
            if env_val is not None:
                kwargs[field.name] = float(env_val)
        return cls(**kwargs)


TOLS = Tolerances.from_env()


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360 wrap."""
    d = abs(a - b)
    if d > 180.0:
        d = 360.0 - d
    return d


def lon_error_arcsec(a: float, b: float) -> float:
    """Longitude error in arcseconds with wrap-around."""
    return angular_diff(a, b) * 3600.0


def generate_test_dates(
    n: int, jd_start: float, jd_end: float, margin: float = 30.0
) -> list[float]:
    """Generate n uniformly-spaced dates within [start+margin, end-margin]."""
    effective_start = jd_start + margin
    effective_end = jd_end - margin
    if n <= 1:
        return [(effective_start + effective_end) / 2.0]
    step = (effective_end - effective_start) / (n - 1)
    return [effective_start + i * step for i in range(n)]


def year_to_jd(year: int) -> float:
    """Convert a year to Julian Day (January 1.0)."""
    return swe_julday(year, 1, 1, 0.0)


# =============================================================================
# COMPARE HELPER
# =============================================================================


class CompareHelper:
    """Execute functions in both Skyfield and LEB mode, saving/restoring state.

    Note: Uses set_calc_mode("skyfield") to prevent LEB auto-discovery from
    interfering with Skyfield reference calculations.
    """

    def __init__(self, leb_path: str):
        self.leb_path = leb_path
        self._saved_leb: str | None = None
        self._saved_mode: str | None = None
        self._saved_reader: Any = None

    def setup(self) -> None:
        """Save current global state."""
        self._saved_leb = ephem.state._LEB_FILE
        self._saved_mode = ephem.state._CALC_MODE
        self._saved_reader = ephem.state._LEB_READER

    def teardown(self) -> None:
        """Restore saved global state."""
        ephem.state._LEB_FILE = self._saved_leb
        ephem.state._LEB_READER = self._saved_reader
        if self._saved_mode is not None:
            ephem.set_calc_mode(self._saved_mode)
        else:
            ephem.set_calc_mode(None)

    def skyfield(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn in forced Skyfield mode (no LEB auto-discovery)."""
        ephem.state._LEB_FILE = None
        ephem.state._LEB_READER = None
        ephem.set_calc_mode("skyfield")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.set_calc_mode(None)

    def leb(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn in LEB mode with explicit file path."""
        ephem.state._LEB_FILE = self.leb_path
        ephem.state._LEB_READER = None
        ephem.set_calc_mode("auto")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.state._LEB_FILE = None
            ephem.state._LEB_READER = None
            ephem.set_calc_mode(None)


# =============================================================================
# FIXTURES
# =============================================================================


_LEB_FILE_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "data", "leb", "ephemeris_medium.leb"
)


@pytest.fixture(scope="session")
def leb_file() -> str:
    """Path to the pre-generated medium tier LEB file."""
    path = os.path.abspath(_LEB_FILE_PATH)
    if not os.path.exists(path):
        pytest.skip(f"LEB file not found: {path}")
    return path


@pytest.fixture
def compare(leb_file: str) -> Generator[CompareHelper, None, None]:
    """CompareHelper instance with automatic state save/restore."""
    helper = CompareHelper(leb_file)
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()


# Date fixtures for different test intensities
_MEDIUM_TIER_START = year_to_jd(1560)
_MEDIUM_TIER_END = year_to_jd(2640)


@pytest.fixture(scope="session")
def test_dates_200() -> list[float]:
    """200 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(200, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


@pytest.fixture(scope="session")
def test_dates_100() -> list[float]:
    """100 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(100, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


@pytest.fixture(scope="session")
def test_dates_50() -> list[float]:
    """50 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(50, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


@pytest.fixture(scope="session")
def test_dates_20() -> list[float]:
    """20 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(20, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


# =============================================================================
# BODY LISTS
# =============================================================================

ICRS_PLANETS = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
    (14, "Earth"),
]

ECLIPTIC_BODIES = [
    (10, "MeanNode"),
    (11, "TrueNode"),
    (12, "MeanApogee"),
    (13, "OscuApogee"),
    (21, "InterpApogee"),
    (22, "InterpPerigee"),
]

ASTEROID_BODIES = [
    (15, "Chiron"),
    (17, "Ceres"),
    (18, "Pallas"),
    (19, "Juno"),
    (20, "Vesta"),
]

HYPOTHETICAL_BODIES = [
    (40, "Cupido"),
    (41, "Hades"),
    (42, "Zeus"),
    (43, "Kronos"),
    (44, "Apollon"),
    (45, "Admetos"),
    (46, "Vulkanus"),
    (47, "Poseidon"),
    (48, "Transpluto"),
]

# Ecliptic body tolerances (per-body)
ECLIPTIC_TOLERANCES = {
    10: {"lon": 0.01, "speed": 0.0001},  # Mean Node
    11: {"lon": 0.5, "speed": 0.01},  # True Node
    12: {"lon": 0.01, "speed": 0.0001},  # Mean Apogee
    13: {"lon": 0.5, "speed": 0.05},  # Oscu Apogee (higher speed variance)
    21: {"lon": 0.5, "speed": 0.01},  # Interp Apogee
    22: {"lon": 0.5, "speed": 0.01},  # Interp Perigee
}

# Formula-based sidereal modes (27 modes + user-defined)
FORMULA_SIDEREAL_MODES = [
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    37,
    38,
    41,
    255,  # user-defined
]

# Star-based sidereal modes (fallback to Skyfield)
STAR_BASED_SIDEREAL_MODES = [17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42]
