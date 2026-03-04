# LEB vs Skyfield Comparison Test Suite: Implementation Plan

> **STATUS: IMPLEMENTATION PLAN**
> This document specifies the complete test suite for validating LEB binary
> ephemeris precision against the Skyfield reference pipeline.

---

## Table of Contents

1. [Motivation and Scope](#1-motivation-and-scope)
2. [Architecture](#2-architecture)
3. [Tolerances](#3-tolerances)
4. [Shared Infrastructure (`conftest.py`)](#4-shared-infrastructure)
5. [Test File Specifications](#5-test-file-specifications)
6. [Excluded Areas](#6-excluded-areas)
7. [Execution and CI](#7-execution-and-ci)
8. [Implementation Order](#8-implementation-order)
9. [File Inventory](#9-file-inventory)

---

## 1. Motivation and Scope

### 1.1 Problem

The existing Swiss vs libephemeris comparison suite (`tests/test_compare_*.py`, 25+
primary files, 100+ total) validates that libephemeris reproduces pyswisseph results.
However, there is no equivalent systematic validation that **LEB mode reproduces
Skyfield mode** across all affected functions.

The existing LEB tests (`tests/test_leb/test_fast_calc.py`, `test_leb_precision.py`)
cover raw `fast_calc_ut()` precision but:

- Use `fast_calc_ut()` directly, not the public API
- Don't test high-level functions that call `swe_calc_ut()` internally
- Don't cover all 30 LEB bodies systematically with all flag combinations
- Don't test indirect callers (crossings, eclipses, rise/transit, stations)

### 1.2 Goal

Create a comprehensive test suite that validates **every function affected by LEB
mode** produces results matching Skyfield mode within strict, configurable tolerances.

### 1.3 LEB Impact Analysis

LEB dispatches in 3 functions:

| Function | File | Line |
|----------|------|------|
| `swe_calc_ut()` | `planets.py` | 772 |
| `swe_calc()` | `planets.py` | 847 |
| `swe_get_ayanamsa_ut()` | `planets.py` | 2276 |

These 3 entry points are called by **62 functions** across 6 files:

| File | Direct callers | Indirect callers | Total |
|------|---------------|-------------------|-------|
| `planets.py` | 3 | 4 | 7 |
| `crossing.py` | 15 | 1 | 16 |
| `eclipse.py` | 11 | 20 | 31 |
| `houses.py` | 3 | 4 | 7 |
| `utils.py` | 1 | 0 | 1 |
| `time_utils.py` | 1 | 0 | 1 |
| **Total** | **34** | **28** | **62** |

### 1.4 Bodies Supported in LEB (30 total)

**Pipeline A -- ICRS Barycentric (16 bodies):**

| ID | Body | Interval | Degree |
|----|------|----------|--------|
| 0 | Sun | 32 days | 13 |
| 1 | Moon | 8 days | 13 |
| 2 | Mercury | 16 days | 15 |
| 3 | Venus | 32 days | 13 |
| 4 | Mars | 32 days | 13 |
| 5 | Jupiter | 64 days | 11 |
| 6 | Saturn | 64 days | 11 |
| 7 | Uranus | 128 days | 9 |
| 8 | Neptune | 128 days | 9 |
| 9 | Pluto | 128 days | 9 |
| 14 | Earth | 8 days | 13 |
| 15 | Chiron | 32 days | 13 |
| 17 | Ceres | 32 days | 13 |
| 18 | Pallas | 32 days | 13 |
| 19 | Juno | 32 days | 13 |
| 20 | Vesta | 32 days | 13 |

**Pipeline B -- Ecliptic Direct (6 bodies):**

| ID | Body | Interval | Degree |
|----|------|----------|--------|
| 10 | Mean Node | 8 days | 13 |
| 11 | True Node | 8 days | 13 |
| 12 | Mean Apogee (Lilith) | 8 days | 13 |
| 13 | Osculating Apogee | 8 days | 13 |
| 21 | Interpolated Apogee | 8 days | 13 |
| 22 | Interpolated Perigee | 8 days | 13 |

**Pipeline C -- Heliocentric Ecliptic (9 bodies):**

| ID | Body | Interval | Degree |
|----|------|----------|--------|
| 40 | Cupido | 64 days | 11 |
| 41 | Hades | 64 days | 11 |
| 42 | Zeus | 64 days | 11 |
| 43 | Kronos | 64 days | 11 |
| 44 | Apollon | 64 days | 11 |
| 45 | Admetos | 64 days | 11 |
| 46 | Vulkanus | 64 days | 11 |
| 47 | Poseidon | 64 days | 11 |
| 48 | Isis/Transpluto | 64 days | 11 |

### 1.5 Flags Handled Natively by LEB

| Flag | Handling |
|------|----------|
| `SEFLG_SPEED` | Central difference velocity (Pipeline A) or stored velocity (B/C) |
| `SEFLG_HELCTR` | Observer = Sun position from LEB |
| `SEFLG_BARYCTR` | Observer = origin (0,0,0) |
| `SEFLG_EQUATORIAL` | Precession-nutation matrix rotation |
| `SEFLG_J2000` | J2000 ecliptic or J2000 equatorial (combined with EQUATORIAL) |
| `SEFLG_SIDEREAL` | Formula-based modes only (27 modes + user-defined mode 255) |
| `SEFLG_TRUEPOS` | Skips light-time correction and aberration |
| `SEFLG_NOABERR` | Skips aberration only |
| `SEFLG_MOSEPH` | Silently stripped |

### 1.6 Flags that Trigger Fallback to Skyfield

| Flag | Reason |
|------|--------|
| `SEFLG_TOPOCTR` | Requires observer location + Earth rotation |
| `SEFLG_XYZ` | Returns Cartesian instead of spherical |
| `SEFLG_RADIANS` | Returns radians instead of degrees |
| `SEFLG_NONUT` | Mean ecliptic (skip nutation) |

Star-based sidereal modes (17, 27-36, 39, 40, 42) also fall back.

---

## 2. Architecture

### 2.1 Directory Layout

```
tests/test_leb/compare/
    __init__.py
    conftest.py                          # Shared infrastructure
    test_compare_leb_planets.py          # 1. ICRS planet positions
    test_compare_leb_lunar.py            # 2. Lunar/ecliptic bodies
    test_compare_leb_asteroids.py        # 3. Main-belt asteroids
    test_compare_leb_hypothetical.py     # 4. Uranian planets + Transpluto
    test_compare_leb_observations.py     # 5. All flag combinations
    test_compare_leb_sidereal.py         # 6. Sidereal mode (formula-based)
    test_compare_leb_ayanamsha.py        # 7. Ayanamsha values (direct API)
    test_compare_leb_velocities.py       # 8. Velocity precision (all bodies)
    test_compare_leb_distances.py        # 9. Distance precision
    test_compare_leb_crossings.py        # 10. Crossing functions
    test_compare_leb_eclipses_solar.py   # 11. Solar eclipses
    test_compare_leb_eclipses_lunar.py   # 12. Lunar eclipses
    test_compare_leb_elongation.py       # 13. Elongation helpers
    test_compare_leb_stations.py         # 14. Stations + retrogrades
    test_compare_leb_houses.py           # 15. Houses (Sunshine system)
    test_compare_leb_gauquelin.py        # 16. Gauquelin sectors
    test_compare_leb_rise_transit.py     # 17. Rise/transit/set
```

### 2.2 Testing Pattern

Every test follows the same pattern:

1. Force **Skyfield mode** via `set_calc_mode("skyfield")` -- compute reference result
2. Activate **LEB mode** via `set_leb_file(path)` + `set_calc_mode("auto")` -- compute LEB result
3. Compare within tolerance
4. Restore original state in `finally` block

```python
def test_example(self, compare):
    """Example pattern."""
    ref = compare.skyfield(ephem.swe_calc_ut, jd, SE_SUN, SEFLG_SPEED)
    leb = compare.leb(ephem.swe_calc_ut, jd, SE_SUN, SEFLG_SPEED)

    lon_err = angular_diff(ref[0][0], leb[0][0]) * 3600.0
    assert lon_err < Tolerances.POSITION_ARCSEC
```

### 2.3 Pytest Markers

All tests use `@pytest.mark.leb_compare` for selective execution. Additionally:

| Marker | Purpose |
|--------|---------|
| `@pytest.mark.leb_compare` | All LEB comparison tests |
| `@pytest.mark.slow` | Tests with 100+ dates or heavy computation |
| `@pytest.mark.parametrize` | Body/flag/date parametrization |

### 2.4 Data Source

Tests use the **pre-generated** `.leb` file at `data/leb/ephemeris_medium.leb`
(de440, 1550-2650). No on-the-fly generation. If the file is missing, the entire
suite is skipped with `pytest.skip()`.

---

## 3. Tolerances

### 3.1 Design Principle

LEB is a Chebyshev polynomial approximation of the **same Skyfield pipeline**.
Unlike Swiss vs libephemeris (different algorithms), LEB vs Skyfield differences
are purely numerical approximation error. Tolerances should be **much stricter**
than the Swiss comparison suite.

All tolerances are defined in a single `Tolerances` class in `conftest.py` for
easy adjustment.

### 3.2 Tolerance Table

#### Core Position/Velocity (Files 1-9)

| Metric | Tolerance | Unit | Notes |
|--------|-----------|------|-------|
| `POSITION_ARCSEC` | 0.1 | arcsec | ICRS planets lon/lat |
| `ECLIPTIC_ARCSEC` | 0.5 | arcsec | Pipeline B bodies (nodes, Lilith) |
| `HYPOTHETICAL_ARCSEC` | 0.5 | arcsec | Pipeline C bodies (Uranians) |
| `EQUATORIAL_ARCSEC` | 0.2 | arcsec | After equatorial transform |
| `J2000_ARCSEC` | 0.2 | arcsec | After J2000 precession |
| `SIDEREAL_ARCSEC` | 0.01 | arcsec | Sidereal positions (formula-based) |
| `DISTANCE_AU` | 1e-6 | AU | Geocentric distance |
| `SPEED_LON_DEG_DAY` | 0.001 | deg/day | Longitude velocity |
| `SPEED_LAT_DEG_DAY` | 0.0001 | deg/day | Latitude velocity |
| `SPEED_DIST_AU_DAY` | 1e-6 | AU/day | Distance velocity |
| `AYANAMSHA_ARCSEC` | 0.001 | arcsec | Direct ayanamsha values |

#### Indirect Callers -- Timing (Files 10-17)

| Metric | Tolerance | Unit | Notes |
|--------|-----------|------|-------|
| `CROSSING_SUN_SEC` | 1.0 | seconds | Solar longitude crossings |
| `CROSSING_MOON_SEC` | 5.0 | seconds | Lunar longitude crossings |
| `CROSSING_MOON_NODE_SEC` | 10.0 | seconds | Lunar node crossings |
| `CROSSING_PLANET_SEC` | 30.0 | seconds | Planet crossings |
| `ECLIPSE_TIMING_SEC` | 1.0 | seconds | Eclipse global/local timing |
| `ECLIPSE_MAGNITUDE` | 1e-4 | ratio | Eclipse magnitude |
| `ECLIPSE_POSITION_DEG` | 0.001 | degrees | Eclipse central line position |
| `ELONGATION_ARCSEC` | 0.001 | arcsec | Elongation angle |
| `STATION_TIMING_SEC` | 1.0 | seconds | Station timing |
| `RISE_TRANSIT_SEC` | 1.0 | seconds | Rise/transit/set timing |
| `HOUSE_CUSP_DEG` | 0.0 | degrees | Non-Sunshine cusps (must be identical) |
| `HOUSE_SUNSHINE_ARCSEC` | 0.001 | arcsec | Sunshine system cusps |
| `GAUQUELIN_SECTOR` | 0.001 | sector | Gauquelin sector number |

### 3.3 Per-Body Ecliptic Tolerances (File 2)

Pipeline B bodies have varying precision due to different orbital characteristics:

| Body | Position (arcsec) | Speed (deg/day) | Notes |
|------|------------------|-----------------|-------|
| Mean Node (10) | 0.01 | 0.0001 | Smooth, linear-ish motion |
| True Node (11) | 0.5 | 0.01 | Oscillatory |
| Mean Apogee (12) | 0.01 | 0.0001 | Smooth |
| Oscu Apogee (13) | 0.5 | 0.01 | Oscillatory |
| Interp Apogee (21) | 0.5 | 0.01 | Interpolated, less smooth |
| Interp Perigee (22) | 0.5 | 0.01 | Interpolated, less smooth |

### 3.4 Tolerance Escalation Strategy

If any test fails at the ultra-strict level:

1. Log the max error, body, JD, and flag combination
2. Check if the error is consistent (systematic) or sporadic (edge case)
3. If systematic: investigate Chebyshev fitting parameters in `leb_format.py`
4. If sporadic (near segment boundaries): relax tolerance for that body only
5. Document any relaxation in the test file with a comment explaining why

---

## 4. Shared Infrastructure

### 4.1 `conftest.py`

```
tests/test_leb/compare/conftest.py
```

#### Fixtures

| Fixture | Scope | Description |
|---------|-------|-------------|
| `leb_file` | `session` | Path to `data/leb/ephemeris_medium.leb`. Skips if missing. |
| `compare` | `function` | `CompareHelper` instance. Saves/restores global state. |
| `test_dates_200` | `session` | 200 uniformly-spaced JDs across 1560-2640 |
| `test_dates_100` | `session` | 100 uniformly-spaced JDs across 1560-2640 |
| `test_dates_50` | `session` | 50 uniformly-spaced JDs across 1560-2640 |
| `test_dates_20` | `session` | 20 uniformly-spaced JDs across 1560-2640 |

#### `CompareHelper` Class

```python
class CompareHelper:
    """Execute functions in both Skyfield and LEB mode, saving/restoring state."""

    def __init__(self, leb_path: str):
        self.leb_path = leb_path
        self._saved_leb = None
        self._saved_mode = None

    def setup(self):
        """Save current global state."""
        self._saved_leb = ephem.state._LEB_FILE
        self._saved_mode = ephem.state._CALC_MODE

    def teardown(self):
        """Restore saved global state."""
        ephem.set_leb_file(self._saved_leb)
        ephem.state._CALC_MODE = self._saved_mode

    def skyfield(self, fn, *args, **kwargs):
        """Call fn in forced Skyfield mode."""
        ephem.set_leb_file(None)
        ephem.set_calc_mode("skyfield")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.set_calc_mode(None)

    def leb(self, fn, *args, **kwargs):
        """Call fn in LEB mode."""
        ephem.set_leb_file(self.leb_path)
        ephem.set_calc_mode("auto")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.set_leb_file(None)
            ephem.set_calc_mode(None)
```

#### Helper Functions

```python
def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360 wrap."""
    d = abs(a - b)
    if d > 180.0:
        d = 360.0 - d
    return d

def lon_error_arcsec(a: float, b: float) -> float:
    """Longitude error in arcseconds with wrap-around."""
    return angular_diff(a, b) * 3600.0

def generate_test_dates(n: int, jd_start: float, jd_end: float,
                        margin: float = 30.0) -> list[float]:
    """Generate n uniformly-spaced dates within [start+margin, end-margin]."""
    ...
```

#### `Tolerances` Class

```python
@dataclasses.dataclass
class Tolerances:
    """Configurable tolerance thresholds for LEB vs Skyfield comparison.

    All values can be overridden via environment variables with the prefix
    LEB_TOL_, e.g. LEB_TOL_POSITION_ARCSEC=0.5
    """
    # Core position
    POSITION_ARCSEC: float = 0.1
    ECLIPTIC_ARCSEC: float = 0.5
    HYPOTHETICAL_ARCSEC: float = 0.5
    EQUATORIAL_ARCSEC: float = 0.2
    J2000_ARCSEC: float = 0.2
    SIDEREAL_ARCSEC: float = 0.01
    DISTANCE_AU: float = 1e-6
    # Velocity
    SPEED_LON_DEG_DAY: float = 0.001
    SPEED_LAT_DEG_DAY: float = 0.0001
    SPEED_DIST_AU_DAY: float = 1e-6
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
    ELONGATION_ARCSEC: float = 0.001
    STATION_TIMING_SEC: float = 1.0
    RISE_TRANSIT_SEC: float = 1.0
    HOUSE_SUNSHINE_ARCSEC: float = 0.001
    GAUQUELIN_SECTOR: float = 0.001

    @classmethod
    def from_env(cls) -> "Tolerances":
        """Load tolerance overrides from environment variables."""
        ...
```

---

## 5. Test File Specifications

### 5.1 `test_compare_leb_planets.py` -- ICRS Planet Positions

**Counterpart:** `tests/test_compare_planets.py`

**Purpose:** Validate ecliptic longitude, latitude, and distance for all 11 ICRS
pipeline planets across the full medium tier range.

| Parameter | Value |
|-----------|-------|
| Bodies | Sun(0), Moon(1), Mercury(2), Venus(3), Mars(4), Jupiter(5), Saturn(6), Uranus(7), Neptune(8), Pluto(9), Earth(14) |
| Dates | 200 uniformly-spaced, 1560-2640 |
| Flags | `SEFLG_SPEED` |
| Tolerances | lon: 0.1", lat: 0.1", dist: 1e-6 AU |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestPlanetLongitude` | 11 parametrized (1 per planet) | Max longitude error across 200 dates |
| `TestPlanetLatitude` | 11 parametrized | Max latitude error across 200 dates |
| `TestPlanetDistance` | 11 parametrized | Max distance error across 200 dates |
| `TestPlanetAllComponents` | 11 parametrized | Combined check of all 6 output components |
| `TestPlanetStatistics` | 1 | Report mean/max/p99 error per planet (informational) |

**Statistical reporting:** Each test class collects `max_err`, `mean_err`,
`worst_jd` and includes them in the assertion message for debugging.

---

### 5.2 `test_compare_leb_lunar.py` -- Lunar/Ecliptic Bodies

**Counterpart:** `tests/test_compare_lunar_nodes_lilith.py`

**Purpose:** Validate all 6 Pipeline B (ecliptic-direct) bodies with per-body
tolerances reflecting their different orbital characteristics.

| Parameter | Value |
|-----------|-------|
| Bodies | Mean Node(10), True Node(11), Mean Apogee(12), Oscu Apogee(13), Interp Apogee(21), Interp Perigee(22) |
| Dates | 200 uniformly-spaced, 1560-2640 |
| Flags | `SEFLG_SPEED` (default) + `SEFLG_EQUATORIAL \| SEFLG_SPEED` |
| Tolerances | Per-body, see Section 3.3 |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestLunarLongitude` | 6 parametrized (1 per body) | Position with per-body tolerance |
| `TestLunarLatitude` | 6 parametrized | Latitude check |
| `TestLunarSpeed` | 6 parametrized | Velocity with per-body tolerance |
| `TestLunarEquatorial` | 6 parametrized | Equatorial transform of ecliptic bodies |
| `TestLunarDistance` | 6 parametrized | Distance (where meaningful) |

---

### 5.3 `test_compare_leb_asteroids.py` -- Main-Belt Asteroids

**Counterpart:** `tests/test_compare_minor_bodies.py`

**Purpose:** Validate the 5 asteroids in LEB (all Pipeline A).

| Parameter | Value |
|-----------|-------|
| Bodies | Chiron(15), Ceres(17), Pallas(18), Juno(19), Vesta(20) |
| Dates | 200 uniformly-spaced, 1560-2640 |
| Flags | `SEFLG_SPEED` |
| Tolerances | lon/lat: 0.1", speed: 0.001 deg/day, dist: 1e-6 AU |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestAsteroidPosition` | 5 parametrized | Longitude + latitude |
| `TestAsteroidSpeed` | 5 parametrized | Velocity |
| `TestAsteroidDistance` | 5 parametrized | Distance |

---

### 5.4 `test_compare_leb_hypothetical.py` -- Uranian Planets

**Counterpart:** `tests/test_compare_hypothetical.py`

**Purpose:** Validate all 9 Pipeline C (heliocentric analytical) bodies.

| Parameter | Value |
|-----------|-------|
| Bodies | Cupido(40)-Poseidon(47), Isis/Transpluto(48) |
| Dates | 200 uniformly-spaced, 1560-2640 |
| Flags | `SEFLG_SPEED` |
| Tolerances | lon: 0.5", lat: 0.5", speed: 0.01 deg/day |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestHypotheticalPosition` | 9 parametrized | Longitude + latitude |
| `TestHypotheticalSpeed` | 9 parametrized | Velocity |

---

### 5.5 `test_compare_leb_observations.py` -- Flag Combinations

**Counterpart:** `tests/test_compare_observations.py`

**Purpose:** Validate all flag combinations supported natively by LEB. This is
the most important test for correctness of the coordinate transform pipeline.

| Parameter | Value |
|-----------|-------|
| Bodies | Sun(0), Moon(1), Mars(4), Jupiter(5) |
| Dates | 50 uniformly-spaced, 1560-2640 |
| Flags | 8 combinations (see below) |
| Tolerances | 0.2" for transformed coordinates |

**Flag Combinations:**

| # | Flags | Description |
|---|-------|-------------|
| 1 | `SEFLG_SPEED` | Default ecliptic |
| 2 | `SEFLG_SPEED \| SEFLG_EQUATORIAL` | Equatorial (RA/Dec) |
| 3 | `SEFLG_SPEED \| SEFLG_J2000` | J2000 ecliptic |
| 4 | `SEFLG_SPEED \| SEFLG_EQUATORIAL \| SEFLG_J2000` | J2000 equatorial (ICRS) |
| 5 | `SEFLG_SPEED \| SEFLG_HELCTR` | Heliocentric |
| 6 | `SEFLG_SPEED \| SEFLG_BARYCTR` | Barycentric |
| 7 | `SEFLG_SPEED \| SEFLG_TRUEPOS` | True (geometric) position |
| 8 | `SEFLG_SPEED \| SEFLG_NOABERR` | No aberration correction |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestFlagCombinations` | 4 bodies x 8 flags = 32 parametrized | All 6 components checked |
| `TestFlagCombinationsExtended` | 4 bodies x compound flags | Helio+equatorial, bary+J2000, etc. |
| `TestFallbackFlags` | 4 tests | Verify `SEFLG_TOPOCTR`, `SEFLG_XYZ`, `SEFLG_RADIANS`, `SEFLG_NONUT` produce identical results (both fallback to Skyfield) |

---

### 5.6 `test_compare_leb_sidereal.py` -- Sidereal Mode

**Counterpart:** `tests/test_compare_sidereal.py`

**Purpose:** Validate sidereal positions for all 27 formula-based ayanamsha modes.

| Parameter | Value |
|-----------|-------|
| Bodies | Sun(0), Moon(1), Mars(4), Jupiter(5) |
| Dates | 20 uniformly-spaced, 1900-2100 (ayanamsha accuracy range) |
| Flags | `SEFLG_SPEED \| SEFLG_SIDEREAL` |
| Ayanamshas | 27 formula-based modes (0-16, 18-26, 37, 38, 41, 255) |
| Tolerances | 0.01" position, 0.001 deg/day speed |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestSiderealPosition` | 4 bodies x 27 modes = 108 parametrized | Sidereal longitude |
| `TestSiderealSpeed` | 4 bodies x 27 modes = 108 parametrized | Sidereal speed (includes precession correction) |
| `TestStarBasedFallback` | 1 | Verify star-based modes (17, 27-36, 39, 40, 42) produce identical results (both fallback to Skyfield) |

---

### 5.7 `test_compare_leb_ayanamsha.py` -- Ayanamsha Values

**Counterpart:** Part of `tests/test_compare_sidereal.py`

**Purpose:** Validate `swe_get_ayanamsa_ut()` LEB dispatch directly.

| Parameter | Value |
|-----------|-------|
| Ayanamshas | 27 formula-based modes |
| Dates | 50 uniformly-spaced, 1560-2640 |
| Tolerances | 0.001" |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestAyanamshaValues` | 27 parametrized | Direct value comparison |
| `TestAyanamshaConsistency` | 4 | Verify `swe_get_ayanamsa_ut(jd, sid_mode=X)` matches `swe_calc_ut(jd, body, SEFLG_SIDEREAL)` offset |

---

### 5.8 `test_compare_leb_velocities.py` -- Velocity Precision

**Counterpart:** Velocity aspects of multiple `test_compare_*` files

**Purpose:** Dedicated, exhaustive velocity validation for all 30 LEB bodies.

| Parameter | Value |
|-----------|-------|
| Bodies | All 30 LEB bodies (16 ICRS + 6 ecliptic + 5 asteroids + 9 helio, filtered by availability) |
| Dates | 100 uniformly-spaced, 1560-2640 |
| Flags | `SEFLG_SPEED` |
| Tolerances | lon: 0.001 deg/day, lat: 0.0001 deg/day, dist: 1e-6 AU/day |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestLongitudeVelocity` | 30 parametrized | `result[3]` comparison |
| `TestLatitudeVelocity` | 30 parametrized | `result[4]` comparison |
| `TestDistanceVelocity` | 30 parametrized (ICRS + asteroids only) | `result[5]` comparison |
| `TestVelocityEquatorial` | 11 parametrized (ICRS planets) | Speed in equatorial coordinates |

---

### 5.9 `test_compare_leb_distances.py` -- Distance Precision

**Purpose:** Dedicated distance validation for geocentric bodies.

| Parameter | Value |
|-----------|-------|
| Bodies | 11 ICRS planets + Chiron(15), Ceres(17) |
| Dates | 100 uniformly-spaced, 1560-2640 |
| Flags | `SEFLG_SPEED` |
| Tolerances | 1e-6 AU |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestGeocentricDistance` | 13 parametrized | `result[2]` comparison |
| `TestHeliocentricDistance` | 5 parametrized (Mars-Saturn) | Distance with `SEFLG_HELCTR` |
| `TestDistanceStatistics` | 1 | Report distribution of errors |

---

### 5.10 `test_compare_leb_crossings.py` -- Crossing Functions

**Counterpart:** `tests/test_compare_crossings.py`

**Purpose:** Validate that crossing functions (iterative solvers using `swe_calc_ut`
internally) produce consistent timing in both modes.

| Parameter | Value |
|-----------|-------|
| Functions | `swe_solcross_ut`, `swe_mooncross_ut`, `swe_mooncross_node_ut`, `swe_cross_ut`, `swe_helio_cross_ut` |
| Tolerances | Sun: 1s, Moon: 5s, Moon node: 10s, Planet: 30s |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestSunCrossings` | 7 parametrized (target longitudes 0-180) | `swe_solcross_ut` |
| `TestMoonCrossings` | 4 parametrized (0, 90, 180, 270) | `swe_mooncross_ut` |
| `TestMoonNodeCrossings` | 2 x 4 = 8 parametrized | Ascending + descending, 4 start dates |
| `TestPlanetCrossings` | 3 planets x 4 longitudes = 12 | `swe_cross_ut` for Mars, Jupiter, Saturn |
| `TestHelioCrossings` | 3 planets x 4 longitudes = 12 | `swe_helio_cross_ut` |
| `TestCrossingPositionVerify` | 7 | Verify planet is at target longitude at crossing time |

---

### 5.11 `test_compare_leb_eclipses_solar.py` -- Solar Eclipses

**Counterpart:** `tests/test_compare_eclipses.py`

**Purpose:** Validate solar eclipse search and circumstances.

| Parameter | Value |
|-----------|-------|
| Functions | `sol_eclipse_when_glob`, `swe_sol_eclipse_when_loc`, `swe_sol_eclipse_where`, `swe_sol_eclipse_how` |
| Date range | 2024-2034 (5-6 eclipses) |
| Locations | Rome, New York, Sydney, Equator |
| Tolerances | timing: 1s, position: 0.001 deg, magnitude: 1e-4 |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestSolarEclipseGlobal` | 5-6 eclipses | Global search timing |
| `TestSolarEclipseLocal` | 4 locations x 3 eclipses | Local circumstances |
| `TestSolarEclipseWhere` | 3 eclipses | Central line coordinates |
| `TestSolarEclipseHow` | 3 eclipses x 4 locations | Magnitude and type |
| `TestSolarEclipseConsistency` | 1 | Type/flag agreement between modes |

---

### 5.12 `test_compare_leb_eclipses_lunar.py` -- Lunar Eclipses

**Counterpart:** `tests/test_compare_eclipses.py`

**Purpose:** Validate lunar eclipse search and circumstances.

| Parameter | Value |
|-----------|-------|
| Functions | `lun_eclipse_when`, `lun_eclipse_when_loc`, `swe_lun_eclipse_how` |
| Date range | 2024-2034 |
| Locations | Rome, New York, Sydney, Equator |
| Tolerances | timing: 1s, magnitude: 1e-4 |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestLunarEclipseGlobal` | 5-6 eclipses | Global search timing |
| `TestLunarEclipseLocal` | 4 locations x 3 eclipses | Local visibility |
| `TestLunarEclipseHow` | 3 eclipses | Magnitude and type |
| `TestLunarEclipseContacts` | 3 eclipses | P1/P4, U1-U4 contact times |

---

### 5.13 `test_compare_leb_elongation.py` -- Elongation Helpers

**Counterpart:** `tests/test_compare_elongation.py`

**Purpose:** Validate elongation functions that call `swe_calc_ut` internally.

| Parameter | Value |
|-----------|-------|
| Functions | `get_elongation_from_sun`, `is_morning_star`, `is_evening_star`, `get_elongation_type` |
| Bodies | Mercury(2), Venus(3), Mars(4), Jupiter(5), Saturn(6) |
| Dates | 50 uniformly-spaced, 2020-2030 |
| Tolerances | angle: 0.001", boolean/type: must be identical |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestElongationAngle` | 5 bodies x 50 dates | Numerical elongation value |
| `TestElongationClassification` | 5 bodies x 50 dates | Morning/evening/type must match exactly |

---

### 5.14 `test_compare_leb_stations.py` -- Stations and Retrogrades

**Counterpart:** New (no direct Swiss counterpart for stations)

**Purpose:** Validate station-finding and retrograde detection.

| Parameter | Value |
|-----------|-------|
| Functions | `swe_find_station_ut`, `is_retrograde`, `get_station_type`, `swe_next_retrograde_ut` |
| Bodies | Mars(4), Jupiter(5), Saturn(6) |
| Date range | 2024-2030 |
| Tolerances | timing: 1s, type: identical, is_retrograde: identical |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestFindStation` | 3 planets x ~3 stations each | Station JD timing |
| `TestIsRetrograde` | 3 planets x 50 dates | Boolean match |
| `TestStationType` | 3 planets x 50 dates | SR/SD/None match |
| `TestNextRetrograde` | 3 planets x 3 start dates | Retrograde period start/end |

---

### 5.15 `test_compare_leb_houses.py` -- House Calculations

**Counterpart:** `tests/test_compare_houses.py`

**Purpose:** Validate that house calculations are unaffected by LEB mode for
non-Sunshine systems, and correct for Sunshine ('I') system which calls
`swe_calc_ut` for Sun declination.

| Parameter | Value |
|-----------|-------|
| Functions | `swe_houses`, `swe_houses_ex` |
| Systems | All 24 house systems |
| Locations | Rome, New York, Sydney, Equator |
| Dates | 3 test dates |
| Tolerances | Sunshine: 0.001", all others: exact match (0.0) |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestNonSunshineIdentical` | 23 systems x 4 locations x 3 dates | Cusps + ASCMC must be bit-for-bit identical |
| `TestSunshinePrecision` | 4 locations x 3 dates | Sunshine cusps within tolerance |
| `TestHousesExSidereal` | 4 sidereal modes x 4 locations | `swe_houses_ex` with sidereal flag |

**Key insight:** For all systems except Sunshine ('I'), `swe_houses` never calls
`swe_calc_ut`, so results MUST be identical regardless of LEB mode. This test
verifies that no unintended side effects leak through.

---

### 5.16 `test_compare_leb_gauquelin.py` -- Gauquelin Sectors

**Counterpart:** Part of `tests/test_compare_houses_ext.py`

**Purpose:** Validate `swe_gauquelin_sector` which calls `swe_calc_ut` for
planet positions.

| Parameter | Value |
|-----------|-------|
| Functions | `swe_gauquelin_sector` |
| Bodies | Sun(0), Moon(1), Mars(4), Jupiter(5) |
| Locations | Rome, New York, Sydney, Equator |
| Dates | 3 test dates |
| Methods | 0 (Placidus), 1 (Koch) |
| Tolerances | sector: 0.001 |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestGauquelinSector` | 4 bodies x 4 locations x 3 dates x 2 methods | Sector number comparison |

---

### 5.17 `test_compare_leb_rise_transit.py` -- Rise/Transit/Set

**Counterpart:** `tests/test_compare_rise_transit.py`

**Purpose:** Validate rise/transit/set timing. These functions use `swe_calc_ut`
in iterative position refinement loops.

| Parameter | Value |
|-----------|-------|
| Functions | `rise_trans` |
| Bodies | Sun(0), Moon(1), Venus(3) |
| Locations | Rome, New York, Sydney, Equator |
| RSMI types | 1 (rise), 2 (set), 4 (transit) |
| Dates | 3 test dates (equinox, solstice, random) |
| Tolerances | 1s timing |

**Test Classes:**

| Class | Tests | Description |
|-------|-------|-------------|
| `TestSunrise` | 4 locations x 3 dates | Rise timing |
| `TestSunset` | 4 locations x 3 dates | Set timing |
| `TestSunTransit` | 4 locations x 3 dates | Transit timing |
| `TestMoonrise` | 4 locations x 3 dates | Moon rise timing |
| `TestMoonset` | 4 locations x 3 dates | Moon set timing |
| `TestVenusRise` | 4 locations x 3 dates | Venus rise timing |

---

## 6. Excluded Areas

These areas have NO LEB impact and are correctly excluded from this suite:

| Area | Reason | Swiss counterpart |
|------|--------|-------------------|
| Coordinate transforms (`cotrans`, `azalt`) | Pure math, no ephemeris | `test_compare_coordinates.py` |
| Time functions (`julday`, `revjul`, `deltat`) | No ephemeris dependency | `test_compare_time.py` |
| Fixed stars (`swe_fixstar_ut`) | Skyfield-only, no LEB dispatch | `test_compare_fixedstars.py` |
| Phenomena (`swe_pheno_ut`) | Uses Skyfield directly, not `swe_calc_ut` | `test_compare_phenomena.py` |
| Heliacal events (`swe_heliacal_ut`) | Uses `swe_pheno_ut` (Skyfield-only) | `test_compare_heliacal.py` |
| Occultations (`swe_lun_occult_*`) | Not using `swe_calc_ut` directly | `test_compare_occultations.py` |
| Planet occultations | Not using `swe_calc_ut` directly | `test_compare_planet_occultations.py` |
| Planetary moons | Requires SPK kernels, not in LEB | `test_compare_planetary_moons.py` |
| Orbital elements (`swe_nod_aps_ut`) | No LEB dispatch | `test_compare_orbital.py` |
| Planetocentric (`swe_calc_pctr`) | No LEB dispatch | `test_compare_calc_pctr.py` |
| Utility functions (`degnorm`, `split_deg`) | Pure math | `test_compare_utilities.py` |
| Benchmark/performance | Different concern | `test_compare_benchmark.py` |

---

## 7. Execution and CI

### 7.1 Pytest Configuration

Register the `leb_compare` marker in `pyproject.toml`:

```toml
[tool.pytest.ini_options]
markers = [
    # ... existing markers ...
    "leb_compare: LEB vs Skyfield comparison tests",
]
```

### 7.2 Task Runner (poethepoet)

Add to `pyproject.toml`:

```toml
[tool.poe.tasks.test_leb_compare]
help = "Run LEB vs Skyfield comparison tests"
cmd = "pytest tests/test_leb/compare/ -v -m leb_compare"

[tool.poe.tasks.test_leb_compare_quick]
help = "Run LEB vs Skyfield comparison (core only, no slow)"
cmd = "pytest tests/test_leb/compare/ -v -m 'leb_compare and not slow' --timeout=120"
```

### 7.3 Execution Examples

```bash
# Full suite
poe test:leb:compare

# Quick (skip slow tests)
poe test:leb:compare:quick

# Single area
pytest tests/test_leb/compare/test_compare_leb_planets.py -v

# Single body
pytest tests/test_leb/compare/test_compare_leb_planets.py -v -k "Sun"

# With relaxed tolerances
LEB_TOL_POSITION_ARCSEC=0.5 pytest tests/test_leb/compare/ -v

# Verbose with statistics
pytest tests/test_leb/compare/ -v -s  # -s to see print output
```

### 7.4 Prerequisites

- `data/leb/ephemeris_medium.leb` must exist (pre-generated, checked in)
- `de440.bsp` ephemeris must be available (for Skyfield reference)
- No `pyswisseph` dependency (unlike the Swiss comparison suite)

---

## 8. Implementation Order

### Phase 1: Infrastructure + Core (files 0-4)

| Step | File | Est. LOC | Priority |
|------|------|----------|----------|
| 1.0 | `__init__.py` | 1 | P0 |
| 1.1 | `conftest.py` | ~150 | P0 |
| 1.2 | `test_compare_leb_planets.py` | ~250 | P0 |
| 1.3 | `test_compare_leb_lunar.py` | ~200 | P0 |
| 1.4 | `test_compare_leb_asteroids.py` | ~150 | P0 |
| 1.5 | `test_compare_leb_hypothetical.py` | ~150 | P0 |

**Validate:** Run Phase 1 tests, verify they pass, tune tolerances if needed.

### Phase 2: Flags + Precision (files 5-9)

| Step | File | Est. LOC | Priority |
|------|------|----------|----------|
| 2.1 | `test_compare_leb_observations.py` | ~300 | P0 |
| 2.2 | `test_compare_leb_sidereal.py` | ~250 | P0 |
| 2.3 | `test_compare_leb_ayanamsha.py` | ~150 | P1 |
| 2.4 | `test_compare_leb_velocities.py` | ~250 | P0 |
| 2.5 | `test_compare_leb_distances.py` | ~150 | P1 |

**Validate:** Run Phase 1+2, ensure all pass.

### Phase 3: Indirect Callers (files 10-17)

| Step | File | Est. LOC | Priority |
|------|------|----------|----------|
| 3.1 | `test_compare_leb_crossings.py` | ~300 | P0 |
| 3.2 | `test_compare_leb_eclipses_solar.py` | ~350 | P1 |
| 3.3 | `test_compare_leb_eclipses_lunar.py` | ~300 | P1 |
| 3.4 | `test_compare_leb_elongation.py` | ~150 | P1 |
| 3.5 | `test_compare_leb_stations.py` | ~250 | P1 |
| 3.6 | `test_compare_leb_houses.py` | ~250 | P1 |
| 3.7 | `test_compare_leb_gauquelin.py` | ~150 | P2 |
| 3.8 | `test_compare_leb_rise_transit.py` | ~200 | P1 |

**Validate:** Full suite green.

### Phase 4: Polish

| Step | Task | Priority |
|------|------|----------|
| 4.1 | Add `poe` task runner entries | P1 |
| 4.2 | Register `leb_compare` pytest marker | P1 |
| 4.3 | Add environment variable tolerance overrides | P2 |
| 4.4 | Review and finalize tolerances based on actual results | P0 |

### Estimated Total

| Metric | Value |
|--------|-------|
| Files | 18 (1 init + 1 conftest + 17 test files) |
| Estimated LOC | ~3,800 |
| Test cases | ~600+ (parametrized) |
| Bodies covered | 30/30 LEB bodies |
| Flag combinations | 8+ native + 4 fallback |
| Functions covered | 62 (34 direct + 28 indirect) |

---

## 9. File Inventory

| # | File | Type | LOC (est.) |
|---|------|------|------------|
| 0 | `tests/test_leb/compare/__init__.py` | Package | 1 |
| 1 | `tests/test_leb/compare/conftest.py` | Infrastructure | 150 |
| 2 | `tests/test_leb/compare/test_compare_leb_planets.py` | Core | 250 |
| 3 | `tests/test_leb/compare/test_compare_leb_lunar.py` | Core | 200 |
| 4 | `tests/test_leb/compare/test_compare_leb_asteroids.py` | Core | 150 |
| 5 | `tests/test_leb/compare/test_compare_leb_hypothetical.py` | Core | 150 |
| 6 | `tests/test_leb/compare/test_compare_leb_observations.py` | Flags | 300 |
| 7 | `tests/test_leb/compare/test_compare_leb_sidereal.py` | Flags | 250 |
| 8 | `tests/test_leb/compare/test_compare_leb_ayanamsha.py` | Precision | 150 |
| 9 | `tests/test_leb/compare/test_compare_leb_velocities.py` | Precision | 250 |
| 10 | `tests/test_leb/compare/test_compare_leb_distances.py` | Precision | 150 |
| 11 | `tests/test_leb/compare/test_compare_leb_crossings.py` | Indirect | 300 |
| 12 | `tests/test_leb/compare/test_compare_leb_eclipses_solar.py` | Indirect | 350 |
| 13 | `tests/test_leb/compare/test_compare_leb_eclipses_lunar.py` | Indirect | 300 |
| 14 | `tests/test_leb/compare/test_compare_leb_elongation.py` | Indirect | 150 |
| 15 | `tests/test_leb/compare/test_compare_leb_stations.py` | Indirect | 250 |
| 16 | `tests/test_leb/compare/test_compare_leb_houses.py` | Indirect | 250 |
| 17 | `tests/test_leb/compare/test_compare_leb_gauquelin.py` | Indirect | 150 |
| 18 | `tests/test_leb/compare/test_compare_leb_rise_transit.py` | Indirect | 200 |
| | **Total** | | **~3,800** |
