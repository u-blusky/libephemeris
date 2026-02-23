# Migration Guide: pyswisseph to libephemeris

This guide helps users migrate from `pyswisseph` (Python bindings to Swiss Ephemeris) to `libephemeris`.

## Overview

`libephemeris` is designed as a **drop-in replacement** for `pyswisseph` in most common use cases. It provides:

- Same function names (with both `swe_` prefixed and unprefixed aliases)
- Same flag constants (`SEFLG_*`, `SE_SIDM_*`, `SE_*`)
- Same return value structures
- Same behavior for global state management

However, there are important differences to be aware of when migrating.

---

## Quick Migration

### Minimal Change

In most cases, you can simply replace your import:

```python
# Before (pyswisseph)
import swisseph as swe

# After (libephemeris)
import libephemeris as swe
```

Your existing code should continue to work for common planetary calculations, house systems, and ayanamshas.

### Constants Import

```python
# Before (pyswisseph)
import swisseph as swe
planet = swe.SUN
flag = swe.FLG_SPEED

# After (libephemeris) - Option 1: Same style
import libephemeris as swe
planet = swe.SE_SUN   # Note: SE_ prefix is used
flag = swe.SEFLG_SPEED

# After (libephemeris) - Option 2: Explicit constants import (recommended)
import libephemeris as swe
from libephemeris.constants import SE_SUN, SEFLG_SPEED
```

---

## API Differences

### Function Names

Most functions are available with **both** the `swe_` prefix and without:

| pyswisseph | libephemeris |
|------------|--------------|
| `swe.calc_ut()` | `swe.calc_ut()` or `swe.swe_calc_ut()` |
| `swe.houses()` | `swe.houses()` or `swe.swe_houses()` |
| `swe.julday()` | `swe.julday()` or `swe.swe_julday()` |
| `swe.set_sid_mode()` | `swe.set_sid_mode()` or `swe.swe_set_sid_mode()` |
| `swe.get_ayanamsa_ut()` | `swe.get_ayanamsa_ut()` or `swe.swe_get_ayanamsa_ut()` |

### Constant Names

libephemeris uses `SE_` and `SEFLG_` prefixes consistently:

| pyswisseph | libephemeris |
|------------|--------------|
| `swe.SUN` | `SE_SUN` (0) |
| `swe.MOON` | `SE_MOON` (1) |
| `swe.FLG_SPEED` | `SEFLG_SPEED` |
| `swe.FLG_SIDEREAL` | `SEFLG_SIDEREAL` |
| `swe.FLG_TOPOCTR` | `SEFLG_TOPOCTR` |
| `swe.SIDM_LAHIRI` | `SE_SIDM_LAHIRI` (1) |

### House Cusp Array Indexing

Both libraries return 12 house cusps, but indexing may differ:

```python
# pyswisseph returns 13 elements (index 0 is unused, cusps are 1-12)
cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, b'P')
cusp1_swe = cusps_swe[1]

# libephemeris returns 12 elements (0-indexed)
cusps, ascmc = ephem.swe_houses(jd, lat, lon, b'P')
cusp1 = cusps[0]  # First house cusp
```

---

## Precision Differences

libephemeris uses NASA JPL DE ephemerides (via Skyfield) instead of Swiss Ephemeris data files. Here are the validated precision differences:

### Planetary Positions

| Component | Max Difference | Notes |
|-----------|---------------|-------|
| Longitude (geocentric) | < 0.001 degrees | Sub-arcsecond accuracy |
| Latitude (geocentric) | < 0.001 degrees | |
| Distance | < 0.0001 AU | |
| Heliocentric longitude | < 0.03 degrees | Slightly relaxed tolerance |

### House Cusps and Angles

| Component | Max Difference | Notes |
|-----------|---------------|-------|
| House cusps | < 0.001 degrees | All 19 house systems |
| Ascendant | < 0.001 degrees | |
| MC | < 0.001 degrees | |
| ARMC | < 0.001 degrees | |
| Vertex | < 0.01 degrees | |

### Ayanamshas

| Component | Max Difference | Notes |
|-----------|---------------|-------|
| Standard ayanamshas | < 0.01 degrees | Fagan-Bradley, Lahiri, etc. |
| Star-based ayanamshas | < 0.06 degrees | True Citra, Galactic Center, etc. |

Star-based ayanamshas have slightly higher tolerance due to differences in star position calculations, precession models, and coordinate transformations between the underlying ephemeris engines.

### Velocities

| Component | Max Difference | Notes |
|-----------|---------------|-------|
| Angular velocity | < 0.01 degrees/day | |
| Radial velocity | < 0.001 AU/day | |

### Lunar Nodes

| Component | Max Difference | Notes |
|-----------|---------------|-------|
| Mean Node | < 0.005 degrees | High precision |
| True Node | < 0.14 degrees | Different oscillation model |

The True Node (osculating node) shows larger differences due to different algorithms for computing the short-period oscillations. Both implementations agree on the mean position but differ on the instantaneous oscillation.

### Lilith (Lunar Apogee)

| Component | Max Difference | Notes |
|-----------|---------------|-------|
| Mean Apogee (Mean Lilith) | < 0.01 degrees | High precision |
| Osculating Apogee (True Lilith) | ~0.015° mean, ~0.065° max | Sub-arcminute precision |

**Note**: True Lilith (SE_OSCU_APOG, body ID 13) now achieves excellent precision (~0.015° mean difference from pyswisseph) through calibrated perturbation corrections applied to osculating orbital elements derived from JPL DE440 state vectors. See `TRUE_LILITH_METHODS.md` for details.

---

## Features Not Yet Implemented

The following features are present in pyswisseph but **not yet fully implemented** in libephemeris:

### Eclipse Functions (Partial)

Eclipse functions are implemented but some return values are not yet calculated:

- **Saros series number**: Returns 0 (not implemented)
- **Inex number**: Returns 0 (not implemented)
- **Sunrise/sunset on central line**: Returns 0 for solar eclipses (not implemented)

Affected functions:
- `sol_eclipse_when_glob()` / `swe_sol_eclipse_when_glob()`
- `sol_eclipse_when_loc()` / `swe_sol_eclipse_when_loc()`
- `lun_eclipse_when()` / `swe_lun_eclipse_when()`
- `lun_occult_when_glob()` / `swe_lun_occult_when_glob()`

### Fixed Star Velocities

Fixed star velocity calculations return 0:

```python
pos, _ = ephem.fixstar_ut("Aldebaran", jd, SEFLG_SPEED)
# pos[3], pos[4], pos[5] are 0.0 (velocities not implemented)
```

### Date Range Limitations

libephemeris uses JPL DE ephemerides with specific date ranges:

| Ephemeris | Date Range | Size | Notes |
|-----------|-----------|------|-------|
| DE440s | 1849-2150 | ~31 MB | Lightweight subset of DE440 |
| DE440 (default) | 1550-2650 | ~128 MB | ICRF 3.0, recommended |
| DE441 | -13200 to 17191 | ~3.4 GB | Extended version of DE440 |

### Precision Tiers

libephemeris organises the current-generation files into three precision tiers:

| Tier | File | Use Case |
|------|------|----------|
| `base` | de440s.bsp | Lightweight, modern-era usage |
| `medium` | de440.bsp | General purpose **(DEFAULT)** |
| `extended` | de441.bsp | Historical/far-future research |

```python
from libephemeris import set_precision_tier

set_precision_tier("extended")   # uses de441.bsp
```

```bash
export LIBEPHEMERIS_PRECISION=extended
```

### Selecting an Ephemeris File Directly

You can also bypass the tier system and select a specific file:

```python
from libephemeris import set_ephemeris_file, set_ephe_path

# Set custom ephemeris file
set_ephemeris_file("de441.bsp")  # Extended range version of DE440

# Or specify a custom directory
set_ephe_path("/path/to/ephemeris/files")
```

You can also select the ephemeris file via the `LIBEPHEMERIS_EPHEMERIS` environment variable:

```bash
export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

### Resolution Priority

Resolution priority (highest to lowest):

1. `LIBEPHEMERIS_EPHEMERIS` environment variable
2. `set_ephemeris_file()` / `set_jpl_file()` programmatic call
3. Precision tier (`LIBEPHEMERIS_PRECISION` env var or `set_precision_tier()`)
4. Default: `de440.bsp` (medium tier)

---

## Thread Safety with EphemerisContext

**This is a major difference from pyswisseph.**

The Swiss Ephemeris (and pyswisseph) uses global state and is **NOT thread-safe**. libephemeris provides the same behavior for the module-level API, but also offers a **thread-safe alternative** via `EphemerisContext`.

### Global API (pyswisseph-compatible, NOT thread-safe)

```python
import libephemeris as swe

# This works like pyswisseph - uses global state
swe.set_topo(12.5, 41.9, 0)  # Rome
swe.set_sid_mode(1)  # Lahiri

pos, _ = swe.calc_ut(2451545.0, swe.SE_SUN, swe.SEFLG_SIDEREAL)
```

**Warning**: Do NOT use the global API in multi-threaded applications. Different threads modifying global state (topo, sidereal mode) will cause race conditions.

### Context API (thread-safe)

For multi-threaded applications, use `EphemerisContext`:

```python
from libephemeris import EphemerisContext, SE_SUN, SE_MOON, SEFLG_SIDEREAL

# Each thread creates its own context with isolated state
ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)    # Rome - isolated to this context
ctx.set_sid_mode(1)            # Lahiri - isolated to this context

# Thread-safe calculations
sun, _ = ctx.calc_ut(2451545.0, SE_SUN, SEFLG_SIDEREAL)
moon, _ = ctx.calc_ut(2451545.0, SE_MOON, SEFLG_SIDEREAL)
cusps, ascmc = ctx.houses(2451545.0, 41.9, 12.5, ord('P'))
```

### Multi-Threading Example

```python
from libephemeris import EphemerisContext, SE_SUN, SE_MOON, SEFLG_SIDEREAL
from concurrent.futures import ThreadPoolExecutor

def calculate_chart(location: dict, jd: float) -> dict:
    """Thread-safe chart calculation function."""
    # Each thread creates its own context
    ctx = EphemerisContext()
    ctx.set_topo(location['lon'], location['lat'], 0)
    ctx.set_sid_mode(1)  # Lahiri

    sun, _ = ctx.calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
    moon, _ = ctx.calc_ut(jd, SE_MOON, SEFLG_SIDEREAL)
    cusps, ascmc = ctx.houses(jd, location['lat'], location['lon'], ord('P'))

    return {
        'location': location['name'],
        'sun': sun[0],
        'moon': moon[0],
        'asc': ascmc[0]
    }

# Different locations
locations = [
    {'name': 'Rome', 'lon': 12.5, 'lat': 41.9},
    {'name': 'London', 'lon': -0.1, 'lat': 51.5},
    {'name': 'Tokyo', 'lon': 139.7, 'lat': 35.7},
    {'name': 'New York', 'lon': -74.0, 'lat': 40.7},
    {'name': 'Sydney', 'lon': 151.2, 'lat': -33.9},
]

# Calculate concurrently - each thread has its own isolated context
jd = 2451545.0  # J2000.0
with ThreadPoolExecutor(max_workers=5) as executor:
    results = list(executor.map(
        lambda loc: calculate_chart(loc, jd),
        locations
    ))

for r in results:
    print(f"{r['location']}: Sun={r['sun']:.2f}, Moon={r['moon']:.2f}, Asc={r['asc']:.2f}")
```

### EphemerisContext Methods

| Method | Description |
|--------|-------------|
| `__init__(ephe_path, ephe_file)` | Create context with optional ephemeris path/file |
| `set_topo(lon, lat, alt)` | Set observer location (isolated to this context) |
| `get_topo()` | Get current observer location |
| `set_sid_mode(mode, t0, ayan_t0)` | Set sidereal mode (isolated to this context) |
| `get_sid_mode(full=False)` | Get current sidereal mode |
| `calc_ut(tjd_ut, ipl, iflag)` | Calculate planetary position (UT) |
| `calc(tjd, ipl, iflag)` | Calculate planetary position (TT/ET) |
| `calc_pctr(tjd_ut, ipl, iplctr, iflag)` | Planet-centric position |
| `houses(tjd_ut, lat, lon, hsys)` | Calculate house cusps and angles |
| `close()` | Class method - close all shared ephemeris resources |

### Resource Sharing

`EphemerisContext` is designed for both thread safety and memory efficiency:

- **Isolated State**: Each context has its own observer location, sidereal mode, and angles cache
- **Shared Resources**: Expensive resources (ephemeris files ~16MB+, timescale data) are shared across all contexts
- **Thread-Safe Loading**: Uses double-checked locking pattern for lazy initialization

```python
# Multiple contexts share the same ephemeris file (memory efficient)
ctx1 = EphemerisContext()
ctx2 = EphemerisContext()
ctx3 = EphemerisContext()

# Each context has isolated state
ctx1.set_sid_mode(1)  # Lahiri
ctx2.set_sid_mode(0)  # Fagan-Bradley
ctx3.set_sid_mode(27) # True Citra

# But they all share the same ephemeris data in memory
```

---

## SEFLG_MOSEPH (Moshier Ephemeris Flag)

The `SEFLG_MOSEPH` flag is accepted for API compatibility but **silently ignored**. All calculations in libephemeris always use JPL DE440/DE441 via Skyfield, regardless of whether `SEFLG_MOSEPH` is passed. Code that previously used `SEFLG_MOSEPH` to select the Moshier semi-analytical ephemeris will continue to work without errors, but will use the JPL ephemeris instead.

```python
# This still works, but SEFLG_MOSEPH is silently ignored:
pos, _ = swe.calc_ut(jd, swe.SE_SUN, swe.SEFLG_MOSEPH | swe.SEFLG_SPEED)
# Equivalent to:
pos, _ = swe.calc_ut(jd, swe.SE_SUN, swe.SEFLG_SPEED)
```

---

## Migration Checklist

- [ ] Replace `import swisseph as swe` with `import libephemeris as swe`
- [ ] Update constant names if using unprefixed versions (`SUN` -> `SE_SUN`)
- [ ] Check house cusp array indexing (0-based in libephemeris)
- [ ] Verify date range is within ephemeris coverage (1550-2650 for DE440)
- [ ] For multi-threaded apps: migrate to `EphemerisContext` API
- [ ] Update tests for relaxed tolerances on star-based ayanamshas (< 0.06 degrees)
- [ ] Handle eclipse functions that return 0 for Saros/Inex numbers
- [ ] Review True Node usage (up to 0.14 degrees difference from pyswisseph)
- [ ] Review True Lilith usage (~0.065° max difference - sub-arcminute precision)

---

## Reporting Issues

If you encounter compatibility issues not covered in this guide, please report them at:

https://github.com/g-battaglia/libephemeris/issues

Include:
1. Your pyswisseph code that doesn't work
2. The expected result from pyswisseph
3. The actual result from libephemeris
4. Python version and libephemeris version
