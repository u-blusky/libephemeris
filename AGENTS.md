## Project Overview

LibEphemeris is a high-precision astronomical ephemeris library for Python, providing Swiss Ephemeris-compatible API using NASA JPL DE440 ephemeris via Skyfield. The library focuses on scientific accuracy, maintainability, and clean code.

## Build, Lint, and Test Commands

### Package Manager
Uses `uv` for dependency management and `poe` (poethepoet) for task running.

### Installation
```bash
uv pip install -e ".[dev]"        # Install with dev dependencies
```

### Testing Commands
```bash
poe test                          # Fast tests only (excludes @pytest.mark.slow)
poe test:full                     # All tests including slow ones

# Run a single test file
pytest tests/test_file.py -v

# Run a single test by name
pytest tests/test_file.py::TestClass::test_method -v
pytest -k "test_name_pattern" -v

# Run tests by marker
pytest -m "unit"                  # Unit tests only
pytest -m "not slow"              # Exclude slow tests

# Coverage
poe coverage                      # With coverage report
```

### Linting and Formatting
```bash
poe lint                          # Run ruff linter with auto-fix
poe format                        # Run ruff formatter
poe typecheck                     # Run mypy type checking
```

### Test Markers (from pytest.ini)
- `slow` - Long-running tests (deselect with `-m "not slow"`)
- `network` - Requires network access
- `unit` / `integration` - Test categories
- `precision` / `comparison` / `edge_case` - Specialized tests

### Important: Never Run the Full Test Suite
**Do NOT run `poe test` or `poe test:full` during development.** The full suite is
too large and slow. Always run targeted tests for the area you changed:

```bash
# Run a single test file
pytest tests/test_lunar/test_elp2000_perigee_perturbations.py -v

# Run a specific test class
pytest tests/test_lunar/test_file.py::TestClassName -v

# Run by marker or keyword
pytest -k "perigee" -v
pytest -m "unit" tests/test_lunar/ -v

# Feature-specific suites (use these instead of poe test)
poe test:lunar              # All lunar tests (no slow)
poe test:lunar:perigee      # Perigee tests only
poe test:lunar:apogee       # Apogee tests only
poe test:lunar:lilith       # Lilith tests only
```

## Code Style Guidelines

### Imports
```python
from __future__ import annotations      # Always first

import math                              # Standard library
from typing import Tuple, TYPE_CHECKING

from skyfield.api import Star            # Third-party

from .constants import SE_SUN, SEFLG_SPEED  # Relative internal imports
```

- Use `from __future__ import annotations` at top of every module
- Group imports: stdlib, third-party, local (separated by blank lines)
- Star imports allowed in `__init__.py` and `constants.py`

### Formatting
- **Line length**: 88 characters (Ruff/Black default)
- **Python version**: 3.9+
- **Formatter**: Ruff (preferred)
- **Quotes**: Double quotes

### Naming Conventions
- **Functions**: `snake_case` - `calc_mean_lunar_node()`, `swe_calc_ut()`
- **Swiss Ephemeris functions**: `swe_` prefix - `swe_calc_ut()`, `swe_houses()`
- **Classes**: `PascalCase` - `EphemerisContext`, `BesselianElements`
- **Constants**: `SCREAMING_SNAKE_CASE` - `SE_SUN`, `SEFLG_SPEED`
- **Private**: Single underscore prefix - `_PLANET_MAP`, `_calc_body()`

### Docstrings (Google-style)
```python
def swe_calc_ut(tjd_ut: float, ipl: int, iflag: int) -> Tuple[...]:
    """Calculate planetary position for Universal Time.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (SE_SUN, SE_MOON, etc.)
        iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.)

    Returns:
        Tuple of (position_tuple, return_flag)

    Raises:
        EphemerisRangeError: If date is outside ephemeris coverage
    """
```

### Error Handling
Use exceptions from `libephemeris/exceptions.py`:
```python
from .exceptions import (
    Error,                    # Base class (pyswisseph compatible)
    CoordinateError,          # Invalid lat/lon
    UnknownBodyError,         # Unknown celestial body
    EphemerisRangeError,      # Date outside range
    PolarCircleError,         # House calc at polar latitudes
)

raise CoordinateError(
    message=f"latitude {lat} is out of valid range",
    coordinate_name="latitude",
    value=lat,
    min_value=-90.0,
    max_value=90.0,
)
```

### Test Structure
```python
import pytest
import libephemeris as ephem
from libephemeris.constants import *

class TestHouseAlgorithms:
    @pytest.mark.unit
    def test_equal_house_formula(self):
        """Test Equal house formula."""
        jd = 2451545.0  # J2000
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("E"))
        assert abs(cusps[0] - ascmc[0]) < 0.001
```

Use fixtures from `tests/conftest.py`: `standard_jd`, `test_locations`, `all_planets`, `default_tolerances`

## Project Structure

```
libephemeris/
├── __init__.py           # Public API exports
├── constants.py          # All SE_* and SEFLG_* constants
├── planets.py            # Core planetary calculations
├── houses.py             # House system calculations
├── lunar.py              # Lunar nodes and Lilith
├── minor_bodies.py       # Asteroids and TNOs
├── eclipse.py            # Eclipse calculations
├── exceptions.py         # Exception hierarchy
├── context.py            # Thread-safe EphemerisContext
├── astrometry.py         # IAU precession/nutation/aberration utilities
├── state.py              # Global state management
├── time_utils.py         # Julian day conversions
├── leb_format.py         # LEB binary format constants, dataclasses, struct helpers
├── leb_reader.py         # LEB mmap reader + Clenshaw polynomial evaluation
└── fast_calc.py          # LEB calculation pipelines (ICRS, ecliptic, heliocentric)
data/
└── leb/                  # Precomputed LEB binary ephemeris files
    └── ephemeris_*.leb   # Per-tier merged files (base, medium, extended)
scripts/
└── generate_leb.py       # LEB binary ephemeris generator (group/merge workflow)
tests/
├── conftest.py           # Shared fixtures and markers
└── test_*/               # Test subdirectories by feature (includes test_leb/)
```

## Scripts & Calibration

The `scripts/` directory contains data generation and calibration tools:

```
scripts/
├── generate_leb.py                     # LEB binary ephemeris generator (group/merge workflow)
├── calibrate_perigee_perturbations.py  # Perigee perturbation series calibration (v2.2)
├── generate_lunar_corrections.py       # Lunar correction table generation
├── generate_planet_centers_spk.py      # Planet centers SPK file generation
├── build_star_catalog.py               # Fixed star catalog builder
├── calibrate_true_node.py              # True node calibration
├── download_spk.py                     # SPK kernel downloader
├── download_max_range_spk.py           # Extended-range SPK downloader
├── update_orbital_elements.py          # Orbital elements updater
├── tier_diagnostic_*.py                # Tier diagnostic tables
├── profile_hot_paths.py                # Performance profiling
└── release_planet_centers.py           # Planet centers release script
```

### Lunar Calibration Workflow

When modifying the perigee perturbation series, follow this order:

1. **Calibrate**: `poe calibrate-perigee` (or `poe calibrate-perigee:quick` for validation)
2. **Apply**: paste new coefficients into `_calc_elp2000_perigee_perturbations()` in `lunar.py`
3. **Regenerate**: `poe generate-lunar-corrections` (regenerates `lunar_corrections.py`)
4. **Test**: `poe test:lunar:perigee`

Note: `generate_lunar_corrections.py` imports the perturbation series directly from
`lunar.py`, so there is no manual sync step needed.

See `docs/interpolated_perigee_methodology.md` for the full methodology.

### Lunar Test Commands
```bash
poe test:lunar              # All lunar tests (no slow)
poe test:lunar:perigee      # Perigee: perturbations + interpolated + osculating
poe test:lunar:apogee       # Apogee: perturbations + interpolated
poe test:lunar:lilith       # Lilith: mean + true (all correction modules)
```

## Key Patterns

### Swiss Ephemeris API Compatibility
```python
# Both work identically:
ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
ephem.calc_ut(jd, SE_SUN, SEFLG_SPEED)
```


### Adding a New Function
1. Implement in appropriate module (`planets.py`, `lunar.py`, etc.)
2. Add Google-style docstring with Args, Returns, Raises
3. Export from `__init__.py` with both `swe_` prefix and alias
4. Write tests in `tests/test_*` directory

### Return Native Python Floats
```python
# Always convert numpy types to native Python floats
return _to_native_floats((lon, lat, dist, dlon, dlat, ddist)), iflag
```

## Important Notes

> `SEFLG_MOSEPH` is accepted for API compatibility but silently ignored. All calculations always use JPL DE440/DE441 via Skyfield.

### Ephemeris File Selection

LibEphemeris supports three precision tiers, each mapping to a different JPL ephemeris:

| Tier | File | Date Range | Size | Use Case |
|------|------|------------|------|----------|
| `base` | de440s.bsp | 1849-2150 | ~31 MB | Lightweight, modern-era usage |
| `medium` | de440.bsp | 1550-2650 | ~128 MB | General purpose **(DEFAULT)** |
| `extended` | de441.bsp | -13200 to +17191 | ~3.1 GB | Historical/far-future research |

**Resolution priority** (highest to lowest):
1. `LIBEPHEMERIS_EPHEMERIS` environment variable
2. `set_ephemeris_file()` / `set_jpl_file()` programmatic call
3. Precision tier (`LIBEPHEMERIS_PRECISION` env var or `set_precision_tier()`)
4. Default: `de440.bsp` (medium tier)

```bash
# Select by tier
export LIBEPHEMERIS_PRECISION=extended   # uses de441.bsp

# Or select ephemeris file directly
export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

```python
# Programmatic equivalents
from libephemeris import set_precision_tier, set_ephemeris_file

set_precision_tier("extended")       # uses de441.bsp
# or
set_ephemeris_file("de441.bsp")      # direct override (takes priority over tier)
```

*KEEP ALWAYS 1:1 Compatibility with PySwissEphemeris!*

### Binary Ephemeris Mode (LEB)

LibEphemeris supports a precomputed binary ephemeris mode (`.leb` files) that
stores Chebyshev polynomial approximations of planetary positions. This provides
~14x speedup over the default Skyfield pipeline.

**Architecture:**
- `libephemeris/leb_format.py` — Binary format constants, dataclasses, struct helpers
- `libephemeris/leb_reader.py` — mmap reader + Clenshaw evaluation (pure Python)
- `libephemeris/fast_calc.py` — Three calculation pipelines (ICRS, ecliptic-direct, heliocentric)
- `scripts/generate_leb.py` — CLI generator for `.leb` files

**Two modes (transparent to callers):**
1. **Skyfield mode (default):** Full Skyfield/JPL pipeline, unchanged
2. **Binary mode:** Read precomputed `.leb` files, automatic silent fallback to Skyfield for unsupported bodies/flags

**Activation:**
```python
from libephemeris import set_leb_file
set_leb_file("/path/to/ephemeris.leb")  # enables binary mode
set_leb_file(None)                       # disables binary mode
```
Or via environment variable: `LIBEPHEMERIS_LEB=/path/to/ephemeris.leb`

**Flags that trigger Skyfield fallback:**
`SEFLG_TOPOCTR`, `SEFLG_XYZ`, `SEFLG_RADIANS`, `SEFLG_NONUT`

**LEB Test Commands:**
```bash
poe test:leb                  # LEB unit tests (excludes slow)
poe test:leb:precision        # Full precision test suite (slow)
poe test:leb:precision:quick  # Precision tests for medium tier only
```

**LEB Generation Commands:**

The generator supports two workflows: full (all bodies at once) and group-based
(generate each body group independently, then merge). The group workflow is
**recommended** — it avoids macOS multiprocessing deadlocks and allows
regenerating a single group without redoing everything.

```bash
# Full generation (all bodies at once)
poe leb:generate:base         # Generate base tier LEB
poe leb:generate:medium       # Generate medium tier LEB
poe leb:generate:extended     # Generate extended tier LEB
poe leb:generate:all          # Generate all three tiers

# Group generation + merge (RECOMMENDED)
poe leb:generate:base:groups      # All 3 groups + merge for base tier
poe leb:generate:medium:groups    # All 3 groups + merge for medium tier
poe leb:generate:extended:groups  # All 3 groups + merge for extended tier

# Individual groups (for targeted regeneration)
poe leb:generate:base:planets     # Planets only (Sun-Pluto, Earth)
poe leb:generate:base:asteroids   # Asteroids only (Chiron, Ceres-Vesta)
poe leb:generate:base:analytical  # Analytical bodies (nodes, Lilith, Uranians)
poe leb:generate:base:merge       # Merge partial files into final .leb
```

**Body groups:**
- `planets` — Sun, Moon, Mercury-Pluto, Earth (body IDs 0-9, 14)
- `asteroids` — Chiron, Ceres, Pallas, Juno, Vesta (body IDs 15, 17-20)
- `analytical` — Mean/True/Osculating Node, Mean/True Lilith, 8 Uranians (body IDs 10-13, 16, 40-47)

**Group workflow details:**
- Each `--group` run produces a partial `.leb` file (e.g. `ephemeris_base_planets.leb`)
- `--merge` combines partial files into one complete file (zero re-computation)
- Partial files are generation-time intermediates only; **runtime always uses a single merged file**
- No auto-discovery or multi-file reader at runtime

**Key implementation details:**
- Zero new runtime dependencies (only `mmap`, `struct`, `math`, `dataclasses`)
- `numpy` only needed by the generator script (already a dev dependency)
- Thread-safe: sidereal params passed explicitly through call chain
- `EphemerisContext` passes sidereal config via keyword-only args (no global swap)
- `LEBReader.__init__` uses try/except around parse to prevent resource leaks
- Pipeline A defers Earth position fetch for `SEFLG_HELCTR`/`SEFLG_BARYCTR`
- Sidereal speed correction subtracts IAU 2006 general precession rate from `dlon`
- Analytical bodies run sequentially in main process (no `ProcessPoolExecutor` — removed due to macOS deadlocks with numpy/BLAS/Accelerate)
- Asteroid generation uses `spktype21` for SPK type 21 files (~36x faster than scalar fallback)

**Reference:** See `docs/LEB_PLAN.md` for the implementation plan and `docs/LEB_GUIDE.md` for the comprehensive technical guide.
