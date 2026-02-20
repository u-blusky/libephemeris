# AGENTS.md - Guide for AI Coding Agents

This document provides essential information for AI agents working on the libephemeris codebase.

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
└── time_utils.py         # Julian day conversions
tests/
├── conftest.py           # Shared fixtures and markers
└── test_*/               # Test subdirectories by feature
```

## Scripts & Calibration

The `scripts/` directory contains data generation and calibration tools:

```
scripts/
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
3. **Sync**: update the matching series in `scripts/generate_lunar_corrections.py`
4. **Regenerate**: `poe generate-lunar-corrections` (regenerates `lunar_corrections.py`)
5. **Test**: `poe test:lunar:perigee`

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

The `LIBEPHEMERIS_EPHEMERIS` environment variable controls which JPL ephemeris file is used. For example, set `LIBEPHEMERIS_EPHEMERIS=de441.bsp` to use DE441 for extended date range coverage. The default is DE440.

*KEEP ALWAYS 1:1 Compatibility with PySwissEphemeris!*
