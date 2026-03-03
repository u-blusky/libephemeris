## Project Overview

LibEphemeris is a high-precision astronomical ephemeris library for Python, providing Swiss Ephemeris-compatible API using NASA JPL DE440 ephemeris via Skyfield. 

VERY IMPORTANT:
*Keep always 1:1 compatibility with PySwissEphemeris.*

## Commands

Uses `uv` for dependencies and `poe` (poethepoet) for task running.

```bash
uv pip install -e ".[dev]"        # Install with dev dependencies
poe lint                          # Ruff linter with auto-fix
poe format                        # Ruff formatter
poe typecheck                     # mypy
```

### Important: Never Run the Full Test Suite

**Do NOT run `poe test` or `poe test:full` during development.** Always run targeted tests:

```bash
pytest tests/test_file.py -v                          # Single file
pytest tests/test_file.py::TestClass::test_method -v  # Single test
pytest -k "pattern" -v                                # By keyword

# Feature-specific suites
poe test:lunar              # All lunar tests (no slow)
poe test:lunar:perigee      # Perigee tests only
poe test:lunar:apogee       # Apogee tests only
poe test:lunar:lilith       # Lilith tests only
poe test:leb                # LEB unit tests (no slow)
```

## Code Style

- `from __future__ import annotations` at top of every module
- Imports grouped: stdlib, third-party, local (relative imports)
- Line length 88, Python 3.9+, double quotes, Ruff formatter
- Google-style docstrings with Args/Returns/Raises
- Naming: `snake_case` functions, `PascalCase` classes, `SCREAMING_SNAKE_CASE` constants, `_underscore` private
- Swiss Ephemeris compatible functions use `swe_` prefix
- Always return native Python floats (not numpy types)
- Exceptions in `libephemeris/exceptions.py`: `Error`, `CoordinateError`, `UnknownBodyError`, `EphemerisRangeError`, `PolarCircleError`

## Ephemeris File Selection

Three precision tiers: `base` (de440s.bsp, 1849-2150), `medium` (de440.bsp, 1550-2650, **default**), `extended` (de441.bsp, -13200 to +17191).

`SEFLG_MOSEPH` is accepted for API compatibility but silently ignored. All calculations always use JPL DE440/DE441 via Skyfield.

## Binary Ephemeris Mode (LEB)

Precomputed `.leb` files with Chebyshev polynomial approximations (~14x speedup). Automatic fallback to Skyfield for unsupported bodies/flags (`SEFLG_TOPOCTR`, `SEFLG_XYZ`, `SEFLG_RADIANS`, `SEFLG_NONUT`).

```python
from libephemeris import set_leb_file
set_leb_file("/path/to/ephemeris.leb")  # or env LIBEPHEMERIS_LEB=...
```

Group workflow is **recommended** for generation (avoids macOS multiprocessing deadlocks):

```bash
poe leb:generate:base:groups      # All 3 groups + merge for base tier
poe leb:generate:medium:groups    # All 3 groups + merge for medium tier
poe leb:generate:extended:groups  # All 3 groups + merge for extended tier
```

Runtime always uses a **single merged file**. See `docs/leb/design.md` and `docs/leb/guide.md` for details.

## Lunar Calibration Workflow

1. `poe calibrate-perigee` (or `poe calibrate-perigee:quick`)
2. Paste coefficients into `_calc_elp2000_perigee_perturbations()` in `lunar.py`
3. `poe generate-lunar-corrections` (regenerates `lunar_corrections.py`)
4. `poe test:lunar:perigee`

See `docs/methodology/interpolated-perigee.md` for the full methodology.
