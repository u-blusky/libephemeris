## Project Overview

LibEphemeris is a pure-Python astronomical ephemeris library with a Swiss Ephemeris-compatible API, powered by NASA JPL DE440/DE441.

**Always maintain 1:1 compatibility with PySwissEphemeris.**

## Commands

```bash
uv pip install -e ".[dev]"              # Install with dev deps

leph code lint                          # Ruff linter with auto-fix
leph code format                        # Ruff formatter
leph code typecheck                     # mypy

leph test skyfield essential            # ~490 tests, ~20s (fast sanity check)
leph test leb-backend unit-fast         # ~5890 tests, ~1 min [RECOMMENDED]
pytest tests/test_file.py -v            # Single file
pytest tests/test_file.py::test_name -v # Single test
```

**Never run the full test suite.** Always use `leph test <subgroup>` or `pytest` on specific files. Full CLI reference (dev, prod, poe shortcuts, TAB completion setup): see `CLI.md`.

## Code Style

- `from __future__ import annotations` at top of every module
- Line length 88, Python 3.9+, double quotes, Ruff formatter
- Google-style docstrings, `snake_case` functions, `PascalCase` classes, `_underscore` private
- `swe_` prefix for Swiss Ephemeris-compatible functions
- Always return native Python floats (not numpy types)

## Architecture

Three tiers: `base` (de440s, 1849-2150), `medium` (de440, 1550-2650, **default**), `extended` (de441, -13200 to +17191).

Four calc modes (`set_calc_mode()` / `LIBEPHEMERIS_MODE`): `auto` (default: LEB->Horizons->Skyfield), `leb`, `horizons`, `skyfield`.

LEB = precomputed Chebyshev polynomials (~14x speedup). Set via `set_leb_file()` or `LIBEPHEMERIS_LEB` env var.
