# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

A pure Python astronomical ephemeris library based on NASA JPL data. Designed as a 1:1 API-compatible drop-in replacement for PySwisseph, strictly focused on scientific precision.

- Drop-in replacement for `pyswisseph`
- Uses NASA JPL ephemerides (DE440/DE441) via Skyfield
- Uses modern IAU standards (ERFA/pyerfa)

> [!WARNING]
> **Pre-Alpha** -- The public API may change without notice.

---

## Why

Swiss Ephemeris is fast and widely used, but its architecture reflects design
decisions from the 1990s. LibEphemeris makes a different trade-off:
**absolute consistency with current NASA JPL ephemerides and IAU standards**,
at the cost of speed and strict numerical agreement with Swiss Ephemeris.

- **Modern ephemeris, no fallback.** All calculations use JPL DE440/DE441 (2021) exclusively. The `SEFLG_MOSEPH` flag is accepted for compatibility but silently ignored -- there is no degradation to 1980s analytical theories.
- **Physical planet positions.** Outer planet positions are automatically corrected from system barycenters to true planet body centers, using JPL satellite ephemerides and analytical moon theories. No extra flags or files required.
- **JPL-grounded lunar apsides.** The interpolated perigee and apogee are derived from actual physical distance-extrema passages in JPL data, rather than from truncated analytical series.
- **Current Delta T model.** Uses Stephenson, Morrison & Hohenkerk (2016) with optional IERS observational data -- the current astronomical standard for historical timescale conversion.
- **Transparent.** Pure Python, fully testable, readable, and documented with explicit astronomical references.
- **Slower than Swiss Ephemeris.** Swiss is C; LibEphemeris is Python (see Performance).

Methodology and rationale: `docs/methodology.md`.
Precision measurements: `docs/PRECISION.md`.

---

## Installation

```bash
pip install libephemeris
```

DE440 (~128 MB) downloads automatically on first use.

Download all data files for your use case (ephemeris + SPK kernels):

```bash
libephemeris download:medium      # recommended for most users
```

See [CLI commands](#cli-commands) for all download tiers.

### Optional extras

```bash
pip install libephemeris[nbody]  # REBOUND/ASSIST n-body integration
pip install libephemeris[stars]  # Star catalog building (astropy)
pip install libephemeris[all]    # Everything
```

**Requirements:** Python 3.9+ &bull; skyfield >= 1.54 &bull; pyerfa >= 2.0

---

## Quick start

```python
import libephemeris as swe
from libephemeris.constants import SE_SUN, SE_MOON, SEFLG_SPEED

jd = swe.julday(2000, 1, 1, 12.0)  # J2000.0

sun, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SPEED)
moon, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SPEED)

print("Sun lon:", sun[0], "deg")
print("Moon lon:", moon[0], "deg")
```

Houses:

```python
jd = swe.julday(2024, 11, 5, 18.0)
cusps, ascmc = swe.houses(jd, 41.9028, 12.4964, b"P")  # Placidus
print("ASC:", ascmc[0], "deg")
print("MC:", ascmc[1], "deg")
```

---

## Calculation flags

Flags are bitmasks that control what is calculated and how the result is returned. `calc_ut()` returns a 6-element tuple `(longitude, latitude, distance, speed_lon, speed_lat, speed_dist)` and a return flag. Combine multiple flags with bitwise OR (`|`).

### Velocity

- `SEFLG_SPEED` — Populates the speed fields (`pos[3]`–`pos[5]`) with daily motion in longitude, latitude, and distance. Without this flag those values are zero. Almost every call should include it.

### Observer

By default the observer is at Earth's center (geocentric).

- `SEFLG_HELCTR` — Heliocentric: moves the observer to the Sun. Distances become heliocentric AU.
- `SEFLG_TOPOCTR` — Topocentric: places the observer on Earth's surface at the position set with `swe_set_topo()`. This matters most for the Moon (up to ~1° parallax).

### Coordinates

By default the output is ecliptic longitude/latitude of date.

- `SEFLG_EQUATORIAL` — Switches to equatorial coordinates: `pos[0]` becomes Right Ascension (0–360°) and `pos[1]` becomes Declination (±90°). Speeds change accordingly.

### Reference frame

By default positions are precessed to the equinox of date.

- `SEFLG_J2000` — Keeps coordinates in the J2000.0 reference frame instead of precessing to the equinox of date.
- `SEFLG_NONUT` — Excludes nutation, giving positions on the mean ecliptic/equator.

### Position corrections

By default positions are apparent (light-time and aberration corrected).

- `SEFLG_TRUEPOS` — Geometric position: no light-time correction. Returns where the body actually is at the instant of calculation.
- `SEFLG_NOABERR` — Astrometric position: light-time corrected but no aberration. Comparable to star catalog positions.
- `SEFLG_ASTROMETRIC` — Convenience shorthand for `SEFLG_NOABERR | SEFLG_NOGDEFL`.

### Sidereal zodiac

- `SEFLG_SIDEREAL` — Subtracts the ayanamsha from ecliptic longitude, returning sidereal rather than tropical positions. Requires a prior `swe_set_sid_mode()` call to select the ayanamsha (Lahiri, Fagan-Bradley, etc.).

### Combining flags

```python
# Heliocentric position with velocity
pos, _ = swe.calc_ut(jd, SE_MARS, SEFLG_SPEED | SEFLG_HELCTR)

# Sidereal equatorial coordinates
pos, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_SIDEREAL)
```

> [!NOTE]
> `SEFLG_MOSEPH` and `SEFLG_SWIEPH` are accepted for API compatibility but
> silently ignored — all calculations always use JPL DE440/DE441 via Skyfield.
> `SEFLG_BARYCTR` is mapped to heliocentric. `SEFLG_XYZ`, `SEFLG_RADIANS`,
> `SEFLG_SPEED3`, and `SEFLG_ICRS` are defined but not yet implemented.

---

## Choose the ephemeris (DE440 vs DE441)

Default ephemeris: **DE440** (1550--2650 CE).

To use **DE441** (extended range -13200 to +17191 CE):

```python
import libephemeris as swe
swe.set_ephemeris_file("de441.bsp")
```

Or via environment variable:

```bash
export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

Dates outside the loaded kernel raise `EphemerisRangeError`.

---

## Outer-planet centers (barycenter vs center)

JPL DE ephemerides provide system barycenters for Jupiter/Saturn/Uranus/Neptune/Pluto.
LibEphemeris corrects these to planet body centers automatically:

1. Uses tier-specific SPK files (`planet_centers_{tier}.bsp`) downloaded by `libephemeris download:<tier>`
2. Falls back to analytical satellite models when SPK coverage is not available:
    - Jupiter: Galilean moon theory (E5/Meeus) — ~0.05 arcsec
    - Saturn: TASS 1.7 (Titan-dominated) — ~0.02 arcsec
    - Uranus: 5 major moons Keplerian — ~0.01 arcsec
    - Neptune: Triton Keplerian — ~0.003 arcsec
    - Pluto: Charon two-body — ~0.15 arcsec

SPK coverage varies by tier:

| Planet  | base (1850-2150) | medium (1550-2650) | extended                            |
| ------- | ---------------- | ------------------ | ----------------------------------- |
| Jupiter | SPK              | SPK                | SPK 1600-2200, fallback outside     |
| Saturn  | SPK              | SPK                | SPK -502 to +4500, fallback outside |
| Uranus  | SPK              | SPK                | SPK full (-12000 to +17000)         |
| Neptune | SPK              | SPK                | SPK full (-12000 to +17000)         |
| Pluto   | SPK              | SPK 1800-2200      | SPK 1800-2200, fallback outside     |

Full technical details are in `docs/PRECISION.md`.

---

## Minor bodies (SPK)

For many asteroids and TNOs, high precision requires JPL SPK kernels.
SPK kernels are downloaded automatically via `astroquery` (a core dependency).

```python
import libephemeris as swe
from libephemeris.constants import SE_CHIRON

swe.set_auto_spk_download(True)
swe.set_spk_cache_dir("./spk_cache")

pos, _ = swe.calc_ut(2460000.0, SE_CHIRON, 0)
print(pos[0])
```

### Optional Dependencies

LibEphemeris has several optional dependencies for enhanced functionality:

| Extra     | Description           | Dependencies        |
| --------- | --------------------- | ------------------- |
| `[stars]` | Star catalog building | `astropy`           |
| `[nbody]` | N-body integration    | `rebound`, `assist` |
| `[all]`   | All optional features | All above           |

```bash
pip install libephemeris[stars]  # For star catalog building via astropy
pip install libephemeris[nbody]  # For REBOUND/ASSIST n-body fallback
pip install libephemeris[all]    # Install all optional dependencies
```

**Note:** `pyerfa` and `astroquery` are required dependencies. `pyerfa` provides IAU 2006/2000A precession-nutation models; `astroquery` enables automatic SPK downloads from JPL Horizons.

---

## Thread safety

The global API is pyswisseph-compatible and uses global mutable state.
For concurrency, use `EphemerisContext`:

```python
from libephemeris import EphemerisContext, SE_SUN

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)
pos, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
```

---

## Docs

- `docs/methodology.md` (computational methodology, comparison with Swiss Ephemeris)
- `docs/PRECISION.md` (scientific models, measured precision, references)
- `docs/migration-guide.md` (pyswisseph -> libephemeris)
- `docs/HOUSE_SYSTEMS.md`
- `docs/AYANAMSHA.md`
- `docs/testing.md`
- `docs/methodology_lunar_apsides.md` (lunar apsides methodology)
- `docs/interpolated_perigee_methodology.md` (perigee calibration method and precision)
- `docs/INTERPOLATED_APOGEE.md` (apogee methodology)
- `docs/TRUE_LILITH_METHODS.md` (True Lilith correction methods)
- `docs/PLANET_CENTERS_SPK.md` (planet centers SPK system)
- `docs/PRECISION_TUNING.md` (precision tuning guide)
- `docs/PYERFA_BENEFITS.md` (pyerfa integration benefits)
- `docs/REBOUND_BENEFITS.md` (REBOUND n-body integration benefits)

---

## Performance

LibEphemeris is pure Python; Swiss Ephemeris is C. Expect LibEphemeris to be slower.
For batch workloads, use `EphemerisContext`, parallelism, and caching.

---

## CLI commands

Download all required data for your precision tier:

| Command                          | Ephemeris  | SPK range | Download |
| -------------------------------- | ---------- | --------- | -------- |
| `libephemeris download:base`     | de440s.bsp | 1850-2150 | ~35 MB   |
| `libephemeris download:medium`   | de440.bsp  | 1900-2100 | ~130 MB  |
| `libephemeris download:extended` | de441.bsp  | 1600-2500 | ~3.3 GB  |

Other commands:

| Command                  | Description                     |
| ------------------------ | ------------------------------- |
| `libephemeris status`    | Show installed data file status |
| `libephemeris --version` | Show version information        |

---

## Development

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris
uv pip install -e ".[dev]"
```

All development tasks use [poethepoet](https://poethepoet.naberhaus.dev/) (`poe`).

### Code quality

| Command         | Description                 |
| --------------- | --------------------------- |
| `poe format`    | Format code with Ruff       |
| `poe lint`      | Lint and auto-fix with Ruff |
| `poe typecheck` | Type-check with mypy        |

### Tests

| Command                   | Description                                               |
| ------------------------- | --------------------------------------------------------- |
| `poe test`                | Fast tests (excludes `@pytest.mark.slow`)                 |
| `poe test:full`           | All tests including slow ones                             |
| `poe test:fast`           | Fast tests, parallel (`-n auto`)                          |
| `poe test:fast:essential` | Parallel essential subset (~670 tests, 1 file per module) |
| `poe test:unit`           | Unit tests only (`tests/`)                                |
| `poe test:unit:fast`      | Unit tests, parallel                                      |
| `poe test:compare`        | Comparison tests vs pyswisseph (`compare_scripts/tests/`) |
| `poe test:compare:fast`   | Comparison tests, parallel                                |
| `poe test:lunar`          | All lunar tests (nodes, Lilith, perigee, apogee)          |
| `poe test:lunar:perigee`  | Perigee tests (perturbations + interpolated + osculating) |
| `poe test:lunar:apogee`   | Apogee tests (perturbations + interpolated)               |
| `poe test:lunar:lilith`   | Lilith tests (mean + true, all correction modules)        |
| `poe coverage`            | Fast tests with coverage report                           |
| `poe coverage:full`       | All tests with coverage report                            |

### Tier diagnostics

Run diagnostic tables showing all celestial bodies with coordinates, velocities, and data source for each precision tier:

| Command             | Description                                     |
| ------------------- | ----------------------------------------------- |
| `poe diag:base`     | Diagnostic for base tier (1850-2150)            |
| `poe diag:medium`   | Diagnostic for medium tier (1550-2650)          |
| `poe diag:extended` | Diagnostic for extended tier (-13200 to +17191) |

### SPK downloads (dev)

Download SPK kernels for minor bodies directly (without the full CLI tier setup):

| Command                     | Description                                           |
| --------------------------- | ----------------------------------------------------- |
| `poe spk:download:base`     | SPK for base tier range (1850-2150)                   |
| `poe spk:download:medium`   | SPK for medium tier range (1900-2100)                 |
| `poe spk:download:extended` | Max-range SPK files (1600-2500, single file per body) |

### Data generation

Generate `planet_centers_*.bsp` files for each precision tier. Requires `spiceypy >= 6.0.0`.

| Command                                | Description                                               | Download |
| -------------------------------------- | --------------------------------------------------------- | -------- |
| `poe generate-planet-centers:base`     | Generate for base tier (1850-2150)                        | ~500 MB  |
| `poe generate-planet-centers:medium`   | Generate for medium tier (1550-2650)                      | ~4 GB    |
| `poe generate-planet-centers:extended` | Generate for extended tier (partial)                      | ~6.5 GB  |
| `poe generate-planet-centers:all`      | Generate all 3 tiers                                      | ~11 GB   |
| `poe generate-lunar-corrections`       | Regenerate lunar correction tables (requires `de441.bsp`) | —        |

### Calibration

Calibrate perturbation series coefficients against JPL DE441 ephemeris. See `docs/interpolated_perigee_methodology.md` for full details.

| Command                       | Description                                                            |
| ----------------------------- | ---------------------------------------------------------------------- |
| `poe calibrate-perigee`       | Full perigee calibration, 1500-2500 CE (~30 min, requires `de441.bsp`) |
| `poe calibrate-perigee:quick` | Quick validation run, 100-year range (~2 min)                          |

After calibration, update coefficients in `lunar.py`, sync `generate_lunar_corrections.py`, and run `poe generate-lunar-corrections` followed by `poe test:lunar:perigee`.

Generated files are saved in the workspace root:

- `planet_centers_base.bsp` (~15-20 MB)
- `planet_centers_medium.bsp` (~40-50 MB)
- `planet_centers_extended.bsp` (~80-100 MB)

---

## License

LGPL-3.0. See `LICENSE`.
