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

Swiss Ephemeris is fast and widely used. LibEphemeris makes a different trade-off:

- **Accuracy-first:** JPL numerical integrations (DE440/DE441) + IAU-standard precession/nutation
- **Transparent:** pure Python, fully testable, and documented with references
- **Slower than SwissEph:** Swiss is C; LibEphemeris is Python (see Performance)

Precision details (models, term counts, measured comparisons, references): `docs/PRECISION.md`.

---

## Installation

```bash
pip install libephemeris
```

DE440 (~128 MB) downloads automatically on first use.

Optional data for outer-planet center corrections:

```bash
libephemeris download-data
```

### Optional extras

```bash
pip install libephemeris[spk]    # Automatic SPK downloads from JPL Horizons
pip install libephemeris[nbody]  # REBOUND/ASSIST n-body integration
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

1. Uses a compact SPK file (`planet_centers.bsp`) downloaded by `libephemeris download-data`
2. Falls back to analytical satellite models when SPK coverage is not available

Full technical details are in `docs/PRECISION.md`.

---

## Minor bodies (SPK)

For many asteroids and TNOs, high precision requires JPL SPK kernels.

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

| Extra | Description | Dependencies |
|-------|-------------|--------------|
| `[spk]` | Automatic SPK downloads from JPL Horizons | `astroquery` |
| `[stars]` | Star catalog access | `astropy` |
| `[nbody]` | N-body integration | `rebound`, `reboundx` |
| `[all]` | All optional features | All above |

```bash
pip install libephemeris[spk]    # For auto-downloading SPK kernels from Horizons
pip install libephemeris[stars]  # For accessing star catalogs via astropy
pip install libephemeris[all]    # Install all optional dependencies
```

**Note:** `pyerfa` is a required dependency and provides IAU 2006/2000A precession-nutation models for high-precision calculations.

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

- `docs/PRECISION.md` (scientific models, measured precision, references)
- `docs/migration-guide.md` (pyswisseph -> libephemeris)
- `docs/HOUSE_SYSTEMS.md`
- `docs/AYANAMSHA.md`
- `docs/testing.md`

---

## Performance

LibEphemeris is pure Python; Swiss Ephemeris is C. Expect LibEphemeris to be slower.
For batch workloads, use `EphemerisContext`, parallelism, and caching.

---

## Development

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris
uv pip install -e ".[dev]"

poe test           # fast tests
poe test:full      # all tests
poe lint           # ruff check --fix
poe format         # ruff format
poe typecheck      # mypy
poe coverage       # pytest with coverage
```

---

## License

LGPL-3.0. See `LICENSE`.
