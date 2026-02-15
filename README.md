# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

Pure-Python astronomical ephemeris library focused on scientific rigor.

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
