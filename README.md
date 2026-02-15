# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

**The most scientifically rigorous astronomical ephemeris library for Python. Drop-in replacement for `pyswisseph`, powered by NASA JPL DE440/DE441 and modern IAU standards.**

> [!WARNING]
> **Pre-Alpha** -- The public API may change without notice. Not yet recommended for production use.

---

## Why LibEphemeris

Swiss Ephemeris has been the industry standard for decades. LibEphemeris is a ground-up replacement that uses newer science, newer data, and a more transparent architecture.

### More precise where it matters

| Component | LibEphemeris | Swiss Ephemeris |
|---|---|---|
| **Ephemeris** | JPL DE440/DE441 (Park et al. 2021) | JPL DE431 (2013) |
| **Nutation** | IAU 2006/2000A -- **1365 terms**, ~0.01 mas | IAU 2006/2000A -- 1365 terms, ~0.01 mas |
| **Precession** | IAU 2006 (Capitaine et al. 2003) | IAU 2006 |
| **Outer planet positions** | **Planet centers** (corrected from barycenters) | System barycenters (no correction by default) |
| **Delta T** | Stephenson et al. 2016 (newer, more accurate for historical dates) | Espenak & Meeus 2006 |
| **True Node** | Geometric **h = r x v** from JPL state vectors | Analytical perturbation series |

DE440 incorporates 8 years more observational data than DE431, including Juno-era Jupiter data and MESSENGER-era Mercury data, and uses the updated ICRF 3.0 reference frame. Lunar Laser Ranging precision: ~1 milliarcsecond.

### Planet center corrections

Swiss Ephemeris returns **system barycenters** for Jupiter, Saturn, Uranus, Neptune, and Pluto -- the center of mass of the planet plus all its moons. LibEphemeris corrects this automatically using bundled SPK segments and analytical satellite theories (Lieske E5, TASS 1.7, Keplerian models for Triton and Charon). No configuration needed.

### Full pyswisseph compatibility

1:1 API compatibility. Swap `import swisseph` for `import libephemeris`. Same functions, same constants, same return types.

### Readable, testable, open

Pure Python. Every algorithm is documented, every model is cited. LGPL-3.0. Thread-safe Context API for concurrent applications.

### The trade-off

LibEphemeris is slower than Swiss Ephemeris. It is pure Python; Swiss Ephemeris is C. For a single natal chart this is imperceptible. For batch processing millions of charts, it matters. A Rust core (Milestone 2) will close this gap.

| Operation | pyswisseph (C) | libephemeris (Python) | Ratio |
|---|---|---|---|
| `calc_ut` (5000 calls) | 0.04s | 4.2s | ~117x |
| `houses` (1200 calls) | 0.004s | 0.13s | ~31x |
| `get_ayanamsa_ut` (1500 calls) | 0.001s | 0.03s | ~29x |

*Apple M1 Pro, Python 3.10.*

---

## Installation

```bash
pip install libephemeris
```

DE440 (~128 MB) is downloaded automatically on first use. For planet center corrections on outer planets:

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

### Planetary positions

```python
import libephemeris as swe
from libephemeris.constants import *

jd = swe.julday(2000, 1, 1, 12.0)  # J2000.0

sun, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SPEED)
print(f"Sun: {sun[0]:.6f}°  speed: {sun[3]:.6f}°/day")

moon, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SPEED)
print(f"Moon: {moon[0]:.6f}°")
```

### Houses and angles

```python
jd = swe.julday(2024, 11, 5, 18.0)
cusps, ascmc = swe.houses(jd, 41.9028, 12.4964, b"P")  # Placidus, Rome

print(f"ASC: {ascmc[0]:.2f}°   MC: {ascmc[1]:.2f}°")
for i in range(1, 13):
    print(f"  House {i:2d}: {cusps[i]:.2f}°")
```

### Sidereal positions

```python
swe.set_sid_mode(SE_SIDM_LAHIRI)
ayanamsha = swe.get_ayanamsa_ut(jd)
sun_sid, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SIDEREAL | SEFLG_SPEED)
print(f"Lahiri ayanamsha: {ayanamsha:.6f}°")
print(f"Sidereal Sun: {sun_sid[0]:.6f}°")
```

### Event finding

```python
# Next Aries ingress after jd
ingress_jd = swe.solcross_ut(0.0, jd, SEFLG_SWIEPH)

# Next time Moon crosses 0° longitude
moon_cross = swe.mooncross_ut(0.0, jd, SEFLG_SWIEPH)
```

### Eclipses

```python
# Next solar eclipse after J2000
ecl_type, times = swe.sol_eclipse_when_glob(jd, SEFLG_SWIEPH)
max_jd = times[0]
y, m, d, h = swe.revjul(max_jd, SE_GREG_CAL)
print(f"Next solar eclipse: {y}-{m:02d}-{d:02d}")
```

### Extended date range

```python
# Use DE441 for historical/future dates (-13200 to +17191 CE)
swe.set_ephemeris_file("de441.bsp")

jd_ancient = swe.julday(-2500, 6, 21, 12.0)  # Summer solstice, 2500 BC
sun, _ = swe.calc_ut(jd_ancient, SE_SUN, SEFLG_SPEED)
```

---

## Scientific foundations

The full calculation pipeline, with models, term counts, and measured precision, is documented in [PRECISION.md](docs/PRECISION.md). Summary below.

### Nutation: IAU 2006/2000A (1365 terms)

All ecliptic coordinates use the IAU 2006/2000A nutation model via the ERFA library (`erfa.nut06a()`): 678 lunisolar + 687 planetary terms, ~0.01--0.05 milliarcsecond precision. This is the highest-precision nutation model adopted by the IAU.

`SEFLG_MOSEPH` is accepted for API compatibility but silently ignored. LibEphemeris never falls back to semi-analytical models -- all calculations use the full JPL numerical integration.

### Precession: IAU 2006

IAU 2006 precession (Capitaine et al. 2003) via `erfa.pmat06()`, Fukushima-Williams four-angle formulation with terms to T^5. Frame bias from ICRS to J2000 is applied via the combined bias-precession-nutation matrix (`erfa.pnm06a()`).

### Planet centers vs barycenters

JPL ephemerides give **system barycenters** for outer planets. LibEphemeris corrects to **body centers** automatically:

1. **SPK planet centers** (<0.001 arcsec): bundled `planet_centers.bsp` with segments for Jupiter (599), Saturn (699), Uranus (799), Neptune (899), Pluto (999) from JPL satellite SPK files.
2. **Analytical satellite theories** (sub-arcsecond fallback): Lieske E5 for Galilean moons (21--50 terms/satellite), TASS 1.7 for Saturn (Vienne & Duriez 1995, ~1500 terms), Keplerian/J2 for Triton (Jacobson 2009), two-body Keplerian for Charon (Brozovic & Jacobson 2024).
3. **Raw barycenters** (last resort if both above fail).

Standard Swiss Ephemeris returns barycenters by default. The `SEFLG_CENTER_BODY` flag exists but is rarely used and requires additional satellite data.

### Delta T: Stephenson et al. 2016

Cubic spline from 720 BC to ~2016 AD (Table S15), parabolic model outside, optional IERS observed values for 1973--present. This is newer and more accurate for historical dates than the Espenak & Meeus (2006) polynomials used by Swiss Ephemeris.

### True Node: geometric from JPL state vectors

Angular momentum **h = r x v** from JPL DE440 Moon position/velocity, refined with 120+ ELP2000-82B perturbation terms. Computes the True Node by definition: the intersection of the instantaneous orbital plane with the ecliptic.

---

## Measured precision vs Swiss Ephemeris

Tested at 100+ random dates within DE440 range (1550--2650 CE):

| Planet | Max difference | Mean difference |
|--------|---------------|-----------------|
| Sun | 0.20" | 0.04" |
| Moon | 3.32" | 0.70" |
| Mercury | 0.32" | 0.05" |
| Venus | 0.33" | 0.08" |
| Mars | 0.58" | 0.06" |
| Jupiter | 0.44" | 0.12" |
| Saturn | 0.51" | 0.13" |
| Uranus | 0.50" | 0.23" |
| Neptune | 1.17" | 0.24" |
| Pluto | 0.75" | 0.26" |

All sub-arcsecond except Moon (~3.3"), which reflects different COB correction pipelines. Both values are well within the intrinsic precision of the JPL ephemeris.

| Component | Precision |
|-----------|-----------|
| Ascendant | <0.001° |
| MC | <0.001° |
| House cusps | 0.001°--0.01° |
| Ayanamshas (formula) | <0.0002° |
| Ayanamshas (star-based) | <0.006° |
| Eclipse timing | <10 seconds |
| Crossing events | <0.1 arcsec |

See [PRECISION.md](docs/PRECISION.md) for the full analysis.

---

## Features

### Celestial bodies

| Category | Bodies |
|---|---|
| **Planets** | Sun, Moon, Mercury -- Pluto (+ Earth for heliocentric) |
| **Lunar points** | Mean Node, True Node, Mean Lilith, True Lilith, Interpolated Apogee/Perigee |
| **Main belt** | Ceres, Pallas, Juno, Vesta (+ any numbered asteroid via SBDB) |
| **Centaurs** | Chiron, Pholus, Nessus, Asbolus, Chariklo |
| **TNOs** | Eris, Sedna, Haumea, Makemake, Orcus, Ixion, Quaoar, Gonggong, Varuna |
| **Near-Earth** | Apophis, Eros, Bennu, Ryugu, Itokawa, Toutatis |
| **Hypothetical** | 8 Hamburg School Uranian planets, Transpluto, Vulcan, White Moon, Proserpina, Waldemath, Pickering/Leverrier/Lowell/Adams Planet X |
| **Fixed stars** | 102 stars from Hipparcos catalog (Royal Stars, Behenian stars, Pleiades, full zodiacal coverage) |
| **Planetary moons** | Galilean (Io, Europa, Ganymede, Callisto), Saturn (Mimas -- Iapetus), Uranus (Miranda -- Oberon), Triton, Phobos/Deimos, Charon |
| **Arabic parts** | Part of Fortune, Spirit, Eros, Faith |
| **Angles** | Ascendant, MC, Descendant, IC, Vertex, Antivertex |

### House systems (24)

`P` Placidus &bull; `K` Koch &bull; `O` Porphyry &bull; `R` Regiomontanus &bull; `C` Campanus &bull; `E` Equal (Asc) &bull; `A`/`D` Equal (MC) &bull; `W` Whole Sign &bull; `M` Morinus &bull; `B` Alcabitius &bull; `T` Polich/Page (Topocentric) &bull; `U` Krusinski &bull; `G` Gauquelin &bull; `V` Vehlow &bull; `X` Meridian &bull; `H` Horizontal &bull; `F` Carter &bull; `S` Sripati &bull; `L` Pullen SD &bull; `Q` Pullen SR &bull; `N` Natural Gradient &bull; `Y` APC &bull; `I`/`i` Sunshine

### Sidereal modes (43 ayanamshas)

Western (Fagan/Bradley, De Luce, Djwhal Khul), Indian (Lahiri, Raman, Krishnamurti, Yukteshwar, Suryasiddhanta, Aryabhata), Babylonian (Kugler 1-3, Huber, ETPSC, Britton), star-based (Aldebaran 15Tau, True Citra, True Revati, True Pushya, True Mula, True Sheoran), galactic (Galactic Center, Galactic Equator variants), historical epochs (Hipparchos, Sassanian, J2000, J1900, B1950), Valens Moon, and user-defined.

### Calculation flags

| Flag | Effect |
|---|---|
| `SEFLG_SPEED` | Include velocities (6-component output) |
| `SEFLG_TOPOCTR` | Topocentric positions (requires `set_topo()`) |
| `SEFLG_HELCTR` | Heliocentric positions |
| `SEFLG_BARYCTR` | Barycentric positions |
| `SEFLG_EQUATORIAL` | Right ascension / declination |
| `SEFLG_XYZ` | Cartesian (X, Y, Z) coordinates |
| `SEFLG_SIDEREAL` | Sidereal positions (requires `set_sid_mode()`) |
| `SEFLG_J2000` | J2000.0 reference frame |
| `SEFLG_NONUT` | No nutation |
| `SEFLG_NOABERR` | No aberration correction |
| `SEFLG_NOGDEFL` | No gravitational deflection |
| `SEFLG_TRUEPOS` | True geometric position (no light-time) |
| `SEFLG_RADIANS` | Output in radians |
| `SEFLG_ICRS` | ICRS reference frame |

### Eclipses and occultations

- **Solar eclipses:** global search, local circumstances, contact times (C1-C4), path width, central line, magnitude, obscuration, Besselian elements.
- **Lunar eclipses:** global search, local circumstances, umbral/penumbral magnitude and contact times, gamma parameter, duration.
- **Occultations:** lunar and planetary occultations of stars.
- **Saros and Inex** series number computation.

### Heliacal events

Schaefer (1990) atmospheric visibility model: Rayleigh/aerosol/ozone extinction, twilight sky brightness, Ptolemaic visibility thresholds, Schaefer contrast model with configurable observer skill levels.

### Rise, set, and transit

`rise_trans()` and `rise_trans_true_hor()` with disc-center, disc-bottom, and true-horizon modes. Configurable refraction and twilight definitions (civil, nautical, astronomical).

### Minor bodies via SPK

Default Keplerian elements give arcminute-level accuracy. For sub-arcsecond precision, download SPK kernels from JPL Horizons:

```python
swe.download_and_register_spk(
    body="Chiron", ipl=SE_CHIRON,
    start="2000-01-01", end="2100-01-01",
)
pos, _ = swe.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)
```

Or enable automatic download: `swe.set_auto_spk_download(True)`.

Strict precision mode (default) prevents accidentally using low-precision Keplerian fallbacks for Chiron, Ceres, Pallas, Juno, and Vesta.

### Planetary phenomena

Phase angle, elongation, apparent diameter, apparent magnitude. Hapke model for Moon, Mallama (2018) for Pluto, IAU standard formulas for planets.

### N-body orbit propagation (optional)

```python
# pip install libephemeris[nbody]
from libephemeris import propagate_orbit_assist
result = propagate_orbit_assist(SE_CHIRON, jd_start, jd_end, step_days=1.0)
```

### Time utilities

Julian Day conversions, Delta T (Stephenson 2016 + optional IERS), UTC/TT/TAI conversions, sidereal time, equation of time, LMT/LAT conversions.

### Coordinate transformations

`cotrans()` for ecliptic-equatorial conversion, `azalt()` for horizontal coordinates, `refrac()` for atmospheric refraction.

---

## Configuring the ephemeris

Default: **DE440** (1550--2650 CE, ~128 MB, auto-downloaded).

```python
swe.set_ephemeris_file("de441.bsp")   # -13200 to +17191 CE, ~3.4 GB
swe.set_ephe_path("/data/jpl-kernels")  # Custom directory
```

Or via environment variable: `export LIBEPHEMERIS_EPHEMERIS=de441.bsp`

Priority: `set_ephemeris_file()` > `LIBEPHEMERIS_EPHEMERIS` env var > default `de440.bsp`.

| Kernel | Date range | Size | Notes |
|---|---|---|---|
| `de440.bsp` | 1550 -- 2650 | ~128 MB | Default. Modern dates, highest data density. |
| `de441.bsp` | -13200 -- +17191 | ~3.4 GB | Same precision as DE440 for historical/future work. |
| `de422.bsp` | -3000 -- +3000 | ~623 MB | Older generation (2009). |
| `de431.bsp` | -13200 -- +17191 | ~3.4 GB | Used by Swiss Ephemeris. Superseded by DE441. |

Dates outside the loaded kernel's range raise `EphemerisRangeError`.

---

## Thread safety

### Global API (pyswisseph-compatible)

```python
import libephemeris as swe

swe.set_topo(12.5, 41.9, 0)
pos, _ = swe.calc_ut(2451545.0, SE_SUN, 0)
```

Uses global state -- not thread-safe, exactly like `pyswisseph`.

### Context API (thread-safe)

```python
from libephemeris import EphemerisContext, SE_SUN

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)
sun, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
```

Each context has isolated state (observer location, sidereal mode, SPK registrations, angle cache) while sharing ephemeris data. Memory overhead per context: ~1 KB.

```python
from concurrent.futures import ThreadPoolExecutor

def calc_sun(location):
    ctx = EphemerisContext()
    ctx.set_topo(location["lon"], location["lat"], 0)
    pos, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
    return {"name": location["name"], "sun": pos[0]}

locations = [
    {"name": "Rome", "lon": 12.5, "lat": 41.9},
    {"name": "London", "lon": -0.1, "lat": 51.5},
    {"name": "Tokyo", "lon": 139.7, "lat": 35.7},
]

with ThreadPoolExecutor(max_workers=3) as pool:
    results = list(pool.map(calc_sun, locations))
```

---

## pyswisseph compatibility

LibEphemeris exports every function with both `swe_` prefix and bare name:

```python
swe.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)  # pyswisseph style
swe.calc_ut(jd, SE_SUN, SEFLG_SPEED)       # short form (identical)
```

Same constants (`SE_SUN`, `SEFLG_SPEED`, `SE_SIDM_LAHIRI`, etc.) and same return types.

**Known differences:**
- Dates outside the loaded JPL kernel raise `EphemerisRangeError` instead of silently degrading to Moshier.
- True Node uses geometric `h = r x v` from JPL state vectors rather than perturbation series (~0.06° practical difference; the geometric method is more rigorous by construction).
- Outer planet positions are corrected from system barycenters to planet body centers.

---

## Exception hierarchy

```
Error (base, pyswisseph-compatible)
├── InputValidationError
│   ├── CoordinateError          # invalid lat/lon
│   └── InvalidBodyError         # body not valid for operation
├── DataNotFoundError
│   ├── UnknownBodyError         # unknown body ID
│   ├── StarNotFoundError        # star not in catalog
│   ├── SPKNotFoundError         # SPK file missing
│   └── SPKRequiredError         # SPK needed in strict mode
├── CalculationError
│   ├── PolarCircleError         # houses at polar latitudes
│   ├── EphemerisRangeError      # date outside kernel range
│   └── ConvergenceError         # algorithm didn't converge
└── ConfigurationError           # missing config (e.g. set_topo)
```

---

## Development

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris
uv pip install -e ".[dev]"

poe test           # fast tests (excludes @pytest.mark.slow)
poe test:full      # all tests
poe lint           # ruff check --fix
poe format         # ruff format
poe typecheck      # mypy
poe coverage       # pytest with coverage
```

### Project layout

```
libephemeris/
├── __init__.py            # Public API (all exports and aliases)
├── constants.py           # SE_*, SEFLG_*, SE_SIDM_* constants
├── planets.py             # Core planetary calculations
├── houses.py              # 24 house system implementations
├── lunar.py               # Lunar nodes, Lilith, interpolated apsides
├── moon_theories/         # Satellite theories for COB corrections
│   ├── galilean.py        # Lieske E5 (Jupiter's Galilean moons)
│   ├── tass17.py          # TASS 1.7 (Saturn's moons, Vienne & Duriez 1995)
│   ├── triton.py          # Triton Keplerian (Jacobson 2009)
│   └── charon.py          # Charon two-body (Brozovic & Jacobson 2024)
├── minor_bodies.py        # Asteroids, TNOs (Keplerian + SBDB lookup)
├── eclipse.py             # Solar/lunar eclipses, Besselian elements
├── crossing.py            # Longitude crossing event search
├── heliacal.py            # Heliacal rising/setting (Schaefer 1990)
├── extinction.py          # Atmospheric extinction and twilight models
├── fixed_stars.py         # 102 Hipparcos stars with proper motion
├── hypothetical.py        # Uranian planets, seorbel.txt parser
├── planetary_moons.py     # Galilean moons, Titan, etc.
├── arabic_parts.py        # Arabic parts (Lots)
├── spk.py                 # SPK kernel download and registration
├── rebound_integration.py # N-body integration (optional)
├── context.py             # Thread-safe EphemerisContext
├── state.py               # Global state management
├── time_utils.py          # Julian Day, Delta T, sidereal time
├── astrometry.py          # IAU 2006 precession, IAU 2006/2000A nutation
├── utils.py               # Coordinate transforms, angle utilities
├── exceptions.py          # Exception hierarchy
├── iers_data.py           # IERS Delta T observed data
└── cli.py                 # CLI entry point
tests/                     # Comprehensive test suite
compare_scripts/           # Swiss Ephemeris comparison tools
docs/
├── PRECISION.md           # Full scientific precision analysis
├── PRECISION_TUNING.md    # Precision configuration guide
└── migration-guide.md     # pyswisseph migration guide
```

---

## Roadmap

1. **Milestone 1 (current):** Pure-Python library, 1:1 pyswisseph drop-in with modern IAU standards.
2. **Milestone 2:** Rust core for C-level performance with the same Python API.

---

## License

**LGPL-3.0.** You can use LibEphemeris as a dependency in proprietary software without releasing your source code. If you modify LibEphemeris itself, those modifications must be released under the LGPL. See [LICENSE](LICENSE).
