# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

**A pure-Python, high-precision astronomical ephemeris library. Drop-in replacement for `pyswisseph`, powered by NASA JPL DE440 via Skyfield.**

> [!WARNING]
> **Pre-Alpha** -- The public API may change without notice. Not yet recommended for production use.

---

## Why LibEphemeris

Swiss Ephemeris is the industry standard, but its C codebase is opaque and hard to maintain. LibEphemeris is a ground-up rewrite in Python that:

- Uses **NASA JPL DE440/DE441** -- the same ephemerides that guide Mars rovers and the James Webb Space Telescope.
- Provides a **1:1 `pyswisseph`-compatible API** so you can swap `import swisseph` for `import libephemeris` with minimal changes.
- Is **fully open source** (LGPL-3.0), readable, and testable.
- Offers a **thread-safe Context API** for concurrent applications.

### Precision at a glance

| Component | vs Swiss Ephemeris | Notes |
|---|---|---|
| Planets | < 1 arcsec | Effectively identical |
| Moon | < 3.5 arcsec | Different lunar theories |
| True Node | ~206 arcsec RMS | Different methodology (both valid) |
| Houses | < 0.001° | Same algorithms |
| Ayanamshas | < 0.06° | 43 systems supported |
| Velocities | < 0.01°/day | Central-difference method |

DE440 uses ICRF 3.0, validated with Lunar Laser Ranging to ~1 milliarcsecond for the Moon. See [PRECISION.md](docs/PRECISION.md) for full tables.

---

## Installation

```bash
pip install libephemeris
```

The DE440 ephemeris file (~128 MB) is downloaded automatically on first use. For optional precision data on outer planets:

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
| **Fixed stars** | 100+ stars with Hipparcos proper motion (Regulus, Spica, Sirius, Aldebaran, Antares, Pleiades cluster, etc.) |
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

### Eclipse and occultation calculations

- **Solar eclipses:** global search, local circumstances, geographic coordinates, contact times (C1-C4), path width, central line, magnitude, obscuration, Besselian elements.
- **Lunar eclipses:** global search, local circumstances, umbral/penumbral magnitude and contact times (P1, U1, U2, U3, U4, P4), gamma parameter, duration.
- **Occultations:** lunar occultations and planetary occultations of stars.
- **Saros and Inex** series number computation.

### Heliacal events and atmospheric modelling

Schaefer (1990) atmospheric visibility model for heliacal rising/setting with Ptolemaic visibility thresholds. Includes atmospheric extinction (Rayleigh, aerosol, ozone, water vapour), twilight sky brightness, and Schaefer contrast threshold model with configurable observer skill levels.

### Rise, set, and transit

`rise_trans()` and `rise_trans_true_hor()` for computing rise, set, upper and lower culmination times. Supports disc-center, disc-bottom, and true-horizon modes with configurable refraction and twilight definitions (civil, nautical, astronomical).

### High-precision minor bodies via SPK

Default Keplerian elements give arcminute-level accuracy. For arcsecond-level precision, download SPK kernels directly from JPL Horizons:

```python
import libephemeris as swe

# One-liner: download + register
swe.download_and_register_spk(
    body="Chiron", ipl=SE_CHIRON,
    start="2000-01-01", end="2100-01-01",
)

# Now calc_ut uses the SPK automatically
pos, _ = swe.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)
```

Or enable automatic download:

```python
swe.set_auto_spk_download(True)
swe.set_spk_cache_dir("./spk_cache")
# First calculation triggers download
pos, _ = swe.calc_ut(2460000.0, SE_CHIRON, 0)
```

### Planetary phenomena

`pheno_ut()` / `pheno()` for computing phase angle, phase, elongation, apparent diameter, and apparent magnitude. Includes Hapke model for Moon, Mallama (2018) for Pluto, IAU standard formulas for planets. Elongation helpers: `get_elongation_from_sun()`, `is_morning_star()`, `is_evening_star()`.

### N-body orbit propagation (optional)

With `pip install libephemeris[nbody]`, use REBOUND and ASSIST for high-fidelity orbit propagation:

```python
from libephemeris import propagate_orbit_assist, compare_with_keplerian
result = propagate_orbit_assist(SE_CHIRON, jd_start, jd_end, step_days=1.0)
```

### Orbital elements and nodal/apsidal data

`get_orbital_elements_ut()` returns full Keplerian elements. `nod_aps_ut()` computes ascending/descending nodes and perihelion/aphelion for any body, with mean, osculating, and barycentric methods.

### Time utilities

Julian Day conversions (`julday`, `revjul`), Delta T (`deltat`, `deltat_ex`), UTC/TT/TAI conversions, sidereal time (`sidtime`, `sidtime0`), equation of time, LMT/LAT conversions, IERS observed Delta T with automatic data download.

### Coordinate transformations

`cotrans()` / `cotrans_sp()` for ecliptic-equatorial conversion, `azalt()` / `azalt_rev()` for horizontal coordinates, `refrac()` / `refrac_extended()` for atmospheric refraction.

---

## Thread safety

LibEphemeris provides two APIs:

### Global API (pyswisseph-compatible)

```python
import libephemeris as swe

swe.set_topo(12.5, 41.9, 0)
pos, _ = swe.calc_ut(2451545.0, SE_SUN, 0)
```

Uses global state -- not thread-safe, exactly like `pyswisseph`.

### Context API (thread-safe)

```python
from libephemeris import EphemerisContext, SE_SUN, SE_MOON

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)
sun, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
moon, _ = ctx.calc_ut(2451545.0, SE_MOON, 0)
```

Each `EphemerisContext` instance has its own isolated state (observer location, sidereal mode, SPK registrations, angle cache) while sharing the expensive ephemeris data across instances. Safe for `ThreadPoolExecutor` and similar patterns:

```python
from concurrent.futures import ThreadPoolExecutor
from libephemeris import EphemerisContext, SE_SUN

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

## Configuring the ephemeris

Default: **DE440** (1550--2650 CE, ~128 MB, auto-downloaded).

```python
# Use a different kernel
swe.set_ephemeris_file("de441.bsp")  # -13200 to +17191 CE, ~3.4 GB

# Use a custom directory
swe.set_ephe_path("/data/jpl-kernels")
```

Or via environment variable:

```bash
export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

Priority: `set_ephemeris_file()` > `LIBEPHEMERIS_EPHEMERIS` env var > default `de440.bsp`.

| Kernel | Date range | Size |
|---|---|---|
| `de440.bsp` | 1550 -- 2650 | ~128 MB |
| `de441.bsp` | -13200 -- +17191 | ~3.4 GB |
| `de422.bsp` | -3000 -- +3000 | ~623 MB |
| `de431.bsp` | -13200 -- +17191 | ~3.4 GB |

> `SEFLG_MOSEPH` is accepted for API compatibility but silently ignored. All calculations always use JPL ephemerides.

### Planet centers vs barycenters

JPL DE440/DE441 provide positions for outer planets (Jupiter, Saturn, Uranus, Neptune, Pluto) as **system barycenters** -- the center of mass of the planet plus all its moons -- not the planet center itself. For Jupiter the difference can reach ~0.01°.

LibEphemeris corrects this **automatically** using a three-tier strategy:

1. **SPK-based planet centers** (default, <0.001 arcsec): a bundled `planet_centers.bsp` file (~25 MB) contains precise center-of-body segments (NAIF IDs 599, 699, 799, 899, 999) extracted from JPL satellite ephemerides. Downloaded with `libephemeris download-data`.

2. **Analytical moon theories** (fallback if `planet_centers.bsp` is missing): computes the barycenter-to-center offset from satellite positions using Lieske E5 for Jupiter's Galilean moons, TASS 1.7 for Saturn (Titan dominates at 96% of moon mass), Keplerian elements with J2 precession for Triton, and a two-body solution for Charon.

3. **Raw barycenters** (last resort): used only if both methods above fail.

No configuration is needed -- the correction is applied transparently to all `calc_ut()` / `calc()` calls. Standard Swiss Ephemeris returns barycenters for outer planets by default and requires the `SEFLG_CENTER_BODY` flag for planet center positions.

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

Catch broad categories (`except CalculationError`) or specific errors (`except PolarCircleError`). All inherit from `Error`, which is compatible with `swisseph.Error`.

---

## pyswisseph compatibility

LibEphemeris exports every function with both `swe_` prefix and bare name:

```python
# These are identical:
swe.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
swe.calc_ut(jd, SE_SUN, SEFLG_SPEED)
```

Same constants (`SE_SUN`, `SEFLG_SPEED`, `SE_SIDM_LAHIRI`, etc.) and same return types. Both `SE_` and bare-name aliases are provided for all constants.

**Known differences:**
- Very long time ranges beyond the loaded JPL kernel raise `EphemerisRangeError` instead of silently degrading.
- True Node uses geometric `h = r x v` from JPL state vectors rather than perturbation series. Both approaches are valid; the difference is ~0.06° for practical purposes.

---

## Performance

LibEphemeris is pure Python; `pyswisseph` is C. The trade-off is readability and maintainability vs raw speed.

| Operation | pyswisseph (C) | libephemeris (Python) | Ratio |
|---|---|---|---|
| `calc_ut` (5000 calls) | 0.04s | 4.2s | ~117x |
| `houses` (1200 calls) | 0.004s | 0.13s | ~31x |
| `get_ayanamsa_ut` (1500 calls) | 0.001s | 0.03s | ~29x |

*Apple M1 Pro, Python 3.10. A future Rust core (Milestone 2) will close this gap.*

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
├── minor_bodies.py        # Asteroids, TNOs (Keplerian + SBDB lookup)
├── eclipse.py             # Solar/lunar eclipses, Besselian elements
├── crossing.py            # Longitude crossing event search
├── heliacal.py            # Heliacal rising/setting
├── extinction.py          # Atmospheric extinction and twilight models
├── fixed_stars.py         # Fixed star catalog with proper motion
├── hypothetical.py        # Uranian planets, seorbel.txt parser
├── planetary_moons.py     # Galilean moons, Titan, etc.
├── arabic_parts.py        # Arabic parts (Lots)
├── spk.py                 # SPK kernel download and registration
├── rebound_integration.py # N-body integration (optional)
├── context.py             # Thread-safe EphemerisContext
├── state.py               # Global state management
├── time_utils.py          # Julian Day, Delta T, sidereal time
├── astrometry.py          # Precession, nutation, aberration (IAU)
├── utils.py               # Coordinate transforms, angle utilities
├── exceptions.py          # Exception hierarchy
├── iers_data.py           # IERS Delta T observed data
└── cli.py                 # CLI entry point
tests/                     # Comprehensive test suite
compare_scripts/           # Swiss Ephemeris comparison tools
```

---

## Roadmap

1. **Milestone 1 (current):** Pure-Python library, 1:1 pyswisseph drop-in.
2. **Milestone 2:** Rust core via [Starfield](https://docs.rs/starfield/latest/starfield/) for C-level performance with the same Python API.

---

## License

**LGPL-3.0.** You can use LibEphemeris as a dependency in proprietary software without releasing your source code. If you modify LibEphemeris itself, those modifications must be released under the LGPL. See [LICENSE](LICENSE).
