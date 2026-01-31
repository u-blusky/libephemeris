# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>


**High-precision astronomical ephemeris library for Python (Swiss Ephemeris compatible, powered by Skyfield and JPL DE ephemerides).**

> [!WARNING] > **Pre-Alpha Release**
>
> LibEphemeris is currently in an early pre-alpha stage. The public API is unstable and may change without notice. Do not rely on it in production yet.

## Project Philosophy & Roadmap

LibEphemeris is born from the need for a **truly open, maintainable, and scientifically rigorous** alternative to the Swiss Ephemeris. While Swiss Ephemeris is the industry standard, its codebase is complex, difficult to maintain, and rooted in older C practices.

**Our Goals:**

1.  **True Open Source**: A codebase that is easy to read, contribute to, and maintain.
2.  **Scientific Precision**: Leveraging modern astronomical data from NASA JPL (via Skyfield/Starfield) to guarantee precision that matches or exceeds current standards.
3.  **Independence**: A completely original implementation, not just a wrapper around the existing C library.

### Roadmap

-   **Milestone 1 (Current)**: **Pure Python Library**

    -   A 1:1 drop-in replacement for `pyswisseph`.
    -   Powered by [Skyfield](https://rhodesmill.org/skyfield/) and NASA JPL DE421+ ephemerides.
    -   Focus on correctness, readability, and higher scientific precision than Swiss Ephemeris implementations.

-   **Milestone 2 (Next Step)**: **Rust Core Rewrite**
    -   Porting the core logic to **Rust** for maximum performance, memory safety, and stability.
    -   Will utilize [Starfield](https://docs.rs/starfield/latest/starfield/) (a Rust port of Skyfield) instead of the Python Skyfield library.
    -   This will provide a high-performance backend while maintaining the easy-to-use Python interface.

---

## Features at a Glance

-   **Planetary positions**: Sun, Moon, all major planets and Pluto.
-   **High precision**: Based on NASA JPL DE440 by default (configurable to other DE files).
-   **Multiple coordinate systems**: Ecliptic, equatorial, J2000 and of-date frames.
-   **Observation modes**: Geocentric, topocentric, heliocentric, barycentric.
-   **Velocities**: Full 6-component state vectors (position + velocity).
-   **House systems (19)**: Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign, Porphyry, Alcabitius, Polich/Page (Topocentric), Morinus, Meridian, Vehlow, Horizontal, Carter, Krusinski, Natural Gradient, and more.
-   **Sidereal zodiac (43 ayanamshas)**: Fagan/Bradley, Lahiri, Raman, Krishnamurti, star-based and historical variants.
-   **Extended points**: Lunar nodes, Lilith (mean and true), interpolated apogee/perigee, major asteroids (Chiron, Pholus, Ceres, Pallas, Juno, Vesta), centaurs (Nessus, Asbolus, Chariklo), TNOs (Orcus, Ixion, Haumea, Quaoar, Makemake, Gonggong, Eris, Sedna), major fixed stars and Arabic parts.
-   **High-precision minor bodies via SPK**: Download SPK kernels from JPL Horizons for arcsecond-level precision on asteroids and TNOs.
-   **Event finding**: Sun/Moon longitude crossings (e.g. ingress), with additional events planned (eclipses, etc.).
-   **Thread safety**: Optional thread-safe `EphemerisContext` API for concurrent calculations.
-   **Swiss Ephemeris compatibility**: Same function names, flags and result structure as `pyswisseph` in most common use cases.

---

## Installation

Using `pip`:

```bash
pip install libephemeris
```

Using [`uv`](https://github.com/astral-sh/uv) (recommended for development):

```bash
uv pip install libephemeris
```

From source:

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris
uv pip install -e .
```

### Requirements

-   Python **3.9+**
-   `skyfield>=1.54`
-   `skyfield-data>=7.0.0`
-   A JPL ephemeris file (DE440 by default, downloaded automatically on first use if not present locally)

### Optional Dependencies

LibEphemeris has several optional dependencies that enhance functionality when installed:

| Package | Install Command | Features Enabled |
|---------|-----------------|------------------|
| **pyerfa** | `pip install pyerfa` | High-precision IAU 2006 precession and obliquity for true lunar node calculations. Falls back to Lieske (1979) formulas when not installed. |
| **astroquery** | `pip install astroquery` | Automatic SPK kernel downloads from JPL Horizons for arcsecond-level precision on asteroids and TNOs. Required for `set_auto_spk_download(True)` and `download_spk()`. |
| **astropy** | `pip install astropy` | Required by the star catalog build script (`scripts/build_star_catalog.py`) for unit handling when querying Hipparcos data. Not needed for normal library usage. |

**Installation with optional dependencies:**

```bash
# Install with high-precision lunar node support
pip install libephemeris[precision]

# Install with automatic SPK download support
pip install libephemeris[spk]

# Install with all optional features
pip install libephemeris[all]

# Install with star catalog building support (for development)
pip install libephemeris[stars]
```

---

## Quick Start

### Basic planetary positions

```python
import libephemeris as ephem
from libephemeris.constants import *

# Julian Day (J2000.0)
jd = ephem.swe_julday(2000, 1, 1, 12.0)

# Sun position (longitude, latitude, distance, and speeds)
sun, flags = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
print(f"Sun longitude: {sun[0]:.6f}°")
print(f"Sun speed: {sun[3]:.6f}°/day")

# Moon position
moon, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH)
print(f"Moon longitude: {moon[0]:.6f}°")
```

### Houses and angles

```python
# Rome coordinates
lat, lon, alt = 41.9028, 12.4964, 0.0
jd = ephem.swe_julday(2024, 11, 5, 18.0)

# Placidus houses
cusps, ascmc = ephem.swe_houses(jd, lat, lon, b"P")

print(f"Ascendant: {ascmc[0]:.2f}°")
print(f"MC:        {ascmc[1]:.2f}°")
print(f"House 1:   {cusps[1]:.2f}°")
print(f"House 10:  {cusps[10]:.2f}°")
```

### Sidereal calculations

```python
# Lahiri ayanamsha
ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

ayanamsha = ephem.swe_get_ayanamsa_ut(jd)
print(f"Lahiri ayanamsha: {ayanamsha:.6f}°")

sun_sid, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)
print(f"Sidereal Sun: {sun_sid[0]:.6f}°")
```

### Longitude crossings

```python
# Next time the Sun enters 0° Aries
next_cross = ephem.swe_solcross_ut(0.0, jd, SEFLG_SWIEPH)
print(f"Next Aries ingress (JD): {next_cross:.6f}")
```

---

## Thread Safety

LibEphemeris provides **two APIs** for different use cases:

### Global API (pyswisseph-compatible, NOT thread-safe)

The traditional Swiss Ephemeris API using global state:

```python
import libephemeris as swe

swe.set_topo(12.5, 41.9, 0)  # Rome
pos, _ = swe.calc_ut(2451545.0, swe.SE_SUN, 0)
```

✅ **100% pyswisseph compatible**  
⚠️ **NOT thread-safe** (matches Swiss Ephemeris behavior)

### Context API (thread-safe)

For multi-threaded applications, use `EphemerisContext`:

```python
from libephemeris import EphemerisContext, SE_SUN, SE_MOON

# Each thread creates its own context
ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)  # Rome

sun, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
moon, _ = ctx.calc_ut(2451545.0, SE_MOON, 0)
```

✅ **Thread-safe**  
✅ **Independent state per context**  
✅ **Shared ephemeris files** (memory efficient)

#### Multi-threading example

```python
from libephemeris import EphemerisContext, SE_SUN
from concurrent.futures import ThreadPoolExecutor

def calculate_chart(location, jd):
    """Thread-safe calculation function."""
    ctx = EphemerisContext()
    ctx.set_topo(location['lon'], location['lat'], 0)

    sun, _ = ctx.calc_ut(jd, SE_SUN, 0)
    return {'location': location['name'], 'sun_lon': sun[0]}

locations = [
    {'name': 'Rome', 'lon': 12.5, 'lat': 41.9},
    {'name': 'London', 'lon': -0.1, 'lat': 51.5},
    {'name': 'Tokyo', 'lon': 139.7, 'lat': 35.7},
]

# Calculate concurrently
with ThreadPoolExecutor(max_workers=3) as executor:
    results = list(executor.map(
        lambda loc: calculate_chart(loc, 2451545.0),
        locations
    ))

for r in results:
    print(f"{r['location']}: Sun at {r['sun_lon']:.2f}°")
```

---

## Configuring Ephemeris Files

By default, LibEphemeris uses the **JPL DE440** kernel (`de440.bsp`), which covers roughly **1550–2650**. The file is:

-   loaded from a local path if already present, or
-   automatically downloaded via Skyfield the first time it is needed.

You can control which ephemeris file is used and where it is loaded from.

### Choosing a different JPL kernel

Use `set_ephemeris_file()` to select a different `.bsp` file:

```python
from libephemeris import set_ephemeris_file

set_ephemeris_file("de441.bsp")  # very long time span, larger file
```

Supported JPL kernels include, for example:

-   `de421.bsp`: 1900–2050 (legacy, ~16 MB)
-   `de422.bsp`: −3000–3000 (~623 MB)
-   `de430.bsp`: 1550–2650 (~128 MB)
-   `de431.bsp`: −13200–17191 (~3.4 GB)
-   `de440.bsp`: 1550–2650 (default, ~128 MB) - recommended for most uses
-   `de441.bsp`: −13200–17191 (~3.4 GB) - for extended historical/future work

If the chosen file is not present locally, Skyfield will attempt to download it.

### Custom ephemeris directory

Use `set_ephe_path()` to point LibEphemeris to a directory containing JPL kernels:

```python
from libephemeris import set_ephe_path

set_ephe_path("/path/to/jpl-kernels")
```

Resolution order for the ephemeris file is:

1. The directory set via `set_ephe_path()`, if any.
2. The project/workspace root.
3. Download via Skyfield.

If you try to compute positions outside the date range covered by the selected kernel, LibEphemeris will raise an exception describing the supported range.

---

## High-Precision Minor Bodies with SPK Kernels

By default, minor bodies (asteroids, TNOs like Chiron, Eris, Sedna) are calculated using **Keplerian orbital elements** with secular perturbations. This provides arcminute-level accuracy, which may not be sufficient for precise astrological or astronomical work.

For **arcsecond-level precision**, you can download **SPK kernels** directly from NASA JPL Horizons and use them for calculations.

### Quick example

```python
import libephemeris as eph

# Download SPK for Chiron covering 2000-2100
spk_path = eph.download_spk(
    body="Chiron",
    start="2000-01-01",
    end="2100-01-01",
    directory="./spk_kernels"
)

# Register the SPK for Chiron calculations
eph.register_spk_body(
    ipl=eph.SE_CHIRON,
    spk_path=spk_path,
    naif_id=eph.NAIF_CHIRON  # 2002060
)

# Now calc_ut automatically uses the SPK kernel
pos, _ = eph.calc_ut(2451545.0, eph.SE_CHIRON, eph.SEFLG_SPEED)
print(f"Chiron (SPK): {pos[0]:.6f}°")
```

### Convenience function

For a one-liner approach:

```python
import libephemeris as eph

# Download and register in one call
eph.download_and_register_spk(
    body="Eris",
    ipl=eph.SE_ERIS,
    start="2000-01-01",
    end="2100-01-01",
)

# Calculations now use SPK automatically
pos, _ = eph.calc_ut(2460000.0, eph.SE_ERIS, 0)
```

### Available functions

| Function | Description |
|----------|-------------|
| `download_spk(body, start, end, ...)` | Download SPK from JPL Horizons |
| `register_spk_body(ipl, spk_path, naif_id)` | Register SPK for a body |
| `unregister_spk_body(ipl)` | Remove SPK registration |
| `download_and_register_spk(...)` | Download and register in one call |
| `list_spk_bodies()` | List all registered SPK bodies |
| `get_spk_body_info(ipl)` | Get SPK info for a body |
| `get_spk_coverage(spk_path)` | Get date range covered by SPK |
| `set_auto_spk_download(enable)` | Enable/disable automatic SPK download |
| `set_spk_cache_dir(path)` | Set SPK cache directory |
| `set_spk_date_padding(days)` | Set date range padding for auto-downloads |

### Automatic SPK Download (Experimental)

LibEphemeris can automatically download and cache SPK kernels on demand:

```python
import libephemeris as eph

# Enable automatic SPK download
eph.set_auto_spk_download(True)
eph.set_spk_cache_dir("./spk_cache")  # Optional: custom cache location

# First calculation triggers automatic download
pos, _ = eph.calc_ut(2460000.0, eph.SE_CHIRON, 0)  # Downloads SPK if needed
```

### NAIF ID constants

SPK kernels use NASA NAIF IDs. LibEphemeris provides constants for common bodies:

```python
NAIF_ASTEROID_OFFSET = 2000000  # asteroid_number + offset = NAIF ID
NAIF_CHIRON = 2002060           # Chiron (2060)
NAIF_CERES = 2000001            # Ceres (1)
NAIF_PALLAS = 2000002           # Pallas (2)
NAIF_JUNO = 2000003             # Juno (3)
NAIF_VESTA = 2000004            # Vesta (4)
NAIF_PHOLUS = 2005145           # Pholus (5145)
NAIF_NESSUS = 2007066           # Nessus (7066)
NAIF_ASBOLUS = 2008405          # Asbolus (8405)
NAIF_CHARIKLO = 2010199         # Chariklo (10199)
NAIF_ERIS = 2136199             # Eris (136199)
NAIF_SEDNA = 2090377            # Sedna (90377)
NAIF_GONGGONG = 2225088         # Gonggong (225088)
```

For any numbered asteroid, the NAIF ID is `asteroid_number + 2000000`.

### Thread-safe usage with EphemerisContext

```python
from libephemeris import EphemerisContext, SE_CHIRON

ctx = EphemerisContext()

# Register SPK in this context only
ctx.register_spk_body(SE_CHIRON, "./chiron.bsp", 2002060)

pos, _ = ctx.calc_ut(2451545.0, SE_CHIRON, 0)
```

### Precision comparison

| Method | Typical Accuracy | Use Case |
|--------|------------------|----------|
| Keplerian (default) | ~1-10 arcminutes | Quick estimates, historical charts |
| SPK kernel | ~1-5 arcseconds | Precise calculations, transit timing |

### Notes

- SPK files can be large (1-50 MB depending on date range)
- Files are cached locally after download
- If a date is outside SPK coverage, an exception is raised (no silent fallback)
- SPK kernels are downloaded from JPL Horizons API (requires internet for download)

---

## Scientific Accuracy and Validation

### Ephemeris data

-   **Source**: NASA JPL DE ephemerides (DE440 by default).
-   **Time span**: 1550–2650 for DE440; extendable by selecting other kernels.
-   **Precision**: Sub-arcsecond accuracy for major planets within the supported range.
-   **Reference frame**: ICRS/J2000.0 (ICRF 3.0 for DE440/DE441).

### Comparison with Swiss Ephemeris

LibEphemeris has been tested against Swiss Ephemeris using an automated test suite.

| Component            | Tests | Pass Rate | Max Difference |
| -------------------- | ----- | --------- | -------------- |
| Planetary positions  | 229   | 100%      | < 0.001°       |
| House systems        | 113   | 100%      | < 0.001°       |
| Ayanamsha values     | 129   | 100%      | < 0.06°        |
| Lunar nodes / Lilith | 40+   | 100%      | < 0.01°        |
| Velocities           | 100   | 100%      | < 0.01°/day    |

These comparisons are implemented in the `tests/` and `compare_scripts/` directories, and are run regularly during development.

### Performance Benchmarks

LibEphemeris is a pure Python implementation, while pyswisseph is a C library. The performance difference reflects this architectural choice, prioritizing readability and maintainability over raw speed. Future Rust core implementation (Milestone 2) will significantly improve performance.

| Operation       | Iterations | pyswisseph (C)       | libephemeris (Python) | Slowdown |
| --------------- | ---------- | -------------------- | --------------------- | -------- |
| calc_ut         | 5000       | 0.04s (~138,000/s)   | 4.2s (~1,180/s)       | ~117x    |
| houses          | 1200       | 0.004s (~296,000/s)  | 0.13s (~9,400/s)      | ~31x     |
| get_ayanamsa_ut | 1500       | 0.001s (~1,300,000/s)| 0.03s (~45,000/s)     | ~29x     |

*Benchmarks run on Apple M1 Pro, Python 3.10. Results may vary depending on hardware and system load.*

**Key observations:**

-   **House calculations** are the fastest operation, with only ~31x slowdown.
-   **Ayanamsa calculations** are similarly efficient at ~29x slowdown.
-   **Planetary positions** (calc_ut) are the most computationally intensive, as they involve JPL ephemeris file lookups via Skyfield.

Run the benchmarks yourself:

```bash
pytest tests/test_precision/test_benchmark_comparison.py -v -s
```

---

## Swiss Ephemeris Compatibility

LibEphemeris aims to behave as a **drop-in replacement** for `pyswisseph` in many scenarios:

-   Same function names (e.g. `swe_calc_ut`, `swe_houses`, `swe_julday`, `swe_revjul`, `swe_get_ayanamsa_ut`).
-   Same integer constants and flags from `libephemeris.constants` (e.g. `SE_SUN`, `SEFLG_SWIEPH`, `SEFLG_SPEED`, `SE_SIDM_LAHIRI`).
-   Similar return types and value ordering.

There are still differences and missing features compared to the full Swiss Ephemeris API, especially around:

-   very long time ranges (beyond the chosen JPL kernel),
-   eclipse and occultation functions,
-   the full minor-planet and fixed-star catalogues.

Please open an issue if you hit a compatibility gap that is important for your use case.

---

## Development

### Project layout

```text
libephemeris/
├── libephemeris/
│   ├── __init__.py
│   ├── constants.py      # Constants and flags
│   ├── context.py        # Thread-safe EphemerisContext
│   ├── planets.py        # Planetary calculations
│   ├── houses.py         # House systems
│   ├── lunar.py          # Nodes and Lilith
│   ├── minor_bodies.py   # Asteroids and TNOs (Keplerian)
│   ├── spk.py            # SPK kernel support for high-precision minor bodies
│   ├── fixed_stars.py    # Fixed stars and points
│   ├── crossing.py       # Longitude crossing events
│   ├── angles.py         # Angle helpers (Asc, MC, etc.)
│   ├── arabic_parts.py   # Arabic parts calculations
│   ├── time_utils.py     # Time conversion helpers
│   └── state.py          # Global state (loader, ephemeris, sidereal mode)
├── tests/                # Comprehensive test suite
├── compare_scripts/      # Swiss Ephemeris comparison tools
└── README.md
```

### Development workflow

Install in editable mode with development dependencies:

```bash
uv pip install -e ".[dev]"
```

Run tests (via `poethepoet` tasks defined in `pyproject.toml`):

```bash
poe test       # run pytest
poe coverage   # run tests with coverage
poe lint       # run Ruff (lint)
poe format     # run Ruff formatter
```

---

## License

LibEphemeris is licensed under the **GNU Lesser General Public License v3.0 (LGPL-3.0)**.

See `LICENSE` for the full text.

**What this means for you:**

1.  **Commercial Use**: You **can** use this library in proprietary/commercial software without releasing your source code, provided you link to the library dynamically (e.g., as a normal Python dependency installed via `pip`).
2.  **Modifications**: If you modify the source code of _LibEphemeris itself_, you **must** release those modifications under the LGPL.
