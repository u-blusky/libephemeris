# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

A high-precision astronomical ephemeris library for Python, powered by NASA JPL DE440/DE441 ephemerides and IAU 2006/2000A standards.

**Drop-in replacement for PySwissEph** - pure Python, no C extensions.

---

## Features

- **NASA JPL DE440/DE441** - the same ephemerides used for spacecraft navigation, via Skyfield
- **IAU 2006/2000A standards** - precession, nutation, and obliquity via the official ERFA library
- **Sub-arcsecond precision** - all planets < 0.75", house cusps < 0.02", verified across 4,400+ test rounds ([full report](docs/PRECISION.md))
- **25 house systems, 43 ayanamsha modes** - all independently verified against pyswisseph
- **Physical planet centers** - outer planets corrected from barycenters using JPL satellite ephemerides
- **4 backends** - Skyfield, LEB (~14x speedup), Horizons API, auto-fallback. All sub-arcsecond
- **15,000+ years** of coverage with extended tier (-13200 to +17191 CE)
- **Drop-in replacement for PySwissEph** - same API, same function names, minimal migration effort ([migration guide](docs/guides/migration-guide.md))
- **Pure Python 3.12+** - no C extensions, runs anywhere

---

## Why LibEphemeris

Swiss Ephemeris is the industry standard for planetary calculations. But its Python binding (pyswisseph) is a C extension - hard to build, hard to debug, tied to a single computation model.

LibEphemeris provides the **same API** with a modern foundation:

- **NASA JPL ephemerides** instead of semi-analytical theory - DE440/DE441 are the latest planetary ephemerides from the Jet Propulsion Laboratory, the same data used for spacecraft navigation.
- **IAU 2006/2000A standards** - precession and nutation computed via the official ERFA library (the open-source implementation of IAU SOFA), not custom routines.
- **Physical planet centers** - Jupiter, Saturn, Uranus, Neptune corrected from system barycenters to actual body centers using JPL satellite ephemerides. Most libraries skip this.
- **Independently verified** - every function cross-validated against pyswisseph, JPL Horizons, and astropy/ERFA. [Precision report with full methodology](docs/PRECISION.md).
- **Pure Python** - readable source, standard debugging, no build toolchain. Runs on any platform, any CI, any serverless environment.

**Switching from pyswisseph?** Your existing code works with minimal changes. [Migration guide](docs/guides/migration-guide.md).

---

## Quick Start

```python
import libephemeris as swe
from libephemeris.constants import SE_SUN, SE_MOON, SEFLG_SPEED

jd = swe.julday(2000, 1, 1, 12.0)  # J2000.0

sun, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SPEED)
moon, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SPEED)

print(f"Sun:  {sun[0]:.4f} deg, speed {sun[3]:.4f} deg/day")
print(f"Moon: {moon[0]:.4f} deg, speed {moon[3]:.4f} deg/day")
```

```python
# House cusps (Placidus, Rome)
cusps, ascmc = swe.houses(swe.julday(2024, 11, 5, 18.0), 41.9028, 12.4964, b"P")
print(f"ASC: {ascmc[0]:.4f}, MC: {ascmc[1]:.4f}")
```

More examples: [Getting Started](docs/guides/getting-started.md)

---

## Verified Precision

Every number independently measured against pyswisseph across 4,400+ comparison rounds. [Full report](docs/PRECISION.md).

| Category | Typical | Max | Scope |
|----------|---------|-----|-------|
| Planets (Sun-Pluto) | 0.04-0.26" | 0.75" | All planets, sub-arcsecond |
| Moon | 0.70" | 3.32" | Different underlying lunar models |
| House cusps | < 0.01" | 0.02" | All 25 systems |
| Fixed stars | < 0.1" | 0.51" | 116 Hipparcos catalog stars |
| Solar eclipses | - | < 6s | Timing accuracy |
| Lunar eclipses | - | < 8s | Timing accuracy |
| Ayanamsha | < 0.001 deg | 0.006 deg | All 43 sidereal modes |

---

## Four Backends, One API

Choose your trade-off between speed, precision, and infrastructure. All backends produce sub-arcsecond agreement through the same `calc_ut()` interface.

| Mode | Backend | Speed | Use case |
|------|---------|-------|----------|
| `"auto"` | LEB -> Horizons -> Skyfield | adaptive | **Default.** Works everywhere, zero config |
| `"skyfield"` | JPL DE440/DE441 via Skyfield | ~120 us | Full precision, offline |
| `"leb"` | Precomputed Chebyshev polynomials | ~5 us | Maximum throughput (14x faster) |
| `"horizons"` | NASA JPL Horizons REST API | ~300 ms | Zero local files, CI/serverless |

```python
from libephemeris import set_calc_mode
set_calc_mode("leb")  # or via env: LIBEPHEMERIS_MODE=leb
```

---

## Installation

```bash
pip install libephemeris
```

Works out of the box - the wheel includes bundled ephemeris data for 1850-2150 (14 bodies). For extended coverage:

```bash
libephemeris download medium       # 1550-2650, ~200 MB (recommended)
libephemeris download extended     # -13200 to +17191 CE, full range
libephemeris status                # Show installed data files
```

**Optional extras:** `pip install libephemeris[stars]` for fixed stars, `[nbody]` for n-body integration, `[all]` for everything. [Details](docs/guides/getting-started.md#optional-extras).

---

## Documentation

- [Getting Started](docs/guides/getting-started.md) - installation, ephemeris tiers, first calculations
- [Migration from PySwissEph](docs/guides/migration-guide.md) - API mapping, flag compatibility, known divergences
- [Precision Report](docs/PRECISION.md) - full methodology, comparison tables, verification process
- [Flag Reference](docs/reference/flags.md) - all supported flags with examples
- [House Systems](docs/reference/house-systems.md) - all 25 systems, verified against pyswisseph
- [Ayanamsha Modes](docs/reference/ayanamsha.md) - 43 sidereal modes
- [LEB Binary Ephemeris](docs/leb/guide.md) - format, generation, LEB2 compression
- [Horizons Backend](docs/architecture/horizons-backend.md) - HTTP client, pipeline, precision
- [Architecture](docs/development/architecture-overview.md) - internal design and data flow
- [Methodology](docs/methodology/overview.md) - planet centers, lunar apsides, pyerfa integration
- [CLI Reference](CLI.md) - full command reference
- [Changelog](CHANGELOG.md) - release history

---

## Contributing

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris && uv pip install -e ".[dev]"
poe test:unit:fast                 # Run unit tests
poe test:compare:skyfield          # Cross-validate vs pyswisseph
```

---

## License

AGPL-3.0. See [LICENSE](LICENSE).
