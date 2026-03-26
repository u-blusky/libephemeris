# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

A high-precision astronomical ephemeris library for Python, powered by NASA JPL data and modern IAU standards. API-compatible drop-in replacement for PySwissEph.

> [!NOTE]
> **Alpha 1.0.0** -- API is stabilizing. Breaking changes will be documented in the [changelog](CHANGELOG.md).

---

## Features

- **NASA JPL ephemerides** -- all calculations use DE440/DE441 (2021) exclusively via Skyfield. No analytical fallbacks.
- **Three calculation backends** -- Skyfield (pure JPL), LEB (precomputed Chebyshev, ~14x speedup), and Horizons API (zero-install, no local files needed).
- **IAU standard transforms** -- nutation (IAU 2006/2000A), precession (IAU 2006), obliquity via the official ERFA library.
- **Physical planet centers** -- outer planets corrected from system barycenters to true body centers using JPL satellite ephemerides.
- **Sub-arcsecond precision** -- all planets < 1", house cusps < 0.02", independently verified against JPL Horizons and astropy/ERFA. [Full report](docs/PRECISION.md).
- **Pure Python** -- fully testable, readable, documented. No C extensions.

---

## Installation

```bash
pip install libephemeris
```

Works immediately -- the library fetches data from NASA JPL Horizons API when no local files are present. For offline use or faster calculations, download ephemeris files:

```bash
libephemeris download:medium      # JPL DE440 (~128 MB), recommended
```

**Requirements:** Python 3.9+ | [Optional extras](docs/getting-started.md#optional-extras): `[nbody]`, `[stars]`, `[all]`

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

More examples: [Getting Started](docs/getting-started.md)

---

## Calculation Backends

LibEphemeris supports three backends with automatic fallback:

```
LEB (~5 us/eval)  -->  Horizons API (~300ms)  -->  Skyfield/DE440 (~120 us/eval)
```

| Mode | Behavior | Use case |
|------|----------|----------|
| `"auto"` **(default)** | LEB if configured, then Horizons if no DE440, then Skyfield | Works everywhere |
| `"leb"` | LEB precomputed only | Maximum performance |
| `"horizons"` | NASA JPL Horizons API only | Zero-install, serverless, CI |
| `"skyfield"` | Skyfield/DE440 only | Offline, maximum precision |

```python
from libephemeris import set_calc_mode

set_calc_mode("horizons")  # Or via env: LIBEPHEMERIS_MODE=horizons
```

All backends produce the same positions (sub-arcsecond agreement). The default `"auto"` mode always works -- with or without local data files.

| Backend | Details |
|---------|---------|
| **Skyfield** | Pure JPL DE440/DE441 via Skyfield. Gold standard. [Ephemeris tiers](docs/getting-started.md#ephemeris-tiers) |
| **LEB** | Precomputed Chebyshev polynomials, ~14x speedup. [LEB Guide](docs/leb/guide.md) |
| **LEB2** | Compressed LEB format (4-10x smaller). Core bodies (~8.7 MB) ship in the PyPI wheel. [LEB2 details](docs/leb/guide.md#13-leb2-compressed-format) |
| **Horizons** | NASA JPL Horizons REST API. No local files needed. [Horizons Guide](docs/horizons-backend.md) |

---

## Precision

Measured across 4,400+ comparison rounds. [Full precision report](docs/PRECISION.md).

| Category | Typical | Max | Notes |
|----------|---------|-----|-------|
| Planets (Sun-Pluto) | 0.04-0.26" | 0.75" | Sub-arcsecond, all planets |
| Moon | 0.70" | 3.32" | Different lunar models |
| House cusps | < 0.01" | 0.02" | All 24 systems tested |
| Fixed stars | < 0.1" | 0.51" | 116 Hipparcos stars |
| Solar eclipses | -- | < 6s | Timing precision |
| Lunar eclipses | -- | < 8s | Timing precision |
| Ayanamsha | < 0.001 deg | 0.006 deg | 43 sidereal modes |

---

## Flags Reference

Flags control what is calculated and how results are returned. Combine with `|`.

| Flag | Effect |
|------|--------|
| `SEFLG_SPEED` | Populate speed fields (almost always needed) |
| `SEFLG_HELCTR` | Heliocentric observer |
| `SEFLG_TOPOCTR` | Topocentric observer (set position with `set_topo()`) |
| `SEFLG_EQUATORIAL` | Output RA/Dec instead of ecliptic lon/lat |
| `SEFLG_J2000` | J2000.0 frame instead of equinox of date |
| `SEFLG_SIDEREAL` | Sidereal zodiac (requires `set_sid_mode()`) |
| `SEFLG_TRUEPOS` | Geometric position (no light-time/aberration) |
| `SEFLG_NOABERR` | Astrometric position (no aberration) |

Full flag reference with examples: [docs/flags.md](docs/flags.md)

---

## Documentation

### Guides
- [Getting Started](docs/getting-started.md) -- installation, ephemeris tiers, first calculations
- [Migration from PySwissEph](docs/guides/migration-guide.md)
- [Precision Tuning](docs/guides/precision-tuning.md)

### Architecture
- [LEB Binary Ephemeris](docs/leb/guide.md) -- format, reader, generation, fast-path pipeline
- [Horizons API Backend](docs/horizons-backend.md) -- HTTP client, pipeline, modes, precision
- [Architecture Overview](docs/development/architecture-overview.md)

### Reference
- [Precision Report](docs/PRECISION.md) -- full measurement data
- [Flag Reference](docs/flags.md)
- [House Systems](docs/reference/house-systems.md)
- [Ayanamsha Modes](docs/reference/ayanamsha.md)

### Methodology
- [Overview](docs/methodology/overview.md)
- [Planet Centers](docs/methodology/planet-centers-spk.md)
- [Lunar Apsides](docs/methodology/lunar-apsides.md)
- [True Lilith](docs/methodology/true-lilith.md)
- [pyerfa Integration](docs/methodology/pyerfa-integration.md)

### Development
- [Testing](docs/development/testing.md) -- test suites, commands, coverage
- [Roadmap](docs/development/roadmap.md)
- [Changelog](CHANGELOG.md)

---

## CLI Commands

```bash
# Download ephemeris data
libephemeris download:base         # DE440s (~35 MB), 1850-2150
libephemeris download:medium       # DE440 (~130 MB), 1550-2650 (recommended)
libephemeris download:extended     # DE441 (~3.3 GB), -13200 to +17191

# Status and version
libephemeris status                # Show installed data files
libephemeris --version
```

---

## Development

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris && uv pip install -e ".[dev]"
```

Key commands ([full list](docs/development/testing.md)):

| Command | Description |
|---------|-------------|
| `poe test:unit:fast` | Unit tests, parallel (~1 min) |
| `poe test:leb2:precision:all` | LEB2 compression precision (all tiers) |
| `poe test:horizons:quick` | Horizons API precision (~15s) |
| `poe test:compare` | Cross-validate vs reference |
| `poe lint` | Ruff linter |
| `poe format` | Ruff formatter |

---

## License

AGPL-3.0. See [LICENSE](LICENSE).
