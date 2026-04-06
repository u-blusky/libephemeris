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

---

## Features

- **NASA JPL ephemerides** -- all calculations use DE440/DE441 (2021) via Skyfield, with Horizons API and precomputed LEB as alternative data sources.
- **Four calculation modes** -- `auto` (default), `skyfield`, `leb` (~14x speedup), `horizons` (zero-install). All produce sub-arcsecond agreement.
- **IAU standard transforms** -- nutation (IAU 2006/2000A), precession (IAU 2006), obliquity via the official ERFA library.
- **Physical planet centers** -- outer planets corrected from system barycenters to true body centers using JPL satellite ephemerides.
- **25 house systems** -- all systems supported by Swiss Ephemeris (26 codes including the A/E equal-house alias), independently verified against pyswisseph.
- **Sub-arcsecond precision** -- all planets < 1", house cusps < 0.02", independently verified against JPL Horizons and astropy/ERFA. [Full report](docs/PRECISION.md).
- **Pure Python** -- fully testable, readable, documented. No C extensions.

---

## Installation

```bash
pip install libephemeris
```

The PyPI wheel includes a bundled LEB2 base-tier core (~10.6 MB, 14 bodies, 1850-2150). With the default `medium` tier, LEB2 files are auto-downloaded on first use (~119 MB). For full offline coverage, download the complete data set for a precision tier:

```bash
libephemeris download medium       # DE440 + planet centers + minor-body SPKs (~200 MB total)
```

**Requirements:** Python 3.12+ | [Optional extras](docs/guides/getting-started.md#optional-extras): `[nbody]`, `[stars]`, `[all]`

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

## Calculation Modes

LibEphemeris supports four calculation modes with automatic fallback:

```
LEB (~5 us/eval)  -->  Horizons API (~300ms)  -->  Skyfield/DE440 (~120 us/eval)
```

| Mode | Behavior | Use case |
|------|----------|----------|
| `"auto"` **(default)** | LEB if configured, then Horizons if no DE440, then Skyfield | Works everywhere |
| `"leb"` | Require LEB (auto-discovered or auto-downloaded if needed); unsupported bodies/flags fall back to Skyfield | Maximum performance |
| `"horizons"` | Prefer Horizons API; unsupported bodies/flags fall back to Skyfield | Zero-install, serverless, CI |
| `"skyfield"` | Always Skyfield/DE440 | Offline, full precision |

```python
from libephemeris import set_calc_mode

set_calc_mode("horizons")  # Or via env: LIBEPHEMERIS_MODE=horizons
```

All modes produce the same positions (sub-arcsecond agreement). The default `"auto"` mode resolves data transparently via bundled LEB2, auto-download, Horizons API, or Skyfield.

| Data source | Details |
|-------------|---------|
| **Skyfield** | Pure JPL DE440/DE441 via Skyfield. Gold standard. [Ephemeris tiers](docs/guides/getting-started.md#ephemeris-tiers) |
| **LEB** | Precomputed Chebyshev polynomials, ~14x speedup. [LEB Guide](docs/leb/guide.md) |
| **LEB2** | Compressed LEB format (4-10x smaller). Base-tier core (~10.6 MB) bundled in wheel; other tiers auto-downloaded. [LEB2 details](docs/leb/guide.md#13-leb2-compressed-format) |
| **Horizons** | NASA JPL Horizons REST API. No local files needed. [Horizons Guide](docs/architecture/horizons-backend.md) |

---

## Precision

Measured across 4,400+ comparison rounds. [Full precision report](docs/PRECISION.md).

| Category | Typical | Max | Notes |
|----------|---------|-----|-------|
| Planets (Sun-Pluto) | 0.04-0.26" | 0.75" | Sub-arcsecond, all planets |
| Moon | 0.70" | 3.32" | Different lunar models |
| House cusps | < 0.01" | 0.02" | All 25 systems tested |
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

Full flag reference with examples: [docs/reference/flags.md](docs/reference/flags.md)

---

## Documentation

### Guides
- [Getting Started](docs/guides/getting-started.md) -- installation, ephemeris tiers, first calculations
- [Migration from PySwissEph](docs/guides/migration-guide.md)
- [Optional Modules](docs/guides/optional-modules.md) -- backends, fallback chain, optional extras
- [Precision Tuning](docs/guides/precision-tuning.md)
- [Computation Tracing](docs/guides/tracing.md) -- discover which backend computed each body

### Architecture
- [LEB Binary Ephemeris](docs/leb/guide.md) -- format, reader, generation, fast-path pipeline
- [Horizons API Backend](docs/architecture/horizons-backend.md) -- HTTP client, pipeline, modes, precision
- [Architecture Overview](docs/development/architecture-overview.md)

### Reference
- [Precision Report](docs/PRECISION.md)
- [Flag Reference](docs/reference/flags.md)
- [Known Divergences](docs/reference/divergences.md)
- [House Systems](docs/reference/house-systems.md)
- [Ayanamsha Modes](docs/reference/ayanamsha.md)

### Methodology
- [Overview](docs/methodology/overview.md)
- [Planet Centers](docs/methodology/planet-centers-spk.md)
- [Lunar Apsides](docs/methodology/lunar-apsides.md)
- [True Lilith](docs/methodology/true-lilith.md)
- [pyerfa Integration](docs/methodology/pyerfa-integration.md)

### Manuals
- [Manuale (IT)](docs/manual/it/) -- introduzione ai calcoli astrologici per principianti, 15 capitoli
- [Manual (EN)](docs/manual/en/) -- beginner's guide to astrological calculations, 15 chapters

### Development
- [CLI Reference](CLI.md) -- full command reference for `libephemeris` and `leph`
- [Roadmap](docs/development/roadmap.md)
- [Changelog](CHANGELOG.md)

---

## CLI Commands

```bash
# Download ephemeris data (DE kernel + planet centers + minor-body SPKs)
libephemeris download base          # 1850-2150 (lightweight)
libephemeris download medium        # 1550-2650 (recommended)
libephemeris download extended      # -13200 to +17191 (full range)

# Status and version
libephemeris status                 # Show installed data files
libephemeris --version
```

---

## Development

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris && uv pip install -e ".[dev]"
```

Key commands ([full list](CLI.md)):

| Command | Description |
|---------|-------------|
| `poe test:skyfield:fast` | Skyfield unit tests, parallel (~1 min) |
| `poe test:leb:fast` | LEB unit tests, parallel (~1 min) |
| `poe test:horizons:core` | Horizons API precision (~15s) |
| `poe test:compare:skyfield` | Cross-validate vs pyswisseph |
| `poe test:houses` | All 25 house systems vs pyswisseph |
| `poe lint` | Ruff linter |
| `poe format` | Ruff formatter |

---

## Performance

### `reset_session()` -- lightweight state reset

Resets only per-calculation state (topo, sidereal mode, angles cache) without closing
file handles or clearing LRU caches. Use between independent calculations to avoid the
full teardown cost of `close()`.

```python
import libephemeris as swe

swe.calc_ut(jd1, swe.SE_SUN, flags)  # First calculation
swe.reset_session()                    # Reset topo/sidereal, keep reader alive
swe.calc_ut(jd2, swe.SE_SUN, flags)  # Reuses LEB reader, timescale, caches
```

**Impact**: Consecutive calculations drop from ~3500ms to ~2ms (1750x speedup).

### LEB2 v2 chunked format

LEB2 files use 10-year temporal chunks instead of monolithic per-body compression.
Only the chunk containing the requested Julian Day is decompressed (~300 KB for Moon
instead of 307 MB). The reader transparently supports both v1 (legacy) and v2 (chunked).

**Impact**: Cold-start decompression drops from 1568ms to 47ms (33x speedup).

---

## License

AGPL-3.0. See [LICENSE](LICENSE).
