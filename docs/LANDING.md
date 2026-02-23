# Landing Page Content

> This document contains the structured content for the LibEphemeris landing page.
> Each section is marked with its role in the page layout. Code blocks, tables,
> and copy are ready to be dropped into the final implementation.

---

## [Header]

- **Logo text:** LibEphemeris
- **Nav links:** Docs | API Reference | Examples | GitHub
- **CTA (pill/badge):** `pip install libephemeris`

---

## [Hero]

### Headline

Astronomical ephemeris calculations for modern Python.

### Subheadline

High-precision planetary positions, house systems, eclipses, and more --
computed directly from NASA JPL DE440/DE441 ephemerides, the same data used
by NASA for interplanetary mission planning. API-compatible with Swiss
Ephemeris. Pure Python, type-safe, and thread-safe.

### Primary CTA

Get Started _(links to documentation)_

### Secondary CTA

View on GitHub _(links to https://github.com/g-battaglia/libephemeris)_

### Install snippet (terminal block)

```bash
pip install libephemeris
```

---

## [Code Comparison]

> Side-by-side blocks showing the drop-in replacement. Visually identical code
> with only the import changed.

### Before (pyswisseph)

```python
import swisseph as swe

swe.set_ephe_path("/usr/share/ephe")
jd = swe.julday(2000, 1, 1, 12.0)
pos, flag = swe.calc_ut(jd, swe.SUN, swe.FLG_SPEED)
cusps, ascmc = swe.houses(jd, 41.9, 12.5, b"P")
```

### After (libephemeris)

```python
import libephemeris as swe

# No ephemeris path needed - JPL kernels managed automatically
jd = swe.julday(2000, 1, 1, 12.0)
pos, flag = swe.calc_ut(jd, swe.SUN, swe.FLG_SPEED)
cusps, ascmc = swe.houses(jd, 41.9, 12.5, b"P")
```

### Caption

Same API. Modern engine. Change the import, keep your code.

---

## [Why LibEphemeris] (4 cards)

### Card 1 — Powered by NASA JPL Data

Every calculation is grounded in the JPL DE440/DE441 planetary ephemerides --
the gold standard in positional astronomy, developed by NASA's Jet Propulsion
Laboratory for interplanetary navigation. No analytical approximations, no
fallback to simplified theories. The same data NASA trusts for mission-critical
trajectory planning.

### Card 2 — Scientific & Astrological Precision

Outer planet positions are automatically corrected from system barycenters to
physical planet centers using dedicated JPL SPK kernels. Lunar apsides are
derived from actual physical distance extrema in JPL data, not truncated
analytical series. Every position reflects where the body truly is in space.

### Card 3 — Current IAU Standards

Precession and nutation follow IAU 2006/2000A models via pyerfa. Delta T uses
Stephenson, Morrison & Hohenkerk (2016) with optional IERS observed data --
the current astronomical standard for historical timescale conversion. Not
locked to models from the 1980s.

### Card 4 — Drop-in Swiss Ephemeris Replacement

Full 1:1 API compatibility with pyswisseph. All `swe_*` functions, `SE_*`
constants, and calculation flags work identically. Migrate existing code by
changing a single import line -- your calculations gain NASA JPL precision
with zero refactoring.

---

## [Features] (grid)

> Each item has a short title and a one-line description. Display as a grid
> (3 columns suggested) or a compact list.

### Primary features

- **Planetary Positions** --
  All 10 planets with geocentric, heliocentric, and topocentric coordinates.
  Longitude, latitude, distance, and daily velocities.

- **19 House Systems** --
  Placidus, Koch, Whole Sign, Equal, Regiomontanus, Campanus, Porphyry,
  Alcabitius, Morinus, and 10 more. Includes polar latitude fallback handling.

- **43 Ayanamsha Modes** --
  Full sidereal calculation support. Lahiri, Fagan-Bradley, Krishnamurti,
  True Chitrapaksha, and 39 others.

- **Solar & Lunar Eclipses** --
  Global and local search, Besselian elements, contact points, path geometry,
  central line, magnitude, obscuration, and Saros/Inex series identification.

- **100+ Fixed Stars** --
  Hipparcos-based catalog with rigorous space motion propagation. Royal stars,
  Pleiades, zodiacal constellations, and more.

- **35+ Minor Bodies** --
  Main-belt asteroids, centaurs (Chiron, Pholus, Nessus), and TNOs (Eris,
  Sedna, Makemake). Keplerian model with perturbation corrections, or
  sub-arcsecond precision via automatic SPK kernel downloads.

- **Lunar Nodes & Lilith** --
  Mean and osculating ascending node. Mean and true Black Moon Lilith.
  Interpolated apogee and perigee calibrated against DE441 physical
  distance extrema.

- **Rise, Set & Transit** --
  Rise and set times, meridian transit, with atmospheric refraction correction.

- **15 Hypothetical Bodies** --
  Eight Hamburg School Uranian planets, Transpluto, Vulcan, White Moon
  (Selena), and more, computed from bundled orbital elements.

### Extended features

- **Heliacal Events** --
  Heliacal rising and setting with Schaefer atmospheric extinction model and
  limiting magnitude calculations.

- **21 Planetary Moons** --
  Galilean moons, Titan, Triton, Phobos, Deimos, Charon, and others via SPK
  kernels and analytical satellite theories.

- **Crossing Events** --
  Sun, Moon, and planet longitude crossings (ingresses) with Newton-Raphson
  and Brent fallback convergence.

- **Occultations** --
  Lunar and planetary occultation search, with local and global circumstances.

- **N-body Integration** --
  Optional REBOUND/ASSIST support for n-body orbit propagation beyond
  Keplerian approximation.

- **IERS Data** --
  Automatic download of observed Delta T, UT1-UTC, and leap second data from
  IERS for sub-second timescale accuracy.

---

## [Precision Tiers]

### Headline

From modern horoscopes to archeoastronomy. Choose your range.

### Subheadline

Three precision tiers map to different JPL kernels. DE440 and DE441 have
identical precision -- DE441 is simply the extended-range version.

### Tier table

| Tier | Kernel | Date range | Size | Use case |
| --- | --- | --- | --- | --- |
| Base | `de440s.bsp` | 1849 -- 2150 | ~31 MB | Lightweight web apps, contemporary calculations |
| **Medium** (default) | `de440.bsp` | 1550 -- 2650 | ~128 MB | General purpose, historical astrology |
| Extended | `de441.bsp` | 13,200 BC -- 17,191 AD | ~3.1 GB | Archeoastronomy, deep historical research |

### Code snippet

```python
from libephemeris import set_precision_tier

set_precision_tier("extended")  # 30,000 years of coverage
```

### CLI snippet

```bash
libephemeris download:medium      # recommended for most users
libephemeris status               # verify installed data
```

---

## [Thread Safety]

### Headline

Built for concurrency.

### Description

The global API mirrors pyswisseph's shared-state model for compatibility. For
concurrent applications, `EphemerisContext` provides fully isolated, thread-safe
instances with independent configuration.

### Code snippet

```python
from libephemeris import EphemerisContext

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)
pos, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
```

---

## [By the Numbers]

> Optional section. Compact stats strip or small grid. No commentary needed --
> the numbers speak for themselves.

| Metric | Value |
| --- | --- |
| Lines of library code | 87,000+ |
| Lines of tests | 136,000+ |
| Releases | 20 (since Jan 2024) |
| Supported Python versions | 3.9 -- 3.13 |
| License | LGPL-3.0 |

---

## [Getting Started]

### Headline

Get up and running in three commands.

### Steps

```bash
# 1. Install
pip install libephemeris

# 2. Download ephemeris data (~128 MB)
libephemeris download:medium

# 3. Verify
libephemeris status
```

### Links (inline or as buttons)

- Quick Start Guide
- Migration from pyswisseph _(docs/migration-guide.md)_
- API Reference _(docs/api_reference.rst)_
- Examples _(examples/)_

---

## [Footer]

### Copyright

2024--2026 Giacomo Battaglia

### License

LGPL-3.0. Not affiliated with Astrodienst AG.

### Links

GitHub | PyPI | Issues | Changelog

### Stability note (small text)

The public API may change between minor versions before 1.0.
