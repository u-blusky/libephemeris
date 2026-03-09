# LEB Implementation Plan: Binary Ephemeris Mode

> **STATUS: HISTORICAL DOCUMENT — ORIGINAL IMPLEMENTATION PLAN**
> This document is the original design plan for the `.leb` binary ephemeris mode.
> It has been **superseded** by the [Technical Guide](guide.md) and
> [Algorithms & Theory](algorithms.md) as authoritative references.
>
> **Key changes since this plan was written:**
> - `COORD_ICRS_BARY_SYSTEM` (type 4) added for outer planets (Jupiter-Pluto)
>   to store system barycenters with runtime COB correction
> - PPN gravitational deflection (Sun/Jupiter/Saturn) added to Pipeline A
> - Bodies 13, 21, 22 tightened to interval=4, degree=15
> - `fast_calc_ut` uses `swe_deltat()` instead of `reader.delta_t()`
> - Generator grew from ~600 to ~3668 lines with vectorized evaluation
> - All 31 bodies achieve <0.001" precision on base and medium tiers
> - Parallelization removed; replaced by group generation workflow
>
> See the [Technical Guide](guide.md) for the current accurate state.

---

## Table of Contents

1. [Goal and Architecture](#1-goal-and-architecture)
2. [Two-Mode Design](#2-two-mode-design)
3. [File Format Specification](#3-file-format-specification)
4. [Phase 0: Format Module (`leb_format.py`)](#4-phase-0-format-module)
5. [Phase 1: Generator (`scripts/generate_leb.py`)](#5-phase-1-generator)
6. [Phase 2: Reader (`leb_reader.py`)](#6-phase-2-reader)
7. [Phase 3: Fast Calculation Pipeline (`fast_calc.py`)](#7-phase-3-fast-calculation-pipeline)
8. [Phase 4: Integration](#8-phase-4-integration)
9. [Phase 5: Test Suite](#9-phase-5-test-suite)
10. [Future: Rust Drop-in Replacement](#10-future-rust-drop-in-replacement)
11. [File Inventory](#11-file-inventory)
12. [Risk Register](#12-risk-register)
13. [Reference: Current Pipeline Details](#13-reference-current-pipeline-details)

---

## 1. Goal and Architecture

### 1.1 Objective

Add a **binary ephemeris mode** to LibEphemeris that delivers 10-20x speedup over
the current Skyfield-based pipeline while maintaining:

- **Identical public API** -- `swe_calc_ut()`, `swe_houses()`, etc. unchanged
- **Identical numerical results** -- sub-arcsecond agreement with Skyfield mode
- **Zero new runtime dependencies** -- uses only stdlib (`mmap`, `struct`, `math`)
- **Transparent fallback** -- bodies not in the `.leb` file automatically use Skyfield
- **Same `.leb` file usable by future Rust port** -- format is language-agnostic

### 1.2 Why This Works

The current bottleneck is Skyfield's `.observe()` call chain. Each `swe_calc_ut()`
with `SEFLG_SPEED` triggers **9-12 Chebyshev polynomial evaluations** inside JPL
kernel segments (3 full Skyfield pipelines: position + 2 for central-difference
velocity). The `.leb` approach:

1. **Pre-fits** Chebyshev polynomials to Skyfield/JPL output (done once, offline)
2. **Evaluates** those polynomials directly with Clenshaw algorithm (~1.5us per eval)
3. **Gets velocity for free** via analytical Chebyshev derivative (eliminates 3x overhead)

This is the same technique JPL uses internally -- the process creates a
purpose-built cache of the results.

### 1.3 Implementation Strategy

**Python-first, Rust later.** The entire implementation is pure Python. Once
validated and stable, the reader + pipeline can be ported to Rust as a drop-in
replacement that reads the same `.leb` file format.

---

## 2. Two-Mode Design

### 2.1 User-Facing API

```python
import libephemeris as ephem

# Mode 1: Current (Skyfield) -- DEFAULT, no changes needed
result, flag = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)

# Mode 2: Binary (.leb) -- activated by a single call
ephem.set_leb_file("/path/to/ephemeris.leb")
result, flag = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)  # same API, ~14x faster

# Disable binary mode (back to Skyfield)
ephem.set_leb_file(None)

# Environment variable alternative
# export LIBEPHEMERIS_LEB=/path/to/ephemeris.leb
```

### 2.2 Internal Dispatch (in `planets.py::swe_calc_ut`)

```
swe_calc_ut(tjd_ut, ipl, iflag)
    |
    +-- reader = state.get_leb_reader()
    |
    +-- if reader is not None:
    |       try:
    |           return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)
    |       except (KeyError, ValueError, LEBRangeError):
    |           pass   # body not in .leb or JD out of range
    |
    +-- # ... existing Skyfield code, COMPLETELY UNCHANGED ...
```

**Fallback is automatic and silent.** If a body is not precomputed in the `.leb`
file (e.g., a custom asteroid), the existing Skyfield pipeline handles it. The
user never needs to know which mode served which body.

### 2.3 Scope of Binary Mode

What the `.leb` file covers:

| Body Category | Bodies | Covered? |
|---|---|---|
| Standard planets | Sun, Moon, Mercury-Pluto, Earth | Yes |
| Main asteroids | Chiron, Ceres, Pallas, Juno, Vesta | Yes |
| Lunar nodes | Mean Node, True Node | Yes |
| Lilith variants | Mean, True/Osculating, Interpolated Apogee/Perigee | Yes |
| Uranian hypotheticals | Cupido-Poseidon, Transpluto | Yes |
| Nutation | dpsi, deps (IAU 2006/2000A) | Yes |
| Delta-T | TT - UT1 | Yes |
| Fixed stars | 102-star catalog | Yes (stored as records, not Chebyshev) |
| Custom asteroids | User-registered SPK bodies | No (fallback to Skyfield) |
| Planetary moons | Galilean, Titan, etc. | No (fallback to Skyfield) |

What **stays in Python regardless** (not affected by `.leb`):

- House system calculations (`houses.py`) -- pure trig math
- Ayanamsha formula-based modes (38 of 43) -- pure polynomial math
- Eclipse search (`eclipse.py`) -- algorithmic, but benefits from faster `swe_calc_ut()`
- Crossing/transit calculations (`crossing.py`) -- same benefit
- Atmospheric extinction/heliacal (`extinction.py`, `heliacal.py`)

---

## 3. File Format Specification

### 3.1 Design Principles

- **O(1) segment lookup** -- integer division, no binary search
- **Memory-mapped** -- `mmap.mmap(READ_ONLY)`, zero-copy reads via `struct.unpack_from`
- **Self-describing** -- header contains all metadata needed to read the file
- **Language-agnostic** -- can be read by Python, Rust, C, or any language with mmap
- **Little-endian** -- matches x86/ARM64 native byte order

### 3.2 Overall Layout

```
Byte offset  Content
-----------  ----------------------------------------
0x0000       File Header (64 bytes)
0x0040       Section Directory (N x 24 bytes)
variable     Section 0: Body Index
variable     Section 1: Chebyshev Coefficient Data (bulk, 95%+ of file)
variable     Section 2: Nutation Chebyshev Data
variable     Section 3: Delta-T Sparse Table
variable     Section 4: Star Catalog
variable     Section 5: Orbital Elements (hypothetical bodies)
```

### 3.3 File Header (64 bytes)

| Offset | Type     | Field              | Description |
|--------|----------|--------------------|-------------|
| 0      | `[u8;4]` | magic              | `b"LEB1"` -- format identifier |
| 4      | `u32`    | version            | `1` -- format version |
| 8      | `u32`    | section_count      | Number of sections in directory |
| 12     | `u32`    | body_count         | Number of bodies in body index |
| 16     | `f64`    | jd_start           | Earliest JD covered by this file |
| 24     | `f64`    | jd_end             | Latest JD covered by this file |
| 32     | `f64`    | generation_epoch   | JD when the file was generated |
| 40     | `u32`    | flags              | Reserved flags (0 for now) |
| 44     | `[u8;20]`| reserved           | Zero-filled, future use |

### 3.4 Section Directory (24 bytes per entry)

| Offset | Type  | Field       | Description |
|--------|-------|-------------|-------------|
| 0      | `u32` | section_id  | Section type (0=body_index, 1=chebyshev, 2=nutation, 3=delta_t, 4=stars, 5=orbital_elements) |
| 4      | `u32` | reserved    | Alignment padding |
| 8      | `u64` | offset      | Byte offset from file start |
| 16     | `u64` | size        | Section size in bytes |

### 3.5 Body Index Entry (48 bytes per body)

| Offset | Type  | Field          | Description |
|--------|-------|----------------|-------------|
| 0      | `i32` | body_id        | SE_SUN (0), SE_MOON (1), ..., SE_CUPIDO (40), etc. |
| 4      | `u32` | coord_type     | 0=ICRS_BARY, 1=ECLIPTIC, 2=HELIO_ECL |
| 8      | `u32` | segment_count  | Number of Chebyshev segments |
| 12     | `f64` | jd_start       | Start JD for this body |
| 20     | `f64` | jd_end         | End JD for this body |
| 28     | `f64` | interval_days  | Duration of each segment in days |
| 36     | `u32` | degree         | Polynomial degree (9-16) |
| 40     | `u32` | components     | Number of components (3 for xyz or lon/lat/dist) |
| 44     | `u64` | data_offset    | Byte offset to first coefficient in Section 1 |

### 3.6 Chebyshev Coefficient Storage

For each segment, coefficients are stored contiguously as `f64` (8 bytes each):

```
[c0_comp0, c1_comp0, ..., cN_comp0, c0_comp1, c1_comp1, ..., cN_comp1, ...]
```

Where `N = degree`. Total bytes per segment: `components * (degree + 1) * 8`.

### 3.7 O(1) Segment Lookup

Given Julian Day `jd` and body entry `b`:

```python
segment_idx = int((jd - b.jd_start) / b.interval_days)
segment_idx = max(0, min(segment_idx, b.segment_count - 1))  # clamp
byte_offset = b.data_offset + segment_idx * b.components * (b.degree + 1) * 8
```

Then read `components * (degree + 1)` doubles from `byte_offset`.

### 3.8 Coordinate Types

Bodies are stored in one of 3 coordinate frames, chosen to minimize
runtime transforms:

| coord_type | Frame | Bodies | Stored Components | Why |
|---|---|---|---|---|
| 0 = ICRS_BARY | ICRS barycentric | Sun-Pluto, Earth, Chiron, Ceres-Vesta | (x, y, z) in AU | One dataset supports ALL output frames (geocentric, helio, equatorial, sidereal, J2000) via runtime transforms |
| 1 = ECLIPTIC | Ecliptic of date | Mean/True Node, Mean/True/Osculating Lilith, Interp. Apogee/Perigee | (lon, lat, dist) in deg/deg/AU | Already final coordinates; minimal runtime work |
| 2 = HELIO_ECL | Heliocentric ecliptic | Uranians (40-47), Transpluto (48) | (lon, lat, dist) in deg/deg/AU | Already in final frame; `_calc_body()` returns heliocentric today |

### 3.9 Nutation Section

Same Chebyshev format as bodies, but with 2 components (dpsi, deps) in radians.
Parameters: interval=16 days, degree=16.

### 3.10 Delta-T Section

Sparse table of `(jd, delta_t_days)` pairs, sampled every 30 days.
Each entry is 2 x f64 = 16 bytes. Reader uses bisection + cubic interpolation.

### 3.11 Star Catalog Section

102 fixed records, each containing:

| Offset | Type  | Field     | Description |
|--------|-------|-----------|-------------|
| 0      | `i32` | star_id   | Internal star ID |
| 4      | `f64` | ra_j2000  | Right ascension at J2000 (degrees) |
| 12     | `f64` | dec_j2000 | Declination at J2000 (degrees) |
| 20     | `f64` | pm_ra     | Proper motion in RA (deg/yr, includes cos(dec)) |
| 28     | `f64` | pm_dec    | Proper motion in Dec (deg/yr) |
| 36     | `f64` | parallax  | Parallax (arcsec) |
| 44     | `f64` | rv        | Radial velocity (km/s) |
| 52     | `f64` | magnitude | Visual magnitude |
| 60     | `[u8;4]` | reserved | Padding to 64 bytes |

Stars don't need Chebyshev precomputation. Proper motion is linear:
`ra(t) = ra_j2000 + pm_ra * (t - J2000)`, computed at runtime with 4 arithmetic ops.

### 3.12 Chebyshev Parameters Per Body

| Body       | Interval (days) | Degree | Bytes/segment | Segments (200yr) | Total    |
|------------|-----------------|--------|---------------|-------------------|----------|
| Moon       | 4               | 13     | 336           | 18,250            | 6.1 MB   |
| Mercury    | 16              | 15     | 384           | 4,563             | 1.8 MB   |
| Venus      | 16              | 13     | 336           | 4,563             | 1.5 MB   |
| Sun (EMB)  | 32              | 13     | 336           | 2,282             | 767 KB   |
| Earth      | 4               | 13     | 336           | 18,250            | 6.1 MB   |
| Mars       | 16              | 13     | 336           | 4,563             | 1.5 MB   |
| Jupiter    | 32              | 13     | 336           | 2,282             | 767 KB   |
| Saturn     | 32              | 13     | 336           | 2,282             | 767 KB   |
| Uranus     | 64              | 13     | 336           | 1,141             | 384 KB   |
| Neptune    | 64              | 13     | 336           | 1,141             | 384 KB   |
| Pluto      | 32              | 13     | 336           | 2,282             | 767 KB   |
| Chiron     | 8               | 13     | 336           | 9,125             | 3.1 MB   |
| Ceres-Vesta (4) | 8          | 13     | 336           | 4 x 9,125        | 12.3 MB  |
| Hypotheticals (9) | 32       | 13     | 336           | 9 x 2,282        | 6.9 MB   |
| Nutation   | 16              | 16     | 272           | 4,563             | 1.2 MB   |
| **Total (200yr)** |          |        |               |                   | **~44 MB** |

For the full DE440 range (1550-2650, 1100 years): ~175 MB.
For the full DE441 range (-13200 to +17191, 30,000 years): ~4.8 GB.

### 3.13 Longitude Wrap-Around

Ecliptic longitude has a discontinuity at 0/360 degrees. For bodies stored in
ecliptic coordinates (coord_type 1 and 2):

1. **Generator:** Before Chebyshev fitting, **unwrap** the longitude series
   (remove 360-degree jumps using `numpy.unwrap` with `period=360`).
2. **Generator:** Fit Chebyshev to the unwrapped (continuous) series.
3. **Reader:** After Clenshaw evaluation, **re-wrap** with `degnorm()` (`% 360.0`).

The generator MUST verify that the unwrapped fit reproduces correct results
across wrap boundaries (test with random dates near 0/360 crossings).

---

## 4. Phase 0: Format Module

### 4.1 File: `libephemeris/leb_format.py` (~150 lines)

This module defines format constants, dataclasses, and serialization helpers.
It is used by both the generator (write) and reader (read).

### 4.2 Contents

```python
from __future__ import annotations

import struct
from dataclasses import dataclass
from typing import List, Optional

# --- Format constants ---

MAGIC = b"LEB1"
VERSION = 1

# Coordinate types
COORD_ICRS_BARY = 0    # ICRS barycentric (x, y, z) in AU
COORD_ECLIPTIC = 1     # Ecliptic of date (lon, lat, dist) in deg/deg/AU
COORD_HELIO_ECL = 2    # Heliocentric ecliptic (lon, lat, dist) in deg/deg/AU

# Section IDs
SECTION_BODY_INDEX = 0
SECTION_CHEBYSHEV = 1
SECTION_NUTATION = 2
SECTION_DELTA_T = 3
SECTION_STARS = 4
SECTION_ORBITAL_ELEMENTS = 5

# Struct formats (little-endian)
HEADER_FMT = "<4sIIIdddI20s"      # 64 bytes
HEADER_SIZE = 64
SECTION_DIR_FMT = "<IIQQ"         # 24 bytes
SECTION_DIR_SIZE = 24
BODY_ENTRY_FMT = "<iIIddddIIQ"   # 48 bytes -- NOTE: needs exact calc
BODY_ENTRY_SIZE = 48
STAR_ENTRY_FMT = "<iddddddd4s"   # 64 bytes
STAR_ENTRY_SIZE = 64

# --- Dataclasses ---

@dataclass
class FileHeader:
    magic: bytes
    version: int
    section_count: int
    body_count: int
    jd_start: float
    jd_end: float
    generation_epoch: float
    flags: int

@dataclass
class SectionEntry:
    section_id: int
    offset: int
    size: int

@dataclass
class BodyEntry:
    body_id: int
    coord_type: int        # COORD_ICRS_BARY, COORD_ECLIPTIC, COORD_HELIO_ECL
    segment_count: int
    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int        # 3 for xyz/lonlatdist
    data_offset: int       # byte offset into chebyshev section

@dataclass
class StarEntry:
    star_id: int
    ra_j2000: float        # degrees
    dec_j2000: float       # degrees
    pm_ra: float           # deg/yr (includes cos(dec) factor)
    pm_dec: float          # deg/yr
    parallax: float        # arcsec
    rv: float              # km/s
    magnitude: float

# --- Body parameter table ---
# Maps body_id -> (interval_days, degree, coord_type, components)

BODY_PARAMS: dict[int, tuple[float, int, int, int]] = {
    # SE_SUN through SE_PLUTO, SE_EARTH: ICRS barycentric
    0:  (32, 13, COORD_ICRS_BARY, 3),   # SE_SUN
    1:  (4,  13, COORD_ICRS_BARY, 3),   # SE_MOON
    2:  (16, 15, COORD_ICRS_BARY, 3),   # SE_MERCURY
    3:  (16, 13, COORD_ICRS_BARY, 3),   # SE_VENUS
    4:  (16, 13, COORD_ICRS_BARY, 3),   # SE_MARS
    5:  (32, 13, COORD_ICRS_BARY, 3),   # SE_JUPITER
    6:  (32, 13, COORD_ICRS_BARY, 3),   # SE_SATURN
    7:  (64, 13, COORD_ICRS_BARY, 3),   # SE_URANUS
    8:  (64, 13, COORD_ICRS_BARY, 3),   # SE_NEPTUNE
    9:  (32, 13, COORD_ICRS_BARY, 3),   # SE_PLUTO
    14: (4,  13, COORD_ICRS_BARY, 3),   # SE_EARTH

    # Lunar nodes/Lilith: ecliptic direct
    10: (8,  13, COORD_ECLIPTIC, 3),    # SE_MEAN_NODE  (lon, 0, 0)
    11: (8,  13, COORD_ECLIPTIC, 3),    # SE_TRUE_NODE  (lon, lat, dist)
    12: (8,  13, COORD_ECLIPTIC, 3),    # SE_MEAN_APOG  (lon, lat, 0)
    13: (8,  13, COORD_ECLIPTIC, 3),    # SE_OSCU_APOG  (lon, lat, dist)
    21: (8,  13, COORD_ECLIPTIC, 3),    # SE_INTP_APOG  (lon, lat, dist)
    22: (8,  13, COORD_ECLIPTIC, 3),    # SE_INTP_PERG  (lon, lat, dist)

    # Main asteroids: ICRS barycentric
    15: (8,  13, COORD_ICRS_BARY, 3),   # SE_CHIRON
    17: (8,  13, COORD_ICRS_BARY, 3),   # SE_CERES
    18: (8,  13, COORD_ICRS_BARY, 3),   # SE_PALLAS
    19: (8,  13, COORD_ICRS_BARY, 3),   # SE_JUNO
    20: (8,  13, COORD_ICRS_BARY, 3),   # SE_VESTA

    # Uranian hypotheticals: heliocentric ecliptic
    40: (32, 13, COORD_HELIO_ECL, 3),   # SE_CUPIDO
    41: (32, 13, COORD_HELIO_ECL, 3),   # SE_HADES
    42: (32, 13, COORD_HELIO_ECL, 3),   # SE_ZEUS
    43: (32, 13, COORD_HELIO_ECL, 3),   # SE_KRONOS
    44: (32, 13, COORD_HELIO_ECL, 3),   # SE_APOLLON
    45: (32, 13, COORD_HELIO_ECL, 3),   # SE_ADMETOS
    46: (32, 13, COORD_HELIO_ECL, 3),   # SE_VULKANUS
    47: (32, 13, COORD_HELIO_ECL, 3),   # SE_POSEIDON
    48: (32, 13, COORD_HELIO_ECL, 3),   # SE_ISIS (Transpluto)
}

# --- Serialization helpers ---

def write_header(buf: bytearray, header: FileHeader) -> None: ...
def read_header(data: bytes) -> FileHeader: ...
def write_section_dir(buf: bytearray, offset: int, entry: SectionEntry) -> None: ...
def read_section_dir(data: bytes, offset: int) -> SectionEntry: ...
def write_body_entry(buf: bytearray, offset: int, entry: BodyEntry) -> None: ...
def read_body_entry(data: bytes, offset: int) -> BodyEntry: ...
def write_star_entry(buf: bytearray, offset: int, entry: StarEntry) -> None: ...
def read_star_entry(data: bytes, offset: int) -> StarEntry: ...
```

### 4.3 Implementation Notes

- All struct operations use `struct.pack_into` / `struct.unpack_from` for
  zero-copy compatibility with mmap buffers.
- The `BODY_PARAMS` dict is the single source of truth for Chebyshev parameters.
  Both the generator and reader reference it.
- Struct format strings must be carefully validated against the offset tables
  in Section 3 to ensure correct alignment. Write a unit test that verifies
  `struct.calcsize(FMT) == EXPECTED_SIZE` for each format.

---

## 5. Phase 1: Generator

### 5.1 File: `scripts/generate_leb.py` (~600 lines)

CLI script that produces a `.leb` file using Skyfield and libephemeris's
existing analytical functions as the data source.

### 5.2 CLI Interface

```bash
# Generate 200-year file (default)
python scripts/generate_leb.py --output ephemeris.leb

# Custom range
python scripts/generate_leb.py --start 1550 --end 2650 --output de440_full.leb

# Generate with verification
python scripts/generate_leb.py --output ephemeris.leb --verify --verify-samples 500

# Parallel generation (per-body)
python scripts/generate_leb.py --output ephemeris.leb --workers 8
```

### 5.3 Generation Functions

#### `generate_body_icrs(body_id, jd_start, jd_end, params) -> coefficients`

For planets (SE_SUN through SE_PLUTO), Earth, Chiron, Ceres-Vesta:

1. For each segment interval:
   a. Compute Chebyshev nodes (degree+1 points) within the interval
   b. At each node, query Skyfield for **ICRS barycentric** position:
      ```python
      # For standard planets: use the JPL kernel segment directly
      planets = get_planets()
      target = planets[body_name]     # e.g., planets["mars"]
      t = ts.tt_jd(jd_node)
      pos = target.at(t).position.au  # [x, y, z] ICRS barycentric
      ```
      Note: For outer planets, the generator must get the **planet center** (not barycenter).
      Use `get_planet_target(planets, name)` which resolves to `_SpkCenterTarget`
      or `_CobCorrectedTarget` as appropriate.
   c. Fit Chebyshev coefficients using `numpy.polynomial.chebyshev.chebfit()`
   d. Verify: evaluate the fit at 10 intermediate points, check error < 0.001"

**Critical: Getting true ICRS barycentric positions.**

The Skyfield `target.at(t).position.au` gives the ICRS barycentric position
directly from the SPK kernel, with no observer, no light-time, no aberration.
This is the desired result -- the raw position vector.

For the Earth, the calculation requires `planets["earth"].at(t).position.au`. This gives the
geocenter position in ICRS barycentric coordinates.

For outer planets (Jupiter-Pluto), the generator needs the planet **center**, not the
system barycenter. The current code in `planets.py:567-618` handles this via
`get_planet_target()` which uses `planet_centers.bsp` or analytical COB
corrections. The generator must use the same logic to get accurate planet
center positions.

#### `generate_body_ecliptic(body_id, jd_start, jd_end, params) -> coefficients`

For lunar nodes, Lilith variants:

1. For each segment interval:
   a. Compute Chebyshev nodes within the interval
   b. At each node, call the existing libephemeris function:
      ```python
      # Examples:
      lon = calc_mean_lunar_node(jd_tt)                  # SE_MEAN_NODE
      lon, lat, dist = calc_true_lunar_node(jd_tt)       # SE_TRUE_NODE
      lon, lat = calc_mean_lilith_with_latitude(jd_tt)   # SE_MEAN_APOG
      lon, lat, dist = calc_true_lilith(jd_tt)           # SE_OSCU_APOG
      lon, lat, dist = calc_interpolated_apogee(jd_tt)   # SE_INTP_APOG
      lon, lat, dist = calc_interpolated_perigee(jd_tt)  # SE_INTP_PERG
      ```
   c. **Unwrap longitude** before fitting (remove 360-degree jumps)
   d. Fit Chebyshev coefficients
   e. Verify after re-wrapping with `degnorm()`

#### `generate_body_helio(body_id, jd_start, jd_end, params) -> coefficients`

For Uranian hypotheticals (SE_CUPIDO through SE_POSEIDON) and Transpluto:

1. Same Chebyshev fitting procedure as ecliptic bodies
2. Calls `hypothetical.calc_uranian_planet(jd_tt, body_id)` or
   `hypothetical.calc_transpluto(jd_tt)` at each node
3. These are already in heliocentric ecliptic coordinates

#### `generate_nutation(jd_start, jd_end) -> coefficients`

1. Interval = 32 days, degree = 16
2. At each Chebyshev node: `dpsi, deps = erfa.nut06a(J2000, jd_tt - J2000)`
3. Stored in radians (native erfa output)

#### `generate_delta_t(jd_start, jd_end) -> table`

1. Sample `swe_deltat(jd)` every 30 days
2. Store as sparse `(jd, delta_t_days)` table
3. Reader interpolates with cubic spline

#### `generate_star_catalog() -> records`

1. Extract 102 star records from `fixed_stars.STAR_CATALOG`
2. Store J2000 positions, proper motions, parallax, RV, magnitude
3. No Chebyshev needed -- proper motion is linear

#### `verify_leb(leb_path, n_samples=500)`

Post-generation validation:

1. For each body in the `.leb` file:
   a. Generate `n_samples` random JDs within the file's range
   b. Evaluate via `LEBReader.eval_body()`
   c. Compare with Skyfield/analytical reference
   d. Report max error in arcseconds
   e. **FAIL if any body exceeds 0.001 arcsecond error**

### 5.4 Chebyshev Fitting Details

The fitting uses `numpy.polynomial.chebyshev.chebfit()`:

```python
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

def fit_segment(func, jd_start, jd_end, degree):
    """Fit a function over [jd_start, jd_end] with Chebyshev polynomials."""
    # Chebyshev nodes (Type I) on [-1, 1]
    nodes = np.cos(np.pi * (np.arange(degree + 1) + 0.5) / (degree + 1))
    # Map to [jd_start, jd_end]
    jd_nodes = 0.5 * (jd_end - jd_start) * nodes + 0.5 * (jd_start + jd_end)
    # Evaluate function at nodes
    values = np.array([func(jd) for jd in jd_nodes])
    # Fit
    coeffs = chebfit(nodes, values, degree)
    return coeffs
```

The Clenshaw evaluation in the reader must use the same domain mapping:
```python
tau = 2.0 * (jd - jd_mid) / (jd_end - jd_start)  # map to [-1, 1]
```

### 5.5 Estimated Generation Times

| Range | Duration | File Size |
|---|---|---|
| 200 years (1900-2100) | ~10 minutes | ~16 MB |
| 1100 years (1550-2650) | ~55 minutes | ~88 MB |
| 30,000 years (DE441 full) | ~25 hours | ~2.4 GB |

Generation is embarrassingly parallel per body. With `--workers 8`, the 200-year
generation should complete in ~2 minutes.

---

## 6. Phase 2: Reader

### 6.1 File: `libephemeris/leb_reader.py` (~400 lines)

### 6.2 Class Interface

```python
class LEBReader:
    """Memory-mapped reader for .leb binary ephemeris files."""

    def __init__(self, path: str) -> None:
        """Open and parse a .leb file.

        Args:
            path: Path to the .leb file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file is not a valid .leb file.
        """

    def eval_body(self, body_id: int, jd: float) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Evaluate a body's position and velocity at a given Julian Day.

        Args:
            body_id: SE_* body constant (e.g., SE_SUN=0, SE_MOON=1).
            jd: Julian Day in TT (Terrestrial Time).

        Returns:
            ((pos_0, pos_1, pos_2), (vel_0, vel_1, vel_2))

            For ICRS_BARY bodies: ((x, y, z), (vx, vy, vz)) in AU and AU/day.
            For ECLIPTIC bodies: ((lon, lat, dist), (dlon, dlat, ddist)) in deg, deg, AU.
            For HELIO_ECL bodies: same as ECLIPTIC.

            Longitude values are already wrapped to [0, 360) for ecliptic bodies.
            Velocity is the analytical Chebyshev derivative (not central difference).

        Raises:
            KeyError: If body_id is not in this .leb file.
            ValueError: If jd is outside the body's coverage range.
        """

    def eval_nutation(self, jd_tt: float) -> tuple[float, float]:
        """Evaluate nutation angles at a given Julian Day.

        Args:
            jd_tt: Julian Day in TT.

        Returns:
            (dpsi, deps) in radians (IAU 2006/2000A).
        """

    def delta_t(self, jd: float) -> float:
        """Get Delta-T (TT - UT1) at a given Julian Day.

        Args:
            jd: Julian Day (UT or TT -- the difference is negligible for lookup).

        Returns:
            Delta-T in days.
        """

    def get_star(self, star_id: int) -> StarEntry:
        """Look up a fixed star record.

        Args:
            star_id: Internal star ID.

        Returns:
            StarEntry with J2000 position, proper motion, etc.

        Raises:
            KeyError: If star_id is not in the catalog.
        """

    def has_body(self, body_id: int) -> bool:
        """Check if a body is available in this .leb file."""

    @property
    def jd_range(self) -> tuple[float, float]:
        """Return (jd_start, jd_end) covered by this file."""

    def close(self) -> None:
        """Close the memory-mapped file."""
```

### 6.3 Clenshaw Algorithm (Pure Python)

The Clenshaw algorithm evaluates a Chebyshev series without forming the
polynomials explicitly. For position and velocity simultaneously:

```python
def _clenshaw_with_derivative(coeffs: tuple[float, ...], tau: float) -> tuple[float, float]:
    """Evaluate Chebyshev series and its derivative at tau in [-1, 1].

    Args:
        coeffs: Chebyshev coefficients (c0, c1, ..., cN).
        tau: Evaluation point in [-1, 1].

    Returns:
        (value, derivative) where derivative is d(value)/d(tau).
    """
    n = len(coeffs) - 1
    if n == 0:
        return coeffs[0], 0.0

    # Forward Clenshaw for value
    b_k1 = 0.0  # b_{k+1}
    b_k2 = 0.0  # b_{k+2}
    for k in range(n, 0, -1):
        b_k = coeffs[k] + 2.0 * tau * b_k1 - b_k2
        b_k2 = b_k1
        b_k1 = b_k
    value = coeffs[0] + tau * b_k1 - b_k2

    # Forward Clenshaw for derivative (using Chebyshev derivative recurrence)
    # T'_n(x) = n * U_{n-1}(x), so we evaluate the derivative series
    d_k1 = 0.0
    d_k2 = 0.0
    for k in range(n, 1, -1):
        d_k = k * coeffs[k] + 2.0 * tau * d_k1 - d_k2
        d_k2 = d_k1
        d_k1 = d_k
    derivative = coeffs[1] + 2.0 * tau * d_k1 - d_k2  # ... needs care with recurrence

    return value, derivative
```

**Why pure Python, not numpy?** For single-point evaluation (the common case),
numpy array creation overhead (~5us) dominates the actual Clenshaw loop (~1.5us
for degree 13). Pure Python `float` operations are faster for scalar work.

**Domain scaling for derivative:** The Clenshaw derivative gives `d/d(tau)`.
To convert to `d/d(jd)`, scale by `2.0 / interval_days`:

```python
velocity = raw_derivative * 2.0 / body.interval_days
```

### 6.4 mmap Usage

```python
import mmap

class LEBReader:
    def __init__(self, path: str):
        self._file = open(path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._header = read_header(self._mm[:HEADER_SIZE])
        # ... parse sections and body index ...
        # Store body entries in a dict for O(1) lookup:
        self._bodies: dict[int, BodyEntry] = {}

    def eval_body(self, body_id: int, jd: float):
        body = self._bodies[body_id]  # KeyError if not present
        # O(1) segment lookup
        seg_idx = int((jd - body.jd_start) / body.interval_days)
        seg_idx = max(0, min(seg_idx, body.segment_count - 1))
        # Read coefficients from mmap (zero-copy)
        n_coeffs = body.components * (body.degree + 1)
        offset = body.data_offset + seg_idx * n_coeffs * 8
        coeffs = struct.unpack_from(f"<{n_coeffs}d", self._mm, offset)
        # Compute tau (map jd to [-1, 1] within segment)
        seg_start = body.jd_start + seg_idx * body.interval_days
        seg_end = seg_start + body.interval_days
        seg_mid = 0.5 * (seg_start + seg_end)
        tau = 2.0 * (jd - seg_mid) / body.interval_days
        # Evaluate each component
        pos = []
        vel = []
        deg1 = body.degree + 1
        for c in range(body.components):
            comp_coeffs = coeffs[c * deg1 : (c + 1) * deg1]
            val, deriv = _clenshaw_with_derivative(comp_coeffs, tau)
            pos.append(val)
            vel.append(deriv * 2.0 / body.interval_days)
        # Wrap longitude for ecliptic bodies
        if body.coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL):
            pos[0] = pos[0] % 360.0
        return tuple(pos), tuple(vel)
```

---

## 7. Phase 3: Fast Calculation Pipeline

### 7.1 File: `libephemeris/fast_calc.py` (~800 lines)

This module reimplements the Skyfield pipeline using `LEBReader` as the data
source. It handles coordinate transforms, light-time correction, aberration,
and flag dispatch.

### 7.2 Entry Point

```python
def fast_calc_ut(
    reader: LEBReader,
    tjd_ut: float,
    ipl: int,
    iflag: int,
) -> tuple[tuple[float, float, float, float, float, float], int]:
    """Fast equivalent of swe_calc_ut() using precomputed .leb data.

    Args:
        reader: An open LEBReader instance.
        tjd_ut: Julian Day in Universal Time (UT1).
        ipl: Planet/body ID (SE_SUN, SE_MOON, etc.).
        iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.).

    Returns:
        Same as swe_calc_ut(): ((lon, lat, dist, dlon, dlat, ddist), iflag)

    Raises:
        KeyError: If body is not in the .leb file (caller should fall back).
        ValueError: If JD is outside the .leb file's range.
    """
```

### 7.3 Pipeline A: ICRS Barycentric Bodies

This handles Sun, Moon, Mercury-Pluto, Earth, Chiron, Ceres-Vesta.

**Step-by-step:**

```
1. DELTA-T CONVERSION
   delta_t = reader.delta_t(tjd_ut)   # TT - UT1 in days
   jd_tt = tjd_ut + delta_t

2. BODY POSITIONS (from precomputed Chebyshev)
   (earth_xyz, earth_vel) = reader.eval_body(SE_EARTH, jd_tt)    # ICRS bary
   (target_xyz, target_vel) = reader.eval_body(ipl, jd_tt)        # ICRS bary

3. OBSERVER SELECTION
   if SEFLG_HELCTR:
       (sun_xyz, sun_vel) = reader.eval_body(SE_SUN, jd_tt)
       observer_xyz, observer_vel = sun_xyz, sun_vel
   elif SEFLG_BARYCTR:
       observer_xyz = (0, 0, 0)
       observer_vel = (0, 0, 0)
   else:  # geocentric (default)
       observer_xyz, observer_vel = earth_xyz, earth_vel

4. GEOMETRIC VECTOR
   geo = target_xyz - observer_xyz    # 3-component subtraction

5. LIGHT-TIME CORRECTION (unless SEFLG_TRUEPOS)
   if not SEFLG_TRUEPOS:
       for _ in range(3):  # 3 fixed-point iterations
           dist = sqrt(geo[0]**2 + geo[1]**2 + geo[2]**2)
           lt = dist / C_LIGHT_AU_DAY   # 173.1446326846693 AU/day
           (retarded_xyz, _) = reader.eval_body(ipl, jd_tt - lt)
           geo = retarded_xyz - observer_xyz

6. ABERRATION (unless SEFLG_NOABERR)
   if not (SEFLG_NOABERR or SEFLG_HELCTR or SEFLG_BARYCTR or SEFLG_TRUEPOS):
       # Annual aberration using Earth velocity
       geo = apply_aberration(geo, earth_vel)
       # Uses astrometry.apply_aberration_to_position()

7. COORDINATE TRANSFORM
   dist = sqrt(geo[0]**2 + geo[1]**2 + geo[2]**2)

   if SEFLG_EQUATORIAL and SEFLG_J2000:
       # ICRS J2000 equatorial -- geo is already in this frame!
       ra = atan2(geo[1], geo[0])  # radians -> degrees
       dec = asin(geo[2] / dist)
       lon_deg = degrees(ra) % 360.0
       lat_deg = degrees(dec)

   elif SEFLG_EQUATORIAL:
       # True equator of date
       # Apply precession + nutation rotation matrix to geo
       geo_eq = precession_nutation_matrix(jd_tt) @ geo  # 3x3 matrix
       ra = atan2(geo_eq[1], geo_eq[0])
       dec = asin(geo_eq[2] / dist)
       lon_deg = degrees(ra) % 360.0
       lat_deg = degrees(dec)

   elif SEFLG_J2000:
       # J2000 ecliptic
       eps = 23.4392911  # fixed J2000 mean obliquity
       xe = geo[0]
       ye = geo[1] * cos(eps_rad) + geo[2] * sin(eps_rad)
       ze = -geo[1] * sin(eps_rad) + geo[2] * cos(eps_rad)
       lon_deg = degrees(atan2(ye, xe)) % 360.0
       lat_deg = degrees(asin(ze / dist))

   else:
       # TRUE ECLIPTIC OF DATE (default) -- most common path
       # Step 1: Precess ICRS -> equatorial of date (precession + nutation)
       pn_matrix = precession_nutation_matrix(jd_tt)
       geo_eq = pn_matrix @ geo
       # Step 2: Rotate equatorial -> ecliptic using true obliquity
       dpsi, deps = reader.eval_nutation(jd_tt)
       eps_mean = mean_obliquity_iau2006(jd_tt)  # polynomial, ~0.1us
       eps_true = eps_mean + deps
       lon_deg, lat_deg = equatorial_to_ecliptic(geo_eq, eps_true, dist)

8. SIDEREAL CORRECTION (if SEFLG_SIDEREAL and not SEFLG_EQUATORIAL)
   if SEFLG_SIDEREAL and not SEFLG_EQUATORIAL:
       aya = _calc_ayanamsa_from_leb(reader, tjd_ut)  # uses precession polynomial
       lon_deg = (lon_deg - aya) % 360.0

9. VELOCITY
    # The analytical Chebyshev derivative from eval_body() is transformed
    # through the same rotation matrices as position (option (a) — rigorous).
    #
    # Pipeline steps for velocity:
    # 1. geo_vel = target_vel - observer_vel (geocentric velocity)
    # 2. Light-time: use velocity at retarded time (jd_tt - lt)
    # 3. Apply precession-nutation matrix to velocity vector
    # 4. Rotate equatorial → ecliptic using true obliquity
    # 5. Convert Cartesian velocity to spherical via
    #    _cartesian_velocity_to_spherical(x,y,z, vx,vy,vz)
    #
    # This is both faster (1 pipeline run instead of 3) and more precise
    # than the previous central-difference approach.

    if SEFLG_SPEED:
        # want_velocity=True causes _pipeline_icrs() to return 6 values
        lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
            reader, jd_tt, ipl, iflag, want_velocity=True
        )
    else:
        lon, lat, dist = _pipeline_icrs(reader, jd_tt, ipl, iflag)
        dlon, dlat, ddist = 0.0, 0.0, 0.0

10. RETURN
     return (lon_deg, lat_deg, dist_au, dlon, dlat, ddist), iflag
```

**Note on velocity (Step 9):** The analytical approach (option (a) from the
original design) was implemented. The Chebyshev derivative velocity vector
is transformed through the same linear rotation matrices already computed
for position. This eliminates the 2 extra pipeline evaluations that the
central-difference approach required, and avoids amplifying Chebyshev
fitting errors into velocity errors.

**Architectural limitation:** For ICRS bodies, the Cartesian-to-spherical
velocity conversion involves division by `r_xy` (distance in the ecliptic
plane). For nearby asteroids (Ceres ~2.7 AU, Pallas ~2.8 AU), the
`1/geocentric_distance` amplification factor produces latitude velocity
errors of 0.19-0.71 deg/day. This is inherent to the coordinate
transformation and cannot be fixed without changing the storage format.

### 7.4 Pipeline B: Ecliptic Direct Bodies

For SE_MEAN_NODE, SE_TRUE_NODE, SE_MEAN_APOG, SE_OSCU_APOG, SE_INTP_APOG,
SE_INTP_PERG:

```
1. (lon, lat, dist), (dlon, dlat, ddist) = reader.eval_body(ipl, jd_tt)
   # Already in ecliptic coordinates (lon wrapped to [0, 360))
   # Velocity from analytical Chebyshev derivative -- NO central difference!

2. if SEFLG_EQUATORIAL:
       lon, lat = cotrans(lon, lat, -true_obliquity)
       dlon, dlat = cotrans_sp velocity transform (analytical)

3. if SEFLG_J2000:
       # Precess from ecliptic-of-date to J2000
       lon, lat = _precess_ecliptic(lon, lat, jd_tt, J2000)

4. if SEFLG_SIDEREAL and not SEFLG_EQUATORIAL:
       aya = _calc_ayanamsa_from_leb(reader, tjd_ut)
       lon = (lon - aya) % 360.0

5. return (lon, lat, dist, dlon, dlat, ddist), iflag
```

**Key advantage:** Velocity is directly from the Chebyshev derivative.
No central difference needed. This makes Pipeline B extremely fast (~5us).

### 7.5 Pipeline C: Heliocentric Bodies

Identical to Pipeline B. Uranian hypotheticals and Transpluto are already
stored in heliocentric ecliptic coordinates, which is what `_calc_body()`
returns for these bodies today. No additional transforms needed.

### 7.6 Precession-Nutation Matrix

The precession-nutation matrix is needed for Pipeline A's coordinate transforms.
Implementation options (in order of preference):

1. **Use `erfa.pnm06a()`** if pyerfa is available (~0.02ms, cached via `cache.py`)
2. **Use `astrometry.py::_precession_nutation_matrix()`** as fallback
3. **Build from .leb nutation data:** Read `(dpsi, deps)` from the .leb file,
   compute mean obliquity via the IAU 2006 polynomial, build the rotation
   matrix manually. This avoids calling pyerfa entirely.

Option 3 is the most self-contained but adds ~100 lines of matrix construction.
For the initial implementation, use option 1/2 (existing code) and consider
option 3 for the Rust port.

### 7.7 Ayanamsha from .leb

For formula-based ayanamsha (38 of 43 modes), we can compute without Skyfield:

```python
def _calc_ayanamsa_from_leb(reader, tjd_ut):
    """Compute ayanamsa using only .leb data (no Skyfield)."""
    delta_t = reader.delta_t(tjd_ut)
    jd_tt = tjd_ut + delta_t
    T = (jd_tt - 2451545.0) / 36525.0

    # IAU 2006 general precession polynomial
    precession_arcsec = (
        5028.796195 * T
        + 1.1054348 * T**2
        + 0.00007964 * T**3
        - 0.000023857 * T**4
        - 0.0000000383 * T**5
    )

    # Get reference offset for the active sidereal mode
    aya_j2000 = AYANAMSHA_OFFSETS[sid_mode]

    mean_aya = aya_j2000 + precession_arcsec / 3600.0

    # True ayanamsa: add nutation in longitude
    dpsi, _ = reader.eval_nutation(jd_tt)
    nutation_deg = math.degrees(dpsi)
    return (mean_aya + nutation_deg) % 360.0
```

For the 5 star-based ayanamsha modes ("True Citra", "True Revati", etc.),
the fallback to Skyfield is needed because they require computing the apparent
position of a reference star. These modes raise `KeyError` in `fast_calc` and
the caller falls through to the Skyfield pipeline.

### 7.8 Flags Support Matrix

| Flag | Pipeline A (planets) | Pipeline B (ecliptic) | Pipeline C (helio) |
|---|---|---|---|
| `SEFLG_SPEED` | Central difference on final coords | Chebyshev derivative (free) | Chebyshev derivative (free) |
| `SEFLG_HELCTR` | Use Sun as observer | N/A (already ecliptic) | N/A (already helio) |
| `SEFLG_BARYCTR` | No observer subtraction | N/A | N/A |
| `SEFLG_EQUATORIAL` | ICRS->equatorial rotation | `cotrans()` | `cotrans()` |
| `SEFLG_J2000` | Fixed obliquity rotation | `_precess_ecliptic()` | `_precess_ecliptic()` |
| `SEFLG_SIDEREAL` | Subtract ayanamsa | Subtract ayanamsa | Subtract ayanamsa |
| `SEFLG_TRUEPOS` | Skip light-time | No effect | No effect |
| `SEFLG_NOABERR` | Skip aberration | No effect | No effect |
| `SEFLG_TOPOCTR` | NOT SUPPORTED (fall back to Skyfield) | NOT SUPPORTED | NOT SUPPORTED |
| `SEFLG_MOSEPH` | Stripped (ignored) | Stripped | Stripped |

**`SEFLG_TOPOCTR` fallback:** Topocentric calculations require the observer's
geographic position and Earth rotation parameters that are not in the `.leb`
file. These always fall through to Skyfield. This is acceptable because
topocentric mode is rarely used in batch calculations.

---

## 8. Phase 4: Integration

### 8.1 Changes to `state.py` (~30 lines added)

Add after existing globals (around line 155):

```python
_LEB_FILE: Optional[str] = None          # Path to .leb file
_LEB_READER: Optional["LEBReader"] = None # Cached LEBReader instance
```

New public functions:

```python
def set_leb_file(filepath: Optional[str]) -> None:
    """Set the .leb file path for binary ephemeris mode.

    Args:
        filepath: Path to a .leb file, or None to disable binary mode.
    """
    global _LEB_FILE, _LEB_READER
    if _LEB_READER is not None:
        _LEB_READER.close()
    _LEB_FILE = filepath
    _LEB_READER = None  # force re-creation on next access


def get_leb_reader() -> Optional["LEBReader"]:
    """Get the active LEBReader, if any.

    Returns:
        LEBReader instance if a .leb file is configured, None otherwise.
    """
    global _LEB_READER
    if _LEB_READER is None:
        path = _LEB_FILE or os.environ.get("LIBEPHEMERIS_LEB")
        if path is not None:
            from .leb_reader import LEBReader
            _LEB_READER = LEBReader(path)
    return _LEB_READER
```

Update `close()` (around line 979-998):

```python
# Add to global declarations at top of close():
global _LEB_FILE, _LEB_READER

# Add to reset block:
if _LEB_READER is not None:
    try:
        _LEB_READER.close()
    except Exception:
        pass
_LEB_FILE = None
_LEB_READER = None
```

### 8.2 Changes to `planets.py` (~10 lines modified)

In `swe_calc_ut()`, add before the existing pipeline (around line 770):

```python
def swe_calc_ut(tjd_ut, ipl, iflag):
    # ... existing SE_ECL_NUT handling (line 766-767) ...
    # ... existing SEFLG_MOSEPH stripping (line 770) ...

    # --- NEW: .leb fast path ---
    from .state import get_leb_reader
    reader = get_leb_reader()
    if reader is not None and ipl != SE_ECL_NUT:
        try:
            from . import fast_calc
            return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)
        except (KeyError, ValueError):
            pass  # body not in .leb or JD out of range, fall through
    # --- END NEW ---

    # ... existing validation, Skyfield pipeline, UNCHANGED ...
```

Similarly in `swe_calc()`:

```python
def swe_calc(tjd, ipl, iflag):
    # --- NEW: .leb fast path ---
    reader = get_leb_reader()
    if reader is not None and ipl != SE_ECL_NUT:
        try:
            from . import fast_calc
            return fast_calc.fast_calc_tt(reader, tjd, ipl, iflag)
        except (KeyError, ValueError):
            pass
    # --- END NEW ---

    # ... existing code UNCHANGED ...
```

### 8.3 Changes to `__init__.py` (~5 lines)

Add exports:

```python
from .state import set_leb_file, get_leb_reader
```

And add to `__all__`:

```python
"set_leb_file",
"get_leb_reader",
```

### 8.4 Changes to `context.py`

Add `.leb` support to `EphemerisContext`:

```python
class EphemerisContext:
    def __init__(self, ...):
        # ... existing init ...
        self._leb_file: Optional[str] = None
        self._leb_reader: Optional["LEBReader"] = None

    def set_leb_file(self, filepath: Optional[str]) -> None:
        if self._leb_reader is not None:
            self._leb_reader.close()
        self._leb_file = filepath
        self._leb_reader = None

    def get_leb_reader(self) -> Optional["LEBReader"]:
        if self._leb_reader is None and self._leb_file is not None:
            from .leb_reader import LEBReader
            self._leb_reader = LEBReader(self._leb_file)
        return self._leb_reader
```

---

## 9. Phase 5: Test Suite

### 9.1 Test Files

| File | LOC est. | What it tests |
|---|---|---|
| `tests/test_leb/test_leb_format.py` | ~200 | Format constants, struct sizes, round-trip serialization |
| `tests/test_leb/test_leb_reader.py` | ~250 | mmap opening, eval_body, eval_nutation, delta_t, edge cases |
| `tests/test_leb/test_fast_calc.py` | ~300 | Full pipeline comparison vs Skyfield for all bodies and flag combos |
| `tests/test_leb/test_generate_leb.py` | ~150 | Generator correctness (small date ranges) |

### 9.2 Key Test Cases

#### Format Tests (`test_leb_format.py`)

```python
class TestStructSizes:
    def test_header_size(self):
        assert struct.calcsize(HEADER_FMT) == HEADER_SIZE  # 64

    def test_body_entry_size(self):
        assert struct.calcsize(BODY_ENTRY_FMT) == BODY_ENTRY_SIZE  # 48

    def test_star_entry_size(self):
        assert struct.calcsize(STAR_ENTRY_FMT) == STAR_ENTRY_SIZE  # 64

class TestRoundTrip:
    def test_header_round_trip(self):
        header = FileHeader(MAGIC, VERSION, 6, 25, 2415020.5, 2488069.5, 2460000.5, 0)
        buf = bytearray(HEADER_SIZE)
        write_header(buf, header)
        recovered = read_header(bytes(buf))
        assert recovered == header

    def test_body_entry_round_trip(self):
        entry = BodyEntry(0, COORD_ICRS_BARY, 2282, 2415020.5, 2488069.5, 32.0, 13, 3, 1024)
        # ... round-trip test ...
```

#### Reader Tests (`test_leb_reader.py`)

```python
class TestClenshaw:
    def test_constant_polynomial(self):
        """Chebyshev [5.0] should evaluate to 5.0 everywhere."""
        val, deriv = _clenshaw_with_derivative((5.0,), 0.0)
        assert val == 5.0
        assert deriv == 0.0

    def test_linear_polynomial(self):
        """Chebyshev [a, b] = a + b*x, derivative = b."""
        val, deriv = _clenshaw_with_derivative((3.0, 2.0), 0.5)
        assert abs(val - 4.0) < 1e-14
        assert abs(deriv - 2.0) < 1e-14

    def test_known_chebyshev(self):
        """Verify against numpy.polynomial.chebyshev.chebval."""
        coeffs = (1.0, 0.5, -0.3, 0.1)
        for tau in [-1.0, -0.5, 0.0, 0.5, 1.0]:
            expected = np.polynomial.chebyshev.chebval(tau, coeffs)
            actual, _ = _clenshaw_with_derivative(coeffs, tau)
            assert abs(actual - expected) < 1e-13

class TestEvalBody:
    def test_sun_position_matches_skyfield(self, leb_reader, standard_jd):
        """Sun ICRS position from .leb matches Skyfield within 0.001 arcsec."""
        pos, vel = leb_reader.eval_body(SE_SUN, standard_jd)
        # Compare with Skyfield reference...
        assert angular_separation < 0.001 / 3600  # 0.001 arcsec in degrees
```

#### Pipeline Tests (`test_fast_calc.py`)

```python
class TestFastCalcVsSkyfield:
    """Compare fast_calc output with Skyfield pipeline for all bodies."""

    TOLERANCE_ARCSEC = 0.01  # 10 milliarcseconds

    @pytest.fixture
    def leb_reader(self):
        # Use a small test .leb file generated in conftest
        return LEBReader("tests/fixtures/test_ephemeris.leb")

    @pytest.mark.parametrize("ipl", [SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS,
        SE_MARS, SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO])
    def test_planet_ecliptic(self, leb_reader, ipl):
        """Planet ecliptic lon/lat matches Skyfield within tolerance."""
        for jd in TEST_DATES:
            fast_result, _ = fast_calc_ut(leb_reader, jd, ipl, SEFLG_SPEED)
            sky_result, _ = swe_calc_ut(jd, ipl, SEFLG_SPEED)  # Skyfield
            assert abs(fast_result[0] - sky_result[0]) < self.TOLERANCE_ARCSEC / 3600
            assert abs(fast_result[1] - sky_result[1]) < self.TOLERANCE_ARCSEC / 3600

    @pytest.mark.parametrize("ipl", [SE_MEAN_NODE, SE_TRUE_NODE,
        SE_MEAN_APOG, SE_OSCU_APOG, SE_INTP_APOG, SE_INTP_PERG])
    def test_lunar_ecliptic(self, leb_reader, ipl):
        """Lunar node/Lilith matches Skyfield within tolerance."""
        # ...

    @pytest.mark.parametrize("flag", [
        0, SEFLG_SPEED, SEFLG_EQUATORIAL, SEFLG_J2000,
        SEFLG_SIDEREAL, SEFLG_HELCTR, SEFLG_NOABERR,
        SEFLG_SPEED | SEFLG_EQUATORIAL,
        SEFLG_SPEED | SEFLG_SIDEREAL,
    ])
    def test_all_flags(self, leb_reader, flag):
        """Test all supported flag combinations."""
        # ...

    def test_fallback_for_unknown_body(self, leb_reader):
        """Bodies not in .leb should raise KeyError (caller handles fallback)."""
        with pytest.raises(KeyError):
            fast_calc_ut(leb_reader, 2451545.0, 99999, 0)

    def test_fallback_for_topocentric(self, leb_reader):
        """SEFLG_TOPOCTR should raise or be handled (falls back to Skyfield)."""
        # ...
```

### 9.3 Test Fixtures

Add to `tests/conftest.py` or `tests/test_leb/conftest.py`:

```python
@pytest.fixture(scope="session")
def test_leb_file(tmp_path_factory):
    """Generate a small .leb file for testing (10-year range)."""
    path = tmp_path_factory.mktemp("leb") / "test.leb"
    # Generate a minimal .leb covering 2020-2030
    generate_leb(
        output=str(path),
        start_year=2020,
        end_year=2030,
        bodies=[SE_SUN, SE_MOON, SE_MARS, SE_MEAN_NODE],  # subset
    )
    return str(path)

@pytest.fixture
def leb_reader(test_leb_file):
    reader = LEBReader(test_leb_file)
    yield reader
    reader.close()
```

### 9.4 Benchmark Tests

```python
@pytest.mark.slow
class TestBenchmark:
    def test_speedup_vs_skyfield(self, leb_reader, benchmark):
        """Verify >=10x speedup for swe_calc_ut with SEFLG_SPEED."""
        jd = 2451545.0

        # Time the .leb path
        t0 = time.perf_counter()
        for _ in range(1000):
            fast_calc_ut(leb_reader, jd, SE_MARS, SEFLG_SPEED)
        leb_time = (time.perf_counter() - t0) / 1000

        # Time the Skyfield path
        t0 = time.perf_counter()
        for _ in range(1000):
            swe_calc_ut(jd, SE_MARS, SEFLG_SPEED)
        sky_time = (time.perf_counter() - t0) / 1000

        speedup = sky_time / leb_time
        assert speedup >= 10, f"Expected >=10x speedup, got {speedup:.1f}x"
```

---

## 10. Future: Rust Drop-in Replacement

### 10.1 Strategy

Once the Python implementation is validated and stable, port the **reader +
pipeline** to Rust. The generator stays in Python (runs once, not perf-critical).

### 10.2 Package Structure

```
libephemeris-engine/           # Separate package
├── Cargo.toml
├── pyproject.toml             # maturin build config
├── src/
│   ├── lib.rs                 # PyO3 module entry
│   ├── leb.rs                 # .leb mmap reader        (~200 lines)
│   ├── chebyshev.rs           # Clenshaw + derivative   (~80 lines)
│   ├── calc.rs                # Pipeline A/B/C          (~400 lines)
│   ├── coords.rs              # Coordinate transforms   (~300 lines)
│   └── pymodule.rs            # PyO3 bindings           (~200 lines)
│   Total: ~1,200 lines Rust
```

### 10.3 Integration

```python
# In planets.py:
try:
    from libephemeris_engine import calc_ut as _native_calc_ut
    _HAS_ENGINE = True
except ImportError:
    _HAS_ENGINE = False

def swe_calc_ut(tjd_ut, ipl, iflag):
    reader = get_leb_reader()
    if _HAS_ENGINE and reader is not None:
        try:
            return _native_calc_ut(reader.path, tjd_ut, ipl, iflag)  # ~1-5us
        except (KeyError, ValueError):
            pass
    if reader is not None:
        try:
            return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)  # ~25-75us
        except (KeyError, ValueError):
            pass
    return _calc_body(...)  # Skyfield, ~350-1050us
```

Three performance tiers:
1. **Rust engine** (~1-5us) -- if `libephemeris-engine` is installed
2. **Python .leb reader** (~25-75us) -- if `.leb` file is configured
3. **Skyfield** (~350-1050us) -- fallback, always available

### 10.4 Rust Dependencies

```toml
[dependencies]
pyo3 = { version = "0.23", features = ["extension-module"] }
memmap2 = "0.9"
byteorder = "1.5"

# Zero math dependencies: all trig uses f64::sin(), f64::cos(), etc.
```

### 10.5 Rust Port Timeline (after Python is stable)

| Phase | Weeks | Deliverable |
|---|---|---|
| 1 | 1-2 | `leb.rs` + `chebyshev.rs` + basic `calc.rs` |
| 2 | 2-3 | `coords.rs` (precession, nutation, aberration, frame rotations) |
| 3 | 3-4 | `pymodule.rs` (PyO3 bindings) + integration tests |
| 4 | 4-5 | Full validation vs Python .leb reader + benchmarks |

---

## 11. File Inventory

### 11.1 New Files

| File | Phase | LOC est. | Purpose |
|---|---|---|---|
| `libephemeris/leb_format.py` | 0 | ~150 | Format constants, dataclasses, struct helpers |
| `scripts/generate_leb.py` | 1 | ~600 | CLI generator (Skyfield -> .leb) |
| `libephemeris/leb_reader.py` | 2 | ~400 | mmap reader + Clenshaw evaluation |
| `libephemeris/fast_calc.py` | 3 | ~800 | Calculation pipelines A/B/C |
| `tests/test_leb/test_leb_format.py` | 5 | ~200 | Format unit tests |
| `tests/test_leb/test_leb_reader.py` | 5 | ~250 | Reader unit tests |
| `tests/test_leb/test_fast_calc.py` | 5 | ~300 | Pipeline comparison tests |
| `tests/test_leb/test_generate_leb.py` | 5 | ~150 | Generator tests |
| `tests/test_leb/conftest.py` | 5 | ~50 | Test fixtures |
| **Total new** | | **~2,900** | |

### 11.2 Modified Files

| File | Phase | Changes |
|---|---|---|
| `libephemeris/state.py` | 4 | +30 lines: `_LEB_FILE`, `_LEB_READER`, `set_leb_file()`, `get_leb_reader()`, `close()` update |
| `libephemeris/planets.py` | 4 | +10 lines: fast path dispatch in `swe_calc_ut()` and `swe_calc()` |
| `libephemeris/__init__.py` | 4 | +5 lines: export `set_leb_file`, `get_leb_reader` |
| `libephemeris/context.py` | 4 | +20 lines: `set_leb_file()` and `get_leb_reader()` on `EphemerisContext` |

### 11.3 No New Runtime Dependencies

The entire implementation uses only Python stdlib modules:
- `mmap` -- memory-mapped file access
- `struct` -- binary serialization/deserialization
- `math` -- trigonometry, sqrt, atan2, asin
- `os` -- file path operations
- `dataclasses` -- format dataclasses

`numpy` is needed only by the generator script (already a dev dependency).

---

## 12. Risk Register

| # | Risk | Severity | Likelihood | Mitigation |
|---|---|---|---|---|
| 1 | **Numerical parity with Skyfield** -- small differences in light-time iteration, aberration, or precession-nutation could cause >0.001" discrepancies | Medium | Medium | Extensive comparison test suite (1000+ random dates per body). Accept 0.01" tolerance initially, tighten later. |
| 2 | **Gravitational deflection** -- Skyfield's `.apparent()` applies gravitational light deflection by the Sun (~0.004" max). The .leb pipeline omits this. | Low | High | Document the omission. Max error is 0.004" for bodies near the Sun limb, negligible for most use cases. Add deflection correction later if needed. |
| 3 | **Longitude wrap-around in Chebyshev** -- incorrect unwrapping near 0/360 boundary produces wrong fits | Medium | Low | Generator verifies every segment crossing 0/360. Test with bodies known to cross frequently (Moon node, fast-moving bodies). |
| 4 | **Star-based ayanamsha modes** (5 of 43) require Skyfield | Low | Certain | These modes automatically fall back to Skyfield. Document which modes are .leb-accelerated. |
| 5 | **SEFLG_TOPOCTR not supported** in .leb mode | Low | Certain | Falls back to Skyfield automatically. Document the limitation. |
| 6 | **Delta-T interpolation error** -- cubic interpolation between 30-day samples may not capture rapid DUT1 changes | Low | Low | Use denser sampling (15 days) near known discontinuities. Verify against `swe_deltat()` at random dates. |
| 7 | **Outer planet center accuracy** -- generator must use planet center positions (not barycenters) for Jupiter-Pluto | Medium | Low | Generator uses `get_planet_target()` which resolves to `_SpkCenterTarget`. Verify center vs barycenter offset is captured. |
| 8 | **mmap compatibility** -- edge cases on Windows (file locking) or 32-bit systems (address space) | Low | Low | Test on Windows CI. The 16 MB file is well within 32-bit address space limits. |

---

## 13. Reference: Current Pipeline Details

This section documents the exact behavior of the current Skyfield-based pipeline,
so that `fast_calc.py` can replicate it faithfully.

### 13.1 `swe_calc_ut()` Entry Point (planets.py:710-781)

```
1. If ipl == SE_ECL_NUT (-1): return _calc_nutation_obliquity()
2. Strip SEFLG_MOSEPH bit
3. Validate JD range (only for bodies that use JPL ephemeris)
4. Create Skyfield Time: t = ts.ut1_jd(tjd_ut)
5. Call _calc_body(t, ipl, iflag)
6. Wrap SkyfieldRangeError -> EphemerisRangeError
```

### 13.2 `_calc_body()` Dispatcher (planets.py:1312-2087)

Routes by body type (in priority order):
1. Planetary moons -> `planetary_moons.calc_moon_position()`
2. Lunar nodes (SE_MEAN_NODE, SE_TRUE_NODE) -> `lunar.calc_*()`
3. South nodes (-SE_MEAN_NODE, -SE_TRUE_NODE) -> north node + 180 deg
4. Lilith (SE_MEAN_APOG, SE_OSCU_APOG) -> `lunar.calc_*()`
5. Interp. apogee/perigee -> `lunar.calc_interpolated_*()`
6. Uranian planets (40-47) -> `hypothetical.calc_uranian_planet()`
7. Transpluto (48) -> `hypothetical.calc_transpluto()`
8. Minor bodies -> SPK -> auto-download -> ASSIST -> Keplerian fallback
9. Fixed stars -> `fixed_stars.calc_fixed_star_position()`
10. Angles (SE_ASCENDANT, etc.) -> `angles.get_angle_value()`
11. Arabic parts -> `arabic_parts.calc_arabic_part_*()`
12. **Standard planets** (SE_SUN-SE_PLUTO, SE_EARTH) -> full Skyfield pipeline

### 13.3 Standard Planet Pipeline (planets.py:1776-2087)

```
1. Target resolution: get_planet_target(planets, target_name)
   - Inner planets: direct SPK segment (e.g., planets["mercury"], NAIF 199)
   - Outer planets: _SpkCenterTarget (barycenter + center offset from planet_centers.bsp)
     or _CobCorrectedTarget (analytical center-of-body correction)
     or raw barycenter (fallback)

2. Observer selection:
   - SEFLG_HELCTR: observer = planets["sun"]
   - SEFLG_BARYCTR: observer = None (SSB origin)
   - SEFLG_TOPOCTR: observer = earth + topos
   - default: observer = planets["earth"]

3. Position computation:
   - SEFLG_TRUEPOS: geometric (instantaneous), no light-time
   - SEFLG_HELCTR/BARYCTR apparent: manual 3-iteration light-time
   - default (geocentric): observer.at(t).observe(target).apparent()
     (includes light-time + aberration + gravitational deflection)

4. Coordinate extraction:
   - SEFLG_EQUATORIAL + J2000: pos.radec() -> ICRS RA/Dec
   - SEFLG_EQUATORIAL: pos.radec(epoch="date") -> true equator of date
   - SEFLG_J2000 ecliptic: manual rotation with eps=23.4392911 deg
   - default: pos.frame_latlon(ecliptic_frame) -> true ecliptic of date

5. Central-difference velocity (planets.py:2028-2070):
   - dt = 1/86400 day (1 second) for planets, 7e-5 day (~6 sec) for Moon
   - Recursive: _calc_body(t_prev, ipl, flags & ~SPEED & ~SIDEREAL)
   - Recursive: _calc_body(t_next, ipl, flags & ~SPEED & ~SIDEREAL)
   - dp = (result_next - result_prev) / (2 * dt)
   - Wrap handling: if dp > 180/(2*dt), subtract 360/(2*dt)

6. Sidereal correction (planets.py:2072-2085):
   - Only if SEFLG_SIDEREAL and not SEFLG_EQUATORIAL
   - lon -= _get_true_ayanamsa(t.ut1)
   - Speed: dlon -= (aya_next - aya_prev) / (2*dt)

7. Return: _to_native_floats((lon, lat, dist, dlon, dlat, ddist)), iflag
```

### 13.4 Non-Planet Body Handling

For lunar nodes, Lilith, hypotheticals in the current pipeline:

```
1. Call analytical function (e.g., calc_mean_lunar_node(jd_tt))
2. Compute velocity via central difference:
   - dt = 0.5 days for lunar nodes/Lilith
   - result_prev = calc_*(jd_tt - dt)
   - result_next = calc_*(jd_tt + dt)
   - speed = (next - prev) / (2 * dt) with wrap handling
3. If SEFLG_SIDEREAL and not SEFLG_EQUATORIAL: lon -= ayanamsa
4. If SEFLG_EQUATORIAL: _maybe_equatorial_convert()
   - Uses cotrans_sp() with true obliquity from cache.get_true_obliquity()
   - Or J2000 obliquity 23.4392911 if SEFLG_J2000
5. Return _to_native_floats(result), iflag
```

### 13.5 Key Constants

```python
C_LIGHT_AU_DAY = 173.1446326846693  # Speed of light in AU/day
J2000 = 2451545.0                    # J2000.0 epoch in JD
OBLIQUITY_J2000 = 23.4392911         # Mean obliquity at J2000 (degrees)

# IAU 2006 precession polynomial (arcseconds/century)
PREC_C1 = 5028.796195
PREC_C2 = 1.1054348
PREC_C3 = 0.00007964
PREC_C4 = -0.000023857
PREC_C5 = -0.0000000383

# IAU 2006 mean obliquity polynomial (arcseconds)
OBLIQUITY_COEFFS = (84381.406, -46.836769, -0.0001831, 0.00200340, -0.000000576, -0.0000000434)
```

### 13.6 Skyfield-Dependent Lunar Functions

| Function | Uses Skyfield? | What it gets from Skyfield |
|---|---|---|
| `calc_mean_lunar_node(jd_tt)` | No | Pure polynomial |
| `calc_true_lunar_node(jd_tt)` | **Yes** | Moon position+velocity via `(moon - earth).at(t)` in ecliptic frame |
| `calc_mean_lilith(jd_tt)` | No | Pure polynomial |
| `calc_mean_lilith_with_latitude(jd_tt)` | No | Pure polynomial (calls mean_lilith + mean_node) |
| `calc_true_lilith(jd_tt)` | **Yes** | Moon position+velocity via `(moon - earth).at(t)` in ecliptic frame |
| `calc_interpolated_apogee(jd_tt)` | **Partial** | Longitude is pure math; latitude+distance from `calc_true_lilith()` |
| `calc_interpolated_perigee(jd_tt)` | **Partial** | Longitude is pure math; latitude+distance from `calc_osculating_perigee()` |

All Skyfield-dependent lunar functions share the same access pattern:
`(moon - earth).at(t).frame_xyz(ecliptic_frame)` + `.frame_xyz_and_velocity()`.
Precomputing the final ecliptic coordinates in the `.leb` file eliminates
all Skyfield calls for these bodies.
