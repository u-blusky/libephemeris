# Performance Analysis & Precomputed Ephemeris (.leb) Design

This document summarizes a comprehensive analysis of the libephemeris codebase,
its performance bottlenecks, and the design for a precomputed binary ephemeris
format (`.leb`) with a Python-first implementation followed by a Rust port via PyO3.

---

## Table of Contents

1. [Codebase Overview](#1-codebase-overview)
2. [Performance Bottleneck Analysis](#2-performance-bottleneck-analysis)
3. [Precomputed Binary Format (.leb)](#3-precomputed-binary-format-leb)
4. [Python-First Implementation Plan](#4-python-first-implementation-plan)
5. [Future Rust Port via PyO3](#5-future-rust-port-via-pyo3)
6. [Appendix: Full Module Inventory](#appendix-full-module-inventory)

---

## 1. Codebase Overview

### Project Metrics

| Metric                      | Value                                  |
|-----------------------------|----------------------------------------|
| Package version             | 0.14.0                                 |
| Python support              | 3.9 - 3.13                             |
| Library source files        | 37 `.py` files in `libephemeris/`      |
| Total library LOC           | ~69,000 lines                          |
| Largest module              | `eclipse.py` (14,240 lines)            |
| Unit test files             | ~150+ in `tests/`                      |
| Comparison test files       | ~120 in `compare_scripts/tests/`       |
| Public API entries (`__all__`)| 519                                  |
| Constants defined           | ~300+ (SE_*, SEFLG_*, SIDM_*, etc.)   |
| Exception classes           | 12 (1 base + 4 categories + 7 specific)|
| Runtime dependencies        | 6                                      |

### Architecture Summary

libephemeris does NOT implement its own orbital integrator. It delegates
all positional astronomy to **Skyfield** (which reads JPL DE440/DE441 SPK
kernels via **jplephem**). The library provides:

- A Swiss Ephemeris-compatible Python API (1:1 with pyswisseph)
- 24 astrological house systems
- 43 sidereal/ayanamsha modes
- Solar/lunar eclipse search and contact point calculation
- Besselian element computation
- Lunar node and Lilith calculations (4 methods)
- Fixed star catalog (102 stars)
- Hypothetical/Uranian body support (8 Hamburg School + Transpluto)
- Minor body/asteroid support with auto-SPK download
- Atmospheric extinction and heliacal visibility models

### Runtime Dependencies

| Dependency     | Purpose                                              |
|----------------|------------------------------------------------------|
| `skyfield`     | Core ephemeris engine (JPL DE440 via Skyfield)       |
| `skyfield-data`| Skyfield data files                                  |
| `astroquery`   | Auto SPK kernel downloads from JPL Horizons          |
| `certifi`      | SSL certificates for HTTPS to JPL                    |
| `spktype21`    | SPK type 21 files from JPL Horizons                  |
| `pyerfa`       | IAU ERFA nutation, obliquity, precession, aberration |

### Module Breakdown by LOC

| Module              | LOC    | Domain                                  |
|---------------------|--------|-----------------------------------------|
| `eclipse.py`        | 14,240 | Solar/lunar eclipses, occultations      |
| `fixed_stars.py`    | 5,224  | 102-star catalog, proper motion         |
| `houses.py`         | 4,928  | 24 house systems                        |
| `planets.py`        | 4,877  | Core planetary calculations             |
| `hypothetical.py`   | 3,651  | Uranian planets, Transpluto             |
| `lunar.py`          | 3,033  | Lunar nodes, Lilith (4 methods)         |
| `minor_bodies.py`   | 2,649  | Asteroids, TNO, Kepler equation         |
| `heliacal.py`       | 2,414  | Heliacal rising/setting                 |
| `crossing.py`       | 1,914  | Sun/Moon/planet crossings               |
| `spk_auto.py`       | 1,879  | Automatic SPK download/cache            |
| `utils.py`          | 1,722  | Angular math, coordinate transforms     |
| `extinction.py`     | 1,623  | Atmospheric extinction model            |
| `state.py`          | 1,492  | Global state management                 |
| `true_node_terms.py`| 1,412  | Lunar node perturbation coefficients    |
| `time_utils.py`     | 1,323  | Julian day, Delta-T, calendar           |
| `spk.py`            | 1,304  | SPK kernel support                      |
| `iers_data.py`      | 1,286  | IERS Delta-T data                       |
| `__init__.py`       | 1,142  | Public API surface (519 exports)        |
| `constants.py`      | 1,071  | All SE_*, SEFLG_*, SIDM_* constants     |
| `exceptions.py`     | 962    | Exception hierarchy                     |
| `context.py`        | 547    | Thread-safe EphemerisContext             |
| `profiling.py`      | 518    | Profiling infrastructure                |
| `astrometry.py`     | 511    | IAU precession/nutation/aberration      |
| `cache.py`          | 155    | LRU cache for nutation/obliquity        |

---

## 2. Performance Bottleneck Analysis

### 2.1 Current Performance Profile

The library is **pure Python** with zero native extensions. No Cython, mypyc,
C extensions, or compiled code exists anywhere in the project.

#### Cost per Operation

| Operation                       | Cost     | Why                                           |
|---------------------------------|----------|-----------------------------------------------|
| `swe_calc_ut()` with `SEFLG_SPEED` | ~1ms  | 3 full Skyfield pipelines (9-12 Chebyshev evals + Python object overhead) |
| `swe_calc_ut()` without speed   | ~350us   | 1 Skyfield pipeline                           |
| `swe_houses()` (Placidus)       | ~500us   | Up to 50 trig iterations + GAST from Skyfield |
| `sol_eclipse_when_glob()`       | ~500ms-2s| 600-1050 Skyfield pipelines (iterative search)|
| Full chart (10 planets + houses) | ~12ms   | 30x `swe_calc_ut` + 1x `swe_houses`          |

### 2.2 Where CPU Time Is Spent (Priority Order)

#### 1. Skyfield `.observe()` calls (dominant cost, ~90%+ of total time)

Each `.observe()` involves light-time iteration (2-3 iterations of JPL
Chebyshev polynomial evaluation via jplephem). This is the single largest
cost factor in every operation.

#### 2. Central difference velocity (3x multiplier)

The velocity computation in `_calc_body()` (planets.py:1909-1910) calls
itself recursively twice for central difference:

```python
result_prev, _ = _calc_body(t_prev, ipl, flags_no_speed)
result_next, _ = _calc_body(t_next, ipl, flags_no_speed)
velocity = (result_next - result_prev) / (2 * dt)
```

This means every `swe_calc_ut()` with `SEFLG_SPEED` costs **3x** the
position-only cost. Each of the 3 calls runs the full Skyfield pipeline:
`.at().observe().apparent()`.

**Total Skyfield calls per `swe_calc_ut()` with SEFLG_SPEED:**

| Step                    | `.at()` | `.observe()` | `.apparent()` |
|-------------------------|---------|--------------|---------------|
| Main position           | 1       | 1            | 1             |
| Velocity t-dt (recursive)| 1     | 1            | 1             |
| Velocity t+dt (recursive)| 1     | 1            | 1             |
| **Total**               | **3**   | **3**        | **3**         |

Each `.observe()` internally iterates light-time (~2-3 iterations), so
the true count of JPL ephemeris segment evaluations is approximately
**9-12 Chebyshev polynomial evaluations** per single `swe_calc_ut()` call.

#### 3. Eclipse search (combinatorial explosion)

`sol_eclipse_when_glob()` call breakdown:

| Component           | `swe_calc_ut()` calls | Direct Skyfield pipelines |
|---------------------|-----------------------|---------------------------|
| Find New Moon       | ~42                   | ~126                      |
| Check eclipse       | ~3                    | ~18                       |
| Refine maximum      | 0                     | ~90-180                   |
| Contact times       | 0                     | ~369-729                  |
| **Total**           | **~45**               | **~600-1050**             |

#### 4. No position caching

When calculating a full chart (10 planets at the same JD), the Earth
observer position `earth.at(t)` is recomputed **30 times** (10 planets x
3 pipelines each) despite being identical every time.

#### 5. `erfa.nut06a()` (mitigated by cache)

Would be expensive (~0.02ms per call with 1365 nutation terms), but the
LRU cache in `cache.py` (128 entries) makes this negligible for same-JD
calculations. Cache hit: ~0.0001ms (200x speedup).

### 2.3 Existing Caching Infrastructure

| Cached function              | Cache size | Underlying call    | Speedup |
|------------------------------|------------|--------------------|---------|
| `get_cached_nutation(jd_tt)` | 128 entries| `erfa.nut06a()`    | 200x    |
| `get_cached_obliquity(jd_tt)`| 128 entries| `erfa.obl06()` + nutation | 200x |

**What is NOT cached:**
- Planet positions (every `swe_calc_ut()` hits Skyfield every time)
- Skyfield Time objects (created fresh each call)
- Earth position for observer (recomputed per planet)

### 2.4 numpy Usage in Hot Paths

All numpy usage in the hot paths is **scalar or small 3-element vector
operations** (e.g., `np.array(offset)` for tuple→array conversion). No
vectorized batch computation (`np.dot`, `np.matmul`, etc.) is used anywhere
in the critical path.

---

## 3. Precomputed Binary Format (.leb)

### 3.1 Core Concept

Instead of calling Skyfield at runtime (slow), precompute all necessary
data using Skyfield once, store as Chebyshev polynomial coefficients in a
compact binary file, and evaluate with a fast Clenshaw algorithm at runtime.

```
             OFFLINE (once)                       RUNTIME (every call)
             ==============                       ====================

Skyfield/DE441 --> Generator --> .leb file --> Reader + Clenshaw eval
(slow, precise)   (minutes)    (16MB-2.4GB)   (~1-25us per call)
```

Key insight: **velocities are free**. The analytical derivative of the
Chebyshev polynomial gives exact velocity without central difference:

```
If P(x) = sum(c_k * T_k(x))          (position)
Then P'(x) = sum(c_k * T'_k(x))      (velocity, T'_k = k * U_{k-1})
```

This eliminates the current 3x overhead from recursive `_calc_body()` calls.

### 3.2 Body Coordinate Categories

Bodies fall into 3 categories with different storage strategies:

| Category          | Bodies                                   | Stored as              | Why                                          |
|-------------------|------------------------------------------|------------------------|----------------------------------------------|
| ICRS barycentric  | Sun, Moon, Mercury-Pluto, Earth, Chiron, Ceres-Vesta | (x,y,z) AU ICRS | Reader applies geocentric, light-time, aberration, precession, nutation at runtime |
| Ecliptic direct   | Mean/True Node, Mean/True/Osculating Lilith, Interpolated Apogee/Perigee | (lon,lat,dist) degrees/AU | Already final coordinates; reader applies sidereal/equatorial if needed |
| Heliocentric ecliptic | Uranians (40-47), Transpluto (48)   | (lon,lat,dist) helio ecliptic | Already in final frame; `_calc_body()` returns heliocentric today too |

**Rationale for ICRS barycentric (planets):** A single dataset of (x,y,z)
barycentric ICRS coordinates supports ALL output combinations (geocentric
ecliptic, equatorial, heliocentric, J2000, sidereal) through coordinate
transforms at runtime. This avoids storing 5x redundant datasets.

**Rationale for ecliptic direct (nodes/Lilith):** These bodies are computed
by analytical formulas (Meeus polynomials, ELP2000 perturbation series,
or Skyfield h=rxv) that already produce ecliptic coordinates. Storing the
final ecliptic values is simpler and the Chebyshev derivative gives direct
dlon/dt velocity without central difference.

### 3.3 Nutation, Delta-T, and Stars

| Data              | Storage method                | Size (200yr) |
|-------------------|-------------------------------|-------------|
| Nutation (dpsi, deps) | Chebyshev segments (interval=32d, degree=16) | ~700 KB |
| Delta-T           | Sparse table (jd, delta_t) every 30 days, cubic interpolation | ~20 KB |
| Star catalog      | 102 records × (id, ra, dec, pm_ra, pm_dec, parallax, rv, mag) | ~8 KB |

Stars do NOT need Chebyshev precomputation. They are 102 records with
J2000 positions and proper motion rates. Proper motion is linear and can
be extrapolated at runtime with 4 arithmetic operations.

### 3.4 File Layout

```
Byte offset  Content
-----------  ----------------------------------------
0x0000       File Header (64 bytes)
0x0040       Section Directory (N x 24 bytes)
variable     Section 0: Body Index
variable     Section 1: Chebyshev Data (bulk, 95%+ of file)
variable     Section 2: Nutation Chebyshev
variable     Section 3: Delta-T Table
variable     Section 4: Star Catalog
variable     Section 5: Orbital Elements (hypotheticals)
```

#### File Header (64 bytes)

| Offset | Type     | Field              |
|--------|----------|--------------------|
| 0      | `[u8;4]` | magic = `b"LEB1"`  |
| 4      | `u32`    | version = 1        |
| 8      | `u32`    | section_count      |
| 12     | `u32`    | body_count         |
| 16     | `f64`    | jd_start           |
| 24     | `f64`    | jd_end             |
| 32     | `f64`    | generation_epoch   |
| 40     | `u32`    | flags              |
| 44     | 20 bytes | reserved           |

#### Body Index Entry (48 bytes per body)

| Offset | Type  | Field          |
|--------|-------|----------------|
| 0      | `i32` | body_id        |
| 4      | `u32` | coord_type     | (0=ICRS_BARY, 1=ECLIPTIC, 2=HELIO_ECL)
| 8      | `u32` | segment_count  |
| 12     | `f64` | jd_start       |
| 20     | `f64` | jd_end         |
| 28     | `f64` | interval_days  |
| 36     | `u32` | degree         |
| 40     | `u32` | components     | (3 for xyz/lonlatdist, 2 for lon+lat, 1 for lon)
| 44     | `u64` | data_offset    |

#### O(1) Segment Lookup

Given a Julian Day `t` and a body `b`, the correct Chebyshev segment is
found with a single integer division (no binary search):

```
idx = floor((t - body.jd_start) / body.interval_days)
offset = body.data_offset + idx * components * (degree+1) * 8
```

### 3.5 Chebyshev Parameters Per Body

| Body       | Interval (days) | Degree | Bytes/segment | Segments (200yr) | Total    |
|------------|-----------------|--------|---------------|-------------------|----------|
| Moon       | 8               | 13     | 336           | 9,125             | 3.1 MB   |
| Mercury    | 16              | 15     | 384           | 4,563             | 1.8 MB   |
| Venus      | 32              | 13     | 336           | 2,282             | 767 KB   |
| Sun (EMB)  | 32              | 13     | 336           | 2,282             | 767 KB   |
| Earth      | 8               | 13     | 336           | 9,125             | 3.1 MB   |
| Mars       | 32              | 13     | 336           | 2,282             | 767 KB   |
| Jupiter    | 64              | 11     | 288           | 1,141             | 329 KB   |
| Saturn     | 64              | 11     | 288           | 1,141             | 329 KB   |
| Uranus     | 128             | 9      | 240           | 571               | 137 KB   |
| Neptune    | 128             | 9      | 240           | 571               | 137 KB   |
| Pluto      | 128             | 9      | 240           | 571               | 137 KB   |
| Chiron     | 32              | 13     | 336           | 2,282             | 767 KB   |
| Ceres-Vesta (4) | 32         | 13     | 336           | 4 x 2,282        | 3.1 MB   |
| Nutation   | 32              | 16     | 272           | 2,282             | 621 KB   |
| **Total (200yr)** |          |        |               |                   | **~16 MB** |

For 30,000 years the file scales linearly to ~2.4 GB.

### 3.6 Longitude Wrap-Around Handling

Ecliptic longitude has a discontinuity at 0/360 degrees. The generator
must handle this for bodies stored in ecliptic coordinates:

1. Before Chebyshev fitting: **unwrap** the longitude series (remove 360-degree jumps)
2. Fit Chebyshev to the unwrapped series
3. After evaluation: **re-wrap** with `degnorm()` (modulo 360)

The generator detects segments that cross 0/360 and verifies the unwrapped
fit produces correct results.

---

## 4. Python-First Implementation Plan

### 4.1 New Files

| File                        | LOC est. | Purpose                                   |
|-----------------------------|----------|-------------------------------------------|
| `libephemeris/leb_format.py`| ~150     | Format constants, struct definitions       |
| `scripts/generate_leb.py`   | ~600     | CLI generator using Skyfield               |
| `libephemeris/leb_reader.py`| ~400     | mmap reader + Clenshaw evaluation          |
| `libephemeris/fast_calc.py` | ~800     | Calculation pipeline using leb_reader      |
| `tests/test_leb_*.py`       | ~650     | Validation and benchmark tests             |
| **Total new code**          | **~2,600**| **Zero new dependencies**                 |

### 4.2 Modified Files

| File                     | Changes                                        |
|--------------------------|------------------------------------------------|
| `libephemeris/state.py`  | Add `set_leb_file()`, `get_leb_reader()`       |
| `libephemeris/planets.py`| Add fallback to `fast_calc` when .leb available|
| `libephemeris/__init__.py`| Export new functions                           |

### 4.3 Phase 1: Format + Generator

**`leb_format.py`** defines:
- `MAGIC = b"LEB1"`, `VERSION = 1`
- Coordinate type constants: `COORD_ICRS_BARY=0`, `COORD_ECLIPTIC=1`, `COORD_HELIO_ECL=2`
- Dataclasses: `FileHeader`, `BodyEntry`, `NutationEntry`, `StarEntry`, `DeltaTEntry`
- Functions: `write_header()`, `write_body_index()`, `read_header()`, `read_body_index()`
- `BODY_PARAMS` dict mapping each body_id to `(interval, degree, coord_type, components)`

**`generate_leb.py`** implements:

1. **`generate_body_icrs(body_id, jd_start, jd_end, params)`** -
   For Sun through Pluto, Earth, Chiron, Ceres-Vesta. Samples barycentric
   ICRS positions at Chebyshev nodes using Skyfield, fits with
   `numpy.polynomial.chebyshev.chebfit()`, verifies error < 0.001 arcsec.

2. **`generate_body_ecliptic(body_id, jd_start, jd_end, params)`** -
   For nodes, Lilith, hypotheticals. Calls the existing libephemeris
   functions (`calc_mean_lunar_node`, `calc_true_lunar_node`,
   `calc_true_lilith`, etc.) at Chebyshev nodes.

3. **`generate_nutation(jd_start, jd_end)`** -
   Samples `erfa.nut06a()` at Chebyshev nodes (interval=32d, degree=16).

4. **`generate_delta_t(jd_start, jd_end)`** -
   Sparse table of `swe_deltat()` every 30 days.

5. **`generate_star_catalog()`** -
   Extracts 102 star records from `STAR_CATALOG` in `fixed_stars.py`.

6. **`verify_leb()`** -
   Post-generation validation vs Skyfield on random dates.

**Estimated generation time (200 years):** ~10 minutes (parallelizable per body).

### 4.4 Phase 2: Reader

**`leb_reader.py`** provides:

```python
class LebReader:
    def __init__(self, path: str)
        # Opens file with mmap.mmap (READ_ONLY)
        # Parses header, body index, nutation/delta-t/star sections

    def eval_body(self, body_id: int, jd: float) -> tuple[tuple, tuple]
        # O(1) lookup: idx = (jd - body.jd_start) / body.interval_days
        # Reads coefficients from mmap (zero-copy via struct.unpack_from)
        # Clenshaw evaluation: position + analytical derivative (velocity)
        # Returns ((x,y,z), (vx,vy,vz)) or ((lon,lat,dist), (dlon,dlat,ddist))

    def eval_nutation(self, jd_tt: float) -> tuple[float, float]
        # Same Chebyshev lookup for (dpsi, deps) in radians

    def delta_t(self, jd: float) -> float
        # Bisect + cubic interpolation in sparse table

    def get_star(self, star_id: int) -> StarEntry
        # Direct lookup in catalog
```

The Clenshaw algorithm is implemented in pure Python (no numpy) for
single-point evaluation, because numpy array creation overhead (~5us)
dominates the ~1.5us Clenshaw loop for degree-13 polynomials.

### 4.5 Phase 3: Fast Calculation Pipeline

**`fast_calc.py`** reimplements the Skyfield pipeline in pure Python:

#### Pipeline A: ICRS Barycentric Bodies (planets, asteroids)

```
1. delta_t = reader.delta_t(jd_ut) -> jd_tt
2. (earth_pos, earth_vel) = reader.eval_body(SE_EARTH, jd_tt)
3. (target_pos, target_vel) = reader.eval_body(ipl, jd_tt)
4. Observer selection (geocentric / heliocentric / barycentric)
5. geo = target - observer (geometric vector)
6. Light-time iteration (3 iterations, ~50 lines):
     dist = |geo|; lt = dist / 173.1446; retarded = eval_body(ipl, jd_tt - lt)
7. Annual aberration (if not SEFLG_NOABERR)
8. Coordinate transform:
   - Ecliptic J2000: fixed obliquity rotation (23.4392911 deg)
   - Ecliptic of date: IAU 2006 precession + nutation from .leb
   - Equatorial of date: precession + nutation (no ecliptic rotation)
9. Sidereal correction (if SEFLG_SIDEREAL): lon -= ayanamsha
10. Velocity: central difference on final coordinates (now ~25us per eval
    instead of ~350us, so 3x overhead is ~75us total instead of ~1050us)
```

#### Pipeline B: Ecliptic Direct Bodies (nodes, Lilith)

```
1. (lon, lat, dist), (dlon, dlat, ddist) = reader.eval_body(ipl, jd_tt)
   # Already in final ecliptic coordinates
2. If SEFLG_EQUATORIAL: cotrans(lon, lat, obliquity)
3. If SEFLG_SIDEREAL: lon -= ayanamsha
4. Velocity from Chebyshev derivative (no central difference needed!)
```

#### Pipeline C: Heliocentric Bodies (hypotheticals)

Identical to Pipeline B. No additional transforms needed.

### 4.6 Phase 4: Integration

Minimal changes to existing code:

```python
# In planets.py swe_calc_ut():
reader = state.get_leb_reader()
if reader is not None:
    try:
        return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)
    except (KeyError, ValueError):
        pass  # body not in .leb, fall through to Skyfield
# ... existing Skyfield code unchanged ...
```

### 4.7 Expected Python-Only Speedup

| Operation                | Skyfield (current) | .leb Python | Speedup    |
|--------------------------|--------------------|-------------|------------|
| `swe_calc_ut()` no speed | ~350us             | ~25us       | **~14x**   |
| `swe_calc_ut()` + speed  | ~1050us            | ~75us       | **~14x**   |
| `swe_calc_ut()` mean node| ~100us             | ~5us        | **~20x**   |
| Eclipse search           | ~1s                | ~100ms      | **~10x**   |
| `swe_houses()` Placidus  | ~500us             | ~400us      | ~1.3x (*)  |

(*) Houses don't benefit much because they're already pure math. The
bottleneck is `t.gast` from Skyfield and obliquity, which still requires
reimplementation of GAST/GMST for full benefit.

### 4.8 Implementation Timeline

| Week | Deliverable                                                    |
|------|----------------------------------------------------------------|
| 1    | `leb_format.py`, `generate_leb.py` (planets only), `leb_reader.py`, validate Sun |
| 2    | Generator (all bodies + nutation + delta-t + stars), `fast_calc.py` Pipeline A |
| 3    | `fast_calc.py` Pipelines B+C, all flags, integration, full test suite + benchmark |

---

## 5. Future Rust Port via PyO3

### 5.1 Strategy

Once the Python implementation is validated, port the **reader + pipeline**
to Rust. The same `.leb` file format works for both Python and Rust. The
generator stays in Python (it runs once, not performance-critical).

### 5.2 Rust Crate Structure

```
engine/
├── Cargo.toml
├── pyproject.toml          # maturin build config
├── src/
│   ├── lib.rs              # PyO3 module entry point
│   ├── leb.rs              # .leb reader, mmap, index        (~200 lines)
│   ├── chebyshev.rs        # Clenshaw eval + derivative      (~80 lines)
│   ├── calc.rs             # Main calc pipeline               (~400 lines)
│   ├── coords.rs           # ICRS/ecliptic/equatorial transforms (~300 lines)
│   ├── precession.rs       # IAU 2006 polynomial              (~150 lines)
│   ├── light_time.rs       # 3-iteration light-time           (~60 lines)
│   ├── aberration.rs       # Annual aberration                (~100 lines)
│   ├── ayanamsha.rs        # 43 sidereal systems              (~800 lines)
│   ├── houses/             # 24 house systems                 (~2000 lines)
│   │   ├── mod.rs
│   │   ├── placidus.rs
│   │   ├── koch.rs
│   │   └── ...
│   ├── eclipse.rs          # Eclipse search + contacts         (~1500 lines)
│   ├── crossing.rs         # Crossings, stations               (~500 lines)
│   ├── lunar.rs            # Lunar nodes, Lilith               (~400 lines)
│   ├── hypothetical.rs     # Kepler equation, hypotheticals    (~300 lines)
│   ├── stars.rs            # Fixed stars + proper motion       (~200 lines)
│   ├── constants.rs        # All SE_*, SEFLG_*                 (~400 lines)
│   ├── error.rs            # Error types                       (~100 lines)
│   └── pymodule.rs         # PyO3 bindings                     (~800 lines)
│
│   Total: ~7,800 lines Rust
```

### 5.3 Rust Dependencies

```toml
[dependencies]
pyo3 = { version = "0.23", features = ["extension-module"] }
memmap2 = "0.9"              # Memory-mapped .leb files
byteorder = "1.5"            # Binary parsing
thiserror = "2"              # Typed errors
parking_lot = "0.12"         # Fast Mutex/RwLock

[dev-dependencies]
approx = "0.5"               # Float comparison in tests
criterion = "0.5"            # Benchmarks
```

Zero math dependencies: all trigonometry uses `f64::sin()`, `f64::cos()`, etc.

### 5.4 Expected Speedup After Rust Port

| Operation                     | Python .leb | Rust .leb  | Total vs Skyfield |
|-------------------------------|-------------|------------|-------------------|
| `swe_calc_ut()` + speed       | ~75us       | ~1-5us     | **200-1000x**     |
| `swe_calc_ut()` mean node     | ~5us        | ~0.1us     | **1000x**         |
| `swe_houses()` Placidus       | ~400us      | ~1-5us     | **100-500x**      |
| Eclipse search                | ~100ms      | ~0.5-2ms   | **500-2000x**     |
| Full chart (10 planets+houses)| ~1ms        | ~20-50us   | **240-600x**      |

### 5.5 Python Integration (Transparent Fallback)

```python
# libephemeris/planets.py
try:
    from libephemeris_engine import Ephemeris as _NativeEphemeris
    _engine = _NativeEphemeris("path/to/ephemeris.leb")
    _HAS_ENGINE = True
except ImportError:
    _HAS_ENGINE = False

def swe_calc_ut(tjd_ut, ipl, iflag):
    if _HAS_ENGINE:
        return _engine.calc_ut(tjd_ut, ipl, iflag)
    reader = state.get_leb_reader()
    if reader is not None:
        return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)
    return _calc_body_skyfield(...)  # original code, unchanged
```

Three tiers of performance:
1. **Rust engine** (~1-5us) - if `libephemeris_engine` is installed
2. **Python .leb reader** (~25-75us) - if .leb file is configured
3. **Skyfield** (~350-1050us) - fallback, always available

### 5.6 Rust Port Timeline (After Python Validation)

| Phase | Weeks | Deliverable                                        |
|-------|-------|----------------------------------------------------|
| 1     | 1-3   | `leb.rs` + `chebyshev.rs` + `calc.rs` (core pipeline) |
| 2     | 3-4   | `coords.rs` + `precession.rs` + `aberration.rs`   |
| 3     | 4-6   | `houses/` (all 24 systems) + `ayanamsha.rs`        |
| 4     | 6-9   | `eclipse.rs` + `crossing.rs`                       |
| 5     | 9-11  | `lunar.rs` + `stars.rs` + `hypothetical.rs`        |
| 6     | 11-13 | `pymodule.rs` (PyO3 bindings) + integration         |
| 7     | 13-15 | Test suite + validation vs pyswisseph + benchmarks  |

---

## Appendix: Full Module Inventory

### A.1 Coordinate Frames by Body Type

This table documents the coordinate frame produced by each body's
calculation function, as found in the source code analysis:

| Body                | Source Function                | Uses Skyfield? | Frame                  | Geocentric? |
|---------------------|--------------------------------|----------------|------------------------|-------------|
| Sun-Pluto           | Skyfield `.at().observe().apparent()` | Yes     | True ecliptic of date  | Yes         |
| `SE_MEAN_NODE`      | `calc_mean_lunar_node()`       | No             | Mean ecliptic of date  | N/A         |
| `SE_TRUE_NODE`      | `calc_true_lunar_node()`       | Yes            | True ecliptic of date  | Yes         |
| `SE_MEAN_APOG`      | `calc_mean_lilith_with_latitude()` | No         | Mean ecliptic of date  | N/A         |
| `SE_OSCU_APOG`      | `calc_true_lilith()`           | Yes            | True ecliptic of date  | Yes         |
| `SE_INTP_APOG`      | `calc_interpolated_apogee()`   | No             | Mean ecliptic of date  | N/A         |
| `SE_INTP_PERG`      | `calc_interpolated_perigee()`  | No             | Mean ecliptic of date  | N/A         |
| Uranians (40-47)     | `calc_uranian_planet()`        | No             | Heliocentric ecliptic  | Helio       |
| Transpluto (48)     | `calc_transpluto()`            | No             | Heliocentric ecliptic  | Helio       |
| Fixed stars         | `calc_fixed_star_position()`   | Yes            | True ecliptic of date  | Yes         |
| Chiron, Ceres-Vesta | SPK kernel via Skyfield        | Yes            | True ecliptic of date  | Yes         |

### A.2 Numerical Algorithms Used

| Algorithm                  | Where Used                            | Parameters                |
|----------------------------|---------------------------------------|---------------------------|
| Clenshaw (Chebyshev eval)  | .leb reader, JPL kernel evaluation    | degree 9-16               |
| Newton-Raphson             | `crossing.py`, `eclipse.py` (New Moon)| tol: 0.001-0.1 arcsec    |
| Brent's method             | `crossing.py` (root finding)          | tol: 1e-6 deg/day        |
| Bisection                  | `eclipse.py` (contact times)          | 60 iterations             |
| Golden section search      | `eclipse.py` (min separation)         | 60 iterations             |
| Fixed-point iteration      | `planets.py` (light-time, 3 iter)     | Convergence in 3 iter     |
| Placidus iteration         | `houses.py`                           | max 50 iter, tol 1e-7 deg|
| Gauquelin iteration        | `houses.py`                           | max 100 iter, tol 1e-8   |
| Pullen SR Newton-Raphson   | `houses.py`                           | max 20 iter, tol 1e-10   |
| Kepler equation            | `minor_bodies.py`, `hypothetical.py`  | Elliptic, hyperbolic, parabolic |
| IAU 2006 precession        | `astrometry.py`, `planets.py`         | 5th degree polynomial     |
| IAU 2006/2000A nutation    | `cache.py` -> `erfa.nut06a()`         | 1365 lunisolar + 687 planetary terms |

### A.3 House Systems Implemented (24)

| Code | Name                    | Method        | Iterative? | Polar safe? |
|------|-------------------------|---------------|------------|-------------|
| P    | Placidus                | Time-based    | Yes (50)   | No          |
| K    | Koch (GOH)              | Birthplace    | No         | No          |
| O    | Porphyry                | Arc trisection| No         | Yes         |
| R    | Regiomontanus           | Equator proj. | No         | Unstable    |
| C    | Campanus                | Prime vertical| No         | Unstable    |
| E/A  | Equal (from Asc)        | Geometric     | No         | Yes         |
| W    | Whole Sign              | Geometric     | No         | Yes         |
| X    | Meridian (Zariel)       | RA-based      | No         | Yes         |
| H    | Horizontal              | Azimuthal     | No         | Unstable    |
| T    | Polich-Page (Topocentric)| Modified Regio| No        | Unstable    |
| B    | Alcabitius              | Arabic        | No         | Unstable    |
| M    | Morinus                 | RA projection | No         | Yes         |
| U    | Krusinski-Pisa          | Multi-step    | No         | Unstable    |
| G    | Gauquelin (36 sectors)  | Iterative     | Yes (100)  | No          |
| V    | Vehlow Equal            | Geometric     | No         | Yes         |
| Y    | APC                     | Parallel circle| No        | Unstable    |
| F    | Carter Poli-Equatorial  | RA-based      | No         | Unstable    |
| N    | Natural Gradient        | Fixed (0=H1)  | No         | Yes         |
| S    | Sripati (Indian)        | Midpoints     | No         | Yes         |
| L    | Pullen SD               | Sinusoidal    | No         | Yes         |
| Q    | Pullen SR               | Newton-Raphson| Yes (20)   | Yes         |
| D    | Equal from MC           | Geometric     | No         | Yes         |
| I/i  | Sunshine (Makransky)    | Sun arc       | No         | Unstable    |

### A.4 Ayanamsha Systems (43)

All 43 systems are implemented in `planets.py::_calc_ayanamsa()` (~375 lines).
They use the IAU 2006 general precession polynomial (5th degree) as the
precession backbone, with system-specific reference epochs and offsets.

Categories:
- **Formula-based** (Lahiri, Fagan-Bradley, Krishnamurti, etc.): `aya = offset + precession_polynomial(T)`
- **Star-based** ("True" modes): compute apparent ecliptic longitude of reference star via Skyfield, subtract target position
- **Galactic-based**: calibrated formulas or galactic pole computation
- **User-defined** (`SE_SIDM_USER`): `aya = ayan_t0 + [p(T_now) - p(T0)]`
