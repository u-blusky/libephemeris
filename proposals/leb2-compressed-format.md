# LEB2: Compressed Binary Ephemeris Format

**Status:** IMPLEMENTED — shipped in v1.0.0a2
**Date:** 2026-03-26
**Author:** Giacomo Battaglia
**Note:** Original proposal. Compression ratio estimates (3:1) were incorrect — actual is 1.2x lossless, 4-10x with mantissa truncation. See `release-notes/v1.0.0a2.md` for final implementation details.

---

## Table of Contents

1. [Motivation](#motivation)
2. [Current State: LEB1 Format](#current-state-leb1-format)
3. [Compression Technique: Byte Shuffling + zstd](#compression-technique-byte-shuffling--zstd)
4. [LEB2 Format Specification](#leb2-format-specification)
5. [Reader Architecture](#reader-architecture)
6. [Integration Strategy](#integration-strategy)
7. [Implementation Plan](#implementation-plan)
8. [Expected Results](#expected-results)
9. [Trade-offs](#trade-offs)
10. [What Does NOT Change](#what-does-not-change)
11. [File Change Summary](#file-change-summary)
12. [Open Items](#open-items)

---

## Motivation

The goal is to **ship `.leb` files inside the PyPI package** so users get
precomputed ephemeris out of the box without a separate download step.

PyPI enforces a **100 MB per-file limit** (can be raised by manual request,
but it's not guaranteed). Current LEB1 file sizes:

| Tier | Merged File | Size | Fits PyPI (100 MB)? |
|------|-------------|------|---------------------|
| Base | `ephemeris_base.leb` | **107 MB** | No (7 MB over) |
| Medium | `ephemeris_medium.leb` | **377 MB** | No |
| Extended | `ephemeris_extended.leb` | **2.8 GB** | No |

Group files for reference:

| Group File | Size |
|------------|------|
| `ephemeris_base_planets.leb` | 31 MB |
| `ephemeris_base_asteroids.leb` | 13 MB |
| `ephemeris_medium_planets.leb` | 117 MB |
| `ephemeris_medium_asteroids.leb` | 73 MB |
| `ephemeris_extended_planets.leb` | 1.0 GB |
| `ephemeris_extended_asteroids.leb` | 127 MB |

Even the base tier — the smallest merged file covering 1849-2150 — exceeds
the limit. With ~3:1 compression it drops to ~35 MB, comfortably under the
threshold and practical for pip install.

---

## Current State: LEB1 Format

### File Layout

```
Offset      Content                          Size
-------------------------------------------------------
0x0000      File Header                      64 bytes
0x0040      Section Directory                N x 24 bytes (N=5 currently)
variable    Section 0: Body Index            body_count x 52 bytes
variable    Section 1: Chebyshev Data        variable (bulk of the file)
variable    Section 2: Nutation Data         header(40B) + segments
variable    Section 3: Delta-T Table         header(8B) + entries(16B each)
variable    Section 4: Star Catalog          n_stars x 64 bytes
(reserved)  Section 5: Orbital Elements      not yet used
```

### File Header (64 bytes)

Struct format: `<4sIIIdddI20s`

| Field | Type | Description |
|-------|------|-------------|
| `magic` | `char[4]` | `b"LEB1"` |
| `version` | `uint32` | `1` |
| `section_count` | `uint32` | Number of sections (currently 5) |
| `body_count` | `uint32` | Number of bodies in the file |
| `jd_start` | `float64` | Start of date coverage (JD TT) |
| `jd_end` | `float64` | End of date coverage (JD TT) |
| `generation_epoch` | `float64` | JD when file was generated |
| `flags` | `uint32` | Reserved (0) |
| `reserved` | `char[20]` | Padding to 64 bytes |

### Section Directory Entry (24 bytes)

Struct format: `<IIQQ`

| Field | Type | Description |
|-------|------|-------------|
| `section_id` | `uint32` | 0-5 (see section constants) |
| `reserved` | `uint32` | Padding |
| `offset` | `uint64` | Absolute byte offset from file start |
| `size` | `uint64` | Section size in bytes |

### Section Type Constants (from `leb_format.py:42-47`)

```python
SECTION_BODY_INDEX      = 0   # Body index entries (one per body)
SECTION_CHEBYSHEV       = 1   # Bulk Chebyshev coefficient data
SECTION_NUTATION        = 2   # Nutation Chebyshev data (dpsi, deps)
SECTION_DELTA_T         = 3   # Delta-T sparse table (TT - UT1)
SECTION_STARS           = 4   # Fixed star catalog
SECTION_ORBITAL_ELEMENTS = 5  # Reserved, not yet used
```

### Body Index Entry (52 bytes)

Struct format: `<iIIdddIIQ`

| Field | Type | Description |
|-------|------|-------------|
| `body_id` | `int32` | SE_* constant (0=Sun, 1=Moon, ...) |
| `coord_type` | `uint32` | 0=ICRS_BARY, 1=ECLIPTIC, 2=HELIO_ECL, 4=ICRS_BARY_SYSTEM |
| `segment_count` | `uint32` | Number of Chebyshev segments |
| `jd_start` | `float64` | Body coverage start (JD TT) |
| `jd_end` | `float64` | Body coverage end (JD TT) |
| `interval_days` | `float64` | Segment width in days |
| `degree` | `uint32` | Chebyshev polynomial degree |
| `components` | `uint32` | Always 3 (x/y/z or lon/lat/dist) |
| `data_offset` | `uint64` | Absolute byte offset to first coefficient |

### Chebyshev Coefficient Layout

Coefficients are stored as contiguous `float64` arrays in **component-major**
order. For a body with `degree=13` and `components=3`, one segment is:

```
[c0_x, c1_x, ..., c13_x,      <- 14 coefficients for component 0
 c0_y, c1_y, ..., c13_y,      <- 14 coefficients for component 1
 c0_z, c1_z, ..., c13_z]      <- 14 coefficients for component 2
```

One segment = `components * (degree + 1) * 8` bytes (e.g., 3 * 14 * 8 = 336 bytes).

Segments for a given body are stored contiguously. The body's `data_offset`
points to the first byte of segment 0.

### Body Parameters (single source of truth, `leb_format.py:151-194`)

| Body ID | Name | Interval (days) | Degree | Coord Type | Segment Size |
|---------|------|-----------------|--------|------------|-------------|
| 0 | Sun | 32 | 13 | ICRS_BARY | 336 B |
| 1 | Moon | 4 | 13 | ICRS_BARY | 336 B |
| 2 | Mercury | 16 | 15 | ICRS_BARY | 384 B |
| 3 | Venus | 16 | 13 | ICRS_BARY | 336 B |
| 4 | Mars | 16 | 13 | ICRS_BARY | 336 B |
| 5 | Jupiter | 32 | 13 | ICRS_BARY_SYSTEM | 336 B |
| 6 | Saturn | 32 | 13 | ICRS_BARY_SYSTEM | 336 B |
| 7 | Uranus | 64 | 13 | ICRS_BARY_SYSTEM | 336 B |
| 8 | Neptune | 64 | 13 | ICRS_BARY_SYSTEM | 336 B |
| 9 | Pluto | 64 | 11 | ICRS_BARY_SYSTEM | 288 B |
| 10 | Mean Node | 8 | 13 | ECLIPTIC | 336 B |
| 11 | True Node | 8 | 13 | ECLIPTIC | 336 B |
| 12 | Mean Apogee | 8 | 13 | ECLIPTIC | 336 B |
| 13 | Oscu Apogee | 4 | 15 | ECLIPTIC | 384 B |
| 14 | Earth | 4 | 13 | ICRS_BARY | 336 B |
| 15 | Chiron | 8 | 13 | ICRS_BARY | 336 B |
| 17-20 | Ceres-Vesta | 8 | 13 | ICRS_BARY | 336 B |
| 21 | Interp Apogee | 4 | 15 | ECLIPTIC | 384 B |
| 22 | Interp Perigee | 4 | 15 | ECLIPTIC | 384 B |
| 40-48 | Uranians/Transpluto | 256 | 7 | HELIO_ECL | 192 B |

**Total: 31 bodies.**

### Reader Hot Path (`leb_reader.py`, ~1.5 us/call)

```python
def eval_body(self, body_id, jd):
    body = self._bodies[body_id]

    # 1. O(1) segment lookup via integer division
    seg_idx = int((jd - body.jd_start) / body.interval_days)
    seg_idx = max(0, min(seg_idx, body.segment_count - 1))

    # 2. Map JD to normalized tau in [-1, 1]
    seg_start = body.jd_start + seg_idx * body.interval_days
    seg_mid = seg_start + 0.5 * body.interval_days
    tau = 2.0 * (jd - seg_mid) / body.interval_days  # clamped to [-1, 1]

    # 3. Zero-copy coefficient read from mmap
    byte_offset = body.data_offset + seg_idx * n_coeffs * 8
    coeffs = struct.unpack_from(f"<{n_coeffs}d", self._mm, byte_offset)

    # 4. Clenshaw evaluation per component (value + analytical derivative)
    for c in range(3):
        comp_coeffs = coeffs[c * deg1 : (c + 1) * deg1]
        val, deriv = _clenshaw_with_derivative(comp_coeffs, tau)
        pos.append(val)
        vel.append(deriv * 2.0 / body.interval_days)

    # 5. Longitude wrapping for ecliptic bodies
    if coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL):
        pos[0] = pos[0] % 360.0

    return tuple(pos), tuple(vel)
```

The key property: **zero-copy reads**. `struct.unpack_from()` reads directly
from the mmap buffer with no intermediate allocations.

### Auxiliary Sections

**Nutation** (`SECTION_NUTATION`): 40-byte header + segments of 272 bytes each
(2 components * 17 coefficients * 8 bytes). Used on nearly every calculation call.

**Delta-T** (`SECTION_DELTA_T`): 8-byte header + entries of 16 bytes each
(jd, delta_t_days pairs). Sampled every 30 days, linearly interpolated.

**Stars** (`SECTION_STARS`): 64 bytes per star. Fixed star catalog with J2000
positions, proper motion, parallax, radial velocity, magnitude.

All three are small (< 1 MB even for the extended tier) and not worth compressing.

---

## Compression Technique: Byte Shuffling + zstd

### The Problem with Naive Compression of Float64

Chebyshev coefficients are `float64` values. In IEEE 754, a float64 is
8 bytes: 1 sign bit, 11 exponent bits, 52 mantissa bits. Adjacent coefficients
in a polynomial may have vastly different mantissa patterns but similar
exponent bytes. Naive compression (zstd on raw bytes) achieves only ~1.5:1
because the mantissa bytes appear essentially random.

### Byte Shuffling

Byte shuffling reorganizes the byte layout to group corresponding byte
positions across multiple floats:

```
Before (raw layout):
  float 0: [B0 B1 B2 B3 B4 B5 B6 B7]
  float 1: [B0 B1 B2 B3 B4 B5 B6 B7]
  float 2: [B0 B1 B2 B3 B4 B5 B6 B7]
  ...

After (shuffled layout):
  Lane 0:  [B0_f0, B0_f1, B0_f2, ...]   <- all byte-0s together
  Lane 1:  [B1_f0, B1_f1, B1_f2, ...]   <- all byte-1s together
  ...
  Lane 7:  [B7_f0, B7_f1, B7_f2, ...]   <- all byte-7s together
```

This is equivalent to treating the data as an (N x 8) matrix and transposing
it to (8 x N). The exponent bytes (lanes 6-7 in little-endian) become highly
repetitive sequences; higher mantissa bytes also cluster. This gives zstd
dramatically more compressible input.

This is the same technique used by HDF5's shuffle filter and Blosc.

### Implementation

```python
import numpy as np

def shuffle_bytes(data: bytes, element_size: int = 8) -> bytes:
    """Transpose byte lanes: (N floats x 8 bytes) -> (8 lanes x N bytes)."""
    arr = np.frombuffer(data, dtype=np.uint8).reshape(-1, element_size)
    return arr.T.ascontiguousarray().tobytes()

def unshuffle_bytes(data: bytes, element_size: int = 8) -> bytes:
    """Inverse transpose: (8 lanes x N bytes) -> (N floats x 8 bytes)."""
    n_elements = len(data) // element_size
    arr = np.frombuffer(data, dtype=np.uint8).reshape(element_size, n_elements)
    return arr.T.ascontiguousarray().tobytes()
```

### Why zstd

- **Decompression speed**: 5-10 GB/s on modern CPUs
- **High compression at level 19**: Best practical ratio for offline generation
- **Widely deployed**: Used by Linux kernel, Chrome, npm, etc.
- **Pure C with excellent Python bindings**: `zstandard` package (~200 KB wheel)
- **No dictionary needed**: Per-body blocks are large enough for good ratios

### Expected Ratios

With byte shuffling + zstd level 19, typical scientific float64 data achieves:

| Scenario | Raw-only zstd | Shuffle + zstd |
|----------|--------------|----------------|
| Smooth trajectories (outer planets) | ~1.5:1 | ~3-4:1 |
| High-frequency bodies (Moon, Earth) | ~1.3:1 | ~2.5-3:1 |
| Weighted average across all bodies | ~1.4:1 | **~3:1** |

---

## LEB2 Format Specification

### File Layout

```
LEB2 File:
+-- File Header (64 bytes)
|     magic = b"LEB2"
|     version = 1
|     compression_method field in flags byte
|
+-- Section Directory (N x 24 bytes)
|     Same structure as LEB1
|
+-- Section 0: Body Index
|     CompressedBodyEntry x body_count (68 bytes each)
|     Same fields as LEB1 BodyEntry (52 bytes) plus:
|       compressed_size: uint64
|       uncompressed_size: uint64
|
+-- Section 6: Compressed Chebyshev Blocks
|     Body 0: [zstd(shuffle(raw_coefficients))]
|     Body 1: [zstd(shuffle(raw_coefficients))]
|     ...
|     Each blob is independently decompressible
|
+-- Section 2: Nutation (uncompressed, same as LEB1)
+-- Section 3: Delta-T (uncompressed, same as LEB1)
+-- Section 4: Star Catalog (uncompressed, same as LEB1)
```

### New Constants

```python
LEB2_MAGIC = b"LEB2"
LEB2_VERSION = 1

SECTION_COMPRESSED_CHEBYSHEV = 6  # New section ID

# Compression methods (stored in header flags)
COMPRESSION_NONE = 0
COMPRESSION_ZSTD_SHUFFLE = 1

# CompressedBodyEntry: original 52 bytes + compressed/uncompressed sizes
COMPRESSED_BODY_ENTRY_FMT = "<iIIdddIIQQQ"
COMPRESSED_BODY_ENTRY_SIZE = 68  # 52 + 16
```

### CompressedBodyEntry (68 bytes)

Struct format: `<iIIdddIIQQQ`

| Field | Type | Description |
|-------|------|-------------|
| `body_id` | `int32` | Same as LEB1 |
| `coord_type` | `uint32` | Same as LEB1 |
| `segment_count` | `uint32` | Same as LEB1 |
| `jd_start` | `float64` | Same as LEB1 |
| `jd_end` | `float64` | Same as LEB1 |
| `interval_days` | `float64` | Same as LEB1 |
| `degree` | `uint32` | Same as LEB1 |
| `components` | `uint32` | Same as LEB1 |
| `data_offset` | `uint64` | Offset to **compressed** blob |
| `compressed_size` | `uint64` | **New**: size of compressed blob |
| `uncompressed_size` | `uint64` | **New**: size of raw coefficients |

### Key Design Decisions

1. **Per-body independent blobs**: Each body's coefficients are compressed as
   a single independent block. Only the bodies actually queried at runtime get
   decompressed. A typical astrology calculation queries ~10-15 bodies out of 31.

2. **Only Chebyshev data is compressed**: Nutation (~50 KB), delta-T (a few KB),
   and star catalog (~few KB) remain uncompressed. They're tiny and nutation is
   accessed on nearly every call — not worth the complexity.

3. **`data_offset` points to the compressed blob**: The reader reads
   `compressed_size` bytes starting at `data_offset`, decompresses to get
   `uncompressed_size` bytes, then indexes into segments identically to LEB1.

4. **Compression level 19**: Since generation is an offline one-time step,
   we maximize compression ratio. Level 19 is zstd's highest practical level
   (~10x slower generation than default, but better ratio). Generation speed
   is not user-facing.

5. **`zstandard` is a required dependency**: Not optional. This avoids a
   fallback code path (zlib gives only ~2:1 and no shuffle support). The
   `zstandard` package is well-maintained, ~200 KB wheel, widely used.

---

## Reader Architecture

### LEB2Reader (new class in `leb2_reader.py`)

```python
class LEB2Reader:
    """Reader for compressed .leb files (LEB2 format).

    Provides the same interface as LEBReader. Decompresses each body's
    coefficients lazily on first access and caches the result for
    subsequent evaluations.
    """

    def __init__(self, path: str) -> None:
        self._path = path
        self._file = open(path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._cache: Dict[int, bytes] = {}  # body_id -> decompressed buffer
        self._parse()

    def _decompress_body(self, body_id: int) -> None:
        """Decompress a body's coefficients on first access."""
        entry = self._bodies[body_id]
        compressed = self._mm[entry.data_offset :
                              entry.data_offset + entry.compressed_size]
        self._cache[body_id] = decompress_block(
            compressed, entry.uncompressed_size
        )

    def eval_body(self, body_id: int, jd: float) -> ...:
        if body_id not in self._cache:
            self._decompress_body(body_id)

        body = self._bodies[body_id]
        # From here: IDENTICAL to LEBReader.eval_body
        # O(1) segment lookup
        # struct.unpack_from(self._cache[body_id], offset) instead of mmap
        # Clenshaw evaluation
        ...

    # eval_nutation, delta_t, get_star: read directly from mmap (uncompressed)
    # has_body, jd_range, path, close: same as LEBReader
```

### Hot Path Comparison

```
LEB1:  jd -> O(1) index -> struct.unpack_from(mmap, offset) -> Clenshaw -> result
LEB2:  jd -> O(1) index -> struct.unpack_from(cache[body_id], offset) -> Clenshaw -> result
                            ^^^^^^^^^^^^^^^^^
                            bytes buffer instead of mmap
                            (populated once per body, lazily)
```

After the first access to a body, the hot path is functionally identical.
`struct.unpack_from` works on any buffer implementing the buffer protocol
(`mmap`, `bytes`, `bytearray`, `memoryview`).

### Decompression Cost (first access per body)

For a typical body in the base tier:
- Moon (largest, 4-day intervals): ~27,000 segments * 336 bytes = ~9 MB raw, ~3 MB compressed
- zstd decompresses at ~5 GB/s -> ~0.6 ms for 3 MB
- Subsequent accesses: zero additional cost

### Thread Safety

Same guarantees as LEBReader:
- `_cache` is populated lazily but the dict assignment is atomic in CPython (GIL)
- Once populated, the `bytes` buffer is immutable
- `struct.unpack_from` and Clenshaw use only local variables

### Code Reuse

The Clenshaw functions (`_clenshaw`, `_deriv_coeffs`, `_clenshaw_with_derivative`)
are module-level in `leb_reader.py`. `LEB2Reader` imports and reuses them
directly. No duplication of the evaluation algorithm.

---

## Integration Strategy

### Shared Interface

`fast_calc.py` references `LEBReader` in 14+ function signatures, but only
through `TYPE_CHECKING` (line 56). It calls these methods:

- `reader.eval_body(body_id, jd)` (11 call sites)
- `reader.eval_nutation(jd_tt)` (implicit via `_nutation_from_leb`)
- `reader.delta_t(jd)` (1 call site)
- `reader.has_body(body_id)` (1 call site)
- `reader.get_star(star_id)` (via star lookup)
- `reader.jd_range` (property)
- `reader.path` (property)
- `reader.close()` (cleanup)

Define a `Protocol` for structural typing:

```python
# In leb_reader.py
from typing import Protocol

class LEBReaderLike(Protocol):
    """Protocol for LEB reader interface (LEB1 and LEB2)."""

    def eval_body(self, body_id: int, jd: float
                  ) -> Tuple[Tuple[float, float, float],
                             Tuple[float, float, float]]: ...
    def eval_nutation(self, jd_tt: float) -> Tuple[float, float]: ...
    def delta_t(self, jd: float) -> float: ...
    def has_body(self, body_id: int) -> bool: ...
    def get_star(self, star_id: int) -> StarEntry: ...

    @property
    def jd_range(self) -> Tuple[float, float]: ...
    @property
    def path(self) -> str: ...

    def close(self) -> None: ...
    def __enter__(self) -> "LEBReaderLike": ...
    def __exit__(self, *args) -> None: ...
```

### Factory Function

```python
# In leb_reader.py
def open_leb(path: str) -> Union[LEBReader, "LEB2Reader"]:
    """Open a .leb file, auto-detecting LEB1 vs LEB2 format."""
    with open(path, "rb") as f:
        magic = f.read(4)
    if magic == MAGIC:        # b"LEB1"
        return LEBReader(path)
    elif magic == LEB2_MAGIC:  # b"LEB2"
        from .leb2_reader import LEB2Reader
        return LEB2Reader(path)
    else:
        raise ValueError(f"Unknown LEB format: {magic!r}")
```

### Changes to `state.py` (lines 332-336)

```python
# Before:
from .leb_reader import LEBReader
_LEB_READER = LEBReader(path)

# After:
from .leb_reader import open_leb
_LEB_READER = open_leb(path)
```

### Changes to `context.py` (lines 212-214)

```python
# Before:
from .leb_reader import LEBReader
self._leb_reader = LEBReader(self._leb_file)

# After:
from .leb_reader import open_leb
self._leb_reader = open_leb(self._leb_file)
```

### Changes to `fast_calc.py` (line 56)

```python
# Before:
if TYPE_CHECKING:
    from .leb_reader import LEBReader

# After:
if TYPE_CHECKING:
    from .leb_reader import LEBReaderLike as LEBReader  # alias for compat
```

This is the minimal change. All 14+ type annotations in `fast_calc.py`
continue to say `reader: "LEBReader"` and work correctly via the alias.

---

## Implementation Plan

### Phase 1: Format Constants + Compression Primitives

**Goal:** Self-contained, fully testable foundation.

**Files to edit:**
- `libephemeris/leb_format.py` — Add `LEB2_MAGIC`, `LEB2_VERSION`,
  `SECTION_COMPRESSED_CHEBYSHEV`, `COMPRESSION_ZSTD_SHUFFLE`,
  `CompressedBodyEntry` dataclass, `COMPRESSED_BODY_ENTRY_FMT/SIZE`,
  and `write_compressed_body_entry` / `read_compressed_body_entry` helpers.

**Files to create:**
- `libephemeris/leb_compression.py` — Four functions:
  - `shuffle_bytes(data, element_size=8) -> bytes`
  - `unshuffle_bytes(data, element_size=8) -> bytes`
  - `compress_block(data, level=19) -> bytes` (shuffle + zstd.compress)
  - `decompress_block(data, uncompressed_size) -> bytes` (zstd.decompress + unshuffle)

**Files to create (tests):**
- `tests/test_leb/test_leb_compression.py`:
  - Round-trip `shuffle_bytes` / `unshuffle_bytes` with known data
  - Round-trip `compress_block` / `decompress_block`
  - Verify compressed output is smaller than input for realistic coefficient data
  - Edge cases: empty data, single float, non-multiple-of-8 lengths

**Dependency:**
- `pyproject.toml` — Add `"zstandard>=0.22.0"` to `dependencies`

### Phase 2: LEB2 Writer

**Goal:** Generate LEB2 files from the same Chebyshev fitting pipeline.

**Files to edit:**
- `scripts/generate_leb.py` — Add `--format leb2` CLI flag.
  The generation pipeline is identical up to the serialization step.
  For LEB2:
  1. After computing coefficients for each body, collect all segments into
     one contiguous `bytes` buffer (same as today)
  2. Per body: `compress_block(raw_coeffs, level=19)` -> compressed blob
  3. Write compressed blobs sequentially into `SECTION_COMPRESSED_CHEBYSHEV`
  4. Write `CompressedBodyEntry` per body with `compressed_size` and
     `uncompressed_size`
  5. Nutation, delta-T, stars sections: unchanged (copied verbatim from LEB1
     pipeline)
  6. File header: `magic=b"LEB2"`, `flags` encodes `COMPRESSION_ZSTD_SHUFFLE`

**Testing strategy:**
- Generate LEB2 for base tier with same parameters
- Assert `compressed_size < uncompressed_size` for every body
- Assert file is parseable (header, section directory, body index)
- Print compression ratio per body and total

### Phase 3: LEB2 Reader

**Goal:** Read LEB2 files with lazy per-body decompression.

**Files to create:**
- `libephemeris/leb2_reader.py` — `LEB2Reader` class:
  - Constructor: open file, mmap, parse header/sections/body index
  - `_decompress_body(body_id)`: read compressed blob from mmap, decompress,
    store in `self._cache`
  - `eval_body(body_id, jd)`: check cache, decompress if needed, then
    identical Clenshaw path using `struct.unpack_from(self._cache[body_id], ...)`
  - `eval_nutation`, `delta_t`, `get_star`, `has_body`, `jd_range`, `path`,
    `close`, `__enter__`, `__exit__`: same as `LEBReader` (nutation/delta-T/stars
    are uncompressed, read from mmap directly)
  - Imports `_clenshaw`, `_deriv_coeffs`, `_clenshaw_with_derivative` from
    `leb_reader` — no code duplication

**Files to create (tests):**
- `tests/test_leb/test_leb2_reader.py`:
  - Generate LEB1 and LEB2 from identical parameters (fixture or pre-built)
  - Assert `LEB2Reader.eval_body() == LEBReader.eval_body()` for all 31 bodies
    at multiple JDs across the full date range (**bit-exact**, lossless compression)
  - Assert `eval_nutation()` matches
  - Assert `delta_t()` matches
  - Assert `get_star()` matches
  - Test context manager protocol
  - Test error cases: invalid magic, truncated file, body out of range

### Phase 4: Integration

**Goal:** `set_leb_file()` transparently handles both formats.

**Files to edit:**
- `libephemeris/leb_reader.py`:
  - Add `LEBReaderLike` Protocol class
  - Add `open_leb(path) -> Union[LEBReader, LEB2Reader]` factory function

- `libephemeris/state.py` (line 334):
  ```python
  # from .leb_reader import LEBReader     # remove
  # _LEB_READER = LEBReader(path)         # remove
  from .leb_reader import open_leb         # add
  _LEB_READER = open_leb(path)            # add
  ```

- `libephemeris/context.py` (lines 212-214):
  ```python
  # from .leb_reader import LEBReader     # remove
  # self._leb_reader = LEBReader(...)     # remove
  from .leb_reader import open_leb         # add
  self._leb_reader = open_leb(...)        # add
  ```

- `libephemeris/download.py` (line 904-906):
  ```python
  # from .leb_reader import LEBReader     # remove
  # reader = LEBReader(filepath)          # remove
  from .leb_reader import open_leb         # add
  reader = open_leb(filepath)             # add
  ```

- `libephemeris/fast_calc.py` (line 56):
  ```python
  if TYPE_CHECKING:
      from .leb_reader import LEBReaderLike as LEBReader
  ```

- Update type annotation for `_LEB_READER` in `state.py` (line 166):
  ```python
  _LEB_READER: Optional["LEBReaderLike"] = None
  ```

- Update type annotation for `_leb_reader` in `context.py` (line 98):
  ```python
  self._leb_reader: Optional["LEBReaderLike"] = None
  ```

**Testing:**
- Run `poe test:unit:leb` with a LEB2 file configured via `LIBEPHEMERIS_LEB`
- All existing LEB tests must pass unchanged with LEB1 files
- Add integration test: `set_leb_file("test.leb2")` -> `calc_ut()` produces
  correct results

### Phase 5: Packaging

**Goal:** Ship a pre-built LEB2 file inside the package.

**Files to edit:**
- `pyproject.toml`:
  - Add `"zstandard>=0.22.0"` to `dependencies`
  - Add LEB2 file to package data (if shipping base tier inside the wheel)
- `MANIFEST.in` — Include LEB2 file

**Optional:**
- Add `poe leb:compress` command that takes an existing LEB1 file and produces
  a LEB2 file (useful for users who already have LEB1 files)
- Add `poe leb:benchmark:compression` command that reports size and ratio

---

## Expected Results

| Tier | LEB1 Size | LEB2 Size (est.) | Ratio | PyPI feasible? |
|------|-----------|-------------------|-------|----------------|
| Base | 107 MB | ~30-35 MB | ~3:1 | Yes |
| Medium | 377 MB | ~100-125 MB | ~3:1 | With limit increase request |
| Extended | 2.8 GB | ~700-900 MB | ~3:1 | No (always external download) |

### Performance Impact

| Metric | LEB1 | LEB2 |
|--------|------|------|
| First access to a body | ~0 (mmap page fault) | ~0.2-1 ms (decompress) |
| Subsequent eval_body() | ~1.5 us | ~1.5 us (identical) |
| Memory per body | 0 (mmap, OS-managed) | uncompressed size in heap |
| Total heap for all 31 bodies | ~0 | ~107 MB (same as file, but in RAM) |
| Total heap for typical query (10 bodies) | ~0 | ~30-50 MB |

In practice, a typical astrological chart calculation queries 10-15 bodies.
Only those bodies get decompressed. The rest stay compressed on disk.

---

## Trade-offs

### Pros
- **Base tier fits in PyPI**: 107 MB -> ~35 MB, well under the 100 MB limit
- **Zero install friction**: `pip install libephemeris` gets fast ephemeris
  out of the box
- **Transparent**: Users don't know or care about the format; same API
- **Lossless**: Bit-exact results, since compression is lossless
- **Lazy decompression**: Only touched bodies pay the cost

### Cons
- **Breaks zero-copy mmap**: The current format reads coefficients directly
  from memory-mapped file with no copies. LEB2 requires a decompression
  buffer per body. This is the fundamental architectural cost.
- **New required dependency**: `zstandard` (~200 KB wheel). Well-maintained
  but still a new dependency in the tree.
- **Higher memory usage**: Decompressed bodies live on the heap instead of
  being OS-managed via mmap page cache. For the base tier (107 MB total),
  this is negligible. For extended (2.8 GB), it matters but only the queried
  bodies are decompressed.
- **Slightly higher first-access latency**: ~0.5-1 ms per body on first query.
  Amortized to zero over subsequent calls.
- **Reader complexity**: New module, new format version, protocol type,
  cache management. More code to maintain.

---

## What Does NOT Change

- **Clenshaw evaluation algorithm**: Reused from `leb_reader.py` module-level
  functions. No duplication.
- **`fast_calc.py` pipeline logic**: Zero changes to the 1595-line calculation
  engine. All 14+ functions continue to work via the same `reader.eval_body()`
  interface.
- **LEB1 format**: Fully preserved. Both formats coexist indefinitely.
  `open_leb()` auto-detects via magic bytes.
- **Generation pipeline**: Same Chebyshev fitting, vectorized evaluation,
  error checking. Only the final serialization step differs.
- **Thread safety**: Both readers are thread-safe for the same reasons.
- **Precision**: Lossless compression means bit-exact identical results.
- **API surface**: `set_leb_file()`, `get_leb_reader()`, `swe_calc_ut()`,
  `calc_ut()` — all unchanged.

---

## File Change Summary

| File | Action | Description |
|------|--------|-------------|
| `libephemeris/leb_format.py` | **Edit** | Add LEB2 constants, `CompressedBodyEntry` dataclass, serialization helpers |
| `libephemeris/leb_compression.py` | **New** | `shuffle_bytes`, `unshuffle_bytes`, `compress_block`, `decompress_block` |
| `libephemeris/leb2_reader.py` | **New** | `LEB2Reader` class with lazy per-body decompression |
| `libephemeris/leb_reader.py` | **Edit** | Add `open_leb()` factory function, `LEBReaderLike` Protocol |
| `libephemeris/state.py` | **Edit** | Use `open_leb()` instead of `LEBReader()` (~2 lines) |
| `libephemeris/context.py` | **Edit** | Use `open_leb()` instead of `LEBReader()` (~2 lines) |
| `libephemeris/download.py` | **Edit** | Use `open_leb()` for validation (~1 line) |
| `libephemeris/fast_calc.py` | **Edit** | Update `TYPE_CHECKING` import to Protocol (~1 line) |
| `scripts/generate_leb.py` | **Edit** | Add `--format leb2` flag, compression serialization step |
| `pyproject.toml` | **Edit** | Add `zstandard>=0.22.0` to dependencies |
| `tests/test_leb/test_leb_compression.py` | **New** | Shuffle/unshuffle and compress/decompress round-trip tests |
| `tests/test_leb/test_leb2_reader.py` | **New** | LEB2 reader correctness vs LEB1, edge cases |

---

## Open Items

1. **Actual compression ratio**: The ~3:1 estimate is based on similar
   scientific datasets. Need to measure with real LEB1 coefficient data before
   committing to the approach. If the actual ratio is < 2:1, the base tier
   might still not fit in 100 MB.

2. **Extended tier strategy**: Even with 3:1 compression, the extended tier
   (~900 MB) won't fit in PyPI. It will always need a separate download
   mechanism. LEB2 still helps (faster downloads, less disk), but it's not
   the primary motivation.

3. **LEB1-to-LEB2 converter**: Should we provide a standalone tool that
   converts existing LEB1 files to LEB2 without regenerating? This is
   straightforward (read all bodies, compress, rewrite) and useful for users
   who already have downloaded LEB1 files.

4. **Partial decompression exploration**: Instead of decompressing an entire
   body's coefficients at once, could we compress per-segment or per-group?
   This would reduce memory usage but increase complexity and likely hurt
   compression ratio. Not recommended unless memory is a proven constraint.
