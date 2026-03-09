# LEB File Compression Plan

## Context

The LEB (LibEphemeris Binary) format stores precomputed Chebyshev polynomial
coefficients for 31 astronomical bodies. Each body's trajectory is broken into
fixed-length time segments, and for each segment 3 coordinates (x, y, z in
ICRS barycentric) are approximated by Chebyshev polynomials of degree 12-14.

Current storage uses raw `float64` (8 bytes per coefficient), with no
compression. This results in large files:

| Tier     | Time range         | File size |
|----------|--------------------|-----------|
| base     | 1849-2150 (301 yr) | 107 MB    |
| medium   | 1550-2650 (1100 yr)| 377 MB    |
| extended | -13200 to +17191   | 2.8 GB    |

**98.3% of the file is Chebyshev coefficient data** (105 MB out of 107 MB for
base tier). The remaining 1.7% is header, body index, delta-T table, nutation
data, and fixed stars.

### Coefficient characteristics

Analysis of the base tier reveals properties that make compression very
effective:

- **13,766,034 total coefficients** across all bodies and segments.
- **Steep amplitude decay:** Chebyshev coefficients decrease rapidly with
  order. The first few coefficients of each coordinate carry the bulk of the
  signal (position in AU), while higher-order coefficients encode fine detail
  and are many orders of magnitude smaller.
- **45.3% of coefficients have |value| < 1e-12**, contributing less than
  0.000001" to the final position. These are effectively noise at our
  precision target.
- **57.5% have |value| < 1e-10.**
- **Dynamic range spans ~16 orders of magnitude** (from ~1e+0 AU down to
  ~1e-16 AU).

Storing all of these as float64 is wasteful: the high-order coefficients need
far fewer bits of precision than the dominant low-order ones.

### Precision constraint

All 31 bodies must maintain **< 0.001 arcsecond** accuracy versus the Skyfield
reference. The current worst case is Moon at ~0.000332". Any compression must
not degrade this.

---

## Proposed Strategy

Three complementary techniques, applied in order during generation:

### 1. Truncation of negligible coefficients

**Idea:** For each coordinate of each segment, determine how many trailing
coefficients can be dropped without exceeding the precision budget.

Currently every coordinate stores exactly `ncoe` coefficients (typically 13 or
15). But the last several are often < 1e-13 in magnitude and contribute
sub-micro-arcsecond to the result. We store only the number of significant
coefficients, and the reader treats the rest as zero.

**Implementation:**
- During generation, after fitting Chebyshev coefficients, walk backwards from
  the last coefficient and drop any with |value| below a body-specific
  threshold (e.g. 1e-14 for planets, 1e-15 for Moon).
- Store the actual count of non-zero coefficients (`ncoeff_stored`) in the
  segment header. The reader zero-fills the rest.
- The Clenshaw evaluation loop already handles variable coefficient counts
  since `neval` can differ from `ncoe`.

**Estimated impact:** 20-30% fewer coefficients to store overall.

### 2. Variable-width integer quantization (main compression)

**Idea:** Instead of 8 bytes per coefficient, encode each as a fixed-point
integer using only as many bytes as its magnitude requires.

For each coordinate within a segment:

1. Compute `rmax = max(|coefficients|)` for the coordinate.
2. Quantize each coefficient relative to `rmax`:
   ```
   integer_value = round(coefficient / rmax * SCALE_FACTOR)
   ```
   where `SCALE_FACTOR` is chosen per-width (e.g. 2^31 for 4-byte, 2^23 for
   3-byte, etc.).
3. Assign each coefficient to a width class based on its relative magnitude:

   | Class | Width   | Precision       | Used for                     |
   |-------|---------|-----------------|------------------------------|
   | 0     | 4 bytes | ~9.3 sig digits | Dominant coefficients (c0-c4)|
   | 1     | 3 bytes | ~7.2 sig digits | Medium coefficients (c5-c7)  |
   | 2     | 2 bytes | ~4.8 sig digits | Small coefficients (c8-c10)  |
   | 3     | 1 byte  | ~2.4 sig digits | Tiny coefficients (c11+)     |

4. Store a compact header (2-4 bytes) per coordinate describing the count of
   coefficients at each width, followed by the packed integers.

**Precision analysis:**

The dominant coefficients (class 0) carry the position signal. At 32-bit
precision with `rmax ~ 1 AU`:
```
resolution = 1 AU / 2^31 ~ 7e-10 AU ~ 0.1 m ~ 0.000001"
```
This is 1000x better than our 0.001" target.

For class 1 (24-bit), coefficients that are already 1e-3 to 1e-5 of rmax:
```
their positional contribution ~ rmax * 1e-4 ~ 1e-4 AU
quantization error ~ 1e-4 / 2^23 ~ 1.2e-11 AU ~ 0.000002"
```
Still well within budget. The small and tiny coefficients contribute so little
to the final position that even 8-bit quantization preserves sub-milliarcsecond
accuracy.

**Average bytes per coefficient:** approximately 2.5 (vs 8 currently), for a
**~69% reduction** in coefficient data.

### 3. Block compression (zlib/zstd) on quantized data

**Idea:** After quantization, the integer data has much more regularity than
raw float64. Apply byte-shuffle and generic compression as a final pass.

- **Byte shuffle:** Rearrange N integers so that all MSBs are together, all
  second bytes together, etc. This groups similar-entropy bytes for better
  compression.
- **zlib level 6** or **zstd level 3** on the shuffled blocks.

On raw float64, standard compression only achieves 12-13% reduction. On
quantized integer data, the gains are much larger because:
- The values are smaller and more regular.
- The width classes create natural patterns.
- Byte-shuffled integers compress extremely well.

**Estimated additional reduction:** 20-30% on top of quantization.

---

## Expected Results

| Stage                       | Base tier  | % of original |
|-----------------------------|------------|---------------|
| Current (raw float64)       | 107 MB     | 100%          |
| After truncation            | ~85 MB     | ~79%          |
| After variable quantization | ~30-35 MB  | ~28-33%       |
| After block compression     | ~25-28 MB  | ~23-26%       |

For medium tier (377 MB) → ~90-100 MB.
For extended tier (2.8 GB) → ~650-750 MB.

---

## File Format Changes

### Current segment layout (v1)

```
[body_index entry]
  segment_offset: uint64       -> points to coefficient data

[segment data at offset]
  coefficients: float64[ncoe * 3]    -> 3 coordinates, ncoe coefficients each
```

### Proposed segment layout (v2)

```
[body_index entry]
  segment_offset: uint64       -> points to compressed segment

[compressed segment at offset]
  compressed_size: uint32      -> total bytes of this segment's compressed data
  raw_data_size: uint32        -> uncompressed size (for pre-allocation)
  
  [zlib/zstd compressed block containing:]
    For each of 3 coordinates:
      rmax: float64            -> normalization factor (max |coefficient|)
      ncoeff: uint8            -> number of stored coefficients (after truncation)
      widths_header: 2-4 bytes -> count of coefficients at each width class
      data: variable           -> packed integers at specified widths
```

The header section (magic, version, body index, delta-T, nutation, stars)
remains unchanged except for a version bump.

### Backward compatibility

- Version field in the file header distinguishes v1 (uncompressed) from v2
  (compressed).
- The reader detects the version and uses the appropriate decode path.
- The Clenshaw evaluation is unchanged: after decoding, coefficients are
  float64 arrays as before.

---

## Implementation Plan

### Phase 1: Quantization + truncation in generator

Files to modify:
- `scripts/generate_leb.py` — Add quantization/packing after Chebyshev fitting
- `libephemeris/leb_format.py` — New format version constant, width class definitions

### Phase 2: Decompression in reader

Files to modify:
- `libephemeris/leb_reader.py` — Add dequantization path, version dispatch

The reader flow becomes:
1. Open file, read header, detect version.
2. For a given body+segment: seek to `segment_offset`.
3. Read `compressed_size`, decompress with zlib/zstd.
4. For each coordinate: read `rmax`, `ncoeff`, widths header, unpack integers.
5. Dequantize: `coefficient = integer_value / SCALE_FACTOR * rmax`
6. Zero-fill remaining coefficients up to `ncoe`.
7. Evaluate Clenshaw as before.

Dequantization cost is negligible (~100ns per segment, vs ~1us for Clenshaw).

### Phase 3: Validation

- Regenerate base tier with compression.
- Run full precision test suite against Skyfield reference.
- Verify all 31 bodies remain < 0.001".
- Benchmark read performance (should be within 10% of current).

### Phase 4: Apply to all tiers

- Regenerate medium and extended tiers.
- Run cross-tier precision tests.

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Quantization noise breaches 0.001" for Moon | Moon has the tightest margin (0.000332" fitting error). Use 4-byte quantization for more coefficients on Moon specifically. The noise floor at 32-bit is ~0.000001", so the combined error stays under 0.0004". |
| Read performance degradation | zlib decompression of a ~200-byte segment takes <1us. The main bottleneck remains Clenshaw evaluation. Benchmark to confirm. |
| Complexity in reader | The dequantization is ~30 lines of code. Version dispatch adds ~10 lines. Total reader change is modest. |
| Format migration | Old .leb files (v1) continue to work. New files are v2. No migration needed — just regenerate. |

## Decision Points

Before implementing, choose:

1. **Compression algorithm:** zlib (universally available, Python stdlib) vs
   zstd (better ratio, requires `zstandard` package). Recommendation: **zlib**
   for zero dependencies.

2. **Granularity of compression blocks:** per-segment (finest random access)
   vs per-body (better compression ratio). Recommendation: **per-segment**,
   since the reader already loads one segment at a time.

3. **Quantization width classes:** 4 classes (4/3/2/1 byte) vs 6 classes
   (adding 4-bit and 2-bit). Recommendation: start with **4 classes**, add
   sub-byte packing only if the extra complexity is justified by the gains.
