# LEB Parameter Tuning Plan

## Context

Each body in the LEB file is defined by a tuple of parameters in
`BODY_PARAMS` (`libephemeris/leb_format.py`):

```python
body_id: (segment_days, ncoe, coord_type, components)
```

- **`segment_days`**: length of each time segment in days. Longer segments
  mean fewer segments and less data, but the polynomial must approximate a
  longer arc.
- **`ncoe`**: number of Chebyshev coefficients per coordinate (polynomial
  degree + 1). More coefficients capture finer detail but cost more storage.

Each segment stores `ncoe × components × 8` bytes of coefficient data
(before compression). Total coefficient data for a body is:

```
ceil(total_days / segment_days) × ncoe × 3 × 8 bytes
```

The current parameters were chosen conservatively to guarantee <0.001"
precision for all bodies across all tiers. Many bodies have enormous precision
margins — some 1000x better than the 0.001" target.

## Goal

Tune `segment_days` and `ncoe` per-body to **minimize file size** while
keeping precision **below 0.001"** for all bodies on all tiers. This is
independent of (and complementary to) the binary compression work.

**Estimated impact:** ~50% reduction in raw coefficient data before any
compression is applied. Combined with variable-width quantization, the total
reduction could reach 85-90%.

## Prerequisites

This work should be done **after** the compression format (variable-width
quantization) is implemented and validated. Reasons:

1. Compression provides guaranteed lossless-equivalent precision — it's the
   safe win. Parameter tuning trades precision margin for size and needs
   careful validation.
2. Once compression is in place, the marginal benefit of parameter tuning can
   be measured precisely on the compressed output.
3. The validation pipeline (generate → test) takes significant time. Doing it
   once at the end is more efficient.

## Current Parameters and Margins

Based on the 2000-point precision sweep against Skyfield (base tier):

| Body | Current params | Worst error | Margin vs 0.001" | Category |
|------|---------------|-------------|-------------------|----------|
| Moon | 4d / 13 | 0.000332" | 3x | **Bottleneck** |
| OscuApog | 4d / 15 | ~0.00005" | 20x | Tight |
| Juno | 8d / 13 | ~0.00005" | 21x | Moderate |
| Vesta | 8d / 13 | ~0.00004" | 27x | Moderate |
| Pallas | 8d / 13 | ~0.00004" | 28x | Moderate |
| Ceres | 8d / 13 | ~0.00003" | 38x | Moderate |
| Venus | 16d / 13 | ~0.000007" | 143x | Comfortable |
| Mercury | 16d / 15 | ~0.000004" | 250x | Comfortable |
| IntpApog | 4d / 15 | ~0.000004" | 250x | Comfortable |
| IntpPerg | 4d / 15 | ~0.000004" | 250x | Comfortable |
| Saturn | 32d / 13 | ~0.000005" | 200x | Large |
| Mars | 16d / 13 | ~0.000003" | 333x | Large |
| Jupiter | 32d / 13 | ~0.000002" | 500x | Large |
| Sun | 32d / 13 | ~0.000001" | 1000x | Huge |
| TrueNode | 8d / 13 | ~0.000001" | 1000x | Huge |
| Chiron | 8d / 13 | ~0.000001" | 1000x | Huge |
| Earth | 4d / 13 | ~0.000001" | >>1000x | Huge |
| Uranus | 64d / 13 | ~0.000001" | >>1000x | Huge |
| Neptune | 64d / 13 | <0.000001" | >>1000x | Huge |
| Pluto | 32d / 13 | <0.000001" | >>1000x | Huge |
| MeanNode | 8d / 13 | ~0" | >>1000x | Analytical |
| MeanApog | 8d / 13 | ~0" | >>1000x | Analytical |
| Uranians | 32d / 13 | ~0" | >>1000x | Analytical |

## Tuning Strategy

### Principle

For each body, find the largest `segment_days` and smallest `ncoe` such that
the Chebyshev fitting error stays below **0.0005"** (half the budget, leaving
margin for compression quantization noise and runtime pipeline errors).

### Body categories

#### 1. Analytical bodies (MeanNode, MeanApog, 9 Uranians) — 11 bodies

These are computed from simple analytical formulas (low-order polynomials in
time). Their trajectories are extremely smooth.

**Target:** `segment_days=256`, `ncoe=7` (or even less).

These bodies currently waste 13 coefficients per coordinate to represent what
is essentially a linear or quadratic function. Degree 6 (7 coefficients) is
more than sufficient. Longer segments (256 days) are fine because the
functions change slowly.

**Expected saving:** ~93% per body. Total: ~18 MB → ~1.3 MB.

#### 2. Smooth outer planets (Uranus, Neptune, Pluto) — 3 bodies

System barycenters with very smooth trajectories. Currently have huge margins.

**Target:** `segment_days=64-128`, `ncoe=11`.

**Expected saving:** ~58% per body.

#### 3. Medium outer planets (Jupiter, Saturn) — 2 bodies

Smooth system barycenters but slightly faster orbital motion.

**Target:** `segment_days=64`, `ncoe=11`.

**Expected saving:** ~58% per body.

#### 4. Inner planets (Sun, Venus, Mars, Earth) — 4 bodies

Faster orbital motion, but still large margins (>100x).

**Target:** `segment_days=32-64`, `ncoe=11`.

Earth is special — currently at 4d/13 with >>1000x margin. Can likely go to
8d/11 or even 16d/11.

**Expected saving:** 40-58% per body.

#### 5. Fast inner planets (Mercury) — 1 body

Fastest orbital motion among planets, but still 250x margin.

**Target:** `segment_days=32`, `ncoe=13`.

**Expected saving:** ~57%.

#### 6. Asteroids (Chiron, Ceres, Pallas, Juno, Vesta) — 5 bodies

Moderate margins (21-1000x). Chiron has huge margin, the main belt asteroids
are tighter.

**Target:** `segment_days=16`, `ncoe=11`.

**Expected saving:** ~58%.

#### 7. Ecliptic bodies (TrueNode, OscuApog, IntpApog, IntpPerg) — 4 bodies

Mixed margins. OscuApog has fast oscillation (~2.6 deg/day), needs short
segments. IntpApog/IntpPerg have 250x margin.

**Targets:**
- TrueNode: `8d/11` (1000x margin)
- OscuApog: `4d/13` (20x margin, keep conservative)
- IntpApog: `8d/13` (250x margin, can relax from 4d/15)
- IntpPerg: `8d/13` (250x margin, can relax from 4d/15)

**Expected saving:** 15-57% per body.

#### 8. Moon — 1 body

**Do not touch.** Only 3x margin. The Moon is the precision bottleneck and
any relaxation risks breaching 0.001".

## Size Projections (Base Tier)

| Category | Bodies | Current | After tuning | Saving |
|----------|--------|---------|--------------|--------|
| Analytical | 11 | 18.2 MB | 1.3 MB | 93% |
| Smooth outer | 3 | 2.1 MB | 0.9 MB | 58% |
| Medium outer | 2 | 2.1 MB | 0.9 MB | 58% |
| Inner planets | 4 | 13.9 MB | 5.9 MB | 58% |
| Mercury | 1 | 2.5 MB | 1.1 MB | 57% |
| Asteroids | 5 | 21.3 MB | 9.0 MB | 58% |
| Ecliptic | 4 | 28.3 MB | 20.7 MB | 27% |
| Moon | 1 | 8.6 MB | 8.6 MB | 0% |
| **Total** | **31** | **~103 MB** | **~48 MB** | **~53%** |

Adding ~2 MB for header, nutation, delta-T, and stars: **~107 MB → ~50 MB**.

## Validation Procedure

For each body, after changing parameters:

1. **Regenerate the base tier** with the new parameters.
2. **Run the precision test** for that specific body:
   ```bash
   pytest tests/test_leb/compare/base/test_base_<body>.py -v
   ```
   Or use `measure_precision.py` for a dense sweep:
   ```bash
   python scripts/measure_precision.py --tier base --body <body_id> --samples 2000
   ```
3. **Verify worst-case error < 0.0005"** (half the budget).
4. **Repeat for medium and extended tiers** — longer time ranges may have
   different worst-case segments.

### Iterative approach

Tune one body category at a time, from largest margin to smallest:

1. Analytical bodies (easiest, biggest win: 17 MB saved)
2. Smooth outer planets
3. Inner planets + Earth
4. Asteroids
5. Ecliptic bodies (most complex, smallest win per-body)

Skip Moon entirely.

### Rollback criterion

If any body's worst-case error exceeds 0.0005" with the proposed parameters,
try an intermediate setting (e.g. keep ncoe=13 but extend segment_days, or
vice versa). If no intermediate works within the precision budget, keep the
original parameters for that body.

## Cross-Tier Considerations

Parameters in `BODY_PARAMS` apply to **all tiers**. A setting that works for
base tier (301 years) may not work for extended tier (30,000+ years) because:

- Outer planet orbits have longer-period perturbations that may require more
  coefficients over longer time spans.
- Ecliptic obliquity changes significantly over millennia, affecting ecliptic
  coordinate bodies.

The extended tier already uses a relaxed `ECLIPTIC_ARCSEC=0.1` tolerance for
ecliptic bodies. Parameter tuning should be validated on all three tiers, but
the extended tier may need body-specific overrides (a future enhancement not
currently supported by the format).

## Relationship to Compression

This optimization is **multiplicative** with binary compression:

| Approach | Base tier |
|----------|-----------|
| Current | 107 MB |
| Parameter tuning only | ~50 MB |
| Compression only (quantization + zlib) | ~25-28 MB |
| **Both combined** | **~12-15 MB** |

The two should be implemented and validated independently, then combined.
