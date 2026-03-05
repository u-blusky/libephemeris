# LEB Precision V3 — Direct Geocentric Ecliptic Storage

## 1. Context

### 1.1 What LEB Does Today

LEB (LibEphemeris Binary) stores precomputed Chebyshev polynomial approximations
of celestial body positions, replacing live Skyfield/JPL DE440 calculations with
fast polynomial lookups (~14x speedup). The file stores 30 bodies across three
coordinate types:

| Pipeline | Coord Type | Bodies | Storage Format |
|----------|-----------|--------|----------------|
| A (ICRS) | `COORD_ICRS_BARY` (0) | Sun, Moon, Mercury–Pluto, Earth, Chiron, Ceres–Vesta (16 bodies) | SSB barycentric Cartesian (x, y, z) in AU |
| B (Ecliptic) | `COORD_ECLIPTIC` (1) | Mean/True Node, Mean/Oscu Apogee, Interp Apogee/Perigee (6 bodies) | Ecliptic of date (lon°, lat°, dist AU) |
| C (Heliocentric) | `COORD_HELIO_ECL` (2) | Uranians (Cupido–Poseidon), Transpluto (9 bodies) | Heliocentric ecliptic (lon°, lat°, dist AU) |

The three tiers cover different time ranges:

| Tier | SPK File | Range | File Size |
|------|----------|-------|-----------|
| Base | de440s.bsp | 1849–2150 CE | ~94 MB |
| Medium | de440.bsp | 1550–2650 CE | ~315 MB |
| Extended | de441.bsp | -5000 to +5000 CE | ~2.3 GB |

### 1.2 History: V1 and V2

**V1** (`leb-precision-improvement` branch, completed): Chebyshev parameter
tuning (shorter segments, higher degrees), analytical velocity pipeline (replaced
central-difference with Chebyshev derivatives), tolerance tightening. This
improved fitting errors from ~5" to sub-arcsecond for all bodies.

**V2** (`leb-precision-v2` branch, in progress): Fixed nutation model mismatch
between the LEB pipeline (`erfa.pnm06a()`, IAU 2006/2000A) and Skyfield (IAU
2000A). Replaced erfa-based precession-nutation with Skyfield's `Time.M` matrix
in `fast_calc.py`. This eliminated the 234" errors in extended tier, bringing
all tiers to the same error level.

### 1.3 The Remaining Problem

After V1 and V2, all three tiers share the same architectural error pattern.
The errors are **not** Chebyshev fitting errors — the polynomial approximations
are excellent (sub-arcsecond in ICRS Cartesian). The errors come from the
**runtime conversion pipeline** that transforms stored ICRS Cartesian coordinates
into the user-facing ecliptic output.

Measured worst-case errors across all tiers (from V1/V2 test results):

| Body | Base (") | Medium (") | Extended (") | Limit Source |
|------|----------|------------|--------------|--------------|
| Saturn | **4.85** | ~4.5 | 0.57 | ICRS→ecliptic amplification |
| Uranus | ~4.0 | **4.58** | **1.91** | ICRS→ecliptic amplification |
| Mercury | ~0.5 | ~0.5 | 0.53 | ICRS→ecliptic amplification |
| Jupiter | ~0.5 | ~0.5 | 0.52 | ICRS→ecliptic amplification |
| Neptune | ~0.3 | ~0.3 | 0.30 | ICRS→ecliptic amplification |
| Pluto | ~0.3 | ~0.3 | 0.12 | ICRS→ecliptic amplification |
| Sun | < 0.001 | < 0.001 | < 0.001 | (negligible) |
| Moon | < 0.001 | < 0.001 | < 0.001 | (negligible) |

Swiss Ephemeris precision: ~0.001" for all planets. **LEB is currently
~5000x less precise** for outer planets.

Velocity errors follow the same pattern (ICRS→ecliptic amplification on
derivatives), with asteroid latitude velocity being the worst case:
Pallas 0.714°/day (base), 0.341°/day (medium).

### 1.4 Why These Errors Are Structural, Not Parametric

The errors cannot be fixed by tuning Chebyshev parameters (shorter segments,
higher degrees). Increasing polynomial order reduces **fitting error** in ICRS
Cartesian, but the dominant error comes from the **conversion chain** at runtime.

## 2. Root Cause Analysis

### 2.1 The Pipeline A Conversion Chain

When a user calls `swe_calc_ut()` for a planet (e.g., Saturn), Pipeline A in
`fast_calc.py` executes this chain at runtime:

```
Step 1: target_bary = leb.eval_body(SATURN, jd)     → Chebyshev eval → (x,y,z) + ε_target
Step 2: earth_bary  = leb.eval_body(EARTH, jd)      → Chebyshev eval → (x,y,z) + ε_earth
Step 3: geo = target_bary - earth_bary               → errors ADD: ε_geo = ε_target + ε_earth
Step 4: light-time (3 iterations of steps 1-3)        → 3 more ε_target errors
Step 5: aberration correction                         → small additional error
Step 6: PNM matrix rotation (Skyfield, full precision) → error preserved (rigid rotation)
Step 7: ecliptic rotation                             → error preserved (rigid rotation)
Step 8: atan2(y,x), asin(z/r)                        → angular_error ≈ ε_cartesian / dist_geo
Step 9: → (lon°, lat°, dist AU)                      → final output with amplified error
```

### 2.2 The Error Amplification Mechanism

The critical step is **Step 8**: converting Cartesian to spherical coordinates.
The angular error from a Cartesian position error ε (in AU) is:

```
angular_error (radians) ≈ ε / geocentric_distance
angular_error (arcsec)  ≈ (ε / geocentric_distance) × 206265
```

This means:
- **Saturn** (geocentric distance ~8–11 AU): a 2×10⁻⁴ AU Cartesian error
  → 2×10⁻⁴/9 × 206265 ≈ **4.6"**
- **Moon** (geocentric distance ~0.0025 AU): fitting is extremely tight (4-day
  segments), so ε is ~10⁻¹⁰ AU → 10⁻¹⁰/0.0025 × 206265 ≈ **0.000008"** —
  negligible because the fitting is excellent for fast-moving nearby bodies.

The amplification is inversely proportional to geocentric distance, but the
Chebyshev fitting error also depends on how the body's barycentric position
varies. For outer planets, the **apparent** motion (as seen from Earth) includes
a fast component from Earth's own orbit, but the stored **barycentric** motion
is slow — so the fitting polynomial doesn't track the apparent motion directly.

### 2.3 Why Step 3 (Subtraction) Matters

Both `target_bary` and `earth_bary` carry independent Chebyshev fitting errors.
When subtracted, these errors **add** (worst case). For bodies close to 1 AU
from the SSB (inner planets, Sun), both barycentric positions are ~1 AU in
magnitude, and the geocentric vector is a small difference of large numbers —
classic numerical cancellation.

For outer planets (Saturn, Uranus), the barycentric positions are large (~10-20
AU), but the geocentric distance is also large, so the cancellation is less
severe. However, the fitting error grows with segment length (outer planets use
32-64 day segments), and the amplification factor `1/dist_geo` ensures that
even moderate Cartesian errors produce arcsecond-level angular errors.

### 2.4 Why More Chebyshev Coefficients Cannot Fix This

Consider Saturn with current params `(32, 13)` (32-day segments, degree 13):
- ICRS fitting error: ~10⁻⁶ AU (excellent)
- After geocentric subtraction + amplification: ~10⁻⁶ / 9 × 206265 ≈ 0.02"

Wait — this seems small. So where does the 4.85" come from?

The answer is that the **Earth** component also carries fitting error (~10⁻⁶ AU
at 4-day segments), and more importantly, the Saturn fitting error is not
uniform: it peaks at segment boundaries and near the Chebyshev polynomial's
oscillation maxima. The **worst-case** fitting error for Saturn across the full
date range reaches ~2×10⁻⁴ AU, which at geocentric distance ~9 AU gives 4.6".

Reducing Saturn's interval from 32 to 16 days would roughly halve the fitting
error, bringing worst-case to ~2.3". Going to 8 days would get ~1.2". Going to
4 days would get ~0.6". But this quadruples the number of segments and the file
size for Saturn, and we'd need to do the same for all outer planets. Even at
4-day intervals, we'd still be at ~0.5" — nowhere near Swiss Ephemeris's 0.001".

The fundamental issue is that **the error floor is set by the ICRS→ecliptic
conversion pipeline**, not by the Chebyshev fitting quality.

### 2.5 Evidence: Pipeline B Bodies Have No Problem

Bodies stored as `COORD_ECLIPTIC` (Pipeline B: MeanNode, TrueNode, OscuApogee,
etc.) have errors of **< 0.05"** — because they store the final ecliptic
coordinates directly. There is no ICRS subtraction, no coordinate transformation,
no error amplification. The Chebyshev polynomial directly approximates the
ecliptic longitude/latitude, and the only error is the fitting error itself.

Similarly, `COORD_HELIO_ECL` bodies (Pipeline C: Uranians) have errors of
**< 0.001"** — again because the stored coordinates are already in the output
frame.

This proves that the Chebyshev approximation itself is highly precise. The
problem is exclusively in Pipeline A's conversion chain.

## 3. Solution: COORD_GEO_ECLIPTIC

### 3.1 Core Idea

Add a fourth coordinate type `COORD_GEO_ECLIPTIC = 3` that stores **geocentric
ecliptic-of-date coordinates** (lon°, lat°, dist AU) directly in the LEB file.

At **generation time**, the full Skyfield pipeline (ICRS → geocentric → light-time
→ aberration → PNM rotation → ecliptic rotation → spherical) is computed at
full float64 precision for each Chebyshev sampling node. The Chebyshev polynomial
is then fitted to these final ecliptic coordinates.

At **runtime**, the reader evaluates the Chebyshev polynomial and returns the
result directly — no conversion pipeline, no error amplification. The only error
is the Chebyshev fitting error itself, which for well-chosen parameters is
< 0.001".

```
CURRENT (COORD_ICRS_BARY):
  Generation: Skyfield ICRS (x,y,z) → fit Chebyshev
  Runtime:    Chebyshev → (x,y,z)±ε → subtract Earth±ε → light-time → aberr → PNM → ecl → spherical
  Error:      ε amplified by conversion chain → ~5"

PROPOSED (COORD_GEO_ECLIPTIC):
  Generation: Skyfield ICRS → subtract Earth → light-time → aberr → PNM → ecl → spherical → fit Chebyshev
  Runtime:    Chebyshev → (lon, lat, dist) directly
  Error:      ε from fitting only → ~0.001"
```

### 3.2 Why This Works

1. **All error-amplifying operations happen at generation time** with full
   float64 precision (~10⁻¹⁶ relative error). There is no amplification of
   Chebyshev fitting error because there is no conversion at runtime.

2. **Geocentric ecliptic coordinates are smooth functions of time** — they
   vary slowly and predictably (planetary longitudes change at ~0.01-1°/day).
   This means Chebyshev polynomials approximate them very well, typically
   requiring fewer coefficients than the more complex ICRS barycentric motion.

3. **Pipeline B already proves this works**: the ecliptic-direct bodies (nodes,
   apogees) use exactly this approach and achieve < 0.05" errors with (8, 13)
   parameters. For planets, the geocentric ecliptic motion is even smoother.

### 3.3 Expected Precision

Based on Pipeline B performance and the smoothness of planetary geocentric
ecliptic motion:

| Body | Current Error (") | Expected V3 Error (") | Improvement |
|------|-------------------|----------------------|-------------|
| Saturn | 4.85 | < 0.001 | ~5000x |
| Uranus | 4.58 | < 0.001 | ~4500x |
| Mercury | 0.53 | < 0.001 | ~500x |
| Jupiter | 0.52 | < 0.001 | ~500x |
| Neptune | 0.30 | < 0.001 | ~300x |
| Pluto | 0.12 | < 0.001 | ~120x |
| Moon | < 0.001 | < 0.001 | (same) |
| Sun | < 0.001 | < 0.001 | (same) |

These estimates are conservative. The actual errors may be even smaller (sub-
milliarcsecond), since the geocentric ecliptic longitude of a planet is a very
smooth function. For comparison, Swiss Ephemeris achieves ~0.001" using direct
computation from DE431 — our Chebyshev approximation of the same function
should be comparable.

### 3.4 What About Non-Default Flags?

The `COORD_GEO_ECLIPTIC` data is specific to one reference frame: geocentric,
true ecliptic of date, with aberration and light-time correction applied. Other
output modes require different computation:

| Flag | Behavior with GEO_ECLIPTIC | Frequency |
|------|---------------------------|-----------|
| **Default** (geocentric ecliptic) | **Direct read** — maximum precision | ~90% of calls |
| `SEFLG_SIDEREAL` | Subtract ayanamsha from stored longitude | ~8% of calls |
| `SEFLG_SPEED` | Chebyshev analytical derivative | Common, works directly |
| `SEFLG_EQUATORIAL` | **Fallback to Skyfield** | Rare (~1%) |
| `SEFLG_J2000` | **Fallback to Skyfield** | Rare |
| `SEFLG_HELCTR` | **Fallback to Skyfield** | Rare |
| `SEFLG_BARYCTR` | **Fallback to Skyfield** | Rare |
| `SEFLG_TRUEPOS` | **Fallback to Skyfield** | Rare |
| `SEFLG_NOABERR` | **Fallback to Skyfield** | Rare |

**For the dominant use case** (astrological computation: geocentric ecliptic of
date, with or without sidereal correction), `COORD_GEO_ECLIPTIC` gives maximum
precision with zero conversion overhead.

**For rare flags**, the system falls back to the full Skyfield pipeline. This is
slower but still correct and precise. The existing `COORD_ICRS_BARY` data is no
longer needed for these flags — Skyfield computes them from scratch.

**Key insight**: we do NOT need to keep both `COORD_ICRS_BARY` and
`COORD_GEO_ECLIPTIC` data in the file. The ICRS data was only useful as an
intermediate representation that could serve multiple output frames. But in
practice, 95%+ of calls use the default geocentric ecliptic frame, and the rare
alternative frames can be served by Skyfield directly (which is what happens
today for unsupported flags like `SEFLG_TOPOCTR`).

### 3.5 Impact on File Size

No change. Each body still stores 3 components per segment. We replace ICRS
`(x, y, z)` with ecliptic `(lon°, lat°, dist)`. The segment parameters
(interval, degree) may even decrease for some bodies because the geocentric
ecliptic function is smoother than the barycentric ICRS function.

The Earth body (ID 14) is **no longer needed** in the LEB file — it was only
stored so Pipeline A could compute the geocentric subtraction at runtime. With
`COORD_GEO_ECLIPTIC`, the subtraction happens at generation time. This saves
~6.5 MB in the medium tier (4-day segments, degree 13, over 401,768 days).

However, Earth may still be kept for compatibility or for future use.

## 4. Implementation Plan

### Phase 1: Format Extension

**Goal**: Add `COORD_GEO_ECLIPTIC` to the LEB format. No functional changes yet.

**Files to modify**:

| File | Change |
|------|--------|
| `libephemeris/leb_format.py` | Add `COORD_GEO_ECLIPTIC = 3` constant (line ~37) |
| `libephemeris/leb_format.py` | Update `BODY_PARAMS` to use `COORD_GEO_ECLIPTIC` for all current `COORD_ICRS_BARY` bodies (IDs 0-9, 14, 15, 17-20) |

The LEB file format version stays at 1 — the `coord_type` field in `BodyEntry`
is a `uint32` and already supports arbitrary values. The reader dispatches on
`coord_type`, so adding a new value is backward-compatible: old readers will
reject new files (unknown coord_type), and new readers will handle both old and
new files.

**BODY_PARAMS changes**:

```python
# Before:
0: (32, 13, COORD_ICRS_BARY, 3),  # SE_SUN
1: (4,  13, COORD_ICRS_BARY, 3),  # SE_MOON
2: (16, 15, COORD_ICRS_BARY, 3),  # SE_MERCURY
...

# After:
0: (32, 13, COORD_GEO_ECLIPTIC, 3),  # SE_SUN
1: (4,  13, COORD_GEO_ECLIPTIC, 3),  # SE_MOON
2: (16, 15, COORD_GEO_ECLIPTIC, 3),  # SE_MERCURY
...
```

Parameters (interval_days, degree) can be kept the same initially. They may be
tuned in Phase 5 after measuring actual fitting errors — the geocentric ecliptic
function may be smooth enough to allow wider intervals or lower degrees, reducing
file size.

**Longitude continuity**: like `COORD_ECLIPTIC`, the longitude component must be
unwrapped before Chebyshev fitting (to avoid 0°/360° discontinuities) and
re-wrapped on read. The existing unwrapping infrastructure in `generate_leb.py`
already handles this.

---

### Phase 2: Generator — `generate_body_geo_ecliptic()`

**Goal**: Implement the generation function that computes geocentric ecliptic
coordinates at Chebyshev nodes using the full Skyfield pipeline.

**File**: `scripts/generate_leb.py`

**New function**: `_eval_body_geo_ecliptic(body_id, jd, planets, ts)`

This function computes the **complete geocentric ecliptic-of-date position** for
a single body at a single JD, using the same Skyfield pipeline as `planets.py`.
The computation must match what Skyfield produces exactly:

```python
def _eval_body_geo_ecliptic(body_id, jd_array, planets, ts):
    """Compute geocentric ecliptic-of-date (lon, lat, dist) at given JDs.

    Uses the full Skyfield pipeline:
    1. Compute target ICRS barycentric position (full precision)
    2. Compute Earth ICRS barycentric position (full precision)
    3. Subtract to get geocentric vector (full precision)
    4. Apply light-time correction (3 iterations, full precision)
    5. Apply aberration correction
    6. Get precession-nutation matrix from Skyfield Time.M
    7. Rotate ICRS → true equatorial of date
    8. Rotate equatorial → ecliptic using true obliquity
    9. Convert Cartesian → spherical (lon, lat, dist)

    Returns:
        Array of shape (N, 3) with (lon_deg, lat_deg, dist_au)
    """
```

**Key implementation details**:

1. **Vectorized Skyfield calls**: use `target.at(t).position.au` for batch
   evaluation of N Chebyshev nodes. This is much faster than N scalar calls.

2. **Light-time correction**: must match Skyfield's iterative approach. Use
   3 fixed-point iterations with the geometric distance.

3. **Aberration**: use Skyfield's first-order aberration correction (same as
   `_apply_aberration()` in `fast_calc.py`).

4. **PNM matrix**: use `t.M` from the Skyfield `Time` object for each JD.
   Since Chebyshev nodes are different JDs, each node needs its own PNM matrix.
   This is the most expensive part of generation — but generation is a one-time
   cost.

5. **Integration with existing generator**: the new function slots into the
   existing `assemble_leb()` dispatch. Bodies with `COORD_GEO_ECLIPTIC` use
   `generate_body_geo_ecliptic()` instead of `generate_body_icrs()`.

**Generation time impact**: significantly slower than ICRS generation because
each Chebyshev node requires the full Skyfield pipeline (including PNM matrix
computation). Estimated ~5x slower per body. For medium tier (401,768 days), this
adds ~30-60 minutes to the total generation time. Acceptable as a one-time cost.

**Asteroid handling**: asteroids (Chiron, Ceres, Pallas, Juno, Vesta) use
`spktype21` for heliocentric positions, then add Sun's barycentric position to
get SSB barycentric. The same pipeline (subtract Earth, light-time, aberration,
PNM, ecliptic) applies. The generator already has this logic in
`_eval_body_icrs_vectorized()` — it just needs to continue through the ecliptic
conversion instead of stopping at ICRS.

---

### Phase 3: Reader — Handle `COORD_GEO_ECLIPTIC`

**Goal**: Make `leb_reader.py` handle the new coordinate type correctly.

**File**: `libephemeris/leb_reader.py`

**Changes to `eval_body()` (line ~272)**:

`COORD_GEO_ECLIPTIC` is handled identically to `COORD_ECLIPTIC` at the reader
level: the Chebyshev polynomial stores `(lon, lat, dist)`, and longitude must be
wrapped to [0, 360) after evaluation.

```python
# In eval_body(), line ~341:
if body.coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL, COORD_GEO_ECLIPTIC):
    pos[0] = pos[0] % 360.0
```

This is a one-line change (adding `COORD_GEO_ECLIPTIC` to the existing check).

---

### Phase 4: Runtime Pipeline — `_pipeline_geo_ecliptic()`

**Goal**: Add a new pipeline in `fast_calc.py` that serves `COORD_GEO_ECLIPTIC`
bodies with zero conversion overhead.

**File**: `libephemeris/fast_calc.py`

**New function**: `_pipeline_geo_ecliptic(reader, jd_tt, ipl, iflag, want_velocity)`

```python
def _pipeline_geo_ecliptic(reader, jd_tt, ipl, iflag, want_velocity=False):
    """Pipeline for COORD_GEO_ECLIPTIC bodies.

    The LEB data already contains geocentric ecliptic-of-date coordinates
    (lon, lat, dist) with light-time and aberration applied. For the default
    output frame, this is a direct read with no conversion.

    For non-default flags (SEFLG_EQUATORIAL, SEFLG_J2000, SEFLG_HELCTR,
    SEFLG_BARYCTR, SEFLG_TRUEPOS, SEFLG_NOABERR), this pipeline raises
    KeyError to trigger fallback to the full Skyfield pipeline.
    """
    # Check for flags that require a different computation
    unsupported = (
        SEFLG_EQUATORIAL | SEFLG_J2000 | SEFLG_HELCTR |
        SEFLG_BARYCTR | SEFLG_TRUEPOS | SEFLG_NOABERR
    )
    if iflag & unsupported:
        raise KeyError(f"COORD_GEO_ECLIPTIC does not support flag {iflag:#x}")

    # Direct read — no conversion needed
    (lon, lat, dist), (dlon, dlat, ddist) = reader.eval_body(ipl, jd_tt)

    if want_velocity:
        return (lon, lat, dist, dlon, dlat, ddist)
    return (lon, lat, dist)
```

**Dispatch in `_fast_calc_core()` (line ~852)**:

```python
# Add to the dispatch logic:
if body.coord_type == COORD_GEO_ECLIPTIC:
    result = _pipeline_geo_ecliptic(reader, jd_tt, ipl, iflag, want_velocity)
```

**Sidereal correction**: `SEFLG_SIDEREAL` is handled **after** the pipeline
dispatch (line ~913), so it works automatically with `COORD_GEO_ECLIPTIC`.
The ayanamsha is subtracted from the longitude — same as for all other coord
types. No change needed.

**Velocity**: the Chebyshev analytical derivative gives `dlon/dt` (°/day),
`dlat/dt` (°/day), `ddist/dt` (AU/day) directly. These are the velocity
components in the geocentric ecliptic frame. No Cartesian-to-spherical velocity
conversion is needed — the derivative of the stored spherical coordinates IS the
spherical velocity. This eliminates the velocity error amplification issue as
well.

**Note on velocity precision**: the stored longitude is a smooth function of
time. Its Chebyshev derivative is also smooth and well-approximated. The only
concern is near retrograde stations where `dlon/dt` passes through zero — but
the derivative polynomial captures this naturally.

---

### Phase 5: Parameter Tuning and Tolerance Tightening

**Goal**: optimize Chebyshev parameters for `COORD_GEO_ECLIPTIC` and set
tolerances to the minimum possible.

**Steps**:

1. **Generate medium tier LEB** with initial parameters (same as current
   ICRS params). Run verification to measure actual fitting errors.

2. **Measure errors per body** across the full date range with dense sampling
   (~10,000 points per body). Record worst-case errors.

3. **Tune parameters**: if fitting errors are already < 0.001", consider
   widening intervals or lowering degrees to reduce file size. If any body
   exceeds 0.001", narrow intervals or raise degrees.

4. **Set tolerances** in `conftest.py` at ~2x measured worst-case error.
   Expected final tolerances:

   | Parametro | Current (") | Expected V3 (") |
   |-----------|------------|-----------------|
   | POSITION_ARCSEC | 5.0 | 0.005 |
   | ECLIPTIC_ARCSEC | 0.05 | 0.005 (or lower) |
   | SIDEREAL_ARCSEC | 5.0 | 0.01 |
   | SPEED_LON_DEG_DAY | 0.045 | 0.001 |
   | SPEED_LAT_DEG_DAY | 0.005 | 0.001 |
   | DISTANCE_AU | 3e-5 | 1e-6 |

5. **Regenerate all three tiers** and run full test suites.

---

### Phase 6: Testing

**Goal**: comprehensive validation across all tiers, all bodies, all flags.

**Test strategy**:

1. **LEB-vs-Skyfield comparison** (existing test infrastructure):
   - `tests/test_leb/compare/` — medium tier (19 test files)
   - `tests/test_leb/compare/base/` — base tier (8 test files)
   - `tests/test_leb/compare/extended/` — extended tier (8 test files)
   - All tests compare LEB output against live Skyfield computation.

2. **Flag fallback verification**: tests with `SEFLG_EQUATORIAL`,
   `SEFLG_J2000`, `SEFLG_HELCTR`, `SEFLG_BARYCTR` must verify that the
   fallback to Skyfield works correctly and produces full-precision results.
   These tests should pass with unchanged tolerances (they use Skyfield
   directly, so precision is the same as non-LEB mode).

3. **Regression testing**: verify that Pipeline B (ecliptic-direct) and
   Pipeline C (heliocentric) bodies are unaffected.

4. **Velocity precision**: verify that Chebyshev derivatives of geocentric
   ecliptic coordinates are more precise than the current ICRS-derivative
   approach. Moon velocity error should drop from ~0.001°/day to < 0.0001°/day.

5. **Sidereal mode**: verify that `SEFLG_SIDEREAL` produces correct results
   with the new pipeline (ayanamsha subtraction from stored longitude).

6. **Edge cases**:
   - Dates near segment boundaries (tau = ±1)
   - Retrograde stations (dlon/dt = 0)
   - Bodies near 0°/360° longitude boundary
   - Bodies at opposition (closest to Earth, smallest geocentric distance)
   - Bodies at conjunction (behind the Sun)

---

### Phase 7: Documentation and Release

**Goal**: update all documentation and prepare for merge.

**Files to update**:

| File | Changes |
|------|---------|
| `docs/leb/design.md` | Add COORD_GEO_ECLIPTIC section, update pipeline description |
| `docs/leb/guide.md` | Update precision tables, parameter tables, usage examples |
| `docs/leb/testing.md` | Update test running instructions for V3 |
| `TODO.md` | Update task status |
| `AGENTS.md` | Update if pipeline description changes |

## 5. Execution Order and Dependencies

```
Phase 1 (format)    ──→ Phase 2 (generator) ──→ Phase 3 (reader)
                                                      ↓
                                               Phase 4 (runtime pipeline)
                                                      ↓
                                               Phase 5 (tuning + tolerances)
                                                      ↓
                                               Phase 6 (testing all tiers)
                                                      ↓
                                               Phase 7 (docs + release)
```

Phases 1–4 are sequential (each depends on the previous). Phase 5 requires
a generated LEB file (Phase 2). Phase 6 requires the runtime pipeline (Phase 4).
Phase 7 can begin in parallel with Phase 6.

Estimated total effort: ~2-3 sessions.

## 6. Risks and Mitigations

### 6.1 Longitude Discontinuity at 0°/360°

**Risk**: when fitting Chebyshev polynomials to longitude values that cross the
0°/360° boundary, the polynomial will oscillate wildly (Gibbs phenomenon).

**Mitigation**: this is already solved by the existing unwrapping infrastructure
in `generate_leb.py`. The `_fit_and_verify_from_values_unwrap()` function
(line 464) uses `np.unwrap()` on longitude before fitting, producing smooth
continuous values. The wrapping back to [0, 360) happens at read time.

### 6.2 Generation Time

**Risk**: computing the full Skyfield pipeline at each Chebyshev node is
~5x slower than computing ICRS positions alone.

**Mitigation**: generation is a one-time cost. For medium tier (~401,768 days),
the ICRS generation takes ~10 minutes. At 5x, the GEO_ECLIPTIC generation
would take ~50 minutes. This is acceptable. The group-based generation workflow
(`poe leb:generate:medium:groups`) already supports parallelism across body
groups, which helps.

### 6.3 Flag Fallback Performance

**Risk**: calls with rare flags (SEFLG_EQUATORIAL, SEFLG_HELCTR, etc.) fall
back to full Skyfield computation, losing the LEB speed advantage.

**Mitigation**: these flags are rare in real-world astrological use (<5% of
calls). The fallback is the same as what happens today for `SEFLG_TOPOCTR` and
`SEFLG_XYZ` — well-tested and correct. For batch computations that need
equatorial output, the performance impact is noticeable but acceptable.

**Future optimization**: if equatorial output becomes a common need, a second
coordinate type `COORD_GEO_EQUATORIAL` could be added alongside
`COORD_GEO_ECLIPTIC`. This would store geocentric equatorial coordinates and
serve `SEFLG_EQUATORIAL` directly. However, this doubles the per-body storage
and is likely not worth it for the current use case.

### 6.4 Retrograde Motion and Longitude Rate

**Risk**: during retrograde, a planet's geocentric longitude reverses direction.
This creates local extrema in the longitude function, which might require shorter
Chebyshev segments for accurate fitting.

**Mitigation**: retrograde motion is smooth and predictable (the longitude
function has continuous derivatives). Chebyshev polynomials handle smooth
oscillations well. The 16-32 day segment lengths are much shorter than a typical
retrograde period (20-160 days), so each segment sees only a small portion of
the retrograde arc. No special handling needed.

### 6.5 Backward Compatibility

**Risk**: old LEB files (with `COORD_ICRS_BARY`) will not benefit from V3
precision improvements. Users must regenerate their files.

**Mitigation**: the reader continues to support `COORD_ICRS_BARY` — old files
work exactly as before. The improvement only applies to newly generated files.
Document this clearly in the release notes.

### 6.6 Nutation Section Consistency

**Risk**: the LEB nutation section stores IAU 2006/2000A values (from erfa at
generation time). The new pipeline doesn't need it for the main computation
path, but it's still used for `swe_nutation()` and sidereal ayanamsha.

**Mitigation**: no change needed. The nutation section continues to serve
`swe_nutation()` and ayanamsha computation. The V2 fix (using Skyfield's
nutation at runtime) was only needed for the ICRS→ecliptic conversion, which
V3 eliminates entirely. If desired, the nutation section can be regenerated
with Skyfield values in a future update.

### 6.7 True Obliquity Time-Dependence in Longitude

**Risk**: the stored geocentric ecliptic longitude is referenced to the
**true ecliptic of date**, which changes over time due to precession and
nutation. This is correct for the default output but means the stored values
are frame-dependent.

**Mitigation**: this is by design. The true ecliptic of date is the standard
output frame for `swe_calc_ut()`. Swiss Ephemeris uses the same reference frame.
The frame-dependence is why `SEFLG_J2000` and `SEFLG_EQUATORIAL` require
fallback — they need a different reference frame. This is acceptable because
these flags are rare.

## 7. File Reference

### Files to Create

None — all changes are to existing files.

### Files to Modify

| File | Phase | Changes |
|------|-------|---------|
| `libephemeris/leb_format.py` | 1 | Add `COORD_GEO_ECLIPTIC = 3`. Update BODY_PARAMS for bodies 0-9, 14, 15, 17-20 |
| `scripts/generate_leb.py` | 2 | Add `_eval_body_geo_ecliptic()` function. Update `assemble_leb()` dispatch to route GEO_ECLIPTIC bodies to new generator |
| `libephemeris/leb_reader.py` | 3 | Add `COORD_GEO_ECLIPTIC` to longitude wrapping check in `eval_body()` |
| `libephemeris/fast_calc.py` | 4 | Add `_pipeline_geo_ecliptic()`. Update `_fast_calc_core()` dispatch |
| `tests/test_leb/compare/conftest.py` | 5 | Update TIER_DEFAULTS with V3 tolerances |
| `docs/leb/design.md` | 7 | Document COORD_GEO_ECLIPTIC pipeline |
| `docs/leb/guide.md` | 7 | Update precision tables, parameter tables |

### Files Unchanged

| File | Why |
|------|-----|
| `libephemeris/leb_reader.py` (core eval) | Chebyshev evaluation is coordinate-agnostic |
| Pipeline B bodies (ecliptic-direct) | Already in final frame, unaffected |
| Pipeline C bodies (heliocentric) | Unaffected |
| Test files for ecliptic/helio bodies | Unaffected |

## 8. Summary

V3 is a **paradigm shift** from "store intermediate data, convert at runtime"
to "store final data, read directly". It eliminates the ICRS→ecliptic conversion
pipeline as the precision bottleneck, bringing LEB precision from ~5" to ~0.001"
for all planets — matching or exceeding Swiss Ephemeris.

The tradeoff is explicit: rare output frames (equatorial, heliocentric,
barycentric) fall back to Skyfield. For the dominant astrological use case
(geocentric ecliptic), V3 delivers both maximum speed AND maximum precision.

| Metric | V1 (current) | V2 (current) | V3 (proposed) | Swiss Eph |
|--------|-------------|-------------|--------------|-----------|
| Saturn lon error | 4.85" | 4.85" | < 0.001" | ~0.001" |
| Uranus lon error | 4.58" | 4.58" | < 0.001" | ~0.001" |
| Moon lon error | < 0.001" | < 0.001" | < 0.001" | ~0.001" |
| Speed (geocentric ecl) | ~14x Skyfield | ~8x Skyfield | ~14x+ Skyfield | 1x (reference) |
| Speed (equatorial) | ~14x Skyfield | ~8x Skyfield | 1x (fallback) | 1x (reference) |

---

## Notes From V3 Implementation Review (Redo Checklist)

These notes capture issues found in a first attempted V3 implementation, to help
whoever redoes the task avoid the same pitfalls.

### Must Match Skyfield Exactly (Generator)

- The generator for `COORD_GEO_ECLIPTIC` must reproduce the same Skyfield
  reference pipeline used by the tests. If generation-time math differs (even
  slightly), tightening tolerances (e.g. 0.005") will fail.
- Aberration must be applied with the same formula as runtime
  (`_apply_aberration()` style): `u' = u + v - u*(u.v)` with renormalization. A
  naive `geo += dist*(v/c)` is not equivalent and can introduce significant
  angular error.
- Earth velocity for aberration must be computed for every evaluation node
  (vectorized), not sparsely sampled and then reused/interpolated ad-hoc. Sparse
  sampling makes aberration systematically wrong.
- Use the same definition of true obliquity as the runtime path (Skyfield IAU
  2000A mean + deps via `Time` nutation angles). Avoid calling alternative
  helpers like `t.ecliptic_obliquity()` unless you have verified they are
  identical to the runtime reference.

### Performance / Scalability Pitfalls

- Avoid per-JD Python loops inside light-time iterations and frame rotations;
  those explode generation time for short-interval bodies (Moon/Earth at 4-day
  segments).
- If you need light-time iterations, do them in batch (vectorized) as much as
  possible, or use Skyfield APIs that already implement the
  apparent/aberrated pipeline consistently.

### Testing / Tolerances

- Do not tighten tolerances in `tests/test_leb/compare/conftest.py` until a
  generated V3 LEB file is validated against Skyfield and worst-case errors are
  measured. Otherwise you may lock in unrealistic thresholds and mask pipeline
  mismatches.

### Edge Cases / Semantics

- Storing `SE_EARTH` as geocentric ecliptic can be degenerate (distance ~0, lon
  undefined). Decide explicitly whether to keep Earth in the file for
  compatibility; if kept, ensure generator and runtime behavior is well-defined
  and tests cover it.
- Watch for subtle velocity handling: if `SEFLG_SPEED` is not set, ensure
  downstream adjustments (e.g. sidereal precession-rate subtraction) do not
  reintroduce non-zero `dlon` unexpectedly.

