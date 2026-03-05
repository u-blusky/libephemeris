# LEB Precision V2 — Nutation Model Alignment

## 1. Problem Statement

### 1.1 Symptom

LEB-vs-Skyfield comparison tests show large position errors that grow with
distance from J2000. The errors are **not** Chebyshev fitting errors (verify_leb
confirms sub-arcsecond fitting precision for all bodies). They are caused by a
**nutation/precession model mismatch** between the LEB pipeline and Skyfield.

Observed worst-case errors by tier:

| Tier | Range | Max T (centuries) | Max Position Error | Dominant Body |
|------|-------|-------------------|--------------------|---------------|
| Base | 1849-2150 | ~1.5 | ~0.04" | Negligible (hidden by pipeline errors) |
| Medium | 1550-2650 | ~6.5 | ~0.2" | Negligible (hidden by pipeline errors) |
| Extended | -5000 to 5000 | ~70 | **234"** (0.065 deg) | Mars truepos at ~4800 BCE |

For base and medium tiers, the mismatch is small (~0.04-0.2") and hidden
beneath the ICRS-to-ecliptic pipeline amplification errors (Saturn 4.85",
Uranus 4.58"). For extended tier, it completely dominates all other error
sources.

### 1.2 Root Cause

The LEB pipeline (`fast_calc.py`) and the Skyfield pipeline (`planets.py`) use
**different nutation and precession models** for the ICRS-to-ecliptic
coordinate transformation:

#### LEB Pipeline (fast_calc.py)

| Step | Function | Model |
|------|----------|-------|
| Precession-nutation matrix (ICRS → true equatorial) | `erfa.pnm06a()` called **live** | IAU 2006/2000A (Fukushima-Williams angles) |
| True obliquity (equatorial → ecliptic) | `reader.eval_nutation()` from LEB Chebyshev | IAU 2006/2000A (from erfa at generation time) |

#### Skyfield Pipeline (planets.py)

| Step | Function | Model |
|------|----------|-------|
| Combined matrix `M = N × P × B` (ICRS → true equatorial) | Skyfield Python implementation | IAU 2000A (classical angles, no 2006 corrections) |
| True obliquity (equatorial → ecliptic) | `t._nutation_angles_radians` from Skyfield | IAU 2000A |

### 1.3 Why They Diverge

The IAU 2006/2000A model (used by erfa) includes corrections not present in
Skyfield's IAU 2000A implementation:

1. **J2 secular variation**: A correction of ~-25.8 mas/century to the nutation
   in longitude, accounting for the secular decrease in Earth's J2 (dynamic
   flattening) due to post-glacial rebound.

2. **Precession formulation**: erfa uses the Fukushima-Williams 4-angle method;
   Skyfield uses the classical Lieske/Capitaine angles. Numerically equivalent
   near J2000 but diverge at extreme epochs.

3. **Precession rate adjustments**: Small differences in the polynomial
   coefficients for general precession.

At T = 70 centuries (4800 BCE), general precession totals ~352,000 arcseconds.
A relative difference of ~0.07% between the two models yields ~250" of
positional error — matching the observed 234" almost exactly.

### 1.4 Proof Points

| Test | Result | Explanation |
|------|--------|-------------|
| `verify_leb` (Chebyshev vs erfa source function) | < 0.03" | Same model on both sides → fitting error only |
| J2000 frame tests (`SEFLG_J2000`) | PASS at 5" | No precession/nutation rotation applied |
| Default ecliptic tests (planet longitude) | 14-81" FAIL | Nutation+precession model mismatch |
| Flag tests (equatorial/helio/bary/truepos/noaberr) | 175-234" FAIL | Full frame rotation mismatch amplified |
| Extended > Medium > Base error growth | Confirmed | More centuries from J2000 = more divergence |

### 1.5 Impact Across Tiers

The fix benefits all tiers, but the visible improvement depends on tier:

- **Base tier**: Errors will drop from ~0.04" to ~0.001" for the nutation
  component. However, the dominant errors (Saturn lat 4.85") are from the
  ICRS-to-ecliptic pipeline amplification (1/geocentric_distance) and are
  unchanged. Net improvement: **marginal**.

- **Medium tier**: Same story — dominant errors are pipeline amplification.
  The nutation component drops from ~0.2" to ~0.001". Net improvement:
  **small but measurable** in flag/equatorial tests.

- **Extended tier**: The nutation mismatch IS the dominant error. Fixing it
  drops errors from **234" to sub-arcsecond** — a ~1000x improvement.
  This makes extended tier precision comparable to base/medium.

## 2. Solution

### 2.1 Approach: Use Skyfield's Rotation Matrix at Runtime

Replace `erfa.pnm06a()` in `_precession_nutation_matrix()` with Skyfield's
own `Time.M` matrix. This ensures the LEB pipeline uses exactly the same
precession-nutation model as the Skyfield reference.

The change is localized to a single function in `fast_calc.py`.

### 2.2 Code Changes

#### fast_calc.py — `_precession_nutation_matrix()`

Current code (line ~242):

```python
def _precession_nutation_matrix(jd_tt: float):
    import erfa
    mat = erfa.pnm06a(J2000, jd_tt - J2000)  # IAU 2006/2000A
    return ((mat[0][0], ...), ...)
```

New code:

```python
def _precession_nutation_matrix(jd_tt: float):
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    return tuple(tuple(row) for row in t.M)  # IAU 2000A, matching Skyfield
```

#### fast_calc.py — True obliquity computation

Current code (line ~575):

```python
dpsi, deps = reader.eval_nutation(jd_tt)  # LEB Chebyshev (IAU 2006/2000A)
eps_mean = _mean_obliquity_iau2006(jd_tt)
eps_true_rad = math.radians(eps_mean) + deps
```

New code:

```python
ts = get_timescale()
t = ts.tt_jd(jd_tt)
dpsi, deps = t._nutation_angles_radians  # Skyfield IAU 2000A
eps_true_rad = t._mean_obliquity_radians + deps
```

Alternatively, since `t.M` is already computed for the PNM matrix step,
the `t` object can be passed down to avoid creating it twice.

### 2.3 LEB Nutation Section

The LEB file's nutation Chebyshev section (generated from `erfa.nut06a()`)
will become **unused** by the main pipeline. However, it should be kept:

1. It is still used by `swe_nutation()` (SE_ECL_NUT queries)
2. It is used for sidereal ayanamsha computation
3. Removing it would change the LEB file format

A future optimization could regenerate the nutation section using Skyfield's
`iau2000a_radians()` for full consistency, but this is optional.

### 2.4 Performance Impact

| Operation | Before (erfa) | After (Skyfield) | Change |
|-----------|---------------|-------------------|--------|
| PNM matrix | ~5 μs | ~50 μs | ~10x slower |
| Full swe_calc_ut() per body | ~70 μs | ~115 μs | ~1.6x slower |
| vs pure Skyfield | 14x faster | ~8x faster | Still very fast |

The performance cost is acceptable. The LEB speedup remains substantial
(~8x vs Skyfield), and the bottleneck in real applications is typically
evaluating many bodies across many dates, not single-call latency.

### 2.5 Delta-T Consideration

A secondary (minor) source of difference is **delta-T computation**:

- **LEB**: Linear interpolation on a sparse table sampled every 30 days
- **Skyfield**: Stephenson-Morrison-Hohenkerk 2016 model with cubic splines

At extreme dates, linear interpolation on 30-day spacing introduces ~0.01s
error in delta-T, translating to ~0.01" in Moon position. This is negligible
compared to the nutation fix, but could be improved later by either:

- Using Skyfield's delta-T at runtime instead of the LEB table
- Switching to Chebyshev fitting for delta-T in the LEB file
- Reducing the sampling interval (e.g., 10 days)

## 3. Implementation Plan

### Phase 1: Runtime Fix (no regeneration needed)

1. Modify `_precession_nutation_matrix()` to use `t.M` from Skyfield
2. Modify true obliquity computation to use `t._nutation_angles_radians`
3. Optionally: pass the Skyfield `Time` object through the pipeline to avoid
   creating it twice (performance optimization)
4. Run `poe lint` and `poe format`

### Phase 2: Test All Tiers

1. Run base tier tests: `poe test:leb:compare:base`
2. Run medium tier tests: `poe test:leb:compare:medium`
3. Run extended tier tests (targeted): `pytest tests/test_leb/compare/extended/ -v`
4. Verify no regressions on base/medium
5. Measure new error levels on all tiers

### Phase 3: Tighten Tolerances

1. Update `TIER_DEFAULTS["extended"]` with measured errors
2. Tighten `POSITION_ARCSEC` from 250" to expected ~5" or less
3. Tighten `EQUATORIAL_ARCSEC`, `SIDEREAL_ARCSEC` similarly
4. Verify 0 failures on all tiers

### Phase 4: Verify and Document

1. Run full test suite for all tiers
2. Update `TODO.md` and `NOTES.md`
3. Update `docs/leb/precision-improvement-plan.md` with V2 results
4. Commit and push

## 4. Expected Results

After the fix, expected error levels for extended tier:

| Body | Before (") | After (expected ") |
|------|------------|-------------------|
| Sun | 14.3 | < 0.5 |
| Moon | 28.4 | < 0.5 |
| Mercury | 19.1 | < 0.5 |
| Mars | 29.0 | < 0.5 |
| Pluto | 81.1 | < 0.5 |
| Flag tests (max) | 234 | < 5 (pipeline limit) |

The remaining errors after the fix will be:
- Chebyshev fitting error (< 0.03", from verify_leb)
- ICRS-to-ecliptic pipeline amplification (body-dependent, same as base/medium)
- Light-time iteration convergence (~0.001")

## 5. Risks

1. **Skyfield API stability**: `t.M` and `t._nutation_angles_radians` are
   semi-private APIs. If Skyfield changes them, the code will break. Mitigation:
   pin Skyfield version and add import-time validation.

2. **Performance regression**: The ~1.6x slowdown per call may be noticeable in
   batch computations. Mitigation: benchmark before/after, consider caching the
   Skyfield Time object when evaluating multiple bodies at the same JD.

3. **LEB nutation section becomes inconsistent**: The nutation section stores
   IAU 2006/2000A values but the pipeline now uses IAU 2000A. The section is
   only used for SE_ECL_NUT and sidereal queries, where the difference is
   <0.001" for typical dates. Mitigation: document this, optionally regenerate
   with Skyfield model later.
