# LEB Precision V3 — Historical: COORD_GEO_ECLIPTIC (Abandoned)

> **Status:** ABANDONED. This approach was fully implemented, tested, and
> found to fail catastrophically at retrograde stations (2–135" Chebyshev
> fitting errors). The code was removed and replaced by
> `COORD_ICRS_BARY_SYSTEM` (type 4) with runtime PPN gravitational deflection
> and relativistic aberration. See `algorithms.md` §16 (Problem 1) for the
> canonical summary.
>
> This document is preserved as a historical record and educational reference
> for why storing geocentric ecliptic coordinates in Chebyshev polynomials
> does not work.

---

## Timeline

| Date | Event | Commit |
|------|-------|--------|
| Design | This document written as a proposal | — |
| Implementation | Full implementation across generator, reader, format, pipeline | `6420874` |
| Tuning | Parameter adjustments, tolerance work | `2113075`, `84954bb`, `07941b5`, `f16ce62` |
| **Failure** | Retrograde cusps cause 2–135" fitting errors | `59c03af` |
| Replacement | `COORD_ICRS_BARY_SYSTEM` (type 4) + PPN deflection + full aberration | `01d33f5` |
| Cleanup | ~700 lines of dead COORD_GEO_ECLIPTIC code removed | `89541cb` |
| Release | Shipped as v0.23.0 | `8fb5220` |

The constant `COORD_GEO_ECLIPTIC = 3` remains in `leb_format.py` as
"reserved, not used" for format compatibility.

---

## 1. The Problem That Motivated This Approach

After V1 (Chebyshev parameter tuning) and V2 (nutation model fix), LEB had
a persistent structural precision problem: **ICRS-to-ecliptic conversion
amplified Chebyshev fitting errors** for outer planets.

Measured worst-case errors:

| Body | Error (") | Root Cause |
|------|-----------|------------|
| Saturn | 4.85 | ICRS→ecliptic amplification |
| Uranus | 4.58 | ICRS→ecliptic amplification |
| Mercury | 0.53 | ICRS→ecliptic amplification |
| Jupiter | 0.52 | ICRS→ecliptic amplification |

The errors came from the runtime conversion pipeline (geocentric subtraction,
light-time iteration, coordinate rotation, Cartesian-to-spherical), not from
the Chebyshev fitting itself. The fitting in ICRS Cartesian was excellent
(sub-arcsecond), but the conversion chain amplified errors by a factor of
`1/geocentric_distance × 206265`.

## 2. The Proposed Solution

Store **geocentric ecliptic-of-date coordinates** (lon°, lat°, dist AU)
directly in the LEB file as `COORD_GEO_ECLIPTIC = 3`. The full Skyfield
pipeline would run at generation time with float64 precision, and the
Chebyshev polynomial would approximate the final ecliptic output directly.
At runtime, no conversion — just polynomial evaluation.

```
COORD_ICRS_BARY (existing):
  Generation: Skyfield ICRS (x,y,z) → fit Chebyshev
  Runtime:    Chebyshev → (x,y,z)±ε → subtract Earth → light-time → PNM → ecl → spherical
  Error:      ε amplified by conversion chain → ~5"

COORD_GEO_ECLIPTIC (proposed):
  Generation: Skyfield → full pipeline → ecliptic (lon, lat, dist) → fit Chebyshev
  Runtime:    Chebyshev → (lon, lat, dist) directly
  Error:      ε from fitting only → predicted <0.001"
```

The rationale was that Pipeline B bodies (TrueNode, OscuApogee, etc.) already
stored ecliptic coordinates directly and achieved <0.05" errors, proving that
Chebyshev approximation of ecliptic coordinates could work. Geocentric
planetary longitudes were expected to be similarly smooth.

## 3. Why It Failed: Retrograde Cusps

### The Fatal Flaw

The assumption that "geocentric ecliptic longitude is a smooth function of
time" was **wrong**. Planets undergo retrograde motion, and at retrograde
stations the ecliptic longitude has a **cusp** — a sharp reversal where the
derivative changes sign abruptly.

While the longitude function is technically C∞ (infinitely differentiable),
the curvature at stations is extreme. Chebyshev polynomials of practical
degree (13–15) cannot follow these sharp reversals without large oscillation
errors (Gibbs-like ringing).

### Measured Failure

When tested across the full date range:

- **Fitting errors at retrograde stations**: 2–135 arcseconds
- **Worst bodies**: outer planets (Saturn, Jupiter) where retrograde arcs
  are most pronounced
- **Pattern**: errors concentrated at the ~2-week window around each
  retrograde station, with near-zero errors elsewhere

The 135" worst case is **135,000x worse** than the 0.001" target. No
parameter tuning could fix this — shorter segments would need to be ~0.5
days to capture the cusp shape, which would increase file size by 60x
and defeat the purpose of precomputation.

### Why Pipeline B Bodies Don't Have This Problem

Pipeline B bodies (TrueNode, MeanNode, OscuApogee, MeanApogee, IntpApogee,
IntpPerigee) are **lunar orbit elements**, not geocentric positions of
distant bodies. They don't undergo retrograde motion — their ecliptic
coordinates are smooth, monotonic, or slowly oscillating functions. The
success of Pipeline B was not evidence that geocentric planetary longitudes
would also work.

### Why ICRS Barycentric Doesn't Have This Problem

In ICRS barycentric coordinates, there is no retrograde motion. A planet's
barycentric position traces a smooth, nearly elliptical orbit. Retrograde
is an artifact of the geocentric perspective (Earth overtaking an outer
planet, or an inner planet lapping Earth). By storing the raw barycentric
position and computing the geocentric view at runtime, the Chebyshev
polynomial only needs to approximate the smooth orbital motion.

## 4. What Actually Shipped

The precision problem was ultimately solved by a different approach:

### COORD_ICRS_BARY_SYSTEM (type 4)

For outer planets (Jupiter–Pluto), the original `COORD_ICRS_BARY` stored
planet-center positions, which included high-frequency center-of-body (COB)
oscillations from inner moons. These oscillations were the main source of
Chebyshev fitting difficulty.

`COORD_ICRS_BARY_SYSTEM` stores the **system barycenter** (smooth) and
applies the COB correction at runtime using `planet_centers.bsp` or
analytical moon theory. This eliminated the high-frequency component from
the stored data, achieving <0.001" fitting errors for all outer planets.

### PPN Gravitational Deflection

Added `_apply_gravitational_deflection()` to the runtime pipeline, modeling
light bending by the Sun, Jupiter, and Saturn using the PPN formula. This
fixed the 3.95" systematic error for Saturn and ~0.002" errors for all
planets.

### Full Relativistic Aberration

The runtime pipeline now applies proper stellar aberration consistent with
Skyfield's `apparent()` method.

### Result

All 31 bodies achieve <0.001" precision across all three tiers, with the
smooth-in-ICRS approach preserved. No ecliptic storage needed.

## 5. Lessons Learned

1. **Smoothness in one frame does not imply smoothness in another.** ICRS
   barycentric motion is inherently smooth; geocentric ecliptic motion has
   frame-dependent singularities (retrograde cusps).

2. **Test with adversarial cases early.** The initial testing used random
   date samples that rarely hit retrograde stations. Dense sampling near
   stations immediately revealed the failure.

3. **Pipeline B success was misleading.** It proved Chebyshev works for
   ecliptic coordinates of *specific bodies* (lunar nodes/apsides), not for
   ecliptic coordinates in general.

4. **The real bottleneck was elsewhere.** The ~5" errors weren't from the
   ICRS→ecliptic conversion amplification alone — they were primarily from
   COB oscillations in the stored ICRS data. Fixing the storage (system
   barycenters) fixed the actual problem without changing the coordinate
   frame.

---

## Appendix: Original Root Cause Analysis

The following sections from the original proposal document the error
amplification mechanism in Pipeline A. This analysis was correct but
led to the wrong solution (changing the storage frame instead of fixing
what was stored in the existing frame).

### The Pipeline A Conversion Chain

When a user calls `swe_calc_ut()` for a planet (e.g., Saturn), Pipeline A
executes this chain at runtime:

```
Step 1: target_bary = leb.eval_body(SATURN, jd)     → Chebyshev eval → (x,y,z) + ε_target
Step 2: earth_bary  = leb.eval_body(EARTH, jd)      → Chebyshev eval → (x,y,z) + ε_earth
Step 3: geo = target_bary - earth_bary               → errors ADD: ε_geo = ε_target + ε_earth
Step 4: light-time (3 iterations of steps 1-3)        → 3 more ε_target errors
Step 5: aberration correction                         → small additional error
Step 6: PNM matrix rotation                           → error preserved (rigid rotation)
Step 7: ecliptic rotation                             → error preserved (rigid rotation)
Step 8: atan2(y,x), asin(z/r)                        → angular_error ≈ ε_cartesian / dist_geo
Step 9: → (lon°, lat°, dist AU)                      → final output with amplified error
```

### The Error Amplification Mechanism

The critical step is Step 8: converting Cartesian to spherical coordinates.
The angular error from a Cartesian position error ε (in AU) is:

```
angular_error (radians) ≈ ε / geocentric_distance
angular_error (arcsec)  ≈ (ε / geocentric_distance) × 206265
```

For Saturn (geocentric distance ~8–11 AU): a 2×10⁻⁴ AU Cartesian error
gives 2×10⁻⁴/9 × 206265 ≈ 4.6". For Moon (geocentric distance ~0.0025 AU):
fitting is extremely tight (4-day segments), so ε is ~10⁻¹⁰ AU, giving
negligible angular error.

### Why Pipeline B Bodies Had No Problem

Bodies stored as `COORD_ECLIPTIC` (Pipeline B) had errors of <0.05" because
they stored final ecliptic coordinates directly — no subtraction, no
coordinate transformation, no error amplification. Similarly,
`COORD_HELIO_ECL` bodies (Pipeline C) had <0.001" errors.

This correctly identified that the ICRS→ecliptic conversion was the error
source, but the solution (store ecliptic directly) failed due to retrograde
cusps — a problem that Pipeline B bodies don't have.
