# LEB Precision V3 — Retrograde Fitting Problem

## Problem Statement

The V3 `COORD_GEO_ECLIPTIC` approach stores precomputed geocentric ecliptic-of-date coordinates (lon°, lat°, dist AU) directly in the `.leb` file. This eliminates the ICRS→ecliptic conversion error that plagued V2 (~0.5–2"), but introduces a **worse problem**: Chebyshev polynomials cannot fit geocentric ecliptic longitude near planetary retrograde stations.

## Root Cause

When a planet goes retrograde, its **geocentric ecliptic longitude** reverses direction, creating a cusp-like local extremum. Chebyshev polynomials (which are smooth global approximants) struggle to fit these sharp reversals, producing large fitting errors concentrated near retrograde stations.

In contrast, **ICRS barycentric coordinates** (the V2 approach) represent the planet's actual spatial trajectory in 3D Cartesian — which is smooth everywhere, even during apparent retrogrades. The retrograde is an artifact of the geocentric perspective, not a physical trajectory feature.

## Measured Fitting Errors (V3 with current BODY_PARAMS)

| Body | Params (interval/degree) | Worst Fitting Error | Target |
|------|--------------------------|---------------------|--------|
| Sun | 16d/15 | 0.0000" | <0.001" ✅ |
| Moon | 4d/13 | 0.0008" | <0.001" ✅ |
| Mercury | 1d/17 | **14.29"** | <0.001" ❌ |
| Venus | 4d/13 | **6.92"** | <0.001" ❌ |
| Mars | 1d/17 | **3.91"** | <0.001" ❌ |
| Jupiter | 0.5d/21 | **30.98"** | <0.001" ❌ |
| Saturn | 2d/15 | **5.68"** | <0.001" ❌ |
| Uranus | 1d/23 | **134.81"** | <0.001" ❌ |
| Neptune | 4d/17 | **9.26"** | <0.001" ❌ |
| Pluto | 8d/13 | **2.34"** | <0.001" ❌ |

Sun and Moon work perfectly because they never go retrograde as seen from Earth.

## Why the Previous Validation Missed This

The fast test script `test_chebyshev_params.py --scan 500` sampled only ~500 out of 50,000+ segments (0.5–1%), randomly. The worst-case segments (near retrograde stations) were never sampled. This produced misleadingly low error figures (<0.001" for all planets).

The full generator reveals the true worst-case when it evaluates all segments.

## What Would Fix It (Brute Force)

Reducing the interval makes the fitting work (e.g., Saturn at 0.5d/23 gives 0.0005"), but at enormous file size cost:

| Body | Params for <0.001" | Estimated Size Impact |
|------|--------------------|-----------------------|
| Saturn | 0.5d/23 | ~126 MB (vs ~18 MB at 2d/15) |
| Jupiter | 0.1d/21 | ~350 MB+ |
| Uranus | 0.1d/23 | ~400 MB+ |

This would push the base tier to **multi-GB**, which is impractical.

## V2 vs V3 Comparison

| Aspect | V2 (ICRS_BARY) | V3 (GEO_ECLIPTIC) |
|--------|----------------|---------------------|
| Chebyshev fitting | <0.001" (smooth 3D trajectory) | 2–135" (retrograde cusps) |
| Runtime conversion | 0.5–2" (ICRS→ecliptic amplification) | 0" (direct read) |
| **Total worst-case** | **0.5–2"** | **2–135"** |
| Sun/Moon | 0.5–2" | **0.0000–0.0008"** |

V3 is dramatically better for Sun and Moon, but dramatically worse for planets 2–9.

## Possible Solutions

1. **Hybrid V2+V3**: Use `COORD_GEO_ECLIPTIC` for Sun and Moon (which benefit enormously), keep `COORD_ICRS_BARY` for planets 2–9 (which are smooth in ICRS).

2. **Adaptive segmentation**: Detect retrograde stations and use finer intervals only near them. Complex to implement but could achieve <0.001" without multi-GB files.

3. **Different basis**: Instead of storing (lon, lat, dist), store geocentric ICRS Cartesian (x, y, z) — which is smooth — and convert to ecliptic at runtime. This is similar to V2 but with the geocentric subtraction already done, avoiding the barycentric→geocentric amplification that caused V2's 0.5–2" error.

4. **Residual correction**: Store V2 ICRS_BARY Chebyshev + a low-order correction polynomial for the ecliptic conversion error.

## Status

- Base tier planets generated with V3 params but fitting errors are unacceptable
- Investigation complete, solution approach TBD
- See `docs/leb/leb_precision_v3.md` for original V3 spec
