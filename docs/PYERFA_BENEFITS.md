# PyERFA Integration Benefits

## Overview

LibEphemeris optionally integrates with [PyERFA](https://github.com/liberfa/pyerfa), the Python wrapper for the IAU SOFA/ERFA library. When installed, PyERFA provides access to the most accurate nutation and precession models available.

## Installation

PyERFA is included as a dependency of libephemeris and is installed automatically:

```bash
pip install libephemeris
```

Or explicitly:

```bash
pip install pyerfa
```

## Precision Comparison

LibEphemeris uses three tiers of nutation models:

| Model | Source | Precision | Terms |
|-------|--------|-----------|-------|
| IAU 2000B | Skyfield (default fallback) | ~1 milliarcsecond (mas) | 77 |
| IAU 2000A | Skyfield/PyERFA | ~0.1 mas | 1365 |
| IAU 2006/2000A (nut06a) | PyERFA only | ~0.01-0.05 mas | 1365 + corrections |

### IAU 2000B vs IAU 2000A

IAU 2000B is a truncated version of IAU 2000A with only 77 nutation terms (vs 1365). The difference is typically ~1 mas, which translates to:
- Moon position error: ~2 seconds of time
- For most astrological applications: negligible

### IAU 2000A vs IAU 2006/2000A

IAU 2006 adds J2 secular variation and obliquity frame corrections to IAU 2000A. The improvement is ~0.01-0.05 mas for dates within a century of J2000. This difference grows slowly with distance from J2000.

## Quantified Benefits

| Time from J2000 | IAU 2000B-2000A diff | IAU 2000A-2006 diff |
|----------------|---------------------|-------------------|
| 0 years | ~0.5 mas | ~0.01 mas |
| 10 years | ~0.5 mas | ~0.02 mas |
| 50 years | ~1.0 mas | ~0.05 mas |
| 100 years | ~1.5 mas | ~0.1 mas |

### Obliquity Models

PyERFA provides IAU 2006 obliquity via `erfa.obl06()`, which differs from the Laskar 1986 polynomial by ~42 mas at J2000. This offset is a known difference between the models.

## PyERFA Functions Available

When PyERFA is installed, the following functions are available via `libephemeris.erfa_nutation`:

- `get_erfa_nutation_nut00a(jd_tt)` — IAU 2000A nutation (dpsi, deps)
- `get_erfa_nutation_nut06a(jd_tt)` — IAU 2006/2000A nutation (dpsi, deps)
- `get_erfa_obliquity_iau2006(jd_tt)` — IAU 2006 mean obliquity
- `get_erfa_pnm06a_matrix(jd_tt)` — Combined precession-nutation-bias matrix
- `get_erfa_nutation_cached(jd_tt)` — Cached nutation with automatic fallback
- `compare_nutation_models(jd_tt)` — Compare all available models
- `has_erfa()` — Check if PyERFA is available

## When to Use PyERFA

### Recommended for:
- **Research requiring highest precision** — sub-milliarcsecond accuracy
- **Long-term Ephemeris**: Calculations spanning centuries where error accumulates
- **Validation against official IERS data** — matching published pole positions
- **Precession-nutation matrix applications** — using `pnm06a` avoids cross-term accumulation errors

### Not necessary for:
- **Typical astrological work** — IAU 2000A is already sufficient
- **Short-term calculations** — differences are negligible within decades of J2000
- **Applications with > 1 arcsecond tolerance** — the nutation model difference is sub-arcsecond

## Summary

PyERFA provides measurable precision improvements through the IAU 2006/2000A nutation model (`nut06a`). The improvement is ~0.01-0.05 mas for dates near J2000, growing to ~0.1 mas at 100 years. For typical astrological work, the built-in Skyfield IAU 2000A model is already sufficient, but PyERFA is recommended for research-grade precision.
