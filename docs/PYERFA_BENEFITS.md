# PyERFA Integration Benefits

## Overview

PyERFA is the Python wrapper for the SOFA/ERFA library, providing IAU-standard implementations of fundamental astronomical algorithms. Libephemeris integrates PyERFA to provide milliarcsecond (mas) precision for critical operations like precession and nutation calculations.

This document describes the precision benefits of using PyERFA integration and when to use it for astronomical calculations.

## Installation

PyERFA is an optional dependency. Install it with:

```bash
pip install pyerfa
```

Or with libephemeris extras:

```bash
pip install libephemeris[precision]
```

To verify PyERFA is available:

```python
from libephemeris.erfa_nutation import has_erfa
if has_erfa():
    print("PyERFA is available for high-precision calculations")
```

## Precision Comparison

Libephemeris provides multiple nutation models with different precision characteristics:

### IAU 2000A vs IAU 2000B vs IAU 2006

| Model | Precision | Use Case |
|-------|-----------|----------|
| IAU 2000A | ~0.2 mas | Highest precision, ~1300 terms |
| IAU 2000B | ~1 mas | High precision, ~80 terms, faster |
| IAU 2006/2000A (nut06a) | ~0.01-0.05 mas | Best precision with improved precession |

The nut06a model combines the IAU 2006 precession with IAU 2000A nutation, providing the best available precision for equinox-based coordinates.

### Error Growth Over Time

The differences between models grow with distance from J2000.0:

| Years from J2000 | IAU 2000B vs 2000A | nut06a vs nut00a |
|------------------|-------------------|------------------|
| 0 | ~0.2 mas | ~0.01 mas |
| 10 | ~0.5 mas | ~0.02 mas |
| 50 | ~2 mas | ~0.05 mas |
| 100 | ~5 mas | ~0.1 mas |

## Quantified Benefits

### Nutation Precision

- **IAU 2000A**: Full-precision nutation model (~0.2 mas accuracy)
- **IAU 2000B**: Truncated model (~1 mas accuracy), 10x faster
- **nut06a**: IAU 2006 precession + IAU 2000A nutation (~0.05 mas accuracy)

### Obliquity Precision

PyERFA provides the IAU 2006 obliquity model with milliarcsecond accuracy:

```python
from libephemeris.erfa_nutation import get_erfa_obliquity_iau2006
obliquity = get_erfa_obliquity_iau2006(jd_tt)  # Returns radians
```

### Transformation Matrix

For complete coordinate transformations:

```python
from libephemeris.erfa_nutation import get_erfa_pnm06a_matrix
# Returns the precession-nutation matrix for J2000 to date transformation
matrix = get_erfa_pnm06a_matrix(jd_tt)
```

## PyERFA Functions Available

Libephemeris exposes the following PyERFA-based functions:

| Function | Description |
|----------|-------------|
| `has_erfa()` | Check if PyERFA is available |
| `get_erfa_nutation_nut00a(jd)` | IAU 2000A nutation (dpsi, deps in radians) |
| `get_erfa_nutation_nut00b(jd)` | IAU 2000B nutation (truncated, faster) |
| `get_erfa_nutation_nut06a(jd)` | IAU 2006/2000A combined nutation |
| `get_erfa_obliquity_iau2006(jd)` | Mean obliquity using IAU 2006 model |
| `get_erfa_pnm06a_matrix(jd)` | Precession-nutation matrix |
| `compare_nutation_models(jd)` | Compare all available models |
| `get_erfa_nutation_cached(jd)` | Cached nutation calculation |

### Cached Nutation

For performance-critical applications:

```python
from libephemeris.erfa_nutation import get_erfa_nutation_cached
dpsi, deps = get_erfa_nutation_cached(jd_tt)
```

This uses an LRU cache to avoid redundant calculations for the same Julian date.

## When to Use PyERFA

### Recommended Use Cases

1. **High-precision coordinate transformations**: When you need milliarcsecond accuracy for precessing coordinates between epochs

2. **Observatory-grade calculations**: Professional astronomical applications requiring IAU-standard algorithms

3. **Validation against external sources**: When comparing results with JPL Horizons or other professional ephemeris services

4. **Long-term extrapolation**: Calculations more than 50 years from J2000.0 benefit most from IAU 2006 improvements

### When the Default is Sufficient

1. **Typical astrological applications**: The default precision (~1 arcsec) is adequate for most astrological chart calculations

2. **Performance-critical code**: The built-in approximation is faster when milliarcsecond precision isn't required

3. **Dates near J2000.0**: Within 50 years of J2000.0, the differences between models are minimal for most applications

## Summary

PyERFA integration provides:

- **Milliarcsecond (mas) precision** for nutation and precession calculations
- **IAU-standard algorithms**: IAU 2000A, IAU 2000B, and IAU 2006/2000A models
- **Cached calculations** for performance optimization
- **Seamless fallback** when PyERFA is not installed

The precision improvement is most significant for:
- Dates far from J2000.0 (before 1950 or after 2050)
- Applications requiring arcsecond or better accuracy
- Professional astronomical calculations
