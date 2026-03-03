# PyERFA Integration

LibEphemeris integrates PyERFA (the Python wrapper for the IAU SOFA/ERFA library) to provide IAU-standard nutation and precession models with milliarcsecond (mas) precision for coordinate transformations.

## Table of Contents

- [Background](#background)
- [Method](#method)
  - [Nutation Models](#nutation-models)
  - [Error Growth Over Time](#error-growth-over-time)
  - [Obliquity](#obliquity)
  - [Precession-Nutation Matrix](#precession-nutation-matrix)
  - [Cached Nutation](#cached-nutation)
- [API Reference](#api-reference)
- [Precision and Validation](#precision-and-validation)
- [When to Use PyERFA](#when-to-use-pyerfa)
- [Installation](#installation)
- [References](#references)

## Background

Precession and nutation are fundamental corrections required for transforming celestial coordinates between epochs. Precession describes the slow, secular drift of Earth's rotational axis over millennia, while nutation describes the shorter-period oscillations superimposed on this drift, driven primarily by the gravitational torques of the Moon and Sun on Earth's equatorial bulge.

The International Astronomical Union (IAU) has adopted progressively refined models for these effects. The IAU 2000A nutation model uses ~1300 terms to achieve ~0.2 mas precision. The IAU 2006 precession model improves upon the IAU 2000 precession by incorporating updated rate corrections derived from VLBI observations. PyERFA provides reference implementations of these models, enabling LibEphemeris to offer observatory-grade accuracy as an optional precision tier.

## Method

### Nutation Models

LibEphemeris provides access to three nutation models via PyERFA, each offering a different precision-performance tradeoff:

| Model | Precision | Terms | Use Case |
|-------|-----------|-------|----------|
| IAU 2000A | ~0.2 mas | ~1300 | Highest precision nutation |
| IAU 2000B | ~1 mas | ~80 | High precision, faster computation |
| IAU 2006/2000A (nut06a) | ~0.01-0.05 mas | ~1300 | Best overall precision with improved precession |

The nut06a model combines the IAU 2006 precession with IAU 2000A nutation, providing the best available precision for equinox-based coordinates. It applies small corrections to the IAU 2000A nutation to account for the updated precession rates.

### Error Growth Over Time

The differences between models grow with distance from J2000.0:

| Years from J2000 | IAU 2000B vs 2000A | nut06a vs nut00a |
|------------------|-------------------|------------------|
| 0 | ~0.2 mas | ~0.01 mas |
| 10 | ~0.5 mas | ~0.02 mas |
| 50 | ~2 mas | ~0.05 mas |
| 100 | ~5 mas | ~0.1 mas |

### Obliquity

PyERFA provides the IAU 2006 obliquity model with milliarcsecond accuracy:

```python
from libephemeris.erfa_nutation import get_erfa_obliquity_iau2006
obliquity = get_erfa_obliquity_iau2006(jd_tt)  # Returns radians
```

### Precession-Nutation Matrix

For complete coordinate transformations from J2000 to a date of interest:

```python
from libephemeris.erfa_nutation import get_erfa_pnm06a_matrix
# Returns the precession-nutation matrix for J2000 to date transformation
matrix = get_erfa_pnm06a_matrix(jd_tt)
```

### Cached Nutation

For performance-critical applications, a cached variant avoids redundant calculations for the same Julian date via an LRU cache:

```python
from libephemeris.erfa_nutation import get_erfa_nutation_cached
dpsi, deps = get_erfa_nutation_cached(jd_tt)
```

## API Reference

LibEphemeris exposes the following PyERFA-based functions:

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

## Precision and Validation

### Nutation Precision

- **IAU 2000A**: Full-precision nutation model (~0.2 mas accuracy)
- **IAU 2000B**: Truncated model (~1 mas accuracy), approximately 10x faster
- **nut06a**: IAU 2006 precession + IAU 2000A nutation (~0.05 mas accuracy)

The precision improvement is most significant for:
- Dates far from J2000.0 (before 1950 or after 2050)
- Applications requiring arcsecond or better accuracy
- Professional astronomical calculations

## When to Use PyERFA

### Recommended Use Cases

1. **High-precision coordinate transformations**: When milliarcsecond accuracy is required for precessing coordinates between epochs.

2. **Observatory-grade calculations**: Professional astronomical applications requiring IAU-standard algorithms.

3. **Validation against external sources**: When comparing results with JPL Horizons or other professional ephemeris services.

4. **Long-term extrapolation**: Calculations more than 50 years from J2000.0 benefit most from IAU 2006 improvements.

### When the Default Is Sufficient

1. **Typical astrological applications**: The default precision (~1 arcsec) is adequate for most astrological chart calculations.

2. **Performance-critical code**: The built-in approximation is faster when milliarcsecond precision is not required.

3. **Dates near J2000.0**: Within 50 years of J2000.0, the differences between models are minimal for most applications.

## Installation

PyERFA is an optional dependency. Install it with:

```bash
pip install pyerfa
```

Or with LibEphemeris extras:

```bash
pip install libephemeris[precision]
```

To verify PyERFA is available:

```python
from libephemeris.erfa_nutation import has_erfa
if has_erfa():
    print("PyERFA is available for high-precision calculations")
```

When PyERFA is not installed, LibEphemeris falls back seamlessly to its built-in approximations.

## References

- IAU SOFA Library: [http://www.iausofa.org/](http://www.iausofa.org/)
- Capitaine, N. et al. "Expressions for IAU 2000 precession quantities" (2003), Astronomy & Astrophysics, 412, 567-586
- Mathews, P.M. et al. "Modeling of nutation and precession: New nutation series for nonrigid Earth and insights into the Earth's interior" (2002), Journal of Geophysical Research, 107(B4)
- IERS Conventions (2010), IERS Technical Note No. 36
