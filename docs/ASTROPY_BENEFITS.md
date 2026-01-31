# Astropy Integration Benefits

This document describes what astropy could provide to supplement libephemeris calculations, and recommendations for when to use it.

## Overview

Astropy is the core Python package for astronomy, providing comprehensive tools for time handling, coordinate transformations, and astronomical constants. LibEphemeris evaluates astropy's `astropy.time` and `astropy.coordinates` modules as optional supplements to its Skyfield-based core.

While Skyfield provides the primary ephemeris calculations in libephemeris, astropy offers specialized features not directly available in Skyfield, particularly Galactic coordinates, extended time scales, and barycentric corrections.

## Installation

```bash
# Install astropy separately (evaluation/development)
pip install astropy

# Or install libephemeris with star catalog building support
pip install libephemeris[stars]
```

Note: Astropy is currently used for the star catalog building script and is being evaluated for additional integrations documented here.

## Feature Comparison

### Time Scales

Astropy provides additional time scales beyond what Skyfield offers:

| Time Scale | Skyfield | Astropy | Use Case |
|------------|----------|---------|----------|
| UTC | Yes | Yes | Civil time |
| UT1 | Yes | Yes | Earth rotation time |
| TT | Yes | Yes | Terrestrial dynamical time |
| TAI | Limited | Yes | Atomic time (basis for TT) |
| TDB | Limited | Yes | Barycentric dynamical time |
| TCG | No | Yes | Geocentric coordinate time |
| TCB | No | Yes | Barycentric coordinate time |
| GPS | No | Yes | GPS navigation time |

### Time Scale Visualization

```
                    Civil Time                  Physical Time
                        |                            |
    GPS ----------+     |     +-------------------- TCB (Barycentric)
                  |     |     |
    UTC ----+-----+-----+     +-------------------- TCG (Geocentric)
            |                 |
    TAI ----+     +-----+-----+-------------------- TDB (Barycentric Dynamical)
            |     |     |
            +-----+-----+----- TT (Terrestrial Time)
                  |
    UT1 ----------+------------ Earth Rotation
```

### Coordinate Frames

Astropy provides coordinate frames not directly available in Skyfield:

| Frame | Skyfield | Astropy | Description |
|-------|----------|---------|-------------|
| ICRS | Yes | Yes | International Celestial Reference System |
| Equatorial (RA/Dec) | Yes | Yes | Right Ascension / Declination |
| Ecliptic | Yes | Yes | Sun-centered plane of solar system |
| Horizontal (Alt/Az) | Yes | Yes | Observer-centered sky coordinates |
| Galactic | No | Yes | Galaxy-centered (l, b) coordinates |
| Supergalactic | No | Yes | Local supercluster centered |
| ITRS | Limited | Yes | Earth-fixed terrestrial frame |
| GCRS | Yes | Yes | Geocentric celestial reference |

## Quantified Benefits

### Time Scale Precision

Delta T (TT - UT1) comparison at J2000.0:

| Source | Delta T | Model |
|--------|---------|-------|
| Skyfield | ~63.8 seconds | Stephenson et al. 2016 |
| Astropy | ~63.8 seconds | IERS data tables |
| Difference | <100 milliseconds | Model differences |

### Coordinate Transform Accuracy

Alt/Az comparison for the same ICRS coordinates:

| Comparison | Typical Difference | Notes |
|------------|-------------------|-------|
| Skyfield vs Astropy Alt | <60 arcseconds | Different nutation models |
| Skyfield vs Astropy Az | <60 arcseconds | Different nutation models |

For most astrological applications, these differences are negligible.

## Astropy Functions Available

LibEphemeris exposes astropy capabilities through the `astropy_integration` module:

### Time Conversions

```python
from libephemeris.astropy_integration import (
    check_astropy_available,
    compare_time_conversions,
    get_extended_time_scales,
    parse_time_string,
)

# Check if astropy is available
if check_astropy_available():
    jd = 2451545.0  # J2000.0
    
    # Compare time handling between Skyfield and astropy
    result = compare_time_conversions(jd)
    print(f"Delta T difference: {result['delta_t_difference_ms']:.2f} ms")
    
    # Get all time scales
    scales = get_extended_time_scales(jd)
    for name, jd_value in scales.items():
        print(f"{name}: {jd_value:.10f}")
    
    # Parse time strings flexibly
    jd = parse_time_string("2000-01-01T12:00:00")
```

### Coordinate Transformations

```python
from libephemeris.astropy_integration import (
    compare_coordinate_transforms,
    icrs_to_galactic,
    galactic_to_icrs,
)

# Compare Alt/Az between Skyfield and astropy
result = compare_coordinate_transforms(
    ra_deg=101.2875,   # Sirius RA
    dec_deg=-16.7161,  # Sirius Dec
    jd_utc=2451545.0,
    lon_deg=12.5,      # Observer longitude
    lat_deg=41.9       # Observer latitude
)
print(f"Alt difference: {result['alt_difference_arcsec']:.2f} arcsec")

# Convert to Galactic coordinates (not available in Skyfield)
galactic_l, galactic_b = icrs_to_galactic(266.417, -29.008)
print(f"Galactic: l={galactic_l:.2f}, b={galactic_b:.2f}")

# Convert back to ICRS
ra, dec = galactic_to_icrs(galactic_l, galactic_b)
```

### Barycentric Corrections

```python
from libephemeris.astropy_integration import get_barycentric_correction

# Calculate barycentric corrections for high-precision work
result = get_barycentric_correction(
    ra_deg=180.0,      # Target RA
    dec_deg=45.0,      # Target Dec
    jd_utc=2451545.0,
    lon_deg=0.0,       # Observer longitude
    lat_deg=50.0       # Observer latitude
)

print(f"BJD (TDB): {result['bjd_tdb']:.6f}")
print(f"Light time: {result['light_time_seconds']:.2f} s")
print(f"RV correction: {result['rv_correction_km_s']:.2f} km/s")
```

### Geodetic Utilities

```python
from libephemeris.astropy_integration import (
    geodetic_to_geocentric,
    geocentric_to_geodetic,
)

# Convert geodetic (WGS84) to geocentric Cartesian
x, y, z = geodetic_to_geocentric(
    lon_deg=12.5,    # Longitude
    lat_deg=41.9,    # Geodetic latitude
    height_m=100.0   # Height above ellipsoid
)
print(f"Geocentric: x={x/1e6:.3f} Mm, y={y/1e6:.3f} Mm, z={z/1e6:.3f} Mm")

# Convert back
lon, lat, height = geocentric_to_geodetic(x, y, z)
```

### Integration Evaluation

```python
from libephemeris.astropy_integration import evaluate_integration_potential

# Get comprehensive evaluation of astropy benefits
evaluation = evaluate_integration_potential()

print("Time features astropy provides:")
for feature in evaluation['time_features']:
    print(f"  - {feature}")

print("\nCoordinate features astropy provides:")
for feature in evaluation['coordinate_features']:
    print(f"  - {feature}")
```

## When to Use Astropy

### Recommended Use Cases

1. **Galactic Coordinate Work**
   - Converting between ICRS and Galactic (l, b)
   - Supergalactic coordinate transformations
   - Not available in Skyfield

2. **Extended Time Scales**
   - TDB (Barycentric Dynamical Time) for solar system barycenter calculations
   - TCG/TCB for relativistic corrections
   - GPS time for navigation applications

3. **Barycentric Corrections**
   - High-precision radial velocity measurements
   - Pulsar timing
   - BJD (Barycentric Julian Date) calculations
   - Light travel time corrections

4. **Flexible Time Parsing**
   - ISO format strings
   - FITS date/time formats
   - Unix timestamps
   - Various astronomical time formats

5. **Geodetic Calculations**
   - WGS84 ellipsoid-based conversions
   - Precise geodetic to geocentric transformations

### Not Necessary For

1. **Standard Astrological Charts**
   - Natal chart calculations
   - Transit predictions
   - Synastry analysis
   
   For these applications, Skyfield provides all required functionality with excellent precision.

2. **Basic Ephemeris Generation**
   - Planetary positions
   - Lunar phases
   - House calculations
   
   These are fully handled by libephemeris's Skyfield-based core.

3. **Standard Ecliptic/Equatorial Coordinates**
   - Skyfield already handles these transformations

4. **Delta T Calculations**
   - Skyfield's Delta T is sufficient for most purposes
   - Differences are typically <1 second

## Performance Considerations

Astropy adds overhead compared to Skyfield-only calculations:

| Operation | Skyfield Only | With Astropy |
|-----------|---------------|--------------|
| Time scale conversion | ~50 microseconds | ~200 microseconds |
| Coordinate transform | ~100 microseconds | ~500 microseconds |
| Galactic conversion | N/A | ~300 microseconds |
| Barycentric correction | N/A | ~1 millisecond |

For batch processing, consider caching or vectorization.

## Compatibility Notes

Skyfield and astropy are designed to work together:

- Both use ICRS as the fundamental reference frame
- Both use IERS data for Earth orientation parameters
- Delta T values may differ slightly due to different models
- Nutation models may have small differences (~0.1 arcsec)
- For most astrological uses, differences are negligible

## Implementation Approach

The recommended approach for astropy integration in libephemeris:

1. **Optional Dependency**: Install via `pip install libephemeris[stars]` or manually
2. **Evaluation Module**: Use `astropy_integration.py` for comparison/testing
3. **Skyfield Primary**: Core calculations remain Skyfield-based for performance
4. **Astropy Supplements**: Use astropy for specialized features not in Skyfield
5. **Graceful Fallback**: Functions check for astropy availability

## Comparison with PyERFA

| Aspect | PyERFA | Astropy |
|--------|--------|---------|
| Focus | Nutation/precession precision | Coordinate frames, time scales |
| Installation | `pip install libephemeris[precision]` | `pip install astropy` |
| Adds | IAU 2006/2000A nutation | Galactic coords, TDB, barycentric |
| Performance | Minimal overhead | Moderate overhead |
| Required | No | No |
| Recommended for | Research, high precision | Galactic work, barycentric corrections |

## Summary

| Aspect | Skyfield Only | With Astropy |
|--------|---------------|--------------|
| Time scales | UTC, UT1, TT | + TAI, TDB, TCG, TCB, GPS |
| Galactic coordinates | No | Yes |
| Barycentric corrections | Manual | Built-in |
| Time string parsing | Manual | Flexible |
| Geodetic (WGS84) | Limited | Full |
| Performance | Optimal | Some overhead |
| Required for astrology | No | No |
| Recommended for research | Sometimes | Yes (specialized cases) |

Astropy provides valuable supplementary capabilities for specialized astronomical work, particularly Galactic coordinates and barycentric corrections. For standard astrological calculations, libephemeris's Skyfield-based core provides all necessary functionality with optimal performance.

## References

1. Astropy Collaboration, 2022, "The Astropy Project: Sustaining and Growing a Community-oriented Open-source Project and the Latest Major Release (v5.0) of the Core Package", ApJ 935, 167
2. Astropy documentation: https://docs.astropy.org/
3. Skyfield documentation: https://rhodesmill.org/skyfield/
4. IERS Conventions 2010, Chapter 5: Transformation between celestial and terrestrial reference systems
5. Galactic coordinate system: https://en.wikipedia.org/wiki/Galactic_coordinate_system
