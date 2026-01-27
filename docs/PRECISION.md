# LibEphemeris Precision Limitations

This document provides detailed information about the precision and limitations of LibEphemeris compared to pyswisseph (Swiss Ephemeris). Use this guide to determine if LibEphemeris is suitable for your use case.

## Summary Table

| Component | Precision vs pyswisseph | Notes |
|-----------|------------------------|-------|
| Sun position | ±1 arcsec | High precision |
| Mercury, Venus positions | ±1 arcsec | High precision |
| Mars position | ±2 arcsec | Slightly relaxed due to perturbation differences |
| Moon position | ±5 arcsec | Relaxed due to lunar theory differences |
| Jupiter, Saturn positions | ±5 arcsec | Outer planets |
| Uranus, Neptune, Pluto | ±5 arcsec | Keplerian model limitations apply to Pluto |
| True Lunar Node | ~0.07° (~260 arcsec) | Osculating orbital elements method |
| House cusps | ±0.001° (~3.6 arcsec) | All 19 house systems |
| Ayanamsha (standard) | ±0.01° | Fagan-Bradley, Lahiri, etc. |
| Ayanamsha (star-based) | ±0.06° | True Citra, Galactic Center, etc. |
| Sun crossings (ingress) | ±0.001 arcsec | Sub-milliarcsecond precision |
| Moon crossings | ±0.05 arcsec | Sub-arcsecond precision |
| Planet crossings | ±0.1 arcsec | Sub-arcsecond precision |
| Polar houses (>66.5°) | Falls back to Porphyry | Placidus/Koch undefined |

---

## Planetary Positions

### Longitude Precision

Tested against pyswisseph at 1000+ random dates within DE421 range (1900-2050):

| Planet | Max Difference | Tolerance |
|--------|---------------|-----------|
| Sun | <1 arcsec | ±1 arcsec |
| Moon | <5 arcsec | ±5 arcsec |
| Mercury | <1 arcsec | ±1 arcsec |
| Venus | <1 arcsec | ±1 arcsec |
| Mars | <2 arcsec | ±2 arcsec |
| Jupiter | <5 arcsec | ±5 arcsec |
| Saturn | <5 arcsec | ±5 arcsec |
| Uranus | <5 arcsec | ±5 arcsec |
| Neptune | <5 arcsec | ±5 arcsec |
| Pluto | <5 arcsec | ±5 arcsec |

**Note on Moon**: The Moon has a relaxed tolerance due to differences in lunar theory implementations between LibEphemeris (Skyfield/JPL DE) and pyswisseph. Both are accurate, but use slightly different models.

**Note on Mars**: Mars requires slightly relaxed tolerance due to differences in perturbation calculations between the implementations.

### Latitude Precision

| Planet | Max Difference |
|--------|---------------|
| Sun | <1 arcsec |
| Moon | <5 arcsec |
| Other planets | <1 arcsec |

### Distance Precision

| Component | Max Difference |
|-----------|---------------|
| Distance (AU) | <0.0001 AU |
| Pluto distance | <0.001 AU (Keplerian model) |

### Velocity Precision

| Component | Max Difference |
|-----------|---------------|
| Angular velocity | <0.01°/day |
| Radial velocity | <0.001 AU/day |

---

## Coordinate Systems

### Geocentric vs Heliocentric

| Coordinate System | Longitude Tolerance |
|-------------------|---------------------|
| Geocentric | ±0.001° |
| Topocentric | ±0.001° |
| Heliocentric | ±0.03° (slightly relaxed) |
| Barycentric | ±0.03° (slightly relaxed) |

---

## House Systems

### Cusp Precision

LibEphemeris supports all 19 Swiss Ephemeris house systems:

| System | Precision | Notes |
|--------|-----------|-------|
| Placidus (P) | ~0.001-0.01° | Iterative, fails at polar latitudes |
| Koch (K) | ~0.001-0.01° | Similar to Placidus limitations |
| Porphyry (O) | ~0.01-0.1 arcsec | Geometric, works everywhere |
| Regiomontanus (R) | ~0.01° | Complex calculation |
| Campanus (C) | ~0.01° | Complex calculation |
| Equal (A/E) | Exact | No astronomical calculations |
| Whole Sign (W) | Exact | No astronomical calculations |
| Morinus (M) | ~0.01-0.1 arcsec | Geometric |
| Meridian (X) | ~0.01-0.1 arcsec | Geometric |
| Polich-Page/Topocentric (T) | ~0.01° | Complex calculation |
| Alcabitus (B) | ~0.01° | Medieval system |
| Vehlow (V) | Exact | Equal variant |
| Horizontal (H) | ~0.01° | Falls back to Porphyry if iteration fails |
| Gauquelin (G) | Approximation | Uses Placidus approximation (not true 36-sector) |
| Krusinski-Pisa (U) | ~0.01° | Uses Porphyry fallback |
| Carter Poli-Equatorial (F) | ~0.01° | |
| APC Houses (Y) | ~0.01° | Uses Porphyry fallback |
| Natural Graduation (N) | ~0.01° | |
| Sripati (S) | ~0.01° | |

### Angle Precision

| Angle | Max Difference |
|-------|---------------|
| Ascendant | <0.001° |
| MC (Midheaven) | <0.001° |
| ARMC | <0.001° |
| Vertex | <0.01° |
| Co-Ascendant | Not implemented (returns 0.0) |
| Polar Ascendant | Not implemented (returns 0.0) |

### Polar Latitude Limitations

**Critical Limitation**: Some house systems fail at polar latitudes where the ecliptic does not properly intersect the horizon.

| Latitude Range | Behavior |
|----------------|----------|
| |lat| < 66.5° | All house systems work normally |
| |lat| > 66.5° | Placidus, Koch may fail |
| |lat| > 85° | Most time-based systems unreliable |
| |lat| = 89.9° | Only Equal, Whole Sign, Porphyry reliable |

**Fallback Behavior**: When Placidus or Koch calculations fail at polar latitudes, LibEphemeris automatically falls back to Porphyry house system.

**Systems that work at all latitudes**:
- Equal (A/E)
- Whole Sign (W)
- Porphyry (O)
- Morinus (M)
- Meridian (X)

---

## Ayanamsha (Sidereal Modes)

### Standard Ayanamshas

LibEphemeris supports all 43 Swiss Ephemeris ayanamsha modes.

| Category | Max Difference | Examples |
|----------|---------------|----------|
| Formula-based | <0.01° | Fagan-Bradley, Lahiri, Raman |
| Epoch-based | <0.01° | J2000, J1900, B1950 |
| Historical | <0.01° | Babylonian variants |

### Star-Based Ayanamshas

Star-based ayanamshas have slightly higher tolerance due to differences in:
- Star position calculations
- Precession models
- Proper motion calculations
- Coordinate transformations

| Ayanamsha | Max Difference |
|-----------|---------------|
| True Citra | <0.06° |
| True Revati | <0.06° |
| True Pushya | <0.06° |
| True Mula | <0.06° |
| True Sheoran | <0.06° |
| Galactic Center 0 Sag | <0.06° |
| Galactic Center Rgilbrand | <0.06° |
| Galactic Center Cochrane | <0.06° |
| Galactic Center Mula Wilhelm | <0.06° |

---

## Crossing Events (Ingress Times)

### Sun Crossings

| Event | Precision |
|-------|-----------|
| Zodiac sign ingress | ±0.001 arcsec (sub-milliarcsecond) |
| Solstices | ±0.001 arcsec |
| Equinoxes | ±0.001 arcsec |
| Time difference vs pyswisseph | <1 second |

### Moon Crossings

| Event | Precision |
|-------|-----------|
| Zodiac sign ingress | ±0.05 arcsec (sub-arcsecond) |
| Node crossing (latitude=0) | ±0.05 arcsec |

### Planet Crossings

| Event | Precision |
|-------|-----------|
| Generic planet crossing | ±0.1 arcsec |
| Heliocentric crossing | ±0.1 arcsec |

**Note on Retrograde Stations**: Near retrograde stations where planetary velocity approaches zero, the algorithm uses Brent's method as a fallback to ensure convergence.

---

## Time Calculations

### Julian Day Precision

| Function | Precision |
|----------|-----------|
| julday/revjul | <1e-10 days (~0.01 microseconds) |
| Delta T | <0.5 seconds |

### Valid Date Ranges

LibEphemeris uses JPL DE ephemerides with specific date ranges:

| Ephemeris | Date Range | Precision |
|-----------|------------|-----------|
| DE421 (default) | 1900-2050 | Full precision |
| DE422 | -3000 to 3000 | Full precision |
| DE430 | 1550-2650 | Full precision |
| DE431 | -13200 to 17191 | Full precision |

**Outside valid range**: Calculations will raise an exception with the supported date range.

---

## Minor Bodies and Extended Points

### Lunar Nodes and Lilith

LibEphemeris uses different calculation models for some lunar points, resulting in varying precision:

| Point | Max Difference vs pyswisseph | Notes |
|-------|------------------------------|-------|
| Mean Node | < 0.01° | High precision |
| True Node | ~0.07° (~260 arcsec) | Osculating orbital elements method |
| Mean Lilith | ~0.1° | Minor formula differences |
| True Lilith | ~5-15° | Different orbital model (see below) |

**Note on True Node**: The True Node calculation uses osculating orbital elements derived from the Moon's instantaneous position and velocity, with IAU 2006 precession correction. Mean error is ~0.02° (86 arcsec), max error ~0.07° (260 arcsec). For most astrological purposes, this precision is sufficient.

**Note on True Lilith**: The True Lilith (osculating apogee) calculation uses the eccentricity vector method derived from JPL DE ephemeris state vectors. This differs fundamentally from Swiss Ephemeris which uses an integrated analytical lunar theory. The difference arises because the "osculating apogee" concept is model-dependent when solar perturbations are significant. Typical error is 5-15°, with occasional differences up to 25°. For applications requiring closer Swiss Ephemeris compatibility, use Mean Lilith instead.

### Asteroids

| Asteroid | Precision | Notes |
|----------|-----------|-------|
| Chiron | ~1° | Keplerian approximation |
| Pholus | ~1° | Keplerian approximation |
| Ceres | ~1° | Keplerian approximation |
| Pallas | ~1° | Keplerian approximation |
| Juno | ~1° | Keplerian approximation |
| Vesta | ~1° | Keplerian approximation |

**Note**: Asteroids use Keplerian (two-body) orbital propagation which does not account for planetary perturbations. Errors grow over time, especially for dates far from the orbital elements epoch (2023.0).

### Trans-Neptunian Objects (TNOs)

| TNO | Precision | Notes |
|-----|-----------|-------|
| Orcus | ~5° | Keplerian approximation |
| Haumea | ~5° | Keplerian approximation |
| Quaoar | ~5° | Keplerian approximation |
| Makemake | ~5° | Keplerian approximation |
| Gonggong | ~5° | Keplerian approximation |
| Eris | ~5° | Keplerian approximation, highly distant |
| Sedna | ~10° | Keplerian approximation, extremely distant |

**Note**: TNOs use simplified Keplerian propagation. For high-precision TNO work, use dedicated ephemeris files.

---

## Fixed Stars

### Position Precision

| Component | Precision |
|-----------|-----------|
| Longitude | ~0.01 arcsec |
| Latitude | ~0.01 arcsec |
| Proper motion | Corrected for stellar curvature |
| Nutation | IAU 2000A model (1365 terms) |

### Velocity Limitation

Fixed star velocity calculations return 0:

```python
pos, _ = ephem.fixstar_ut("Aldebaran", jd, SEFLG_SPEED)
# pos[3], pos[4], pos[5] are 0.0 (velocities not implemented)
```

---

## Eclipse Functions

### Implemented Features

| Function | Precision |
|----------|-----------|
| Solar eclipse timing | ~1-2 minutes |
| Lunar eclipse timing | ~1-2 minutes |
| Occultation timing | ~1-2 minutes |
| Eclipse magnitude | ~0.01 |

### Not Implemented

The following return 0 or placeholder values:

- Saros series number
- Inex number
- Sunrise/sunset on central line for solar eclipses

---

## Performance Considerations

### Calculation Speed

LibEphemeris (pure Python with Skyfield) is approximately 10-100x slower than pyswisseph (C library). For most astrological applications, this is negligible. For batch processing of millions of charts, consider:

1. Using multi-threading with `EphemerisContext`
2. Caching frequently-used calculations
3. Using pyswisseph for performance-critical paths

### Memory Usage

| Resource | Approximate Size |
|----------|-----------------|
| DE421 ephemeris | ~16 MB |
| DE431 ephemeris | ~3.4 GB |
| Timescale data | ~2 MB |
| Per EphemerisContext | ~1 KB |

Ephemeris files are shared across all `EphemerisContext` instances for memory efficiency.

---

## Summary of Use Case Suitability

### Suitable For

- Natal chart calculations (standard latitudes)
- Transit calculations
- Synastry and composite charts
- Progression calculations
- Most astrological software applications
- Research requiring reproducible, documented precision

### Use With Caution

- Polar latitude locations (>66.5°) - use Equal or Whole Sign houses
- Star-based ayanamshas - expect ±0.06° variance
- Asteroid/TNO positions - Keplerian approximation only
- Very old historical dates (before 1900 with DE421)

### Not Suitable For

- Observatory-grade positional astronomy (use NOVAS or SOFA)
- Spacecraft navigation (use SPICE)
- Precise eclipse predictions for scientific purposes
- High-precision Saros/Inex series analysis

---

## References

1. NASA JPL DE Ephemeris documentation
2. Swiss Ephemeris documentation and source code
3. Meeus, "Astronomical Algorithms" 2nd Edition
4. IERS Conventions 2003 (nutation models)
5. IAU SOFA Library documentation
6. Skyfield documentation (Rhodes Mill)
