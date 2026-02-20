# Interpolated Perigee: Methodological Differences

## Overview

The interpolated lunar perigee computed by LibEphemeris differs from Swiss Ephemeris by approximately 10-11 degrees RMS. This document explains why this difference exists and what it represents.

## Background: What is Interpolated Perigee?

The Moon's orbit around Earth is continuously perturbed by the Sun. These perturbations cause:

1. The orbital axis (line of apsides) to precess with a period of ~8.85 years
2. The perigee (point of closest approach) to **oscillate** around this mean precession

These oscillations are substantial: approximately ±25 degrees from the mean position.

The term "interpolated perigee" refers to a **smoothed** version of the instantaneous (osculating) perigee position, with short-period oscillations removed. However, there are different valid approaches to defining and computing this smoothing.

## Two Methodological Approaches

### Swiss Ephemeris Approach (Moshier/ELP2000-82B)

**Method:** Semi-analytical perturbation theory

**How it works:**
- Uses the ELP2000-82B lunar theory developed by Chapront-Touzé and Chapront
- Separates perturbations **analytically** based on their physical origin
- Removes "short-period" perturbations (periods less than a few months)
- Retains "long-period" perturbations that represent the true apsidal motion

**Physical interpretation:**
> The perigee position after removing rapid oscillations caused by instantaneous Sun-Moon-Earth geometry, while preserving perturbations intrinsic to the orbital dynamics.

**Characteristics:**
- Perturbation-aware smoothing
- Distinguishes between different perturbation sources
- Based on analytical lunar theory
- Standard in astrological software

### LibEphemeris Approach (JPL DE441 + Quadratic Regression)

**Method:** Numerical smoothing of osculating elements

**How it works:**
- Computes instantaneous (osculating) perigee from JPL DE441 Moon positions
- The osculating perigee oscillates violently (±30° within days)
- Applies quadratic polynomial regression over 9 samples in a ±28 day window
- The result is a numerically smoothed perigee position

**Physical interpretation:**
> The local trend of the perigee position, with high-frequency components removed via uniform mathematical smoothing.

**Characteristics:**
- Perturbation-agnostic smoothing
- Uses actual high-precision JPL Moon positions
- Based on numerical curve fitting
- Covers extended date range (DE441: -13200 to +17191)

## Quantitative Difference

| Metric | Value |
|--------|-------|
| RMS difference | ~10-11 degrees |
| Maximum difference | ~18-20 degrees |
| Range of agreement | Varies with lunar phase |

The difference is not an error in either implementation, but a consequence of applying different smoothing methodologies to the same underlying phenomenon.

## Why the Difference Occurs

The two approaches differ in **what is removed** and **what is retained**:

| Perturbation | Approx. Period | SE retains? | JPL-smooth retains? |
|--------------|----------------|-------------|---------------------|
| Evection (2D-2M') | ~31.8 days | Yes (dominant) | Yes (but attenuated) |
| Annual equation | ~1 year | Yes | Partially |
| Monthly variation | ~14.77 days | Removed | Partially removed |
| Higher harmonics | Variable | Theory-dependent | Uniform smoothing |

- **SE approach:** Selectively removes perturbations based on physical understanding of their origin and period
- **LibEphemeris approach:** Applies uniform frequency-based smoothing without distinguishing perturbation types

## Analogy

Consider audio signal processing with bass and treble components:

- **SE approach:** Uses an intelligent filter that knows which frequencies represent "signal" vs "noise" based on the source characteristics
- **LibEphemeris approach:** Uses a uniform low-pass filter that removes everything above a certain frequency, affecting both signal and noise together

Both produce a smoothed result, but the outputs differ because the filtering logic differs.

## Implications

### For Astronomy

Both approaches are scientifically valid. The choice depends on the application:

- **SE-style:** Better when studying long-term orbital dynamics and apsidal precession
- **JPL-style:** Better when needing consistency with actual Moon positions over extended time ranges

### For Astrology

Swiss Ephemeris is the **de facto standard** in astrological software. Astrological interpretations and ephemeris tables worldwide use SE's definition of interpolated perigee. This means:

- Software using LibEphemeris will show interpolated perigee positions that differ from other astrological software
- The difference is systematic and predictable, not random error
- Historical or comparative astrological work may require awareness of this difference

### For Software Compatibility

Libraries claiming "Swiss Ephemeris compatibility" should note that this specific function produces different results. The ~11° RMS difference exceeds typical precision requirements for astrological applications.

## Date Range Considerations

| Ephemeris | Valid Range |
|-----------|-------------|
| Swiss Ephemeris | ~-5400 to +5400 |
| JPL DE441 | -13200 to +17191 |

The LibEphemeris approach provides perigee calculations over a significantly extended range, but this comes at the cost of divergence from the SE standard within SE's valid range.

## References

- Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris", Astronomy & Astrophysics, 190, 342-352
- Chapront-Touzé, M. & Chapront, J. (1991). "Lunar Tables and Programs from 4000 B.C. to A.D. 8000", Willmann-Bell
- Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441", The Astronomical Journal, 161:105
- Moshier, S.L. (implementation of ELP2000-82B in Swiss Ephemeris)
