# Interpolated Lunar Apogee and Perigee

This document provides comprehensive documentation on the interpolated lunar apogee
(SE_INTP_APOG) and perigee (SE_INTP_PERG) calculations in libephemeris, explaining
how they differ from the osculating (true) apogee/perigee, and when to use each variant.

## Table of Contents

1. [Overview](#overview)
2. [The Problem with Osculating Apsides](#the-problem-with-osculating-apsides)
3. [The Interpolated Solution](#the-interpolated-solution)
4. [Implementation Details](#implementation-details)
5. [When to Use Each Variant](#when-to-use-each-variant)
6. [API Usage](#api-usage)
7. [Precision and Limitations](#precision-and-limitations)
8. [References](#references)

## Overview

The lunar apsides (apogee and perigee) represent the farthest and nearest points of
the Moon's orbit from Earth. In astrological tradition, the lunar apogee is known as
**Black Moon Lilith** (among other names), while the perigee is sometimes called
**Priapus** or **White Moon Selena** in some systems.

libephemeris provides three variants for each apsidal point:

| Body ID | Constant | Description |
|---------|----------|-------------|
| 12 | `SE_MEAN_APOG` | Mean Lunar Apogee (Mean Lilith) |
| 13 | `SE_OSCU_APOG` | Osculating Lunar Apogee (True Lilith) |
| 21 | `SE_INTP_APOG` | Interpolated Lunar Apogee (Natural Lilith) |
| 22 | `SE_INTP_PERG` | Interpolated Lunar Perigee (Natural Priapus) |

## The Problem with Osculating Apsides

### What is the Osculating Apogee?

The osculating apogee is computed from the Moon's instantaneous position and velocity
at any given moment. It represents the apogee of the hypothetical Keplerian ellipse
(two-body orbit) that would result if all gravitational perturbations suddenly ceased.

### The Two-Body Approximation Problem

The Moon's orbit is strongly perturbed by the Sun's gravity. When we compute "osculating"
orbital elements, we are forcing a three-body problem (Earth-Moon-Sun) into a two-body
framework (Earth-Moon only). This creates significant artifacts.

As the Swiss Ephemeris documentation explains:

> "The solar perturbation results in gigantic monthly oscillations in the ephemeris
> of the osculating apsides (the amplitude is **30 degrees**). These oscillations
> have to be considered an artifact of the insufficient model, they do not really
> show a motion of the apsides."

> "The lunar orbit is far from being an ellipse!"

### Characteristics of the Osculating Apogee

1. **Large oscillations**: The osculating apogee oscillates **+/- 30 degrees** around
   the mean apogee position with a period of approximately one synodic month (~29.5 days).

2. **Artificial motion**: These oscillations do not represent actual motion of the
   apsidal line. They are mathematical artifacts of the two-body approximation.

3. **True twice per month**: The osculating apogee only represents the "true" position
   (where the Moon is actually at maximum or minimum distance) twice per month:
   - When in conjunction with the Moon (Moon at maximum distance)
   - When in opposition to the Moon (Moon at minimum distance from Earth)

### Example: Osculating vs Mean Apogee

Over a typical month, you might observe:

```
Date          Mean Lilith    True Lilith    Difference
2025-01-01    115.23°        127.45°        +12.22°
2025-01-08    116.01°        98.31°         -17.70°
2025-01-15    116.79°        138.62°        +21.83°
2025-01-22    117.57°        95.24°         -22.33°
2025-01-29    118.35°        131.15°        +12.80°
```

The True Lilith (osculating) oscillates wildly while Mean Lilith moves smoothly.

## The Interpolated Solution

### Philosophy

The interpolated (or "natural") apogee attempts to find a middle ground between:
- The **mean apogee**: Too smooth, ignores real apsidal variations
- The **osculating apogee**: Too noisy, dominated by spurious oscillations

The goal is to smooth out the spurious 30-degree oscillations while preserving
the genuine apsidal motion.

### What Swiss Ephemeris Research Discovered

Through analysis of actual lunar apogee and perigee passages (when the Moon is at
its farthest/nearest distance from Earth), Swiss Ephemeris researchers found that:

1. **Apogee oscillates ~5 degrees from mean position** (vs. 30 degrees for osculating)
2. **Perigee oscillates ~25 degrees from mean position** (due to asymmetric solar effects)
3. **Apogee and perigee are not exactly opposite** - they deviate from 180 degrees apart,
   being only roughly opposite when the Sun is in conjunction with one of them or at
   a 90-degree angle

### The "Natural" Apsides

The interpolated apsides represent the "natural" position of the apsidal line -
where the actual lunar distance extremes occur. This is derived by smoothing the
osculating elements to remove the artificial oscillations while retaining the
physical variations.

## Implementation Details

### Swiss Ephemeris Approach

Swiss Ephemeris (from version 1.70) uses an interpolation method derived from
analytical lunar theory (Moshier's lunar ephemeris):

> "Conventional interpolation algorithms do not work well in the case of the lunar
> apsides. The supporting points are too far away from each other in order to
> provide a good interpolation, the error estimation is greater than 1 degree for
> the perigee. Therefore, [we] derived an 'interpolation method' from the analytical
> lunar theory which we have in the form of Moshier's lunar ephemeris."

### libephemeris Approach

libephemeris uses an **analytical perturbation series** based on ELP2000-82B lunar
theory to compute the interpolated apsides. This approach adds perturbation
corrections to the mean apogee (Mean Lilith) position.

#### Algorithm Overview

**For Interpolated Apogee:**
```
1. Calculate Mean Lilith (mean lunar apogee) longitude
2. Calculate Julian centuries from J2000.0
3. Compute fundamental lunar arguments (D, M, M', F)
4. Apply ELP2000-82B perturbation series:
   - Dominant term: +4.53° × sin(2D - 2M')
   - Additional periodic terms from lunar theory
5. Normalize result to [0°, 360°)
```

**For Interpolated Perigee:**
```
1. Calculate Mean Lilith + 180° (mean perigee longitude)
2. Calculate Julian centuries from J2000.0
3. Compute fundamental lunar arguments (D, M, M', F)
4. Apply calibrated perigee perturbation series:
   - Dominant term: -22.2° × sin(2D - 2M')  (opposite sign to apogee!)
   - ~15 additional terms fitted to Swiss Ephemeris data
5. Normalize result to [0°, 360°)
```

#### Key Insight: Asymmetric Perturbations

The apogee and perigee experience very different perturbation amplitudes:

| Apside | Oscillation from Mean | Dominant Coefficient |
|--------|----------------------|----------------------|
| Apogee | ~5° | +4.53° |
| Perigee | ~25° | -22.2° |

This asymmetry arises because:
1. Solar perturbations affect the perigee more strongly
2. The perturbation series terms have opposite signs for apogee vs perigee
3. The Moon spends more time near apogee (slower orbital velocity)

#### Perturbation Series Details

The perturbation series uses fundamental lunar arguments:
- **D**: Mean elongation of Moon from Sun
- **M**: Sun's mean anomaly
- **M'**: Moon's mean anomaly
- **F**: Moon's argument of latitude

**Apogee perturbation (excerpt):**
```python
# Dominant term (2D - 2M' argument)
delta = 4.5306 * sin(2*D - 2*M')  # degrees

# Additional terms
delta += 0.4193 * sin(2*D - M')
delta += 0.1320 * sin(2*M')
# ... ~10 more terms
```

**Perigee perturbation (calibrated to Swiss Ephemeris):**
```python
# Dominant term (opposite sign!)
delta = -22.2018 * sin(2*D - 2*M')  # degrees

# Additional terms fitted to SE data
delta += 1.5335 * E * sin(2*D - M)
delta += 1.1813 * sin(4*D - 2*M')
# ... ~15 more terms
```

#### Coefficient Calibration

The perigee perturbation coefficients were derived by:
1. Sampling 500 Swiss Ephemeris data points across a wide date range
2. Computing the residual (SE perigee - mean perigee - 180°)
3. Using least-squares fitting to find optimal coefficients
4. Validating against independent test dates

This approach achieves good agreement with Swiss Ephemeris while using
a computationally efficient analytical formula.

## When to Use Each Variant

### Mean Lilith (SE_MEAN_APOG) - `calc_mean_lilith`

**Use when:**
- You need smooth, predictable motion
- Historical or traditional astrological work
- Compatibility with older software that only supported mean values
- Applications where short-term oscillations are not relevant

**Characteristics:**
- Moves at approximately 40.7 degrees/year (prograde)
- No short-term oscillations
- Does not represent actual apsidal position on any given day

### True Lilith (SE_OSCU_APOG) - `calc_true_lilith`

**Use when:**
- You specifically need the instantaneous osculating elements
- Research into lunar orbital mechanics
- Understanding the theoretical two-body apogee

**Characteristics:**
- Large oscillations (+/- 30 degrees from mean)
- Period of approximately one synodic month
- Only "true" twice per month
- Most of the oscillation is artificial (two-body approximation artifact)

**Caution:** For most astrological purposes, the osculating apogee is misleading
because the 30-degree oscillations do not represent real apsidal motion.

### Interpolated Lilith (SE_INTP_APOG) - `calc_interpolated_apogee`

**Use when:**
- You want the "actual" apogee position (where Moon reaches maximum distance)
- Modern astrological work seeking physical accuracy
- Applications where you want real variations without spurious noise
- Transit or timing work related to lunar distance

**Characteristics:**
- Oscillates approximately +/- 5 degrees from mean (for apogee)
- Smooth, continuous motion
- Represents the genuine apsidal line orientation
- Available in libephemeris and Swiss Ephemeris (version 1.70+)

### Interpolated Perigee (SE_INTP_PERG) - `calc_interpolated_perigee`

**Use when:**
- You need the "actual" perigee position (where Moon reaches minimum distance)
- Supermoon calculations (perigee near Full Moon)
- Studying lunar distance variations

**Characteristics:**
- Oscillates approximately +/- 25 degrees from mean (larger than apogee)
- Not exactly opposite to interpolated apogee (can deviate from 180 degrees)
- Represents genuine perigee position

### Decision Matrix

| Requirement | Recommended |
|-------------|-------------|
| Smooth ephemeris for natal charts | Mean Lilith |
| Traditional/historical astrology | Mean Lilith |
| Physical accuracy | Interpolated |
| Supermoon calculations | Interpolated Perigee |
| Lunar distance extremes timing | Interpolated |
| Swiss Ephemeris compatibility (basic) | Mean Lilith |
| Swiss Ephemeris compatibility (advanced) | Interpolated |
| Research into orbital mechanics | Osculating |

## API Usage

### Using swe_calc_ut / swe_calc

```python
import libephemeris as swe

jd = 2460676.5  # 2025-01-01

# Mean Lilith
mean_pos, _ = swe.swe_calc_ut(jd, swe.SE_MEAN_APOG, swe.SEFLG_SPEED)
print(f"Mean Lilith: {mean_pos[0]:.4f} deg, speed: {mean_pos[3]:.4f} deg/day")

# Osculating (True) Lilith
oscu_pos, _ = swe.swe_calc_ut(jd, swe.SE_OSCU_APOG, swe.SEFLG_SPEED)
print(f"True Lilith: {oscu_pos[0]:.4f} deg, speed: {oscu_pos[3]:.4f} deg/day")

# Interpolated (Natural) Lilith
intp_pos, _ = swe.swe_calc_ut(jd, swe.SE_INTP_APOG, swe.SEFLG_SPEED)
print(f"Interpolated Lilith: {intp_pos[0]:.4f} deg, speed: {intp_pos[3]:.4f} deg/day")

# Interpolated Perigee
perg_pos, _ = swe.swe_calc_ut(jd, swe.SE_INTP_PERG, swe.SEFLG_SPEED)
print(f"Interpolated Perigee: {perg_pos[0]:.4f} deg, speed: {perg_pos[3]:.4f} deg/day")
```

### Using Direct Lunar Functions

```python
from libephemeris import lunar

jd_tt = 2460676.5  # 2025-01-01 (TT)

# Mean Lilith
mean_lon = lunar.calc_mean_lilith(jd_tt)
print(f"Mean Lilith: {mean_lon:.4f} deg")

# True Lilith (returns longitude, latitude, eccentricity)
oscu_lon, oscu_lat, oscu_ecc = lunar.calc_true_lilith(jd_tt)
print(f"True Lilith: {oscu_lon:.4f} deg, lat: {oscu_lat:.4f} deg, ecc: {oscu_ecc:.5f}")

# Interpolated Apogee
intp_lon, intp_lat, intp_ecc = lunar.calc_interpolated_apogee(jd_tt)
print(f"Interpolated Apogee: {intp_lon:.4f} deg")

# Interpolated Perigee
perg_lon, perg_lat, perg_ecc = lunar.calc_interpolated_perigee(jd_tt)
print(f"Interpolated Perigee: {perg_lon:.4f} deg")
```

### Return Values

All lunar apside functions return:
- **Longitude**: Ecliptic longitude in degrees [0, 360)
- **Latitude**: Ecliptic latitude in degrees (typically small, < 5 degrees)
- **Eccentricity/Distance**: Orbital eccentricity (~0.055 for Moon)

When using `swe_calc_ut` with `SEFLG_SPEED`, velocity is also calculated:
- **Speed (longitude)**: Daily motion in degrees/day
- **Speed (latitude)**: Daily change in latitude
- **Speed (distance)**: Daily change in eccentricity

## Precision and Limitations

### Comparison with Swiss Ephemeris

| Variant | Mean Error | Max Error |
|---------|------------|-----------|
| Mean Lilith | ~0.003° (~12") | ~0.005° (~18") |
| True Lilith | ~5° | ~15° |
| Interpolated Apogee | ~1.1° | ~3.3° |
| Interpolated Perigee | ~2.3° | ~8° |

**Note on True Lilith differences:** The ~5 degree mean differences arise because:
1. libephemeris computes osculating elements from JPL DE state vectors
2. Swiss Ephemeris uses integrated analytical lunar theory
3. The osculating apogee concept is inherently model-dependent for strongly perturbed orbits

**Note on Interpolated differences:** The remaining differences arise from:
1. Different underlying osculating apogee calculations feeding into the series
2. Coefficient calibration was done on a finite sample of dates
3. Swiss Ephemeris uses Moshier's full analytical lunar theory; libephemeris uses
   a reduced perturbation series with calibrated coefficients

### Smoothness Comparison

Testing over 100 consecutive days shows variance in daily motion:

| Variant | Daily Motion Variance |
|---------|----------------------|
| Osculating Apogee | ~38 deg²/day² |
| Interpolated Apogee | ~0.6 deg²/day² |

The interpolated apogee is approximately **60 times smoother** than the osculating.

### Apogee-Perigee Relationship

**Important:** Both Swiss Ephemeris and libephemeris compute apogee and perigee
using **independent** perturbation series. This means they are not constrained to
be exactly 180° apart.

Per Swiss Ephemeris documentation (section 2.2.4):
> "Apogee and perigee are not exactly opposite - they are only roughly opposite
> when the Sun is in conjunction with one of them or at a 90-degree angle."

The deviation from 180° can be up to **28 degrees** in extreme cases.

| Implementation | Apogee-Perigee Separation |
|----------------|---------------------------|
| libephemeris | ~180° ± 15° (physically correct) |
| Swiss Ephemeris | ~180° ± 28° (varies with lunar/solar geometry) |

This is expected physical behavior, not a limitation.

### Valid Date Range

The analytical perturbation series approach has no inherent date range limitations
beyond the validity of the underlying lunar theory.

For practical purposes, the implementation is valid for:
- DE421 (default): 1900 - 2050 (full precision)
- Extended range: 3000 BCE - 3000 CE (reduced precision)

## References

### Swiss Ephemeris Documentation

1. Section 2.2.3 "The Osculating Apogee" - Explains why the osculating apogee oscillates
2. Section 2.2.4 "The Interpolated or Natural Apogee and Perigee" - Describes the interpolation approach

### Academic References

1. Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs from 4000 B.C. to A.D. 8000" (1991), Willmann-Bell
2. Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Willmann-Bell, Chapter 47
3. Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
4. Brouwer, D. & Clemence, G.M. "Methods of Celestial Mechanics" (1961)

### Astrological References

1. Santoni, Francis. "Ephemerides de la lune noire vraie 1910-2010" (Editions St. Michel, 1993)
2. Koch, Dieter. "Was ist Lilith und welche Ephemeride ist richtig", Meridian 1/95
3. Koch, Dieter & Rindgen, Bernhard. "Lilith und Priapus" (Frankfurt/Main, 2000)

### Related Documentation

- [TRUE_LILITH_METHODS.md](TRUE_LILITH_METHODS.md) - Comparison of True Lilith calculation methods
- [PRECISION.md](PRECISION.md) - Overall precision limitations of libephemeris
