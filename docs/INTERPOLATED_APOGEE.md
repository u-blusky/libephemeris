# Interpolated Lunar Apogee and Perigee

## Overview

The **interpolated apogee** (SE_INTP_APOG, body 21) and **interpolated perigee**
(SE_INTP_PERG, body 22) are smoothed versions of the osculating lunar apogee and
perigee that filter out short-period oscillations while preserving the underlying
apsidal motion.

## Three Variants of Lunar Apogee

LibEphemeris computes three distinct apogee variants:

### 1. Mean Apogee (SE_MEAN_APOG, body 12)

The time-averaged apogee position computed from analytical polynomial formulae.
Moves smoothly at the mean apsidal precession rate (~40.7°/year prograde).

### 2. Osculating Apogee / True Lilith (SE_OSCU_APOG, body 13)

The instantaneous apogee direction derived from the Moon's osculating orbital
elements. Oscillates approximately **30 degrees** from the mean position due to
solar perturbations, evection, and other short-period effects. The oscillation
amplitude of ~30° makes the osculating apogee highly variable on timescales of
weeks to months.

### 3. Interpolated Apogee (SE_INTP_APOG, body 21)

A smoothed apogee that oscillates approximately **5 degrees** from the mean
position — substantially less than the 30° oscillation of the osculating variant.
Computed by sampling the osculating apogee over a 56-day window (9 equally-spaced
samples) and applying a linear fit to extract the smoothed trend.

## Algorithm

The interpolated apogee is computed as follows:

1. **Sample** the osculating apogee (True Lilith) at 9 equally-spaced points
   across a 56-day window centered on the target date
2. **Unwrap** the longitude values to handle 360°/0° crossings
3. **Fit** a linear regression to the unwrapped longitudes
4. **Evaluate** the fitted line at the target date to get the smoothed longitude
5. **Latitude and distance** are taken from the osculating apogee at the target date

The 56-day window spans approximately two anomalistic months, which effectively
averages out the dominant short-period oscillations while preserving the secular
apsidal precession and medium-period perturbations.

## Interpolated Perigee (SE_INTP_PERG, body 22)

The interpolated perigee uses the same smoothing algorithm but samples the
osculating perigee (computed directly from the eccentricity vector, not by
adding 180° to the apogee). The interpolated perigee longitude is approximately
opposite to the interpolated apogee (within ~2-3°), and its latitude has the
opposite sign.

## Comparison with pyswisseph

The interpolated apogee/perigee implementation differs from pyswisseph's Swiss
Ephemeris due to fundamental algorithmic differences:

- LibEphemeris uses JPL DE440/DE441 ephemeris via Skyfield for the osculating
  elements, while Swiss Ephemeris uses its own semi-analytical lunar theory
- The ELP2000-82B perturbation terms used internally produce slightly different
  osculating positions, which propagate to the interpolated result
- Typical longitude divergence: ~1° (up to ~1.5° in worst cases)
- This is classified as an inherent engine difference in the hyper-validation

## API Usage

```python
import libephemeris

# Interpolated Apogee (body 21)
result = libephemeris.swe_calc_ut(2451545.0, 21)
lon, lat, dist = result[0][0], result[0][1], result[0][2]

# Interpolated Perigee (body 22)
result = libephemeris.swe_calc_ut(2451545.0, 22)
lon, lat, dist = result[0][0], result[0][1], result[0][2]

# Direct function calls
from libephemeris.lunar import calc_interpolated_apogee, calc_interpolated_perigee

lon, lat, dist = calc_interpolated_apogee(2451545.0)
lon, lat, dist = calc_interpolated_perigee(2451545.0)
```

## References

- Chapront, J. et al. (2002) "A new determination of lunar orbital parameters",
  A&A 387, 700-708
- ELP2000-82B lunar theory (Chapront-Touzé & Chapront, 1983)
- Simon, J.L. et al. (1994) "Numerical expressions for precession formulae",
  A&A 282, 663-683
- Swiss Ephemeris documentation: "Interpolated Lunar Apogee and Perigee"
