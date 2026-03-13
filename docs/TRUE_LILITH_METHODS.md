# True Lilith Calculation Methods

## Overview

LibEphemeris supports three lunar apogee calculation methods, each with different precision characteristics and use cases. This document describes the methods, their precision, and recommendations for method selection.

## Methods Comparison

### Mean Lilith (Mean Lunar Apogee)

Mean Lilith uses the mean orbital elements of the Moon to calculate the apogee position. It provides a smooth, continuously advancing position.

- **Precision**: ~12 arcsec mean difference (~0.003 degrees)
- **Calculation**: Based on mean lunar orbital elements
- **Constant**: `SE_MEAN_APOG`
- **Use case**: Traditional astrological interpretations

### True Lilith (Osculating Lunar Apogee)

True Lilith calculates the osculating (instantaneous) apogee of the lunar orbit using the **eccentricity vector** method. The function `calc_true_lilith` computes the direction of the osculating apogee from the instantaneous orbital elements derived from the Moon's position and velocity.

- **Precision**: ~52 arcsec mean difference, ~235 arcsec maximum difference
- **In degrees**: ~0.015 degrees mean, ~0.065 degrees maximum
- **Calculation**: Eccentricity vector from instantaneous position and velocity
- **Constant**: `SE_OSCU_APOG`
- **Use case**: Precise astrological work requiring the actual osculating apogee

This achieves **sub-arcminute** precision for mean difference, which is sufficient for most astrological applications.

### Interpolated Apogee

The interpolated apogee provides a smoothed version of the osculating apogee, removing short-period oscillations.

- **Precision**: ~1.1° difference from reference
- **Constant**: `SE_INTP_APOG`
- **Use case**: Research applications requiring smoothed osculating values

## Precision Measurement Methodology

Precision was measured using 500 random date samples across the range 1900-2050. Each sample compared libephemeris output against pyswisseph using the `SE_OSCU_APOG` body constant.

## Recommendation

For typical astrological use:

1. **Mean Lilith** (`SE_MEAN_APOG`) is recommended for traditional interpretations. It has the best precision (~0.003°) and produces smooth positions.

2. **True Lilith** (`SE_OSCU_APOG`) is recommended when the osculating apogee position is needed. The sub-arcminute mean precision is adequate for chart interpretation.

3. **Interpolated Apogee** (`SE_INTP_APOG`) is available for research requiring smoothed osculating values, but has lower precision (~1.1°).

## Technical Details

### Eccentricity Vector Method

The `calc_true_lilith` function implements the eccentricity vector method:

1. Obtain the Moon's geocentric position and velocity from DE440
2. Calculate the instantaneous orbital elements (semi-major axis, eccentricity)
3. Compute the eccentricity vector **e** = (v × h) / μ - r/|r|
4. The direction of **e** gives the osculating perigee; the apogee is 180° opposite
5. Project onto the ecliptic plane for the ecliptic longitude

This method is independent of analytical lunar theories and works directly from the numerical ephemeris data.

### pyswisseph Compatibility

The True Lilith calculation is compatible with pyswisseph's `SE_OSCU_APOG` body constant. Both libraries compute the osculating apogee, but use different internal methods, resulting in the documented precision differences.
