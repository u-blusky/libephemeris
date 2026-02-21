# Interpolated Perigee Methodology

## Overview

LibEphemeris computes the interpolated lunar perigee using a three-layer architecture calibrated against JPL DE441 ephemeris:

1. **Mean perigee** — polynomial model for the secular apsidal precession
2. **Harmonic perturbation series** — 61-term trigonometric series capturing periodic oscillations
3. **Residual correction table** — precomputed corrections absorbing remaining model error

The combined result achieves ~1.5 degree RMS agreement with Swiss Ephemeris (down from ~10-11 degrees in v1), while being grounded entirely in JPL ephemeris data rather than analytical lunar theory.

## Background: What is Interpolated Perigee?

The Moon's orbit around Earth is continuously perturbed by the Sun. These perturbations cause:

1. The orbital axis (line of apsides) to precess with a period of ~8.85 years
2. The perigee (point of closest approach) to **oscillate** around this mean precession

These oscillations are substantial: approximately ±25 degrees from the mean position.

The term "interpolated perigee" refers to a **smoothed** version of the instantaneous (osculating) perigee position, with short-period oscillations removed. There are different valid approaches to defining and computing this smoothing.

## Calibration Method: Passage-Interpolated Harmonic Fitting (v2.2)

### The Problem

The lunar perigee oscillates violently in its osculating (instantaneous) form — ±30 degrees within days. Extracting clean perturbation coefficients from this noisy signal requires care. Two earlier approaches failed:

- **v2.0 (raw osculating fit):** Fitting a harmonic series directly to ~365K daily osculating perigee samples. Failed because osculating noise is *correlated* with the design matrix terms — both are functions of the same Delaunay arguments (D, M, M', F). This produced nonsensical coefficients (9.18 deg passage RMS).

- **v2.1 (passage-only fit):** Fitting only at perigee passages (Earth-Moon distance minima), where the Moon's longitude *is* the perigee longitude. Failed because at every passage, the mean anomaly M' ≈ 0 by definition, making the design matrix degenerate (condition number 3.7×10⁵).

### The Solution (v2.2)

The v2.2 method combines the strengths of both approaches:

1. **Find perigee passages** — Locate all Earth-Moon distance minima over a 1000-year span [1500, 2500] CE using JPL DE441. At each passage, the Moon's ecliptic longitude is an unambiguous measure of the perigee longitude (12,958 passages found).

2. **Cubic spline interpolation** — Fit a cubic spline through the (time, perigee_longitude) passage data, with angle unwrapping to handle the 0°/360° discontinuity. This creates a smooth, continuous perigee longitude function.

3. **Daily resampling** — Sample the spline at daily intervals to produce ~365K clean data points with full coverage of all Delaunay argument combinations (unlike passage-only data where M' ≈ 0 always).

4. **Harmonic least-squares fit** — Construct a 120-term design matrix of candidate trigonometric terms from lunar theory and fit via ordinary least squares. Terms with |coefficient| < 0.001 degrees are discarded, yielding the final 61-term series.

### Why This Works

The spline interpolation acts as a physically-motivated smoothing: it passes exactly through the ground-truth passage points while providing well-conditioned data at arbitrary intermediate times. The daily samples have full M' coverage (unlike passages) and no correlated osculating noise (unlike raw data).

## Architecture

### Layer 1: Mean Perigee (Polynomial)

The mean perigee position uses a standard polynomial in Julian centuries T from J2000:

```
ω₀ = 83.3532° + 4069.0137° × T + ...  (higher-order terms)
```

This captures the secular apsidal precession (~40.7° per year). The mean perigee is 180° from the mean apogee (Mean Lilith).

### Layer 2: Perturbation Series (61 terms)

The perturbation series `_calc_elp2000_perigee_perturbations()` adds periodic corrections organized by physical origin:

| Category | Terms | Dominant Coefficient |
|----------|-------|---------------------|
| Primary evection sin(kD - kM') | 13 | sin(2D-2M') = **-22.21°** |
| Evection phase cos(kD - kM') | 7 | cos(2D-2M') = -0.075° |
| Solar anomaly coupling (E×sin) | 9 | E×sin(4D-4M'-M) = +0.53° |
| Solar double coupling (E²×sin) | 4 | E²×sin(2D-2M'+2M) = +0.071° |
| Lunar anomaly (sin kM') | 2 | sin(M') = +0.011° |
| Latitude coupling (F terms) | 3 | sin(2F-2M') = +0.17° |
| Cross-coupling (D,M' combos) | 7 | sin(6D-5M') = -0.45° |
| Solar-latitude cross | 3 | E×sin(2F-2D-M) = +0.010° |
| Higher-order evection-solar | 3 | E×sin(8D-8M'-M) = +0.038° |
| Secular (T×trig) | 3 | T×cos(2D-2M') = -0.004° |
| Cosine phase corrections | 2 | cos(2F-2M') = +0.022° |
| Sun-Moon anomaly coupling | 2 | E×sin(M-M') = -0.002° |
| Polynomial corrections | 4 | const = -0.175° |

The polynomial correction terms (const, T, T², T³) absorb secular drift between the mean perigee model and the JPL-calibrated series.

### Layer 3: Residual Correction Table

After applying the perturbation series, residual errors (~2 deg RMS near J2000, ~8 deg over the full range) are absorbed by a precomputed correction table (`lunar_corrections.py`, 15,195 entries). Linear interpolation between table entries brings the final error to < 0.1 degrees.

The correction table is generated by `scripts/generate_lunar_corrections.py`, which:
1. Computes JPL ground-truth perigee via ±28 day quadratic regression of osculating elements
2. Computes model perigee (mean + perturbation series)
3. Stores the difference at regular time intervals

## Key Coefficient Changes (v1 → v2.2)

The v1 coefficients were severely attenuated because v1 used quadratic regression smoothing of osculating data, which damped the harmonic amplitudes by ~50%:

| Term | v1 Coefficient | v2.2 Coefficient | Change |
|------|---------------|-------------------|--------|
| sin(2D-2M') | -9.62° | **-22.21°** | ×2.3 |
| sin(4D-4M') | +0.94° | **+6.45°** | ×6.9 |
| sin(6D-6M') | -0.13° | **-2.28°** | ×17.5 |
| sin(2D-M') | +2.61° | **-0.04°** | sign flip |
| sin(M') | +0.72° | **+0.01°** | ×0.01 |
| E×sin(2D-2M'-M) | -0.33° | **-0.97°** | ×3.0 |

The dramatic increase in the dominant sin(2D-2M') term (evection) from -9.62° to -22.21° is the most significant change. The old value was roughly half the true amplitude due to signal attenuation by the quadratic regression window.

## Precision

| Metric | v1 | v2.2 |
|--------|-----|------|
| Perturbation series RMS (near J2000) | ~10-11° | ~2° |
| Perturbation series RMS (full range) | ~10-11° | ~8° |
| After correction table (vs SE) | ~10-11° | **~1.5°** |
| After correction table (vs JPL) | N/A | **< 0.1°** |

The remaining ~1.5° difference vs Swiss Ephemeris is expected: it reflects the fundamental methodological difference between LibEphemeris (JPL-grounded numerical approach) and Swiss Ephemeris (ELP2000-82B analytical theory). See "Methodological Differences" below.

## Methodological Differences vs Swiss Ephemeris

### Swiss Ephemeris Approach (Moshier/ELP2000-82B)

**Method:** Semi-analytical perturbation theory

- Uses the ELP2000-82B lunar theory developed by Chapront-Touzé and Chapront
- Separates perturbations **analytically** based on their physical origin
- Removes "short-period" perturbations (periods less than a few months)
- Retains "long-period" perturbations that represent the true apsidal motion

### LibEphemeris Approach (JPL DE441 + Passage-Interpolated Fitting)

**Method:** Empirical harmonic series calibrated against numerical ephemeris

- Identifies perigee passages from JPL DE441 Moon positions (physical ground truth)
- Spline interpolation between passages creates a smooth perigee longitude function
- Harmonic series coefficients are fit empirically to this function
- Residual correction table absorbs remaining model error

### Why ~1.5° Residual Difference Exists

The two approaches differ in **what is smoothed** and **how**:

| Aspect | Swiss Ephemeris | LibEphemeris |
|--------|----------------|--------------|
| Ground truth | Analytical lunar theory | JPL numerical ephemeris |
| Smoothing | Theory-based term selection | Passage-interpolated spline |
| Perturbation definition | Physics-derived | Empirically fitted |
| Extended range | ~-5400 to +5400 | -13200 to +17191 (DE441) |

The ~1.5° difference is a consequence of different smoothing philosophies applied to the same underlying phenomenon, not an error in either implementation.

## Analogy

Consider audio signal processing with bass and treble components:

- **SE approach:** Uses an intelligent filter that knows which frequencies represent "signal" vs "noise" based on the source characteristics
- **LibEphemeris approach:** Identifies the signal at known clean points (passage times), interpolates between them, then fits a harmonic model to the result

Both produce a smoothed result, but the outputs differ because the smoothing methodology differs.

## Implications

### For Astronomy

Both approaches are scientifically valid. The choice depends on the application:

- **SE-style:** Better when studying long-term orbital dynamics within analytical theory frameworks
- **JPL-style:** Better when needing consistency with actual Moon positions and extended time ranges

### For Astrology

Swiss Ephemeris is the **de facto standard** in astrological software. The ~1.5° RMS difference between LibEphemeris and SE for interpolated perigee is within astrological orb tolerances for most applications, though users should be aware of it.

### For Software Compatibility

The v2.2 recalibration reduced the SE discrepancy from ~10-11° to ~1.5°, bringing LibEphemeris much closer to practical compatibility with Swiss Ephemeris for this function.

## Date Range

| Ephemeris | Valid Range |
|-----------|-------------|
| Swiss Ephemeris | ~-5400 to +5400 |
| JPL DE440 | 1550 to 2650 |
| JPL DE441 | -13200 to +17191 |

Set `LIBEPHEMERIS_EPHEMERIS=de441.bsp` for extended date range coverage.

## Calibration Reproducibility

The calibration can be reproduced using:

```bash
LIBEPHEMERIS_EPHEMERIS=de441.bsp python scripts/calibrate_perigee_perturbations.py \
    --start-year 1500 --end-year 2500 --output /tmp/perigee_v22_full.json
```

The script outputs calibrated coefficients as Python code ready to paste into `lunar.py`. The correction table can be regenerated with:

```bash
python scripts/generate_lunar_corrections.py
```

## References

- Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris", Astronomy & Astrophysics, 190, 342-352
- Chapront-Touzé, M. & Chapront, J. (1991). "Lunar Tables and Programs from 4000 B.C. to A.D. 8000", Willmann-Bell
- Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441", The Astronomical Journal, 161:105
- Moshier, S.L. "Ephemeris Calculation Software" (1989-2023). Available at http://www.moshier.net — independent C implementation of ELP2000-82B used as reference for the analytical series coefficients.
