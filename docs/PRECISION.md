# Precision Summary

Compact precision summary for LibEphemeris vs pyswisseph 2.10.03.
For full details, models, and methodology see [reference/precision.md](reference/precision.md)
and [reference/swisseph-comparison.md](reference/swisseph-comparison.md).

## Planetary Positions (geocentric ecliptic, 1550–2650 CE)

| Body | Mean Diff | Max Diff | Notes |
|------|-----------|----------|-------|
| Sun | 0.04" | 0.20" | DE440 vs DE431 |
| Moon | 0.70" | 3.32" | Numerical vs analytical lunar theory |
| Mercury | 0.05" | 0.32" | |
| Venus | 0.08" | 0.33" | |
| Mars | 0.06" | 0.58" | |
| Jupiter | 0.12" | 0.44" | Includes COB correction |
| Saturn | 0.13" | 0.51" | Includes COB correction |
| Uranus | 0.23" | 0.50" | |
| Neptune | 0.24" | 1.17" | |
| Pluto | 0.26" | 0.75" | Includes COB correction |

All planets sub-arcsecond. Moon ~3" max reflects different lunar models (JPL DE440 numerical integration vs ELP/MPP02 + DE431).

## Velocities

| Component | Max Diff |
|-----------|----------|
| Longitude speed | < 0.003°/day |
| Latitude speed | < 0.004°/day |
| Distance speed | < 0.0001 AU/day |

## Lunar Points

| Point | Max Diff | Independent Verification |
|-------|----------|--------------------------|
| Mean Node | < 0.001° | — |
| True Node | < 0.01" | Verified vs JPL Horizons to machine precision |
| Mean Lilith | < 0.015" (lon) | Latitude ~20" systematic (different node formulas) |
| True Lilith | < 0.5" | Both libraries ~240" from Horizons (inherent two-body limit) |
| Interpolated Apogee | ~0.36° | Genuine algorithm difference (JPL DE440 vs ELP2000-82B perturbation series) |
| Interpolated Perigee | ~2.6° | JPL DE440 physical passages vs truncated ELP2000-82B perturbation series |

## House Cusps

< 0.02" for all 24 supported house systems, tested at 11 global locations.
Iterative systems (Placidus, Koch) use 10⁻⁷° convergence threshold.

## Fixed Stars

116 stars from Hipparcos catalog with van Leeuwen 2007 proper motions.
Max difference: 0.51" (Rigil Kentaurus — nearest star, parallax not modeled).
98% of 101 comparable stars within 0.5". Two catalog bugs found and fixed
(Algedi wrong component, Asellus Borealis wrong HIP number).

## Ayanamsha

- Standard modes (Lahiri, Fagan-Bradley, Raman): < 0.0002°
- Star-based modes (True Citra, True Revati): < 0.006°

## Eclipses

- Solar eclipse timing: < 6 seconds
- Lunar eclipse timing: < 8 seconds
- Rise/set/transit: < 30 seconds

## Delta T

- Modern (1900–2025): < 1 second
- Historical (< 1700): up to ~187 seconds (different models: SMH 2016 vs E&M 2006)
- Future (> 2050): grows with extrapolation divergence

## Minor Bodies

- With SPK kernels: sub-arcsecond (matching JPL Horizons)
- Keplerian fallback: ~10–30" for main belt asteroids near epoch, degrees over decades

## Hypothetical Planets

Uranian hypothetical planets: < 1" (Keplerian from published elements).

## Heliocentric / Barycentric / Equatorial / XYZ

| Mode | Max Diff |
|------|----------|
| Heliocentric | < 0.0004° (1.1") |
| Barycentric (non-Sun) | < 0.001° |
| Equatorial RA/Dec | < 0.0005° (1.7") |
| XYZ Cartesian | < 0.00005 AU |

## Hyper-Validation

4400+ comparison rounds across 29 API sections:
**3947 PASS, 441 KNOWN, 0 FAIL, 12 SKIP**.
All divergences documented in [divergences.md](divergences.md).
