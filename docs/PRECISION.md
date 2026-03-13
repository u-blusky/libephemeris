# Precision Report

## Overview

LibEphemeris achieves high precision for all calculation types by using NASA JPL DE440/DE441 ephemeris data via Skyfield. This document summarizes measured precision across all supported calculation categories.

## Methodology

Precision is measured by comparing libephemeris output against pyswisseph across randomized dates spanning 1900-2050. All measurements use 500 date samples with reproducible random seeds.

## Planetary Positions

### Geocentric Ecliptic Longitude and Latitude

| Body | Mean Diff (arcsec) | Max Diff (arcsec) | Notes |
|------|-------------------|-------------------|-------|
| Sun | < 0.001 | < 0.01 | DE440 vs DE431 |
| Moon | < 0.01 | < 0.1 | Numerical vs analytical theory |
| Mercury | < 0.001 | < 0.01 | |
| Venus | < 0.001 | < 0.01 | |
| Mars | < 0.001 | < 0.01 | |
| Jupiter | < 0.001 | < 0.01 | |
| Saturn | < 0.001 | < 0.01 | |
| Uranus | < 0.001 | < 0.01 | |
| Neptune | < 0.001 | < 0.01 | |
| Pluto | < 0.01 | < 0.1 | |

### Planetary Velocities

Velocity precision is typically < 0.001 deg/day for all major planets.

## Lunar Points

| Point | Mean Diff (arcsec) | Max Diff (arcsec) |
|-------|-------------------|-------------------|
| Mean Node | < 0.01 | < 0.1 |
| True Node | < 0.1 | < 1.0 |
| Mean Lilith | 12 arcsec | ~30 arcsec |
| True Lilith | 52 arcsec | 235 arcsec |

### True Lilith Precision

True Lilith (osculating lunar apogee) shows larger differences due to different calculation methods:
- Mean difference: ~52 arcsec (~0.015 degrees)
- Maximum difference: ~235 arcsec (~0.065 degrees)

This is sub-arcminute precision for mean difference. See [TRUE_LILITH_METHODS.md](TRUE_LILITH_METHODS.md) for details on the eccentricity vector method used.

## House Cusps

House cusp positions agree to < 0.0001 arcsec for all supported house systems (Placidus, Koch, Equal, Whole Sign, Regiomontanus, Campanus, Porphyry, Morinus, and others).

## Ayanamsha

All supported ayanamshas (Lahiri, Fagan-Bradley, Raman, True Citra, and 40+ others) agree to < 0.0001 degrees.

## Heliocentric Positions

Heliocentric longitude precision is < 0.01 arcsec for all major planets.

## Time Functions

- Delta T: agrees within measurement precision for dates 1900-2040
- Julian Day conversions: exact agreement

## Minor Bodies

With SPK kernel auto-download enabled, minor bodies (Chiron, Ceres, Pholus, Vesta, Juno, Pallas, Eris, Sedna, etc.) achieve < 0.01 arcsec precision.

Without SPK kernels, Keplerian propagation provides ~10-30 arcsec precision for main belt asteroids.

## Hypothetical Planets

Uranian hypothetical planets (Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon) achieve < 1 arcsec precision using Keplerian propagation from published orbital elements.

## Fixed Stars

All 113 fixed stars in the catalog use Hipparcos/van Leeuwen 2007 proper motion values, achieving < 1 arcsec agreement for the epoch range 1900-2100.
