# REBOUND Integration Benefits

## Overview

LibEphemeris optionally integrates with [REBOUND](https://rebound.readthedocs.io/), an N-body integrator, and [ASSIST](https://github.com/matthewholman/assist), which provides high-accuracy asteroid ephemerides using JPL planetary ephemeris data. When installed, these packages provide improved precision for minor body (asteroid) calculations beyond simple Keplerian propagation.

## Installation

REBOUND and ASSIST are optional dependencies:

```bash
pip install libephemeris[nbody]
```

Or explicitly:

```bash
pip install rebound assist
```

## Precision Comparison

Minor body position precision depends on the propagation method:

| Method | Precision | Best For |
|--------|-----------|----------|
| Keplerian (default) | ~10-30 arcsec (main belt) | Quick calculations, short propagation |
| REBOUND IAS15 (2-body) | ~1-10 arcsec | Medium-term propagation |
| REBOUND + ASSIST | < 0.1 arcsec | High-precision research |
| SPK kernels (JPL) | < 0.01 arcsec | Best available precision |

### Keplerian Propagation

The default Keplerian method uses published orbital elements with secular perturbation corrections. Precision degrades with propagation time from the epoch.

### REBOUND N-body Integration

REBOUND provides symplectic and high-order integrators for N-body problems. The IAS15 integrator (15th-order Gauss-Radau) provides machine-precision accuracy for the gravitational problem.

### ASSIST Ephemeris-Driven Integration

ASSIST uses JPL planetary ephemeris data to drive the integration, accounting for gravitational perturbations from all major planets, the Moon, and selected asteroids.

## Quantified Benefits

### Propagation Error Growth

| Propagation Time | Keplerian | REBOUND (2-body) | REBOUND + ASSIST |
|-----------------|-----------|-------------------|------------------|
| 30 days | ~1 arcsec | ~0.1 arcsec | < 0.01 arcsec |
| 1 year | ~10 arcsec | ~1 arcsec | < 0.1 arcsec |
| 10 years | ~100 arcsec | ~10 arcsec | < 1 arcsec |
| 100 years | ~1000 arcsec | ~100 arcsec | < 10 arcsec |

Note: SPK kernel auto-download (when available) bypasses propagation entirely and provides the best precision.

## REBOUND Functions Available

When REBOUND is installed, the following functions are available via `libephemeris.rebound_integration`:

- `check_rebound_available()` — Check if REBOUND is installed
- `check_assist_available()` — Check if ASSIST is installed
- `get_rebound_version()` — Get REBOUND version string
- `get_assist_version()` — Get ASSIST version string
- `elements_to_rebound_orbit(elements, jd)` — Convert orbital elements to REBOUND format
- `propagate_orbit_rebound(elements, jd_start, jd_end)` — Propagate orbit using REBOUND
- `propagate_trajectory(elements, jd_start, jd_end, num_points)` — Generate trajectory points
- `compare_with_keplerian(elements, jd)` — Compare REBOUND vs Keplerian results

### Integrators

REBOUND supports multiple integrators:

- **IAS15** — 15th-order Gauss-Radau, adaptive step size, machine precision (default)
- **WHFast** — Symplectic Wisdom-Holman, fast for long-term integrations
- **MERCURIUS** — Hybrid symplectic/IAS15, handles close encounters
- **TRACE** — Time-reversible hybrid integrator
- **Leapfrog** — Simple 2nd-order symplectic integrator

### Data Types

- `ReboundIntegrator` — Enum of available integrators
- `AssistEphemConfig` — Configuration for ASSIST ephemeris data directory
- `PropagationResult` — Result dataclass with (x, y, z, vx, vy, vz, jd_tt) and properties for `distance`, `ecliptic_lon`, `ecliptic_lat`, and `to_tuple()`

## When to Use REBOUND/ASSIST

### Recommended for:
- **Long-term asteroid tracking** — propagation over years to decades
- **Centaur and TNO calculations** — objects with significant perturbations (Chiron, Eris)
- **Close encounter scenarios** — MERCURIUS handles gravitational close approaches
- **Research applications** — when sub-arcsecond precision is needed without SPK kernels

### Not necessary for:
- **Short-term calculations** — Keplerian is adequate for < 30 days
- **Main planets** — DE440/DE441 ephemeris provides direct positions
- **When SPK kernels are available** — SPK auto-download provides better precision

## Summary

REBOUND and ASSIST provide significant precision improvements for minor body calculations, especially for long propagation times and objects with complex orbital dynamics (centaurs, TNOs). For most use cases, the SPK kernel auto-download feature provides even better precision. REBOUND is most valuable when SPK kernels are unavailable and high-precision propagation is needed.
