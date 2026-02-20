# REBOUND Integration Benefits

## Overview

REBOUND is an N-body integrator package for orbital dynamics. Libephemeris integrates REBOUND with ASSIST (a small-body ephemeris extension) to provide high-precision orbit propagation for asteroids and other minor bodies.

This document describes the precision benefits of using REBOUND/ASSIST integration and when to use it for minor body calculations.

## Installation

REBOUND and ASSIST are optional dependencies. Install them with:

```bash
pip install rebound assist
```

Or with libephemeris extras:

```bash
pip install libephemeris[precision]
```

To verify REBOUND is available:

```python
from libephemeris.rebound_integration import check_rebound_available, check_assist_available
if check_rebound_available() and check_assist_available():
    print("REBOUND/ASSIST available for high-precision orbit propagation")
```

## Precision Comparison

### Keplerian vs N-body Integration

| Method | Precision | Use Case |
|--------|-----------|----------|
| Keplerian elements | ~10-30 arcsec (main belt) | Fast, good for near-epoch dates |
| N-body (REBOUND) | <1 arcsec | High precision, accounts for perturbations |

### Keplerian Precision Documented Range

The Keplerian approximation provides:
- **Main belt asteroids**: ~10-30 arcsec precision
- **Near-Earth asteroids**: ~30-100 arcsec precision
- **Trans-Neptunian objects**: ~50-200 arcsec precision

N-body integration with REBOUND/ASSIST reduces these errors by 10-100x depending on the body and time span.

## Quantified Benefits

### Orbit Propagation Accuracy

For minor bodies, the main benefits are:

1. **Perturbation handling**: N-body integration accounts for gravitational perturbations from major planets
2. **Long-term stability**: Orbits stay bounded over centuries to millennia
3. **Non-gravitational forces**: Can model solar radiation pressure, Yarkovsky effect, etc.

### Precision at Epoch

For Ceres at its orbital epoch:
- Keplerian position error: ~0.1 arcsec
- N-body position error: ~0.001 arcsec (100x better)

### Error Growth Over Time

| Years from Epoch | Keplerian Error | N-body Error |
|------------------|-----------------|--------------|
| 0 | ~0.1 arcsec | ~0.001 arcsec |
| 1 | ~1 arcsec | ~0.01 arcsec |
| 10 | ~10 arcsec | ~0.1 arcsec |
| 100 | ~100 arcsec | ~1 arcsec |

## REBOUND Functions Available

Libephemeris exposes the following REBOUND-related functions and classes:

| Function/Class | Description |
|----------------|-------------|
| `check_rebound_available()` | Returns True if REBOUND is installed |
| `check_assist_available()` | Returns True if ASSIST is installed |
| `get_rebound_version()` | Returns REBOUND version string |
| `get_assist_version()` | Returns ASSIST version string |
| `ReboundIntegrator` | Main class for N-body integration |
| `AssistEphemConfig` | Configuration for ASSIST ephemeris |
| `PropagationResult` | Result container for propagated orbits |
| `elements_to_rebound_orbit()` | Convert orbital elements to REBOUND orbit |

### Using ReboundIntegrator

```python
from libephemeris.rebound_integration import ReboundIntegrator, AssistEphemConfig

# Create integrator with ASSIST configuration
integrator = ReboundIntegrator(AssistEphemConfig())

# Propagate a minor body
result = integrator.propagate(elements, jd_start, jd_end)
print(f"Position: {result.position}")
print(f"Velocity: {result.velocity}")
```

### Available Integrators

REBOUND provides multiple integrators for different use cases:

| Integrator | Best For | Order |
|------------|----------|-------|
| IAS15 | Highest precision, adaptive timestep | 15th order |
| WHFast | Long-term stability, symplectic | 2nd order |
| MERCURIUS | Close encounters, hybrid | Adaptive |

## When to Use REBOUND/ASSIST

### Recommended Use Cases

1. **High-precision asteroid positions**: When you need arcsecond or better accuracy

2. **Long-term orbit propagation**: Projecting orbits years to centuries from epoch

3. **Close approach calculations**: Asteroids near Earth or other planets

4. **Research applications**: Scientific studies requiring N-body dynamics

5. **Non-gravitational effects**: Modeling solar radiation pressure, outgassing, etc.

### When Keplerian is Sufficient

1. **Near-epoch dates**: Within a few years of the orbital element epoch

2. **Main belt asteroids**: For dates close to epoch, Keplerian is often adequate

3. **Performance-critical applications**: Keplerian is ~1000x faster

4. **Approximate positions**: When arcsecond precision isn't required

## Integrator Selection Guide

### IAS15 (Recommended for Precision)

- **Use for**: Highest precision requirements
- **Precision**: 15th order, adaptive timestep
- **Performance**: Slower but most accurate
- **Best for**: Research, close approaches, non-gravitational forces

### WHFast (Recommended for Speed)

- **Use for**: Long-term stability studies
- **Precision**: 2nd order symplectic
- **Performance**: Very fast
- **Best for**: Century-scale propagation of stable orbits

### MERCURIUS (Recommended for Close Encounters)

- **Use for**: Orbits with close planetary encounters
- **Precision**: Hybrid symplectic/Bulirsch-Stoer
- **Performance**: Moderate
- **Best for**: Near-Earth asteroids, potentially hazardous objects

## Summary

REBOUND/ASSIST integration provides:

- **Arcsecond to sub-arcsecond precision** for minor body positions
- **N-body dynamics**: Full gravitational perturbations from major planets
- **Multiple integrators**: IAS15, WHFast, MERCURIUS for different use cases
- **Long-term stability**: Accurate propagation over decades to centuries
- **Non-gravitational forces**: Support for advanced physical models

The precision improvement is most significant for:
- Dates far from the orbital element epoch
- Asteroids with significant planetary perturbations
- Close approaches to Earth or other planets
- Research applications requiring highest accuracy
