# Getting Started

## Installation

```bash
pip install libephemeris
```

The library works immediately after installation. When no local ephemeris files are present, it fetches data from NASA JPL Horizons API transparently.

For offline use or faster calculations:

```bash
libephemeris download:medium      # DE440 (~128 MB), 1550-2650 — recommended
libephemeris download:base        # DE440s (~35 MB), 1850-2150 — lightweight
libephemeris download:extended    # DE441 (~3.3 GB), -13200 to +17191 — full range
```

### Optional Extras

```bash
pip install libephemeris[nbody]   # REBOUND/ASSIST n-body integration for TNOs
pip install libephemeris[stars]   # Star catalog building (astropy)
pip install libephemeris[all]     # Everything
```

**Requirements:** Python 3.9+ | skyfield >= 1.54 | pyerfa >= 2.0

## Ephemeris Tiers

| Tier | JPL File | Date Range | Size | Use Case |
|------|----------|------------|------|----------|
| `base` | de440s.bsp | 1849-2150 | ~31 MB | Modern-era, lightweight |
| `medium` | de440.bsp | 1550-2650 | ~128 MB | General purpose **(default)** |
| `extended` | de441.bsp | -13200 to +17191 | ~3.1 GB | Historical/far-future |

DE440 and DE441 have identical precision -- DE441 is the extended-range version.

```python
import libephemeris as swe

# Select by tier
swe.set_precision_tier("extended")

# Or by file directly
swe.set_ephemeris_file("de441.bsp")

# Or via environment variable
# export LIBEPHEMERIS_PRECISION=extended
# export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

## First Calculation

```python
import libephemeris as swe
from libephemeris.constants import SE_SUN, SE_MOON, SE_MARS, SEFLG_SPEED

jd = swe.julday(2024, 3, 26, 12.0)

# Planet positions with velocity
sun, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SPEED)
moon, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SPEED)
mars, _ = swe.calc_ut(jd, SE_MARS, SEFLG_SPEED)

# Returns: (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
print(f"Sun:  {sun[0]:.4f} deg, speed {sun[3]:.4f} deg/day")
print(f"Moon: {moon[0]:.4f} deg, speed {moon[3]:.4f} deg/day")
print(f"Mars: {mars[0]:.4f} deg, speed {mars[3]:.4f} deg/day")
```

## House Cusps

```python
jd = swe.julday(2024, 11, 5, 18.0)

# Placidus houses for Rome (41.9N, 12.5E)
cusps, ascmc = swe.houses(jd, 41.9028, 12.4964, b"P")

print(f"ASC: {ascmc[0]:.4f}")
print(f"MC:  {ascmc[1]:.4f}")
for i, cusp in enumerate(cusps[1:13], 1):
    print(f"House {i:2d}: {cusp:.4f}")
```

## Thread Safety

The global API uses mutable global state (compatible with PySwissEph).
For concurrent calculations, use `EphemerisContext`:

```python
from libephemeris import EphemerisContext, SE_SUN

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)
pos, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)
```

Each context has independent state (sidereal mode, topographic position, LEB reader).

## Minor Bodies

For asteroids and TNOs, JPL SPK kernels provide high precision.
Auto-download is supported:

```python
swe.set_auto_spk_download(True)
pos, _ = swe.calc_ut(2460000.0, swe.SE_CHIRON, 0)
```

Fallback chain: SPK kernel > auto-download > REBOUND/ASSIST > Keplerian propagation.

See [REBOUND integration](methodology/rebound-integration.md) for n-body details.
