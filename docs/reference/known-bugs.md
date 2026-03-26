# Known Bugs and Limitations

Documented issues in LibEphemeris. These are tracked and will be addressed
in future releases where possible.

## Active Bugs

### 1. True Node distance tolerance

**Severity:** Low (distance value only, longitude is correct)
**Affected:** `swe_calc_ut(jd, SE_TRUE_NODE, ...)` — the distance component.

**Symptom:** Distance value differs from reference by ~2.3e-4 AU. Longitude
is correct within ~4" (expected difference due to different lunar models).

**Root cause:** Different osculating orbit computation from Moon state vectors.
The distance represents the osculating distance at the node, which is sensitive
to the exact eccentricity vector computation.

**Impact:** Negligible for astrological use. Longitude is the primary output
for True Node calculations.

## Fixed Bugs (v1.0.0a3)

### Sun heliocentric (FIXED)

Previously, `swe_calc(jd, SE_SUN, SEFLG_HELCTR)` returned ~126° longitude
(garbage from `arctan2` of near-zero values). Now correctly returns (0,0,0,0,0,0).

### Uranian geocentric bodies (FIXED)

Previously, `swe_calc(jd, 40, SEFLG_SPEED)` without `SEFLG_HELCTR` raised
`UnknownBodyError`. Geocentric conversion path added for all 8 Uranians (40-47).

## Horizons Backend Limitations

These are inherent to the Horizons API approach, not bugs:

| Limitation | Details |
|------------|---------|
| **Moon velocity** | ~0.001 deg/day difference vs analytical derivative (numerical differentiation) |
| **Heliocentric offset** | ~0.01-0.03" systematic difference (Horizons Sun center vs Skyfield SSB offset) |
| **Requires internet** | No offline fallback (use LEB or Skyfield for offline) |
| **Latency** | ~300-600ms first chart, cached afterward |
| **Unsupported bodies** | True Node (11), Oscu Apogee (13), Interp Apogee/Perigee (21,22), fixed stars, planetary moons — fall through to Skyfield |

Full details: [Horizons Backend Guide](../architecture/horizons-backend.md)
