# Flag Reference

Flags are bitmasks that control what is calculated and how results are returned.
`calc_ut()` / `calc()` return a 6-element tuple `(longitude, latitude, distance,
speed_lon, speed_lat, speed_dist)` and a return flag. Combine flags with `|`.

## Velocity

| Flag | Effect |
|------|--------|
| `SEFLG_SPEED` | Populate speed fields (pos[3]-pos[5]) with daily motion. Without this, speeds are zero. Almost every call should include it. |

## Observer

By default the observer is at Earth's center (geocentric).

| Flag | Effect |
|------|--------|
| `SEFLG_HELCTR` | Heliocentric: observer at the Sun. Distances become heliocentric AU. |
| `SEFLG_BARYCTR` | Barycentric: observer at the solar system barycenter. |
| `SEFLG_TOPOCTR` | Topocentric: observer on Earth's surface. Set position with `set_topo(lon, lat, alt)`. Matters most for the Moon (~1 deg parallax). |

## Coordinates

By default output is ecliptic longitude/latitude of date.

| Flag | Effect |
|------|--------|
| `SEFLG_EQUATORIAL` | Right Ascension (0-360 deg) and Declination (+/-90 deg) instead of ecliptic coordinates. |
| `SEFLG_XYZ` | Cartesian (x, y, z) in AU instead of spherical coordinates. |
| `SEFLG_RADIANS` | Angles in radians instead of degrees. |

## Reference Frame

By default positions are precessed to the equinox of date.

| Flag | Effect |
|------|--------|
| `SEFLG_J2000` | J2000.0 reference frame (no precession to date). |
| `SEFLG_NONUT` | Mean ecliptic/equator (no nutation applied). |
| `SEFLG_ICRS` | ICRS frame (no precession, no nutation, no frame bias). |

## Position Corrections

By default positions are apparent (light-time + aberration corrected).

| Flag | Effect |
|------|--------|
| `SEFLG_TRUEPOS` | Geometric position: no light-time correction. |
| `SEFLG_NOABERR` | Astrometric: light-time corrected, no aberration. |
| `SEFLG_NOGDEFL` | Skip gravitational deflection correction. |
| `SEFLG_ASTROMETRIC` | Shorthand for `SEFLG_NOABERR \| SEFLG_NOGDEFL`. |

## Sidereal Zodiac

| Flag | Effect |
|------|--------|
| `SEFLG_SIDEREAL` | Subtract ayanamsha from ecliptic longitude. Requires prior `set_sid_mode()` call to select ayanamsha (Lahiri, Fagan-Bradley, etc.). |

## Compatibility

| Flag | Effect |
|------|--------|
| `SEFLG_MOSEPH` | Accepted for API compatibility, silently ignored. All calculations always use JPL DE440/DE441. |
| `SEFLG_SWIEPH` | Same -- accepted and ignored. |
| `SEFLG_SPEED3` | Converted to `SEFLG_SPEED` internally. |

## Examples

```python
import libephemeris as swe
from libephemeris.constants import *

jd = swe.julday(2024, 3, 26, 12.0)

# Default: geocentric ecliptic apparent position with speed
pos, _ = swe.calc_ut(jd, SE_MARS, SEFLG_SPEED)

# Heliocentric with speed
pos, _ = swe.calc_ut(jd, SE_MARS, SEFLG_SPEED | SEFLG_HELCTR)

# Sidereal equatorial
swe.set_sid_mode(SE_SIDM_LAHIRI)
pos, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_EQUATORIAL)

# J2000 ecliptic, no aberration
pos, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SPEED | SEFLG_J2000 | SEFLG_NOABERR)

# Topocentric (Rome)
swe.set_topo(12.4964, 41.9028, 0)
pos, _ = swe.calc_ut(jd, SE_MOON, SEFLG_SPEED | SEFLG_TOPOCTR)
```
