# Ayanamsha Definitions in LibEphemeris

## Overview

An **ayanamsha** is the angular difference between the tropical zodiac and the
sidereal zodiac. The tropical zodiac is defined by the vernal equinox (0° Aries),
while the sidereal zodiac is defined by fixed stars. Due to the precession of
the equinoxes, the two zodiacs slowly drift apart at approximately 50.3
arcseconds per year.

LibEphemeris supports all 43 ayanamsha modes defined in the Swiss Ephemeris
reference API, accessible via `swe_set_sid_mode()` and `swe_get_ayanamsa_ut()`.

## Ayanamsha Modes

### Traditional / Indian Systems

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 0 | `SE_SIDM_FAGAN_BRADLEY` | 24.7403 | Western sidereal (Fagan/Bradley) |
| 1 | `SE_SIDM_LAHIRI` | 23.857 | Official Indian government standard |
| 2 | `SE_SIDM_DELUCE` | 27.8158 | Robert DeLuce |
| 3 | `SE_SIDM_RAMAN` | 22.4108 | B.V. Raman |
| 4 | `SE_SIDM_USHASHASHI` | 20.0575 | Usha/Shashi |
| 5 | `SE_SIDM_KRISHNAMURTI` | 23.7602 | Krishnamurti Paddhati |
| 6 | `SE_SIDM_DJWHAL_KHUL` | 28.3597 | Djwhal Khul (Alice Bailey) |
| 7 | `SE_SIDM_YUKTESHWAR` | 22.4788 | Sri Yukteshwar |
| 8 | `SE_SIDM_JN_BHASIN` | 22.7621 | J.N. Bhasin |

### Babylonian Systems

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 9 | `SE_SIDM_BABYL_KUGLER1` | 23.5336 | Babylonian (Kugler 1) |
| 10 | `SE_SIDM_BABYL_KUGLER2` | 24.9336 | Babylonian (Kugler 2) |
| 11 | `SE_SIDM_BABYL_KUGLER3` | 25.7836 | Babylonian (Kugler 3) |
| 12 | `SE_SIDM_BABYL_HUBER` | 24.7336 | Babylonian (Huber) |
| 13 | `SE_SIDM_BABYL_ETPSC` | 24.5225 | Babylonian (Eta Piscium) |
| 14 | `SE_SIDM_ALDEBARAN_15TAU` | 24.76 | Aldebaran at 15° Taurus |
| 28 | `SE_SIDM_BABYL_BRITTON` | 24.6158 | Babylonian (Britton) |

### Historical / Classical

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 15 | `SE_SIDM_HIPPARCHOS` | 20.2478 | Hipparchos |
| 16 | `SE_SIDM_SASSANIAN` | 19.9930 | Sassanian |

### Reference Epoch Systems

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 17 | `SE_SIDM_GALCENT_0SAG` | ~26.85 | Galactic Center at 0° Sagittarius |
| 18 | `SE_SIDM_J2000` | 0.0000 | J2000.0 ecliptic (no precession applied) |
| 19 | `SE_SIDM_J1900` | 1.3966 | J1900.0 ecliptic |
| 20 | `SE_SIDM_B1950` | 0.6984 | B1950.0 ecliptic |

### Suryasiddhanta / Aryabhata

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 21 | `SE_SIDM_SURYASIDDHANTA` | 20.8951 | Suryasiddhanta |
| 22 | `SE_SIDM_SURYASIDDHANTA_MSUN` | 20.6804 | Suryasiddhanta (mean Sun) |
| 23 | `SE_SIDM_ARYABHATA` | 20.8951 | Aryabhata |
| 24 | `SE_SIDM_ARYABHATA_MSUN` | 20.6574 | Aryabhata (mean Sun) |
| 36 | `SE_SIDM_ARYABHATA_522` | 20.5758 | Aryabhata 522 CE |
| 25 | `SE_SIDM_SS_REVATI` | 20.1034 | SS Revati |
| 26 | `SE_SIDM_SS_CITRA` | 23.0058 | SS Citra |

### True Star-Based (Proper Motion Corrected)

These ayanamshas are defined by the true position of a reference star,
corrected for proper motion and precession in real time. Unlike fixed
ayanamshas, they track the actual stellar position.

| Mode | Constant | J2000 Value (°) | Reference Star |
|------|----------|-----------------|----------------|
| 27 | `SE_SIDM_TRUE_CITRA` | ~23.86 | Spica (α Virginis) at 0° Libra |
| 28 | `SE_SIDM_TRUE_REVATI` | ~20.05 | ζ Piscium at 0° Aries |
| 29 | `SE_SIDM_TRUE_PUSHYA` | ~22.72 | δ Cancri at 16° Cancer |
| 34 | `SE_SIDM_TRUE_MULA` | ~24.59 | λ Scorpii at 0° Sagittarius |
| 38 | `SE_SIDM_TRUE_SHEORAN` | ~25.23 | True Sheoran |
| 39 | `SE_SIDM_GALCENT_COCHRANE` | ~-3.15 | Galactic Center (Cochrane, 270° ecliptic) |
| 41 | `SE_SIDM_VALENS_MOON` | ~22.80 | Valens Moon |

### Galactic Center Based

These systems reference the Galactic Center (Sagittarius A*) position in
the ecliptic coordinate system.

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 17 | `SE_SIDM_GALCENT_0SAG` | ~26.85 | Galactic Center at 0° Sagittarius |
| 30 | `SE_SIDM_GALCENT_RGILBRAND` | ~22.47 | Galactic Center (Gil Brand) |
| 35 | `SE_SIDM_GALCENT_MULA_WILHELM` | ~20.04 | Galactic Center Mula (Wilhelm) |
| 39 | `SE_SIDM_GALCENT_COCHRANE` | ~-3.15 | Galactic Center (Cochrane) |

### Galactic Equator Based

| Mode | Constant | J2000 Value (°) | Description |
|------|----------|-----------------|-------------|
| 31 | `SE_SIDM_GALEQU_IAU1958` | ~30.11 | Galactic Equator (IAU 1958) |
| 32 | `SE_SIDM_GALEQU_TRUE` | ~30.11 | Galactic Equator (true) |
| 33 | `SE_SIDM_GALEQU_MULA` | ~23.51 | Galactic Equator Mula |
| 40 | `SE_SIDM_GALEQU_FIORENZA` | 25.0000 | Galactic Equator (Fiorenza) |
| 37 | `SE_SIDM_GALALIGN_MARDYKS` | ~30.11 | Galactic Alignment (Mardyks) |

### User-Defined

| Mode | Constant | Description |
|------|----------|-------------|
| 255 | `SE_SIDM_USER` | User-defined ayanamsha via `swe_set_sid_mode(SE_SIDM_USER, t0, ayan_t0)` |

The `SE_SIDM_USER` mode allows defining a custom ayanamsha by specifying:
- `t0`: Reference epoch as Julian Day
- `ayan_t0`: Ayanamsha value at the reference epoch

The library then applies standard precession to compute the ayanamsha at
any other date.

## Usage

```python
import libephemeris as leph
from libephemeris.constants import SE_SIDM_LAHIRI, SE_SIDM_USER

# Set a standard ayanamsha mode
leph.swe_set_sid_mode(SE_SIDM_LAHIRI)
ayan = leph.swe_get_ayanamsa_ut(2451545.0)  # J2000
print(f"Lahiri ayanamsha at J2000: {ayan:.4f}°")

# Use a custom ayanamsha
leph.swe_set_sid_mode(SE_SIDM_USER, 2451545.0, 23.5)
ayan = leph.swe_get_ayanamsa_ut(2460000.0)
print(f"Custom ayanamsha: {ayan:.4f}°")

# Compute sidereal positions
leph.swe_set_sid_mode(SE_SIDM_LAHIRI)
from libephemeris.constants import SE_SUN, SEFLG_SIDEREAL
pos, _ = leph.swe_calc_ut(2451545.0, SE_SUN, SEFLG_SIDEREAL)
print(f"Sun sidereal longitude: {pos[0]:.4f}°")
```

## Precession Model

All ayanamsha calculations use the IAU precession model via Skyfield,
matching the Swiss Ephemeris reference implementation. The precession
rate is approximately 50.29 arcseconds per year (varying slightly over
centuries due to the non-uniform precession model).

## Compatibility

All 43 ayanamsha modes produce results compatible with pyswisseph 2.10.
Star-based ("True") ayanamshas may show small differences (< 1°) due to
different proper motion catalogs and stellar position algorithms. See
`docs/divergences.md` for details.
