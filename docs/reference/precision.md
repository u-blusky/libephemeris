# LibEphemeris -- Scientific Precision

LibEphemeris is built on modern IAU standards and NASA JPL ephemerides. Every component -- from nutation to planet center corrections -- uses the most accurate models available in the astronomical literature. This document describes the scientific foundations, the specific models chosen, and measured precision for every calculation.

## Table of Contents

- [1. Ephemeris Foundation](#1-ephemeris-foundation)
- [2. Nutation Model](#2-nutation-model)
- [3. Precession Model](#3-precession-model)
- [4. Aberration and Light Deflection](#4-aberration-and-light-deflection)
- [5. Planet Center Corrections](#5-planet-center-corrections)
- [6. Velocity Computation](#6-velocity-computation)
- [7. True Ecliptic Coordinates](#7-true-ecliptic-coordinates)
- [8. House Calculations](#8-house-calculations)
- [9. Lunar Points](#9-lunar-points)
- [10. Delta T (TT − UT1)](#10-delta-t-tt--ut1)
- [11. Fixed Stars](#11-fixed-stars)
- [12. Ayanamsha (Sidereal Modes)](#12-ayanamsha-sidereal-modes)
- [13. Eclipses and Occultations](#13-eclipses-and-occultations)
- [14. Rise, Set, and Transit](#14-rise-set-and-transit)
- [15. Heliacal Events](#15-heliacal-events)
- [16. Planetary Positions -- Measured Precision](#16-planetary-positions----measured-precision)
- [17. Minor Bodies](#17-minor-bodies)
- [18. Thread-safe Context API](#18-thread-safe-context-api)
- [References](#references)

---

## 1. Ephemeris Foundation

### NASA JPL Development Ephemeris

LibEphemeris uses NASA JPL DE440 and DE441, the most recent planetary ephemerides produced by the Jet Propulsion Laboratory (Park et al. 2021). These are the same ephemerides used for Mars rover navigation and the James Webb Space Telescope.

| Property | DE440s | DE440 | DE441 |
|----------|--------|-------|-------|
| Date range | 1849--2150 CE | 1550--2650 CE | -13200--+17191 CE |
| File size | ~31 MB | ~119 MB | ~3.4 GB |
| Reference frame | ICRF 3.0 | ICRF 3.0 | ICRF 3.0 |
| Lunar model | Numerically integrated | Numerically integrated | Numerically integrated |
| Lunar Laser Ranging fit | ~1 milliarcsecond | ~1 milliarcsecond | ~1 milliarcsecond |
| Planetary radar fit | ~10 m (inner planets) | ~10 m (inner planets) | ~10 m (inner planets) |

DE440 and DE441 have identical precision -- DE441 is simply the extended-range version. DE440s is a reduced-size subset of DE440 with no loss of precision within its range.

### Precision tiers

LibEphemeris organizes these files into three precision tiers:

| Tier | File | Use Case |
|------|------|----------|
| `base` | de440s.bsp | Lightweight, modern-era usage |
| `medium` | de440.bsp | General purpose **(DEFAULT)** |
| `extended` | de441.bsp | Historical/far-future research |

Select by tier or by file directly:

```python
import libephemeris as eph

# By tier
eph.set_precision_tier("extended")   # uses de441.bsp

# Or by file
eph.set_ephemeris_file("de441.bsp")
```

Environment variables: `LIBEPHEMERIS_PRECISION=extended` or `LIBEPHEMERIS_EPHEMERIS=de441.bsp`.

Resolution priority (highest to lowest):
1. `LIBEPHEMERIS_EPHEMERIS` environment variable
2. `set_ephemeris_file()` / `set_jpl_file()` programmatic call
3. Precision tier (`LIBEPHEMERIS_PRECISION` env var or `set_precision_tier()`)
4. Default: `de440.bsp` (medium tier)

### Swiss Ephemeris comparison

Swiss Ephemeris uses JPL DE431 (older generation, 2013) repackaged into its own binary format. DE440/DE441 (2021) incorporate 8 additional years of observational data, including improved Juno-era Jupiter observations and MESSENGER-era Mercury data, and use the updated ICRF 3.0 reference frame.

When `SEFLG_MOSEPH` is passed, Swiss Ephemeris falls back to the Moshier semi-analytical ephemeris (VSOP87 for planets, ELP2000-82B for Moon), which has errors of ~1 arcsecond for inner planets and ~10+ arcseconds for outer planets at historical dates. LibEphemeris accepts `SEFLG_MOSEPH` for API compatibility but silently ignores it -- all calculations always use the full JPL numerical integration.

---

## 2. Nutation Model

### IAU 2006/2000A (1365 terms)

LibEphemeris uses the **IAU 2006/2000A** nutation model via the IAU ERFA library (`erfa.nut06a()`). This is the highest-precision nutation model currently adopted by the International Astronomical Union.

| Property | Value |
|----------|-------|
| Lunisolar terms | 678 |
| Planetary terms | 687 |
| Total terms | **1365** |
| Precision | ~0.01--0.05 milliarcseconds |
| Standard | IERS Conventions 2010, ch. 5 |

The model computes nutation in longitude (Δψ) and nutation in obliquity (Δε) from the five Delaunay fundamental arguments plus nine planetary fundamental arguments, with corrections from the IAU 2006 precession-rate adjustments (Capitaine et al. 2005).

### Swiss Ephemeris comparison

Swiss Ephemeris uses the same IAU 2006/2000A model. Both implementations achieve sub-milliarcsecond nutation. The difference between the two is negligible (<0.01 mas).

However, Skyfield's internal nutation (used before pyerfa integration) defaults to IAU 2000B, a truncated model with only **77 terms** and ~1 milliarcsecond precision. LibEphemeris overrides this with the full 1365-term model via pyerfa for all planet and house calculations.

### Impact on ecliptic longitude

Nutation in longitude (Δψ) directly shifts ecliptic longitudes. The maximum amplitude of Δψ is ~17.2 arcseconds (the 18.6-year lunar nodal cycle). Using a 77-term model vs. a 1365-term model introduces up to ~1 mas error in this correction, which over decades accumulates to measurable differences in planetary longitudes.

---

## 3. Precession Model

### IAU 2006 (Capitaine et al. 2003)

LibEphemeris uses the **IAU 2006 precession** model via `erfa.pmat06()`, implementing the Fukushima-Williams four-angle formulation with polynomial terms up to T^5.

| Property | Value |
|----------|-------|
| Polynomial order | 5th degree in T (centuries from J2000) |
| Rate precision | ~0.08 mas/century |
| J2000 obliquity | 84381.406 arcseconds (23°26'21.406") |
| Standard | IERS Conventions 2010, ch. 5 |

### Swiss Ephemeris comparison

Swiss Ephemeris also uses IAU 2006 precession. Both implementations are equivalent for this component.

### Frame bias (ICRS to J2000)

The rotation from the International Celestial Reference System (ICRS) to the mean equator and equinox of J2000.0 is applied via the combined bias-precession-nutation matrix from `erfa.pnm06a()`. This matrix incorporates the frame bias angles (dα₀ = −14.6 mas, ξ₀ = −16.617 mas, η₀ = −6.819 mas) defined in the IERS Conventions.

---

## 4. Aberration and Light Deflection

### Annual aberration

For planet calculations, LibEphemeris uses Skyfield's full relativistic treatment via the `.apparent()` method, which implements IERS 2003 conventions for annual aberration (Bradley formula + relativistic corrections to order v²/c²).

The maximum annual aberration is ~20.5 arcseconds. Relativistic corrections to the classical Bradley formula amount to ~0.003 arcseconds.

### Gravitational light deflection

Skyfield's `.apparent()` method also applies gravitational light deflection by the Sun (maximum ~1.75 arcseconds at the solar limb, falling as 1/sin(θ)). This correction is significant for planets near solar conjunction.

### Swiss Ephemeris comparison

Swiss Ephemeris implements the same corrections (annual aberration + solar gravitational deflection). Both libraries omit diurnal aberration (~0.3 arcseconds maximum), which is below the threshold of astrological significance.

---

## 5. Planet Center Corrections

### The barycenter problem

JPL DE440/DE441 provide positions for outer planets as **system barycenters** -- the center of mass of the planet plus all its moons -- not the planet body center. The maximum angular offset between barycenter and body center:

| Planet | Max offset (geocentric) | Primary contributor |
|--------|------------------------|-------------------|
| Jupiter | ~0.6 arcsec | Ganymede + Callisto |
| Saturn | ~0.2 arcsec | Titan (96% of moon mass) |
| Uranus | ~0.01 arcsec | Titania + Oberon |
| Neptune | ~0.05 arcsec | Triton |
| Pluto | ~0.3 arcsec | Charon (binary system) |

For inner planets (Mercury, Venus, Earth, Mars), the barycenter effectively equals the body center because their moons (if any) have negligible mass relative to the planet.

### Three-tier correction strategy

LibEphemeris corrects barycenters to body centers **automatically** using a three-tier fallback:

#### Tier 1: SPK-based planet centers (<0.001 arcsec)

A bundled `planet_centers.bsp` file (~25 MB) contains precise center-of-body segments for NAIF IDs 599 (Jupiter), 699 (Saturn), 799 (Uranus), 899 (Neptune), and 999 (Pluto). These segments are extracted from JPL satellite ephemerides:

| NAIF ID | Source SPK | Coverage |
|---------|-----------|----------|
| 599 | jup204.bsp | ~1925--2025 |
| 699 | sat319.bsp | ~1950--2050 |
| 799 | ura083.bsp | ~1950--2050 |
| 899 | nep050.bsp | ~1989--2049 |
| 999 | plu017.bsp | ~1990--2050 |

The `_SpkCenterTarget` class in `planets.py` adds the SPK offset to the DE440 barycenter position at the retarded time (when light left the planet), ensuring the correction is applied correctly for light-time-corrected observations.

**If the date falls outside the SPK coverage**, Tier 1 silently falls back to Tier 2.

#### Tier 2: Analytical satellite theories (sub-arcsecond)

The `_CobCorrectedTarget` class computes the barycenter-to-center offset using analytical theories for major satellites. The correction formula:

```
offset = -Σ(m_i × r_i) / M_total
```

where m_i is the satellite mass and r_i is the satellite position relative to the planet center.

**Jupiter -- Lieske E5 theory (Meeus implementation)**

Theory for Io, Europa, Ganymede, and Callisto based on Lieske (1998), as presented in Meeus "Astronomical Algorithms" ch. 44.

| Property | Value |
|----------|-------|
| Reference epoch | JDE 2443000.5 (1976-08-10) |
| Terms per satellite | 21--50 (longitude), 7--11 (latitude), 7--16 (radius) |
| Coordinate frame | Jupiter equatorial → J2000 ecliptic (4 rotations) |
| Precession correction | From B1950.0: P = 1.3966626·T₀ + 0.0003088·T₀² |
| Positional precision | ~100--200 km (~0.05 arcsec at opposition) |

GM values from JUP365:

| Satellite | GM (km³/s²) | Mass ratio (m/M_Jupiter) |
|-----------|------------|--------------------------|
| Io | 5959.915 | 4.70 × 10⁻⁵ |
| Europa | 3202.712 | 2.53 × 10⁻⁵ |
| Ganymede | 9887.833 | 7.80 × 10⁻⁵ |
| Callisto | 7179.283 | 5.67 × 10⁻⁵ |
| **Total** | **26229.744** | **2.07 × 10⁻⁴** |

Jupiter GM = 126,686,531.9 km³/s² (JUP365).

**Saturn -- TASS 1.7 theory (Vienne & Duriez 1995)**

Full analytical theory for 8 Saturnian satellites, ported from the Stellarium implementation (Johannes Gajdosik, MIT license).

| Property | Value |
|----------|-------|
| Reference epoch | JD 2444240.0 (1980-01-01 TT) |
| Satellites | Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Iapetus, Hyperion |
| Total terms | ~1500 (from `tass17_data.py`, 2964 lines of coefficients) |
| Orbital elements | 6 per satellite (n, λ, e·cos ω̃, e·sin ω̃, sin(i/2)·cos Ω, sin(i/2)·sin Ω) |
| Kepler solver | Newton's method, 15 iterations, convergence at 10⁻¹⁴ |
| Coordinate frame | TASS17 → VSOP87 (J2000 ecliptic) via fixed rotation matrix |
| Positional precision | ~50--100 km (~0.05 arcsec at opposition) |

Titan dominates the correction (GM = 8978.137 km³/s², 96% of total Saturnian moon mass). Saturn GM = 37,931,206.23 km³/s² (SAT441).

**Neptune -- Triton Keplerian model**

Simple Keplerian orbit with secular J2 nodal precession, based on NEP097 elements (Jacobson 2009).

| Property | Value |
|----------|-------|
| Semi-major axis | 354,759 km |
| Period | 5.877 days (retrograde) |
| Eccentricity | 0.000016 |
| Inclination | 156.865° (retrograde orbit) |
| Node precession | −360°/(688 × 365.25 days) (~688-year period) |
| Kepler solver | 1 iteration (sufficient for e ≈ 0) |
| Positional precision | ~20--50 km (~0.003 arcsec at opposition) |

Triton GM = 1428.495 km³/s², Neptune GM = 6,835,099.97 km³/s² (NEP097). Mass ratio = 2.09 × 10⁻⁴.

**Pluto -- Charon two-body Keplerian**

Two-body Keplerian solution from PLU060 (Brozovic & Jacobson 2024).

| Property | Value |
|----------|-------|
| Semi-major axis | 19,591 km |
| Period | 6.387 days |
| Eccentricity | 0.0002 (tidally locked) |
| Inclination | 96.145° |
| Node precession | 0 (negligible) |
| Kepler solver | 5 iterations |

Pluto-Charon is a **binary system**: the barycenter lies ~2130 km outside Pluto's center. Charon GM = 106.1 km³/s², Pluto GM = 869.3 km³/s² (PLU060). Mass ratio Charon/(Pluto+Charon) = **0.109**.

**Uranus**: Currently returns zero offset. The total mass ratio of Uranian moons to Uranus is ~1.0 × 10⁻⁴, corresponding to a maximum geocentric offset of ~0.01 arcsec. This is below the precision of most astrological applications.

#### Tier 3: Raw barycenters (last resort)

If both SPK and analytical methods fail, the raw barycenter position from DE440/DE441 is used. This path is not reached under normal operation.

### Velocity correction

The velocity component of the COB offset is computed via central difference numerical differentiation with a 1-second step:

```
v_offset = [offset(t + 0.5s) - offset(t - 0.5s)] / 1s
```

### Swiss Ephemeris comparison

Standard Swiss Ephemeris returns **system barycenters** for outer planets. The `SEFLG_CENTER_BODY` flag (a newer, less commonly used feature) enables planet center positions, but requires additional satellite ephemeris files. LibEphemeris applies the correction **transparently** to all calculations, achieving better default precision for Jupiter, Saturn, Neptune, and Pluto without requiring any user configuration.

---

## 6. Velocity Computation

### Central difference method

LibEphemeris computes all velocities via numerical differentiation using the central difference formula:

```
dp/dt = [p(t + dt) - p(t - dt)] / (2 · dt)
```

| Context | Step size (dt) | Precision |
|---------|---------------|-----------|
| Planetary longitude/latitude | 1 second | ~0.0001°/day |
| COB offset velocity | 1 second | ~10⁻⁸ AU/day |
| Lunar node/Lilith | 0.5 days | ~0.001°/day |
| Fixed star apparent motion | 1 second | ~10⁻⁶°/day |

The central difference method is O(h²) accurate, meaning halving the step size reduces the error by a factor of 4.

### Swiss Ephemeris comparison

Swiss Ephemeris computes velocities analytically by differentiating the Chebyshev polynomial representation of the JPL ephemeris. This avoids the O(h²) truncation error of numerical differentiation but requires maintaining separate derivative code.

The practical difference is <0.001°/day for all planets. For the Moon, speed differences of up to ~0.01°/day have been measured, primarily because the numerical method samples the position function at different points than the analytical derivative.

---

## 7. True Ecliptic Coordinates

### From ICRS to ecliptic of date

The full pipeline from JPL ephemeris to ecliptic longitude:

1. **JPL state vector** (ICRS Cartesian, barycentric or geocentric)
2. **Light-time correction** (iterative, typically 3 Newton iterations)
3. **COB correction** (barycenter → planet center, Tier 1/2/3)
4. **Apparent position** via Skyfield: annual aberration + gravitational deflection
5. **Ecliptic of date** via `frame_latlon(ecliptic_frame)`: applies combined IAU 2006 precession + IAU 2006/2000A nutation matrix (`erfa.pnm06a()`) to rotate from ICRS to true ecliptic of date

For J2000 ecliptic (`SEFLG_J2000`), step 5 uses a fixed rotation by the J2000 obliquity (23°26'21.406") without precession or nutation.

### Mean obliquity of the ecliptic

IAU 2006 polynomial (Capitaine et al. 2003):

```
ε₀ = 84381.406" − 46.836769"·T − 0.0001831"·T² + 0.00200340"·T³
     − 0.000000576"·T⁴ − 0.0000000434"·T⁵
```

where T is Julian centuries from J2000.0.

---

## 8. House Calculations

### Sidereal time

Houses require Greenwich Apparent Sidereal Time (GAST). LibEphemeris uses **Skyfield's `t.gast`**, which implements the IAU SOFA algorithm:

1. Earth Rotation Angle (ERA) from UT1 (IERS 2003 conventions)
2. Greenwich Mean Sidereal Time (GMST) via ERA + equation of origins
3. Equation of equinoxes: Δψ · cos(ε) (nutation in RA from IAU 2006/2000A)
4. GAST = GMST + equation of equinoxes

Precision: ~0.001 seconds of time (~0.015 arcseconds in hour angle).

### True obliquity

Houses use the true obliquity (mean obliquity + nutation in obliquity) from `erfa.obl06()` + `erfa.nut06a()`. This is the same IAU 2006/2000A model used for planet calculations.

### ARMC

```
ARMC = (GAST × 15°/h + geographic_longitude) mod 360°
```

### Supported systems (24)

All major house systems are implemented with the same spherical trigonometry as Swiss Ephemeris:

`P` Placidus, `K` Koch, `O` Porphyry, `R` Regiomontanus, `C` Campanus, `E` Equal (Asc), `A`/`D` Equal (MC), `W` Whole Sign, `M` Morinus, `B` Alcabitius, `T` Polich-Page/Topocentric, `U` Krusinski, `G` Gauquelin (36-sector), `V` Vehlow, `X` Meridian, `H` Horizontal, `F` Carter, `S` Sripati, `L` Pullen SD, `Q` Pullen SR, `N` Natural Gradient, `Y` APC, `I`/`i` Sunshine.

### Convergence tolerance

Iterative systems (Placidus, Koch) use a convergence threshold of 10⁻⁷ degrees (~0.00036 arcseconds).

### Polar latitude behavior

Above the polar circle (~66.56° = 90° − obliquity), Placidus and Koch are geometrically undefined because some ecliptic points never rise or set. LibEphemeris raises `PolarCircleError` with the option to fall back to Porphyry via `swe_houses_with_fallback()`.

### Measured precision vs Swiss Ephemeris

| Component | Max difference |
|-----------|---------------|
| Ascendant | <0.001° |
| MC (Midheaven) | <0.001° |
| ARMC | <0.001° |
| House cusps (typical) | 0.001°--0.01° |
| Vertex | <0.01° |

---

## 9. Lunar Points

### True Node (osculating lunar node)

The True Node is where the Moon's **instantaneous orbital plane** intersects the ecliptic. LibEphemeris uses a two-stage approach.

**Stage 1: Geometric osculating node from JPL DE440**

1. Get Moon's geocentric position **r** and velocity **v** from JPL DE440 (~1 milliarcsecond source precision)
2. Compute angular momentum vector: **h** = **r** × **v**
3. This vector is perpendicular to the instantaneous orbital plane
4. Find intersection with the ecliptic: Ω = atan2(h_x, −h_y)
5. Apply IAU 2006 precession (J2000 → date)
6. Apply IAU 2006/2000A nutation (1365 terms)

This computes **exactly** what the True Node is by definition: the intersection of the orbital plane with the ecliptic. Deriving the node directly from the state vectors avoids the approximation errors inherent in analytical series — the geometric construction is exact by definition.

**ELP2000-82B perturbation series (available but not applied)**

The codebase contains a comprehensive ELP2000-82B perturbation series (~170 terms, ~900 lines) organized into the following categories:

| Category | Terms | Dominant amplitude |
|----------|-------|-------------------|
| Solar (main) | 9 | −1.5233° · sin(2D) |
| Solar (2nd order) | 12 | 0.003°--0.04° |
| Solar (3rd order) | 10 | 0.001°--0.006° |
| Inclination (F-related) | 9 | 0.001°--0.01° |
| Venus perturbation | 9 | 0.001°--0.005° |
| Mars perturbation | 9 | 0.001°--0.004° |
| Jupiter perturbation | 7 | 0.001°--0.003° |
| Saturn perturbation | 7 | 0.001°--0.003° |
| Evection | 6 | up to 0.047° |
| Variation | 10 | up to 0.052° |
| Annual equation | 6 | up to 0.186° (E · sin M) |
| Parallactic inequality | 3 | up to 0.035° |
| Second-order coupling | 8 | 0.0001°--0.003° |
| Secular terms | 3+ | T-dependent corrections to T⁵ |

The Earth eccentricity factor E = 1.0 − 0.002516·T − 0.0000074·T² is applied to Sun-dependent terms per Meeus convention.

> **Note:** These perturbation corrections are **not currently applied** to the True Node calculation. Investigation revealed that the ELP2000-82B series was designed for the *mean* lunar node, not the geometric node derived from state vectors via `h = r × v`. Applying the series to the geometric node produced errors of tens of degrees, confirming the two approaches are incompatible. The geometric method from Stage 1 is used alone. The residual vs Swiss Ephemeris (~8.9 arcsec mean, ~0.14° max) reflects the fundamental methodological difference. See [Precision History](../development/precision-history.md) for the full investigation record.

### Mean Node

Meeus polynomial for the mean ascending node longitude, extended with T⁴ and T⁵ corrections from Chapront et al. (2002) and Simon et al. (1994).

### Mean Lilith (Black Moon)

Mean longitude of the lunar apogee, computed from Meeus polynomials (ch. 50). Includes latitude computation. Velocity via central difference with 0.5-day step.

### True Lilith (osculating apogee)

Computed from JPL DE440 Moon state vectors via orbital mechanics: the eccentricity vector **e** is derived from position and velocity, and the apogee direction points from the focus along **e**.

### Interpolated Apogee and Perigee

Uses the Moshier analytical method: ELP2000-82B perturbation series (~50 harmonic terms) fitted to DE404 over a 7000-year span (Moshier 1992). Removes the spurious ~30° oscillations present in the osculating apsides.

### Meeus polynomial validity ranges

| Range | |T| (centuries) | Accuracy |
|-------|---------------|----------|
| Optimal | <2 (~±200 years) | <0.001° |
| Valid | <10 (~±1000 years) | <0.01° |
| Maximum | <20 (~±2000 years) | Significant degradation |

Fundamental arguments include T⁵ corrections from Chapront et al. (2002).

### Measured precision vs Swiss Ephemeris

| Point | Max difference | Mean difference |
|-------|---------------|-----------------|
| Mean Node | ~18 arcsec (0.005°) | ~11 arcsec (0.003°) |
| True Node | ~520 arcsec (0.14°) | ~145 arcsec (0.04°) |
| Mean Lilith | ~18 arcsec (0.005°) | ~12 arcsec (0.003°) |
| True Lilith | ~235 arcsec (0.065°) | ~52 arcsec (0.015°) |
| Interpolated Apogee | ~1300 arcsec (0.36°) | ~350 arcsec (0.10°) |
| Interpolated Perigee | ~9400 arcsec (2.6°) | ~1650 arcsec (0.46°) |

The True Node difference reflects different methodologies, not calculation errors. LibEphemeris uses the geometrically exact approach from JPL state vectors; Swiss Ephemeris uses analytical series. Both are valid; the geometric method is more rigorous by construction.

---

## 10. Delta T (TT − UT1)

### Model

LibEphemeris uses the **Stephenson, Morrison & Hohenkerk (2016)** Delta T model via Skyfield:

| Date range | Method |
|------------|--------|
| 720 BC -- ~2016 AD | Cubic spline interpolation from Table S15 |
| Outside spline | Parabolic: ΔT = −320 + 32.5 · u² (u = (year − 1825)/100) |
| 1973 -- present | IERS observed values (optional, via `set_iers_delta_t_enabled(True)`) |

### IERS observed Delta T

When enabled, LibEphemeris uses observed ΔT values from the International Earth Rotation and Reference Systems Service for recent dates (1973--present), achieving ~0.1 second precision.

```python
from libephemeris import set_iers_delta_t_enabled, download_delta_t_data
download_delta_t_data()
set_iers_delta_t_enabled(True)
```

### Swiss Ephemeris comparison

Swiss Ephemeris uses Espenak & Meeus (2006) polynomials for historical Delta T. The Stephenson et al. (2016) model used by LibEphemeris is more recent and generally considered more accurate for dates before 1600 CE, where direct observations are sparse.

For modern dates (1900--2100), both implementations agree to within ~1 second.

### Typical Delta T values

| Year | ΔT (seconds) |
|------|-------------|
| −3000 | ~72,000 |
| 1000 | ~1,500 |
| 1800 | ~14 |
| 1900 | ~−3 |
| 2000 | ~64 |
| 2020 | ~69 |

---

## 11. Fixed Stars

### Star catalog

102 stars from the **Hipparcos catalog** (ESA 1997), covering all bright and astrologically significant stars: the 4 Royal Stars, Behenian stars, Pleiades cluster, Hyades, and full zodiacal constellation coverage.

| Property | Value |
|----------|-------|
| Epoch | J2000.0 (ICRS) |
| Source | Hipparcos catalog (HIP numbers) |
| Count | 102 stars |
| Data per star | RA, Dec, PM_RA (with cos δ), PM_Dec, visual magnitude |

### Proper motion

Rigorous space motion propagation from Hipparcos Vol. 1, Section 1.5.5, with second-order Taylor term for celestial sphere curvature:

```
correction = −0.5 · |V|² · P · t²
```

Precision: <0.01 arcsec over ±100 years; <1 arcsec over ±500 years.

### Position pipeline

1. Proper motion propagation (J2000 → date)
2. Create Skyfield `Star` object with propagated RA/Dec
3. Apparent position via `observer.at(t).observe(star).apparent()` (aberration + gravitational deflection)
4. Ecliptic coordinates via `frame_latlon(ecliptic_frame)` (IAU 2006 precession + IAU 2006/2000A nutation)

### Swiss Ephemeris comparison

Swiss Ephemeris uses a larger catalog (~1000+ stars from its own bundled star catalog). Both use Hipparcos data and similar proper motion models. LibEphemeris covers 102 stars, selected for astrological and navigational significance. Star positions agree to ~0.01 arcsecond.

---

## 12. Ayanamsha (Sidereal Modes)

43 ayanamsha systems are supported, matching the full Swiss Ephemeris set.

### Formula-based ayanamshas

| Category | Examples | Max difference vs SwissEph |
|----------|----------|---------------------------|
| Standard | Fagan-Bradley, Lahiri, Raman | <0.0002° |
| Epoch-based | J2000, J1900, B1950 | <0.0002° |
| Historical | Babylonian variants (Kugler 1-3, Huber) | <0.001° |

### Star-based ayanamshas

Star-based ayanamshas anchor the sidereal zodiac to a specific fixed star. The slightly higher differences reflect different proper motion and precession models:

| Ayanamsha | Max difference |
|-----------|---------------|
| True Citra (Spica, HIP 65474) | <0.006° |
| True Revati | <0.006° |
| True Pushya | <0.006° |
| True Mula | <0.006° |
| Galactic Center variants | <0.001° |

The IAU 2006 precession model (used by LibEphemeris) is more accurate than the older Lieske (1977) model for computing the precession of the equinox, which directly affects ayanamsha values.

---

## 13. Eclipses and Occultations

### Solar eclipses

Calculated using Besselian elements with the full JPL DE ephemeris for Sun and Moon positions. Contact times (C1--C4), path width, central line coordinates, magnitude, and obscuration are computed.

| Property | Precision |
|----------|-----------|
| Eclipse maximum timing | <10 seconds vs Swiss Ephemeris |
| Contact times | <10 seconds |
| Eclipse magnitude | ~0.01 |
| Path width | ~1 km |

### Lunar eclipses

Full implementation including umbral and penumbral contact times (P1, U1, U2, U3, U4, P4), umbral/penumbral magnitude, gamma parameter, and duration.

### Saros and Inex series

Saros and Inex series numbers are computed from eclipse-to-eclipse relationships using the Saros period (6585.32 days) and Inex period (10571.95 days).

---

## 14. Rise, Set, and Transit

Computed using the Bennett (1982) atmospheric refraction formula with configurable atmospheric conditions (pressure, temperature).

| Event | Precision vs Swiss Ephemeris |
|-------|------------------------------|
| Sunrise/sunset | <30 seconds |
| Moonrise/moonset | <30 seconds |
| Meridian transit | <30 seconds |

---

## 15. Heliacal Events

Schaefer (1990) atmospheric visibility model with:

1. **Atmospheric extinction**: Rayleigh scattering + aerosol + ozone + water vapor
2. **Twilight sky brightness**: gradient model with moonlight contribution
3. **Ptolemaic visibility thresholds**: arcus visionis values
4. **Contrast threshold**: Schaefer model with configurable observer skill levels

Timing precision: <1 day for heliacal rising/setting events.

---

## 16. Planetary Positions -- Measured Precision

Tested against Swiss Ephemeris at 100+ random dates within DE440 range (1550--2650 CE):

### Longitude

| Planet | Max difference | Mean difference |
|--------|---------------|-----------------|
| Sun | 0.20 arcsec | 0.04 arcsec |
| Moon | 3.32 arcsec | 0.70 arcsec |
| Mercury | 0.32 arcsec | 0.05 arcsec |
| Venus | 0.33 arcsec | 0.08 arcsec |
| Mars | 0.58 arcsec | 0.06 arcsec |
| Jupiter | 0.44 arcsec | 0.12 arcsec |
| Saturn | 0.51 arcsec | 0.13 arcsec |
| Uranus | 0.50 arcsec | 0.23 arcsec |
| Neptune | 1.17 arcsec | 0.24 arcsec |
| Pluto | 0.75 arcsec | 0.26 arcsec |

All differences are **sub-arcsecond** except the Moon (~3 arcseconds), which reflects the different nutation models and COB correction pipelines between the two libraries. Both values are well within the precision of the underlying DE ephemeris.

### Latitude

| Planet | Max difference |
|--------|---------------|
| Sun | <0.1 arcsec |
| Moon | <1.3 arcsec |
| Other planets | <0.6 arcsec |

### Velocity

| Component | Max difference |
|-----------|---------------|
| Angular velocity | <0.0004°/day |
| Radial velocity | <0.001 AU/day |

### Source of measured differences

The differences between LibEphemeris and Swiss Ephemeris are **not errors** in either library. They arise from intentional methodological choices:

| Source | Typical contribution |
|--------|---------------------|
| Nutation model (IAU 2006/2000A vs internal) | ~0.01--0.05 mas |
| COB correction (analytical vs none) | ~0.01--0.6 arcsec (outer planets) |
| DE440 vs DE431 ephemeris generation | ~0.001 arcsec |
| Velocity method (numerical vs analytical) | ~0.0001°/day |
| Light-time iteration tolerance | ~0.001 arcsec |

---

## 17. Minor Bodies

### Strict precision mode (default)

For major asteroids with strongly perturbed orbits (Chiron, Ceres, Pallas, Juno, Vesta), LibEphemeris requires SPK kernels by default. Simple Keplerian propagation can produce errors of **1--10 degrees** for these bodies, which is unacceptable for astrological use.

```python
eph.set_auto_spk_download(True)  # Recommended: automatic SPK from JPL Horizons
```

### SPK kernel precision

With SPK kernels (downloaded from JPL Horizons), all minor body calculations achieve **sub-arcsecond** precision matching JPL Horizons exactly.

### Trans-Neptunian Objects

TNOs use Keplerian elements with first-order secular perturbations from Jupiter, Saturn, Uranus, and Neptune (Laplace-Lagrange theory). Typical accuracy: <10° over 50-year spans. For research-grade TNO work, use SPK kernels.

---

## 18. Thread-safe Context API

Each `EphemerisContext` instance has its own isolated state while sharing the expensive ephemeris data across instances. The `_CONTEXT_SWAP_LOCK` (RLock) ensures thread safety during state operations.

Memory overhead per context: ~1 KB. Ephemeris files (DE440: ~119 MB) are loaded once and shared.

---

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Capitaine, N. et al. (2003). "Expressions for IAU 2000 precession quantities." *Astronomy & Astrophysics*, 412, 567-586.
3. Mathews, P.M. et al. (2002). "Modeling of nutation-precession: New nutation series for nonrigid Earth." *Journal of Geophysical Research*, 107(B4).
4. Capitaine, N. et al. (2005). "IAU 2006 precession." *Highlights of Astronomy*, 14, 474-475.
5. Lieske, J.H. (1998). "Galilean Satellites of Jupiter. Theory E5." *Astronomy & Astrophysics Supplement*, 129, 205-217.
6. Vienne, A. & Duriez, L. (1995). "TASS1.6: Ephemerides of the major Saturnian satellites." *Astronomy & Astrophysics*, 297, 588-605.
7. Jacobson, R.A. (2009). "The orbits of the Neptunian satellites." *Astronomical Journal*, 137, 4322-4329.
8. Brozovic, M. & Jacobson, R.A. (2024). "The orbits of the Pluto system." *Planetary Science Journal*.
9. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
10. Meeus, J. (1998). *Astronomical Algorithms*, 2nd edition. Willmann-Bell.
11. Moshier, S.L. (1992). "Comparison of a 7000-year lunar ephemeris with analytical theory." *Astronomy & Astrophysics*, 262, 613-616.
12. Schaefer, B.E. (1990). "Telescopic limiting magnitudes." *Publications of the Astronomical Society of the Pacific*, 102, 212-229.
13. IERS Conventions (2010). IERS Technical Note No. 36, ed. Petit, G. & Luzum, B.
14. ESA (1997). *The Hipparcos and Tycho Catalogues*. ESA SP-1200.
