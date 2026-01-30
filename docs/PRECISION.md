# LibEphemeris Precision Limitations

This document provides detailed information about the precision and limitations of LibEphemeris compared to pyswisseph (Swiss Ephemeris). Use this guide to determine if LibEphemeris is suitable for your use case.

## Summary Table

| Component | Precision vs pyswisseph | Notes |
|-----------|------------------------|-------|
| Sun position | ±0.2 arcsec | High precision (max 0.2", mean 0.04") |
| Mercury, Venus positions | ±0.4 arcsec | High precision |
| Mars position | ±0.6 arcsec | High precision |
| Moon position | ±3.5 arcsec | Relaxed due to lunar theory differences |
| Jupiter, Saturn positions | ±0.6 arcsec | Outer planets |
| Uranus, Neptune, Pluto | ±1.2 arcsec | High precision with DE440 |
| Mean Lunar Node | ~0.005° (~18 arcsec) | High precision formula |
| True Lunar Node | ~0.14° (~520 arcsec) | Osculating orbital elements method |
| Mean Lilith | ~0.12° (~430 arcsec) | Minor formula differences |
| House cusps | ±0.001° (~3.6 arcsec) | All 19 house systems |
| Ayanamsha (standard) | ±0.0002° | Fagan-Bradley, Lahiri, Raman (high precision) |
| Ayanamsha (star-based) | ±0.006-0.06° | True Citra (0.006°), others (0.06°) |
| Sun crossings (ingress) | ±0.001 arcsec | Sub-milliarcsecond precision |
| Moon crossings | ±0.05 arcsec | Sub-arcsecond precision |
| Planet crossings | ±0.1 arcsec | Sub-arcsecond precision |
| Polar houses (>66.5°) | Falls back to Porphyry | Placidus/Koch undefined |
| Planetary moons | ±1 arcsec | Requires SPK kernel registration |

---

## Planetary Positions

### Longitude Precision

Tested against pyswisseph at 100+ random dates within DE440 range (1550-2650):

| Planet | Max Difference | Mean Difference | Tolerance |
|--------|---------------|-----------------|-----------|
| Sun | 0.20 arcsec | 0.04 arcsec | ±0.3 arcsec |
| Moon | 3.32 arcsec | 0.70 arcsec | ±3.5 arcsec |
| Mercury | 0.32 arcsec | 0.05 arcsec | ±0.4 arcsec |
| Venus | 0.33 arcsec | 0.08 arcsec | ±0.4 arcsec |
| Mars | 0.58 arcsec | 0.06 arcsec | ±0.6 arcsec |
| Jupiter | 0.44 arcsec | 0.12 arcsec | ±0.5 arcsec |
| Saturn | 0.51 arcsec | 0.13 arcsec | ±0.6 arcsec |
| Uranus | 0.50 arcsec | 0.23 arcsec | ±0.6 arcsec |
| Neptune | 1.17 arcsec | 0.24 arcsec | ±1.2 arcsec |
| Pluto | 0.75 arcsec | 0.26 arcsec | ±0.8 arcsec |

**Note on Moon**: The Moon has a relaxed tolerance due to differences in lunar theory implementations between LibEphemeris (Skyfield/JPL DE440) and pyswisseph. Both are accurate, but use slightly different models.

**Note on DE440**: LibEphemeris now uses DE440 as the default ephemeris, which provides improved accuracy for outer planets with 12+ years more observational data than DE421, and uses the ICRF 3.0 reference frame.

### Latitude Precision

| Planet | Max Difference |
|--------|---------------|
| Sun | <0.1 arcsec |
| Moon | <1.3 arcsec |
| Other planets | <0.6 arcsec |

### Distance Precision

| Component | Max Difference |
|-----------|---------------|
| Distance (AU) | <0.0001 AU |

### Velocity Precision

| Component | Max Difference |
|-----------|---------------|
| Angular velocity | <0.0004°/day |
| Radial velocity | <0.001 AU/day |

---

## Coordinate Systems

### Geocentric vs Heliocentric

| Coordinate System | Longitude Tolerance |
|-------------------|---------------------|
| Geocentric | ±0.001° |
| Topocentric | ±0.001° |
| Heliocentric | ±0.03° (slightly relaxed) |
| Barycentric | ±0.03° (slightly relaxed) |

---

## House Systems

### Cusp Precision

LibEphemeris supports all 19 Swiss Ephemeris house systems:

| System | Precision | Notes |
|--------|-----------|-------|
| Placidus (P) | ~0.001-0.01° | Iterative, fails at polar latitudes |
| Koch (K) | ~0.001-0.01° | Similar to Placidus limitations |
| Porphyry (O) | ~0.01-0.1 arcsec | Geometric, works everywhere |
| Regiomontanus (R) | ~0.01° | Complex calculation |
| Campanus (C) | ~0.01° | Complex calculation |
| Equal (A/E) | Exact | No astronomical calculations |
| Whole Sign (W) | Exact | No astronomical calculations |
| Morinus (M) | ~0.01-0.1 arcsec | Geometric |
| Meridian (X) | ~0.01-0.1 arcsec | Geometric |
| Polich-Page/Topocentric (T) | ~0.01° | Complex calculation |
| Alcabitus (B) | ~0.01° | Medieval system |
| Vehlow (V) | Exact | Equal variant |
| Horizontal (H) | ~0.01° | Falls back to Porphyry if iteration fails |
| Gauquelin (G) | ~0.01° | True 36-sector algorithm (18 diurnal + 18 nocturnal) |
| Krusinski-Pisa (U) | ~0.01° | Uses Porphyry fallback |
| Carter Poli-Equatorial (F) | ~0.01° | |
| APC Houses (Y) | ~0.01° | Uses Porphyry fallback |
| Natural Graduation (N) | ~0.01° | |
| Sripati (S) | ~0.01° | |

### Angle Precision

| Angle | Max Difference |
|-------|---------------|
| Ascendant | <0.001° |
| MC (Midheaven) | <0.001° |
| ARMC | <0.001° |
| Vertex | <0.01° |
| Equatorial Ascendant | <0.001° |
| Co-Ascendant (Koch) | <0.001° |
| Co-Ascendant (Munkasey) | <0.001° |
| Polar Ascendant | <0.001° |

### Polar Latitude Limitations

**Critical Limitation**: Some house systems fail at polar latitudes where the ecliptic does not properly intersect the horizon.

| Latitude Range | Behavior |
|----------------|----------|
| |lat| < 66.5° | All house systems work normally |
| |lat| > 66.5° | Placidus, Koch may fail |
| |lat| > 85° | Most time-based systems unreliable |
| |lat| = 89.9° | Only Equal, Whole Sign, Porphyry reliable |

**Fallback Behavior**: When Placidus or Koch calculations fail at polar latitudes, LibEphemeris automatically falls back to Porphyry house system.

**Systems that work at all latitudes**:
- Equal (A/E)
- Whole Sign (W)
- Porphyry (O)
- Morinus (M)
- Meridian (X)

---

## Ayanamsha (Sidereal Modes)

### Standard Ayanamshas

LibEphemeris supports all 43 Swiss Ephemeris ayanamsha modes.

| Category | Max Difference | Examples |
|----------|---------------|----------|
| Formula-based | <0.0002° | Fagan-Bradley, Lahiri, Raman |
| Epoch-based | <0.0002° | J2000, J1900, B1950 |
| Historical | <0.001° | Babylonian variants |

### Star-Based Ayanamshas

Star-based ayanamshas have slightly higher tolerance due to differences in:
- Star position calculations
- Precession models
- Proper motion calculations
- Coordinate transformations

| Ayanamsha | Max Difference |
|-----------|---------------|
| True Citra | <0.006° |
| True Revati | <0.06° |
| True Pushya | <0.06° |
| True Mula | <0.06° |
| True Sheoran | <0.06° |
| Galactic Center 0 Sag | <0.06° |
| Galactic Center Rgilbrand | <0.06° |
| Galactic Center Cochrane | <0.06° |
| Galactic Center Mula Wilhelm | <0.06° |

**Note on True Citra**: Uses high-precision Hipparcos coordinates for Spica
(HIP 65474) with full proper motion correction including parallax and radial
velocity. This achieves ~0.006° precision vs Swiss Ephemeris, significantly
better than other star-based ayanamshas.

---

## Crossing Events (Ingress Times)

### Sun Crossings

| Event | Precision |
|-------|-----------|
| Zodiac sign ingress | ±0.001 arcsec (sub-milliarcsecond) |
| Solstices | ±0.001 arcsec |
| Equinoxes | ±0.001 arcsec |
| Time difference vs pyswisseph | <1 second |

### Moon Crossings

| Event | Precision |
|-------|-----------|
| Zodiac sign ingress | ±0.05 arcsec (sub-arcsecond) |
| Node crossing (latitude=0) | ±0.05 arcsec |

### Planet Crossings

| Event | Precision |
|-------|-----------|
| Generic planet crossing | ±0.1 arcsec |
| Heliocentric crossing | ±0.1 arcsec |

**Note on Retrograde Stations**: Near retrograde stations where planetary velocity approaches zero, the algorithm uses Brent's method as a fallback to ensure convergence.

---

## Time Calculations

### Julian Day Precision

| Function | Precision |
|----------|-----------| 
| julday/revjul | <1e-10 days (~0.01 microseconds) |
| Delta T (with IERS) | <0.1 seconds (1973-present) |
| Delta T (modeled) | <3.6 seconds (historical/future dates) |

### IERS Delta T Support

LibEphemeris can optionally download observed Delta T values from IERS 
(International Earth Rotation and Reference Systems Service) for high-precision 
calculations on recent dates (1973-present):

```python
from libephemeris import set_iers_delta_t_enabled, download_delta_t_data
download_delta_t_data()  # Download IERS data
set_iers_delta_t_enabled(True)  # Enable IERS Delta T
```

### TAI Time Scale Support

LibEphemeris supports TAI (International Atomic Time) conversions:
- `utc_to_tai_jd()`: Convert UTC date/time to TAI Julian Day
- `tai_jd_to_utc()`: Convert TAI Julian Day to UTC date/time
- `tt_to_tai_jd()`: Convert TT to TAI (fixed 32.184s offset)
- `tai_to_tt_jd()`: Convert TAI to TT

TAI is a continuous time scale without leap seconds, forming the basis for 
TT (Terrestrial Time) where TT = TAI + 32.184 seconds exactly.

### Valid Date Ranges

LibEphemeris uses JPL DE ephemerides with specific date ranges:

| Ephemeris | Date Range | Precision |
|-----------|------------|-----------| 
| DE440 (default) | 1550-2650 | Full precision, ICRF 3.0 |
| DE421 | 1900-2050 | Full precision |
| DE422 | -3000 to 3000 | Full precision |
| DE430 | 1550-2650 | Full precision |
| DE431 | -13200 to 17191 | Full precision |
| DE441 | -13200 to 17191 | Full precision |

**Outside valid range**: Calculations will raise an exception with the supported date range.

---

## Minor Bodies and Extended Points

### Lunar Nodes and Lilith

LibEphemeris uses different calculation models for some lunar points, resulting in varying precision:

| Point | Max Difference vs pyswisseph | Mean Difference | Notes |
|-------|------------------------------|-----------------|-------|
| Mean Node | ~0.005° (~18 arcsec) | ~0.003° (~11 arcsec) | High precision |
| True Node | ~0.14° (~520 arcsec) | ~0.04° (~145 arcsec) | Osculating orbital elements method |
| Mean Lilith | ~0.12° (~430 arcsec) | ~0.07° (~250 arcsec) | Minor formula differences |
| True Lilith | ~5-15° | - | Different orbital model (see below) |
| Interpolated Apogee | ~8-10° | - | Smoothed osculating apogee (see below) |
| Interpolated Perigee | ~10-15° | - | Smoothed osculating perigee (see below) |

**Note on True Node**: The True Node calculation uses osculating orbital elements derived from the Moon's instantaneous position and velocity, with IAU 2006 precession correction. Mean error is ~0.04° (145 arcsec), max error ~0.14° (520 arcsec). For most astrological purposes, this precision is sufficient.

**Note on Mean Node**: The Mean Node now achieves high precision (~0.005° max), significantly better than the True Node calculation.

**Note on True Lilith**: The True Lilith (osculating apogee) calculation uses the eccentricity vector method derived from JPL DE ephemeris state vectors. This differs fundamentally from Swiss Ephemeris which uses an integrated analytical lunar theory. The difference arises because the "osculating apogee" concept is model-dependent when solar perturbations are significant. Typical error is 5-15°, with occasional differences up to 25°. For applications requiring closer Swiss Ephemeris compatibility, use Mean Lilith instead.

**Note on Interpolated Apogee/Perigee**: The interpolated apogee (SE_INTP_APOG) and perigee (SE_INTP_PERG) smooth out the spurious 30-degree oscillations in the osculating apsides using polynomial regression over a 56-day window. While Swiss Ephemeris uses an analytical method derived from Moshier's lunar theory, libephemeris uses least-squares fitting to achieve smoothing. The interpolated apogee oscillates ~5° from the mean (vs. 30° for osculating), making it significantly smoother. However, differences of 8-10° from Swiss Ephemeris persist due to the different underlying osculating calculations. See [INTERPOLATED_APOGEE.md](INTERPOLATED_APOGEE.md) for comprehensive documentation.

### Asteroids

| Asteroid | Precision | Notes |
|----------|-----------|-------|
| Chiron | ~1° | Keplerian approximation |
| Pholus | ~1° | Keplerian approximation |
| Ceres | ~1° | Keplerian approximation |
| Pallas | ~1° | Keplerian approximation |
| Juno | ~1° | Keplerian approximation |
| Vesta | ~1° | Keplerian approximation |

**Note**: Asteroids use Keplerian (two-body) orbital propagation which does not account for planetary perturbations. Errors grow over time, especially for dates far from the orbital elements epoch (2023.0).

### Trans-Neptunian Objects (TNOs)

LibEphemeris now includes first-order secular perturbations from Jupiter, Saturn,
Uranus, and Neptune for all TNOs, providing significantly improved accuracy over
pure Keplerian propagation.

| TNO | Precision | Notes |
|-----|-----------|-------|
| Eris | <10° | Keplerian + secular perturbations, a~68 AU |
| Makemake | <10° | Keplerian + secular perturbations, a~45 AU |
| Ixion | <10° | Plutino (2:3 resonance with Neptune), a~39 AU |
| Orcus | <10° | Plutino (2:3 resonance with Neptune), a~39 AU |
| Haumea | <10° | Keplerian + secular perturbations |
| Quaoar | <10° | Keplerian + secular perturbations |
| Sedna | ~10° | Extremely distant (a~550 AU), minimal perturbations |

#### Secular Perturbation Model

The TNO calculation model applies Laplace-Lagrange secular perturbation theory:

- **Jupiter & Saturn**: Dominant perturbations for all bodies
- **Uranus**: Significant for TNOs, comparable to Saturn's influence at large distances
- **Neptune**: Critical for TNOs, especially plutinos in 2:3 mean motion resonance

Secular perturbations correct for:
- Perihelion precession (dω/dt)
- Node regression (dΩ/dt)

Example perturbation rates (arcsec/year):
| TNO | dω/dt | dΩ/dt |
|-----|-------|-------|
| Eris (68 AU) | +0.08 | -0.02 |
| Makemake (45 AU) | +0.35 | -0.15 |
| Ixion (39 AU, plutino) | +0.94 | -0.41 |
| Orcus (39 AU, plutino) | +0.93 | -0.41 |

#### Resonance Detection

LibEphemeris can detect bodies in Neptune mean motion resonances:

```python
from libephemeris.minor_bodies import detect_mean_motion_resonance, MINOR_BODY_ELEMENTS
from libephemeris.constants import SE_IXION

result = detect_mean_motion_resonance(MINOR_BODY_ELEMENTS[SE_IXION])
print(result.is_resonant)  # True
print(result.resonance.name)  # "plutino"
```

**Note on Resonant Bodies**: Plutinos like Ixion and Orcus are in 2:3 mean motion
resonance with Neptune. While secular perturbation theory is applied, resonant
dynamics are not fully captured. For research-grade precision, use SPK kernels.

#### Validation (2000-2050)

TNO positions were validated against pyswisseph over a 50-year span:

- Maximum longitude error: <10° (within tolerance for astrological use)
- Mean longitude error: <5° (typical case)
- Latitude error: <5°

**Note**: For high-precision TNO work requiring sub-arcminute accuracy, use
dedicated SPK ephemeris files from JPL Horizons.

---

## Planetary Moons

LibEphemeris supports calculating positions of planetary moons (Galilean moons,
Titan, Triton, Phobos, Deimos, Charon, etc.) using JPL satellite SPK files.

### Supported Moons

| System | Moons | SPK File |
|--------|-------|----------|
| Jupiter | Io, Europa, Ganymede, Callisto | jup365.bsp |
| Saturn | Titan, Enceladus, Mimas, Rhea, Dione, Tethys, Iapetus | sat441.bsp |
| Uranus | Titania, Oberon, Ariel, Umbriel, Miranda | ura116.bsp |
| Neptune | Triton | nep102.bsp |
| Mars | Phobos, Deimos | mar097.bsp |
| Pluto | Charon | plu058.bsp |

### Precision

| Component | Precision |
|-----------|-----------|
| Longitude | ±1 arcsec |
| Latitude | ±1 arcsec |
| Distance | <0.001 AU |

### Usage

```python
import libephemeris as eph
eph.register_moon_spk('jup365.bsp')
pos, _ = eph.calc_ut(2451545.0, eph.SE_MOON_IO, eph.SEFLG_SPEED)
```

**Note**: SPK kernels must be downloaded separately from JPL Horizons and
registered using `register_moon_spk()` before calculating moon positions.

---

## Fixed Stars

### Position Precision

| Component | Precision |
|-----------|-----------|
| Longitude | ~0.01 arcsec |
| Latitude | ~0.01 arcsec |
| Proper motion | Corrected for stellar curvature |
| Nutation | IAU 2000A model (1365 terms) |

### Velocity Limitation

Fixed star velocity calculations return 0:

```python
pos, _ = ephem.fixstar_ut("Aldebaran", jd, SEFLG_SPEED)
# pos[3], pos[4], pos[5] are 0.0 (velocities not implemented)
```

---

## Eclipse Functions

### Implemented Features

| Function | Precision |
|----------|-----------| 
| Solar eclipse timing | < 10 seconds |
| Lunar eclipse timing | < 10 seconds |
| Occultation timing | ~1-2 minutes |
| Eclipse magnitude | ~0.01 |

### High-Precision Eclipse Timing

Solar and lunar eclipse times are calculated using proper Besselian elements,
achieving timing precision of better than 10 seconds compared to Swiss Ephemeris.
This is accomplished by:

1. **Refined Eclipse Maximum**: The time of maximum eclipse is found by minimizing
   the gamma parameter (shadow axis distance from Earth center) using golden
   section search, achieving sub-second precision.

2. **Accurate Contact Times**: First/fourth contacts (penumbral limits) and
   second/third contacts (umbral limits) are calculated by finding when the
   shadow cone boundaries cross Earth's limb using Besselian element geometry.

3. **Full Ephemeris Precision**: All calculations use high-precision JPL DE
   ephemeris positions for the Moon and Sun.

### Not Implemented

The following return 0 or placeholder values:

- Saros series number
- Inex number
- Sunrise/sunset on central line for solar eclipses

---

## Performance Considerations

### Calculation Speed

LibEphemeris (pure Python with Skyfield) is approximately 10-100x slower than pyswisseph (C library). For most astrological applications, this is negligible. For batch processing of millions of charts, consider:

1. Using multi-threading with `EphemerisContext`
2. Caching frequently-used calculations
3. Using pyswisseph for performance-critical paths

### Memory Usage

| Resource | Approximate Size |
|----------|-----------------| 
| DE440 ephemeris | ~128 MB |
| DE421 ephemeris | ~16 MB |
| DE431 ephemeris | ~3.4 GB |
| Timescale data | ~2 MB |
| Per EphemerisContext | ~1 KB |

Ephemeris files are shared across all `EphemerisContext` instances for memory efficiency.

---

## Summary of Use Case Suitability

### Suitable For

- Natal chart calculations (standard latitudes)
- Transit calculations
- Synastry and composite charts
- Progression calculations
- Most astrological software applications
- Research requiring reproducible, documented precision
- Planetary moon calculations (with SPK kernels)

### Use With Caution

- Polar latitude locations (>66.5°) - use Equal or Whole Sign houses
- Star-based ayanamshas - expect ±0.006° for True Citra, ±0.06° for others
- Asteroid/TNO positions - Keplerian approximation only
- True Node calculations - ~0.14° precision vs Swiss Ephemeris
- Very old historical dates (before 1550 with DE440)

### Not Suitable For

- Observatory-grade positional astronomy (use NOVAS or SOFA)
- Spacecraft navigation (use SPICE)
- Precise eclipse predictions for scientific purposes
- High-precision Saros/Inex series analysis

---

## Precision Validation

LibEphemeris includes comprehensive precision validation tooling to verify accuracy
against pyswisseph reference values.

### Precision Report Generator

Generate detailed precision reports using the included script:

```bash
# Generate text report
python compare_scripts/generate_precision_report.py

# Generate JSON report for automation
python compare_scripts/generate_precision_report.py --json

# Generate CSV report for analysis
python compare_scripts/generate_precision_report.py --csv

# Increase test samples for higher confidence
python compare_scripts/generate_precision_report.py -n 500
```

The report includes:
- Maximum, mean, and standard deviation for each calculation type
- P95 and P99 percentile values
- Pass/fail status against documented tolerances
- Support for JSON and CSV output formats

### Precision Test Suite

Run the precision documentation tests to verify all claims in this document:

```bash
pytest tests/test_precision/test_precision_docs.py -v
```

---

## References

1. NASA JPL DE Ephemeris documentation
2. Swiss Ephemeris documentation and source code
3. Meeus, "Astronomical Algorithms" 2nd Edition
4. IERS Conventions 2003 (nutation models)
5. IAU SOFA Library documentation
6. Skyfield documentation (Rhodes Mill)
