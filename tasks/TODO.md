# Exhaustive Test Plan — LibEphemeris v1.0.0a3+

Comprehensive test checklist for full validation across all backends, body types,
flag combinations, house systems, and precision tiers. Designed to verify 1:1
compatibility with the reference implementation and internal consistency across
the Skyfield, LEB, and Horizons calculation paths.

**Context:** LibEphemeris has 3 backends (Skyfield/DE440, LEB precomputed Chebyshev,
Horizons API), 4 calc modes (auto/skyfield/leb/horizons), supports 40+ body types,
10+ flag combos, 24 house systems, 43 ayanamsha modes, eclipses, rise/set,
heliacal visibility, and more. This plan ensures no regressions after the
tech-debt cleanup (except narrowing, thread safety, dedup, Skyfield TypeError fix).

---

## 1. Cross-Backend Position Accuracy

### 1.1 Skyfield vs Reference (pyswisseph) — per body per flag

Bodies: SE_SUN(0), SE_MOON(1), SE_MERCURY(2), SE_VENUS(3), SE_MARS(4),
SE_JUPITER(5), SE_SATURN(6), SE_URANUS(7), SE_NEPTUNE(8), SE_PLUTO(9),
SE_EARTH(14), SE_MEAN_NODE(10), SE_TRUE_NODE(11), SE_MEAN_APOG(12),
SE_OSCU_APOG(13), SE_CHIRON(15), SE_CERES(17), SE_PALLAS(18), SE_JUNO(19),
SE_VESTA(20), SE_INTP_APOG(21), SE_INTP_PERG(22)

Flags:
- [ ] SEFLG_SPEED (default)
- [ ] SEFLG_SPEED | SEFLG_SIDEREAL (Lahiri)
- [ ] SEFLG_SPEED | SEFLG_EQUATORIAL
- [ ] SEFLG_SPEED | SEFLG_J2000
- [ ] SEFLG_SPEED | SEFLG_NOABERR
- [ ] SEFLG_SPEED | SEFLG_HELCTR
- [ ] SEFLG_SPEED | SEFLG_BARYCTR
- [ ] SEFLG_SPEED | SEFLG_TRUEPOS
- [ ] SEFLG_SPEED | SEFLG_NONUT
- [ ] SEFLG_SPEED | SEFLG_NOGDEFL
- [ ] SEFLG_SPEED | SEFLG_XYZ
- [ ] SEFLG_SPEED | SEFLG_RADIANS
- [ ] SEFLG_SPEED | SEFLG_TOPOCTR (requires set_topo first)

Dates: 500 random JDs in [2415020, 2488069] (1900-2100)

**Commands:**
```bash
LIBEPHEMERIS_COMPARE_MODE=skyfield poe test:compare
```

### 1.2 LEB vs Skyfield — all 14 core bodies

- [ ] Default flags, 500 dates, threshold < 0.005"
- [ ] Sidereal flags, 200 dates
- [ ] Equatorial flags, 200 dates
- [ ] J2000 flags, 200 dates
- [ ] Heliocentric flags, 200 dates

**Commands:**
```bash
poe test:unit:leb:fast
poe test:leb2:precision:base
poe test:leb2:precision:medium
```

### 1.3 LEB vs Reference (pyswisseph)

- [ ] All compare tests in LEB mode, all bodies, no @slow
- [ ] Same with @slow included

**Commands:**
```bash
LIBEPHEMERIS_COMPARE_MODE=leb poe test:compare:leb
LIBEPHEMERIS_COMPARE_MODE=leb poe test:compare:leb:full
```

### 1.4 Horizons vs Skyfield — 13 bodies x 6 flags

- [ ] 200 dates geocentric, threshold < 0.001"
- [ ] 50 dates heliocentric
- [ ] Cross-validate Horizons vs LEB2

**Commands:**
```bash
poe test:horizons
poe test:horizons:quick
poe test:horizons:vs:leb
```

### 1.5 Horizons vs Reference (pyswisseph)

- [ ] All compare tests in Horizons mode

**Commands:**
```bash
poe test:compare:horizons
```

---

## 2. House Systems (24 systems)

Systems: P K O R C E W X M H T B G I i A L N Q Y F U D J

### 2.1 Per-system accuracy vs reference

- [ ] 10 dates x 6 locations (equator, mid-lat, high-lat, polar) per system
- [ ] Cusps 1-12 all in [0, 360)
- [ ] ASC and MC in [0, 360)
- [ ] Vertex, East Point computed

### 2.2 Sidereal houses

- [ ] All systems with SEFLG_SIDEREAL (Lahiri, Fagan-Bradley, Raman)
- [ ] swe_houses_ex() with sidereal flag

### 2.3 Polar latitudes

- [ ] Placidus/Koch fallback to Porphyry above 66.5 deg
- [ ] Equal/Whole Sign work at all latitudes
- [ ] PolarCircleError raised where expected

### 2.4 swe_house_pos()

- [ ] Planet in correct house for known chart
- [ ] All 24 systems

**Commands:**
```bash
pytest compare_scripts/tests/test_all_houses.py -v
pytest compare_scripts/tests/test_compare_houses.py -v
pytest compare_scripts/tests/test_compare_houses_ext.py -v
```

---

## 3. Ayanamsha Modes (43 modes)

- [ ] All 43 modes produce different sidereal offsets
- [ ] swe_get_ayanamsa_ut() returns finite value for all modes
- [ ] Sidereal longitude = tropical - ayanamsha (within tolerance)
- [ ] Custom ayanamsha (mode 255) with user-defined t0, ayan_t0
- [ ] Star-based ayanamshas (modes 17, 27, 28, 29) use star positions

**Commands:**
```bash
pytest compare_scripts/tests/test_ayanamsha.py -v
pytest compare_scripts/tests/test_ayanamsha_all_modes.py -v
pytest compare_scripts/tests/test_compare_sidereal.py -v
pytest compare_scripts/tests/test_compare_sidereal_regression.py -v
```

---

## 4. Eclipse Calculations

### 4.1 Solar eclipses

- [ ] swe_sol_eclipse_when_glob() finds eclipses for known years
- [ ] swe_sol_eclipse_when_loc() at specific location
- [ ] swe_sol_eclipse_where() geographic coordinates
- [ ] Eclipse type matches reference (total, annular, partial)
- [ ] Timing within 6 seconds of reference

### 4.2 Lunar eclipses

- [ ] swe_lun_eclipse_when() finds eclipses for known years
- [ ] swe_lun_eclipse_when_loc() at specific location
- [ ] Eclipse type matches reference
- [ ] Timing within 8 seconds of reference

### 4.3 Planet occultations

- [ ] swe_lun_occult_when_glob() for major bodies
- [ ] swe_lun_occult_when_loc() at specific location

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_eclipses.py -v
pytest compare_scripts/tests/test_compare_occultations.py -v
pytest compare_scripts/tests/test_compare_planet_occultations.py -v
```

---

## 5. Rise/Set/Transit

- [ ] swe_rise_trans() for Sun, Moon, planets
- [ ] Multiple locations (equator, mid-lat, polar)
- [ ] Sun rise/set times within 1 minute of reference
- [ ] Moon rise/set times within 2 minutes
- [ ] Transit times (meridian crossing)
- [ ] Circumpolar handling at high latitudes

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_rise_transit.py -v
```

---

## 6. Heliacal Visibility

- [ ] swe_heliacal_ut() for major planets and stars
- [ ] swe_vis_limit_mag() visual magnitude calculation
- [ ] Different atmospheric conditions
- [ ] Heliacal rising/setting types (ACRONYCHAL, COSMICAL, EVENING, MORNING)

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_heliacal.py -v
```

---

## 7. Fixed Stars

- [ ] swe_fixstar_ut() for 50+ known stars
- [ ] swe_fixstar_mag() magnitude
- [ ] Proper motion over 100+ years
- [ ] All stars from Hipparcos catalog accessible
- [ ] Star-based ayanamsha stars (Spica, Revati, etc.)

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_fixedstars.py -v
```

---

## 8. Coordinate Transforms

- [ ] swe_cotrans() ecliptic <-> equatorial round-trip
- [ ] swe_cotrans_sp() with speed components
- [ ] Multiple obliquity values (23.0, 23.4393, 24.0)
- [ ] Edge cases: poles, equinoxes, solstices

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_coordinates.py -v
```

---

## 9. Time Functions

- [ ] julday/revjul round-trip for 1000 random dates (-5000 to +5000)
- [ ] deltat() returns positive value for modern dates
- [ ] deltat_ex() returns source information
- [ ] sidtime() matches reference
- [ ] utc_to_jd / jdet_to_utc / jdut1_to_utc round-trip
- [ ] date_conversion validation
- [ ] day_of_week correctness

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_time.py -v
```

---

## 10. Lunar Calculations

### 10.1 Nodes (Mean and True)

- [ ] Mean Node smooth retrograde motion
- [ ] True Node oscillation around mean
- [ ] Distance values match reference

### 10.2 Apsides (Perigee/Apogee/Lilith)

- [ ] Mean Apogee (Black Moon Lilith)
- [ ] Osculating Apogee
- [ ] Interpolated Apogee (SE_INTP_APOG)
- [ ] Interpolated Perigee (SE_INTP_PERG)
- [ ] Perigee/apogee correction terms

### 10.3 Lunar eclipses from lunar.py directly

- [ ] Moon at perigee and apogee dates
- [ ] Distance extrema match physical values

**Commands:**
```bash
poe test:lunar
poe test:lunar:perigee
poe test:lunar:apogee
poe test:lunar:lilith
pytest compare_scripts/tests/test_compare_lunar.py -v
pytest compare_scripts/tests/test_compare_lunar_nodes_lilith.py -v
```

---

## 11. Minor Bodies

- [ ] Chiron, Ceres, Pallas, Juno, Vesta positions
- [ ] SPK-based precision for available bodies
- [ ] Keplerian fallback for bodies without SPK
- [ ] ASSIST n-body fallback (if installed)
- [ ] Auto SPK download (if network available)

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_minor_bodies.py -v
pytest compare_scripts/tests/test_compare_hypothetical.py -v
```

---

## 12. Uranian / Hypothetical Bodies

- [ ] Cupido through Poseidon (40-47) heliocentric
- [ ] Cupido through Poseidon (40-47) geocentric
- [ ] Transpluto / Isis (48)
- [ ] All with SEFLG_SPEED, SEFLG_SIDEREAL, SEFLG_EQUATORIAL

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_hypothetical.py -v
```

---

## 13. Heliocentric and Barycentric

- [ ] All planets heliocentric vs reference
- [ ] All planets barycentric vs reference
- [ ] Sun heliocentric = (0,0,0)
- [ ] Earth heliocentric = opposite of Sun geocentric

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_helio_bary.py -v
```

---

## 14. Observations (azalt, elongation)

- [ ] swe_azalt() for Sun, Moon, planets at known locations
- [ ] swe_azalt_rev() reverse azimuth to ecliptic
- [ ] swe_pheno_ut() phase angle, elongation, magnitude
- [ ] Solar elongation for inner planets

**Commands:**
```bash
pytest compare_scripts/tests/test_azalt.py -v
pytest compare_scripts/tests/test_azalt_rev.py -v
pytest compare_scripts/tests/test_compare_elongation.py -v
pytest compare_scripts/tests/test_compare_observations.py -v
```

---

## 15. Crossings and Phenomena

- [ ] swe_solcross_ut() solstice/equinox times
- [ ] swe_mooncross_ut() lunar crossing times
- [ ] swe_crossings() general body crossings

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_crossings.py -v
```

---

## 16. LEB-Specific Tests

### 16.1 LEB1/LEB2 Reader

- [ ] open_leb() auto-detects LEB1 vs LEB2
- [ ] LEB2Reader lazy decompression
- [ ] CompositeLEBReader multi-file dispatch
- [ ] eval_body() for all 14 core bodies
- [ ] eval_nutation() matches Skyfield
- [ ] delta_t() matches Skyfield
- [ ] get_star() for fixed stars in LEB

### 16.2 LEB2 Compression

- [ ] shuffle/unshuffle round-trip
- [ ] compress/decompress round-trip
- [ ] Precision vs LEB1 < 0.001" for all bodies

### 16.3 LEB Precision per tier

- [ ] Base tier: all core bodies, max error < 0.003"
- [ ] Medium tier: all core bodies
- [ ] Extended tier: all core bodies

**Commands:**
```bash
poe test:leb2
poe test:leb2:precision:base
poe test:leb2:precision:medium
poe test:leb2:precision:all
poe test:unit:leb:fast
```

---

## 17. Horizons-Specific Tests

### 17.1 HTTP Client

- [ ] LRU cache hit/miss
- [ ] Parallel fetch (multiple bodies same JD)
- [ ] Retry with exponential backoff
- [ ] Timeout handling
- [ ] Cache clear and shutdown

### 17.2 Pipeline

- [ ] Light-time iteration
- [ ] Gravitational deflection (Sun, Jupiter, Saturn)
- [ ] Stellar aberration
- [ ] Frame rotation (ICRS -> ecliptic of date)
- [ ] Analytical bodies (Mean Node, Mean Apogee) — no HTTP

### 17.3 Unsupported bodies fallback

- [ ] True Node -> Skyfield fallback
- [ ] Oscu Apogee -> Skyfield fallback
- [ ] SEFLG_TOPOCTR -> Skyfield fallback

**Commands:**
```bash
poe test:horizons
poe test:horizons:quick
poe test:horizons:vs:leb
```

---

## 18. State Management

### 18.1 Setter/getter consistency

- [ ] set_calc_mode -> get_calc_mode round-trip
- [ ] set_topo -> get_topo round-trip
- [ ] set_sid_mode -> get_sid_mode round-trip (all 43 modes)
- [ ] set_ephe_path -> get_ephe_path
- [ ] set_precision_tier -> get_precision_tier
- [ ] close() resets all state
- [ ] Environment variables override programmatic settings

### 18.2 Thread safety

- [ ] Concurrent calc_ut with different sidereal modes (via EphemerisContext)
- [ ] Concurrent calc_ut with same JD, different bodies
- [ ] close() during in-flight calculation
- [ ] LEB reader concurrent access

---

## 19. EphemerisContext (Thread-Safe API)

- [ ] ctx.calc_ut() matches global calc_ut()
- [ ] ctx.set_topo() isolated per context
- [ ] ctx.set_sid_mode() isolated per context
- [ ] ctx.set_leb_file() per-context LEB reader
- [ ] Multiple contexts concurrent, different settings
- [ ] Context manager (__enter__ / __exit__)

---

## 20. Edge Cases

- [ ] JD at exact tier boundaries (start/end of DE440 range)
- [ ] JD = 0 (historical date)
- [ ] Body ID = -1 (SE_ECL_NUT nutation)
- [ ] Body ID = 999 (unknown -> error)
- [ ] NaN / Inf JD -> error
- [ ] Latitude = 90 (North Pole houses)
- [ ] Latitude = -90 (South Pole)
- [ ] Longitude = 180 / -180
- [ ] Very fast Moon motion near perigee
- [ ] Retrograde planets (Mars, Mercury)

---

## 21. Orbital Elements

- [ ] swe_nod_aps_ut() for all planets
- [ ] swe_orbit_max_min_true_distance() for all planets
- [ ] Ascending/descending node
- [ ] Perihelion/aphelion

**Commands:**
```bash
pytest compare_scripts/tests/test_compare_orbital.py -v
```

---

## 22. Validation Suite (automated)

- [ ] run_validation_suite.py --quick (8/8 PASS)
- [ ] run_validation_suite.py full (8/8 PASS)

**Commands:**
```bash
python scripts/run_validation_suite.py --quick
python scripts/run_validation_suite.py
```

---

## 23. Golden Regression

- [ ] test_golden_regression.py matches regenerated golden file
- [ ] No unexpected position changes

**Commands:**
```bash
pytest tests/test_golden_regression.py -v
```

---

## 24. CLI and Download

- [ ] `libephemeris status` works
- [ ] `libephemeris --version` shows correct version
- [ ] `libephemeris download:base` (if network)
- [ ] download_leb2_for_tier() downloads from GitHub Release

---

## Execution Commands Summary

```bash
# Unit tests (no slow)
poe test:unit:fast

# Unit tests in LEB mode
poe test:unit:leb:fast

# LEB2 precision (all tiers)
poe test:leb2:precision:all

# Horizons precision
poe test:horizons

# Compare vs reference (Skyfield mode)
poe test:compare

# Compare vs reference (LEB mode)
poe test:compare:leb

# Compare vs reference (Horizons mode)
poe test:compare:horizons

# Lunar tests
poe test:lunar

# Validation suite
python scripts/run_validation_suite.py

# Golden regression
pytest tests/test_golden_regression.py -v

# Full compare (all tests including slow)
poe test:compare:full
```
