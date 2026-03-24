# Validation Plan v2 — Next-Phase Quality Assurance

> Created: March 2026
> Status: Active
> Prerequisite: Hyper-validation v1 complete (4400 rounds, 0 FAIL)

## Overview

With pyswisseph 1:1 API compatibility verified and the full test suite passing
(4066 tests, 0 FAIL), this plan covers the next phase of quality assurance:
independent cross-validation against authoritative sources, LEB accuracy
verification, and robustness testing.

---

## 1. JPL Horizons Cross-Validation

**Goal:** Verify libephemeris positions against NASA JPL Horizons (the
authoritative source) independently of pyswisseph.

### 1.1 Planetary Positions (10 bodies × 50 dates)

- [ ] Query Horizons for geocentric apparent ecliptic lon/lat for Sun–Pluto
- [ ] Date range: 1900–2100, including solstices, equinoxes, and random dates
- [ ] Expected tolerance: < 0.5" for all planets, < 1" for Moon
- [ ] Automate via `astroquery.jplhorizons` or Horizons batch API
- [ ] Save results as reproducible JSON for regression

### 1.2 Lunar Node and Lilith

- [ ] True Node: verify < 0.01" vs Horizons (already spot-checked at 24 dates)
- [ ] Extend to 100 dates spanning 1800–2200
- [ ] Mean Node: verify polynomial accuracy vs Horizons osculating elements

### 1.3 Minor Bodies with SPK

- [ ] Chiron, Ceres, Pallas, Juno, Vesta: 20 dates each
- [ ] Expected: < 0.01" with SPK kernels
- [ ] Verify fallback cascade: SPK → ASSIST → Keplerian

### 1.4 Outer Planet COB Corrections

- [ ] Jupiter, Saturn, Neptune, Pluto body-center positions vs Horizons
- [ ] Verify analytical fallback accuracy outside SPK coverage
- [ ] 20 dates per planet, including dates near/outside SPK boundaries

**Deliverable:** `scripts/horizons_cross_validate.py` with JSON report

---

## 2. LEB Accuracy Sweep

**Goal:** Verify that LEB (binary ephemeris) mode reproduces Skyfield mode
to within documented tolerances across the full parameter space.

### 2.1 Position Accuracy

- [ ] All 31 LEB bodies at 1000 random dates within tier range
- [ ] Tolerance: < 1" for planets, < 5" for ecliptic bodies
- [ ] Compare LEB vs Skyfield (not vs pyswisseph)
- [ ] Test all three tiers: base, medium, extended

### 2.2 Flag Combinations

- [ ] Verify LEB fallback triggers correctly for unsupported flags:
      SEFLG_TOPOCTR, SEFLG_XYZ, SEFLG_RADIANS, SEFLG_NONUT
- [ ] Verify transparent fallback produces identical results to pure Skyfield
- [ ] Test SEFLG_SPEED, SEFLG_EQUATORIAL, SEFLG_J2000, SEFLG_SIDEREAL in LEB

### 2.3 Boundary Conditions

- [ ] Test dates at exact LEB segment boundaries (Chebyshev interval edges)
- [ ] Test dates at per-body range boundaries (asteroid coverage limits)
- [ ] Test graceful fallback when body is outside LEB range

### 2.4 Performance Regression

- [ ] Benchmark LEB vs Skyfield: target ≥ 10x speedup
- [ ] Measure cold-start (first read) vs warm (cached) performance
- [ ] Profile memory usage with all three tier files loaded

**Deliverable:** `scripts/leb_accuracy_sweep.py` with per-body error report

---

## 3. Property-Based Testing

**Goal:** Use hypothesis-style random testing to find edge cases that
structured tests miss.

### 3.1 Coordinate Transform Round-Trips

- [ ] ecliptic → equatorial → ecliptic: verify identity within 0.001"
- [ ] ecliptic → XYZ → ecliptic: verify identity
- [ ] degrees → radians → degrees: verify identity
- [ ] Random JD, random body, random flag combinations

### 3.2 API Contract Testing

- [ ] `calc_ut(jd, body, flags)` always returns (6-tuple, int)
- [ ] `houses(jd, lat, lon, system)` always returns (tuple, tuple)
- [ ] `fixstar(name, jd, flags)` always returns (6-tuple, str)
- [ ] No function raises unhandled exceptions for valid inputs
- [ ] All functions accept SEFLG_MOSEPH without error

### 3.3 Monotonicity Properties

- [ ] Sun longitude is monotonically increasing (mod 360) over 1 day
- [ ] Mean Node longitude is monotonically decreasing over 1 day
- [ ] Julian day conversion is monotonic: julday(y,m,d) < julday(y,m,d+1)
- [ ] Delta T is monotonically increasing for dates > 1972

### 3.4 Symmetry Properties

- [ ] `calc_ut(jd, body, SEFLG_HELCTR)` for Sun returns (0, 0, 0, ...)
- [ ] `pheno_ut(jd, SE_SUN)` returns phase_angle = 0
- [ ] `houses(jd, lat, lon)` cusps are in ascending order (mod 360)

**Deliverable:** `tests/test_property_based.py` using hypothesis library

---

## 4. Fuzz Testing

**Goal:** Verify robustness against malformed, extreme, and boundary inputs.

### 4.1 Extreme Julian Dates

- [ ] JD = 0 (4713 BC)
- [ ] JD = -1e6 (far past, outside all ephemeris ranges)
- [ ] JD = 1e8 (far future)
- [ ] JD = NaN, Inf, -Inf
- [ ] Verify: either valid result or clean exception (no crash/hang)

### 4.2 Invalid Body IDs

- [ ] Negative body IDs
- [ ] Body IDs > 100000 (SE_AST_OFFSET range)
- [ ] Body ID = SE_AST_OFFSET + 0 (edge case)
- [ ] Verify: UnknownBodyError or clean error

### 4.3 Extreme Geographic Coordinates

- [ ] Latitude = ±90° (poles)
- [ ] Latitude = ±91° (invalid)
- [ ] Longitude = ±180°, ±360°, ±720°
- [ ] Altitude = -1000m, 0m, 100000m
- [ ] Verify house calculations handle all gracefully

### 4.4 Flag Exhaustion

- [ ] All 2^14 combinations of the 14 SEFLG flags
- [ ] Verify no crashes, hangs, or unhandled exceptions
- [ ] Identify and document any invalid/conflicting flag combinations

**Deliverable:** `tests/test_fuzz.py`

---

## 5. Eclipse Event Validation

**Goal:** Cross-validate eclipse predictions against NASA eclipse catalogs.

### 5.1 Solar Eclipses (Espenak Canon)

- [ ] 20 historical solar eclipses (1900–2025) from NASA Eclipse Website
- [ ] Compare: type (total/annular/partial), maximum time, magnitude
- [ ] Tolerance: < 60 seconds timing, correct type classification

### 5.2 Lunar Eclipses

- [ ] 20 historical lunar eclipses from NASA catalog
- [ ] Compare: type, contact times (P1, U1, U2, U3, U4, P4), magnitude
- [ ] Verify gamma sign and magnitude

### 5.3 Future Eclipses

- [ ] Next 10 solar eclipses (2026–2035)
- [ ] Next 10 lunar eclipses (2026–2035)
- [ ] Compare against Horizons-derived circumstances

**Deliverable:** `tests/test_eclipse_catalog_validation.py`

---

## 6. Thread Safety and Concurrency

**Goal:** Verify that `EphemerisContext` is truly thread-safe under load.

### 6.1 Concurrent Context Stress Test

- [ ] 50 threads, each with own EphemerisContext
- [ ] Each thread computes 100 random positions
- [ ] Verify all results match single-threaded baseline
- [ ] No deadlocks, no data races, no crashes

### 6.2 Global State Isolation

- [ ] Thread A sets sid_mode(LAHIRI), Thread B sets sid_mode(FAGAN_BRADLEY)
- [ ] Verify each thread gets its own ayanamsha
- [ ] Test with `set_topo()`, `set_sid_mode()`, `set_jpl_file()`

### 6.3 LEB + Skyfield Mixed Mode

- [ ] Some contexts use LEB, others force Skyfield
- [ ] Verify no interference between modes
- [ ] Test mode switching mid-computation

**Deliverable:** `tests/test_concurrency_stress.py`

---

## 7. Regression Test Infrastructure

**Goal:** Establish infrastructure to prevent regressions.

### 7.1 Golden File Tests

- [ ] Generate reference outputs for 100 representative calculations
- [ ] Store as JSON golden files in `tests/golden/`
- [ ] CI runs golden file comparison on every commit
- [ ] Any difference triggers explicit review

### 7.2 Performance Benchmarks

- [ ] Establish baseline timings for key operations
- [ ] `calc_ut` (Skyfield): target < 200 us
- [ ] `calc_ut` (LEB): target < 15 us
- [ ] `houses`: target < 500 us
- [ ] CI alerts on > 20% regression

### 7.3 Hyper-Validation CI Integration

- [ ] Run hyper-validation (fast sections only) as CI check
- [ ] Sections A–AC excluding J (occultations) and U (heliacal)
- [ ] Target: < 10 minutes runtime
- [ ] 0 FAIL gate on merge to main

**Deliverable:** CI configuration, `tests/golden/` directory

---

## Priority Order

1. **JPL Horizons Cross-Validation** (§1) — highest value, independent ground truth
2. **LEB Accuracy Sweep** (§2) — user-facing feature needs verification
3. **Property-Based Testing** (§3) — catches systematic blind spots
4. **Eclipse Catalog Validation** (§5) — high-visibility feature
5. **Fuzz Testing** (§4) — robustness
6. **Thread Safety** (§6) — important for production use
7. **Regression Infrastructure** (§7) — long-term maintenance

## Success Criteria

- All §1 tests pass with documented tolerances
- LEB accuracy confirmed < 1" for all planetary bodies
- No crashes from any fuzz input
- Eclipse timing within 60 seconds of NASA catalog
- Thread safety verified under 50-thread stress
- Golden file infrastructure in CI
