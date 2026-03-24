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

- [x] Query Horizons for geocentric apparent ecliptic lon/lat for Sun–Pluto
- [x] Date range: 1900–2100, including solstices, equinoxes, and random dates
- [x] Expected tolerance: < 0.5" for all planets, < 1" for Moon (core era 1972–2040; era-adaptive tolerances for historical/extrapolated dates)
- [x] Automate via `astroquery.jplhorizons` or Horizons batch API
- [x] Save results as reproducible JSON for regression

### 1.2 Lunar Node and Lilith

- [x] True Node: validated via pyswisseph hyper-validation (4400 rounds, 0 FAIL)
- [x] Mean Node: validated via pyswisseph hyper-validation
- [x] Note: These are derived astrological quantities not available as Horizons body positions

### 1.3 Minor Bodies with SPK

- [x] Chiron, Ceres, Pallas, Juno, Vesta: 20 dates each
- [x] All pass within 0.12" with SPK kernels (max 0.11")
- [x] Verify fallback cascade: SPK → ASSIST → Keplerian

### 1.4 Outer Planet COB Corrections

- [x] Jupiter, Saturn, Neptune, Pluto body-center positions vs Horizons
- [x] All pass within 0.10" (core era max 0.10")
- [x] 20 dates per planet, including dates near/outside SPK boundaries

**Deliverable:** `scripts/horizons_cross_validate.py` with JSON report ✅

---

## 2. LEB Accuracy Sweep

**Goal:** Verify that LEB (binary ephemeris) mode reproduces Skyfield mode
to within documented tolerances across the full parameter space.

### 2.1 Position Accuracy

- [x] All 31 LEB bodies at 1000 random dates within tier range
- [x] Tolerance: < 1" for planets, < 5" for ecliptic bodies (all pass, max 0.000352")
- [x] Compare LEB vs Skyfield (not vs pyswisseph)
- [x] Test all three tiers: base, medium, extended (medium verified)

### 2.2 Flag Combinations

- [x] Verify LEB fallback triggers correctly for unsupported flags:
      SEFLG_TOPOCTR, SEFLG_XYZ, SEFLG_RADIANS, SEFLG_NONUT
- [x] Verify transparent fallback produces identical results to pure Skyfield
- [x] Test SEFLG_SPEED, SEFLG_EQUATORIAL, SEFLG_J2000, SEFLG_SIDEREAL in LEB

### 2.3 Boundary Conditions

- [x] Test dates at exact LEB segment boundaries (Chebyshev interval edges)
- [x] Test dates at per-body range boundaries (asteroid coverage limits)
- [x] Test graceful fallback when body is outside LEB range

### 2.4 Performance Regression

- [x] Benchmark LEB vs Skyfield: ≥8x speedup for JPL-backed planets (achieved 10–23x)
- [x] Measure cold-start (first read) vs warm (cached) performance
- [x] Polynomial bodies (MeanNode) exempt from speedup target (Skyfield already ~50µs)

**Deliverable:** `scripts/leb_accuracy_sweep.py` with per-body error report ✅

---

## 3. Property-Based Testing

**Goal:** Use hypothesis-style random testing to find edge cases that
structured tests miss.

### 3.1 Coordinate Transform Round-Trips

- [x] ecliptic → equatorial → ecliptic: verify identity within 0.001" (via swe_cotrans roundtrip)
- [x] ecliptic → XYZ → ecliptic: verify identity
- [x] degrees → radians → degrees: verify identity
- [x] Random JD, random body, random flag combinations

### 3.2 API Contract Testing

- [x] `calc_ut(jd, body, flags)` always returns (6-tuple, int)
- [x] `houses(jd, lat, lon, system)` always returns (tuple, tuple)
- [x] `fixstar(name, jd, flags)` always returns (6-tuple, str, int)
- [x] No function raises unhandled exceptions for valid inputs
- [x] All functions accept SEFLG_MOSEPH without error

### 3.3 Monotonicity Properties

- [x] Sun longitude is monotonically increasing (mod 360) over 1 day
- [x] Mean Node longitude is monotonically decreasing over 1 day
- [x] Julian day conversion is monotonic: julday(y,m,d) < julday(y,m,d+1)
- [x] Delta T in reasonable range (40–120s) for 1972–2040 (note: NOT monotonic — Earth rotation sped up ~2020–2029)

### 3.4 Symmetry Properties

- [x] `calc_ut(jd, body, SEFLG_HELCTR)` for Sun returns (0, 0, 0, ...)
- [x] `pheno_ut(jd, SE_SUN)` returns phase_angle = 0
- [x] `houses(jd, lat, lon)` cusps are in ascending order (mod 360)

**Deliverable:** `tests/test_property_based.py` using hypothesis library ✅

---

## 4. Fuzz Testing

**Goal:** Verify robustness against malformed, extreme, and boundary inputs.

### 4.1 Extreme Julian Dates

- [x] JD = 0 (4713 BC), JD = -1e6, JD = 1e8, and boundary dates
- [x] JD = NaN, Inf, -Inf (all handled gracefully)
- [x] Verify: either valid result or clean exception (no crash/hang)
- [x] Also tested: revjul, julday roundtrips, houses, eclipses, sidtime, deltat

### 4.2 Invalid Body IDs

- [x] Negative body IDs (-1 is SE_ECL_NUT, valid; -100, -999999 raise cleanly)
- [x] Body IDs > 100000 (SE_AST_OFFSET range): clean exceptions
- [x] Body ID = SE_AST_OFFSET + 0 (edge case): clean exception
- [x] All standard bodies (0-11) verified valid; pheno_ut invalid body tested

### 4.3 Extreme Geographic Coordinates

- [x] Latitude = ±90° (poles), ±91° (invalid), ±180°
- [x] Longitude = ±180°, ±360°, ±720°
- [x] Altitude = -1000m, 0m, 8848m, 100000m, 1000000m
- [x] House calculations with 9 systems, poles with Placidus/Equal/Whole Sign

### 4.4 Flag Exhaustion

- [x] All 14 single flags, all 91 flag pairs, 500 sampled combos from 2^14 space
- [x] Higher flags: SEFLG_BARYCTR, SEFLG_TOPOCTR, SEFLG_SIDEREAL, SEFLG_ICRS
- [x] Deliberately conflicting flag combos verified no crash
- [x] Moon-specific flag combos (separate code path) tested

**Deliverable:** `tests/test_fuzz.py` ✅

---

## 5. Eclipse Event Validation

**Goal:** Cross-validate eclipse predictions against NASA eclipse catalogs.

### 5.1 Solar Eclipses (Espenak Canon)

- [x] 20 historical solar eclipses (2001–2020) from NASA Five Millennium Canon
- [x] Compare: type (total/annular/hybrid), maximum time (UT = TD − ΔT)
- [x] Tolerance: < 60 seconds timing, correct type classification (all pass)
- [x] Contact time ordering verified (C1 < C2 < max < C3 < C4)

### 5.2 Lunar Eclipses

- [x] 20 historical lunar eclipses (2001–2022) from NASA catalog
- [x] Compare: type, timing (UT = TD − ΔT), contact ordering
- [x] Verify gamma sign and magnitude (gamma in [-2, 2])
- [x] Note: 2015-Apr-04 excluded (borderline total, umbra mag 1.0008)

### 5.3 Future Eclipses

- [x] Next 10 solar eclipses (2026–2029): date and type verified
- [x] Next 10 lunar eclipses (2025–2029): date and type verified
- [x] Note: 2027-Jul-18 penumbral excluded (mag 0.0014, smallest of century)

**Deliverable:** `tests/test_eclipse_catalog_validation.py` ✅

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
