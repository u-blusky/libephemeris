# Precision V3 — Final Plan

## Context

LibEphemeris is a pure-Python astronomical ephemeris library using NASA JPL DE440
via Skyfield. It aims for 1:1 API compatibility with PySwissEphemeris while being
built on modern IAU standards and JPL data.

We have completed a deep, exhaustive precision audit across 9 sessions with 15
commits on `dev/precision-v3`. All 4 phases of the original plan (Critical,
High Priority, Keplerian, Low Priority) are complete.

### What's already been verified (Session 9 triangulation)

**True Node (longitude):**
- 3-source comparison: libephemeris vs pyswisseph vs JPL Horizons
- We match JPL Horizons to **<0.01"** across 24 test dates (1950-2050)
- We match SE to **<0.01"** (machine precision) in 200-date dense comparison
- Old tolerance was 0.15° — tightened to 0.001° (commit `284bed7`)

**True Lilith / Osculating Apogee (longitude):**
- Both us AND SE differ from Horizons by ~240" (~4 arcmin)
- This is NOT a bug — it's the two-body approximation in the eccentricity vector formula
- We match SE to **mean 0.1", max 0.5"** in 200-date dense comparison
- Old tolerance was 0.1° — tightened to 0.001° (commit `284bed7`)

### What remains unverified

Several areas still have tolerances > 0.001° where we haven't determined who is
correct (us or SE). These need triangulated cross-checks against independent
sources before we can document them as advantages or fix them as bugs.

---

## BLOCK A — Cross-Checks (Who Is Right?)

For every area with tolerance > 0.001°, we perform a 3-source triangulation:
**us vs SE vs independent source** (JPL Horizons, astropy/ERFA, Gaia/SIMBAD).

### A1. J2000 Frame — systematic ~0.004° offset (nodes/Lilith) ✅ COMPLETED — WE ARE MORE ACCURATE

- **Symptom:** ~14-22" systematic offset when `SEFLG_J2000` is used
- **Ground truth:** JPL Horizons "ecliptic J2000" output + astropy ICRS-to-ecliptic-J2000
- **Results:**
  - SE applies a spurious ~14" offset at J2000 epoch that should not exist
  - Our IAU 2006 precession produces the correct J2000 frame transformation
  - Verified independently against JPL Horizons and astropy/ERFA
- **Conclusion:** Documented as a selling point — our J2000 frame is more accurate than SE's

### A2. Mean Lilith (MEAN_APOG) — 0.01° (~36") ✅ COMPLETED — NOT a bug

- **Symptom:** Consistent ~15" mean difference in mean Lilith longitude
- **Ground truth:** Chapront et al. (2002) published mean orbital elements of lunar perigee
- **Results:**
  - Longitude: < 0.015" within 1950–2050 (essentially zero). Grows to ~1.9" at year 2600 due to minor higher-order polynomial coefficient difference (~0.05–0.07"/cy²). Negligible.
  - Latitude: Systematic ~15–20" sinusoidal difference. Traced to phase offset in the node calculation used for orbital-plane-to-ecliptic projection (`lat = i × sin(apogee − node)`). Both libraries use slightly different analytical node formulas.
  - Speed: < 0.00005°/day — perfect agreement.
- **Conclusion:** No action needed. Differences are inherent to the different analytical approaches (DE404-fitted Chapront coefficients vs SE internal coefficients). Mean Lilith is a mathematical point, not a physical body — sub-arcsecond longitude agreement is excellent.

### A3. Equatorial Coordinates — TRUE_NODE 0.2°, OSCU_APOG 0.15° ✅ COMPLETED — <0.2"

- **Symptom:** Large tolerances in ecliptic-to-equatorial transformation
- **Ground truth:** Independent calculation with pyerfa `erfa.ecl2equ()`
- **Results:**
  - Tolerances tightened from 0.2° to <0.2" (<0.001°) — sub-arcsecond agreement
  - Verified independently against pyerfa/ERFA obliquity calculations
  - Old tolerances were massively over-conservative
- **Conclusion:** No bug. Tightened in commit `f852000`

### A4. Latitude nodes/Lilith — TRUE_NODE 0.5°, OSCU_APOG 1.5° ✅ COMPLETED — <0.02"

- **Symptom:** Very large differences in latitude
- **Ground truth:** JPL Horizons ecliptic latitude
- **Results:**
  - Tolerances tightened from 1.5° to <0.02" (<0.001°) — sub-arcsecond agreement
  - Old tolerances were massively over-conservative (by ~5400x)
  - No bug in eccentricity vector latitude calculation
- **Conclusion:** No bug. Tightened in commit `f852000`

### A5. Velocity nodes/Lilith — TRUE_NODE 0.2°/day, OSCU_APOG 1.0°/day ✅ COMPLETED — <0.015°/day

- **Symptom:** Large tolerance in velocity values
- **Ground truth:** Numerical differentiation of Horizons positions (t-0.5 and t+0.5)
- **Results:**
  - Tolerances tightened from 1.0°/day to <0.015°/day
  - Sub-arcsecond velocity agreement for most bodies
  - True Lilith velocity tolerance <0.05°/day (slightly larger due to osculating element sensitivity)
- **Conclusion:** No bug. Tightened in commit `f852000`

### A6. Fixed Stars — 0.002° (7.2") ✅ COMPLETED

- **Symptom:** Position differences for fixed stars
- **Ground truth:** van Leeuwen 2007 (new Hipparcos reduction, I/311/hip2) via VizieR TAP + SIMBAD
- **Result:** IMPROVED — 99 stars updated to van Leeuwen 2007 proper motions, 2 catalog bugs fixed
- **Findings:**
  - Proper motions updated from original Hipparcos 1997 to van Leeuwen 2007 for 99/116 stars
  - Algedi: RA/Dec were for wrong component (Alpha-1 instead of Alpha-2 Cap) — fixed
  - Asellus Borealis: HIP number was wrong (43103=Iota Cnc, corrected to 42806=Gamma Cnc) — fixed
  - 5 stars resolve to different physical components vs SE (Menkar, Algedi, Algieba, Albireo, Almach)
  - Independent astropy verification: both us and SE agree with astropy to <0.5" for all 10 principal stars
  - Sirius 0.24" and Fomalhaut 0.11" offsets traced to annual parallax (not modeled)
  - Post-update: 100% of 101 comparable stars within 0.002° of SE, 98% within 0.5"

### A7. Star-based Sidereal Ayanamshas — 0.1° ✅ COMPLETED — up to 19", expected

- **Depends on:** A6 results. If our star positions are more precise, ayanamsha is too
- **Results:**
  - Up to 19" difference for star-dependent ayanamsha modes (Lahiri, Raman, etc.)
  - Differences are expected and directly traceable to different star catalog positions (van Leeuwen 2007 vs original Hipparcos 1997)
  - Formula-based ayanamshas (not star-dependent) agree to <0.001°
- **Conclusion:** No action needed. Star-based ayanamshas inherit the star position differences from A6.

### A8. Heliocentric/Barycentric Planets — 0.03° ✅ COMPLETED

- **Symptom:** Relaxed tolerance in heliocentric/barycentric frames
- **Ground truth:** JPL Horizons heliocentric output, raw Skyfield/DE440
- **Test:** 10 bodies × 10 dates × 4 modes (helio, bary, equatorial, XYZ) = 503 tests
- **Results:**
  - Heliocentric (Mercury–Pluto): < 0.0004° (1.1") — sub-arcsecond, no issues
  - Barycentric (Moon–Pluto): < 0.001° — sub-arcsecond, no issues
  - Barycentric Sun: up to 0.04° angular, but actual 3D offset is only ~120 km (0.017% of solar radius). Angular amplification from tiny Sun-SSB distance (~0.001–0.009 AU). Verified LE is closer to raw Skyfield/DE440 at J2000.
  - Equatorial: < 0.0005° (1.7") — sub-arcsecond, no issues
  - XYZ Cartesian: < 0.00005 AU — sub-arcsecond angular at all distances
- **Bug found:** SEFLG_XYZ and SEFLG_RADIANS not preserved in retflag — fixed in planets.py
- **Test file:** `compare_scripts/tests/test_compare_helio_bary.py` (503 tests, all pass)

**Execution order: A4 > A1 > A3 > A5 > A6 > A2 > A8 > A7**

---

## BLOCK B — Fix `fixed_stars.py` (Silently Ignored Flags) ✅ COMPLETED

### B1. Problem (FIXED in commit `f852000`)

`fixed_stars.py` originally only handled 4 of 19 SEFLG flags (`SPEED`, `NOABERR`, `NOGDEFL`,
`EQUATORIAL`). All others were **silently ignored**. Most critical: `SEFLG_SIDEREAL`
returned tropical coordinates without any warning.

**9 new flags implemented:** SIDEREAL, J2000, NONUT, XYZ, RADIANS, TRUEPOS, MOSEPH, SPEED3, TOPOCTR.

### B2. Implementation (DONE)

Applied to all 4 public functions (`swe_fixstar_ut`, `swe_fixstar`, `swe_fixstar2_ut`,
`swe_fixstar2`). Shared logic extracted into helpers.

| Flag | Status | Implementation |
|------|--------|----------------|
| `SEFLG_SIDEREAL` | ✅ Done | `lon = (lon - ayanamsha) % 360` after calc, like planets.py |
| `SEFLG_J2000` | ✅ Done | Use `ecliptic_J2000_frame` instead of `ecliptic_frame` |
| `SEFLG_NONUT` | ✅ Done | Use mean ecliptic frame (no nutation) |
| `SEFLG_XYZ` | ✅ Done | Spherical-to-Cartesian (reuse pattern from planets.py) |
| `SEFLG_RADIANS` | ✅ Done | Degrees-to-radians (reuse pattern from planets.py) |
| `SEFLG_ICRS` | N/A | Already handled before this work |
| `SEFLG_TRUEPOS` | ✅ Done | Geometric position (use `astrometric` directly) |
| `SEFLG_MOSEPH` | ✅ Done | Silent strip like planets.py |
| `SEFLG_SPEED3` | ✅ Done | Convert to `SEFLG_SPEED` |
| `SEFLG_TOPOCTR` | ✅ Done | Accept silently (stars at infinite distance) |

### B3. Tests (DONE)

Comparison tests added in `test_compare_fixedstars.py`:
- `SEFLG_SIDEREAL` with multiple ayanamshas
- `SEFLG_J2000`
- `SEFLG_NONUT`
- `SEFLG_XYZ` + `SEFLG_RADIANS`

---

## BLOCK C — Documentation of Superiority

### C1. Update `docs/reference/swisseph-comparison.md`

Add section **"7. Independent Verification Results"**:
- True Node triangulation results (match Horizons to <0.01")
- True Lilith triangulation results
- Results from A1-A8 cross-checks
- For each area: table "LibEphemeris vs SE vs Independent Source"

### C2. Update `docs/reference/precision.md`

For each section where we are more precise:
- Add independent source used for verification
- Add concrete numbers (our deviations vs Horizons in arcseconds)

### C3. Update README.md — "Why" section

Add bullet points with **verifiable claims**:
- "Independently verified precision against JPL Horizons"
- Concrete numbers for lunar nodes, Lilith
- Any additional wins from A1-A8

Tone: **factual, verifiable, with references**. Never "we're better" — instead
"here are the numbers, verified against independent sources."

### C4. Update investigation plan

Update `plans/precision-investigation-true-node-lilith.md` with results.

---

## Execution Order

| Step | Block | Description | Status |
|------|-------|-------------|--------|
| 1 | A4 | Latitude cross-check (possible bug, highest priority) | ✅ Done (Session 10) — <0.02", tightened |
| 2 | A1 | J2000 frame cross-check (selling point) | ✅ Done (Session 10) — WE ARE MORE ACCURATE |
| 3 | A3 | Equatorial coordinates cross-check | ✅ Done (Session 10) — <0.2", tightened |
| 4 | A5 | Velocity cross-check | ✅ Done (Session 10) — <0.015°/day, tightened |
| 5 | A6 | Fixed stars cross-check | ✅ Done (Session 11) — 99 stars updated, 2 bugs fixed |
| 6 | B1/B2/B3 | Fix fixed_stars.py flags + tests | ✅ Done (Session 10) — 9 flags implemented |
| 7 | A2, A7, A8 | Remaining cross-checks | ✅ Done (Sessions 11-12) — all verified |
| 8 | C1-C4 | Documentation updates | ✅ Partially done — swisseph-comparison.md updated |
| 9 | Commit | Final commit on dev/precision-v3 | ✅ All merged to main (5 commits ahead of origin) |
