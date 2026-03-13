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

### A1. J2000 Frame — systematic ~0.004° offset (nodes/Lilith)

- **Symptom:** ~14-22" systematic offset when `SEFLG_J2000` is used
- **Ground truth:** JPL Horizons "ecliptic J2000" output + astropy ICRS-to-ecliptic-J2000
- **Hypothesis:** Different precession model in frame transformation (IAU 2006 vs older)
- **Test:** Compute True Node / Lilith in J2000 ecliptic for 10 dates, compare with Horizons
- **If we win:** Document "IAU 2006 precession produces more accurate J2000 frame transformation"
- **Priority:** HIGH (potential selling point)

### A2. Mean Lilith (MEAN_APOG) — 0.01° (~36")

- **Symptom:** Consistent ~15" mean difference in mean Lilith longitude
- **Ground truth:** Chapront et al. (2002) published mean orbital elements of lunar perigee
- **Hypothesis:** Different coefficients in the analytical mean argument of perigee formula
- **Test:** Compare formula coefficients against the original publication
- **If we win:** Document "Uses updated Chapront mean elements"
- **Priority:** LOW (small difference, analytical formula)

### A3. Equatorial Coordinates — TRUE_NODE 0.2°, OSCU_APOG 0.15°

- **Symptom:** Large tolerances in ecliptic-to-equatorial transformation
- **Ground truth:** Independent calculation with pyerfa `erfa.ecl2equ()`
- **Hypothesis:** Different obliquity value (true vs mean, nutation model)
- **Test:** Compute RA/Dec independently with pyerfa and compare
- **If we win:** Document "Equatorial transformation uses IAU 2006/2000A obliquity via ERFA"
- **Priority:** HIGH (these are suspiciously large)

### A4. Latitude nodes/Lilith — TRUE_NODE 0.5°, OSCU_APOG 1.5°

- **Symptom:** Very large differences in latitude
- **Ground truth:** JPL Horizons ecliptic latitude
- **Hypothesis:** Could be a real bug in the eccentricity vector latitude calculation
- **Test:** Compare with Horizons for 20 dates
- **If it's our bug:** Fix it
- **If we win:** Document
- **Priority:** HIGHEST (1.5° could be a real bug)

### A5. Velocity nodes/Lilith — TRUE_NODE 0.2°/day, OSCU_APOG 1.0°/day

- **Symptom:** Large tolerance in velocity values
- **Ground truth:** Numerical differentiation of Horizons positions (t-0.5 and t+0.5)
- **Test:** 10 dates, compare numerical velocities vs ours vs SE
- **Priority:** MEDIUM (velocities are secondary to positions)

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

### A7. Star-based Sidereal Ayanamshas — 0.1°

- **Depends on:** A6 results. If our star positions are more precise, ayanamsha is too
- **Test:** Compute Lahiri ayanamsha (Spica-based) independently
- **Priority:** LOW

### A8. Heliocentric/Barycentric Planets — 0.03°

- **Symptom:** Relaxed tolerance in heliocentric/barycentric frames
- **Ground truth:** JPL Horizons heliocentric output
- **Test:** 5 planets x 5 dates
- **Priority:** LOW

**Execution order: A4 > A1 > A3 > A5 > A6 > A2 > A8 > A7**

---

## BLOCK B — Fix `fixed_stars.py` (Silently Ignored Flags)

### B1. Problem

`fixed_stars.py` only handles 4 of 19 SEFLG flags (`SPEED`, `NOABERR`, `NOGDEFL`,
`EQUATORIAL`). All others are **silently ignored**. Most critical: `SEFLG_SIDEREAL`
returns tropical coordinates without any warning.

### B2. Implementation

Apply to all 4 public functions (`swe_fixstar_ut`, `swe_fixstar`, `swe_fixstar2_ut`,
`swe_fixstar2`). Extract shared logic into a helper.

| Flag | Priority | Implementation |
|------|----------|----------------|
| `SEFLG_SIDEREAL` | **Critical** | `lon = (lon - ayanamsha) % 360` after calc, like planets.py |
| `SEFLG_J2000` | High | Use `ecliptic_J2000_frame` instead of `ecliptic_frame` |
| `SEFLG_NONUT` | High | Use mean ecliptic frame (no nutation) |
| `SEFLG_XYZ` | Medium | Spherical-to-Cartesian (reuse pattern from planets.py) |
| `SEFLG_RADIANS` | Medium | Degrees-to-radians (reuse pattern from planets.py) |
| `SEFLG_ICRS` | Medium | ICRS frame without frame bias |
| `SEFLG_TRUEPOS` | Low | Geometric position (use `astrometric` directly) |
| `SEFLG_MOSEPH` | Low | Silent strip like planets.py |
| `SEFLG_SPEED3` | Low | Convert to `SEFLG_SPEED` |
| `SEFLG_TOPOCTR` | Low | Accept silently (stars at infinite distance) |

### B3. Tests

Add comparison tests in `test_compare_fixedstars.py`:
- `SEFLG_SIDEREAL` with at least 3 ayanamshas
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

| Step | Block | Description |
|------|-------|-------------|
| 1 | A4 | Latitude cross-check (possible bug, highest priority) |
| 2 | A1 | J2000 frame cross-check (selling point) |
| 3 | A3 | Equatorial coordinates cross-check |
| 4 | A5 | Velocity cross-check |
| 5 | A6 | Fixed stars cross-check |
| 6 | B1/B2/B3 | Fix fixed_stars.py flags + tests |
| 7 | A2, A7, A8 | Remaining cross-checks |
| 8 | C1-C4 | Documentation updates |
| 9 | Commit | Final commit on dev/precision-v3 |
