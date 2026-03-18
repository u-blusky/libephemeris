# Extended Tier LEB — Validation Plan V2

## Overview

This is the successor to `EXTENDED_VALIDATION_PLAN.md` (V1), which was fully executed across 8 phases / 38 rounds. V2 is designed to be **100× more thorough**: it targets every gap found during V1 execution, adds systematic coverage of all body×flag×date×ayanamsha combinations, introduces automated regression detection, and defines reproducible acceptance criteria with measured baselines.

### V1 Results Summary (Baseline)

| Metric | V1 Result |
|--------|-----------|
| Phases completed | 8/8 |
| Rounds completed | 38/38 |
| Bugs found & fixed | 4 sidereal bugs |
| Bodies validated | 31 |
| Date range | -5000 to +5000 CE |
| Position precision (LEB vs Skyfield) | Moon <0.00034", planets <0.000016" |
| Position headroom | 60× to 5,000,000× |

### V2 Goals

1. **Exhaustive combinatorial coverage**: Every body × every flag combo × every ayanamsha mode × date strata
2. **Regression-proof**: Every bug fix has a dedicated regression test with tight tolerance
3. **Statistical rigor**: P50/P95/P99/max for every dimension, not just spot checks
4. **Cross-pipeline**: Systematic testing of all 3 pipelines (ICRS, Ecliptic, Heliocentric)
5. **Cross-tier**: Extended vs Medium vs Base in overlap zones
6. **Temporal stratification**: Ancient (-5000 to -1000), Classical (-1000 to 500), Medieval (500-1500), Modern (1500-2500), Future (2500-5000)
7. **Three-way everywhere**: LEB ↔ Skyfield ↔ pyswisseph for all testable scenarios
8. **Automated**: All rounds produce machine-readable output (JSON/CSV) for trend tracking

---

## Phase 1: Exhaustive LEB vs Skyfield Precision Profiling (12 rounds)

### Round 1.1 — Grand Sweep: 31 Bodies × 1000 Dates × 6 Components

- **Bodies**: All 31 LEB bodies
- **Dates**: 1000 uniformly distributed across full range (-5000 to +5000 CE)
- **Flags**: `SEFLG_SPEED` (ecliptic of date)
- **Components**: lon, lat, dist, dlon, dlat, ddist
- **Metrics**: Per body per component: N, mean, median, P50, P75, P90, P95, P99, P99.9, max, RMS, std
- **Stratification**: Report separately for 5 temporal strata
- **Output**: CSV with columns: body, component, stratum, metric, value
- **Acceptance**: All positions within TOLS_EXT thresholds with ≥10× headroom for planets

### Round 1.2 — Moon Intensive: 5000 Dates, Segment-Aware

- **Body**: Moon only
- **Dates**: 5000 dates uniformly distributed across full range
- **Additional**: 500 dates at segment boundaries (±0.001 day from boundary)
- **Components**: All 6
- **Metrics**: Same as 1.1 plus segment-boundary specific statistics
- **Stratification**: By millennium, by segment size (should all be 4-day)
- **Acceptance**: Position <0.0005" P99, speed <0.001 d/d P99

### Round 1.3 — Segment Boundary Continuity: All Bodies × 500 Boundaries

- **Bodies**: All 31
- **Method**: For each body, sample 500 segment boundaries uniformly. At each boundary, evaluate at (edge - 0.1s), (edge), (edge + 0.1s). Measure:
  - Position jump: |pos(edge+ε) - pos(edge-ε)| vs linear prediction from speed
  - Speed jump: |speed(edge+ε) - speed(edge-ε)|
- **Acceptance**: Position jump <0.0001", speed jump <0.01 d/d for all bodies except OscuApog

### Round 1.4 — Segment Interior Precision Profile

- **Bodies**: Moon, Mercury, Venus, Mars, Jupiter, Saturn, TrueNode, OscuApog, Chiron, Ceres
- **Method**: For 200 randomly chosen segments per body, evaluate at 11 Chebyshev-distributed points within each segment (0%, 5%, 15%, 25%, 35%, 50%, 65%, 75%, 85%, 95%, 100%)
- **Metrics**: Error distribution at each relative position
- **Acceptance**: Edge error should not exceed 3× midpoint error

### Round 1.5 — Velocity: All 3 Components × All Bodies × 500 Dates

- **Bodies**: All 31
- **Dates**: 500 uniformly distributed
- **Metrics**: Per body per speed component: full statistical profile
- **Additional**: Identify any body/date with max/P99 ratio >50 (indicates isolated spike)
- **Spike analysis**: For each spike found, determine if it's a Skyfield finite-difference artifact or real LEB issue (check if position is also affected)
- **Acceptance**: lon_speed <0.05 d/d, lat_speed <0.005 d/d, dist_speed <1.2e-4 AU/d

### Round 1.6 — Edge-of-Range Stress Test: First/Last 200 Years

- **Bodies**: Sun, Moon, Mercury, Mars, Jupiter, Pluto, MeanNode, TrueNode, Chiron, Cupido
- **Dates**: 100 dates in JD [-105152.5, -32322.5] (first 200 yr) + 100 in JD [3474442.5, 3547272.5] (last 200 yr)
- **Comparison**: Error magnitude at edges vs range center (from Round 1.1)
- **Acceptance**: ≤5× degradation vs center for all bodies

### Round 1.7 — Distance Precision: Geocentric + Heliocentric + Barycentric

- **Geocentric**: All 31 bodies, 200 dates
- **Heliocentric**: Mercury–Pluto, Chiron, Ceres–Vesta (SEFLG_HELCTR), 200 dates
- **Barycentric**: Mercury–Pluto (SEFLG_BARYCTR), 100 dates
- **Metrics**: |dDist| in AU, relative error |dDist|/|dist|, per body
- **Acceptance**: Geocentric <5e-6 AU, heliocentric <5e-5 AU (inner planet amplification), barycentric <5e-6 AU

### Round 1.8 — Acceleration Consistency: 200 Dates × 10 Bodies

- **Bodies**: Moon, Mercury, Venus, Mars, Jupiter, Saturn, TrueNode, OscuApog, Chiron, Ceres
- **Method**: At each date, compute speed at t-0.5d, t, t+0.5d. Compute acceleration. Compare LEB vs Skyfield.
- **Acceptance**: Acceleration difference <0.001 d/d² for planets, <0.01 d/d² for Moon

### Round 1.9 — Sub-Arcsecond Precision Map: 10,000 Dates × Key Bodies

- **Bodies**: Sun, Moon, Mars, Jupiter, Saturn
- **Dates**: 10,000 uniformly distributed across full range
- **Metrics**: Create 2D heatmap of |error| vs JD for each body
- **Purpose**: Identify if there are any localized pockets of degraded precision
- **Output**: PNG heatmaps + CSV data
- **Acceptance**: No localized pockets exceeding 10× the global P99

### Round 1.10 — Ecliptic Body Deep Dive: MeanNode/MeanApog Polynomial Degradation

- **Bodies**: MeanNode (10), MeanApog (12)
- **Dates**: 2000 dates across full range, concentrated at extremes (±40 centuries from J2000)
- **Focus**: Meeus polynomial degradation quantification
- **Metrics**: Error vs |centuries from J2000| regression analysis
- **Acceptance**: Degradation rate documented, <0.01" within ±20 centuries

### Round 1.11 — TrueNode/OscuApog Oscillation Fidelity

- **Bodies**: TrueNode (11), OscuApog (13)
- **Dates**: 5000 dates at 1-day intervals for 2000-2015 CE
- **Method**: Compare oscillation pattern (TrueNode around MeanNode, OscuApog around MeanApog)
- **Metrics**: Oscillation amplitude, frequency, phase accuracy LEB vs Skyfield
- **Acceptance**: Amplitude error <0.001", no phase drift >0.01 day over 15 years

### Round 1.12 — IntpApog/IntpPerg: Full SPK Coverage

- **Bodies**: IntpApog (21), IntpPerg (22)
- **Dates**: 500 dates within SE coverage (~-3000 to +2900 CE)
- **Focus**: All 6 components at high density
- **Additional**: Verify behavior at SE coverage boundaries (graceful fallback)
- **Acceptance**: Within established tolerances

---

## Phase 2: Three-Way LEB ↔ Skyfield ↔ SE (12 rounds)

### Round 2.1 — Pipeline A Three-Way: Modern Era Dense (1800-2200)

- **Bodies**: All 11 ICRS planets (Sun–Pluto + Earth)
- **Dates**: 200 dates from 1800-2200 CE
- **Flags**: `SEFLG_SPEED`
- **Metrics**: |LEB-SE|, |Sky-SE|, excess = |LEB-SE| - |Sky-SE|, for all 6 components
- **Acceptance**: excess <0.002" for position, <0.001 d/d for speed

### Round 2.2 — Pipeline A Three-Way: Historical/Future (-3000 to +3000)

- **Bodies**: All 11 ICRS planets
- **Dates**: 100 dates from -3000 to +3000 CE
- **Metrics**: Same as 2.1
- **Acceptance**: excess <0.01" for planets, <0.1" for Moon

### Round 2.3 — Pipeline A Three-Way: Extreme Era (-5000 to -3000 and +3000 to +5000)

- **Bodies**: Sun, Moon, Mars, Jupiter, Saturn, Pluto (all available in SE)
- **Dates**: 50 dates in each extreme era
- **Focus**: Verify LEB doesn't add error beyond the DE431↔DE441 baseline
- **Acceptance**: excess <0.1" (baseline is 20-300" at extremes)

### Round 2.4 — Pipeline B Three-Way: Nodes and Apogees

- **Bodies**: MeanNode, TrueNode, MeanApog, OscuApog, IntpApog, IntpPerg
- **Dates**: 100 dates from -3000 to +2900 CE
- **Metrics**: Triple comparison for all 6 components
- **Acceptance**: |LEB-Sky| <0.001" for Mean, <5" for True, <100" for Oscu

### Round 2.5 — Pipeline C Three-Way: Uranians and Transpluto

- **Bodies**: Cupido(40), Hades(41), Zeus(42), Kronos(43), Apollon(44), Admetos(45), Vulkanus(46), Poseidon(47), Transpluto(48)
- **Dates**: 50 dates from -3000 to +3000 CE
- **Metrics**: Triple comparison
- **Acceptance**: |LEB-Sky| <0.001", excess <0.001"

### Round 2.6 — Asteroids Three-Way: SPK Coverage Zone

- **Bodies**: Chiron(15), Ceres(17), Pallas(18), Juno(19), Vesta(20)
- **Dates**: 100 dates within SPK coverage (1900-2100 CE)
- **Metrics**: Triple comparison for position + speed
- **Acceptance**: |LEB-Sky| <0.001", excess <0.01"

### Round 2.7 — Sidereal Three-Way: Modern Era (all flag combos)

- **Bodies**: Sun, Moon, Mars, Jupiter, MeanNode, TrueNode, OscuApog, Cupido
- **Dates**: 20 dates from 1800-2200 CE
- **Flags**: 8 sidereal combinations:
  - SIDEREAL, SID+EQ, SID+J2K, SID+EQ+J2K
  - (same 4 without SIDEREAL as control)
- **Ayanamsha**: Lahiri(1), Fagan-Bradley(0), Raman(3)
- **Metrics**: Triple comparison. Compute:
  - excess for each sidereal flag combo
  - |LEB_sid - LEB_trop| - |SE_sid - SE_trop| (sidereal offset consistency)
- **Acceptance**: excess <0.002" for all bodies in modern era

### Round 2.8 — Sidereal Three-Way: Extended Era

- **Bodies**: Sun, Moon, Mars, MeanNode, TrueNode, Cupido
- **Dates**: 20 dates from -3000 to +3000 CE
- **Flags**: 4 sidereal combinations + 4 controls
- **Ayanamsha**: Lahiri(1)
- **Acceptance**: excess <0.05"

### Round 2.9 — Houses Three-Way

- **House systems**: Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign, Porphyrius, Alcabitius, Morinus, Topocentric, Azimuthal/Horizontal
- **Locations**: Rome(41.9°N), NYC(40.7°N), Tokyo(35.7°N), Sydney(-33.9°S), Reykjavik(64.1°N), Equator(0°), Tromso(69.6°N), Cape Town(-33.9°S), Buenos Aires(-34.6°S), High altitude(27.99°N, 8848m)
- **Dates**: 20 dates from 1800-2200 CE
- **Modes**: Tropical + Sidereal(Lahiri)
- **Metrics**: max |cusp_diff| across 12 cusps + ASC + MC
- **Acceptance**: LEB vs Skyfield <0.001", LEB vs SE <0.01"

### Round 2.10 — Crossing Functions Three-Way

- **Sun crossings**: 0°, 30°, 60°, 90°, 120°, 150°, 180°, 210°, 240°, 270°, 300°, 330° (all 12 sign ingresses) from 1900-2100 CE
- **Moon crossings**: 0°, 90°, 180°, 270° from 2000-2025 CE
- **Node crossings**: Moon node crossings from 2000-2025 CE
- **Helio crossings**: Mars, Jupiter 0° from 2000-2025 CE
- **Metrics**: |dt| in seconds between each pair
- **Acceptance**: LEB vs Skyfield <0.1s, LEB vs SE <0.5s

### Round 2.11 — Eclipse Functions Three-Way

- **Solar eclipses**: Find 30 eclipses (2000-2030) via all three systems
- **Lunar eclipses**: Find 30 eclipses (2000-2030) via all three systems
- **Metrics**: Maximum time difference, eclipse type classification agreement, magnitude difference
- **Acceptance**: Time <1s, type classification 100% match, magnitude <0.001

### Round 2.12 — Rise/Transit/Set Three-Way

- **Bodies**: Sun, Moon, Mars, Jupiter, Venus
- **Locations**: Rome, NYC, Tokyo, Sydney, Equator
- **Dates**: 20 dates from 2000-2025 CE
- **Events**: Rise, Transit, Set
- **Metrics**: |dt| in seconds
- **Acceptance**: <1s for Sun, <5s for Moon, <30s for planets

---

## Phase 3: Sidereal Bug Fix Regression Suite (8 rounds)

This phase exists solely to guard against regressions of the 4 sidereal bugs fixed in the leb/precision branch. Each round targets the specific error signature of one bug.

### Round 3.1 — Bug 1 Regression: Pipeline A SID+EQ (Mean Equator)

- **Bug**: Used nutation matrix instead of mean equator precession for SID+EQ
- **Error signature**: ~36" RA error for SID+EQ, ~0.3" for SID+J2K
- **Bodies**: Sun, Moon, Mars, Jupiter, Saturn (Pipeline A)
- **Dates**: 20 dates from 1800-2200 CE + 10 extreme dates
- **Flags**: SID+EQ, SID+J2K, SID+EQ+J2K
- **Ayanamsha**: Lahiri(1), Fagan-Bradley(0), Raman(3)
- **Reference**: pyswisseph
- **Tolerance**: <0.01° (catches ~36" error)
- **Three-way**: Also compare LEB vs Skyfield (<0.001")

### Round 3.2 — Bug 2 Regression: Pipeline B/C SID+EQ (dpsi Handling)

- **Bug**: Wrong dpsi nutation handling for ecliptic bodies in SID+EQ
- **Error signature**: ~10-20" for MeanNode SID+EQ
- **Bodies**: MeanNode, TrueNode, MeanApog, OscuApog, IntpApog, IntpPerg (all Pipeline B)
- **Dates**: 20 dates from 1800-2200 CE
- **Flags**: SID+EQ, SID+EQ+J2K
- **Ayanamsha**: Lahiri(1), Fagan-Bradley(0), Raman(3)
- **Reference**: pyswisseph
- **Tolerance**: <0.015° for primary, <5.5° for IntpPerg (known algorithm diff)
- **Additional**: Verify MeanNode/MeanApog skip dpsi, TrueNode/OscuApog subtract dpsi

### Round 3.3 — Bug 3 Regression: J2000 Suppression for Non-Mean Bodies

- **Bug**: J2K flag was applied to TrueNode/OscuApog/IntpApog/IntpPerg when sidereal
- **Error signature**: SID+J2K gives different values than SID for suppressed bodies
- **Bodies**: All 6 Pipeline B bodies
- **Dates**: 20 dates from 1800-2200 CE
- **Tests**:
  1. For TrueNode, OscuApog, IntpApog, IntpPerg: verify `|SID - SID+J2K| < 0.0001°`
  2. For MeanNode, MeanApog: verify `|SID - SID+J2K| > 0.001°` (J2K IS applied)
  3. Same checks on pyswisseph (reference), Skyfield, and LEB
- **Reference**: pyswisseph behavior
- **Acceptance**: 100% behavioral match with pyswisseph

### Round 3.4 — Bug 4a Regression: Frame Bias in Precession Matrix

- **Bug**: `_get_precession_matrix()` used `t.P` including ICRS frame bias (~17 mas)
- **Error signature**: ~17 mas systematic at J2000, growing with time
- **Bodies**: Sun, Moon, Mars, Jupiter (Pipeline A)
- **Dates**: 20 dates spanning 500 years (tight tolerance catches 17 mas)
- **Flags**: SID+EQ
- **Method**: Compare LEB SID+EQ vs Skyfield SID+EQ
- **Tolerance**: <0.05" (3× the 17 mas error)
- **Additional**: Verify error does NOT grow systematically with |T - J2000|

### Round 3.5 — Bug 4b Regression: SID+J2K Precession Order (MeanNode/MeanApog)

- **Bug**: Precession applied before ayanamsha subtraction (non-commutative, up to 28" at extreme dates)
- **Error signature**: Error grows with distance from J2000
- **Bodies**: MeanNode, MeanApog
- **Dates**: 30 dates from -3000 to +3000 CE (stress non-commutativity)
- **Flags**: SID+J2K, SID+EQ+J2K
- **Reference**: pyswisseph
- **Tolerance**: <0.015° (catches 28" error)
- **Regression indicator**: Plot error vs |centuries from J2000| — should be flat, not growing

### Round 3.6 — Cross-Mode Sidereal Consistency

- **Purpose**: Verify sidereal corrections are identical in LEB and Skyfield paths
- **Bodies**: Sun, Moon, Mars, MeanNode, TrueNode, Cupido
- **Dates**: 20 dates from 1800-2200 CE
- **Method**: For each body/date/ayanamsha:
  1. Compute tropical lon via LEB and Skyfield
  2. Compute sidereal lon via LEB and Skyfield
  3. Compute ayanamsha via LEB and Skyfield
  4. Verify: (tropical - sidereal) mod 360 ≈ ayanamsha
  5. Verify: |LEB_ayanamsha - Sky_ayanamsha| < 1e-10°
- **Ayanamsha**: All 29 formula-based modes
- **Acceptance**: <1e-8° for tropical-sidereal offset consistency

### Round 3.7 — Sidereal Speed Consistency

- **Purpose**: Verify sidereal speed correction is identical in LEB and Skyfield
- **Bodies**: Sun, Moon, Mars, Jupiter
- **Dates**: 50 dates from 1800-2200 CE
- **Method**: For each:
  1. Compute tropical speed and sidereal speed via LEB and Skyfield
  2. Speed difference should equal precession rate (~0.0000382 d/d)
  3. LEB and Skyfield should agree on sidereal speed within 1e-5 d/d
- **Acceptance**: Speed offset consistency <1e-5 d/d

### Round 3.8 — All Bug Fixes Combined: 4-Way Stress Test

- **Purpose**: Apply all sidereal flag combos simultaneously across all pipelines
- **Bodies**: Sun(A), Moon(A), Mars(A), MeanNode(B), TrueNode(B), OscuApog(B), MeanApog(B), IntpApog(B), Cupido(C)
- **Dates**: 10 dates (5 modern, 5 extended era)
- **Flags**: All 16 combinations of SIDEREAL × EQUATORIAL × J2000 × (3 ayanamsha modes)
- **Reference**: pyswisseph + Skyfield (both)
- **Acceptance**: All match within established tolerances

---

## Phase 4: Coordinate System & Flag Edge Cases (8 rounds)

### Round 4.1 — 0°/360° Longitude Wrap-Around: All Bodies

- **Bodies**: All 31
- **Method**: For each body, find 20 dates where longitude is near 0° (within 1°). Evaluate at 100 points (6-minute intervals) around each crossing.
- **Acceptance**: No jumps >0.001° between consecutive evaluations

### Round 4.2 — Latitude Extrema and Sign Changes

- **Bodies**: Moon(±5.3°), Mercury(±7°), Venus(±3.4°), Mars(±1.8°), Pluto(±17°), Pallas(high inclination)
- **Method**: For each body, find 20 zero-crossings and 20 extrema. Evaluate LEB vs Skyfield at 50 points around each.
- **Acceptance**: |dLat| <0.001" everywhere

### Round 4.3 — Ecliptic↔Equatorial Round-Trip

- **Bodies**: 15 bodies (Sun, Moon, Mercury–Saturn, Earth, MeanNode, TrueNode, Chiron, Cupido)
- **Dates**: 50 dates
- **Method**: Compute ecliptic, manually convert to equatorial, compare with equatorial flag output
- **Acceptance**: <1e-10° (floating-point noise only)

### Round 4.4 — J2000 Precession Round-Trip

- **Bodies**: 15 bodies
- **Dates**: 50 dates
- **Method**: Compute ecliptic of date, manually precess to J2000, compare with J2K flag output
- **Additional**: Verify J2K suppression for TrueNode/OscuApog/IntpApog/IntpPerg
- **Acceptance**: <1e-8° for Pipeline A, exact match for suppressed bodies

### Round 4.5 — Flag Orthogonality: All 64+ Combinations

- **Bodies**: Sun, Moon, Mars, Jupiter, MeanNode, TrueNode, OscuApog, Cupido
- **Dates**: 10 dates (2 per era)
- **Flags**: All valid combinations of: SPEED, EQUATORIAL, J2000, SIDEREAL, TRUEPOS, NOABERR, NOGDEFL, HELCTR, BARYCTR
- **Check**: No NaN, no infinity, lon in [0,360), lat in [-90,90], speed in reasonable range
- **Acceptance**: Zero NaN/infinity/out-of-range across all combinations

### Round 4.6 — TRUEPOS and NOABERR: Difference Quantification

- **Bodies**: Sun, Moon, Mars, Jupiter, Saturn
- **Dates**: 100 dates from 1800-2200 CE
- **Method**: Compare default vs TRUEPOS, default vs NOABERR, default vs NOGDEFL. Quantify the correction magnitude.
- **Acceptance**: LEB vs Skyfield difference <0.001" for each flag mode

### Round 4.7 — Heliocentric and Barycentric: All Bodies

- **Bodies (helio)**: Mercury–Pluto, Chiron, Ceres–Vesta (not Sun, Moon)
- **Bodies (bary)**: Mercury–Pluto, Earth
- **Dates**: 100 dates
- **Components**: All 6
- **Acceptance**: Position <0.02", speed <0.05 d/d (HELCTR amplification factor documented)

### Round 4.8 — Sidereal + Non-Standard Flag Combinations

- **Purpose**: Test sidereal combined with every other flag
- **Bodies**: Sun, Moon, Mars, MeanNode, TrueNode
- **Flags**: SIDEREAL × {EQUATORIAL, J2000, TRUEPOS, NOABERR, HELCTR, BARYCTR} (12 combos)
- **Dates**: 10 dates
- **Acceptance**: LEB vs Skyfield <TOLS_EXT for each combo

---

## Phase 5: Numerical Stability & Physical Consistency (8 rounds)

### Round 5.1 — Nutation/Obliquity/Delta-T at 30 Epoch Points

- **Dates**: -13000, -10000, -7500, -5000, -3000, -2000, -1000, -500, 0, 500, 1000, 1500, 1700, 1800, 1900, 1950, 2000, 2025, 2050, 2100, 2200, 2500, 3000, 5000, 7500, 10000, 12000, 15000, 17000 CE
- **Metrics**: True obliquity, mean obliquity, dpsi, deps, Delta-T
- **Three-way**: LEB vs Skyfield vs SE
- **Acceptance**: LEB = Skyfield exactly; LEB vs SE within documented model differences

### Round 5.2 — Planetary Distance Extrema: 50 Perigees + 50 Apogees

- **Bodies**: Moon (perigee/apogee), Mars (opposition/conjunction), Venus (inferior/superior conjunction)
- **Method**: Find 50 perigees and 50 apogees for Moon (2000-2025 CE). Evaluate all 6 components.
- **Acceptance**: Distance <1e-8 AU at perigee, position <0.001"

### Round 5.3 — Monotonicity: MeanNode and MeanApog Across Full Range

- **Bodies**: MeanNode(10), MeanApog(12)
- **Dates**: 20,000 dates (1 per ~6 months equivalent)
- **Checks**: MeanNode always retrograde, MeanApog always prograde, no wrap-around anomalies
- **Acceptance**: Zero violations

### Round 5.4 — Sun-Earth Heliocentric Consistency

- **Dates**: 200 dates across full range
- **Verify**: Sun_geo_lon ≈ (Earth_helio_lon + 180) mod 360, lat/dist agreement
- **Acceptance**: <1e-8° position, <1e-12 AU distance

### Round 5.5 — TrueNode/OscuApog Oscillation Amplitude Bounds

- **Bodies**: TrueNode(11), OscuApog(13)
- **Dates**: 5000 dates across full range
- **Checks**: |TrueNode - MeanNode| <2.5°, no sudden jumps >0.5° in 10-day intervals
- **Acceptance**: Zero violations

### Round 5.6 — Speed Sign Consistency

- **Bodies**: All 31
- **Dates**: 1000 dates
- **Checks**: For each body, verify speed sign is physically plausible (e.g., MeanNode always negative, Sun generally positive, retrograde planets can be negative)
- **Additional**: Verify no speed exceeds physical limits (Moon <16°/day, Mercury <2.5°/day, etc.)
- **Acceptance**: Zero unphysical speeds

### Round 5.7 — Distance Positivity and Physical Bounds

- **Bodies**: All 31
- **Dates**: 1000 dates
- **Checks**: dist > 0 for all bodies except nodes (which may be 0), dist within physical bounds (Sun ~0.98-1.02 AU, Moon ~0.0024-0.0027 AU, etc.)
- **Acceptance**: Zero violations

### Round 5.8 — Longitude Ordering Consistency

- **Bodies**: Sun, Moon, Mars, Jupiter
- **Dates**: 10,000 consecutive hours from J2000
- **Checks**: Day-to-day longitude change should be <max physical speed, no random jumps
- **Acceptance**: Zero violations

---

## Phase 6: Cross-Tier Consistency (4 rounds)

### Round 6.1 — Extended vs Medium: Overlap Zone (1560-2640)

- **Bodies**: All 31
- **Dates**: 200 dates in overlap
- **Metrics**: |ext_LEB - med_LEB| per body per component
- **Acceptance**: Differences match known DE440↔DE441 baseline

### Round 6.2 — Extended vs Base: Overlap Zone (1860-2140)

- **Bodies**: All 31
- **Dates**: 200 dates in overlap
- **Same analysis as 6.1 but with de440s (base tier BSP)
- **Acceptance**: Consistent with 6.1 (base uses de440s, medium uses de440)

### Round 6.3 — Cross-Tier Sidereal Consistency

- **Bodies**: Sun, Moon, Mars, MeanNode, Cupido
- **Dates**: 50 dates in overlap
- **Flags**: SIDEREAL, SID+EQ, SID+J2K (3 combos)
- **Method**: Cross-tier sidereal difference should equal cross-tier tropical difference
- **Acceptance**: Delta <0.001"

### Round 6.4 — Cross-Tier Speed Comparison

- **Bodies**: Moon, Mercury, Mars, TrueNode, OscuApog
- **Dates**: 200 dates in overlap
- **All 3 speed components
- **Acceptance**: Consistent with DE440↔DE441 baseline

---

## Phase 7: Exhaustive Ayanamsha Coverage (4 rounds)

### Round 7.1 — All 43 Modes × 10 Bodies × 10 Dates: LEB vs Skyfield

- **Bodies**: Sun, Moon, Mars, Jupiter, Saturn, MeanNode, TrueNode, OscuApog, Chiron, Cupido
- **Dates**: 10 dates (1900, 1920, 1940, 1960, 1980, 2000, 2020, 2040, 2060, 2080)
- **Ayanamsha**: All 43 modes (0-42) + user-defined (255)
- **Flags**: SIDEREAL
- **Acceptance**: LEB = Skyfield within 0.001" for formula-based, identical for star-based (both fallback)

### Round 7.2 — Ayanamsha Stability at Extreme Dates

- **Bodies**: Sun, Moon
- **Dates**: -5000, -3000, -1000, 0, 1000, 2000, 3000, 5000 CE
- **Ayanamsha**: All 43 modes
- **Checks**: No NaN, no infinity, values in [0, 360)
- **Acceptance**: 100% valid output

### Round 7.3 — Ayanamsha Rate Consistency

- **Purpose**: Verify d(ayanamsha)/dt is consistent between LEB and Skyfield
- **Method**: For each mode, compute ayanamsha at t and t+1day, derive rate
- **Dates**: 20 dates from 1800-2200 CE
- **Acceptance**: Rate agreement <1e-8 d/d

### Round 7.4 — User-Defined Ayanamsha: Custom t0 and ayan_t0

- **Purpose**: Verify user-defined sidereal mode works correctly on LEB
- **Method**: Set sid_mode=255 with various t0/ayan_t0 values
- **Test cases**: 10 different t0/ayan_t0 combinations
- **Bodies**: Sun, Moon, Mars
- **Acceptance**: LEB = Skyfield exactly

---

## Phase 8: Pathological & Degenerate Cases (6 rounds)

### Round 8.1 — Polar House Systems

- **Latitudes**: 60°, 63°, 66°, 66.5° (Arctic), 67°, 70°, 75°, 80°, 85°, 89°, -66.5°, -70°, -80°, -89°
- **House systems**: Placidus, Koch, Regiomontanus, Campanus, Topocentric
- **Dates**: Summer solstice, winter solstice, equinox (×3 years: 1900, 2000, 2100)
- **Checks**: No crash, no NaN, cusps in [0,360)
- **Acceptance**: Zero crashes/NaN

### Round 8.2 — Bodies at Maximum Speed

- **Moon**: Find 50 dates with highest |speed| in 2000-2025 CE
- **Mercury**: Find 50 dates with highest |speed| in 2000-2025 CE
- **Evaluate**: All 6 components LEB vs Skyfield at each
- **Acceptance**: Position <0.001" even at maximum speed

### Round 8.3 — Retrograde Stations: Zero-Crossing Precision

- **Bodies**: Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
- **Find**: 20 stations per body (2000-2025 CE)
- **Evaluate**: At station ± [0.001, 0.01, 0.1, 1] days
- **Three-way**: LEB vs Skyfield vs SE at each point
- **Acceptance**: Position <0.001", speed sign agreement when |speed| >0.001 d/d

### Round 8.4 — Conjunctions, Oppositions, Squares

- **Sun-Moon**: 30 conjunctions (new moons) + 30 oppositions (full moons) in 2000-2025 CE
- **Sun-Mars**: 10 oppositions
- **Venus-Jupiter**: 5 conjunctions
- **Method**: Compute angular difference via LEB and Skyfield, compare
- **Acceptance**: Angular difference error <0.001"

### Round 8.5 — Near-Simultaneity: Two Bodies at Same Longitude

- **Purpose**: When two bodies have nearly identical longitude, verify both are precisely positioned
- **Method**: Find 20 planetary conjunctions. At each, verify both body longitudes match Skyfield
- **Acceptance**: Each body <0.001"

### Round 8.6 — JD Extremes: First and Last Valid Segments

- **Bodies**: All 31
- **Dates**: JD at (start + 1 day), (start + 7 days), (end - 7 days), (end - 1 day)
- **Checks**: Valid output, no crash, within tolerance
- **Acceptance**: All valid, position within 5× normal tolerance

---

## Phase 9: Automated Regression Detection (4 rounds)

### Round 9.1 — Golden Value Snapshot: 100 Reference Points

- **Purpose**: Create a set of 100 "golden" LEB evaluations that can be recomputed after any code change
- **Method**: Select 100 body/date/flag combinations covering all pipelines, all flag combos, all eras
- **Output**: JSON file with exact results (18 decimal places)
- **Usage**: After any code change, recompute and diff against golden values
- **Acceptance**: Bit-exact match (any change indicates a regression or intended modification)

### Round 9.2 — Sidereal Golden Values: 50 Reference Points

- **Same as 9.1 but focused on sidereal calculations
- **Cover**: All 4 bug fix scenarios, all pipeline types, multiple ayanamshas
- **Output**: JSON file
- **Acceptance**: Bit-exact match

### Round 9.3 — Speed Golden Values: 50 Reference Points

- **Same as 9.1 but focused on velocity calculations
- **Cover**: All speed components, all flag combos, edge cases (max speed, station)
- **Output**: JSON file

### Round 9.4 — Cross-Mode Invariant Checks

- **Purpose**: Define invariants that must always hold regardless of code changes
- **Examples**:
  - `|tropical_lon - sidereal_lon - ayanamsha| < 1e-10°`
  - `MeanNode_speed < 0` always
  - `|Sun_geo_lon - (Earth_helio_lon + 180) mod 360| < 1e-8°`
  - `TrueNode_SID == TrueNode_SID_J2K` exactly
- **Output**: Executable test script
- **Acceptance**: All invariants hold

---

## Phase 10: Full Statistical Report (4 rounds)

### Round 10.1 — Grand Matrix: 31 Bodies × 1000 Dates × 8 Flag Combos

- **Bodies**: All 31
- **Dates**: 1000 uniformly distributed
- **Flags**: default, EQUATORIAL, J2000, EQ+J2K, SIDEREAL, SID+EQ, SID+J2K, SID+EQ+J2K
- **Metrics**: Per body per flag: N, mean, P50, P95, P99, max of |dLon|, |dLat|, |dDist|
- **Stratification**: By era (5 strata)
- **Output**: Grand table (31 × 8 × 5 = 1,240 cells)

### Round 10.2 — Precision Degradation Curve: Error vs Time

- **Bodies**: Sun, Moon, Mars, Jupiter, MeanNode
- **Method**: Compute error at 500 dates. Plot error vs |T - J2000|.
- **Output**: One plot per body showing degradation (or lack thereof)
- **Purpose**: Quantify long-term precision stability

### Round 10.3 — Cross-Tier Precision Comparison Table

- **Method**: At 200 overlap dates, compute precision for all 3 tiers
- **Output**: Side-by-side table: base vs medium vs extended
- **Acceptance**: Extended comparable to medium in overlap zone

### Round 10.4 — Final Acceptance Matrix

- **Output**: Matrix of 31 bodies × all phases/rounds with PASS/WARN/FAIL
- **Criteria**:
  - PASS: Within tolerance with ≥5× headroom
  - WARN: Within tolerance with <5× headroom
  - FAIL: Exceeds tolerance
- **Verdict**: Overall PASS requires zero FAIL cells
- **Known limitations**: Documented and excluded from PASS/FAIL

---

## Test File Mapping

| Phase/Round | Automated Test File | System |
|-------------|-------------------|--------|
| Phase 3 (Bug regression) | `tests/test_leb/test_fast_calc.py::TestSiderealRegressionBug*` | LEB vs Skyfield |
| Phase 3 (Bug regression) | `compare_scripts/tests/test_compare_sidereal_regression.py` | Swiss vs libephemeris |
| Phase 1, 4 (Sidereal) | `tests/test_leb/compare/extended/test_extended_sidereal.py` | LEB vs Skyfield |
| Phase 1 (Asteroids) | `tests/test_leb/compare/extended/test_extended_asteroids.py` | LEB vs Skyfield |
| Phase 1 (Distances) | `tests/test_leb/compare/extended/test_extended_distances.py` | LEB vs Skyfield |
| Phase 1 (Flag velocity) | `tests/test_leb/compare/extended/test_extended_flags.py` | LEB vs Skyfield |
| Phase 1 (Lunar distance) | `tests/test_leb/compare/extended/test_extended_lunar.py` | LEB vs Skyfield |
| Phase 1 (Planets) | `tests/test_leb/compare/extended/test_extended_planets.py` | LEB vs Skyfield |
| Phase 1 (Velocities) | `tests/test_leb/compare/extended/test_extended_velocities.py` | LEB vs Skyfield |
| Phase 2 (Nodes/Lilith) | `compare_scripts/tests/test_compare_lunar_nodes_lilith.py` | Swiss vs libephemeris |
| Phase 2 (Sidereal) | `compare_scripts/tests/test_compare_sidereal.py` | Swiss vs libephemeris |
| Phase 6 (Cross-tier) | `tests/test_leb/compare/crosstier/` | LEB vs LEB |

---

## Tolerances Reference (Updated from V1)

| Body Category | Position vs Skyfield (LEB error) | Position vs SE (model+LEB) | Speed vs Skyfield | Headroom (measured) |
|---|---|---|---|---|
| Sun | <0.001" | <0.005" | <0.001 d/d | 500× |
| Moon | <0.500" | <0.050" | <0.002 d/d | 1,470× |
| Mercury–Neptune | <0.001" | <0.030" | <0.001 d/d | >60× |
| Pluto | <0.025" | <0.100" | <0.001 d/d | >25,000× |
| Earth | <0.001" | <0.005" | <0.001 d/d | >1,000× |
| MeanNode/MeanApog | <0.005" | ~0.1" (model) | <0.001 d/d | 1.7× |
| TrueNode | <5" (model) | ~2-5" (model) | <0.01 d/d | >5M× |
| OscuApog | <100" (model) | ~57-90" (model) | <0.05 d/d | enormous |
| IntpApog/IntpPerg | <5000" (model) | ~200-4600" (model) | <1.5 d/d | N/A |
| Chiron | <0.001" | <0.005" | <0.001 d/d | 20× |
| Ceres–Vesta | <0.001" | <0.100" | <0.001 d/d | 20× |
| Uranians (40-47) | <35" (model) | ~5-33" (model) | <0.001 d/d | enormous |
| Transpluto (48) | <10" (model) | ~1-7" (model) | <0.001 d/d | enormous |

---

## Known Limitations (Not Bugs)

1. **Chiron SPK range**: ~660-4600 CE only
2. **IntpApog/IntpPerg SE range**: ~-3000 to +2900 CE only in pyswisseph
3. **MeanNode/MeanApog Meeus polynomial**: Degraded precision beyond ±20 centuries from J2000 (~0.003")
4. **Neptune/Uranus Skyfield speed artifact**: Isolated finite-difference spikes
5. **Transpluto velocity -180 d/d artifact**: Heliocentric coordinate system for ~55 AU body
6. **HELCTR error amplification**: 8-11× for inner planets (inherent)
7. **DE431↔DE441 baseline**: Moon ~90" at 0 CE, ~280" at 2650 CE
8. **Pholus**: Not in LEB (no SPK data)
9. **Sidereal houses architectural diff**: SE uses GMST+mean_obliquity, LE uses GAST+true_obliquity (~5-16" offset)
10. **IAU 2006 vs IAU 1976 precession**: ~8-18" for J2000 precession vs pyswisseph
11. **Nutation polynomial degradation**: ~0.003" beyond ±20 centuries
12. **LEB sidereal speed**: Uses mean precession rate (analytical) vs Skyfield true ayanamsha rate (finite-diff) — ~5e-5 °/day cross-mode delta
13. **Asteroid SPK coverage**: ~1900-2100 CE only; outside = Keplerian fallback (catastrophically wrong)

---

## Execution Summary

| Metric | V1 | V2 |
|--------|------|------|
| Phases | 8 | 10 |
| Rounds | 38 | 70 |
| Bodies tested | 31 | 31 |
| Date range | -5000 to +5000 | -5000 to +5000 |
| Total body×date evaluations (est.) | ~50,000 | ~5,000,000 |
| Flag combinations | ~8 | ~64+ |
| Ayanamsha modes tested | 27 | 44 |
| Automated regression tests (new) | 0 | ~537 |
| Three-way comparisons | ~300 | ~30,000 |
| Golden reference values | 0 | 200 |
| Estimated execution time | 45-60 min | 8-12 hours |
