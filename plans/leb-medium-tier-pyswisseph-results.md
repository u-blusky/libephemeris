# LEB vs pyswisseph — Medium Tier Results

## Summary

LEB mode (medium tier, de440.bsp, 1550-2650) tested against pyswisseph via `compare_scripts/tests/`.

**Total: 3,987 passed, 230 skipped, 20 pre-existing failures, 2 files timed out (>20 min each)**

## Per-file Results

| File | Passed | Skipped | Failed | Notes |
|------|--------|---------|--------|-------|
| test_compare_planets.py | 563 | 3 | **1** | Moon daily motion 0.001048 vs 0.001 tol |
| test_compare_houses.py | 337 | 0 | 0 | |
| test_compare_crossings.py | 27 | 0 | 0 | |
| test_compare_sidereal.py | 201 | 0 | 0 | |
| test_compare_lunar.py | 31 | 0 | 0 | |
| test_compare_lunar_nodes_lilith.py | 1628 | 0 | 0 | |
| test_compare_hypothetical.py | 90 | 0 | 0 | |
| test_compare_helio_bary.py | 500 | 3 | 0 | Sun helio skipped |
| test_compare_eclipses.py | 21 | 0 | 0 | |
| test_compare_phenomena.py | 6 | 34 | 0 | swe_pheno_ut not impl |
| test_compare_observations.py | 31 | 0 | 0 | |
| test_compare_elongation.py | 418 | 35 | 0 | pheno cross-val skipped |
| test_compare_coordinates.py | 280 | 0 | 0 | |
| test_compare_time.py | 140 | 0 | 0 | 1 xpassed |
| test_compare_utilities.py | 77 | 0 | 0 | |
| test_compare_fixedstars.py | 287 | 27 | 0 | star name mismatches |
| test_compare_nogdefl.py | 59 | 0 | 0 | |
| test_compare_rise_transit.py | 10 | 0 | 0 | |
| test_compare_minor_bodies.py | 28 | 96 | 0 | SPK files missing |
| test_compare_houses_ext.py | 47 | 5 | 0 | Gauquelin skipped |
| test_compare_occultations.py | 12 | 0 | 0 | |
| test_compare_calc_pctr.py | 94 | 0 | 0 | |
| test_compare_planetary_moons.py | 29 | 29 | 0 | SPK files missing |
| test_compare_benchmark.py | 0 | 16 | 0 | benchmarks skipped |
| test_compare_orbital.py | 19 | 0 | **20** | **Pre-existing** (same in Skyfield mode) |
| test_compare_heliacal.py | — | — | — | Timed out (>20 min) |
| test_compare_planet_occultations.py | — | — | — | Timed out (>20 min) |

## LEB-specific Issues

### 1. Moon Daily Motion (marginal, 1 test)
- `TestInnerPlanets::test_inner_planet_daily_motion[1-Moon]`
- Error: 0.001048 deg/day vs tolerance 0.001 deg/day
- Passes in Skyfield mode — LEB Chebyshev approximation amplifies speed error
- Decision pending: loosen tolerance or accept as known limitation

### 2. Pre-existing orbital failures (NOT LEB-related)
- `test_compare_orbital.py` — 20 failures identical in both Skyfield and LEB mode
- `swe_nod_aps_ut()` apse calculations differ by degrees from pyswisseph
- Known libephemeris vs SE algorithmic difference, unrelated to LEB

### 3. Timed-out files
- `test_compare_heliacal.py` (40 tests) — heliacal rise/set is iterative, very slow
- `test_compare_planet_occultations.py` (45 tests) — occultation search is iterative, very slow
- Both use iterative search algorithms that call swe_calc thousands of times
- Need dedicated long-running validation or per-test execution
