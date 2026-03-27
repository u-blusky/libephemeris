# Verification Status

Last updated: 2026-03-27

## Bug Fixes (this session)

| Bug | File | Description | Status |
|-----|------|-------------|--------|
| BUG-006 | `horizons_backend.py` | `_calc_uranian`: swapped args to `calc_uranian_planet(body_id, jd_tt)` + wrong 3-value unpack from 6-tuple return. Removed redundant finite-difference speed computation. | Fixed (`463bef2`) |

## Pytest Tests Added (this session)

Committed in `7888fc5`. 107 new test functions across 5 files:

| File | Tests | Coverage |
|------|-------|----------|
| `tests/test_planets/test_wrap_and_invalid_inputs.py` | 23 | 360 wrap-around, ECL_NUT special body, NaN/Inf inputs |
| `tests/test_coords/test_azalt_refrac_comprehensive.py` | 43 | azalt, azalt_rev round-trip, refrac, refrac_extended |
| `tests/test_time/test_tai_time_functions.py` | 22 | UTC-TAI-UTC, TT-TAI-TT round-trips, leap seconds |
| `tests/test_planets/test_horizons_analytical.py` | 19 | Horizons offline (Mean Node, Mean Apogee, Uranians), calc_mode |
| `tests/test_houses/test_house_pos_placement.py` | fix | Ascendant boundary edge case (12.999 vs 1.0 wrap) |

## Wave 1 Standalone Verification Scripts

All in `tasks/scripts/`. Run with `.venv/bin/python tasks/scripts/<name>.py`.

| Script | Checks | Passed | Failed | Time | Sections covered |
|--------|--------|--------|--------|------|------------------|
| `verify_positions.py` | 20,650 | 20,650 | 0 | 3.7s | 1.1, 1.2, 1.6, 23.4 |
| `verify_houses.py` | 44,870 | 44,854 | 16 | 0.6s | 3.1-3.7 |
| `verify_flags.py` | 77,800 | 77,620 | 180 | 10.5s | 2.1-2.6 |
| `verify_ayanamsha_eclipses_stars.py` | 3,143 | 3,143 | 0 | 18.8s | 4, 5, 6, 7, 9 |
| `verify_time_coords_nodes_pheno.py` | 5,579 | 5,579 | 0 | 0.4s | 10, 11, 12, 17 |
| `verify_batch12.py` | 171 | 171 | 0 | <1s | BUG-006 + Batch 12 spot checks |
| **Total** | **152,213** | **152,017** | **196** | **~34s** | |

## Open Issues

### 1. Sunshine/Makransky house system ('i') at high latitudes
- **16 failures** at latitude 64N (Reykjavik): cusps 2, 6, 8, 12 differ by 30-100 degrees from pyswisseph
- Works correctly at all other tested latitudes (0, 41.9, 51.5, -33.9, 35.7)
- Needs investigation in `houses.py`

### 2. XYZ+HELCTR Moon (flag combo script)
- **180 failures** reported in `verify_flags.py` for Moon with SEFLG_HELCTR | SEFLG_XYZ
- Could not reproduce standalone (Moon helio XYZ returns finite values in isolated test)
- May be a script issue rather than a library bug; needs re-examination

## Sections Not Yet Covered by Verification Scripts

From the 26-section plan (~490K checks target):

| Section | Description | Status |
|---------|-------------|--------|
| 1.3 | LEB vs Skyfield (14 bodies x 500 dates) | Not started |
| 1.4 | LEB vs Skyfield flag modes | Not started |
| 1.5 | Horizons vs Skyfield (requires HTTP) | Not started |
| 1.7 | LEB2 vs LEB1 compression integrity | Not started |
| 1.8 | Triple cross-validation (Skyfield/LEB/Horizons) | Not started |
| 8 | Heliacal visibility (vis_limit_mag) | Not started (slow API) |
| 13 | Uranians deep (sidereal, helio) | Partial (BUG-006 fix covers basic) |
| 14 | Asteroids deep (SE_AST_OFFSET mapping) | Not started |
| 15 | Planetary moons | Not started |
| 16 | Orbital elements deep | Not started |
| 18 | Crossings and stations deep | Partial (existing tests) |
| 19 | LEB-specific (reader API, fast_calc, compression) | Not started |
| 20 | Horizons-specific (HTTP pipeline, fallback) | Partial (analytical only) |
| 21 | State management deep | Not started |
| 22 | EphemerisContext concurrency | Not started |
| 23 | Edge cases deep (date boundaries, special bodies) | Partial |
| 24 | Utility functions deep | Not started |
| 25 | Arabic parts deep | Not started |
| 26 | Golden regression snapshot | Exists already in repo |

## Previously Fixed Bugs (prior sessions)

| Bug | Description |
|-----|-------------|
| BUG-001 | Lunar perigee perturbation coefficients |
| BUG-002 | (details in tasks/results/bugs.md) |
| BUG-003 | (details in tasks/results/bugs.md) |
| BUG-004 | (details in tasks/results/bugs.md) |
| BUG-005 | (details in tasks/results/bugs.md) |

## Total Test Inventory

77 untracked test files with ~1055 test functions across 17 directories, plus the 5 Wave 1 verification scripts covering 152K checks.
