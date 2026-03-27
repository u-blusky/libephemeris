# Verification Status

Last updated: 2026-03-27

## Bug Fixes

| Bug | Commit | File | Description |
|-----|--------|------|-------------|
| BUG-006 | `463bef2` | `horizons_backend.py` | `_calc_uranian`: swapped args to `calc_uranian_planet(body_id, jd_tt)` + wrong 3-value unpack from 6-tuple. Removed redundant finite-difference speed. |
| BUG-008 | `29150e0` | `planets.py` | Skyfield `reify` descriptor for `t.P` stores result as `t.precession_matrix`, shadowing the method. Replaced `t.P` with `t.precession_matrix()` in 3 locations + replaced `pos.frame_latlon(mean_equator_and_equinox_of_date)` with direct matrix multiplication to avoid Skyfield internal `t.P` access. |

## Standalone Verification Scripts

All in `tasks/scripts/`. Run with `.venv/bin/python tasks/scripts/<name>.py`.

| Script | Checks | Passed | Failed | Time | Sections |
|--------|--------|--------|--------|------|----------|
| `verify_positions.py` | 20,650 | 20,650 | 0 | 3.7s | 1.1, 1.2, 1.6, 23.4 |
| `verify_houses.py` | 44,870 | 44,854 | 16 | 0.6s | 3.1-3.7 |
| `verify_flags.py` | 77,200 | 77,200 | 0 | 18.6s | 2.1-2.6 |
| `verify_ayanamsha_eclipses_stars.py` | 3,143 | 3,143 | 0 | 18.8s | 4, 5, 6, 7, 9 |
| `verify_time_coords_nodes_pheno.py` | 5,579 | 5,579 | 0 | 0.4s | 10, 11, 12, 17 |
| `verify_leb_uranians_elements.py` | 17,393 | 17,393 | 0 | 11.6s | 1.3, 1.4, 13, 14, 16, 18 |
| `verify_state_utils_arabic.py` | 8,480 | 8,480 | 0 | 0.5s | 21, 22, 23, 24, 25, 26 |
| `verify_batch12.py` | 171 | 171 | 0 | <1s | BUG-006 + Batch 12 |
| **Total** | **177,486** | **177,470** | **16** | **~55s** | |

The 16 remaining failures are all in the Sunshine/Makransky house system ('i') at latitude 64N (known limitation).

## Pytest Test Files Added

### Batch 12 (committed `7888fc5`)

| File | Tests | Coverage |
|------|-------|----------|
| `test_wrap_and_invalid_inputs.py` | 23 | 360 wrap-around, ECL_NUT, NaN/Inf |
| `test_azalt_refrac_comprehensive.py` | 43 | azalt, azalt_rev, refrac, refrac_extended |
| `test_tai_time_functions.py` | 22 | UTC-TAI, TT-TAI round-trips, leap seconds |
| `test_horizons_analytical.py` | 19 | Horizons offline, calc_mode, Uranians |
| `test_house_pos_placement.py` | fix | Ascendant boundary edge case |

### Wave 2 (committed `2e07aba`)

| File | Tests | Coverage |
|------|-------|----------|
| `test_flag_pair_robustness.py` | 940 | 30 flag pairs, XYZ/RADIANS coherence |
| `test_positions_vs_swisseph.py` | 340 | 10 bodies, geocentric/helio/eq/sidereal vs pyswisseph |
| `test_time_vs_swisseph.py` | 139 | julday, revjul, deltat, sidtime, utc_to_jd |
| `test_ayanamsha_vs_swisseph.py` | 337 | 43 ayanamsha modes + sidereal positions |

## Verification Plan Coverage

| Section | Description | Status |
|---------|-------------|--------|
| 1.1 | Skyfield vs pyswisseph (22 bodies × 100 dates) | Done (verify_positions) |
| 1.2 | Flag variants (10 flags × 10 bodies × 50 dates) | Done (verify_positions) |
| 1.3 | LEB vs Skyfield (14 bodies × 100 dates) | Done (verify_leb_uranians_elements) |
| 1.4 | LEB flag modes (4 flags × 14 bodies × 50 dates) | Done (verify_leb_uranians_elements) |
| 1.5 | Horizons vs Skyfield | Skipped (requires HTTP) |
| 1.6 | Heliocentric positions | Done (verify_positions) |
| 1.7 | LEB2 vs LEB1 | Not started |
| 1.8 | Triple cross-validation | Not started |
| 2.1-2.6 | Flag combinations (13 × 22 × 50) | Done (verify_flags) |
| 3.1-3.7 | House systems (24 systems × 6 locations) | Done (verify_houses) |
| 4 | Ayanamsha (43 modes) | Done (verify_ayanamsha_eclipses_stars) |
| 5-6 | Solar/lunar eclipses | Done (verify_ayanamsha_eclipses_stars) |
| 7 | Rise/set/transit | Done (verify_ayanamsha_eclipses_stars) |
| 8 | Heliacal visibility | Skipped (slow) |
| 9 | Fixed stars | Done (verify_ayanamsha_eclipses_stars) |
| 10 | Coordinate transforms | Done (verify_time_coords_nodes_pheno) |
| 11 | Time functions | Done (verify_time_coords_nodes_pheno) |
| 12 | Lunar nodes/apsides | Done (verify_time_coords_nodes_pheno) |
| 13 | Uranians deep | Done (verify_leb_uranians_elements) |
| 14 | Asteroids deep | Done (verify_leb_uranians_elements) |
| 15 | Planetary moons | Skipped (requires SPK kernels) |
| 16 | Orbital elements | Done (verify_leb_uranians_elements) |
| 17 | Phenomena | Done (verify_time_coords_nodes_pheno) |
| 18 | Crossings/stations | Done (verify_leb_uranians_elements) |
| 19 | LEB-specific | Not started |
| 20 | Horizons-specific | Partial (analytical only) |
| 21-22 | State/context | Done (verify_state_utils_arabic) |
| 23 | Edge cases | Done (verify_state_utils_arabic) |
| 24 | Utility functions | Done (verify_state_utils_arabic) |
| 25 | Arabic parts | Done (verify_state_utils_arabic) |
| 26 | Golden regression | Exists in repo already |

## Known Limitations

### Sunshine/Makransky house system ('i') at high latitudes
- 16 failures at latitude 64N: cusps 2, 6 differ by 30-100 degrees from pyswisseph
- Root cause: both 'I' (Treindl) and 'i' (Makransky) dispatch to the same `_houses_sunshine()` function, but 'i' needs a separate Makransky algorithm for meridian distance > 90
- Divergence starts when 2/3 × NSA > 90 (NSA > 135, roughly lat > 58 in winter)
- Clean-room Makransky implementation needed (cannot reference Swiss Ephemeris C source)

## Total Inventory

- **81 new test files** with ~2,800 test functions
- **8 verification scripts** covering 177,486 checks in ~55 seconds
- **2 bugs fixed** (BUG-006, BUG-008)
- **1 known limitation** documented (Sunshine 'i' at high latitudes)
