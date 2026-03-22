# Hyper-Validation Plan: libephemeris vs pyswisseph

## Objective

Exhaustive 1:1 comparison of **every callable API function** between libephemeris
and pyswisseph 2.10.03, with **1000+ individual test rounds** covering:
- All 101 pyswisseph functions
- All 23 celestial bodies (0–22)
- All 24 house systems
- All 47 ayanamsa modes
- All flag combinations
- Multiple Julian dates spanning the full ephemeris range
- Multiple geographic locations
- Edge cases (poles, equator, sign boundaries, negative values)

## Test Matrix Summary

| # | Category | Function(s) | Axes | Rounds |
|---|----------|-------------|------|--------|
| A | Planetary positions | `calc_ut` | 23 bodies × 10 JDs × 5 flags | 1,150 |
| B | Houses (JD-based) | `houses`, `houses_ex`, `houses_ex2` | 24 sys × 5 locs × 4 JDs | 480 |
| C | Houses (ARMC-based) | `houses_armc`, `houses_armc_ex2` | 24 sys × 4 ARMC × 3 lats | 288 |
| D | Fixed stars | `fixstar2_ut`, `fixstar2_mag` | 20 stars × 3 JDs × 2 flags | 120 |
| E | Ayanamsa | `get_ayanamsa_ex_ut` | 47 modes × 5 JDs | 235 |
| F | Degree splitting | `split_deg` | 20 values × 16 flag combos | 320 |
| G | Nodes & apsides | `nod_aps_ut` | 12 bodies × 3 methods × 3 JDs | 108 |
| H | Solar eclipses | `sol_eclipse_*` | 10 start JDs × 4 funcs | 40 |
| I | Lunar eclipses | `lun_eclipse_*` | 10 start JDs × 3 funcs | 30 |
| J | Occultations | `lun_occult_*` | 5 bodies × 4 JDs | 20 |
| K | Utility math | `degnorm`, `radnorm`, `difdeg2n`, etc. | 12 funcs × 20 values | 240 |
| L | Rise/set/transit | `rise_trans`, `rise_trans_true_hor` | 5 bodies × 3 locs × 4 types | 60 |
| M | Phenomena | `pheno_ut` | 10 bodies × 3 JDs | 30 |
| N | Time functions | `julday`, `revjul`, `deltat`, `sidtime`, etc. | 6 funcs × 20 values | 120 |
| O | House position | `house_pos` | 20 test cases | 20 |
| P | Coord transforms | `cotrans`, `cotrans_sp` | 20 values × 2 funcs | 40 |
| Q | Refraction | `refrac`, `refrac_extended` | 15 alts × 2 dirs | 30 |
| R | Horizontal coords | `azalt`, `azalt_rev` | 15 cases × 2 funcs | 30 |
| S | Orbital elements | `get_orbital_elements` | 10 bodies × 2 JDs | 20 |
| T | Crossings | `solcross_ut`, `mooncross_ut`, `mooncross_node_ut` | 20 cases | 20 |
| U | Heliacal | `heliacal_ut`, `heliacal_pheno_ut` | 10 cases | 10 |
| V | String formatting | `cs2lonlatstr`, `cs2timestr`, `get_planet_name`, etc. | 50 cases | 50 |
| W | SE_AST_OFFSET | `calc_ut` with 10001–10020 | 6 bodies × 3 JDs × 2 flags | 36 |
| X | Sidereal positions | `calc_ut` + `set_sid_mode` | 10 bodies × 5 aya × 2 JDs | 100 |
| Y | Gauquelin sectors | `gauquelin_sector` | 5 bodies × 3 locs | 15 |
| Z | Planet names/consts | `get_planet_name`, `house_name`, `get_ayanamsa_name` | 80 cases | 80 |
| AA | ET/UT conversions | `utc_to_jd`, `jdet_to_utc`, `jdut1_to_utc` | 20 cases | 20 |
| AB | Delta T extended | `deltat_ex` | 20 JDs | 20 |
| AC | Misc utilities | `day_of_week`, `date_conversion`, `d2l`, `csnorm`, etc. | 60 cases | 60 |
| **TOTAL** | | | | **3,782** |

## Tolerances

| Metric | Tolerance | Notes |
|--------|-----------|-------|
| Longitude/Latitude | 0.001" (arcsecond) | ~2.78e-7 degrees |
| Distance | 1e-7 AU | |
| Speed | 0.001"/day | ~2.78e-7 deg/day |
| House cusps | 0.001" | |
| House cusp speeds | 0.01 deg/day | Numerical vs analytical |
| Ayanamsa | 0.001" | |
| Time (JD) | 1e-8 days | ~0.86ms |
| Sidereal time | 1e-8 hours | |
| String output | Exact match | |
| Integer output | Exact match | |
| retflag | Exact match or documented divergence | |

## Known Accepted Divergences (from Phase 1 verification)

These will be flagged as KNOWN, not FAIL:
1. `version` string (intentional)
2. `d2l` negative values (unsigned overflow in pyswisseph)
3. `cs2degstr()` (pyswisseph segfaults)
4. `deltat` for JD > ~2488070 (future dates >2050): up to 3.4s
5. `refrac`/`azalt`: up to 15"
6. `rise_trans_true_hor` negative horhgt: ~100s
7. `heliacal_ut`: up to ~2 days
8. `mooncross_node_ut`: up to ~69s
9. `nod_aps_ut` planetary nodes: 20–700"
10. `sol_eclipse_where` geopos[2:9]: bonus data
11. Exotic ayanamsa modes: up to 25"
12. Fixed star magnitudes: catalog version diff
13. IntpApog/IntpPerg: 374"/8414"
14. Uranian hypotheticals: ~39"
15. Transpluto (body 9): ~41" lon
16. houses_ex2 cusp speeds Koch/Porphyry: algorithmic diff
17-23. Eclipse/occultation retflag diffs, heliacal pheno diffs

## Execution

Run: `.venv/bin/python3 scripts/hyper_validate.py`

Output: JSON report + console summary with PASS/FAIL/KNOWN counts.
