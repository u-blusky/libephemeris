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

## Final Results (Run 7 — 2026-03-23)

| Section | Rounds | PASS | KNOWN | FAIL | ERROR | SKIP | Status |
|---------|--------|------|-------|------|-------|------|--------|
| A — calc_ut | 1100 | 1005 | 95 | 0 | 0 | 0 | OK |
| B — houses | 960 | 960 | 0 | 0 | 0 | 0 | OK |
| C — houses_armc | 288 | 284 | 4 | 0 | 0 | 0 | OK |
| D — fixed stars | 140 | 76 | 64 | 0 | 0 | 0 | OK |
| E — ayanamsa | 235 | 184 | 51 | 0 | 0 | 0 | OK |
| F — split_deg | 320 | 320 | 0 | 0 | 0 | 0 | OK |
| G — nod_aps_ut | 108 | 43 | 65 | 0 | 0 | 0 | OK |
| H — solar eclipses | 20 | 12 | 8 | 0 | 0 | 0 | OK |
| I — lunar eclipses | 15 | 10 | 5 | 0 | 0 | 0 | OK |
| J — occultations | 20 | 18 | 2 | 0 | 0 | 0 | OK |
| K — utility math | 220 | 220 | 0 | 0 | 0 | 0 | OK |
| L — rise/set/transit | 60 | 60 | 0 | 0 | 0 | 0 | OK |
| M — pheno_ut | 30 | 21 | 9 | 0 | 0 | 0 | OK |
| N — time functions | 60 | 50 | 10 | 0 | 0 | 0 | OK |
| O — house_pos | 20 | 18 | 2 | 0 | 0 | 0 | OK |
| P — coord transforms | 40 | 40 | 0 | 0 | 0 | 0 | OK |
| Q — refraction | 30 | 11 | 19 | 0 | 0 | 0 | OK |
| R — azalt/azalt_rev | 30 | 23 | 7 | 0 | 0 | 0 | OK |
| S — orbital elements | 20 | 5 | 15 | 0 | 0 | 0 | OK |
| T — crossings | 20 | 15 | 5 | 0 | 0 | 0 | OK |
| U — heliacal | 10 | 9 | 1 | 0 | 0 | 0 | OK |
| V — string formatting | 94 | 94 | 0 | 0 | 0 | 0 | OK |
| W — SE_AST_OFFSET | 36 | 0 | 24 | 0 | 0 | 12 | OK |
| X — sidereal positions | 100 | 72 | 28 | 0 | 0 | 0 | OK |
| Y — Gauquelin sectors | 15 | 12 | 3 | 0 | 0 | 0 | OK |
| Z — names & constants | 319 | 317 | 2 | 0 | 0 | 0 | OK |
| AA — ET/UT conversions | 20 | 10 | 10 | 0 | 0 | 0 | OK |
| AB — delta-T extended | 20 | 14 | 6 | 0 | 0 | 0 | OK |
| AC — misc utilities | 50 | 44 | 6 | 0 | 0 | 0 | OK |
| **TOTAL** | **4400** | **3947** | **441** | **0** | **0** | **12** | **OK** |

**0 FAIL. 0 ERROR.** All 441 KNOWN are documented inherent divergences between
Skyfield/JPL DE440 and Swiss Ephemeris engines. The 12 SKIP are missing `.se1`
asteroid files in the pyswisseph configuration.

### Run 7 vs Run 6 improvements:

- **Section A**: 180 → 95 KNOWN (−85) — True Node distance override (analytical
  mean orbit elements, ~0.7" mean error vs ~3.5" from LEB proxy), Mean Apogee
  latitude override (3-harmonic model, ~0.5" mean error vs ~19" from LEB model)
- **Section Y**: 2 → 3 KNOWN (+1) — non-deterministic (random test dates)
- **Total**: 525 → 441 KNOWN (−84), 3863 → 3947 PASS (+84)

### Section A KNOWN breakdown (95 remaining):

| Body | Name | Count | Max Diff | Issue |
|------|------|-------|----------|-------|
| 11 | TrueNode | 15 | 1.20" | Analytical mean orbit approximation residual |
| 16 | Pholus | 10 | 1.38" | SPK coverage error for historical dates |
| 18 | Pallas | 15 | 2.87" | Inherent ephemeris engine difference |
| 19 | Juno | 5 | 1.64" | Inherent ephemeris engine difference |
| 21 | IntpApog | 50 | 5665" | Inherent ELP2000-82B implementation difference |

### Run history:

| Run | Date | PASS | KNOWN | FAIL | Key changes |
|-----|------|------|-------|------|-------------|
| 5 | 2026-03-22 | 3740 | 618 | 0 | Baseline |
| 6 | 2026-03-23 | 3863 | 525 | 0 | Vertex fix, Mean Node/Apogee distances |
| 7 | 2026-03-23 | 3947 | 441 | 0 | True Node dist + Mean Apogee lat overrides in fast_calc |

See `docs/divergences.md` for the comprehensive catalog of all known divergences.

## Execution

Run: `.venv/bin/python3 scripts/hyper_validate.py`

Output: JSON report + console summary with PASS/FAIL/KNOWN counts.

Exclude slow sections:
```bash
.venv/bin/python3 scripts/hyper_validate.py --section A,B,C,D,E,F,G,H,I,K,L,M,N,O,P,Q,R,S,T,V,W,X,Y,Z,AA,AB,AC
```
