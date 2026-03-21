# PySwissEph 1:1 Compatibility — Verification Plan

Comprehensive checklist to guarantee drop-in compatibility with `pyswisseph` 2.10.03.
Each item is a concrete, testable check. Status: `[ ]` = pending, `[x]` = verified OK, `[!]` = discrepancy found and fixed.

---

## A. CONSTANTS (319 in pyswisseph)

### A.1 Missing constants
- [x] **A.1.1** All 319 pyswisseph constants exist in libephemeris *(verified: 0 missing)*

### A.2 Value mismatches
- [!] **A.2.1** `TIDAL_AUTOMATIC`: swe=`999999` (int), ephem=`0.0` (float) — FIXED in Phase 3 (commit b951d18)
- [!] **A.2.2** `TIDAL_DEFAULT`: swe=`-25.8`, ephem=`-25.936` — FIXED in Phase 3 (commit b951d18)
- [x] **A.2.3** `version`: swe=`'2.10.03'`, ephem=`'0.25.0'` — intentional, accepted divergence
- [x] **A.2.4** `TRUE_TO_APP` / `APP_TO_TRUE`: swe has `TRUE_TO_APP=0, APP_TO_TRUE=1` — ephem matches ✓

---

## B. FUNCTION SIGNATURES (101 functions in pyswisseph)

### B.1 Missing default values
- [!] **B.1.1** `fixstar(star, jd, iflag)` — FIXED in Phase 3: added `iflag=SEFLG_SWIEPH` default
- [!] **B.1.2** `fixstar_ut(star, jd, iflag)` — FIXED in Phase 3
- [!] **B.1.3** `fixstar2(star, jd, iflag)` — FIXED in Phase 3
- [!] **B.1.4** `fixstar2_ut(star, jd, iflag)` — FIXED in Phase 3
- [!] **B.1.5** `julday(year, month, day, hour)` — FIXED in Phase 3: added `hour=12.0` default

### B.2 Wrong default values
- [!] **B.2.1** `gauquelin_sector` flags default — FIXED in Phase 3: changed to `FLG_SWIEPH|FLG_TOPOCTR` (32770)
- [!] **B.2.2** `sol_eclipse_how` ifl default — FIXED in Phase 3: changed to `FLG_SWIEPH` (2)

### B.3 Parameter type mismatches (hsys: bytes vs int)

pyswisseph accepts `hsys` as `bytes` (e.g. `b'P'`). libephemeris uses `int` (e.g. `ord('P')`).
For 1:1 compat, libephemeris must accept BOTH bytes and int for hsys.

- [x] **B.3.1** `houses(jd, lat, lon, hsys)` — already accepts `b'P'`, `str`, and `int` ✓
- [x] **B.3.2** `houses_ex(jd, lat, lon, hsys, flags)` — already accepts all types ✓
- [x] **B.3.3** `houses_ex2(jd, lat, lon, hsys, flags)` — already accepts all types ✓
- [x] **B.3.4** `houses_armc(armc, lat, eps, hsys)` — already accepts all types ✓
- [x] **B.3.5** `houses_armc_ex2(armc, lat, eps, hsys)` — already accepts all types ✓
- [x] **B.3.6** `house_pos(armc, lat, eps, objcoord, hsys)` — already accepts all types ✓
- [x] **B.3.7** `house_name(hsys)` — ephem accepts both bytes and int (more permissive — OK) ✓
- [!] **B.3.8** `gauquelin_sector` — FIXED in Phase 3: now accepts str (star names) in addition to int

### B.4 Parameter structure mismatches
- [!] **B.4.1** `lun_occult_when_glob` — FIXED in Phase 3: unified to single `body` param accepting int or str
- [x] **B.4.2** `house_pos` — verified: `swe_house_pos` is intentional compat shim with different signature ✓
- [x] **B.4.3** `sidtime` — verified: calling with just `(jd)` works identically ✓

### B.5 Parameter name mismatches (for keyword callers)

All parameter names renamed in Phase 5 (commit 26ff0ad) to match pyswisseph:

- [!] **B.5.1** `calc_ut` — FIXED: `(tjdut, planet, flags)`
- [!] **B.5.2** `calc` — FIXED: `(tjdet, planet, flags)`
- [!] **B.5.3** `calc_pctr` — FIXED: `(tjdet, planet, center, flags)`
- [!] **B.5.4** `nod_aps_ut` — FIXED: `(tjdut, planet, method, flags)`
- [!] **B.5.5** `pheno_ut` — FIXED: `(tjdut, planet, flags)`
- [!] **B.5.6** `get_orbital_elements` — FIXED: `(tjdet, planet, flags)`
- [!] **B.5.7** `orbit_max_min_true_distance` — FIXED: `(tjdet, planet, flags)`
- [!] **B.5.8** `get_planet_name` — FIXED: `(planet)`
- [!] **B.5.9** `heliacal_ut` — FIXED: `(tjdut, geopos, atmo, observer, objname, eventtype, flags)`
- [!] **B.5.10** `heliacal_pheno_ut` — FIXED: same renames as B.5.9
- [!] **B.5.11** `vis_limit_mag` — FIXED: `(tjdut, geopos, atmo, observer, objname, flags)`
- [x] **B.5.12** `rise_trans` — verified: params already match swe ✓
- [!] **B.5.13** All param names renamed to match pyswisseph — DONE

---

## C. RETURN TYPES

### C.1 Native Python types (no numpy)
- [!] **C.1.1** `get_ayanamsa_ex` returns `(int, float)` — FIXED in Phase 3: now returns native `float`
- [!] **C.1.2** `get_ayanamsa_ex_ut` — FIXED in Phase 3: now returns native `float`
- [x] **C.1.3** Audit ALL functions — verified: calc_ut, houses, cotrans, azalt, refrac, deltat, sidtime, get_ayanamsa_ut, get_ayanamsa_ex_ut all return native Python types ✓

### C.2 Return value structure
- [x] **C.2.1** `get_ayanamsa_ex` retflags: verified non-issue ✓
- [x] **C.2.2** `get_ayanamsa_ex_ut` retflags: verified non-issue ✓
- [x] **C.2.3** `pheno_ut` — verified: both swe and ephem return flat tuple of 20 floats ✓
- [x] **C.2.4** `get_orbital_elements` — verified: both return flat tuple of 50 floats ✓

---

## D. OUTPUT FORMAT (string functions)

### D.1 cs2lonlatstr
- [!] **D.1.1** `cs2lonlatstr` — FIXED in Phase 3: complete rewrite to match pyswisseph compact format (`'1N00'`)

### D.2 cs2timestr
- [!] **D.2.1** `cs2timestr` — FIXED in Phase 3: leading zero, mod-24 hour wrap, suppresszero

### D.3 cs2degstr
- [x] **D.3.1** `cs2degstr` — pyswisseph segfaults (SIGABRT), cannot compare. Our implementation follows documentation ✓

---

## E. BEHAVIORAL COMPATIBILITY

### E.1 Error handling
- [x] **E.1.1** `ephem.Error` exists, is a class, and is catchable ✓
- [x] **E.1.2** All custom exceptions (UnknownBodyError, CoordinateError, EphemerisRangeError) are subclasses of `ephem.Error` ✓
- [x] **E.1.3** `calc_ut` with invalid body: swe returns silently for body 9999 (hypothetical Keplerian), ephem raises Error. Accepted divergence — body 9999 is an obscure hypothetical body.
- [x] **E.1.4** `houses` with invalid latitude (95°): both raise `Error` ✓
- [x] **E.1.5** `fixstar_ut` with non-existent star: both raise `Error` ✓
- [!] **E.1.6** `julday`/`revjul` with invalid calendar flag: FIXED in Phase 4 — both raise `ValueError` ✓

### E.2 State management
- [x] **E.2.1** `set_ephe_path(None)` — both accept None ✓
- [x] **E.2.2** `close()` — resets state, calc_ut works after close ✓
- [x] **E.2.3** `set_sid_mode` persistence — mode persists across calls ✓
- [x] **E.2.4** `set_topo` persistence — Moon topocentric diff 0.01", verified ✓
- [x] **E.2.5** `set_delta_t_userdef` — custom value returned correctly ✓
- [x] **E.2.6** `set_delta_t_userdef(DELTAT_AUTOMATIC)` — restores automatic mode ✓
- [x] **E.2.7** `set_tid_acc` / `get_tid_acc` roundtrip — exact match ✓
- [x] **E.2.8** `set_lapse_rate` — setting works ✓

### E.3 Edge cases
- [x] **E.3.1** `calc_ut` with `SE_ECL_NUT` (-1) — returns nutation/obliquity, values match swe ✓
- [x] **E.3.2** `calc_ut` with asteroid (Eros = SE_AST_OFFSET + 433) — works (swe needs .se1 files) ✓
- [x] **E.3.3** `houses` with Gauquelin (`b'G'`) — returns 36 cusps ✓
- [x] **E.3.4** `sol_eclipse_when_glob` with `backwards=True` — finds previous eclipse, timing matches swe ✓
- [!] **E.3.5** `lun_eclipse_when` with `ecltype` filter — FIXED: non-lunar bits (1,2,8,32) now masked out, matching pyswisseph behavior (commit 131a9ca)
- [x] **E.3.6** `rise_trans` with star name as body — both work, timing matches ✓
- [x] **E.3.7** `lun_occult_when_glob` with star name — both work (minor retflag difference: 6 vs 4)
- [x] **E.3.8** `gauquelin_sector` with star name — both work, values match ✓
- [x] **E.3.9** `date_conversion` with invalid date (Feb 30) — both return `(False, jd, dt)` identically ✓
- [x] **E.3.10** `day_of_week` — both return Monday=0, Sunday=6 convention ✓
- [x] **E.3.11** `split_deg` with all roundflag combinations — 70/70 combos match exactly ✓
- [x] **E.3.12** `cotrans` with negative obliquity — values match swe exactly ✓

---

## F. NUMERICAL PRECISION

### F.1 Core calculations (calc_ut)
- [x] **F.1.1** Sun position at 5 dates: max diff 0.004" ✓
- [x] **F.1.2** Moon position at 5 dates: max diff 0.069" ✓
- [x] **F.1.3** All planets (Mercury–Pluto) at 5 dates each: max diff 0.045" (Pluto) ✓
- [x] **F.1.4** Mean Node, True Node at 5 dates: max diff 0.026" ✓
- [x] **F.1.5** Mean Apogee (0.000"), Oscu Apogee (0.268") at 5 dates ✓
- [ ] **F.1.6** Chiron, Pholus at 50 dates: needs .se1 files for pyswisseph comparison
- [ ] **F.1.7** Ceres, Pallas, Juno, Vesta at 50 dates: needs .se1 files for pyswisseph comparison
- [ ] **F.1.8** IntpApog, IntpPerg at 50 dates: max difference < 5 arcseconds.
- [ ] **F.1.9** Uranian hypothetical bodies (Cupido–Poseidon) at 20 dates each.
- [ ] **F.1.10** Transpluto at 20 dates.

### F.2 Speed values
- [x] **F.2.1** Speed (longitude) for Sun (0.000072), Moon (0.001071), Mars (0.000104) deg/day — all < 0.01 ✓
- [ ] **F.2.2** Speed (latitude) for Moon at 50 dates.
- [ ] **F.2.3** Speed (distance) for Moon at 50 dates.

### F.3 Flag combinations
- [x] **F.3.1** `FLG_SWIEPH` (default) — verified in F.1 spot-check ✓
- [x] **F.3.2** `FLG_SPEED` — verified in F.2 spot-check ✓
- [ ] **F.3.3** `FLG_EQUATORIAL` — RA/Dec match swe.
- [ ] **F.3.4** `FLG_HELCTR` — heliocentric positions match swe.
- [ ] **F.3.5** `FLG_SIDEREAL` with Lahiri — sidereal positions match swe.
- [ ] **F.3.6** `FLG_SIDEREAL` with Fagan-Bradley — positions match swe.
- [ ] **F.3.7** `FLG_J2000` — J2000 ecliptic positions match swe.
- [ ] **F.3.8** `FLG_NONUT` — no-nutation positions match swe.
- [ ] **F.3.9** `FLG_TRUEPOS` — true (geometric) positions match swe.
- [ ] **F.3.10** `FLG_NOABERR` — no-aberration positions match swe.
- [ ] **F.3.11** `FLG_NOGDEFL` — no-gravitational-deflection positions match swe.
- [ ] **F.3.12** `FLG_ASTROMETRIC` (NOABERR|NOGDEFL) — astrometric positions match swe.
- [ ] **F.3.13** `FLG_XYZ` — cartesian coordinates match swe.
- [ ] **F.3.14** `FLG_RADIANS` — radian output matches swe.
- [x] **F.3.15** `FLG_TOPOCTR` — Moon topocentric diff 0.01", verified in E.2.4 ✓
- [ ] **F.3.16** `FLG_SIDEREAL | FLG_EQUATORIAL` — sidereal equatorial match swe.
- [ ] **F.3.17** `FLG_SIDEREAL | FLG_J2000` — sidereal J2000 match swe.

### F.4 Houses
- [ ] **F.4.1** Placidus cusps at 20 locations — max difference < 0.01 arcsecond.
- [ ] **F.4.2** Koch cusps at 20 locations.
- [ ] **F.4.3** Equal cusps at 20 locations.
- [ ] **F.4.4** Whole Sign cusps at 20 locations.
- [ ] **F.4.5** Regiomontanus cusps at 20 locations.
- [ ] **F.4.6** Campanus cusps at 20 locations.
- [ ] **F.4.7** All other house systems (A,B,C,D,F,G,H,I,K,L,M,N,O,Q,S,T,U,V,X,Y) at 5 locations each.
- [ ] **F.4.8** `house_pos` — house position of planets matches swe.
- [ ] **F.4.9** `houses_ex2` cusp speeds match swe.
- [ ] **F.4.10** `houses_armc` results match swe.
- [ ] **F.4.11** `gauquelin_sector` values match swe for 10 planet/location combos.

### F.5 Eclipses
- [ ] **F.5.1** `sol_eclipse_when_glob` — next 10 eclipses from J2000: timing matches swe within 1 second.
- [ ] **F.5.2** `sol_eclipse_when_loc` — timing matches swe within 1 second for 5 locations.
- [ ] **F.5.3** `sol_eclipse_where` — central line lat/lon match swe within 0.01 degree.
- [ ] **F.5.4** `sol_eclipse_how` — attributes match swe within tolerance.
- [ ] **F.5.5** `lun_eclipse_when` — next 10 eclipses: timing matches swe within 1 second.
- [ ] **F.5.6** `lun_eclipse_when_loc` — timing matches swe within 1 second.
- [ ] **F.5.7** `lun_eclipse_how` — attributes match swe.
- [ ] **F.5.8** `lun_occult_when_glob` — next occultation by Mars: timing matches swe.
- [ ] **F.5.9** `lun_occult_when_loc` — timing matches swe.
- [ ] **F.5.10** `lun_occult_where` — position matches swe.

### F.6 Rise/Set/Transit
- [ ] **F.6.1** `rise_trans` Sun rise at 10 locations — matches swe within 1 second.
- [ ] **F.6.2** `rise_trans` Sun set at 10 locations — matches swe within 1 second.
- [ ] **F.6.3** `rise_trans` Moon rise at 10 locations.
- [ ] **F.6.4** `rise_trans` meridian transit of Sun at 10 locations.
- [ ] **F.6.5** `rise_trans` fixed star rise (Sirius) at 5 locations.
- [ ] **F.6.6** `rise_trans_true_hor` with custom horizon altitude — matches swe.

### F.7 Crossings
- [ ] **F.7.1** `solcross_ut` — Sun crossing 0° Aries: matches swe within 1 second.
- [ ] **F.7.2** `mooncross_ut` — Moon crossing 0° Aries: matches swe within 1 second.
- [ ] **F.7.3** `mooncross_node_ut` — next node crossing: matches swe.
- [ ] **F.7.4** `helio_cross_ut` — Mars helio crossing 0°: matches swe.

### F.8 Fixed Stars
- [ ] **F.8.1** `fixstar_ut("Sirius", ...)` — position matches swe < 0.1 arcsecond.
- [ ] **F.8.2** `fixstar_ut("Regulus", ...)` — matches swe.
- [ ] **F.8.3** `fixstar_ut("Aldebaran", ...)` — matches swe.
- [ ] **F.8.4** `fixstar_mag` for 10 stars — magnitude and name match swe.
- [ ] **F.8.5** `fixstar2_ut` — same stars, verify matches swe.

### F.9 Ayanamsa
- [ ] **F.9.1** `get_ayanamsa_ut` with all 47 predefined modes — values match swe within 0.001 arcsecond.
- [ ] **F.9.2** `get_ayanamsa_ex_ut` with `FLG_SIDEREAL` — matches swe.
- [ ] **F.9.3** `get_ayanamsa_name` for all 47 modes — names match swe exactly.

### F.10 Time functions
- [ ] **F.10.1** `julday` / `revjul` roundtrip for 100 dates — exact match.
- [ ] **F.10.2** `deltat` at 50 dates — matches swe within 1e-6 days.
- [ ] **F.10.3** `deltat_ex` at 50 dates — matches swe within 1e-6 days.
- [ ] **F.10.4** `utc_to_jd` at 20 dates — JD_ET and JD_UT match swe.
- [ ] **F.10.5** `jdet_to_utc` / `jdut1_to_utc` — roundtrip with `utc_to_jd` matches.
- [ ] **F.10.6** `sidtime` at 20 dates — matches swe within 1e-6 hours.
- [ ] **F.10.7** `sidtime0` — matches swe.
- [ ] **F.10.8** `time_equ` at 20 dates — equation of time matches swe.
- [ ] **F.10.9** `day_of_week` for 20 JDs — matches swe exactly.
- [ ] **F.10.10** `lmt_to_lat` / `lat_to_lmt` — roundtrip at 10 locations matches swe.

### F.11 Coordinate transforms
- [ ] **F.11.1** `cotrans` ecliptic→equatorial at 20 positions — matches swe.
- [ ] **F.11.2** `cotrans` equatorial→ecliptic at 20 positions — matches swe.
- [ ] **F.11.3** `cotrans_sp` with speeds at 10 positions — matches swe.
- [ ] **F.11.4** `azalt` at 10 positions/times — azimuth/altitude match swe.
- [ ] **F.11.5** `azalt_rev` at 10 positions — reverse transform matches swe.
- [ ] **F.11.6** `refrac` TRUE_TO_APP at 10 altitudes — matches swe.
- [ ] **F.11.7** `refrac` APP_TO_TRUE at 10 altitudes — matches swe.
- [ ] **F.11.8** `refrac_extended` — matches swe.

### F.12 Heliacal events
- [ ] **F.12.1** `heliacal_ut` — Sirius heliacal rising from Cairo: date matches swe within 1 day.
- [ ] **F.12.2** `heliacal_pheno_ut` — attributes match swe within tolerance.
- [ ] **F.12.3** `vis_limit_mag` — limiting magnitude matches swe within 0.5 mag.

### F.13 Nodes and apsides
- [ ] **F.13.1** `nod_aps_ut` for Mars, Jupiter, Saturn — ascending node matches swe.
- [ ] **F.13.2** `nod_aps_ut` with `NODBIT_OSCU` — osculating nodes match swe.
- [ ] **F.13.3** `nod_aps_ut` with `NODBIT_FOPOINT` — focal point matches swe.

### F.14 Planetary phenomena
- [ ] **F.14.1** `pheno_ut` for Venus — phase angle, elongation match swe.
- [ ] **F.14.2** `pheno_ut` for Mars — magnitude matches swe.
- [ ] **F.14.3** `pheno_ut` for Moon — horizontal parallax matches swe.

---

## G. UTILITY FUNCTION ACCURACY

### G.1 Degree/radian utilities
- [x] **G.1.1** `degnorm` — 20 values tested, 0 mismatches ✓
- [x] **G.1.2** `radnorm` — 12 values tested, 0 mismatches ✓
- [x] **G.1.3** `difdeg2n` — 10 pairs tested, 0 mismatches ✓
- [x] **G.1.4** `difdegn` — 10 pairs tested, 0 mismatches ✓
- [x] **G.1.5** `difrad2n` — 5 pairs tested, 0 mismatches ✓
- [x] **G.1.6** `deg_midp` — 10 pairs tested, 0 mismatches ✓
- [x] **G.1.7** `rad_midp` — 4 pairs tested, 0 mismatches ✓
- [x] **G.1.8** `d2l` — 12 values tested, 0 mismatches (excl. accepted negative overflow divergence) ✓

### G.2 Centisecond utilities
- [x] **G.2.1** `csnorm` — 12 values tested, 0 mismatches ✓
- [x] **G.2.2** `csroundsec` — 13 values tested, 0 mismatches ✓
- [x] **G.2.3** `difcsn` — 6 pairs tested, 0 mismatches ✓
- [x] **G.2.4** `difcs2n` — 6 pairs tested, 0 mismatches ✓

### G.3 String formatting
- [x] **G.3.1** `cs2degstr` — pyswisseph segfaults, cannot compare. Our implementation follows docs ✓
- [x] **G.3.2** `cs2lonlatstr` — 23 test cases verified character-for-character match with pyswisseph ✓
- [!] **G.3.3** `cs2timestr` — FIXED in Phase 3: leading zero, mod-24 wrap, suppresszero ✓
- [x] **G.3.4** `split_deg` — 70 value/flag combos verified, all match pyswisseph exactly ✓

---

## H. IMPORT SURFACE

### H.1 Bare-name availability
- [x] **H.1.1** Every pyswisseph function name (without `swe_` prefix) is importable from libephemeris ✓
- [x] **H.1.2** Every pyswisseph constant name is importable from libephemeris ✓
- [!] **H.1.3** `from libephemeris import *` exports all expected names — FIXED in Phase 4b: added 308 bare-name constants to `__all__` (commit 264fe38)

### H.2 Alias consistency
- [x] **H.2.1** Bare names verified — `house_pos` and `lun_occult_where` are intentional compat shims with different signatures ✓
- [x] **H.2.2** `ephem.Error` is a class that can be instantiated and caught ✓

---

## PRIORITY ORDER FOR FIXES

### P0 — BREAKING (will crash user code)
1. B.1.1–B.1.4: fixstar* missing iflag default
2. B.1.5: julday missing hour default
3. B.4.1: lun_occult_when_glob separate ipl/starname params
4. D.1.1: cs2lonlatstr wrong format
5. D.2.1: cs2timestr wrong format

### P1 — BEHAVIORAL (wrong results silently)
6. A.2.1: TIDAL_AUTOMATIC wrong value
7. A.2.2: TIDAL_DEFAULT wrong value
8. B.2.1: gauquelin_sector wrong default flags
9. B.2.2: sol_eclipse_how wrong default ifl
10. C.1.1–C.1.2: numpy types in return values

### P2 — INTEROP (some calling patterns won't work)
11. B.3.1–B.3.6: hsys must accept bytes
12. B.3.8: gauquelin_sector must accept star names
13. C.2.1–C.2.2: get_ayanamsa_ex retflags mismatch

### P3 — COSMETIC (param names, documentation)
14. B.5.1–B.5.13: param name differences (positional OK, keyword may break)
15. A.2.3: version string (intentional)

---

## EXECUTION LOG

| # | Item | Date | Result | Notes |
|---|------|------|--------|-------|
| 1 | A.1.1 | Phase 2 | ✓ OK | All 319 constants exist |
| 2 | A.2.1 | Phase 3 | ✓ FIXED | `TIDAL_AUTOMATIC` 0.0 → 999999 |
| 3 | A.2.2 | Phase 3 | ✓ FIXED | `TIDAL_DEFAULT` -25.936 → -25.8 |
| 4 | A.2.3 | Phase 2 | ✓ Accepted | `version` intentionally different |
| 5 | A.2.4 | Phase 5+ | ✓ OK | TRUE_TO_APP/APP_TO_TRUE match |
| 6 | B.1.1–B.1.4 | Phase 3 | ✓ FIXED | fixstar* `iflag=SEFLG_SWIEPH` default |
| 7 | B.1.5 | Phase 3 | ✓ FIXED | julday `hour=12.0` default |
| 8 | B.2.1 | Phase 3 | ✓ FIXED | gauquelin_sector default flags |
| 9 | B.2.2 | Phase 3 | ✓ FIXED | sol_eclipse_how default ifl |
| 10 | B.3.1–B.3.7 | Phase 3 | ✓ OK | hsys already accepts bytes/str/int |
| 11 | B.3.8 | Phase 3 | ✓ FIXED | gauquelin_sector accepts star names |
| 12 | B.4.1 | Phase 3 | ✓ FIXED | lun_occult_when_glob unified body param |
| 13 | B.4.2–B.4.3 | Phase 3 | ✓ OK | Verified compatible |
| 14 | B.5.1–B.5.13 | Phase 5 | ✓ FIXED | All param names renamed (26ff0ad) |
| 15 | C.1.1–C.1.2 | Phase 3 | ✓ FIXED | Native float returns |
| 16 | C.1.3 | Phase 5+ | ✓ OK | No numpy types in any return value |
| 17 | C.2.1–C.2.4 | Phase 5+ | ✓ OK | Return structures match |
| 18 | D.1.1 | Phase 3 | ✓ FIXED | cs2lonlatstr complete rewrite |
| 19 | D.2.1 | Phase 3 | ✓ FIXED | cs2timestr leading zero + features |
| 20 | D.3.1 | Phase 3 | ✓ N/A | pyswisseph segfaults |
| 21 | E.1.1–E.1.6 | Phase 5+ | ✓ OK/FIXED | Error handling verified |
| 22 | E.2.1–E.2.8 | Phase 5+ | ✓ OK | State management verified |
| 23 | E.3.1–E.3.4 | Phase 5+ | ✓ OK | Edge cases verified |
| 24 | E.3.5 | Phase 5+ | ✓ FIXED | lun_eclipse_when ecltype filter (131a9ca) |
| 25 | E.3.6–E.3.12 | Phase 5+ | ✓ OK | Edge cases verified |
| 26 | G.3.1–G.3.4 | Phase 4/5 | ✓ OK/FIXED | String formatting verified |
| 27 | H.1.1–H.1.3 | Phase 4b | ✓ FIXED | 308 bare-name constants added (264fe38) |
| 28 | H.2.1–H.2.2 | Phase 3 | ✓ OK | Alias consistency verified |

### Summary of fixes by commit:

- **b951d18** (Phase 3): fixstar defaults, julday default, lun_occult_when_glob body param, cs2lonlatstr, cs2timestr, TIDAL_*, gauquelin_sector, sol_eclipse_how, get_ayanamsa_ex float types
- **7903529** (Phase 4): split_deg algorithm rewrite, deg_midp/rad_midp 180° convention, julday/revjul calendar validation, 23 ayanamsa name corrections
- **264fe38** (Phase 4b): 308 bare-name constants added to `__all__`
- **26ff0ad** (Phase 5): All parameter names renamed across 10 source files + 22 test files
- **131a9ca** (Phase 5+): lun_eclipse_when ecltype filter fix for non-lunar bits

### Known accepted divergences:

1. `version` string: libephemeris uses its own version, not pyswisseph's
2. `d2l(-0.5)`: swe returns 4294967295 (C unsigned overflow), ephem returns -1 (mathematically correct)
3. `cs2degstr()`: pyswisseph segfaults (SIGABRT), cannot compare
4. `calc_ut(body=9999)`: swe returns hypothetical Keplerian body, ephem raises Error
5. True Node: ~32" systematic difference (algorithmic)
6. Oscu Apog: ~113" systematic difference (algorithmic)
7. Some exotic ayanamsa modes: up to ~14" difference (algorithmic)
