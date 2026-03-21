# PySwissEph 1:1 Compatibility — Verification Plan

Comprehensive checklist to guarantee drop-in compatibility with `pyswisseph` 2.10.03.
Each item is a concrete, testable check. Status: `[ ]` = pending, `[x]` = verified OK, `[!]` = discrepancy found and fixed.

---

## A. CONSTANTS (319 in pyswisseph)

### A.1 Missing constants
- [ ] **A.1.1** All 319 pyswisseph constants exist in libephemeris *(verified: 0 missing)*

### A.2 Value mismatches
- [ ] **A.2.1** `TIDAL_AUTOMATIC`: swe=`999999` (int), ephem=`0.0` (float) — must be `999999`
- [ ] **A.2.2** `TIDAL_DEFAULT`: swe=`-25.8`, ephem=`-25.936` — must match swe value `-25.8`
- [ ] **A.2.3** `version`: swe=`'2.10.03'`, ephem=`'0.25.0'` — intentional, document as accepted divergence
- [ ] **A.2.4** `TRUE_TO_APP` / `APP_TO_TRUE`: swe has `TRUE_TO_APP=0, APP_TO_TRUE=1`, verify ephem matches

---

## B. FUNCTION SIGNATURES (101 functions in pyswisseph)

### B.1 Missing default values
- [ ] **B.1.1** `fixstar(star, jd, iflag)` — ephem has NO default for `iflag`, swe defaults to `FLG_SWIEPH` (2)
- [ ] **B.1.2** `fixstar_ut(star, jd, iflag)` — same: needs `iflag=FLG_SWIEPH` default
- [ ] **B.1.3** `fixstar2(star, jd, iflag)` — same
- [ ] **B.1.4** `fixstar2_ut(star, jd, iflag)` — same
- [ ] **B.1.5** `julday(year, month, day, hour)` — ephem has NO default for `hour`, swe defaults to `12.0`

### B.2 Wrong default values
- [ ] **B.2.1** `gauquelin_sector` flags default: ephem=`0`, swe=`FLG_SWIEPH|FLG_TOPOCTR` (32770) — must match swe
- [ ] **B.2.2** `sol_eclipse_how` ifl default: ephem=`0`, swe=`FLG_SWIEPH` (2) — must be 2

### B.3 Parameter type mismatches (hsys: bytes vs int)

pyswisseph accepts `hsys` as `bytes` (e.g. `b'P'`). libephemeris uses `int` (e.g. `ord('P')`).
For 1:1 compat, libephemeris must accept BOTH bytes and int for hsys.

- [ ] **B.3.1** `houses(jd, lat, lon, hsys)` — must accept `b'P'` as hsys
- [ ] **B.3.2** `houses_ex(jd, lat, lon, hsys, flags)` — must accept `b'P'`
- [ ] **B.3.3** `houses_ex2(jd, lat, lon, hsys, flags)` — must accept `b'P'`
- [ ] **B.3.4** `houses_armc(armc, lat, eps, hsys)` — must accept `b'P'`
- [ ] **B.3.5** `houses_armc_ex2(armc, lat, eps, hsys)` — must accept `b'P'`
- [ ] **B.3.6** `house_pos(armc, lat, eps, objcoord, hsys)` — must accept `b'P'` as hsys
- [ ] **B.3.7** `house_name(hsys)` — swe ONLY accepts bytes, NOT int. ephem accepts both. Verify ephem handles both. *(ephem is more permissive — OK)*
- [ ] **B.3.8** `gauquelin_sector` — swe accepts int OR str for body param, ephem only accepts int. Must accept str (star names)

### B.4 Parameter structure mismatches
- [ ] **B.4.1** `lun_occult_when_glob` — ephem has separate `ipl` + `starname` params, swe has single `body` param accepting int or str. Must unify to single `body` param.
- [ ] **B.4.2** `house_pos` — swe signature: `(armc, geolat, eps, objcoord_seq, hsys=b'P')`. Verify ephem accepts sequence for objcoord in 4th position with hsys in 5th.
- [ ] **B.4.3** `sidtime` — ephem has extra optional params `(jd, longitude=0.0, obliquity=None, nutation=None)`, swe has only `(jd)`. Extra params are OK (backward compat), but verify calling with just `(jd)` works identically.

### B.5 Parameter name mismatches (for keyword callers)

Even if positional calling works, keyword callers could break if param names differ.

- [ ] **B.5.1** `calc_ut` — swe: `(tjdut, planet, flags)`, ephem: `(tjd_ut, ipl, iflag)`. Names differ but positional OK. Document.
- [ ] **B.5.2** `calc` — swe: `(tjdet, planet, flags)`, ephem: `(tjd, ipl, iflag)`. Same.
- [ ] **B.5.3** `calc_pctr` — swe: `(tjd, planet, center, flags)`, ephem: `(tjd_ut, ipl, iplctr, iflag)`.
- [ ] **B.5.4** `nod_aps_ut` — swe: `(tjdut, planet, method, flags)`, ephem: `(tjd_ut, ipl, method, iflag)`.
- [ ] **B.5.5** `pheno_ut` — swe: `(tjdut, planet, flags)`, ephem: `(tjd_ut, ipl, iflag)`.
- [ ] **B.5.6** `get_orbital_elements` — swe: `(tjdet, planet, flags)`, ephem: `(tjd_et, ipl, iflag)`.
- [ ] **B.5.7** `orbit_max_min_true_distance` — swe: `(tjdet, planet, flags)`, ephem: `(tjd_ut, ipl, iflag)`.
- [ ] **B.5.8** `get_planet_name` — swe: `(planet)`, ephem: `(planet_id)`.
- [ ] **B.5.9** `heliacal_ut` — swe: `(tjdut, geopos, atmo, observer, objname, eventtype, flags)`, ephem: `(jd_start, geopos, datm, dobs, object_name, event_type, hel_flag)`. All names differ.
- [ ] **B.5.10** `heliacal_pheno_ut` — same name differences as heliacal_ut.
- [ ] **B.5.11** `vis_limit_mag` — swe: `(tjdut, geopos, atmo, observer, objname, flags)`, ephem: `(jd, geopos, atmo, observer, objname, flags)`.
- [ ] **B.5.12** `rise_trans` — swe: `(tjdut, body, rsmi, geopos, ...)`, ephem: `(tjdut, body, rsmi, geopos, ...)`. Verify match.
- [ ] **B.5.13** Document all param name differences and decide: accept divergence (positional-only is fine) or rename to match.

---

## C. RETURN TYPES

### C.1 Native Python types (no numpy)
- [ ] **C.1.1** `get_ayanamsa_ex` returns `(int, float)` — ephem returns `(int, np.float64)`. Must return native `float`.
- [ ] **C.1.2** `get_ayanamsa_ex_ut` — same: returns `np.float64` instead of native `float`.
- [ ] **C.1.3** Audit ALL functions that might return numpy types instead of native Python types. Check: `calc_ut`, `calc`, `houses`, `azalt`, `cotrans`, `refrac`, etc.

### C.2 Return value structure
- [ ] **C.2.1** `get_ayanamsa_ex` retflags: swe returns `2` (FLG_SWIEPH), ephem returns `0`. Investigate if this is meaningful.
- [ ] **C.2.2** `get_ayanamsa_ex_ut` retflags: same issue.
- [ ] **C.2.3** `pheno_ut` — swe returns flat tuple of 20 floats (not wrapped in outer tuple). Verify ephem matches.
- [ ] **C.2.4** `get_orbital_elements` — swe returns flat tuple of 50 floats. Verify ephem matches.

---

## D. OUTPUT FORMAT (string functions)

### D.1 cs2lonlatstr
- [ ] **D.1.1** `cs2lonlatstr(360000, b'N', b'S')` — swe: `'1N00'`, ephem: `'  1° 0' 0" b'N''`. FORMAT COMPLETELY WRONG. Must match swe format exactly.

### D.2 cs2timestr
- [ ] **D.2.1** `cs2timestr(360000, b':', False)` — swe: `'01:00:00'`, ephem: `' 1:00:00'`. Leading zero vs leading space. Must match swe.

### D.3 cs2degstr
- [ ] **D.3.1** Compare `cs2degstr` output for several values — verify format matches swe exactly.

---

## E. BEHAVIORAL COMPATIBILITY

### E.1 Error handling
- [ ] **E.1.1** `swe.Error` exception class exists and is catchable — verify `ephem.Error` exists and works the same way.
- [ ] **E.1.2** Functions that raise `swe.Error` in pyswisseph should raise `ephem.Error` in libephemeris (not `ValueError`, `RuntimeError`, etc.).
- [ ] **E.1.3** Test: `calc_ut` with invalid body number — both should raise `Error`.
- [ ] **E.1.4** Test: `houses` with invalid latitude (e.g., 95°) — compare error behavior.
- [ ] **E.1.5** Test: `fixstar_ut` with non-existent star name — compare error behavior.
- [ ] **E.1.6** Test: `julday` / `revjul` with invalid calendar flag — swe raises `ValueError`, verify ephem does too.

### E.2 State management
- [ ] **E.2.1** `set_ephe_path(None)` — swe accepts None (resets to empty), verify ephem does too.
- [ ] **E.2.2** `close()` — verify it resets all internal state like swe does.
- [ ] **E.2.3** `set_sid_mode` persistence — set mode, call `get_ayanamsa_ut`, verify mode persists.
- [ ] **E.2.4** `set_topo` persistence — set topo, call `calc_ut` with `FLG_TOPOCTR`, verify topo is used.
- [ ] **E.2.5** `set_delta_t_userdef` — set custom value, call `deltat`, verify custom value is returned.
- [ ] **E.2.6** `set_delta_t_userdef(DELTAT_AUTOMATIC)` — verify it restores automatic mode.
- [ ] **E.2.7** `set_tid_acc` / `get_tid_acc` roundtrip — set value, get it back, verify match.
- [ ] **E.2.8** `set_lapse_rate` / reset — verify setting and resetting works.

### E.3 Edge cases
- [ ] **E.3.1** `calc_ut` with body `SE_ECL_NUT` (-1) — must return nutation/obliquity, not error.
- [ ] **E.3.2** `calc_ut` with asteroid numbers (e.g., `SE_AST_OFFSET + 433` for Eros) — verify works.
- [ ] **E.3.3** `houses` with Gauquelin system (`ord('G')` or `b'G'`) — must return 36 cusps, not 12.
- [ ] **E.3.4** `sol_eclipse_when_glob` with `backwards=True` — verify finds previous eclipse.
- [ ] **E.3.5** `lun_eclipse_when` with `ecltype=ECL_TOTAL` filter — verify only total eclipses returned.
- [ ] **E.3.6** `rise_trans` with `body` as star name string — verify works like swe.
- [ ] **E.3.7** `lun_occult_when_glob` with star name as body — verify works (currently broken: separate ipl/starname params).
- [ ] **E.3.8** `gauquelin_sector` with star name as body — verify works (currently typed as int only).
- [ ] **E.3.9** `date_conversion` with invalid date (e.g., Feb 30) — swe returns `(False, jd, dt)`, verify ephem matches.
- [ ] **E.3.10** `day_of_week` — swe returns Monday=0, verify ephem matches (not Sunday=0).
- [ ] **E.3.11** `split_deg` with all roundflag combinations — verify output matches swe.
- [ ] **E.3.12** `cotrans` with negative obliquity (ecliptic→equatorial vs equatorial→ecliptic) — verify convention matches swe.

---

## F. NUMERICAL PRECISION

### F.1 Core calculations (calc_ut)
- [ ] **F.1.1** Sun position at 100 random dates: max difference vs swe < 1 arcsecond.
- [ ] **F.1.2** Moon position at 100 random dates: max difference vs swe < 1 arcsecond.
- [ ] **F.1.3** All planets (Mercury–Pluto) at 50 dates each: max difference < 1 arcsecond.
- [ ] **F.1.4** Mean Node, True Node at 50 dates: max difference < 1 arcsecond.
- [ ] **F.1.5** Mean Apogee, Oscu Apogee at 50 dates: max difference < 1 arcsecond.
- [ ] **F.1.6** Chiron, Pholus at 50 dates: max difference < 1 arcsecond.
- [ ] **F.1.7** Ceres, Pallas, Juno, Vesta at 50 dates: max difference < 1 arcsecond.
- [ ] **F.1.8** IntpApog, IntpPerg at 50 dates: max difference < 5 arcseconds.
- [ ] **F.1.9** Uranian hypothetical bodies (Cupido–Poseidon) at 20 dates each.
- [ ] **F.1.10** Transpluto at 20 dates.

### F.2 Speed values
- [ ] **F.2.1** Speed (longitude) for Sun, Moon, planets at 50 dates: max difference < 0.01 deg/day.
- [ ] **F.2.2** Speed (latitude) for Moon at 50 dates.
- [ ] **F.2.3** Speed (distance) for Moon at 50 dates.

### F.3 Flag combinations
- [ ] **F.3.1** `FLG_SWIEPH` (default) — positions match swe.
- [ ] **F.3.2** `FLG_SPEED` — speeds match swe.
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
- [ ] **F.3.15** `FLG_TOPOCTR` — topocentric positions match swe (after `set_topo`).
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
- [ ] **G.1.1** `degnorm` — 50 values including negatives and > 360 — matches swe exactly.
- [ ] **G.1.2** `radnorm` — 50 values — matches swe exactly.
- [ ] **G.1.3** `difdeg2n` — 30 pairs — matches swe exactly.
- [ ] **G.1.4** `difdegn` — 30 pairs — matches swe exactly.
- [ ] **G.1.5** `difrad2n` — 30 pairs — matches swe exactly.
- [ ] **G.1.6** `deg_midp` — 20 pairs — matches swe.
- [ ] **G.1.7** `rad_midp` — 20 pairs — matches swe.
- [ ] **G.1.8** `d2l` — 20 values — matches swe exactly.

### G.2 Centisecond utilities
- [ ] **G.2.1** `csnorm` — 20 values — matches swe exactly.
- [ ] **G.2.2** `csroundsec` — 20 values — matches swe exactly.
- [ ] **G.2.3** `difcsn` — 20 pairs — matches swe exactly.
- [ ] **G.2.4** `difcs2n` — 20 pairs — matches swe exactly.

### G.3 String formatting
- [ ] **G.3.1** `cs2degstr` — 20 values — output matches swe character-for-character.
- [ ] **G.3.2** `cs2lonlatstr` — 20 values with various plus/minus chars — matches swe.
- [ ] **G.3.3** `cs2timestr` — 20 values with various separators — matches swe.
- [ ] **G.3.4** `split_deg` — 30 values with all roundflag combos — matches swe.

---

## H. IMPORT SURFACE

### H.1 Bare-name availability
- [ ] **H.1.1** Every pyswisseph function name (without `swe_` prefix) is importable from libephemeris.
- [ ] **H.1.2** Every pyswisseph constant name is importable from libephemeris.
- [ ] **H.1.3** `from libephemeris import *` exports all expected names.

### H.2 Alias consistency
- [ ] **H.2.1** For each bare name, verify `ephem.foo is ephem.swe_foo` (they must be the same object, not copies).
- [ ] **H.2.2** `ephem.Error` is a class that can be instantiated and caught.

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
| | | | | |
