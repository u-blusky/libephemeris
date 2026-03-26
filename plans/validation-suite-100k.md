# Validation Suite: 100K+ Test Checklist

Comprehensive one-time validation checklist. Each section generates thousands
of individual test cases. Run with a script that iterates over the parameter
space and marks pass/fail.

**Target:** 100,000+ individual comparisons across all backends, bodies, dates, and flag combinations.

---

## 1. Position Accuracy: All Backends x All Bodies x Dense Date Sampling

### 1.1 Skyfield vs Horizons (geocentric, ecliptic of date)

For each body in `[0,1,2,3,4,5,6,7,8,9,14,15,17,18,19,20]` (16 bodies):
For each date in 1000 uniform random JDs in [2415020.5, 2488069.5] (1900-2100):
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` Skyfield vs Horizons: lon diff < 0.001"
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` Skyfield vs Horizons: lat diff < 0.001"
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` Skyfield vs Horizons: dist diff < 1e-8 AU

**Subtotal: 16 x 1000 x 3 = 48,000 checks**

### 1.2 Skyfield vs LEB2 (geocentric, ecliptic of date)

For each body in LEB2 core `[0,1,2,3,4,5,6,7,8,9,10,11,12,14]` (14 bodies):
For each date in 1000 uniform random JDs in [2396759, 2506330] (1850-2150):
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` Skyfield vs LEB2: lon diff < 0.001"
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` Skyfield vs LEB2: lat diff < 0.001"
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` Skyfield vs LEB2: dist diff < 1e-8 AU

**Subtotal: 14 x 1000 x 3 = 42,000 checks**

### 1.3 LEB2 vs Horizons (cross-validation)

For each body in `[0,1,2,3,4,5,6,7,8,9,14,15,17]` (13 bodies):
For each date in 500 random JDs in [2415020.5, 2488069.5]:
- [ ] lon diff < 0.001"
- [ ] lat diff < 0.001"

**Subtotal: 13 x 500 x 2 = 13,000 checks**

---

## 2. Flag Combinations

### 2.1 All flags x all backends x 10 bodies x 100 dates

Flags: `SEFLG_SPEED`, `SEFLG_SPEED|SEFLG_SIDEREAL`, `SEFLG_SPEED|SEFLG_EQUATORIAL`,
`SEFLG_SPEED|SEFLG_J2000`, `SEFLG_SPEED|SEFLG_NOABERR`, `SEFLG_SPEED|SEFLG_NOGDEFL`,
`SEFLG_SPEED|SEFLG_TRUEPOS`, `SEFLG_SPEED|SEFLG_HELCTR`, `SEFLG_SPEED|SEFLG_BARYCTR`,
`SEFLG_SPEED|SEFLG_NONUT`, `SEFLG_SPEED|SEFLG_ICRS`, `SEFLG_SPEED|SEFLG_XYZ`,
`SEFLG_SPEED|SEFLG_RADIANS` (13 flag combos)

Bodies: Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Pluto, Earth, Chiron (10)
Dates: 100 random JDs in 1900-2100
Backends: Skyfield, LEB2, Horizons (3)

For each (backend, flag, body, date):
- [ ] No crash/exception
- [ ] Return tuple has 6 elements
- [ ] lon/RA is finite and in valid range
- [ ] lat/dec is finite and in valid range [-90, 90]
- [ ] dist is positive
- [ ] speed values are finite

**Subtotal: 3 x 13 x 10 x 100 x 6 = 234,000 checks**

### 2.2 Multi-flag combinations (combinatorial)

For each pair of compatible flags from the above list (78 pairs):
For each of 5 bodies x 20 dates:
- [ ] No crash
- [ ] Valid output

**Subtotal: 78 x 5 x 20 = 7,800 checks**

---

## 3. Velocity Accuracy

### 3.1 Speed vs finite difference (all backends)

For each body in [0,1,2,3,4,5,6,7,8,9,14] (11 bodies):
For each date in 200 random JDs:
For each backend (Skyfield, LEB2, Horizons):
- [ ] Compute pos at jd and jd+0.001
- [ ] Numerical speed = (pos2 - pos1) / 0.001
- [ ] Speed from SEFLG_SPEED agrees within 0.01 deg/day

**Subtotal: 11 x 200 x 3 = 6,600 checks**

---

## 4. House Calculations

### 4.1 All house systems x multiple locations x dates

House systems: P, K, O, R, C, E, W, X, M, H, T, B, G, I, i, A, L, N, Q, Y, F, U, D, J (24)
Locations: (0,0), (12.5,41.9), (139.7,35.7), (-73.9,40.7), (0,66.5), (0,-66.5) (6)
Dates: 50 random JDs in 1900-2100

For each (system, location, date):
- [ ] `swe.houses()` returns 13 cusps (cusps[1]-cusps[12])
- [ ] `ascmc[0]` (ASC) is in [0, 360)
- [ ] `ascmc[1]` (MC) is in [0, 360)
- [ ] All cusps are in [0, 360)
- [ ] No crash

**Subtotal: 24 x 6 x 50 x 5 = 36,000 checks**

---

## 5. Sidereal Modes

### 5.1 All 43 ayanamsha modes

For each mode 0-42:
For each date in 20 random JDs:
For each body in [0, 1, 4] (Sun, Moon, Mars):
- [ ] `set_sid_mode(mode)` + `calc_ut(jd, body, SEFLG_SPEED|SEFLG_SIDEREAL)` no crash
- [ ] Sidereal lon differs from tropical lon
- [ ] Sidereal lon is in [0, 360)

**Subtotal: 43 x 20 x 3 x 3 = 7,740 checks**

---

## 6. Edge Cases

### 6.1 Boundary dates

For each body in [0,1,2,3,4,5,6,7,8,9,14]:
- [ ] calc_ut at tier start (JD 2396758.5) — no crash
- [ ] calc_ut at tier end (JD 2506331.5) — no crash
- [ ] calc_ut at J2000.0 (JD 2451545.0) — no crash
- [ ] calc_ut at year 1900 Jan 1 — no crash
- [ ] calc_ut at year 2100 Dec 31 — no crash

**Subtotal: 11 x 5 = 55 checks**

### 6.2 Out-of-range dates

For each body in [0,1,2,14]:
- [ ] calc_ut at JD 1000000 — raises EphemerisRangeError or falls back
- [ ] calc_ut at JD 5000000 — raises EphemerisRangeError or falls back
- [ ] calc_ut at JD 0 — handled gracefully

**Subtotal: 4 x 3 = 12 checks**

### 6.3 Invalid inputs

- [ ] calc_ut with body_id=999 — raises or returns error
- [ ] calc_ut with body_id=-1 — handled
- [ ] calc_ut with NaN jd — handled
- [ ] calc_ut with Inf jd — handled
- [ ] houses with invalid system byte — handled
- [ ] houses with lat=91 — handled
- [ ] set_calc_mode("invalid") — raises ValueError
- [ ] set_sid_mode(999) — handled

**Subtotal: 8 checks**

---

## 7. LEB2 Compression Integrity

### 7.1 Round-trip compression

For each body in all 31 bodies:
- [ ] compress_body() -> decompress_body() round-trip preserves coefficients within truncation tolerance
- [ ] Decompressed size matches original
- [ ] Compressed size < original size

**Subtotal: 31 x 3 = 93 checks**

### 7.2 Reader consistency

For each body in LEB2 core (14 bodies):
For each date in 500 random JDs:
- [ ] LEB2Reader.eval_body() matches LEBReader.eval_body() within 0.001"
- [ ] eval_nutation() identical (no compression)
- [ ] delta_t() identical (no compression)

**Subtotal: 14 x 500 + 500 + 500 = 8,000 checks**

### 7.3 CompositeLEBReader dispatch

For each body in all 31 bodies (if companion files present):
For each date in 100 random JDs:
- [ ] CompositeLEBReader dispatches to correct file
- [ ] Result matches individual reader

**Subtotal: 31 x 100 = 3,100 checks**

### 7.4 open_leb() factory

- [ ] LEB1 file opens as LEBReader
- [ ] LEB2 file opens as LEB2Reader
- [ ] Invalid magic raises ValueError
- [ ] Non-existent file raises FileNotFoundError

**Subtotal: 4 checks**

---

## 8. Horizons Backend Specific

### 8.1 Body command mapping

For each body in _HORIZONS_COMMAND (17 bodies):
- [ ] Command string is valid Horizons syntax
- [ ] fetch_state_vector returns 6 components
- [ ] Positions are finite and in plausible AU range

**Subtotal: 17 x 3 = 51 checks**

### 8.2 Cache behavior

- [ ] Same (jd, body) returns identical results (cache hit)
- [ ] Different jd returns different results
- [ ] Cache size doesn't exceed max_cache_size
- [ ] clear_cache() empties the cache
- [ ] shutdown() is idempotent

**Subtotal: 5 checks**

### 8.3 Analytical bodies (no HTTP)

For each body in [10, 12] (Mean Node, Mean Apogee):
For each date in 100 random JDs:
- [ ] horizons_calc_ut() returns result without HTTP call
- [ ] Result matches Skyfield within 0.01"

**Subtotal: 2 x 100 x 2 = 400 checks**

### 8.4 Unsupported bodies fallback

For each body in [11, 13, 21, 22] (True Node, Oscu Apogee, Interp Apogee/Perigee):
- [ ] horizons_calc_ut() raises KeyError
- [ ] swe_calc_ut in auto mode falls through to Skyfield
- [ ] Result is correct (matches pure Skyfield)

**Subtotal: 4 x 3 = 12 checks**

### 8.5 Error handling

- [ ] Timeout: ConnectionError after retries with helpful message
- [ ] Invalid body: KeyError triggers fallback
- [ ] SEFLG_TOPOCTR: KeyError triggers fallback
- [ ] Sun heliocentric: returns (0,0,0,0,0,0)

**Subtotal: 4 checks**

---

## 9. Ecliptic/Equatorial Coordinate Transforms

### 9.1 cotrans round-trip

For each obliquity in [23.0, 23.4393, 24.0]:
For each (lon, lat) in 100 random points:
- [ ] ecl->eq->ecl round-trip: lon matches within 1e-10 deg
- [ ] ecl->eq->ecl round-trip: lat matches within 1e-10 deg

**Subtotal: 3 x 100 x 2 = 600 checks**

---

## 10. Fixed Stars

### 10.1 Star catalog consistency

For each star in first 50 Hipparcos stars:
- [ ] fixstar_ut() returns valid position (finite lon, lat)
- [ ] fixstar_mag() returns valid magnitude
- [ ] Proper motion applied over 100 years gives plausible shift

**Subtotal: 50 x 3 = 150 checks**

---

## 11. Julian Day Conversions

### 11.1 julday/revjul round-trip

For each of 1000 random dates (year -5000 to +5000):
- [ ] revjul(julday(y,m,d,h)) == (y,m,d,h) within 1e-6

**Subtotal: 1000 checks**

### 11.2 Delta-T consistency

For each of 200 JDs in [1900, 2100]:
- [ ] swe_deltat(jd) returns positive value
- [ ] swe_deltat(jd) < 0.01 days (~14 min, reasonable for modern era)
- [ ] Value is smooth (no discontinuities between adjacent JDs)

**Subtotal: 200 x 3 = 600 checks**

---

## 12. Eclipse Computations

### 12.1 Solar eclipse search

For each year in [2000, 2005, 2010, 2015, 2020, 2025]:
- [ ] sol_eclipse_when_glob() finds at least 1 eclipse
- [ ] Eclipse JD is within the year
- [ ] Eclipse type is valid (1-4)

**Subtotal: 6 x 3 = 18 checks**

### 12.2 Lunar eclipse search

For each year in [2000, 2005, 2010, 2015, 2020, 2025]:
- [ ] lun_eclipse_when() finds at least 1 eclipse
- [ ] Eclipse JD is within the year

**Subtotal: 6 x 2 = 12 checks**

---

## 13. Performance Benchmarks

### 13.1 Backend speed verification

- [ ] LEB eval_body: < 20 us/call (median over 10000 calls)
- [ ] Skyfield calc_ut: < 500 us/call (median over 1000 calls)
- [ ] LEB2 eval_body (after first decompress): < 20 us/call
- [ ] Horizons (cached): < 1 ms/call
- [ ] Houses: < 200 us/call

**Subtotal: 5 checks**

---

## 14. Cross-Backend Consistency Matrix

### 14.1 Triple comparison (Skyfield vs LEB2 vs Horizons)

For each body in [0,1,2,3,4,5,6,7,8,9,14] (11 bodies):
For each date in 100 random JDs:
- [ ] |Skyfield - LEB2| < 0.005" (lon)
- [ ] |Skyfield - Horizons| < 0.005" (lon)
- [ ] |LEB2 - Horizons| < 0.005" (lon)
- [ ] |Skyfield - LEB2| < 0.005" (lat)
- [ ] |Skyfield - Horizons| < 0.005" (lat)
- [ ] |LEB2 - Horizons| < 0.005" (lat)

**Subtotal: 11 x 100 x 6 = 6,600 checks**

---

## Summary

| Section | Checks |
|---------|--------|
| 1. Position accuracy (3 backend pairs) | 103,000 |
| 2. Flag combinations | 241,800 |
| 3. Velocity accuracy | 6,600 |
| 4. House calculations | 36,000 |
| 5. Sidereal modes | 7,740 |
| 6. Edge cases | 75 |
| 7. LEB2 compression | 11,197 |
| 8. Horizons backend | 472 |
| 9. Coordinate transforms | 600 |
| 10. Fixed stars | 150 |
| 11. Julian day / Delta-T | 1,600 |
| 12. Eclipses | 30 |
| 13. Performance | 5 |
| 14. Cross-backend consistency | 6,600 |
| **TOTAL** | **~415,869** |
