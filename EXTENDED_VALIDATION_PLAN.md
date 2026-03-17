# Extended Tier LEB — Validation Plan

## Context

The extended tier LEB (`ephemeris_extended.leb`) covers ~10,000 years (JD -105152.5 to 3547272.5, approximately -13200 to +17191 CE) and was generated from DE441 via Skyfield. This plan validates it against two independent references:

- **Skyfield + DE441**: the generation source — any difference is a real LEB Chebyshev error
- **pyswisseph + SE files (DE431-based)**: the compatibility target — differences include a known DE431↔DE441 baseline that grows with distance from J2000

### Available Resources

| Resource | Path | Range | Notes |
|---|---|---|---|
| Extended LEB | `data/leb/ephemeris_extended.leb` | JD -105152.5 to 3547272.5 | 2.8 GB, 31 bodies, generated from DE441 |
| Medium LEB | `data/leb/ephemeris_medium.leb` | JD 2287185.5 to 2688952.5 | 377 MB, 31 bodies, generated from DE440 |
| DE441 BSP | `/Volumes/data/libephemeris/de441.bsp` | Full extended range | 3.1 GB |
| DE440 BSP | `/Volumes/data/libephemeris/de440.bsp` | ~1550-2650 CE | 114 MB |
| SE .se1 files | `swisseph/ephe/` | Full extended range | DE431-based, full planet + moon + asteroid coverage |
| pyswisseph | pip installed | 2.10.03 | `swe.set_ephe_path('swisseph/ephe')` |

### Body Availability in pyswisseph (SE)

| Bodies | Range in SE |
|---|---|
| Sun–Pluto, Earth, Nodes, MeanApog, OscuApog, Ceres–Vesta, Uranians, Transpluto | Full extended range |
| Chiron | ~660–4600 CE (JD ~1967639 to ~3419402) |
| IntpApog / IntpPerg | ~-3000–2900 CE (JD ~625076 to ~2817969) |

### Measured Baseline: SE(DE431) vs Skyfield(DE441)

| Epoch | Sun | Moon | Jupiter | MeanNode |
|---|---|---|---|---|
| -3000 CE | 21" | 217" | 14" | 8.2" |
| 0 CE | 5" | 88" | 0.7" | 0.24" |
| 1900 CE | 0.0006" | 0.009" | 0.027" | 0.11" |
| 2000 CE | 0.001" | 0.0006" | 0.02" | 0.03" |
| 2650 CE | 20" | 283" | 4.5" | 0.41" |

### Setup for All Tests

```python
import libephemeris as ephem
from libephemeris.state import set_calc_mode, set_leb_file, set_ephe_path, set_ephemeris_file
import swisseph as swe

# LEB
set_leb_file('data/leb/ephemeris_extended.leb')

# Skyfield with DE441 (the LEB source)
set_ephe_path('/Volumes/data/libephemeris')
set_ephemeris_file('de441.bsp')

# pyswisseph with full SE files
swe.set_ephe_path('/Users/giacomo/dev/libephemeris/swisseph/ephe')
```

---

## Phase 1: LEB vs Skyfield(DE441) — Chebyshev Accuracy (8 rounds)

Every difference here is a real LEB approximation error.

### Round 1.1 — Full-Range Sweep: All 31 Bodies

**Goal**: Establish the baseline precision profile of the extended LEB across its entire 10,000-year range.

- **Bodies**: All 31 (Sun–Pluto, Earth, Nodes, Apogees, Chiron, Ceres–Vesta, Uranians, Transpluto)
- **Dates**: 200 uniformly distributed dates from JD -105152.5 to 3547272.5
- **Flags**: `SEFLG_SPEED` (default ecliptic of date)
- **Metrics**: For each body compute P50, P95, P99, max of |dLon|, |dLat|, |dDist|, |dLonSpd|
- **Acceptance**: Moon < 0.5", all others < 0.001"
- **Output**: Table per body with statistics, identify any body/date outliers

### Round 1.2 — Moon Deep Dive: 1000 Dates

**Goal**: Moon is the hardest body for Chebyshev (fastest, 4-day segments). Stress-test across the full range.

- **Body**: Moon only
- **Dates**: 1000 uniformly distributed dates, full range
- **Flags**: `SEFLG_SPEED`
- **Metrics**: All 6 components (lon, lat, dist, lon_speed, lat_speed, dist_speed)
- **Additional**: Histogram of errors by millennium (does precision degrade at extremes?)
- **Acceptance**: Position < 0.5", speed < 0.002 d/d

### Round 1.3 — Segment Boundary Continuity

**Goal**: Verify no discontinuities at Chebyshev polynomial edges for all bodies.

- **Bodies**: All 31
- **Method**: At each segment boundary, compare position at (edge - 1sec) and (edge + 1sec) to linear prediction from speed. If the jump exceeds threshold, the polynomials don't connect smoothly.
- **Sample**: 300 boundaries per body (uniformly sampled from all segments)
- **Acceptance**: < 0.001" jump for all bodies (Moon allowed up to 0.001")
- **Additional**: Report max speed discontinuity (|speed_after - speed_before|) as well

### Round 1.4 — Segment Midpoint vs Edge Accuracy

**Goal**: Chebyshev polynomials are most accurate at Chebyshev nodes (near center) and potentially less accurate at segment edges. Verify this.

- **Bodies**: Moon, Mercury, Venus, Mars, Jupiter, TrueNode, OscuApog
- **Method**: For 100 randomly chosen segments per body, evaluate LEB vs Skyfield at segment start, 25%, 50%, 75%, and end. Compare error distribution.
- **Acceptance**: Edge error should not exceed 2× midpoint error systematically
- **Output**: Table showing mean error at each relative position within segment

### Round 1.5 — Velocity: All 3 Components

**Goal**: Velocity precision for all bodies and all 3 speed components (lon_speed, lat_speed, dist_speed).

- **Bodies**: All 31
- **Dates**: 300 uniformly distributed dates, full range
- **Flags**: `SEFLG_SPEED`
- **Metrics**: P50, P95, P99, max for each speed component per body
- **Acceptance**: lon_speed < 0.05 d/d, lat_speed < 0.05 d/d, dist_speed < 1e-4 AU/d
- **Additional**: Flag any body with max/P99 ratio > 100 (indicates isolated spike — investigate if Skyfield artifact vs real LEB issue)

### Round 1.6 — Stress Test: First and Last 100 Years

**Goal**: LEB precision at the very edges of its range, where polynomial fitting may be less constrained.

- **Bodies**: Sun, Moon, Mercury, Mars, Jupiter, Pluto, MeanNode, Cupido
- **Dates**: 50 dates in first 100 years (JD -105152.5 to -68697.5) + 50 in last 100 years (JD 3510817.5 to 3547272.5)
- **Flags**: `SEFLG_SPEED`
- **Metrics**: Same as Round 1.1
- **Additional**: Compare error magnitude at edges vs range center — quantify degradation factor
- **Acceptance**: No body should show > 5× degradation vs center precision

### Round 1.7 — Distance Precision: Geocentric and Heliocentric

**Goal**: Verify distance component (3rd return value) across entire range, both geocentric and heliocentric.

- **Bodies (geocentric)**: All 31
- **Bodies (heliocentric)**: Mercury–Pluto, Chiron, Ceres–Vesta
- **Dates**: 100 dates, full range
- **Metrics**: |dDist| in AU, relative error |dDist|/|dist|
- **Acceptance**: Geocentric < 5e-6 AU, heliocentric < 5e-5 AU

### Round 1.8 — Acceleration Consistency

**Goal**: Verify that the second derivative (acceleration) implied by consecutive speed evaluations is physically smooth — no unphysical jerks.

- **Bodies**: Moon, Mercury, Mars, Jupiter, TrueNode
- **Method**: At 200 dates, compute speed at t-1day, t, t+1day. Compute acceleration = (speed_after - speed_before) / 2. Compare LEB vs Skyfield acceleration.
- **Acceptance**: Acceleration difference < 0.001 d/d² for planets, < 0.01 d/d² for Moon

---

## Phase 2: Three-Way LEB ↔ Skyfield ↔ Swiss Ephemeris (8 rounds)

The gold rule: `|LEB - SE|` must never significantly exceed `|Skyfield - SE|`. If it does, there's a LEB-specific bug beyond the inherent Skyfield↔SE model difference.

Define: `excess = |LEB - SE| - |Sky - SE|`. If excess > threshold, flag as potential bug.

### Round 2.1 — Pipeline A Three-Way: Modern Era (1800-2200 CE)

**Goal**: In the modern era, DE431≈DE441 so SE≈Skyfield. Any LEB↔SE difference should be pure Chebyshev error and should be tiny.

- **Bodies**: Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Earth
- **Dates**: 50 dates from 1800-2200 CE
- **Flags**: `SEFLG_SPEED`
- **Metrics**: For each point: `|LEB-SE|`, `|Sky-SE|`, `excess = |LEB-SE| - |Sky-SE|`
- **Acceptance**: excess < 0.002" for all bodies, excess < 0.01" for Moon
- **Output**: Table with all three values per body per date

### Round 2.2 — Pipeline A Three-Way: Historical Era (-3000 to +3000 CE)

**Goal**: Wider era where DE431↔DE441 baseline grows. Verify LEB doesn't add error on top of the baseline.

- **Bodies**: Same as 2.1
- **Dates**: 30 dates from -3000 to +2900 CE (pyswisseph range)
- **Metrics**: Same triple comparison. For each body/date compute excess.
- **Acceptance**: excess < 0.01" for all planets, excess < 0.1" for Moon (Moon has larger Chebyshev error)

### Round 2.3 — Pipeline B Three-Way: Nodes and Apogees

**Goal**: Pipeline B bodies use analytical formulas (Meeus). Verify LEB matches Skyfield exactly, and both have the same model diff vs SE.

- **Bodies**: MeanNode, TrueNode, MeanApog, OscuApog, IntpApog, IntpPerg
- **Dates**: 30 dates from -3000 to +2900 CE (IntpApog range limit)
- **Flags**: `SEFLG_SPEED`
- **Additional**: Compute the model diff `|Skyfield - SE|` for each body and verify it's stable (doesn't suddenly jump). If `|LEB - SE|` is significantly different from `|Sky - SE|` at the same date, that's a LEB bug.
- **Acceptance**: `|LEB - Sky|` < 0.001" for MeanNode/MeanApog, < 5" for TrueNode

### Round 2.4 — Pipeline C Three-Way: Uranians and Transpluto

**Goal**: These use simple polynomial formulas. All three systems should agree perfectly (within floating-point).

- **Bodies**: Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon, Transpluto
- **Dates**: 30 dates from -3000 to +2900 CE
- **Flags**: `SEFLG_SPEED`
- **Acceptance**: `|LEB - SE|` < 0.1" AND `|LEB - Sky|` < 0.001". Since all use the same formula, larger differences indicate a bug.

### Round 2.5 — Sidereal Three-Way: All Flag Combinations

**Goal**: Verify all 3 sidereal bug fixes work correctly on extended tier by comparing all three systems.

- **Bodies**: Sun, Moon, Mars, MeanNode, TrueNode, OscuApog, Cupido
- **Dates**: 10 dates from 1800-2200 CE (where SE≈DE441, maximum sensitivity to sidereal bugs)
- **Flags**: All 8 sidereal combinations:
  - `SEFLG_SIDEREAL`
  - `SEFLG_SIDEREAL | SEFLG_EQUATORIAL`
  - `SEFLG_SIDEREAL | SEFLG_J2000`
  - `SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_J2000`
  - (same 4 without SIDEREAL as control)
- **Ayanamsha**: Lahiri (1) and Fagan-Bradley (0) — test at least 2 modes
- **Metrics**: Three-way comparison for each. Verify:
  - `|LEB_sid - SE_sid|` ≈ `|Sky_sid - SE_sid|` (no LEB sidereal bug)
  - `|LEB_sid - LEB_trop| - |SE_sid - SE_trop|` < 0.001" (sidereal offset applied identically)
- **Acceptance**: excess < 0.002" for all bodies in modern era

### Round 2.6 — Sidereal Three-Way: Extended Era

**Goal**: Same as 2.5 but at historical/future dates where DE431↔DE441 baseline is larger. Verify sidereal bugs don't emerge at extreme dates.

- **Bodies**: Sun, Moon, Mars, MeanNode, TrueNode, Cupido
- **Dates**: 10 dates spread from -3000 to +2900 CE
- **Flags**: 4 sidereal combinations + 4 non-sidereal controls
- **Ayanamsha**: Lahiri (1)
- **Acceptance**: excess < 0.05" (relaxed due to baseline divergence, but still no LEB-specific error)

### Round 2.7 — Houses Three-Way

**Goal**: Verify house cusps computed via LEB match both Skyfield and SE.

- **House systems**: Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign, Porphyrius, Alcabitius, Morinus
- **Locations**: Rome (41.9°N), NYC (40.7°N), Tokyo (35.7°N), Sydney (33.9°S), Reykjavik (64.1°N), Equator (0°), Tromso (69.6°N near polar)
- **Dates**: 10 dates from 1800-2200 CE
- **Metrics**: For each: max |cusp_diff| across all 12 cusps + ASC + MC, three-way
- **Additional**: Test with `SEFLG_SIDEREAL` (Lahiri) as well — 2 modes per date/location
- **Acceptance**: LEB vs Skyfield < 0.001", LEB vs SE < 0.01" (modern era)

### Round 2.8 — Crossing Functions Three-Way

**Goal**: Verify swe_solcross_ut and swe_mooncross_ut produce consistent crossing times.

- **Sun crossings**: 0°, 90°, 180°, 270° (equinoxes/solstices) from 1900-2100 CE, three-way
- **Moon crossings**: 0°, 90°, 180°, 270° from 2000-2025 CE, three-way
- **Node crossings**: swe_solcross_node_ut if available
- **Metrics**: |dt| in seconds between each pair
- **Acceptance**: LEB vs Skyfield < 0.1s, LEB vs SE < 0.5s

---

## Phase 3: Coordinate System & Flag Edge Cases (6 rounds)

### Round 3.1 — 0°/360° Longitude Wrap-Around

**Goal**: Bodies crossing 0° Aries should not produce discontinuities or wrap-around errors.

- **Method**: For each body, find dates where longitude is near 0° (within 1°). Evaluate at 100 points around the crossing. Verify no jumps > 0.01° between consecutive evaluations spaced 1 hour apart.
- **Bodies**: Sun, Moon, Mars, Jupiter (all cross 0° regularly)
- **Additional**: Also test at 180° (opposition point) and near 90°/270°
- **Acceptance**: No jumps > 0.01° between hourly evaluations

### Round 3.2 — Latitude Extrema and Sign Changes

**Goal**: Verify latitude precision near zero crossings and at extrema.

- **Bodies**: Moon (up to ±5.3°), Mercury (up to ±7°), Pluto (up to ±17°), Pallas (high inclination)
- **Method**: For each body, find dates where latitude crosses zero and where it reaches extrema. Evaluate LEB vs Skyfield at 50 points around each event.
- **Acceptance**: |dLat| < 0.001" everywhere, including at zero crossings

### Round 3.3 — Ecliptic↔Equatorial Round-Trip Consistency

**Goal**: Computing ecliptic, converting to equatorial, and comparing with direct equatorial output should give identical results.

- **Method**: For 10 bodies × 20 dates:
  1. Compute with `SEFLG_SPEED` → ecliptic (lon, lat)
  2. Compute with `SEFLG_SPEED | SEFLG_EQUATORIAL` → equatorial (RA, dec)
  3. Manually convert (1) to equatorial using obliquity
  4. Compare (2) and (3) — must agree to machine precision
- **Acceptance**: < 1e-10° difference (floating-point noise only)

### Round 3.4 — J2000 Precession Round-Trip

**Goal**: Position in ecliptic of date should, when precessed to J2000, match the J2000 flag output.

- **Method**: For 10 bodies × 20 dates:
  1. Compute with `SEFLG_SPEED` (ecliptic of date)
  2. Compute with `SEFLG_SPEED | SEFLG_J2000` (J2000 ecliptic)
  3. Manually precess (1) to J2000 using `_precess_ecliptic`
  4. Compare (2) and (3)
- **Bodies**: Sun, Moon, Mars, Jupiter, MeanNode, Cupido
- **Additional test**: For bodies where SID+J2K suppresses J2000 (TrueNode, OscuApog, IntpApog, IntpPerg), verify that SID output equals SID+J2K output exactly.
- **Acceptance**: < 1e-8° for Pipeline A, exact match for J2K-suppressed bodies

### Round 3.5 — Sidereal↔Tropical Offset Consistency

**Goal**: The difference between sidereal and tropical output should exactly equal the ayanamsha, for all bodies.

- **Method**: For 10 bodies × 20 dates × 3 ayanamsha modes:
  1. Compute tropical longitude
  2. Compute sidereal longitude
  3. Compute ayanamsha via `swe_get_ayanamsa_ut`
  4. Verify: `(tropical - sidereal) mod 360 ≈ ayanamsha`
- **Tolerance**: < 1e-8° (should be exact to floating-point)
- **Additional**: Verify this also holds for latitude (should be unchanged) and speed (should differ by ayanamsha rate)

### Round 3.6 — Flag Orthogonality: All Meaningful Combinations

**Goal**: Verify all flag combinations produce valid, non-NaN, non-infinite results.

- **Bodies**: Sun, Moon, Mars, Jupiter, MeanNode, TrueNode, Cupido
- **Dates**: 5 dates (one per era: -5000, -1000, 2000, +5000, +10000 CE)
- **Flags**: All 32+ meaningful combinations of:
  - `SEFLG_SPEED`
  - `SEFLG_EQUATORIAL`
  - `SEFLG_J2000`
  - `SEFLG_SIDEREAL`
  - `SEFLG_TRUEPOS`
  - `SEFLG_NOABERR`
  - `SEFLG_NOGDEFL`
  - `SEFLG_HELCTR` (where valid)
  - `SEFLG_BARYCTR` (where valid)
- **Check**: No NaN, no infinity, lon in [0, 360), lat in [-90, 90], speed in reasonable range
- **Acceptance**: Zero NaN/infinity across all combinations

---

## Phase 4: Numerical Stability & Physical Consistency (5 rounds)

### Round 4.1 — Nutation, Obliquity, Delta-T at Extreme Dates

**Goal**: These auxiliary quantities feed into all coordinate transforms. Verify LEB, Skyfield, and SE agree.

- **Dates**: -13000, -10000, -5000, -2000, -1000, 0, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000, 17000 CE
- **Metrics**:
  - True obliquity: LEB vs Skyfield vs SE (via `swe_calc_ut` nutation)
  - Mean obliquity: same
  - Nutation in longitude (dpsi): same
  - Nutation in obliquity (deps): same
  - Delta-T: LEB vs Skyfield vs SE
- **Acceptance**: LEB = Skyfield exactly. LEB vs SE: obliquity/nutation within 0.01" in modern era, may diverge at extremes (different models).

### Round 4.2 — Planetary Distance Extrema: Perigee and Apogee

**Goal**: At orbital extrema, errors may be amplified. Verify precision at closest/farthest approach.

- **Bodies**: Moon (perigee ~356,500 km, apogee ~406,700 km), Mars (opposition), Venus (inferior conjunction)
- **Method**: Find 20 perigees and 20 apogees for Moon (2000-2025 CE). At each, evaluate distance LEB vs Skyfield.
- **Additional**: At Moon perigee, also evaluate longitude/latitude (highest speed = hardest for Chebyshev).
- **Acceptance**: Distance < 1e-8 AU at perigee, position < 0.001"

### Round 4.3 — Monotonicity of Mean Bodies

**Goal**: MeanNode longitude should decrease monotonically (retrograde). MeanApog should increase monotonically. Verify across 10,000 years.

- **Bodies**: MeanNode (ID 10), MeanApog (ID 12)
- **Method**: Evaluate at 10,000 dates (1 per year equivalent). Check:
  - MeanNode: lon(t+1) < lon(t) (allowing for 360° wrap)
  - MeanApog: lon(t+1) > lon(t) (allowing for 360° wrap)
  - Speed sign: MeanNode speed always negative, MeanApog speed always positive
- **Acceptance**: Zero violations

### Round 4.4 — Sun-Earth Heliocentric Consistency

**Goal**: The geocentric Sun position should be exactly 180° opposite to the heliocentric Earth position (with matching latitude/distance).

- **Method**: At 50 dates:
  1. Compute Sun geocentric: `swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)`
  2. Compute Earth heliocentric: `swe_calc_ut(jd, SE_EARTH, SEFLG_SPEED | SEFLG_HELCTR)`
  3. Verify: `Sun_lon ≈ (Earth_lon + 180) mod 360`, `Sun_lat ≈ -Earth_lat`, `Sun_dist ≈ Earth_dist`
- **Acceptance**: < 1e-8° for position, < 1e-12 AU for distance

### Round 4.5 — TrueNode Oscillation Pattern

**Goal**: TrueNode oscillates around MeanNode with a known ~18.6-year period and ~1.5° amplitude. Verify the pattern is physically plausible across the full range.

- **Method**: At 1000 dates across full range, compute both MeanNode and TrueNode. Compute `TrueNode - MeanNode`. Verify:
  - Amplitude stays within ±2.5° (known theoretical bound)
  - No sudden jumps > 0.5° between consecutive evaluations (spaced ~10 days)
- **Acceptance**: Zero violations of amplitude bound, zero jumps > 0.5°

---

## Phase 5: Cross-Tier Consistency (3 rounds)

### Round 5.1 — Extended vs Medium LEB: Overlap Zone

**Goal**: In the overlap zone (1550-2650 CE), the two LEB files were generated from different BSPs (DE441 vs DE440). Verify differences match the known DE440↔DE441 baseline.

- **Bodies**: All 31
- **Dates**: 100 dates uniformly in the overlap zone
- **Flags**: `SEFLG_SPEED`
- **Method**: Evaluate both `ephemeris_extended.leb` and `ephemeris_medium.leb` at same dates. Compute differences.
- **Acceptance**:
  - Moon: < 0.5" (known de440↔de441 divergence)
  - All others: < 0.01"
  - If any body shows much larger difference, it indicates a generation bug

### Round 5.2 — Extended vs Medium LEB: Sidereal in Overlap

**Goal**: Same as 5.1 but with sidereal flags. Since sidereal correction is applied post-LEB evaluation, the differences should be identical to the tropical case.

- **Bodies**: Sun, Moon, Mars, Jupiter, MeanNode, TrueNode, Cupido
- **Dates**: 30 dates in overlap zone
- **Flags**: `SEFLG_SPEED | SEFLG_SIDEREAL` (Lahiri) + `SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_EQUATORIAL`
- **Acceptance**: Cross-tier difference in sidereal ≈ cross-tier difference in tropical (< 0.001" delta)

### Round 5.3 — Extended vs Medium LEB: Speed Comparison

**Goal**: Speed precision should be comparable between tiers in the overlap zone.

- **Bodies**: Moon, Mercury, Mars, TrueNode, OscuApog
- **Dates**: 100 dates in overlap
- **Metrics**: |speed_ext - speed_med| for all 3 speed components
- **Acceptance**: Differences consistent with de440↔de441 baseline (< 0.01 d/d for most bodies)

---

## Phase 6: Exhaustive Ayanamsha & Sidereal Coverage (3 rounds)

### Round 6.1 — All 43 Ayanamsha Modes: Full Body Sweep

**Goal**: Verify every ayanamsha mode produces correct sidereal output on extended tier.

- **Bodies**: Sun, Moon, Mars, Jupiter, Saturn, MeanNode, TrueNode, OscuApog, Cupido
- **Dates**: 5 dates (1900, 1950, 2000, 2050, 2100 CE)
- **Ayanamsha**: All 43 modes (0-42)
- **Flags**: `SEFLG_SPEED | SEFLG_SIDEREAL`
- **Method**: LEB vs Skyfield comparison
- **Acceptance**: 0.000000" for all (ayanamsha is subtracted identically)

### Round 6.2 — Ayanamsha at Extreme Dates

**Goal**: Some ayanamsha modes (True Citra, True Revati, etc.) depend on star positions that may behave differently at extreme dates. Verify stability.

- **Bodies**: Sun, Moon
- **Dates**: -5000, -1000, 0, 1000, 2000, 5000, 10000 CE
- **Ayanamsha**: All 43 modes
- **Check**: No NaN, no infinity, values in plausible range (0-360°). LEB vs Skyfield agreement.
- **Acceptance**: Zero NaN/infinity, LEB = Skyfield within 0.001"

### Round 6.3 — SID+J2000 Suppression Verification (Pipeline B)

**Goal**: Explicitly verify that the J2000-suppression fix for TrueNode/OscuApog/IntpApog/IntpPerg works on extended tier.

- **Bodies**: MeanNode (should precess to J2K), TrueNode (should NOT), MeanApog (should precess), OscuApog (should NOT), IntpApog (should NOT), IntpPerg (should NOT)
- **Dates**: 10 dates from 1800-2200 CE
- **Flags**: `SEFLG_SIDEREAL` and `SEFLG_SIDEREAL | SEFLG_J2000`
- **Verification**:
  - For suppressed bodies: `LEB(SID) == LEB(SID+J2K)` exactly (diff = 0.0000")
  - For non-suppressed bodies: `LEB(SID) != LEB(SID+J2K)` (large precession difference)
  - Same check on Skyfield path
  - Same check on pyswisseph
- **Acceptance**: Zero discrepancy with pyswisseph behavior

---

## Phase 7: Pathological & Degenerate Cases (4 rounds)

### Round 7.1 — High-Latitude Houses: Polar Edge Cases

**Goal**: House systems (especially Placidus, Koch) have singularities near the Arctic/Antarctic circles. Verify no crashes or NaN.

- **Latitudes**: 60°, 63°, 66°, 66.5° (Arctic circle), 67°, 70°, 80°, 89° (near pole)
- **Also**: -60°, -66.5°, -70°, -80°, -89° (southern)
- **House systems**: Placidus, Koch, Regiomontanus, Campanus
- **Dates**: Summer solstice, winter solstice, equinox at 2000 CE
- **Check**: No crash, no NaN, cusps in [0, 360), cusps ordered monotonically
- **Acceptance**: Zero crashes/NaN. Values may be degenerate (e.g., intercepted signs) but must be finite.

### Round 7.2 — Bodies at Maximum Geocentric Speed

**Goal**: When Moon is at maximum speed (~15.4°/day), Chebyshev approximation is hardest. Verify precision at speed extrema.

- **Method**: Scan Moon speed over 2000-2025 CE. Find 20 dates with highest |speed|. Evaluate LEB vs Skyfield at those dates.
- **Additional**: Same for Mercury (max ~2.2°/day during direct motion far from station)
- **Acceptance**: Position < 0.001", speed < 0.002 d/d even at maximum speed

### Round 7.3 — Retrograde Stations: Speed Zero-Crossing Precision

**Goal**: At retrograde stations, speed crosses zero. This is where the speed polynomial is most stressed.

- **Bodies**: Mercury, Venus, Mars, Jupiter, Saturn
- **Method**: Find 50 stations (2000-2025 CE) via sign-change detection in speed. At each station ± [0.001, 0.01, 0.1, 1] days, compare LEB vs Skyfield.
- **Additional three-way**: At each station, also compare with pyswisseph.
- **Acceptance**: Position < 0.001", speed < 0.001 d/d (away from exact zero), sign agreement when |speed| > 0.001 d/d

### Round 7.4 — Very Close Conjunctions and Oppositions

**Goal**: When two bodies have nearly the same longitude, angular differences should be preserved by LEB.

- **Method**: Find 20 Sun-Moon conjunctions (new moons) and 20 Sun-Moon oppositions (full moons) in 2000-2025. At each, compute:
  - `Sun_lon - Moon_lon` via LEB
  - `Sun_lon - Moon_lon` via Skyfield
  - The two should agree to < 0.001"
- **Additional**: Find a few Mercury-Venus conjunctions and verify
- **Acceptance**: Angular difference error < 0.001"

---

## Phase 8: Full Statistical Report (3 rounds)

### Round 8.1 — Massive Sweep: 500 Dates × 31 Bodies × 5 Flag Combos

**Goal**: The definitive statistical profile of the extended LEB.

- **Bodies**: All 31
- **Dates**: 500 uniformly distributed across full range
- **Flags**: default, EQUATORIAL, J2000, SIDEREAL, SIDEREAL+EQUATORIAL
- **Metrics**: Per body per flag combo: N, mean, P50, P95, P99, max of |dLon|
- **Output**: Grand table suitable for documentation
- **Additional**: Group statistics by era (pre-CE, 0-1500 CE, 1500-2500 CE, post-2500 CE) to show if precision varies by epoch

### Round 8.2 — Precision Comparison: Extended vs Medium on Overlap Zone

**Goal**: Side-by-side comparison of precision statistics on the same dates using both tiers.

- **Method**: Pick 200 dates in the overlap zone. For each body and date, compute:
  - `|ext_LEB - Skyfield(de441)|`
  - `|med_LEB - Skyfield(de440)|`
- **Output**: Table showing that the two tiers have comparable precision where they overlap
- **Acceptance**: Extended tier precision should be within 2× of medium tier precision for all bodies

### Round 8.3 — Final Acceptance Summary

**Goal**: Single table with pass/fail for every body, every phase.

- **Output**: Matrix of 31 bodies × all rounds, with pass/warn/fail per cell
- **Output**: Overall verdict: PASS if zero FAIL cells (WARN is acceptable)
- **Output**: Known limitations list with explanations

---

## Tolerances Reference

| Body Category | Position vs Skyfield | Position vs SE (modern) | Speed vs Skyfield |
|---|---|---|---|
| Moon | < 0.5" | < 0.05" | < 0.002 d/d |
| Sun, Inner Planets | < 0.001" | < 0.005" | < 0.001 d/d |
| Outer Planets (Jup–Nep) | < 0.001" | < 0.03" | < 0.001 d/d |
| Pluto | < 0.025" | < 0.1" | < 0.001 d/d |
| MeanNode, MeanApog | < 0.001" | model diff ~0.1" | < 0.001 d/d |
| TrueNode | < 5" (model diff) | model diff ~2-5" | < 0.01 d/d |
| OscuApog | < 100" (model diff) | model diff ~57-90" | < 0.05 d/d |
| IntpApog/IntpPerg | < 5000" (model diff) | model diff ~200-4600" | < 1.5 d/d |
| Uranians (40-47) | < 35" (model diff) | model diff ~5-33" | < 0.001 d/d |
| Transpluto (48) | < 10" (model diff) | model diff ~1-7" | < 0.001 d/d |
| Chiron | < 0.001" | < 0.005" | < 0.001 d/d |
| Ceres–Vesta | < 0.001" | < 0.1" | < 0.001 d/d |

---

## Known Limitations (Not Bugs)

1. **Chiron SPK range**: ~660-4600 CE only. ERR outside this range is expected.
2. **IntpApog/IntpPerg SE range**: ~-3000 to +2900 CE only in pyswisseph.
3. **MeanNode/MeanApog Meeus warning**: Polynomial extrapolation beyond ±20 centuries from J2000 has degraded precision. Warning is emitted by libephemeris.
4. **Neptune/Uranus Skyfield speed artifact**: Isolated speed spikes due to finite-difference calculation in Skyfield reference path. Position is unaffected.
5. **Transpluto velocity -180 d/d**: Artifact of heliocentric coordinate system for a nearly stationary body at 55 AU. Not a real motion.
6. **HELCTR error amplification**: Heliocentric positions subtract Earth from body, amplifying Chebyshev error by ~8-11× for inner planets. This is inherent to the derivation, not a bug.
7. **DE431↔DE441 baseline**: pyswisseph uses DE431-based SE files, LEB uses DE441. Moon diverges by ~90" at 0 CE and ~280" at 2650 CE. This is not a libephemeris error.
8. **Pholus**: Not included in LEB (no SPK data). Falls back to Skyfield at runtime.

---

## Execution Notes

- Configure Skyfield to use DE441 via `set_ephe_path('/Volumes/data/libephemeris')` + `set_ephemeris_file('de441.bsp')` for all extended tier validation.
- Configure pyswisseph to use full SE files via `swe.set_ephe_path('swisseph/ephe')`.
- Rounds within each phase can be parallelized in pairs.
- Estimated total execution time: 45-60 minutes.
- Each round should produce a self-contained summary that can be appended to a final report.
