# LEB Precision Improvement Plan

## 1. Context

### 1.1 What is LEB

LEB (LibEphemeris Binary) is a precomputed binary ephemeris format that replaces
live Skyfield/JPL DE440 calculations with Chebyshev polynomial lookups. It
provides a ~14x speedup for astrological calculations.

The file stores 30 celestial bodies as Chebyshev polynomial segments. Each body
has two key parameters:

- **interval_days**: how many days each polynomial segment covers
- **degree**: the polynomial degree (degree N = N+1 coefficients)

Shorter intervals and higher degrees = better precision but larger file.

### 1.2 The three pipelines

When a user calls `swe_calc_ut()` in LEB mode, the code in `fast_calc.py`
dispatches to one of three pipelines:

| Pipeline | Bodies | Storage | Count |
|----------|--------|---------|-------|
| A (ICRS) | Sun, Moon, Mercury-Pluto, Earth, Chiron, Ceres-Vesta | Barycentric Cartesian (x,y,z) in AU | 16 bodies |
| B (Ecliptic) | Mean/True Node, Mean/Oscu Apogee, Interp Apogee/Perigee | Ecliptic of date (lon, lat, dist) | 6 bodies |
| C (Heliocentric) | Uranians (Cupido-Poseidon), Transpluto | Heliocentric ecliptic (lon, lat, dist) | 9 bodies |

Pipeline A is the most complex: it must transform ICRS barycentric positions
through geocentric subtraction, light-time correction, aberration,
precession-nutation, and obliquity rotation to produce ecliptic coordinates.

### 1.3 How velocity works today

This is the core of the problem. The LEB reader (`leb_reader.py:337`) already
computes the **analytical Chebyshev derivative** at every call:

```python
# leb_reader.py, eval_body(), line 337
val, deriv = _clenshaw_with_derivative(comp_coeffs, tau)
pos.append(val)
vel.append(deriv * scale)  # scale = 2.0 / interval_days
```

So `eval_body()` returns `(position, velocity)` where velocity is the exact
mathematical derivative of the fitted polynomial — free and precise.

But Pipeline A in `fast_calc.py` **throws this velocity away** and recomputes
it using numerical central difference:

```python
# fast_calc.py, _fast_calc_core(), lines 781-803
lon, lat, dist = _pipeline_icrs(reader, jd_tt, ipl, iflag)  # position

# Velocity via central difference
lon_prev, lat_prev, dist_prev = _pipeline_icrs(reader, jd_tt - dt, ipl, ...)
lon_next, lat_next, dist_next = _pipeline_icrs(reader, jd_tt + dt, ipl, ...)
dlon = (lon_next - lon_prev) / (2.0 * dt)
```

This means:
- The good velocity (analytical derivative) is computed then **discarded**
- The bad velocity (numerical difference of two approximate positions) is used
- Each velocity costs **2 extra full pipeline runs** (each involving multiple
  `eval_body()` calls for light-time iterations)
- The numerical difference **amplifies** the Chebyshev position error into a
  larger velocity error

Pipelines B and C do NOT have this problem — they use the Chebyshev derivative
directly (`leb_reader.py:514`).

### 1.4 Current BODY_PARAMS

Source: `libephemeris/leb_format.py:149-185`

| ID | Body | Interval | Degree | Coord Type |
|----|------|----------|--------|------------|
| 0 | Sun | 32 | 13 | ICRS_BARY |
| 1 | Moon | 8 | 13 | ICRS_BARY |
| 2 | Mercury | 16 | 15 | ICRS_BARY |
| 3 | Venus | 32 | 13 | ICRS_BARY |
| 4 | Mars | 32 | 13 | ICRS_BARY |
| 5 | Jupiter | 64 | 11 | ICRS_BARY |
| 6 | Saturn | 64 | 11 | ICRS_BARY |
| 7 | Uranus | 128 | 9 | ICRS_BARY |
| 8 | Neptune | 128 | 9 | ICRS_BARY |
| 9 | Pluto | 128 | 9 | ICRS_BARY |
| 14 | Earth | 8 | 13 | ICRS_BARY |
| 10-13,21-22 | Lunar/Lilith | 8 | 13 | ECLIPTIC |
| 15,17-20 | Asteroids | 32 | 13 | ICRS_BARY |
| 40-48 | Uranians+Transpluto | 64 | 11 | HELIO_ECL |

Nutation: interval=32 days, degree=16, 2 components (dpsi, deps).

---

## 2. Problems

### 2.1 Chebyshev fit too coarse for outer planets

Uranus, Neptune, and Pluto use `(128, 9)`: 128-day intervals with degree-9
polynomials (10 coefficients). The position error reaches ~4.8 arcseconds for
Uranus. This is the worst precision of any body in the file.

The rationale was "slow-moving planets need less resolution", which is true for
the orbital motion itself, but misses the fact that these bodies are stored in
**ICRS barycentric** coordinates — their apparent position (as seen from Earth)
can vary faster due to Earth's own motion and parallax effects.

### 2.2 Moon and Earth intervals could be shorter

Moon and Earth use `(8, 13)`. The Moon's position error is ~0.05 arcseconds,
which is good but not optimal. More importantly, **Earth is the observer** for
all geocentric calculations — any error in Earth's position propagates to the
geocentric vector of every other body. Halving the interval for both would
improve precision for all 30 bodies, not just Moon and Earth.

### 2.3 Velocity error from numerical differentiation (Pipeline A)

The central difference amplifies Chebyshev position errors. For the Moon
(fast-moving, ~13 deg/day), this produces velocity errors up to ~14"/day
(~0.004 deg/day). The error comes from differentiating the *approximation
error* itself:

```
v_error ≈ d(position_error)/dt
```

Since the position error oscillates like the (degree+1)-th Chebyshev polynomial,
its derivative oscillates with much larger amplitude. The faster the body moves,
the worse this gets.

This affects all 16 ICRS bodies (Pipeline A), but Moon is worst because:
- It moves fastest (~13 deg/day)
- Its Chebyshev segments are relatively short (8 days) → steeper error
  oscillation

### 2.4 Nutation precision affects latitude

The nutation Chebyshev (32-day intervals, degree 16) introduces small obliquity
errors. These primarily affect the ecliptic latitude of all bodies. The effect
is amplified for nearby planets (Venus at ~0.26 AU, Mars at ~0.37 AU) due to
the 1/distance factor in angular error conversion.

### 2.5 Asteroid SPK coverage gaps (base tier)

The base tier LEB file covers 1850-2150 (de440s.bsp), but the asteroid SPK
files have **shorter coverage**. During LEB generation, the generator falls
back to Keplerian orbital elements for dates outside SPK coverage. This
produces catastrophically wrong data — errors of 7,000-14,000 arcseconds in
position and 0.03-0.05 AU in distance.

The LEB format already stores **per-body date ranges** via `BodyEntry.jd_start`
and `BodyEntry.jd_end` (`leb_format.py:108-109`). At runtime, `eval_body()`
raises `ValueError` for out-of-range dates and `fast_calc.py` falls back to
Skyfield. However, the **base tier compare tests** sample dates across the
full 1860-2140 range without checking per-body coverage, causing 17 test
failures.

This is documented in `docs/leb/guide.md` section 10.4.

#### 2.5.1 Observed failures from base tier test suite

Source: `poe test:leb:compare:base -n 8` (log: `log_leb_base.log`)

Result: **387 passed, 17 failed** out of 404 tests.

**All 17 failures are asteroids. No planet, lunar body, Uranian, hypothetical,
sidereal, or flag-combination test fails.**

**Position failures** (5 tests in `test_base_asteroids.py`):

| Asteroid | Max error | Worst JD | ~Year | Tolerance |
|----------|-----------|----------|-------|-----------|
| Chiron (15) | 7,646" | 2502306.7 | ~2139 | 5.0" |
| Ceres (17) | 11,761" | 2402149.7 | ~1865 | 5.0" |
| Pallas (18) | 10,222" | 2500597.5 | ~2134 | 5.0" |
| Juno (19) | 14,034" | 2401124.2 | ~1862 | 5.0" |
| Vesta (20) | 12,008" | 2401807.8 | ~1864 | 5.0" |

**Distance failures** (7 tests across 2 files):

| Asteroid | Max error | Worst JD | ~Year | File |
|----------|-----------|----------|-------|------|
| Pallas (18) | 0.033 AU | 2401124.2 | ~1862 | test_base_asteroids |
| Juno (19) | 0.044 AU | 2402149.7 | ~1865 | test_base_asteroids |
| Vesta (20) | 0.036 AU | 2409328.2 | ~1884 | test_base_asteroids |
| Chiron (15) | 0.038 AU | 2401466.0 | ~1863 | test_base_asteroids |
| Ceres (17) | 0.044 AU | 2402491.5 | ~1866 | test_base_asteroids |
| Chiron (15) | 0.039 AU | 2401812.4 | ~1864 | test_base_distances |
| Ceres (17) | 0.046 AU | 2402498.4 | ~1866 | test_base_distances |

**Distance velocity failures** (5 tests in `test_base_velocities.py`):

| Asteroid | Max error | Worst JD | ~Year |
|----------|-----------|----------|-------|
| Ceres (17) | 6.44e-4 AU/day | 2499218.7 | ~2130 |
| Chiron (15) | 5.93e-4 AU/day | 2402498.4 | ~1866 |
| Vesta (20) | 7.51e-4 AU/day | 2401812.4 | ~1864 |
| Juno (19) | 7.03e-4 AU/day | 2479325.9 | ~2076 |
| Pallas (18) | 5.61e-4 AU/day | 2408672.0 | ~1882 |

Key observations:
- Worst dates cluster around ~1862-1884 and ~2130-2139 — the boundaries of SPK
  coverage within the base tier range.
- Juno distance-velocity failure at JD 2479325.9 (~2076) indicates a **gap or
  early termination** in Juno's SPK file, not just a boundary issue.
- Longitude velocity and latitude velocity for asteroids **all pass** — the
  Keplerian fallback gives a reasonable orbital direction but wrong position
  along the orbit, so angular rates are approximately correct while absolute
  position and distance are catastrophically wrong.

### 2.6 Medium tier crossing test failures (known, separate issue)

The medium tier compare suite has 12 remaining failures, **none of which are
LEB precision problems**:
- 9 crossing tests (Jupiter/Saturn): `swe_cross_ut()` throws
  `RuntimeError: Planet crossing search diverged` — a bug in the crossing
  solver, not in LEB.
- 2 crossing tests (Saturn helio): same divergence issue.
- 1 solar eclipse test (Sydney): no eclipse found at location.

These are tracked separately and do not affect the precision improvement plan.

### 2.7 Summary of current errors vs targets

| Metric | Current worst | Acceptable target |
|--------|---------------|-------------------|
| Position (outer planets) | ~4.8" (Uranus) | <0.5" |
| Position (Moon) | ~0.05" | <0.01" |
| Position (Venus lat) | ~0.52" | <0.3" |
| Velocity (Moon) | ~14"/day | <0.1"/day |
| Velocity (all ICRS) | amplified | analytical |
| SEFLG_SPEED cost | 3 pipeline calls | 1 pipeline call |
| Asteroids (base tier) | 7,000-14,000" | N/A (test fix) |

---

## 3. Implementation Plan

### Phase 0: Fix asteroid test date ranges (base tier)

**Goal**: eliminate the 17 base tier test failures caused by asteroid SPK
coverage gaps. This is a **test-only fix** — no production code changes.

**Root cause**: the base tier tests sample dates across 1860-2140, but asteroid
SPK files cover a narrower range. The LEB generator stored Keplerian-fallback
data for uncovered dates, producing catastrophic errors (7,000-14,000"). The
test `try/except` blocks catch `KeyError`/`ValueError` from `eval_body()` when
dates are out of the LEB body range, but the Keplerian data is stored *inside*
the body's range — it's bad data that doesn't raise exceptions.

**Solution**: read per-body date ranges from the LEB file and filter test dates
to only include dates within each asteroid's valid SPK coverage. The
`LEBReader` exposes `_bodies[body_id].jd_start` / `jd_end` via the `BodyEntry`
dataclass (`leb_format.py:102-113`). At runtime, `eval_body()` already raises
`ValueError` for out-of-range dates (`leb_reader.py:301-303`).

**Files to modify**:

| File | Change |
|------|--------|
| `tests/test_leb/compare/base/conftest.py` | Add `get_body_date_range()` helper that opens the LEB file and reads `BodyEntry.jd_start`/`jd_end` for a given body_id |
| `tests/test_leb/compare/base/test_base_asteroids.py` | Filter `base_dates_300` to body's valid range before iterating |
| `tests/test_leb/compare/base/test_base_distances.py` | Filter `base_dates_150` for Chiron (15) and Ceres (17) in `DISTANCE_BODIES` |
| `tests/test_leb/compare/base/test_base_velocities.py` | Filter `base_dates_150` for `ASTEROID_BODIES` in `TestBaseDistanceVelocity` |

**Implementation**:

```python
# In tests/test_leb/compare/base/conftest.py:

from libephemeris.leb_reader import LEBReader

def get_body_date_range(leb_path: str, body_id: int) -> tuple[float, float]:
    """Read the valid date range for a body from the LEB file header.

    Returns:
        (jd_start, jd_end) for the given body.
    """
    reader = LEBReader(leb_path)
    body = reader._bodies[body_id]
    return (body.jd_start, body.jd_end)


def filter_dates_for_body(
    dates: list[float], leb_path: str, body_id: int, margin: float = 30.0
) -> list[float]:
    """Filter test dates to only include those within a body's LEB coverage.

    Args:
        dates: list of JD values to filter.
        leb_path: path to the LEB file.
        body_id: SE body ID.
        margin: days of margin inside the boundaries.

    Returns:
        Filtered list of JDs within [jd_start + margin, jd_end - margin].
    """
    jd_start, jd_end = get_body_date_range(leb_path, body_id)
    return [jd for jd in dates if jd_start + margin <= jd <= jd_end - margin]
```

Then in each test method, before the loop:

```python
# Example in test_base_asteroids.py:
def test_position(self, compare, base_dates_300, body_id, body_name):
    filtered_dates = filter_dates_for_body(
        base_dates_300, compare.leb_path, body_id
    )
    if not filtered_dates:
        pytest.skip(f"{body_name}: no test dates within SPK coverage")

    for jd in filtered_dates:
        # ... existing test logic
```

**Expected result**: all 17 failures become passes (dates outside SPK coverage
are skipped). The tests still validate precision for dates within coverage.

**Why not `xfail`**: marking tests as `xfail` would hide real regressions.
Filtering dates is better because it tests the *valid* data and ignores the
*structurally invalid* data.

---

### Phase 1: Chebyshev parameter tuning

**Goal**: reduce fitting error for the worst bodies.

**Files to modify**:
- `libephemeris/leb_format.py` (BODY_PARAMS table, lines 149-185)
- `scripts/generate_leb.py` (NUTATION_INTERVAL constant, line 182)

**Changes**:

| Body | From | To | Rationale |
|------|------|----|-----------|
| Uranus (7) | (128, 9) | (64, 13) | Position from ~4.8" to <0.5" |
| Neptune (8) | (128, 9) | (64, 13) | Same problem as Uranus |
| Pluto (9) | (128, 9) | (64, 13) | Same problem as Uranus |
| Moon (1) | (8, 13) | (4, 13) | Position from ~0.05" to ~0.003" |
| Earth (14) | (8, 13) | (4, 13) | Earth is the observer: improves ALL geocentric positions |
| Nutation | (32, 16) | (16, 16) | Halves obliquity error → helps Venus/Mars latitude |

**What stays the same**:
- Sun (0): `(32, 13)` — already <0.07", good enough
- Mercury (2): `(16, 15)` — already <0.01", excellent
- Venus (3), Mars (4): `(32, 13)` — their fit is perfect (1e-12 AU); the
  remaining error is pipeline-level, not fitting-level
- Jupiter (5), Saturn (6): `(64, 11)` — already <0.05"
- All ecliptic bodies (10-13, 21-22): `(8, 13)` — already <0.5"
- All asteroids (15, 17-20): `(32, 13)` — precision depends on SPK quality
- All Uranians (40-48): `(64, 11)` — analytical, already <0.001"

**File size impact** (medium tier, 1550-2650, ~401,768 days):

| Body | Before | After | Delta |
|------|--------|-------|-------|
| Uranus | 753 KB | 2,109 KB | +1,356 KB |
| Neptune | 753 KB | 2,109 KB | +1,356 KB |
| Pluto | 753 KB | 2,109 KB | +1,356 KB |
| Moon | 3,289 KB | 6,578 KB | +3,289 KB |
| Earth | 3,289 KB | 6,578 KB | +3,289 KB |
| Nutation | ~620 KB | ~1,240 KB | +620 KB |
| **Total delta** | | | **+11,266 KB (~11 MB)** |
| **New file size** | ~88 MB | **~99 MB** | **+12.5%** |

**Verification**: after regeneration, run `generate_leb.py --verify` which
compares every body against its reference source and reports max error in
arcseconds. All bodies should show PASS (< 1").

**Implementation steps**:

1. Edit `libephemeris/leb_format.py`:
   - Line 152: `1: (8, 13, COORD_ICRS_BARY, 3),` → `1: (4, 13, COORD_ICRS_BARY, 3),`
   - Line 158: `7: (128, 9, COORD_ICRS_BARY, 3),` → `7: (64, 13, COORD_ICRS_BARY, 3),`
   - Line 159: `8: (128, 9, COORD_ICRS_BARY, 3),` → `8: (64, 13, COORD_ICRS_BARY, 3),`
   - Line 160: `9: (128, 9, COORD_ICRS_BARY, 3),` → `9: (64, 13, COORD_ICRS_BARY, 3),`
   - Line 161: `14: (8, 13, COORD_ICRS_BARY, 3),` → `14: (4, 13, COORD_ICRS_BARY, 3),`

2. Edit `scripts/generate_leb.py`:
   - Line 182: `NUTATION_INTERVAL = 32.0` → `NUTATION_INTERVAL = 16.0`

3. Regenerate the LEB file:
   ```bash
   poe leb:generate:medium:groups
   ```

4. Verify:
   ```bash
   python scripts/generate_leb.py --tier medium --verify --verify-samples 1000
   ```

5. Update documentation:
   - `docs/leb/guide.md` section 9.1 (BODY_PARAMS table)
   - `docs/leb/guide.md` section 9.2 (parameter rationale)
   - `docs/leb/guide.md` section 10 (precision tables)
   - `docs/leb/guide.md` Appendix A (file size estimation examples)

---

### Phase 2: Analytical velocity for Pipeline A

**Goal**: replace the numerical central difference with the analytical Chebyshev
derivative, transformed through the same rotation matrices already computed for
position.

**File to modify**: `libephemeris/fast_calc.py`

#### 2.1 What we remove

The central difference block in `_fast_calc_core()` (lines 781-803):

```python
# REMOVE THIS:
if iflag & SEFLG_SPEED:
    dt = 1.0 / 86400.0
    flags_no_speed = iflag & ~SEFLG_SPEED
    lon_prev, lat_prev, dist_prev = _pipeline_icrs(reader, jd_tt - dt, ipl, flags_no_speed)
    lon_next, lat_next, dist_next = _pipeline_icrs(reader, jd_tt + dt, ipl, flags_no_speed)
    dlon = (lon_next - lon_prev) / (2.0 * dt)
    # ... wrap handling, dlat, ddist
```

#### 2.2 What we keep

Everything about position computation in `_pipeline_icrs()` stays the same:
- `eval_body()` calls (already return velocity — we just stop ignoring it)
- Geocentric vector computation
- Light-time iteration
- Aberration
- Precession-nutation matrix
- Obliquity rotation

#### 2.3 What we add

Modify `_pipeline_icrs()` to optionally return velocity alongside position.
The velocity vector goes through the same transformations as position:

```
1. geo_vel = target_vel - observer_vel                    (geocentric velocity)
2. Apply light-time velocity correction                   (retarded velocity)
3. Apply aberration velocity correction                   (aberrated velocity)
4. Transform through precession-nutation matrix           (equatorial of date)
5. Rotate equatorial → ecliptic using true obliquity      (ecliptic of date)
6. Convert Cartesian velocity to spherical velocity       (dlon/dt, dlat/dt, ddist/dt)
```

Step 6 (Cartesian to spherical velocity) uses the standard formula:

```python
def _cartesian_velocity_to_spherical(
    x, y, z, vx, vy, vz
) -> Tuple[float, float, float]:
    """Convert Cartesian velocity to spherical velocity.

    Given position (x,y,z) and velocity (vx,vy,vz), compute:
        dlon/dt  = (x*vy - y*vx) / (x² + y²)          [rad/unit_time]
        dlat/dt  = (vz*(x²+y²) - z*(x*vx+y*vy)) / (r²*sqrt(x²+y²))  [rad/unit_time]
        ddist/dt = (x*vx + y*vy + z*vz) / r            [AU/day]
    """
    r_xy_sq = x*x + y*y
    r_sq = r_xy_sq + z*z
    r = math.sqrt(r_sq)
    r_xy = math.sqrt(r_xy_sq)

    if r == 0.0 or r_xy == 0.0:
        return (0.0, 0.0, 0.0)

    dlon_rad = (x * vy - y * vx) / r_xy_sq        # rad/day
    dlat_rad = (vz * r_xy_sq - z * (x*vx + y*vy)) / (r_sq * r_xy)  # rad/day
    ddist = (x*vx + y*vy + z*vz) / r               # AU/day

    dlon_deg = math.degrees(dlon_rad)               # deg/day
    dlat_deg = math.degrees(dlat_rad)               # deg/day

    return (dlon_deg, dlat_deg, ddist)
```

#### 2.4 Detailed implementation

**Step 1**: Add `_cartesian_velocity_to_spherical()` helper function.

**Step 2**: Modify `_pipeline_icrs()` signature:

```python
# BEFORE:
def _pipeline_icrs(reader, jd_tt, ipl, iflag) -> Tuple[float, float, float]:
    # returns (lon, lat, dist)

# AFTER:
def _pipeline_icrs(reader, jd_tt, ipl, iflag, want_velocity=False)
    -> Union[Tuple[float, float, float], Tuple[float, float, float, float, float, float]]:
    # returns (lon, lat, dist) or (lon, lat, dist, dlon, dlat, ddist)
```

**Step 3**: Inside `_pipeline_icrs()`, when `want_velocity=True`:

```python
# 1. We already have target_vel and observer_vel from eval_body() calls
target_pos, target_vel = reader.eval_body(ipl, jd_tt)
earth_pos, earth_vel = reader.eval_body(SE_EARTH, jd_tt)

# 2. Geocentric velocity
geo_vel = _vec3_sub(target_vel, earth_vel)

# 3. Light-time velocity correction
#    At the retarded time, we also get the retarded velocity:
retarded_pos, retarded_vel = reader.eval_body(ipl, jd_tt - lt)
geo_vel = _vec3_sub(retarded_vel, earth_vel)
#    (We use the velocity at retarded time, same as position)

# 4. Aberration velocity: for first-order aberration, the velocity correction
#    is negligible (aberration depends on Earth velocity which changes slowly).
#    We skip the velocity component of aberration — the position aberration
#    is already applied to geo, and the velocity of aberration is ~1e-8 deg/day.

# 5. Apply same rotation matrix to velocity vector
geo_vel_eq = _mat3_vec3(pn_mat, geo_vel)        # same pn_mat as position
ecl_vel = _rotate_equatorial_to_ecliptic(        # same eps_true_rad
    geo_vel_eq[0], geo_vel_eq[1], geo_vel_eq[2], eps_true_rad
)

# 6. Convert to spherical velocity
dlon, dlat, ddist = _cartesian_velocity_to_spherical(
    ecl[0], ecl[1], ecl[2],     # position (already computed)
    ecl_vel[0], ecl_vel[1], ecl_vel[2]  # velocity
)
```

**Step 4**: Modify `_fast_calc_core()` to use the new approach:

```python
# BEFORE (lines 777-803):
lon, lat, dist = _pipeline_icrs(reader, jd_tt, ipl, iflag)
if iflag & SEFLG_SPEED:
    # ... central difference (2 extra pipeline calls)

# AFTER:
if iflag & SEFLG_SPEED:
    lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
        reader, jd_tt, ipl, iflag, want_velocity=True
    )
else:
    lon, lat, dist = _pipeline_icrs(reader, jd_tt, ipl, iflag)
    dlon, dlat, ddist = 0.0, 0.0, 0.0
```

**Step 5**: Handle all coordinate frame branches in `_pipeline_icrs()`:

The velocity transform must follow the same branch as position:

| Frame | Position transform | Velocity transform |
|-------|-------------------|-------------------|
| True ecliptic of date (default) | ICRS → precess-nutate → equatorial → ecliptic rotation | Same matrices applied to velocity vector |
| Equatorial of date | ICRS → precess-nutate | Same matrix applied to velocity vector |
| J2000 equatorial | ICRS (identity) | Velocity is already in ICRS (identity) |
| J2000 ecliptic | ICRS → obliquity rotation | Same obliquity rotation on velocity vector |
| Heliocentric | Sun as observer | Sun velocity as observer velocity |
| Barycentric | Origin as observer | Zero observer velocity |

All branches use linear transforms (rotation matrices), which apply identically
to position and velocity vectors. No special handling needed — just apply the
same matrix to both.

#### 2.5 Edge cases

**Light-time velocity**: the velocity at the retarded time `t - lt` differs
from the velocity at `t`. We use `retarded_vel` (the velocity at the time the
light was emitted), consistent with using `retarded_pos` for position. The
light-time itself depends on distance and changes over time, but this
second-order effect (d(lt)/dt ≈ v_radial/c ≈ 1e-4) is negligible for
astrological precision.

**Aberration velocity**: the aberration correction depends on Earth's velocity,
which changes at ~0.017 deg/day² (Earth's orbital acceleration). Over the
1-second central difference interval this was ~1e-12 deg — negligible. With
the analytical approach we simply skip the velocity component of aberration.
The position aberration is already applied and its effect is captured in the
spherical velocity conversion.

**SEFLG_TRUEPOS, SEFLG_NOABERR**: these flags skip light-time or aberration
for position. The velocity path should respect the same flags (skip the
corresponding velocity corrections).

#### 2.6 Expected results

| Metric | Before | After |
|--------|--------|-------|
| Moon velocity error | ~14"/day | <0.1"/day (est. from Chebyshev derivative precision) |
| All ICRS body velocity | numerically amplified | analytically precise |
| Pipeline A + SPEED cost | 3 pipeline runs | 1 pipeline run |
| Pipeline A + SPEED time | ~8 μs | ~3-4 μs |

#### 2.7 Testing strategy

The existing test infrastructure is sufficient to validate this change:

1. **`tests/test_leb/test_leb_precision.py`**: compares `fast_calc_ut()` vs
   `swe_calc_ut()` for speed components. Current tolerance is 0.01 deg/day.
   After Phase 2, Moon speed should be well within this.

2. **`tests/test_leb/compare/test_compare_leb_velocities.py`**: compares LEB
   velocity vs Skyfield velocity for all 30 bodies. Current tolerance is
   0.05 deg/day. Should tighten after Phase 2.

3. **New targeted test**: add a Moon velocity precision test that samples
   ~1000 dates and asserts max speed error < 0.001 deg/day (~3.6"/day).

4. **Regression check**: run existing compare suite to verify no position
   regressions from the refactored `_pipeline_icrs()`.

---

### Phase 3: Tighten tolerances and update tests

**Goal**: adjust test tolerances to reflect the new precision, ensuring we
don't regress.

**Files to modify**:
- `tests/test_leb/test_leb_precision.py` (standalone precision tests)
- `tests/test_leb/compare/conftest.py` (compare suite tolerances)

**Tolerance changes** (after Phase 1 + 2):

| Metric | Before | After | Rationale |
|--------|--------|-------|-----------|
| `POSITION_ARCSEC` (ICRS planets) | 5.0" | 1.0" | Uranus now <0.5", others <<1" |
| `ECLIPTIC_ARCSEC` | 0.5" | 0.5" | Unchanged (already tight) |
| `SPEED_LON_DEG_DAY` | 0.05 | 0.005 | Analytical velocity ~10x better |
| `SPEED_LAT_DEG_DAY` | 0.05 | 0.005 | Same |
| `EQUATORIAL_ARCSEC` | 5.0" | 1.0" | Follows position improvement |
| `J2000_ARCSEC` | 5.0" | 1.0" | Follows position improvement |
| `SIDEREAL_ARCSEC` | 20.0" | 5.0" | Conservative (some modes have larger errors) |

**Note**: these are conservative initial targets. After measuring actual errors
post-implementation, they can be tightened further. The principle is: set
tolerances at ~2x the measured worst-case error, so tests catch regressions
without being brittle.

**Additional test work**:
- Add the 6 missing test classes from the compare implementation plan
- Unskip station tests (import from `libephemeris.crossing` directly)
- Triage the 12 failing crossing tests (mark as `xfail` with clear reason —
  see section 2.6)

---

### Phase 4: Regenerate, validate, document

**Goal**: produce the new `.leb` files and validate everything end-to-end.

**Steps**:

1. **Regenerate medium tier**:
   ```bash
   poe leb:generate:medium:groups
   ```
   This runs planets, asteroids, and analytical groups as separate subprocesses
   then merges them. Expected time: 10-30 minutes depending on hardware.

2. **Post-generation verification**:
   ```bash
   python scripts/generate_leb.py --tier medium --verify --verify-samples 1000
   ```
   All 30 bodies should show PASS. Record the per-body max errors.

3. **Run precision tests**:
   ```bash
   poe test:leb:precision:quick   # standalone (generates fresh LEB)
   poe test:leb:compare           # compare suite (uses data/leb/ephemeris_medium.leb)
   ```

4. **Record and compare error before/after** for every body. Create a table
   showing the improvement.

5. **Update documentation**:
   - `docs/leb/guide.md`: update all precision tables, parameter rationale,
     file size estimates, Moon speed discussion (section 10.3 — no longer
     relevant after Phase 2)
   - `docs/leb/design.md`: update section 7.3 to document the analytical
     velocity approach (option (a) from the original design)

6. **Regenerate other tiers** (if applicable):
   ```bash
   poe leb:generate:base:groups      # base tier (de440s, 1850-2150)
   poe leb:generate:extended:groups  # extended tier (de441, -5000 to 5000)
   ```

---

## 4. Execution Order and Dependencies

```
Phase 0 (asteroid test fix)  ──→  (independent, can run anytime)

Phase 1 (parameters)  ──→  Phase 4a (regenerate)  ──→  Phase 3 (tolerances)
                                                              ↓
Phase 2 (velocity)  ───────────────────────────────→  Phase 4b (validate)
```

Phase 0 is **fully independent** — it only touches test files and can be done
at any point. It should be done first to unblock a clean base tier test run.

Phase 1 and Phase 2 are **independent** — they can be implemented in parallel
or in either order. However:

- Phase 4 (regeneration) requires Phase 1 to be done first
- Phase 3 (tolerances) should be done last, after measuring actual errors
- Phase 2 should be validated with both old and new `.leb` files

**Recommended execution order**: Phase 0 → Phase 1 → Phase 2 → Phase 4 → Phase 3

---

## 5. Risks and Mitigations

### 5.1 Velocity transform correctness

**Risk**: the Cartesian-to-spherical velocity conversion or the rotation matrix
application to velocity vectors could have subtle sign errors or coordinate
convention mismatches.

**Mitigation**: extensive comparison with the existing central-difference
results (which are correct within their numerical precision). Before removing
the central difference code, run both methods in parallel and compare. Keep the
central difference code as a commented-out fallback until validation is complete.

### 5.2 Light-time velocity interaction

**Risk**: using `retarded_vel` instead of the velocity at the nominal time may
introduce small inconsistencies with Skyfield's approach.

**Mitigation**: this is the physically correct choice (velocity at emission
time), matching what Swiss Ephemeris does. The difference vs velocity at
reception time is ~v²/c ≈ 1e-8 for the Moon, negligible.

### 5.3 LEB file backwards compatibility

**Risk**: changing BODY_PARAMS means old `.leb` files become incompatible with
the new code (the reader uses BODY_PARAMS to know interval/degree).

**Mitigation**: the interval and degree are stored **per-body in the file
header** (`BodyEntry` in `leb_format.py`), not read from BODY_PARAMS at
runtime. The reader uses the values from the file. So old files will continue
to work — they'll just have the old (less precise) parameters. BODY_PARAMS
is only used by the **generator** to decide what parameters to use when
creating new files.

### 5.4 File size growth

**Risk**: ~11 MB growth (88 → 99 MB) might concern users.

**Mitigation**: 12.5% growth for a measurable precision improvement is a
good tradeoff. The file fits comfortably in memory. If size is critical,
the base tier (300 years) grows proportionally less.

### 5.5 Venus/Mars latitude residual error

**Risk**: even after all improvements, Venus latitude may still show ~0.3-0.4"
error due to the inherent 1/distance amplification in the ICRS→ecliptic
pipeline.

**Mitigation**: this is an architectural constraint, not a bug. Documenting it
clearly in the precision tables is sufficient. To fully eliminate it would
require storing ecliptic coordinates directly (losing the ability to serve
multiple output frames from one dataset), which is not worth the tradeoff.

### 5.6 Asteroid date filtering may reduce test coverage

**Risk**: filtering dates for asteroids means fewer dates are tested per
asteroid. If the SPK coverage is very narrow, some asteroids might have very
few test dates.

**Mitigation**: the filter logs the number of dates tested. If fewer than 10
dates remain after filtering, the test should still run (not skip) — even a
small number of dates validates precision within the coverage range. The
`pytest.skip()` only triggers when zero dates remain. The real validation for
asteroid precision is the medium tier, which has much wider SPK coverage.

---

## 6. File Reference

### Modified files

| File | Phase | Changes |
|------|-------|---------|
| `tests/test_leb/compare/base/conftest.py` | 0 | Add `get_body_date_range()` and `filter_dates_for_body()` helpers |
| `tests/test_leb/compare/base/test_base_asteroids.py` | 0 | Filter dates per asteroid body before iteration |
| `tests/test_leb/compare/base/test_base_distances.py` | 0 | Filter dates for Chiron/Ceres in geocentric distance test |
| `tests/test_leb/compare/base/test_base_velocities.py` | 0 | Filter dates for asteroids in distance velocity test |
| `libephemeris/leb_format.py` | 1 | Update BODY_PARAMS for Moon, Earth, Uranus, Neptune, Pluto |
| `scripts/generate_leb.py` | 1 | Update NUTATION_INTERVAL from 32 to 16 |
| `libephemeris/fast_calc.py` | 2 | Add `_cartesian_velocity_to_spherical()`, modify `_pipeline_icrs()` to return velocity, remove central difference block |
| `tests/test_leb/test_leb_precision.py` | 3 | Update tolerance constants |
| `tests/test_leb/compare/conftest.py` | 3 | Update Tolerances dataclass defaults |
| `docs/leb/guide.md` | 4 | Update precision tables, parameter tables, rationale, file size examples |
| `docs/leb/design.md` | 4 | Update velocity computation section |

### Read-only reference files

| File | Relevance |
|------|-----------|
| `libephemeris/leb_reader.py:82-143` | Chebyshev derivative computation (`_clenshaw_with_derivative`) |
| `libephemeris/leb_reader.py:272-315` | `eval_body()` — already returns (pos, vel), raises `ValueError` for out-of-range JD |
| `libephemeris/leb_format.py:102-113` | `BodyEntry` dataclass — has `jd_start`, `jd_end` fields |
| `libephemeris/fast_calc.py:409-495` | `_pipeline_icrs()` — current implementation |
| `libephemeris/fast_calc.py:777-803` | `_fast_calc_core()` — central difference block to remove |
| `libephemeris/fast_calc.py:186-230` | `_precession_nutation_matrix()` and `_mat3_vec3()` — reuse for velocity |
| `libephemeris/fast_calc.py:110-128` | `_rotate_equatorial_to_ecliptic()` — reuse for velocity |
| `libephemeris/fast_calc.py:138-179` | `_apply_aberration()` — reference for aberration velocity decision |
| `scripts/generate_leb.py:182-184` | Nutation generation constants |
| `scripts/generate_leb.py:2404-2428` | Error reporting during generation |
| `scripts/generate_leb.py:2691-2855` | Post-generation verification (`verify_leb()`) |
| `docs/leb/guide.md:section 10.4` | Documents asteroid SPK coverage limitations |
