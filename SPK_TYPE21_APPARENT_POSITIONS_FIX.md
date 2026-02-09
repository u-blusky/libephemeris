# Bug Report: SPK Type 21 Returns J2000 Geometric Instead of Apparent Positions

**Bug ID:** `LIBEPH-SPK21-PRECESSION`  
**Severity:** High  
**Affected Version:** 0.9.0  
**Status:** ✅ **FIXED** in commit `b7b13ac`  
**Fixed Version:** 0.9.1  
**Date:** 2026-02-09  
**Author:** AI Assistant  

---

## Executive Summary

libephemeris correctly calculates **geometric J2000 positions** from SPK Type 21 files, but astrology requires **apparent positions** (equinox of date with aberration). This causes a discrepancy of ~0.33° (~20 arcminutes) compared to Swiss Ephemeris and JPL Horizons for bodies like Chiron.

---

## Problem Description

### Symptoms

When calculating Chiron's position using SPK Type 21 kernel:

| Source | Longitude (2024-01-15 00:00 UTC) | Reference Frame |
|--------|----------------------------------|-----------------|
| libephemeris SPK | 15.2662° | J2000 geometric |
| Swiss Eph (J2000+NOABERR) | 15.2652° | J2000 geometric |
| Swiss Eph (default) | 15.5988° | Apparent (of-date) |
| JPL Horizons Observer | 15.6157° | Apparent observer |
| **Difference** | **~0.33° (1258 arcsec)** | |

### Root Cause Analysis

The SPK file contains coordinates in the **ICRF/J2000 equatorial frame**. The current implementation in `_calc_type21_position()`:

1. ✅ Correctly converts ICRS → ecliptic J2000 (rotation by obliquity)
2. ✅ Correctly converts heliocentric → geocentric (subtracting Earth position)
3. ❌ **Does NOT apply precession** from J2000 equinox to equinox of date
4. ❌ **Does NOT apply aberration** (annual aberration due to Earth's orbital motion)
5. ❌ **Does NOT apply light-time correction** (position at retarded time)

### Why This Matters

- Astrology traditionally uses **apparent positions** (what you would observe in the sky)
- Swiss Ephemeris default output is apparent positions
- The ~0.33° error is significant for astrological interpretation (Chiron moves ~0.02°/day)
- This affects ALL bodies calculated via SPK Type 21 (Chiron, Pholus, TNOs, etc.)

---

## Evidence

### Test Script Output

```python
# JD 2460323.5 = 2024-01-15 00:00 UTC

# Swiss Eph with different flags:
J2000 + NOABERR + NOGDEFL      lon=15.265165°
J2000 + NOABERR                lon=15.265166°
J2000 only                     lon=15.264405°
NOABERR only                   lon=15.599585°
Default (apparent)             lon=15.598825°
NONUT                          lon=15.600086°

# libephemeris SPK Type 21 result:
libephemeris SPK               lon=15.266175°  # Matches J2000+NOABERR!
```

### Breakdown of the ~0.33° Difference

| Correction | Contribution | Applied? |
|------------|-------------|----------|
| Precession (J2000 → 2024) | ~0.334° | ❌ No |
| Annual Aberration | ~0.005° | ❌ No |
| Nutation | ~0.003° | ❌ No |
| Light-time | ~0.0003° | ❌ No |
| **Total** | **~0.34°** | |

---

## Affected Code

### File: `libephemeris/spk.py`

### Function: `_calc_type21_position()` (line 818-930)

```python
def _calc_type21_position(
    kernel,
    naif_id: int,
    t,
    iflag: int,
) -> Optional[tuple]:
    """
    Calculate body position using SPK type 21 kernel.
    
    CURRENT ISSUE: Returns J2000 geometric positions instead of apparent.
    """
    # ... current implementation ...
    
    # Line 861: Converts to ecliptic J2000 (correct)
    pos_ecl = np.array(_icrs_to_ecliptic_j2000(*pos_au))
    
    # Line 887: Computes geocentric (correct)
    pos_geo = pos_ecl - earth_helio_ecl
    
    # Line 901: Converts to spherical (correct)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    
    # MISSING: Precession from J2000 to date
    # MISSING: Aberration correction
    # MISSING: Light-time correction
    # MISSING: Nutation (optional)
```

---

## Available Resources

The required transformation functions **already exist** in libephemeris:

### Precession

| Function | Location | Description |
|----------|----------|-------------|
| `precess_from_j2000(lon, lat, jd_tt)` | `moshier/precession.py:677` | Precess ecliptic coords J2000→date |
| `precession_matrix_j2000_to_date(jd_tt)` | `moshier/precession.py:436` | 3x3 precession matrix |

### Aberration

| Function | Location | Description |
|----------|----------|-------------|
| `apply_aberration_to_position(pos, earth_vel)` | `moshier/utils.py:481` | Apply aberration to position vector |
| `annual_aberration_cartesian(target_dir, earth_vel)` | `moshier/utils.py:423` | Cartesian aberration correction |

### Nutation

| Function | Location | Description |
|----------|----------|-------------|
| `nutation_angles(jd_tt)` | `moshier/precession.py:335` | IAU 2000B nutation Δψ, Δε |

### Constants

| Constant | Value | Location |
|----------|-------|----------|
| `C_LIGHT_AU_DAY` | 173.1446326846693 | `moshier/utils.py:24` |
| `SEFLG_J2000` | 2048 | `constants.py:575` |
| `SEFLG_NOABERR` | 1024 | `constants.py:574` |
| `SEFLG_NONUT` | 64 | `constants.py:570` |
| `SEFLG_TRUEPOS` | 16 | `constants.py:568` |

---

## Implementation Plan

### Overview

Transform the current calculation flow:

```
CURRENT:
  SPK ICRS → ecliptic J2000 → geocentric J2000 → spherical

REQUIRED:
  SPK ICRS → ecliptic J2000 → light-time → geocentric J2000 
           → aberration → precession → nutation → spherical
```

### Step 1: Add Light-Time Correction

**Location:** `_calc_type21_position()`, after line 851

**Implementation:**
```python
from .constants import SEFLG_TRUEPOS

# Constants
C_LIGHT_AU_DAY = 173.1446326846693

# Get initial position
pos_km, vel_km = kernel.compute_type21(10, naif_id, jd)

# Apply iterative light-time correction (unless SEFLG_TRUEPOS)
if not (iflag & SEFLG_TRUEPOS):
    for _ in range(3):  # 3 iterations is standard
        pos_au = np.array(pos_km) / AU_KM
        pos_ecl = np.array(_icrs_to_ecliptic_j2000(*pos_au))
        
        # Get Earth position for geocentric distance
        earth_helio_ecl = ...  # existing code
        pos_geo = pos_ecl - earth_helio_ecl
        
        dist_au = np.linalg.norm(pos_geo)
        light_time_days = dist_au / C_LIGHT_AU_DAY
        jd_retarded = jd - light_time_days
        
        pos_km, vel_km = kernel.compute_type21(10, naif_id, jd_retarded)
```

**Respects flag:** `SEFLG_TRUEPOS` skips light-time correction

### Step 2: Add Aberration Correction

**Location:** After computing geocentric position, before converting to spherical

**Implementation:**
```python
from .constants import SEFLG_NOABERR

# Apply aberration (unless SEFLG_NOABERR)
if not (iflag & SEFLG_NOABERR):
    # Get Earth's velocity in ecliptic J2000
    earth_vel_icrs = earth.at(t).velocity.au_per_d - sun.at(t).velocity.au_per_d
    earth_vel_ecl = np.array(_icrs_to_ecliptic_j2000(*earth_vel_icrs))
    
    # Aberration formula: apparent = geometric + (v_earth / c) × direction
    # Using the relativistic aberration approximation
    c_au_day = 173.1446326846693
    direction = pos_geo / np.linalg.norm(pos_geo)
    
    # Bradley's classical aberration
    aberration = np.cross(earth_vel_ecl, np.cross(direction, earth_vel_ecl)) / c_au_day
    pos_geo = pos_geo + aberration * np.linalg.norm(pos_geo)
```

**Alternative:** Use existing `apply_aberration_to_position()` from `moshier/utils.py`

**Respects flag:** `SEFLG_NOABERR` skips aberration

### Step 3: Add Precession from J2000 to Date

**Location:** After computing spherical coordinates (lon, lat)

**Implementation:**
```python
from .constants import SEFLG_J2000
from .moshier.precession import precess_from_j2000

# Convert to spherical (existing code)
lon = math.degrees(math.atan2(y, x)) % 360.0
lat = math.degrees(math.asin(z / r))

# Apply precession from J2000 to equinox of date (unless SEFLG_J2000)
if not (iflag & SEFLG_J2000):
    jd_tt = t.tt  # Need TT for precession
    lon, lat = precess_from_j2000(lon, lat, jd_tt)
```

**Respects flag:** `SEFLG_J2000` returns J2000 positions

### Step 4: Add Nutation (Optional)

**Location:** After precession

**Implementation:**
```python
from .constants import SEFLG_NONUT
from .moshier.precession import nutation_angles

# Apply nutation (unless SEFLG_NONUT or SEFLG_J2000)
if not (iflag & SEFLG_NONUT) and not (iflag & SEFLG_J2000):
    delta_psi, delta_eps = nutation_angles(jd_tt)
    lon += delta_psi  # Nutation in longitude (degrees)
```

**Respects flag:** `SEFLG_NONUT` skips nutation

### Step 5: Update Speed Calculations

The velocity/speed calculations (lines 905-930) also need to be updated to reflect the precessed frame. This may require:

1. Precessing the velocity vector along with position
2. Or computing speed as finite difference in the apparent frame

---

## Testing Plan

### New Test File: `tests/test_spk/test_type21_apparent.py`

```python
"""
Tests for SPK Type 21 apparent position calculations.

Verifies that positions match Swiss Ephemeris within tolerance
for various flag combinations.
"""

import pytest
import libephemeris as eph
import swisseph as swe

# Test dates
JD_2024_01_15 = 2460323.5
JD_2000_01_01 = 2451545.0
JD_1990_01_01 = 2447892.5

class TestChironApparentPositions:
    """Test Chiron apparent positions match Swiss Ephemeris."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Ensure Chiron SPK is registered."""
        eph.download_and_register_spk("2060", eph.SE_CHIRON, "1800-01-01", "2200-12-31")
        swe.set_ephe_path("/path/to/sweph")
    
    @pytest.mark.parametrize("jd", [JD_2024_01_15, JD_2000_01_01, JD_1990_01_01])
    def test_apparent_matches_swisseph(self, jd):
        """Default (apparent) position should match Swiss Eph within 1 arcsec."""
        lib_result = eph.swe_calc_ut(jd, eph.SE_CHIRON, 0)
        swe_result = swe.calc_ut(jd, swe.CHIRON)
        
        diff_lon = abs(lib_result[0][0] - swe_result[0][0]) * 3600
        diff_lat = abs(lib_result[0][1] - swe_result[0][1]) * 3600
        
        assert diff_lon < 1.0, f"Longitude diff {diff_lon}\" exceeds 1 arcsec"
        assert diff_lat < 1.0, f"Latitude diff {diff_lat}\" exceeds 1 arcsec"
    
    def test_j2000_flag_returns_geometric(self):
        """SEFLG_J2000 should return J2000 geometric position."""
        jd = JD_2024_01_15
        
        lib_result = eph.swe_calc_ut(jd, eph.SE_CHIRON, eph.SEFLG_J2000)
        swe_result = swe.calc_ut(jd, swe.CHIRON, swe.FLG_J2000)
        
        diff_lon = abs(lib_result[0][0] - swe_result[0][0]) * 3600
        assert diff_lon < 5.0, f"J2000 longitude diff {diff_lon}\" exceeds 5 arcsec"
    
    def test_noaberr_flag(self):
        """SEFLG_NOABERR should skip aberration correction."""
        jd = JD_2024_01_15
        
        lib_result = eph.swe_calc_ut(jd, eph.SE_CHIRON, eph.SEFLG_NOABERR)
        swe_result = swe.calc_ut(jd, swe.CHIRON, swe.FLG_NOABERR)
        
        diff_lon = abs(lib_result[0][0] - swe_result[0][0]) * 3600
        assert diff_lon < 2.0, f"NOABERR longitude diff {diff_lon}\" exceeds 2 arcsec"
    
    def test_apparent_differs_from_j2000(self):
        """Apparent position should differ from J2000 by ~0.33° in 2024."""
        jd = JD_2024_01_15
        
        apparent = eph.swe_calc_ut(jd, eph.SE_CHIRON, 0)
        geometric = eph.swe_calc_ut(jd, eph.SE_CHIRON, eph.SEFLG_J2000)
        
        diff_deg = abs(apparent[0][0] - geometric[0][0])
        
        # Precession from J2000 to 2024 is about 0.33°
        assert 0.30 < diff_deg < 0.40, f"Apparent-J2000 diff {diff_deg}° not in expected range"
```

### Integration Tests

1. Run full kerykeion test suite after fix
2. Verify Chiron positions in chart SVGs match expected
3. Test with Pholus and other SPK bodies

---

## Verification Matrix

| Test Case | Current | Expected | Flag |
|-----------|---------|----------|------|
| Chiron default | 15.266° | 15.599° | 0 |
| Chiron J2000 | 15.266° | 15.264° | SEFLG_J2000 |
| Chiron NOABERR | 15.266° | 15.600° | SEFLG_NOABERR |
| Chiron J2000+NOABERR | 15.266° | 15.265° | SEFLG_J2000 \| SEFLG_NOABERR |

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Regression in existing tests | Medium | High | Run full test suite before merge |
| Performance impact | Low | Low | Light-time iteration adds 2 SPK lookups |
| Incorrect precession | Low | High | Validate against Swiss Eph for multiple dates |
| Flag handling errors | Medium | Medium | Test each flag combination |

---

## Estimated Effort

| Phase | Time |
|-------|------|
| Implementation | 3-4 hours |
| Unit Testing | 1-2 hours |
| Integration Testing | 1 hour |
| Documentation | 30 min |
| **Total** | **5-7 hours** |

---

## Dependencies

1. `moshier/precession.py` functions must work correctly (✅ verified)
2. `moshier/utils.py` aberration functions must work correctly (needs verification)
3. kerykeion tests must pass after fix

---

## Appendix: Reference Implementation

### How Swiss Ephemeris Handles This

From Swiss Ephemeris source (`sweph.c`):

1. **Light-time:** Iterative correction (3 iterations)
2. **Aberration:** Annual aberration using Earth velocity
3. **Precession:** IAU 2006 precession from J2000 to date
4. **Nutation:** IAU 2000B nutation series
5. **Frame rotation:** Ecliptic of date

### How Skyfield Handles This

Skyfield's `.apparent()` method applies:
1. Light-time correction (iterative)
2. Gravitational deflection (Sun)
3. Aberration (annual + diurnal)

But Skyfield returns ICRS positions, not ecliptic of date.

---

## Changelog Entry

```markdown
## [0.9.1] - 2026-02-XX

### Fixed
- SPK Type 21 calculations now return apparent positions (equinox of date)
  instead of J2000 geometric positions, matching Swiss Ephemeris behavior
- Added light-time correction for SPK Type 21 bodies
- Added precession from J2000 to equinox of date
- Added annual aberration correction
- Respects SEFLG_J2000, SEFLG_NOABERR, SEFLG_NONUT, SEFLG_TRUEPOS flags
```
