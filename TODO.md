# TODO - Improvements and Technical Debt

Items identified during development (March 2026). Grouped by category.

---

## Completed

- ~~#6~~ PyPI Packaging — `base_core.leb` bundled in wheel, auto-discovered
- ~~#7~~ Download command — `download_leb2_for_tier()` in download.py, 12 files in DATA_FILES
- ~~#8~~ LEB1 regenerated — Pluto 64d/deg11, Uranians 256d/deg7
- ~~#9~~ GitHub Release data-v2 — 12 LEB2 files published
- ~~#11a~~ Uranian geocentric — added geocentric path in `planets.py`

---

## Open

### 10. Run validation suite (415K checks)

**Priority: medium | Effort: CPU time**

Execute `plans/validation-suite-100k.md`. Create a script, run, document results.

### 11b. Sun heliocentric bug

**Priority: medium | Effort: small**

`swe_calc(jd, SE_SUN, SEFLG_HELCTR)` via Skyfield path returns ~180 deg error.
The Sun heliocentric should be (0,0,0) or Earth's heliocentric position reversed.
Documented in `proposals/leb-optimization-findings.md`.

**Files:** `libephemeris/planets.py` (_calc_body Sun + SEFLG_HELCTR path)

### 11c. True Node distance tolerance

**Priority: low | Effort: small**

True Node distance exceeds DISTANCE_AU tolerance in LEB comparison tests.
Likely a rounding/precision issue in the eccentricity vector computation.

---

## Technical Debt

### 1. Thread Safety

**Priority: high | Risk: high | Effort: large**

~30 mutable global variables in `state.py` without locks. Requires architectural
decision: lock-per-variable vs context-only.

**Files:** `state.py`, `context.py`

### 2. Broad `except Exception` (101 sites)

**Priority: medium | Risk: low | Effort: medium**

101 bare `except Exception:` blocks. 24 in `eclipse.py`, ~10 in `houses.py`,
~8 in `heliacal.py`. Most are intentional graceful degradation.

### 3. `MeeusRangeError` — dead code

**Priority: low | Risk: medium | Effort: small**

Defined and exported but never raised. Either remove or fix hierarchy
(`ValueError` -> should be `EphemerisRangeError`).

### 4. `houses.py` deduplication

**Priority: low | Risk: medium | Effort: medium**

Repeated calculation blocks across house systems.

### 5. `constants.py` missing `__all__`

**Priority: low | Risk: high | Effort: small**

51 files do `from .constants import *`. Adding `__all__` requires full audit.
