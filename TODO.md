# TODO - Improvements and Technical Debt

Items identified during development (March 2026). Grouped by category.

---

## Completed

- ~~#6~~ PyPI Packaging — `base_core.leb` bundled in wheel, auto-discovered
- ~~#7~~ Download command — `download_leb2_for_tier()` in download.py, 12 files in DATA_FILES
- ~~#8~~ LEB1 regenerated — Pluto 64d/deg11, Uranians 256d/deg7
- ~~#9~~ GitHub Release data-v2 — 12 LEB2 files published
- ~~#10~~ Validation suite — `scripts/run_validation_suite.py` created and run (7/8 PASS, 1 preexisting XYZ numpy bug)
- ~~#11a~~ Uranian geocentric — added geocentric path in `planets.py`
- ~~#11b~~ Sun heliocentric — returns (0,0,0) correctly now

---

## Open

### 11c. True Node distance tolerance

**Priority: low | Effort: unclear**

True Node distance differs from reference by ~2.3e-4 AU. Longitude is correct (~4").
Different osculating orbit computation from Moon state vectors.
Documented in `docs/reference/known-bugs.md`.

---

## Technical Debt

### 1. Thread Safety

**Priority: high | Risk: high | Effort: large**

~30 mutable global variables in `state.py` without locks. Requires architectural
decision: lock-per-variable vs context-only.

**Files:** `state.py`, `context.py`

### 2. Broad `except Exception` (~92 sites remaining)

**Priority: medium | Risk: low | Effort: medium**

`state.py` done (9 -> 0). Remaining: 24 in `eclipse.py`, ~8 in `planets.py`,
~3 in `houses.py`, ~8 in `heliacal.py`. Most are intentional graceful
degradation in numerical iteration loops.

### ~~3. `MeeusRangeError` — dead code~~ DONE

Marked as deprecated in docstring, kept for backward compat.

### 4. `houses.py` deduplication

**Priority: low | Risk: medium | Effort: medium**

Repeated calculation blocks across house systems.

### 5. `constants.py` missing `__all__`

**Priority: low | Risk: high | Effort: small**

51 files do `from .constants import *`. Adding `__all__` requires full audit.
