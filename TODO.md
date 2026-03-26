# TODO - Improvements and Technical Debt

Items identified during development (March 2026). Grouped by category.

---

## Next Up (v1.0.0a4+)

### ~~6. PyPI Packaging — ship `base_core.leb` in the wheel~~ DONE

Completed: `libephemeris/data/leb2/base_core.leb` included in package-data,
`_discover_leb_file()` falls back to bundled file.

### 7. Update download command for LEB2 modular files

**Priority: high | Effort: medium**

The `libephemeris download:leb:*` commands still download monolithic LEB1 files. Update to download LEB2 modular files (core, asteroids, apogee, uranians) from GitHub Releases.

- Update `libephemeris/download.py` DATA_FILES with LEB2 file URLs
- Support `libephemeris download:leb2:base:core`, `download:leb2:base:asteroids`, etc.
- Keep backward compat: old `download:leb:base` still works (downloads LEB1)

### ~~8. Regenerate LEB1 base with optimized BODY_PARAMS~~ DONE

Already regenerated — LEB1 has Pluto 64d/deg11, Uranians 256d/deg7.

### 9. Publish LEB2 non-core files to GitHub Releases

**Priority: medium | Effort: small**

LEB2 files for asteroids, apogee, uranians (all 3 tiers) need to be published as GitHub Release assets so users can download them.

- Use `gh release create` or update existing release
- Update download URLs in download.py

### 10. Run validation suite (415K checks)

**Priority: medium | Effort: small (CPU time)**

Execute `plans/validation-suite-100k.md` — the 415K-check validation plan covering positions, flags, velocity, houses, sidereal modes, edge cases, compression, Horizons, transforms, stars, eclipses, performance, and cross-backend consistency.

- Create a script that implements the checklist
- Run and document results

### 11. Known bugs from LEB optimization findings

**Priority: medium | Effort: varies**

Documented in `proposals/leb-optimization-findings.md`:

- ~~Uranian geocentric bodies (40-47)~~ FIXED — added geocentric path in `planets.py`.
- Sun heliocentric — `swe_calc(jd, SE_SUN, SEFLG_SPEED | SEFLG_HELCTR)` returns ~180 deg error. Known Skyfield path bug.
- True Node distance — exceeds DISTANCE_AU tolerance in comparison tests.

---

## Technical Debt (from code quality audit)

---

## 1. Thread Safety (architectural)

**Priority: high | Risk: high | Effort: large**

~30 mutable global variables in `state.py` are read/written without locks.
The current `_CONTEXT_SWAP_LOCK` serialises context-based calculations but
does not protect individual setter/getter pairs.

### What needs to happen

- Introduce a `threading.Lock` (or `RLock`) around every setter/getter in
  `state.py` (`set_topo`, `set_sid_mode`, `set_ephe_path`, `get_planets`,
  `get_timescale`, etc.).
- Alternatively, move to a fully context-based architecture where each
  `EphemerisContext` instance holds its own copy of the state, eliminating
  shared globals entirely.
- The `close()` function also needs a lock to avoid racing with in-flight
  calculations.

### Why it was skipped

Requires an architectural decision (lock-per-variable vs context-only) and
risks breaking the pyswisseph-compatible global API that many callers rely on.

### Files involved

- `libephemeris/state.py` (primary)
- `libephemeris/context.py`

---

## 2. Broad `except Exception` Blocks (101 sites)

**Priority: medium | Risk: low | Effort: medium**

101 `except Exception:` blocks across the codebase silently swallow errors.
24 are in `eclipse.py` alone. During the audit, logging was added to the 5
most critical ones (velocity fallbacks in `planets.py`, state init paths in
`state.py`). The remaining ~96 still silently discard exceptions.

### What needs to happen

- Audit each remaining site and decide:
  - Expected fallback path -> `logger.debug()`
  - Genuinely unexpected error -> `logger.warning()` or narrower exception type
- Highest-value files to address first:
  - `eclipse.py` (24 sites)
  - `houses.py` (~10 sites)
  - `heliacal.py` (~8 sites)
- Consider narrowing `except Exception` to specific exception types where the
  failure mode is well-understood.

### Why it was skipped

Most of these are in complex numerical calculation paths where the broad catch
is arguably intentional (graceful degradation). Changing them requires
understanding each individual call site's failure modes.

### Files involved

- `libephemeris/eclipse.py` (24 sites)
- `libephemeris/houses.py` (~10 sites)
- `libephemeris/heliacal.py` (~8 sites)
- `libephemeris/crossing.py`, `libephemeris/hypothetical.py`, and others

---

## 3. `MeeusRangeError` Exception Hierarchy

**Priority: low | Risk: medium | Effort: small**

`MeeusRangeError` inherits from `ValueError` but semantically should inherit
from `EphemerisRangeError` (which inherits from `Error`, not `ValueError`).
Additionally, `MeeusRangeError` is dead code: it is defined and exported in
`__all__` but never raised or caught anywhere in the library.

### What needs to happen

- Decide: either remove `MeeusRangeError` entirely (it is unused), or fix its
  hierarchy to inherit from `EphemerisRangeError`.
- If changing the hierarchy: update `tests/test_meeus_validity.py` which
  asserts `issubclass(MeeusRangeError, ValueError)`.

### Why it was skipped

The test explicitly asserts the current hierarchy, and `MeeusRangeError` is in
the public `__all__`. Changing it could break downstream code that catches
`ValueError` expecting to catch `MeeusRangeError`.

### Files involved

- `libephemeris/exceptions.py`
- `tests/test_meeus_validity.py`

---

## 4. `houses.py` Code Deduplication

**Priority: low | Risk: medium | Effort: medium**

`houses.py` contains repeated calculation blocks across different house system
implementations. Similar to the sidereal correction deduplication done in
`planets.py` (Fix 6), common patterns could be extracted into helper functions.

### What needs to happen

- Identify repeated blocks (pole handling, cusp normalization, etc.).
- Extract into shared helpers.
- Requires running the full house system test suite to verify no regressions,
  which is slow.

### Why it was skipped

Too risky without being able to run the full test suite. House system
calculations are sensitive to floating-point order of operations.

### Files involved

- `libephemeris/houses.py`

---

## 5. `constants.py` Missing `__all__`

**Priority: low | Risk: high | Effort: small**

`constants.py` does not define `__all__`, so `from .constants import *` exports
every name in the module. Adding `__all__` would make the public API explicit.

### What needs to happen

- Define `__all__` in `constants.py` listing all public constants.
- Verify that all 51 files doing `from .constants import *` still work.

### Why it was skipped

51 files use `from .constants import *`. Any constant accidentally omitted from
`__all__` would silently break those imports at runtime. Requires a complete
audit of which constants are actually used.

### Files involved

- `libephemeris/constants.py`
- 51 files that do `from .constants import *` or `from libephemeris.constants import *`
