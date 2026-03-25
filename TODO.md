# TODO - Improvements and Technical Debt

Items identified during the code quality audit (March 2026) that were
deliberately skipped because they require deeper architectural work or carry
backward-compatibility risk.

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
