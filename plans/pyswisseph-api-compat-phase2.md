# PySwissEph API Compatibility — Phase 2

## Context

Phase 1 (commit `4305239`) fixed fixstar return types, rise_trans signatures, and added
missing constants. A follow-up batch (`b39b3cb`..`b6d8795`) by another agent fixed all
~100 missing constant aliases and several eclipse wrapper signatures.

This plan covers the **remaining** signature / return-type / default-value mismatches
discovered during the full re-audit on 2026-03-21.

---

## P0 — Bug (we introduced)

### 1. `fixstar_mag` returns bare `float` instead of `(float, str)`

**File:** `libephemeris/fixed_stars.py:5520`

**Current:**
```python
def swe_fixstar_mag(star_name: str) -> float:
    ...
    return _STAR_MAGNITUDES[star_id]
```

**pyswisseph:** `fixstar_mag(star)` → `(float mag, str stnam)` (confirmed via docstring
and identical to `fixstar2_mag`)

**Fix:** Change return to `(magnitude, star_name_out)` like `swe_fixstar2_mag` does.
Need `_format_star_name` or equivalent for the canonical name. Use `_canonical_name`
from `_resolve_star_id` to build `"Name,Nomenclature"` format.

**Callers to update:**
- `libephemeris/eclipse.py` — grep for `fixstar_mag` calls
- `libephemeris/heliacal.py` — grep for `fixstar_mag` calls
- `tests/test_fixed_stars/test_fixstar_mag.py` — tests that unpack bare float
- Any compare_scripts that call `fixstar_mag`

---

## P1 — Missing default values

### 2. `calc_ut` / `calc` — `iflag` is required

**File:** `libephemeris/planets.py:794,905`

**Current:**
```python
def swe_calc_ut(tjd_ut: float, ipl: int, iflag: int) -> ...:
def swe_calc(tjd: float, ipl: int, iflag: int) -> ...:
```

**pyswisseph:** `flags=FLG_SWIEPH|FLG_SPEED` (= `2|256` = `258`)

**Fix:** Change `iflag: int` → `iflag: int = SEFLG_SWIEPH | SEFLG_SPEED` on both.

**Callers:** No callers break (adding a default is backwards-compatible).

### 3. `calc_pctr` — `iflag` is required

**File:** `libephemeris/planets.py:986`

**Current:** `def swe_calc_pctr(tjd_ut: float, ipl: int, iplctr: int, iflag: int)`

**pyswisseph:** `flags=FLG_SWIEPH|FLG_SPEED`

**Fix:** `iflag: int = SEFLG_SWIEPH | SEFLG_SPEED`

### 4. `pheno` / `pheno_ut` — `iflag` is required

**File:** `libephemeris/planets.py:5002,5042`

**Current:** `def swe_pheno_ut(tjd_ut: float, ipl: int, iflag: int)`

**pyswisseph:** `flags=FLG_SWIEPH`

**Fix:** `iflag: int = SEFLG_SWIEPH`

### 5. `set_topo` — `alt` is required

**File:** `libephemeris/state.py:668`

**Current:** `def set_topo(lon: float, lat: float, alt: float)`

**pyswisseph:** `set_topo(lon, lat, alt=0.0)`

**Fix:** `alt: float = 0.0`

---

## P2 — Wrong parameter order / missing params

### 6. `nod_aps_ut` — swapped param order + no defaults

**File:** `libephemeris/planets.py:3827`

**Current:**
```python
def swe_nod_aps_ut(tjd_ut: float, ipl: int, iflag: int, method: int)
```

**pyswisseph:** `nod_aps_ut(tjdut, planet, method=NODBIT_MEAN, flags=FLG_SWIEPH|FLG_SPEED)`
(method 3rd, flags 4th, both with defaults)

**Fix:** Change to:
```python
def swe_nod_aps_ut(
    tjd_ut: float, ipl: int,
    method: int = SE_NODBIT_MEAN,
    iflag: int = SEFLG_SWIEPH | SEFLG_SPEED,
)
```

**Callers to update:** Every call site that passes `iflag` as 3rd positional and `method`
as 4th positional must be swapped. Search:
- `libephemeris/` for `swe_nod_aps_ut(` and `nod_aps_ut(`
- `tests/` for same
- `compare_scripts/` for same
- `examples/` for same

### 7. `nod_aps` — wrong default for `iflag`, `method` has no default

**File:** `libephemeris/planets.py:3876`

**Current:**
```python
def swe_nod_aps(tjd_et: float, ipl: int, method: int, iflag: int = SEFLG_SPEED)
```

**pyswisseph:** `nod_aps(tjdet, planet, method=NODBIT_MEAN, flags=FLG_SWIEPH|FLG_SPEED)`

**Fix:** Change to:
```python
def swe_nod_aps(
    tjd_et: float, ipl: int,
    method: int = SE_NODBIT_MEAN,
    iflag: int = SEFLG_SWIEPH | SEFLG_SPEED,
)
```

**Callers:** Search for `swe_nod_aps(` — existing callers that pass `method` positionally
still work.

### 8. `cs2timestr` — missing `sep` and `suppresszero` params

**File:** `libephemeris/utils.py:1186`

**Current:** `def cs2timestr(cs: int) -> str`

**pyswisseph:** `cs2timestr(cs, sep, suppresszero=False)` where `sep` is a 1-byte
separator character. Returns format like `"12h34m56s"` with the sep char.

**Fix:** Add `sep` and `suppresszero` params:
```python
def cs2timestr(cs: int, sep: str = ":", suppresszero: bool = False) -> str:
```
- `sep` replaces the `:` separator (pyswisseph uses bytes like `b':'`, we accept str)
- `suppresszero` suppresses trailing zero components

The current implementation hardcodes `:` separator in the f-string. Change to use `sep`.
If `suppresszero` is True and seconds==0, omit `:00`; if minutes also==0, omit `:00:00`.

**Note:** pyswisseph accepts `bytes` for sep. We should accept both `str` and `bytes`:
```python
if isinstance(sep, bytes):
    sep = sep.decode("ascii")
```

### 9. `houses` / `houses_ex` / `houses_ex2` — `hsys` has no default

**Files:**
- `libephemeris/houses.py:510` — `swe_houses(tjdut, lat, lon, hsys, iflag=0)`
- `libephemeris/houses.py:1514` — `swe_houses_ex(tjdut, lat, lon, hsys, flags=0)`
- `libephemeris/houses.py:1627` — `swe_houses_ex2(tjdut, lat, lon, hsys, flags=0)`

**pyswisseph:**
- `houses(tjdut, lat, lon, hsys=b'P')` — no `flags` param at all
- `houses_ex(tjdut, lat, lon, hsys=b'P', flags=0)`
- `houses_ex2(tjdut, lat, lon, hsys=b'P', flags=0)`

**Fix:**
- `swe_houses`: add default `hsys: int = ord("P")`. The `iflag` param is a libephemeris
  extension (pyswisseph `houses` has no flags). Keep it but with default `0` so
  `houses(jd, lat, lon)` works.
- `swe_houses_ex`: add default `hsys: int = ord("P")`
- `swe_houses_ex2`: add default `hsys: int = ord("P")`

**Note on `hsys` type:** pyswisseph uses `bytes` (`b'P'`). Our `hsys` is `int` (`ord('P')`).
Our functions already handle both via `_normalize_hsys()` or similar. The default should
be `ord("P")` since our type annotation is `int`. Users passing `b'P'` already works if
the function converts it.

### 10. `houses_armc` — missing `ascmc9` param

**File:** `libephemeris/houses.py:1081`

**Current:** `def swe_houses_armc(armc, lat, eps, hsys)`

**pyswisseph:** `houses_armc(armc, lat, eps, hsys=b'P', ascmc9=0.0)`

**Fix:**
```python
def swe_houses_armc(armc: float, lat: float, eps: float,
                    hsys: int = ord("P"), ascmc9: float = 0.0)
```
The `ascmc9` param is used for Sunshine house system; it can be accepted and stored but
may not affect our implementation if Sunshine doesn't use it. Add the param for
signature compatibility; if unused internally, that's fine.

### 11. `houses_armc_ex2` — has `flags` but pyswisseph has `ascmc9`

**File:** `libephemeris/houses.py:1372`

**Current:** `def swe_houses_armc_ex2(armc, lat, eps, hsys, flags=0)`

**pyswisseph:** `houses_armc_ex2(armc, lat, eps, hsys=b'P', ascmc9=0.0)`

**Fix:** pyswisseph's `houses_armc_ex2` does NOT have a `flags` param — it has `ascmc9`.
Change to:
```python
def swe_houses_armc_ex2(armc: float, lat: float, eps: float,
                        hsys: int = ord("P"), ascmc9: float = 0.0)
```
Since our implementation uses `flags` internally for speed calculation, we can keep
`flags` as an internal detail or accept it as a keyword-only arg for backwards compat:
```python
def swe_houses_armc_ex2(armc: float, lat: float, eps: float,
                        hsys: int = ord("P"), ascmc9: float = 0.0,
                        *, flags: int = 0)
```

### 12. `helio_cross` / `helio_cross_ut` — missing `backwards` param

**Files:**
- `libephemeris/crossing.py:1065` — `swe_helio_cross_ut(planet_id, x2cross, jd_ut, flag=SEFLG_SWIEPH)`
- `libephemeris/crossing.py:1228` — `swe_helio_cross(planet_id, x2cross, jd_tt, flag=SEFLG_SWIEPH)`

**pyswisseph:** `helio_cross(planet, x2cross, tjdet, flags=FLG_SWIEPH, backwards=False)`

**Fix:** Add `backwards: bool = False` param. When `True`, search backward in time
(negate the search step).

**Implementation:** The internal search loop needs to reverse direction. Read the current
implementation to see how the search step works, then negate it when `backwards=True`.

### 13. `date_conversion` / `swe_date_conversion` — missing defaults

**Files:**
- `libephemeris/time_utils.py:293` — `date_conversion(year, month, day, hour, calendar)`
- `libephemeris/__init__.py:592` — `swe_date_conversion(year, month, day, hour, calendar)`

**pyswisseph:** `date_conversion(year, month, day, hour=12.0, cal=b'g')`

**Fix — `swe_date_conversion`:** Add defaults:
```python
def swe_date_conversion(year: int, month: int, day: int,
                        hour: float = 12.0, calendar: "str | bytes" = b"g")
```

**Fix — `date_conversion`:** This is a different function (converts between calendars,
different return type). It is NOT the pyswisseph-compatible one. The bare `date_conversion`
should be an alias of `swe_date_conversion` for pyswisseph compatibility. Currently
`__init__.py` imports `date_conversion` from `time_utils` (different function) AND defines
`swe_date_conversion` inline. We need to make `date_conversion = swe_date_conversion` in
`__init__.py`, and keep the old one available under a different name if needed.

### 14. `swe_sol_eclipse_where` — `ifl` has no default

**File:** `libephemeris/eclipse.py:2704`

**Current:** `def swe_sol_eclipse_where(tjd_ut: float, ifl: int)`

**pyswisseph:** `sol_eclipse_where(tjdut, flags=FLG_SWIEPH)`

**Fix:** `ifl: int = SEFLG_SWIEPH`

### 15. `vis_limit_mag` — returns `int` as first element, pyswisseph returns `float`

**File:** `libephemeris/heliacal.py:2354` (and duplicate in `eclipse.py:9036`)

**Current:** `-> Tuple[int, Tuple[float, ...]]`

**pyswisseph:** Returns `(float res, tuple[10] dret)` — first element is float
(values -2, 0, 1, 2)

**Fix:** Change return type annotation to `Tuple[float, Tuple[float, ...]]` and ensure
the first element is returned as `float(result)` not `int`.

---

## P3 — Bare-name wrappers should alias swe_ versions

In pyswisseph, `sol_eclipse_how` IS `swe_sol_eclipse_how` — same function. In
libephemeris, the bare names are **different** functions with "Pythonic" signatures
(separate lat/lon/alt instead of geopos tuple, etc.).

### Strategy

Make bare names into aliases of the `swe_` versions. The old "Pythonic" functions can be
kept under different names (e.g. `sol_eclipse_how_geo`, or just internal) if any internal
code depends on them, but the PUBLIC bare names must match pyswisseph.

### Functions to reassign in `__init__.py`:

| Bare name | Currently points to | Should point to |
|-----------|--------------------|-----------------| 
| `sol_eclipse_how` | `eclipse.sol_eclipse_how` (Pythonic) | `eclipse.swe_sol_eclipse_how` |
| `sol_eclipse_when_glob` | `eclipse.sol_eclipse_when_glob` (Pythonic) | `eclipse.swe_sol_eclipse_when_glob` |
| `sol_eclipse_when_loc` | `eclipse.sol_eclipse_when_loc` (Pythonic) | `eclipse.swe_sol_eclipse_when_loc` |
| `sol_eclipse_where` | `eclipse.sol_eclipse_where` (Pythonic) | `eclipse.swe_sol_eclipse_where` |
| `lun_eclipse_how` | `eclipse.lun_eclipse_how` (Pythonic) | `eclipse.swe_lun_eclipse_how` |
| `lun_eclipse_when` | `eclipse.lun_eclipse_when` (Pythonic) | `eclipse.swe_lun_eclipse_when` |
| `lun_eclipse_when_loc` | `eclipse.lun_eclipse_when_loc` (Pythonic) | `eclipse.swe_lun_eclipse_when_loc` |
| `lun_occult_when_loc` | `eclipse.lun_occult_when_loc` (Pythonic) | `eclipse.swe_lun_occult_when_loc` |
| `heliacal_ut` | `heliacal.heliacal_ut` (Pythonic) | `heliacal.swe_heliacal_ut` |
| `heliacal_pheno_ut` | `heliacal.heliacal_pheno_ut` (Pythonic) | `heliacal.swe_heliacal_pheno_ut` |
| `gauquelin_sector` | `houses.gauquelin_sector` (Pythonic) | `houses.swe_gauquelin_sector` |

**Impact:** Any internal code or tests that call the bare names with the old Pythonic
signatures will break. Need to grep for all call sites.

**Callers to update:**
- `tests/` — any test using `sol_eclipse_how(jd, lat, lon, alt)` style
- `compare_scripts/` — same
- `examples/` — same
- Internal library code (e.g., `eclipse.py` itself, `heliacal.py`)

**Note:** The old Pythonic functions should remain in eclipse.py/heliacal.py/houses.py
as internal helpers (rename with `_` prefix or `_pythonic` suffix), but NOT be exported
as the public bare names.

---

## P4 — Minor

### 16. `version` attribute

**File:** `libephemeris/__init__.py:758`

**Current:** Has `__version__ = "0.25.0"` and `swe_version()` function, but no `version`
attribute.

**pyswisseph:** `swe.version` returns `'2.10.03'` (module-level string attribute)

**Fix:** After `__version__` definition, add:
```python
version = __version__
```

---

## Execution Order

1. **P0** (fixstar_mag) — fix the bug we introduced
2. **P1** (defaults) — backwards-compatible, no callers break
3. **P4** (version) — trivial
4. **P2.14** (sol_eclipse_where default) — backwards-compatible
5. **P2.15** (vis_limit_mag return type) — may need test updates
6. **P2.8** (cs2timestr) — new params with defaults, backwards-compatible
7. **P2.9** (houses hsys defaults) — backwards-compatible
8. **P2.10-11** (houses_armc) — backwards-compatible if ascmc9 is optional
9. **P2.6-7** (nod_aps_ut order) — **BREAKING**: callers must be updated
10. **P2.12** (helio_cross backwards) — new param with default, backwards-compatible
11. **P2.13** (date_conversion) — **BREAKING**: bare name changes return type
12. **P3** (bare-name aliases) — **BREAKING**: callers using Pythonic signatures break

---

## Test Strategy

After each breaking change:
1. Run the specific test files that cover the changed functions
2. Fix any broken call sites
3. Run `pytest tests/ -x --timeout=30` to catch regressions (skip slow tests)

For P3 (bare-name aliases), the impact is large. Search for all callers first:
```bash
grep -rn 'sol_eclipse_how\|sol_eclipse_when_glob\|sol_eclipse_when_loc' tests/ examples/ compare_scripts/ --include='*.py'
grep -rn 'lun_eclipse_how\|lun_eclipse_when\b\|lun_eclipse_when_loc' tests/ examples/ compare_scripts/ --include='*.py'
grep -rn 'heliacal_ut\|heliacal_pheno_ut' tests/ examples/ compare_scripts/ --include='*.py'
grep -rn 'gauquelin_sector' tests/ examples/ compare_scripts/ --include='*.py'
```

---

## Files to Modify (Summary)

### Library files:
- `libephemeris/fixed_stars.py` — fixstar_mag return type
- `libephemeris/planets.py` — calc/calc_ut/calc_pctr/pheno/pheno_ut defaults, nod_aps/nod_aps_ut order+defaults
- `libephemeris/state.py` — set_topo alt default
- `libephemeris/utils.py` — cs2timestr sep+suppresszero
- `libephemeris/houses.py` — houses/houses_ex/houses_ex2 hsys default, houses_armc ascmc9, houses_armc_ex2 ascmc9
- `libephemeris/crossing.py` — helio_cross/helio_cross_ut backwards param
- `libephemeris/eclipse.py` — swe_sol_eclipse_where ifl default, vis_limit_mag return float
- `libephemeris/heliacal.py` — vis_limit_mag return float
- `libephemeris/__init__.py` — version attribute, bare-name alias reassignments, date_conversion alias, swe_date_conversion defaults

### Test/script files (to be updated for breaking changes):
- Tests calling `nod_aps_ut` with old param order
- Tests calling bare-name eclipse/heliacal/gauquelin with Pythonic signatures
- Tests checking `fixstar_mag` bare float return
- Tests checking `vis_limit_mag` int return
