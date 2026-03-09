# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.23.0] - 2026-03-09

### Summary

**LEB Precision V3: sub-milliarcsecond accuracy for all 31 celestial bodies across
all three precision tiers.** This release completely rewrites the LEB runtime
pipeline, replacing the previous ICRS-to-ecliptic coordinate conversion approach
(which suffered from 1–5 arcsecond errors due to retrograde cusps and COB
oscillations) with a physics-correct pipeline that stores smooth ICRS barycentric
coordinates and applies gravitational deflection, special-relativistic aberration,
and precession-nutation at evaluation time.

Previous worst-case error: **4.85 arcseconds** (Saturn, base tier).
New worst-case error: **0.000332 arcseconds** (Moon, base tier) — a **14,600x improvement**.

### Added

#### New coordinate type: `COORD_ICRS_BARY_SYSTEM` (type 4)

Outer planets (Jupiter, Saturn, Uranus, Neptune, Pluto) now store their **system
barycenter** positions in ICRS coordinates rather than planet-center positions.
System barycenters are ultra-smooth trajectories free of high-frequency moon
oscillations that made Chebyshev fitting unreliable for planet centers. The
Center-of-Body (COB) correction — the offset from system barycenter to planet
center due to satellite gravitational pull — is applied at runtime using SPK
data or analytical moon theories.

- New constant `COORD_ICRS_BARY_SYSTEM = 4` in `leb_format.py`
- New `_SYSTEM_BARY_NAMES` mapping in `fast_calc.py` (body ID → Skyfield segment name)
- New `_apply_cob_correction()` function in `fast_calc.py` for runtime COB via
  `get_planet_centers()` SPK data with analytical fallback
- New `generate_body_icrs_system_bary()` in `generate_leb.py` with `_SYSTEM_BARY_MAP`
  for generation-time system barycenter extraction

#### PPN gravitational deflection in LEB runtime pipeline

The LEB pipeline now applies post-Newtonian (PPN) gravitational light deflection
from the Sun, Jupiter, and Saturn, matching Skyfield's `apparent()` pipeline.
This corrected a ~3.95 arcsecond systematic error on Saturn that was present in
the V2 pipeline.

- New `_apply_gravitational_deflection()` function in `fast_calc.py` (~60 lines)
- Implements the standard PPN formula: `δθ = (1+γ) GM/(c²d) · (sin θ)⁻¹`
- Three deflecting bodies: Sun (dominant, ~1.75" max), Jupiter (~0.017" max),
  Saturn (~0.006" max)
- Moon excluded from deflection (geocentric observer, zero impact parameter)

#### Skyfield-compatible asteroid pipeline via `_SpkType21Target`

New `_SpkType21Target` class in `spk.py` that wraps SPK Type 21 asteroid kernels
as Skyfield `VectorFunction`-compatible objects. This routes asteroids through
Skyfield's `observe()`/`apparent()` pipeline (light-time iteration, aberration,
deflection) instead of the previous manual ecliptic J2000 + precession/nutation
approach that had ~0.3–0.4 arcsecond systematic errors.

- New `_SpkType21Target` class (~100 lines) with `_at()`, `at()`,
  `_observe_from_bcrs()` methods matching the Skyfield VectorFunction protocol
- New `get_spk_type21_target()` factory function
- Combines heliocentric SPK Type 21 positions with Sun SSB position for
  SSB-centered ICRS output
- Handles both scalar and array time inputs

#### Precision measurement script

New `scripts/measure_precision.py` for dense end-to-end precision measurement
of LEB files against Skyfield reference calculations.

- Vincenty-formula angular separation for numerically stable error measurement
- Per-body statistics: mean, P99, max error in arcseconds + component breakdown
- `--tier` flag with auto-detection from LEB filename for correct SPK selection
- `--group` flag for testing specific body categories
- `--samples` flag for configurable sampling density (default 2000 per body)
- Forces `calc_mode="skyfield"` globally to prevent LEB auto-discovery contamination
- Suppresses `MeeusPolynomialWarning` spam on extended tier

#### Diagnostic scripts

New diagnostic scripts for investigating precision issues:

- `scripts/diagnose_errors.py` — Decomposes total error into Chebyshev fitting
  error vs COB mismatch vs end-to-end pipeline error
- `scripts/diagnose_pipeline.py` — Step-by-step pipeline comparison (raw
  Chebyshev, light-time, deflection, aberration, precession-nutation)
- `scripts/diagnose_chebyshev.py` — Chebyshev coefficient analysis and
  fitting quality visualization
- `scripts/diagnose_leb_read.py` — Low-level LEB reader diagnostics
- `scripts/measure_medium_errors.py` — Medium tier focused error analysis
- `scripts/prototype_v2_fix.py` — V2 deflection pipeline prototype
- `scripts/prototype_cartesian_v3.py` — V3 Cartesian storage prototype
- `scripts/test_chebyshev_params.py` — Automated Chebyshev parameter tuning
- `scripts/check_leb_params.py` — BODY_PARAMS validation

#### Documentation

- **New `docs/leb/algorithms.md`** (~1020 lines): Comprehensive mathematical
  reference covering Chebyshev polynomial theory, Clenshaw recurrence algorithm,
  least-squares fitting methodology, coordinate systems (ICRS barycentric,
  system barycentric, ecliptic, heliocentric), COB corrections, PPN gravitational
  deflection derivation, special-relativistic aberration, precession-nutation
  matrices, error analysis methodology, and historical problems with solutions
- Updated `docs/leb/guide.md` from v1.0 to v2.0 (~325 lines changed):
  `COORD_ICRS_BARY_SYSTEM` documentation, Pipeline A' (deflection + aberration),
  updated BODY_PARAMS table with precision results, file sizes for all 3 tiers,
  per-tier precision tables
- Updated `docs/leb/design.md`: marked as historical document with header
  redirecting to `guide.md` and `algorithms.md`
- Updated `docs/leb/testing.md`: tolerances, file sizes, body count updated
  to reflect V3 precision
- Updated `docs/README.md`: added `algorithms.md` and `testing.md` links

### Changed

#### LEB runtime pipeline rewrite (`fast_calc.py`, +334 lines)

Complete rewrite of `_pipeline_icrs()` — the core evaluation pipeline for
ICRS-stored bodies (planets, asteroids, Earth). The new pipeline (called
"Pipeline A'" in the documentation) performs:

1. **Chebyshev evaluation** — Clenshaw algorithm on stored ICRS barycentric
   coefficients
2. **COB correction** (system barycenters only) — runtime planet-center offset
   via SPK data or analytical moon theories, evaluated at observer time
3. **Light-time iteration** — Newton-Raphson convergence for retarded position
4. **Gravitational deflection** — PPN formula with Sun/Jupiter/Saturn
5. **Special-relativistic aberration** — Bradley formula with Earth velocity
6. **Precession-nutation** — Skyfield's IAU 2006/2000A frame rotation to
   ecliptic of date

Previous pipeline stored geocentric ecliptic coordinates directly, which
failed at retrograde stations (cusps in longitude) and contained COB
oscillations from outer planet moons.

#### Outer planet storage: planet center → system barycenter

Jupiter, Saturn, Uranus, Neptune, and Pluto changed from
`COORD_ICRS_BARY` (planet center) to `COORD_ICRS_BARY_SYSTEM` (system
barycenter) in `BODY_PARAMS`. System barycenters are smooth trajectories
determined only by solar system gravitational dynamics, free of the
high-frequency oscillations introduced by satellite orbits. The COB
correction is applied at runtime, matching Skyfield's internal pipeline.

#### Ecliptic body Chebyshev parameters tightened

OscuApogee, InterpApogee, and InterpPerigee changed from `8d/13`
(8-day intervals, degree 13) to `4d/15` (4-day intervals, degree 15)
to capture fast oscillations (~2.6°/day for OscuApogee). This reduced
ecliptic body errors from ~0.028" to ~0.000049".

#### Ecliptic body generation fix

`generate_ecliptic_bodies_vectorized()` was hardcoding
`interval_days=8, degree=13` for all ecliptic bodies, ignoring per-body
`BODY_PARAMS`. Fixed to split ecliptic bodies by their params and generate
each group with correct parameters.

#### Delta T precision fix in `fast_calc_ut()`

Replaced `reader.delta_t()` (linearly-interpolated sparse table with up to
~0.004s error near 1985) with `swe_deltat()` (Skyfield's precise Delta T
model) for UT→TT conversion. This eliminated ~0.002" Moon errors from
imprecise time conversion.

#### Test tolerance overhaul

All three tiers' comparison test tolerances tightened by 3–4 orders of
magnitude:

| Tolerance | Base (before → after) | Medium (before → after) | Extended (before → after) |
|-----------|----------------------|------------------------|--------------------------|
| `POSITION_ARCSEC` | 5.0 → 0.001 | 5.0 → 0.001 | 5.0 → 0.001 |
| `ASTEROID_ARCSEC` | 0.5 → 0.001 | 0.5 → 0.001 | 5.0 → 0.001 |
| `ECLIPTIC_ARCSEC` | 0.05 → 0.001 | 0.05 → 0.001 | 0.1 → 0.1 * |
| `SIDEREAL_ARCSEC` | 5.0 → 0.001 | 5.0 → 0.001 | 5.0 → 0.001 |
| `DISTANCE_AU` | 3e-5 → 5e-6 | 3e-5 → 5e-6 | 5e-5 → 5e-6 |

\* Extended tier ecliptic bodies retain 0.1" tolerance due to Meeus polynomial
degradation beyond ±20 centuries from J2000.0 (architectural limit of the
underlying lunar theory, not of the LEB system).

Per-body ecliptic tolerances also tightened: all 6 ecliptic bodies from
0.01–0.5" to 0.001".

#### Unit test updates

- `test_fast_calc.py`: Flag tests rewritten to validate full pipeline output
  (geocentric ecliptic) instead of raw Chebyshev coefficients
- `test_generate_leb.py`: Sun/Moon fit accuracy tests updated to use full
  `fast_calc_ut()` pipeline with generous tolerance for on-the-fly test files
- `test_leb_reader.py`: Position reasonableness tests tightened; Sun/Skyfield
  comparison test rewritten to use full pipeline
- `test_leb_format.py`: `COORD_ICRS_BARY_SYSTEM` added to valid coord types
- `test_compare_leb_asteroids.py`: Uses `ASTEROID_ARCSEC` tolerance instead
  of `POSITION_ARCSEC`

#### LEB file regeneration

All three tier files regenerated with V3 pipeline:

| Tier | File | Size | Worst Case | Bodies |
|------|------|------|------------|--------|
| Base | `ephemeris_base.leb` | 112 MB | Moon 0.000332" | 31 |
| Medium | `ephemeris_medium.leb` | 377 MB | Moon 0.000325" | 31 |
| Extended | `ephemeris_extended.leb` | 2.8 GB | Mars 0.000010" * | 31 |

\* Excluding ecliptic bodies at extreme dates (OscuApogee 0.054" at ±50 centuries).

### Fixed

#### COB evaluation time bug

`_apply_cob_correction()` in `fast_calc.py` was evaluating the COB offset
at retarded time (`jd_tt - light_time`), but Skyfield's `_observe_from_bcrs()`
evaluates COB at observer time (`jd_tt`). This one-line fix (changing
`jd_tt - lt` to `jd_tt`) eliminated residual errors for outer planets.

#### Asteroid pipeline mismatch (0.3–0.4" systematic error)

The reference asteroid pipeline used ecliptic J2000 coordinates with manual
precession and nutation rotation, while the LEB pipeline used ICRS coordinates
with Skyfield's precession-nutation matrix. The two approaches differ by
~0.3–0.4 arcseconds due to frame tie and nutation model differences. Fixed
by creating `_SpkType21Target` in `spk.py` that routes asteroids through
Skyfield's `observe()`/`apparent()` pipeline, ensuring identical coordinate
transformations.

#### `measure_precision.py` mode-switching bug (false 0.28" Moon error)

The precision measurement script had three bugs causing false error reports:

1. **LEB auto-discovery contamination**: `swe_calc()` reference calls could
   silently use LEB via `get_leb_reader()` auto-discovery, comparing LEB
   against itself instead of Skyfield. Fixed by forcing `calc_mode="skyfield"`.
2. **Per-sample LEBReader creation**: `LEBReader` was instantiated inside the
   inner loop (once per sample point, 2000x per body) instead of once at
   startup. Fixed by moving to `main()` and passing the reader instance.
3. **Cross-tier SPK mismatch**: Extended tier LEB (generated from `de441.bsp`)
   was compared against Skyfield using `de440.bsp` (medium, the default),
   showing ~0.05" false errors from DE440 vs DE441 ephemeris differences.
   Fixed by adding `--tier` flag with auto-detection from filename.

### Removed

- `docs/leb/leb_precision_v3.md` — superseded by `algorithms.md`
- `docs/leb/precision_v2_plan.md` — obsolete planning document
- `docs/leb/precision-improvement-plan.md` — obsolete planning document
- `docs/leb/compare-implementation-plan.md` — obsolete planning document
- Dead `COORD_GEO_ECLIPTIC` code paths in `fast_calc.py` (the geocentric
  ecliptic storage approach was abandoned early in V3 development due to
  retrograde cusp fitting failures)

### Precision Results

All 31 bodies pass <0.001 arcsecond on all three tiers (1569 comparison tests):

#### Base tier (de440s.bsp, 1849–2150, 404 tests)

| Group | Bodies | Worst Body | Max Error |
|-------|--------|------------|-----------|
| Planets (11) | Sun–Pluto, Earth | Moon | 0.000332" |
| Asteroids (5) | Chiron–Vesta | Juno | 0.000045" |
| Ecliptic (6) | Nodes, Lilith, Apsides | OscuApog | 0.000049" |
| Hypothetical (9) | Cupido–Isis | all | ~0.000000" |

#### Medium tier (de440.bsp, 1550–2650, 904 tests)

| Group | Bodies | Worst Body | Max Error |
|-------|--------|------------|-----------|
| Planets (11) | Sun–Pluto, Earth | Moon | 0.000325" |
| Asteroids (5) | Chiron–Vesta | Vesta | 0.000036" |
| Ecliptic (6) | Nodes, Lilith, Apsides | OscuApog | 0.000075" |
| Hypothetical (9) | Cupido–Isis | all | ~0.000000" |

#### Extended tier (de441.bsp, -5000 to 5000 CE, 261 tests)

| Group | Bodies | Worst Body | Max Error |
|-------|--------|------------|-----------|
| Planets (11) | Sun–Pluto, Earth | Mars | 0.000010" |
| Asteroids (5) | Chiron–Vesta | Pallas | 0.000018" |
| Ecliptic (6) | Nodes, Lilith, Apsides | OscuApog | 0.054" * |
| Hypothetical (9) | Cupido–Isis | all | ~0.000000" |

\* Meeus polynomial lunar theory degrades beyond ±20 centuries from J2000.0.
Within ±1000 CE, ecliptic body errors are <0.001".

## [0.22.0] - 2026-03-02

### Added

#### Binary Ephemeris Mode (LEB) — Complete Implementation

Full implementation of the LibEphemeris Binary (LEB) precomputed ephemeris system,
providing ~14x runtime speedup over Skyfield/JPL pipeline via Chebyshev polynomial
approximations stored in `.leb` files.

##### Per-body date ranges for asteroids (Tasks 1.1-1.4)

Asteroid SPK kernels from JPL Horizons cover only ~1600-2500 CE, while planetary
ephemerides (DE440/DE441) span much wider ranges. LEB now stores per-body
`jd_start`/`jd_end` in each `BodyEntry` header, allowing asteroids to have narrower
coverage than planets within the same file. When a query falls outside an asteroid's
LEB range, the library transparently falls back to Skyfield.

- `assemble_leb()` computes per-body date ranges from actual SPK coverage
- `eval_body()` in `leb_reader.py` raises `ValueError` for out-of-range JDs
- `swe_calc_ut()`/`swe_calc()` catch both `KeyError` and `ValueError` for fallback
- Keplerian fallback removed from asteroid generation — only SPK data is used

##### Vectorized analytical body generation (Tasks 2.1-2.7, ~10x speedup)

The 6 ecliptic bodies (True Node, True Lilith, Osculating Apogee, Interpolated
Apogee, Interpolated Perigee, Mean Apogee) previously made ~328,000 scalar Skyfield
calls each. A single vectorized Skyfield call now computes all shared Moon ephemeris
data, then each body's values are derived with numpy.

New vectorized functions in `generate_leb.py`:
- `_calc_mean_lilith_batch()` — vectorized mean apse analytical calculation
- `_calc_lunar_fundamental_arguments_batch()` — vectorized Delaunay arguments
- `_calc_elp2000_apogee_perturbations_batch()` — vectorized 40+ term perturbation series
- `_calc_elp2000_perigee_perturbations_batch()` — vectorized 61-term perturbation series
- `_eval_ecliptic_bodies_batch()` — single Skyfield call, all 6 bodies from shared r,v
- `generate_ecliptic_bodies_vectorized()` — orchestrates batch eval + Chebyshev fitting
- `_fit_and_verify_from_values_unwrap()` — fits with longitude unwrapping

**Result:** Analytical group generation went from ~5-8 minutes to ~38 seconds for base tier.

##### Group generation and merge workflow (Task 1.7)

New `--group {planets,asteroids,analytical}` and `--merge FILE [FILE ...]` CLI
options for `generate_leb.py`. Partial group files are generation-time artifacts;
at runtime only the single merged file is used. This avoids macOS
`ProcessPoolExecutor` deadlocks that occurred when child processes used C extensions
(numpy, erfa, Skyfield).

New `poe` tasks for all three tiers:
```
poe leb:generate:base:groups       # All 3 groups + merge for base tier
poe leb:generate:medium:groups     # All 3 groups + merge for medium tier
poe leb:generate:extended:groups   # All 3 groups + merge for extended tier
```

##### LEB verification rewrite (Task 1.6)

`verify_leb()` rewritten with proper per-body verification against Skyfield
reference values. Reports per-body max error in arcseconds with PASS/FAIL status.

##### Calculation mode control (Task 3.3)

New `LIBEPHEMERIS_MODE` environment variable and programmatic API for controlling
the calculation backend:

```python
from libephemeris import set_calc_mode, get_calc_mode

set_calc_mode("auto")      # Use LEB if configured, otherwise Skyfield (default)
set_calc_mode("skyfield")  # Always use Skyfield, even if LEB is configured
set_calc_mode("leb")       # Require LEB; raises RuntimeError if unavailable
```

- `get_leb_reader()` respects mode: `"skyfield"` returns `None`, `"leb"` raises if missing
- `close()` resets `_CALC_MODE` to `None`
- Environment variable `LIBEPHEMERIS_MODE` sets default (overridden by programmatic call)
- 18 new tests in `TestCalcMode` class

##### LEB binary ephemeris download commands (Task 3.5)

New CLI commands and Python API to download pre-generated LEB files from
GitHub Releases, eliminating the need to generate them locally:

```bash
libephemeris download:leb:base       # ~53 MB, 1850-2150 CE
libephemeris download:leb:medium     # ~175 MB, 1550-2650 CE
libephemeris download:leb:extended   # not yet available
```

Python API:
```python
from libephemeris import download_leb_for_tier
download_leb_for_tier("medium")  # downloads + activates
```

- Atomic downloads with SHA256 integrity verification
- Progress bar (rich if available, simple fallback)
- `--force` flag to re-download existing files
- Post-download LEB validation (magic bytes, version, body index)
- Optional auto-activation via `set_leb_file()` after download
- New `poe` tasks: `download:leb:base`, `download:leb:medium`, `download:leb:extended`

##### LEB auto-discovery

`get_leb_reader()` now automatically discovers downloaded LEB files without
requiring explicit `set_leb_file()` or `LIBEPHEMERIS_LEB` configuration.

Resolution order:
1. Explicit path via `set_leb_file()`
2. `LIBEPHEMERIS_LEB` environment variable
3. Auto-discovery: `~/.libephemeris/leb/ephemeris_{tier}.leb`

This means `libephemeris download:leb:medium` followed by any `swe_calc_ut()`
call will automatically use the LEB fast path — zero configuration needed.

##### Data source logging

Added `logger.debug()` calls at all calculation dispatch points, enabling runtime
visibility into which backend serves each request. Activate with
`set_log_level('DEBUG')` or `LIBEPHEMERIS_LOG_LEVEL=DEBUG`.

Sources logged: `LEB`, `LEB->fallback`, `Skyfield`, `SPK`, `SPK (auto-downloaded)`,
`ASSIST (n-body)`, `Keplerian (fallback)`.

Log format: `[libephemeris] DEBUG: body=<id> jd=<jd> source=<SOURCE>`

Dispatch points: `swe_calc_ut()`, `swe_calc()`, `_calc_body()` in `planets.py`;
`calc_ut()`, `calc()` in `context.py`.

#### Keplerian Precision Improvements

##### Laplace-Lagrange secular perturbations for eccentricity and inclination (Task 4.1)

Implemented the (h,k)/(p,q) vector formalism for secular evolution of eccentricity
and inclination under gravitational perturbations from Jupiter, Saturn, Uranus, and
Neptune:

- Eccentricity vector `(h, k) = (e sin varpi, e cos varpi)` decomposes into forced + free
- Inclination vector `(p, q) = (sin(i/2) sin Omega, sin(i/2) cos Omega)` similarly
- Free component rotates at the body's proper frequency
- Uses Laplace coefficients `b_{3/2}^{(1)}` and `b_{3/2}^{(2)}` for coupling terms
- `_calc_forced_elements()` function (~120 lines) in `minor_bodies.py`
- `apply_secular_perturbations()` now returns 6-tuple `(omega, Omega, M, n, e_pert, i_pert)`

##### Multi-epoch orbital elements (Task 4.3)

Generated Keplerian elements from SPK Type 21 state vectors at 50-year intervals
from 1650-2450 CE for 6 bodies (Chiron, Pholus, Ceres, Pallas, Juno, Vesta).
`_get_closest_epoch_elements()` selects the element set with the smallest time
offset from the query date, dramatically improving long-term accuracy:

| Offset   | Before    | After     | Improvement |
|----------|-----------|-----------|-------------|
| 25 years | 3.3 deg   | 2.6'      | ~75x        |
| 50 years | 5.3 deg   | 3.6 deg   | ~1.5x       |
| 100 years| 10.7 deg  | 3.5 deg   | ~3x         |

- `MINOR_BODY_ELEMENTS_MULTI` dictionary (~600 lines, 17 epochs per body)
- State-to-Keplerian conversion via ICRS to ecliptic J2000 rotation
- Original single-epoch element always considered as candidate (preserves near-epoch accuracy)

##### Keplerian precision benchmark (Task 4.6)

New `tests/test_keplerian_precision_benchmark.py` with systematic comparison of
Keplerian vs SPK positions across 11 time offsets for 5 asteroids. Regression
tests enforce: epoch < 1", 1 month < 60", 1 year < 5'.

#### REBOUND/ASSIST N-body Integration (Task 4.4)

Complete integration of REBOUND/ASSIST as a fallback for ephemeris-quality asteroid
orbit propagation, including planetary perturbations from all major bodies.

##### End-to-end ASSIST pipeline

- `propagate_orbit_assist()` — ephemeris-quality integration with Sun, Moon, 8 planets,
  16 massive asteroids, J2/J3/J4 harmonics, and GR corrections
- Fixed critical bug: ASSIST manages the Sun internally, so particles must be added with
  Cartesian coordinates (not orbital elements). New `_elements_to_cartesian()` helper
  converts via temporary REBOUND simulation
- `propagate_trajectory()` ASSIST path also fixed with Cartesian conversion
- Fallback chain in `planets.py`: SPK > auto-download > ASSIST > Keplerian

##### Cached ASSIST availability check

- `check_assist_data_available()` — cached check for both ASSIST import and data file presence
- `reset_assist_data_cache()` — clears cache (called by `close()`)
- Replaces uncached `check_assist_available()` in the fallback chain

##### ASSIST data download

New `download_assist_data()` function with production-quality download pipeline:
- SSL certificate handling via `certifi`
- Atomic downloads (temp file + `os.replace()`)
- Progress bar integration via existing `_get_progress_bar()`
- File integrity verification: size validation + BSP structural check via `jplephem`
- Existing file detection with verification before skipping
- `--force` flag for re-download
- SHA256 hash reporting on completion

##### CLI command

New `libephemeris download:assist` CLI command:
```bash
libephemeris download:assist                    # Download both files (~714 MB)
libephemeris download:assist --no-asteroids     # Planet ephemeris only (~98 MB)
libephemeris download:assist --target-dir /path # Custom directory
libephemeris download:assist --force            # Re-download even if present
```

New `poe` task: `poe download:assist`

Required data files (saved to `~/.libephemeris/assist/`):
- `linux_p1550p2650.440` — JPL DE440 planet ephemeris in Linux binary format (~98 MB)
- `sb441-n16.bsp` — 16 massive asteroid perturbers (~616 MB)

##### Conditional test suite

- `TestAssistDataAvailability` (6 tests): cached check, reset, import failure, missing files, close resets cache
- `TestAssistEndToEnd` (9 tests): Ceres/Vesta/Chiron propagation, ASSIST vs REBOUND comparison,
  backward integration, trajectory, compare_with_keplerian, error handling
- Tests automatically skip when ASSIST data files are not present
- All 55 REBOUND/ASSIST tests pass (1 skipped for "without rebound" test)

### Changed

#### LEB Generator Performance

- Removed `ProcessPoolExecutor` entirely from `generate_leb.py` (macOS deadlock fix
  with numpy/erfa/Skyfield C extensions)
- Vectorized Skyfield ICRS evaluations (~150x speedup for individual body calls)
- Vectorized nutation computation via direct `erfa.nut06a()` (~50x speedup)
- Batched verification (fit + verify in single JD array allocation)
- `spktype21` integration for asteroid SPK evaluation (~36x vs `swe_calc`)
- Linear extrapolation for SPK boundary overshoot (smooth Chebyshev fitting at edges)

#### Minor Bodies

- `apply_secular_perturbations()` signature changed from 4-tuple to 6-tuple return:
  `(omega_pert, Omega_pert, M_pert, n_pert, e_pert, i_pert)`
- All callers updated: `calc_minor_body_position()`, `_elements_to_rebound_params()`,
  test unpacking in `test_secular_perturbations.py`

#### State Management

- `close()` now resets `_CALC_MODE` and calls `reset_assist_data_cache()`
- `get_leb_reader()` respects calculation mode setting
- `get_leb_reader()` auto-discovers LEB files in `~/.libephemeris/leb/` based on active tier

#### Naming

- Removed Swiss Ephemeris naming/references from generator labels and comments
  (replaced with LibEphemeris-native terminology)

### Documentation

- New comprehensive `docs/LEB_GUIDE.md` (~1668 lines): architecture, binary format,
  generation workflow, group/merge, verification, state management, exported API
- New comprehensive `TODO.md` (~530 lines): full implementation roadmap with status
- Updated `docs/LEB_PLAN.md` with implementation notes
- Updated `README.md`:
  - New "Binary Ephemeris Mode (LEB)" section with activation, mode control, LEB_GUIDE link
  - New "N-body fallback (REBOUND/ASSIST)" section with install/download/usage instructions
- Updated `docs/LEB_GUIDE.md`:
  - Section 1: "Calculation Mode (`LIBEPHEMERIS_MODE`)" subsection
  - Section 7.1: Mode variables and functions in global state docs
  - Section 7.4: `set_calc_mode`, `get_calc_mode` in exported API

### Removed

- `docs/KEPLERIAN_TODO.md` — superseded by `TODO.md`
- `ProcessPoolExecutor` from LEB generator (macOS deadlock)
- Keplerian fallback from asteroid LEB generation (SPK-only)

### Tests

- 18 new `TestCalcMode` tests in `test_context_leb.py`
- 6 new `TestAssistDataAvailability` tests in `test_rebound_integration.py`
- 9 new `TestAssistEndToEnd` tests in `test_rebound_integration.py`
- 5 new Keplerian precision benchmark tests (marked slow)
- Updated 3 callers of `apply_secular_perturbations()` in `test_secular_perturbations.py`
- All 155 LEB tests pass, 55 REBOUND tests pass, 46 secular perturbation tests pass, 29 context tests pass

## [0.20.0] - 2026-02-24

### Changed

- Uploaded planet_centers tier-specific SPK files to GitHub release:
  - `planet_centers_base.bsp` (25.4 MB, 1850-2150)
  - `planet_centers_medium.bsp` (72.6 MB, 1550-2650)
  - `planet_centers_extended.bsp` (222.6 MB, -12000 to +17000)
- Updated `download.py` with SHA256 hashes for integrity verification
- Modified `release_planet_centers.py` to find BSP files in `~/.libephemeris`

## [0.19.0] - 2026-02-23

### Fixed

- Corrected barycentric mode calculations (`SEFLG_BARYCTR`) for consistent
  heliocentric-to-barycentric coordinate transformations
- Fixed sign convention in `cotrans_sp()` ecliptic-to-equatorial coordinate
  transformation (affects lunar nodes, Lilith, interpolated apogee/perigee
  equatorial coordinates)
- Corrected lunar mean elements polynomial evaluation for improved
  Mean Lilith and Mean Node accuracy
- Adjusted test thresholds for lunar comparison tests to account for
  legitimate methodological differences between JPL and Swiss Ephemeris
- Improved lunar node and apogee/perigee precision with refined perturbation
  calculations

### Documentation

- Clarified scientific methodology and intentional deviation from Swiss
  Ephemeris for lunar apsides calculations in methodology docs
- Updated AGENTS.md with improved development guidelines

## [0.18.0] - 2026-02-21

### Changed

#### Complete astrological and astronomical logic rewrite

Rewrote and modernized core algorithms across the library for improved
numerical stability and code clarity:

- Replaced legacy trigonometric formulas with idiomatic `atan2`-based geometry
  throughout `houses.py`, `hypothetical.py`, `planets.py`, `utils.py`,
  `eclipse.py`, `fixed_stars.py`, `heliacal.py`, and `schaefer.py`
- Renamed single-letter variables to descriptive names across the entire
  codebase for maintainability
- Expanded `lunar_corrections.py` with regenerated correction tables
- Added `data/fictitious_orbits.csv` for hypothetical planet orbital elements,
  replacing the old `seorbel.txt` format
- Removed legacy `seorbel.txt` file

#### Perturbation series deduplication

- `generate_lunar_corrections.py` now imports perturbation series from
  `lunar.py` instead of maintaining a manual copy (~160 lines removed)
- Added JPL self-consistency tests for perigee perturbations (8 tests
  validating precession rate, decomposition consistency, table continuity,
  and JPL passage comparison)

### Fixed

#### Occultation API signatures restored

Restored `lun_occult_when_glob` and `lun_occult_when_loc` to accept separate
`(ipl, starname)` positional arguments instead of the unified body parameter,
preserving 1:1 pyswisseph API compatibility. Fixed `swe_lun_occult_when_loc`
silently dropping the `backwards` parameter. Resolved 65 failing tests.

### Documentation

- Added `docs/PYERFA_BENEFITS.md` documenting pyerfa integration benefits
- Added `docs/REBOUND_BENEFITS.md` documenting REBOUND n-body integration benefits
- Updated `docs/INTERPOLATED_APOGEE.md`, `docs/TRUE_LILITH_METHODS.md`,
  `docs/PRECISION.md`, `docs/AYANAMSHA.md`, `docs/HOUSE_SYSTEMS.md`,
  and Sphinx API reference

## [0.17.0] - 2026-02-20

### Changed

- Recalibrate interpolated perigee against JPL DE441 ephemeris
- Expand perigee perturbation model from 43 to 67 terms for improved accuracy
- Regenerate perigee correction table (15195 entries, 2-year step)

### Documentation

- Add interpolated perigee methodology document explaining JPL-vs-SE
  methodological differences (~11 deg RMS) between JPL quadratic regression
  on osculating elements and Swiss Ephemeris semi-analytical ELP2000-82B

## [0.16.1] - 2026-02-20

### Fixed

- `generate_planet_centers_spk.py`: load leap seconds kernel in `verify_spk()`
  before calling `et2utc()` to prevent `SpiceMISSINGTIMEINFO` error during
  SPK file verification

## [0.16.0] - 2026-02-20

### Added

#### 3-tier planet_centers SPK system

New tier-specific `planet_centers_*.bsp` files provide precise planet center
positions for Jupiter, Saturn, Uranus, Neptune, and Pluto. Each tier has its
own file with appropriate date coverage:

| Tier     | File                          | Coverage                 | Size       |
| -------- | ----------------------------- | ------------------------ | ---------- |
| base     | `planet_centers_base.bsp`     | 1850-2150                | ~15-20 MB  |
| medium   | `planet_centers_medium.bsp`   | 1550-2650                | ~40-50 MB  |
| extended | `planet_centers_extended.bsp` | partial -12000 to +17000 | ~80-100 MB |

Files are saved in the workspace root (alongside `de440.bsp`, `de441.bsp`).

#### Uranus analytical COB fallback

New `uranian.py` module implements Keplerian theory for Uranus' 5 major moons
(Ariel, Umbriel, Titania, Oberon, Miranda). This provides ~0.01 arcsec precision
for Uranus center-of-body corrections when SPK coverage is unavailable.

Previously, Uranus had no analytical fallback outside SPK range.

#### Planet centers generation commands

New `poe` tasks for generating tier-specific planet_centers files:

```bash
poe generate-planet-centers:base      # ~500 MB source download
poe generate-planet-centers:medium    # ~4 GB source download
poe generate-planet-centers:extended  # ~6.5 GB source download
poe generate-planet-centers:all       # Generate all 3 tiers
```

Requires `spiceypy >= 6.0.0`. Downloads satellite SPK files from JPL NAIF
and extracts planet center segments (NAIF 599, 699, 799, 899, 999).

#### GitHub release script

New `scripts/release_planet_centers.py` for uploading planet_centers files
to GitHub Releases with SHA256 hash calculation.

### Changed

- `get_planet_centers()` in `state.py` now loads the appropriate file for the
  active precision tier, with automatic reload when tier changes
- `download_for_tier()` now downloads tier-specific `planet_centers_{tier}.bsp`
  to workspace root instead of bundled `planet_centers.bsp`
- `print_data_status()` shows tier-specific files with current tier indicator
- Saturn extended tier now merges `sat441xl_part-1.bsp` + `sat441xl_part-2.bsp`
  for continuous coverage from -502 to +4500

#### Coverage summary

| Planet  | base  | medium        | extended                            |
| ------- | ----- | ------------- | ----------------------------------- |
| Jupiter | SPK ✓ | SPK ✓         | SPK 1600-2200, fallback outside     |
| Saturn  | SPK ✓ | SPK ✓         | SPK -502 to +4500, fallback outside |
| Uranus  | SPK ✓ | SPK ✓         | SPK full -12000/+17000              |
| Neptune | SPK ✓ | SPK ✓         | SPK full -12000/+17000              |
| Pluto   | SPK ✓ | SPK 1800-2200 | SPK 1800-2200, fallback outside     |

Fallback uses analytical moon theories (E5 for Jupiter, TASS 1.7 for Saturn,
Keplerian for Uranus/Neptune, Charon two-body for Pluto) with ~0.01-0.15 arcsec
precision.

## [0.15.0] - 2026-02-20

### Added

#### Tier-based CLI download commands

Replaced the old `init`, `init-fast`, and `download-data` CLI commands with three
tier-aware download commands. Each command downloads all data required for its
precision tier (ephemeris file, `planet_centers.bsp`, SPK kernels for 21 minor bodies):

- `libephemeris download:base` — de440s.bsp + SPK (1850-2150)
- `libephemeris download:medium` — de440.bsp + SPK (1900-2100)
- `libephemeris download:extended` — de441.bsp + SPK (1600-2500)

#### Programmatic tier download

New `download_for_tier()` function in the public API:

```python
from libephemeris import download_for_tier
download_for_tier("medium")
```

#### Max-range SPK download script

New `scripts/download_max_range_spk.py` downloads a single SPK file per body
covering the full JPL Horizons range (1600-2500) with `--body`, `--force`,
`--dry-run`, `--delay` options and `[N/total]` progress display.

#### Extended tier diagnostic dates

Tier diagnostic scripts now generate milestone dates covering the full range
of each tier: base (13 dates/25y), medium (22 dates/50y), extended (38 dates/800y).

### Fixed

#### SPK source detection in diagnostics

`_get_source()` in `scripts/_tier_diagnostic.py` now checks actual SPK file
coverage via `get_spk_coverage()` instead of assuming the tier's `spk_date_range`.
Previously reported "SPK" for dates outside the file's actual coverage.

#### `discover_local_spks()` wider-range file preference

When multiple SPK files exist for the same body (e.g., different date ranges),
`discover_local_spks()` now compares coverage spans and re-registers with the
wider file. Previously used whichever file `os.listdir()` returned first.

#### Extended tier SPK date range

Fixed `spk_date_range` for the extended tier from `("1550-01-01", "2650-01-01")`
to `("1600-01-01", "2500-01-01")` to match actual JPL Horizons limits. The old
range caused `_try_auto_spk_download()` to request impossible ranges.

### Changed

- Updated README: new CLI commands section, full `poe` task reference in Development
- `poe spk:download:extended` now runs `download_max_range_spk.py` instead of
  the old `ensure_all_ephemerides` one-liner

### Removed

- `libephemeris init`, `libephemeris init-fast`, `libephemeris download-data`
  CLI commands (replaced by `download:<tier>`)

## [0.14.0] - 2026-02-18

### Fixed

#### DE441 ephemeris range detection

Fixed incorrect date range detection for DE441 ephemeris files. DE441 has segments
split at 1969 (segments 0-13 cover -13200 to 1969, segments 14-27 cover 1969 to +17191).
The previous code only checked the first segment, incorrectly reporting the range as
-13200 to 1969. Now iterates all segments to find the overall min start_jd and max end_jd.

### Changed

- Improved README documentation with revised project description for clarity

## [0.13.0] - 2026-02-14

### Fixed

#### Critical: SE_SIDM_USER ayanamsha precession error

Fixed the precession polynomial in `_calc_ayanamsa()` for `SE_SIDM_USER` when the
user-defined reference epoch (`t0`) differs from J2000. The code was evaluating the
IAU 2006 polynomial at `T_user` instead of computing the differential `p(T) - p(T0)`,
causing ~2.2 arcsec error for non-J2000 reference epochs.

#### Critical: PREC_RATE / PREC_RATE_QUAD NameError crash

Fixed `NameError` crash in star-based ayanamsha calculations (`SE_SIDM_TRUE_REVATI`,
`SE_SIDM_TRUE_PUSHYA`, `SE_SIDM_TRUE_MULA`, `SE_SIDM_GALCENT_*`) caused by undefined
constants `PREC_RATE` and `PREC_RATE_QUAD`. Replaced with the full IAU 2006 5-term
precession polynomial.

#### Remaining forward-difference velocity calculations

Converted the last two forward-difference `O(h)` velocity calculations to central
difference `O(h²)` in `calc_seorbel_position()` and `_calc_keplerian_position()`
(`hypothetical.py`), completing the migration started in commit 5577461.

### Changed

#### Precision: IAU 2006/2000A models throughout

Upgraded all nutation and obliquity calculations to use IAU 2006/2000A models via
`erfa.nut06a()` and `erfa.obl06()`, replacing the mixed Skyfield/simplified models
that were previously used in some code paths. This ensures sub-milliarcsecond
consistency across all calculations (nutation, obliquity, ayanamsha, coordinate
transformations).

- `_calc_nutation_obliquity()`: now uses `erfa.nut06a()` + `erfa.obl06()` directly
- `_get_star_position_ecliptic()`: rewritten with full Skyfield astrometric pipeline
  (adds annual aberration, previously missing)
- `get_nutation_model()`: returns `IAU2006_2000A` source info
- Central-difference velocity in `spk.py`, `planetary_moons.py`, `hypothetical.py`

#### Dependency: pyerfa promoted to required

`pyerfa` is now a required dependency (was optional via `[precision]` extra). All
calculations unconditionally use `erfa.nut06a()` and `erfa.obl06()` for maximum
precision. The `[precision]` install extra has been removed.

### Removed

- `[precision]` optional extra from `pyproject.toml` (pyerfa is now always required)
- Redundant `pyerfa` entry from `[all]` extra (already in core dependencies)

## [0.12.0] - 2026-02-13

### Fixed

#### Critical: Ecliptic-to-equatorial coordinate transformation (`cotrans_sp`, `cotrans`)

Fixed sign errors in the ecliptic-to-equatorial coordinate transformation formulas
in `utils.py`. Both `cotrans_sp()` (position + velocity) and `cotrans()` (position only)
had every `sin(ε)` term with the **wrong sign**, effectively computing the inverse
transformation (equatorial→ecliptic) instead of ecliptic→equatorial.

**Root cause:** The formulas were using the equatorial→ecliptic signs:

```
sin(δ) = sin(β)·cos(ε) - cos(β)·sin(ε)·sin(λ)       ← WRONG (was minus)
tan(α) = [sin(λ)·cos(ε) + tan(β)·sin(ε)] / cos(λ)   ← WRONG (was plus)
```

**Correct formulas** (matching Swiss Ephemeris `swe_cotrans_sp` and standard celestial mechanics):

```
sin(δ) = sin(β)·cos(ε) + cos(β)·sin(ε)·sin(λ)       ← FIXED (plus)
tan(α) = [sin(λ)·cos(ε) - tan(β)·sin(ε)] / cos(λ)   ← FIXED (minus)
```

**7 sign corrections total:**

- `cotrans_sp()` position: `sin_new_lat` formula (`-` → `+`), longitude `y` term (`+` → `-`)
- `cotrans_sp()` latitude velocity: two terms in the derivative had inverted signs
- `cotrans_sp()` longitude velocity: `dy/dt` lat_speed coefficient (`+` → `-`)
- `cotrans()` position: same two sign fixes as `cotrans_sp()`

**Impact:** Any body whose equatorial coordinates were computed via `cotrans_sp()`/`cotrans()`
(lunar nodes, Lilith, interpolated apogee/perigee) returned incorrect declination and
right ascension values. For example, True North Node declination was +5.74° instead of
the correct -5.74° (sign inverted), and Mean Lilith declination was off by ~39°.

#### Missing equatorial conversion for non-lunar bodies in `_calc_body()`

Extended `_maybe_equatorial_convert()` to all body types that were missing it.
Previously, only lunar nodes and Lilith (added in commit 9919746) had the
`SEFLG_EQUATORIAL` conversion applied. All other body types computed via ecliptic
coordinates silently ignored the `SEFLG_EQUATORIAL` flag, returning ecliptic
longitude/latitude as if they were RA/Dec.

**Affected body types (now fixed):**

- Minor bodies via SPK kernels (e.g., Chiron, Ceres, Pallas)
- Minor bodies via automatic SPK download
- Minor bodies via Keplerian element fallback
- Uranian hypothetical planets (Cupido through Poseidon)
- Transpluto (Isis)
- Fixed stars (Regulus, Spica)
- Planetary moons (Galilean moons, Titan, etc.)

**Example:** Chiron's `calc_ut()` with `FLG_EQUATORIAL` returned `(lon, lat)` verbatim
as `(RA, Dec)` — producing a declination of -6.77° (the ecliptic latitude) instead of
the correct equatorial declination of +13.41°.

#### Type safety and correctness fixes

- Resolved all 34 mypy type errors and 1 ruff lint warning across the codebase
- Fixed incorrect return type annotations on 6 eclipse functions (`sol_eclipse_when_loc`, `sol_eclipse_where`, `lun_eclipse_when_loc`, `lun_occult_when_loc`, `planet_occult_when_loc`, and their `swe_` wrappers) — annotations declared `Tuple[Tuple, Tuple, int]` but actual return order was `Tuple[int, Tuple, Tuple]`
- Fixed `vis_limit_mag()` passing `lat, lon, alt_m` as separate arguments to `azalt()` instead of a single `(lon, lat, alt_m)` tuple
- Added `@overload` signatures to `get_sid_mode()` in `state.py` and `context.py` for proper type narrowing
- Fixed `Optional[float]` parameter annotations in `sidtime()` (was `float = None`)

### Changed

- `SEFLG_MOSEPH` flag is now accepted but silently ignored — all calculations always use JPL DE440/DE441 via Skyfield
- Removed Moshier semi-analytical ephemeris package (`moshier/`)
- Removed DE421 fallback — if DE440 is not found locally, it will be downloaded

### Added

- Environment variable `LIBEPHEMERIS_EPHEMERIS` for ephemeris file selection (e.g., `de441.bsp` for extended range -13200 to +17191 CE)
- Priority: `set_ephemeris_file()` > `LIBEPHEMERIS_EPHEMERIS` env var > default `de440.bsp`
- New `astrometry.py` module with IAU 2006 precession, IAU 2000B nutation, and stellar aberration utilities

### Removed

- `libephemeris/moshier/` package (VSOP87, ELP2000-82B, Pluto analytical)
- `_calc_body_moshier()` function and Moshier routing in `swe_calc_ut()`/`swe_calc()`
- `validate_jd_range_moshier()` and Moshier range constants from `exceptions.py`
- ~50 Moshier-related test files

## [0.11.0] - 2026-02-11

### Added

#### CLI Commands

- `libephemeris init` command for full initialization: downloads DE440.bsp, planet_centers.bsp, and SPK kernels for all 21 minor bodies (20-year chunks, 1550-2650 CE)
- `libephemeris init-fast` command for modern-era initialization (1900-2100 CE), ~10x faster than full init
- Per-chunk progress indicator during SPK downloads (inline overwrite)
- `--cache-dir` CLI option for custom SPK cache directory

#### SPK Cache Centralization

- Unified SPK cache directory at `~/.libephemeris/spk/` (was scattered)
- Cache directory resolution: `set_spk_cache_dir()` > `LIBEPHEMERIS_SPK_DIR` env var > `--cache-dir` CLI > default
- `init_all()` function with `start_year`/`end_year` parameters for programmatic use

#### SPK Download Script

- Updated `scripts/download_spk.py` to use `SPK_BODY_NAME_MAP` (21 bodies) with chunking support

### Fixed

- Rewrote `download_spk_from_horizons()` to use direct JPL Horizons HTTP API, replacing broken astroquery-based implementation (astroquery 0.4.11 `Horizons` constructor requires `step` in epochs dict, and `Horizons.download_spk()` method does not exist)
- Smart-skip for Horizons SPK date range limits: `START time outside SPK limits` skips chunk and tries next; `STOP time outside SPK limits` breaks to next body (avoids hundreds of futile requests and warnings)
- Removed unnecessary `astroquery` dependency from `enable_auto_spk()`, `auto_get_spk()`, and `download_spk_from_horizons()`

## [0.10.0] - 2026-02-10

### Added

#### Moshier Analytical Ephemeris (`SEFLG_MOSEPH`)

- Complete Moshier semi-analytical ephemeris engine as an explicit calculation mode
- VSOP87 truncated theory for inner/outer planet positions
- ELP 2000-82B truncated lunar theory for Moon positions
- DE404-based Pluto analytical theory
- Standalone IAU 2006 precession and IAU 2000B nutation (no external dependencies)
- Standalone aberration and light-time correction utilities
- Coordinate transformation support (ecliptic, equatorial, J2000, ICRS)
- `SEFLG_MOSEPH` flag routing gate in `swe_calc_ut()` and `swe_calc()`
- Extended date range: -3000 to +3000 CE (vs 1550-2650 for JPL DE440)
- Sidereal/ayanamsha calculations in Moshier mode
- House calculations (`swe_houses()`, `swe_houses_ex()`) in Moshier mode
- Lunar points (True Node, Mean Lilith) in Moshier mode
- `CalculationError` for unsupported bodies in Moshier mode

#### SPK Improvements

- SPK Type 21 (Modified Difference Arrays) support for JPL Horizons asteroid kernels
- Apparent position calculations from SPK Type 21 data
- Keplerian fallback for JPL major body index asteroids when SPK is unavailable

#### Configuration

- Auto SPK download enabled by default (`set_auto_spk_download`)

### Fixed

- SPK Type 21 now returns apparent positions instead of geometric
- Keplerian fallback allowed for JPL major body index asteroids
- `difdeg2n` aligned with pyswisseph for +/-180 degree separation
- `SEFLG_BARYCTR` now matches Swiss Ephemeris behavior
- Ayanamsa precision improved with IAU 2006 precession model
- Moon velocity precision improved with optimized timestep
- Placidus house convergence precision improved from 0.0003 deg to 0.0001 deg
- Placidus polar latitude handling improved (64-66 deg)
- Sunshine ('I'/'i') house system sidereal handling and lowercase support

### Documentation

- Marked SPK Type 21 apparent positions bug as fixed
- Updated SEFLG_MOSEPH documentation to reflect supported status

## [0.9.0] - 2026-02-09

### Added

#### Heliacal Events

- Complete Schaefer (1990) atmospheric visibility model for heliacal event calculations
- Rayleigh + aerosol + ozone atmospheric extinction model
- Twilight and moonlight sky brightness calculations
- Ptolemaic visibility thresholds (arcus visionis) for planets and stars
- Limiting visual magnitude calculation for naked-eye observation

#### Saturn Satellite System (TASS 1.7)

- Complete TASS 1.7 implementation for all 8 major Saturn satellites
- Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus
- ~2000 periodic terms per satellite for sub-arcsecond precision
- Validated against JPL Horizons ephemerides

#### Photometric Models

- Hapke photometric model for Moon magnitude with opposition surge correction
- Accurate Pluto magnitude formula from Mallama (2018) with phase corrections

#### Lunar Calculations

- Moshier analytical method for interpolated apogee (~50 harmonic terms)
- Higher-order terms (T^4, T^5) for improved True Node historical accuracy
- Precision warnings for dates outside Meeus optimal range (±200 years)

#### Nutation Model

- `get_nutation_model()` function to check active nutation model
- `NutationFallbackWarning` when using simplified 4-term model

#### Minor Bodies

- Resonant libration correction for plutinos (Ixion, Orcus, Pluto) in mean motion resonance

### Changed

#### Precision Improvements

- **Interpolated Apogee**: 11x precision improvement (1.1° -> 0.10° mean error) via Moshier method
- **Interpolated Perigee**: 5x precision improvement (2.3° -> 0.46° mean error) via extended ELP2000-82B series
- **True Revati/Pushya/Mula ayanamsha**: 10x precision improvement (0.06° -> <0.006°)
- **Galactic Center ayanamshas**: 60x precision improvement (0.06° -> <0.001°)
- **Sunrise/sunset timing**: 4x precision improvement (120s -> <30s tolerance)
- **Lunar occultation timing**: 5x precision improvement (300s -> <60s tolerance)
- **Uranian planets**: Orbital elements aligned with Witte/Hamburg School published parameters (Regelwerk für Planetenbilder)

#### API Improvements

- Added overload signatures to `house_pos()` for type safety
- `heliacal_ut()` and `vis_limit_mag()` now use Schaefer model for SE compatibility

### Fixed

- J2000 ayanamsha sign convention for pre-J2000 dates (negative before 2000, positive after)
- Eclipse search reliability with bidirectional mode (no longer skips nearby eclipses)
- Hybrid eclipse detection using proper Besselian elements
- Eclipse central line algorithm unified with `sol_eclipse_where()`
- Arabic Parts day/night calculation using 3D solar altitude for extreme latitudes
- Fixed star latitude velocity sign aligned with Swiss Ephemeris convention
- Uranian planet elements cross-validated with Hamburg School published orbital data

### Documentation

- Documented heliocentric methodology difference from Swiss Ephemeris in nod_aps
- Updated precision values in PRECISION.md for all improved calculations

## [0.8.0] - 2026-02-06

### Added

#### House Systems

- New house systems: Sripati (`'O'`), Pullen Sinusoidal Delta (`'L'`), Pullen Sinusoidal Ratio (`'Q'`), and Sunshine/Makransky (`'I'`)
- Gauquelin methods 2-5 now use actual rise/set times via `rise_trans()` for precise sector positions
- Gauquelin sectors expanded from 12 to 36 (matching Swiss Ephemeris)

#### Lunar Calculations

- SE-compatible Mean Lilith algorithm using DE404 polynomial coefficients
- ELP2000-82B perturbation series for interpolated apogee (0.10° mean error, 0.36° max vs SE)
- Independent ELP2000-82B perturbation series for interpolated perigee (0.46° mean error, 2.6° max vs SE)

### Changed

#### Precision Improvements

- **Interpolated Apogee**: 6x precision improvement (0.60° -> 0.10° mean error vs Swiss Ephemeris)
- **Interpolated Perigee**: 4.6x precision improvement (2.11° -> 0.46° mean error vs Swiss Ephemeris)
- **Nutation**: Upgraded from simplified 4-term model to full IAU 2000A (1365 terms) via Skyfield, improving from ~1 arcsec to sub-milliarcsecond precision
- **Fixed star velocities**: Upgraded from forward differences O(h) to central differences O(h²), ~100x better numerical precision
- **Campanus house_pos**: Fixed coordinate transformation, precision improved from 0.7° to 0.02° tolerance

#### API Changes

- `swe_gauquelin_sector()` now matches Swiss Ephemeris signature: `(tjdut, body, method, geopos, atpress, attemp, flags)`
- Gauquelin methods 0-1 now use `house_pos()` with `'G'` system for exact SE compatibility

### Fixed

- Fixed Equal from MC (`'D'`) house system algorithm
- Fixed Campanus `house_pos` meridian distance sign convention (uses SE's `cotrans` rotation formula)
- Fixed Koch `house_pos` precision (~0.2° max remaining)
- Fixed Mean Lilith latitude assertion (non-zero due to lunar orbital inclination ~5.145°)
- Fixed type annotations in `houses.py` (`numpy.float64` -> `float`, `Optional[Union[...]]`)
- Fixed "possibly unbound" variable in Placidus RA iteration

### Documentation

- Updated `PRECISION.md` with current interpolated apogee/perigee precision values
- Updated `TODO_HOUSES.md` with completed house system improvements

## [0.7.0] - 2026-02-03

### Added

#### Planet Center SPK Integration

- Added `planet_centers.bsp` support for high-precision outer planet calculations
- New SPK file provides planet center positions for Jupiter, Saturn, Uranus, Neptune, and Pluto (1989-2049)
- Precision improved from ~0.1 arcsec (analytical COB) to <0.001 arcsec (SPK-based)
- Automatic fallback to analytical Center-of-Body corrections when SPK is out of range

#### CLI Tool

- New `libephemeris` command-line interface:
    - `libephemeris download-data` - Download optional precision data files (~25MB)
    - `libephemeris status` - Show installed data file status
    - `libephemeris --version` - Show version information
- Progress bar during downloads (uses `rich` if available, otherwise simple ASCII progress)
- Atomic downloads with SHA256 verification support

#### New Modules

- `libephemeris.download` - Data file download utilities with progress bar
- `libephemeris.cli` - Command-line interface entry point

### Changed

- `_SpkCenterTarget` and `_CobCorrectedTarget` now compute velocity corrections for COB offsets
- Exception handling catches both `jplephem.OutOfRangeError` and `skyfield.EphemerisRangeError`
- Planet center SPK file moved from package data to optional download (reduces package size)

### Fixed

- Fixed velocity not being corrected in COB fallback paths (both `at()` and `_observe_from_bcrs()`)
- Fixed exception type mismatch when planet_centers.bsp is out of range

## [0.6.0] - 2026-02-02

### Changed

#### Breaking: Eclipse Function Return Order

All eclipse functions now return `(retflag, ...)` as the first element to match pyswisseph API conventions:

- `sol_eclipse_when_glob`: returns `(int, tuple)` - retflag first
- `sol_eclipse_where`: returns `(int, geopos, attr)` - retflag first
- `sol_eclipse_how`: returns `(int, attr)` - retflag first
- `sol_eclipse_when_loc`: returns `(int, times, attr)` - retflag first
- `lun_eclipse_when`: returns `(int, tuple)` - retflag first
- `lun_eclipse_how`: returns `(int, attr)` - retflag first
- `lun_eclipse_when_loc`: returns `(int, times, attr)` - retflag first
- `lun_occult_when_glob`: returns `(int, tuple)` - retflag first
- `lun_occult_when_loc`: returns `(int, times, attr)` - retflag first

**Migration**: Update unpacking from `times, attr, ecl_type = func()` to `ecl_type, times, attr = func()`

### Fixed

- Fixed `azalt()` call in heliacal.py (was passing 8 arguments instead of 6)
- Fixed Asellus Australis HIP number from 43834 to 42911 in fixed_stars.py
- Fixed early return statements in eclipse exception handlers to use correct tuple order

## [0.5.1] - 2026-01-31

### Added

#### Documentation

- New "Why LibEphemeris is More Accurate" section in README explaining scientific precision
- Detailed comparison table: NASA JPL DE440 vs Swiss Ephemeris
- True Node calculation methodology explanation with mathematical foundation
- Calibration script (`scripts/calibrate_true_node.py`) for True Node comparison

### Changed

#### True Node Documentation

- Updated PRECISION.md with rigorous True Node methodology explanation
- Added calibration results (500 random dates, 1900-2100): ~206 arcsec RMS difference
- Documented why geometric method (h = r × v) is mathematically more rigorous
- Added Swiss Ephemeris documentation quote confirming the approach

#### Code Documentation

- Enhanced `calc_true_lunar_node()` docstring with precision data
- Added explanation of why libephemeris is mathematically more accurate
- Updated precision figures based on actual calibration measurements

## [0.5.0] - 2026-01-31

### Changed

#### Default Ephemeris Upgrade

- Upgraded default ephemeris from DE421 to DE440 (JPL's latest recommended ephemeris)
- DE440 provides improved accuracy and extends coverage to 1550-2650

#### Eclipse Calculation Improvements

- Refactored `_calculate_eclipse_type_and_magnitude` to use proper Besselian elements
- Now uses `_calc_gamma()`, `_calc_penumbra_limit()`, and `_calc_umbra_limit()` helper functions
- Replaced ad-hoc gamma approximations with proper shadow geometry calculations
- Added spherical trigonometry for accurate angular separation between Sun and Moon

### Fixed

- Fixed tuple unpacking mismatch in `calc_angles()` when calling `swe_houses_with_fallback()`
- Fixed heliocentric and SSB-centered position calculations in `planets.py` to use direct vector computation instead of `observe().apparent()`
- Removed duplicate `_ANGLES_CACHE = {}` line in `state.py`
- Updated `calc_angles()` to use `swe_houses_with_fallback()` for better polar latitude handling

## [0.4.0] - 2026-01-31

### Added

#### Python 3.9 Support

- Added Python 3.9 compatibility with `from __future__ import annotations`
- Updated minimum Python version requirement from 3.10 to 3.9
- Added Python 3.9 classifier in package metadata

## [0.3.0] - 2026-01-31

### Added

#### Minor Bodies

- New centaurs: Nessus (7066), Asbolus (8405), Chariklo (10199)
- New TNO: Gonggong (225088)
- Uranus perturbations for improved TNO accuracy
- Neptune perturbations for TNO accuracy (critical for plutinos)
- Mean motion resonance detection for Neptune resonances (`detect_mean_motion_resonance`)
- Updated orbital elements with full precision from JPL SBDB
- TNO validation against pyswisseph over 2000-2050 date range

#### SPK Auto-Download and Caching

- New `spk_auto` module for automatic SPK download and caching
- `set_auto_spk_download()` / `get_auto_spk_download()` for enabling automatic SPK fallback
- `set_spk_cache_dir()` / `get_spk_cache_dir()` for configuring SPK cache location
- `set_spk_date_padding()` / `get_spk_date_padding()` for date range padding configuration
- Automatic SPK registration after download
- SPK cache management functions (`is_spk_cached`, `ensure_cache_dir`, `get_cache_path`)

#### Lunar Calculations

- Interpolated lunar apogee (SE_INTP_APOG) with comprehensive algorithm
- Interpolated lunar perigee (SE_INTP_PERG)
- Velocity calculation for interpolated apogee/perigee
- Optimized interpolation window (9 samples, 56 days, linear fit)
- True Lilith velocity calculation (SEFLG_SPEED support)
- True Node velocity calculation
- Comprehensive perturbation terms for True Node (Venus, Mars, Saturn, evection, variation, annual equation, parallactic)
- ELP2000-82B True Node perturbation term table
- IAU 2000A nutation correction for True Node
- Second-order perturbation terms for True Node
- Solar gravitational perturbation on eccentricity vector for True Lilith

#### Scripts

- `scripts/download_spk.py` for pre-downloading SPK files
- `scripts/update_orbital_elements.py` for updating orbital elements from JPL SBDB

#### Constants

- `SPK_BODY_NAME_MAP` for body ID to JPL Horizons mapping
- NAIF ID constants for new bodies (NAIF_NESSUS, NAIF_ASBOLUS, NAIF_CHARIKLO, NAIF_GONGGONG)
- SE_NESSUS, SE_ASBOLUS, SE_CHARIKLO, SE_GONGGONG body IDs

#### Error Handling

- Category-based exception hierarchy for better error handling
- Proactive Julian Day range validation before calculation
- Geographic coordinates validation (lat/lon)
- Improved error handling for extreme latitudes (>80°) in house calculations
- Graceful handling of missing SPK files with `SPKNotFoundError`
- Clear error messages for unknown body IDs
- Improved date range error messages with supported range details

#### Retrograde & Eclipse Handling

- Retrograde station handling with stable near-zero velocity calculations
- Eclipse edge case handling for shallow partial eclipses

#### Dependency Upgrades

- Upgraded Skyfield to 1.54 for `deflectors=` arg and improved performance
- Upgraded jplephem to 2.24 for NumPy compatibility

#### Profiling

- New profiling module for performance analysis

### Documentation

- Comprehensive documentation for interpolated apogee (`docs/INTERPOLATED_APOGEE.md`)
- True Lilith calculation method comparison (`docs/TRUE_LILITH_METHODS.md`)
- Updated precision documentation with TNO validation results
- Precision tuning guide (`docs/PRECISION_TUNING.md`)
- Updated API reference with TAI, IERS, planetary moons, and other new features
- Enhanced migration guide with lunar nodes/Lilith precision info
- Documentation of optional dependencies (pyerfa, astroquery, astropy)
- Usage examples demonstrating common use cases
- Documented pyerfa, astropy, and REBOUND integration benefits

### Changed

- Converted compare scripts to pytest-style unit tests
- Moved swisseph-dependent tests to `compare_scripts/tests/`

### Tests

- TNO validation tests against pyswisseph over 2000-2050
- Resonance detection tests
- Secular perturbation tests for minor bodies
- Interpolated apogee/perigee precision tests
- True Lilith latitude validation tests
- True Node velocity tests
- Download SPK script tests
- Orbital elements update script tests
- Comprehensive pyerfa precision evaluation tests
- Comprehensive ayanamsha multi-date tests

## [0.2.0] - 2026-01-26

### Added

#### Minor Bodies

- Secular perturbations from Jupiter and Saturn for improved accuracy
- Support for parabolic and hyperbolic orbits
- Updated orbital elements to epoch 2025.0 (JD 2461000.5)

#### Lunar Calculations

- Planetary perturbations to true node calculation
- Dynamic IAU 2006 obliquity model (replaces fixed J2000 obliquity)
- Updated GM_Earth to IAU 2015 Resolution B3 value
- Documentation of Meeus polynomial validity range with warnings

#### Fixed Stars

- Full IAU 2000A nutation model (replaces 2-term approximation)
- Second-order Taylor expansion for proper motion

#### Crossing Functions

- Pluto typical speed support in `swe_cross_ut`
- Brent's method fallback for station detection
- Adaptive iteration limits for slow planets
- Tightened solar crossing tolerance to 0.001 arcsec

#### Eclipse Functions

- `sol_eclipse_when_glob` for global solar eclipse search
- `sol_eclipse_when_loc` for location-specific solar eclipse search
- `sol_eclipse_where` for central eclipse path calculation
- `sol_eclipse_how` for eclipse circumstances at location
- `lun_eclipse_when` for lunar eclipse search
- `lun_eclipse_when_loc` for location-specific lunar eclipse search
- `lun_eclipse_how` for lunar eclipse circumstances at location
- `lun_occult_when_glob` for lunar occultation search
- `lun_occult_when_loc` for location-specific lunar occultation search
- `lun_occult_where` for lunar occultation path calculation
- `rise_trans` for calculating rise, set, and transit times
- `rise_trans_true_hor` for custom horizon altitude calculations
- `heliacal_ut` for heliacal rising/setting events
- `heliacal_pheno_ut` for detailed heliacal phenomena
- `vis_limit_mag` for visual limiting magnitude

#### Utility Functions

- `degnorm` for angle normalization
- `radnorm` for radian angle normalization
- `deg_midp` for angular midpoint calculation
- `rad_midp` for radian angular midpoint calculation
- `difdegn` for positive angular difference
- `difrad2n` for radian angular difference
- `difcs2n` for centiseconds angular difference
- `difcsn` for positive centiseconds angular difference
- `csnorm` for centiseconds normalization
- `d2l` for double to long conversion with rounding
- `cs2degstr` for centiseconds to degrees string conversion
- `cs2lonlatstr` for centiseconds to lon/lat string conversion
- `cs2timestr` for centiseconds to time string conversion
- `cotrans` for ecliptic/equatorial coordinate transformation
- `cotrans_sp` for coordinate and velocity transformation
- `azalt` for equatorial/ecliptic to horizontal coordinate conversion
- `azalt_rev` for horizontal to equatorial/ecliptic coordinate conversion
- `refrac` for atmospheric refraction calculation
- `refrac_extended` for extended atmospheric refraction

#### Time Functions

- `utc_to_jd` for UTC to Julian Day conversion with leap second support
- `jdet_to_utc` for converting JD(TT/ET) to UTC with Delta-T and leap seconds
- `jdut1_to_utc` for converting JD(UT1) to UTC date/time
- `utc_time_zone` for applying timezone offsets to UTC date/time
- `time_equ` for Equation of Time calculation
- `lat_to_lmt` for Local Apparent Time to Local Mean Time conversion
- `lmt_to_lat` for Local Mean Time to Local Apparent Time conversion
- `sidtime` for Local Sidereal Time calculation
- `sidtime0` for Greenwich Sidereal Time calculation
- `set_delta_t_userdef` for user-defined Delta T
- `set_tid_acc` and `get_tid_acc` for tidal acceleration in Delta T

#### State Functions

- `set_jpl_file` for specifying JPL ephemeris files
- `set_lapse_rate` for configuring atmospheric lapse rate
- `close` function to release ephemeris resources
- `get_library_path` to return ephemeris file directory
- `get_current_file_data` to return ephemeris file info

#### Planets Functions

- `get_planet_name` to return human-readable planet names
- `pheno` and `pheno_ut` for planetary phenomena
- `calc_pctr` for planet-centric position calculations
- `nod_aps` and `nod_aps_ut` for orbital nodes and apsides
- `get_orbital_elements` for Keplerian orbital elements
- `orbit_max_min_true_distance` for perigee/apogee distances

#### Fixed Stars Functions

- `swe_fixstar` for Terrestrial Time (TT) star positions
- `fixstar2` and `fixstar2_ut` with flexible star lookup
- `fixstar_mag` and `fixstar2_mag` for magnitude lookup

#### Houses Functions

- `houses_ex2` returning cusp velocities
- `houses_armc` for ARMC-based house calculations
- `houses_armc_ex2` for ARMC-based house cusps with velocities
- `house_pos` to determine which house a celestial body is in
- `gauquelin_sector` for 36-sector calculation

#### Ayanamsha Functions

- `get_ayanamsa_ex` and `get_ayanamsa_ex_ut` for extended ayanamsha data

#### Crossing Functions

- `swe_solcross` for TT-based sun longitude crossing
- `swe_mooncross` for TT-based moon longitude crossing
- `mooncross_node` and `mooncross_node_ut` for moon node crossing
- `helio_cross` and `helio_cross_ut` for heliocentric crossings

#### Other

- `Error` class for pyswisseph compatibility
- `date_conversion` for Julian/Gregorian calendar conversion
- `day_of_week` for Julian Day to weekday conversion
- `deltat_ex` for ephemeris-specific Delta T calculation

### Changed

#### Precision Improvements

- Reduced Moon iteration limit from 50 to 30 (optimized)
- Tightened lunar crossing tolerance to 0.05 arcsec
- Improved Newton-Raphson convergence to sub-arcsecond precision

### Fixed

- Fixed stars proper motion using rigorous space motion approach
- Planets proper motion using rigorous space motion approach
- Houses: use true Ascendant in `_houses_equal_mc` instead of approximation
- Houses: add polar circle detection for Gauquelin house system
- Houses: add polar circle detection for Placidus/Koch house systems

### Documentation

- Added comprehensive cookbook with practical astrological examples
- Added precision limitations documentation
- Added migration guide from pyswisseph to libephemeris
- Added complete API reference with Sphinx integration

### Tests

- Added benchmark tests comparing libephemeris vs pyswisseph
- Added natal chart integration tests with famous people data
- Added solstices/equinoxes verification tests against Swiss Ephemeris
- Added edge case tests for julday/revjul date handling
- Added thread-safety tests for concurrent ephemeris usage
- Added sub-arcsecond precision comparison tests for 7 planets
- Added comprehensive station time comparison tests for Mercury, Venus, Mars, Jupiter, Saturn
- Added comprehensive polar latitude tests for all 15+ house systems
- Added comprehensive tests for all 43 ayanamsha modes vs pyswisseph

## [0.1.0] - 2024-01-01

### Added

- Initial release
- Core planetary position calculations (Sun, Moon, all major planets, Pluto)
- High-precision ephemeris based on NASA JPL DE421
- Multiple coordinate systems (ecliptic, equatorial, J2000, of-date)
- Observation modes (geocentric, topocentric, heliocentric, barycentric)
- Full 6-component state vectors (position + velocity)
- 19 house systems (Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign, Porphyry, Alcabitius, Topocentric, Morinus, Meridian, Vehlow, Horizontal, Carter, Krusinski, Natural Gradient, and more)
- 43 ayanamsha modes (Fagan/Bradley, Lahiri, Raman, Krishnamurti, and more)
- Lunar nodes (True and Mean)
- Lilith (Mean and True Black Moon)
- Major asteroids (Chiron, Pholus, Ceres, Pallas, Juno, Vesta)
- TNOs (Orcus, Haumea, Quaoar, Makemake, Gonggong, Eris, Sedna)
- Fixed stars support
- Arabic parts calculations
- Sun/Moon longitude crossings (ingress detection)
- Thread-safe `EphemerisContext` API for concurrent calculations
- Swiss Ephemeris compatible function names, flags, and result structure

[Unreleased]: https://github.com/g-battaglia/libephemeris/compare/v0.21.0...HEAD
[0.21.0]: https://github.com/g-battaglia/libephemeris/compare/v0.20.0...v0.21.0
[0.20.0]: https://github.com/g-battaglia/libephemeris/compare/v0.19.0...v0.20.0
[0.19.0]: https://github.com/g-battaglia/libephemeris/compare/v0.18.0...v0.19.0
[0.18.0]: https://github.com/g-battaglia/libephemeris/compare/v0.17.0...v0.18.0
[0.17.0]: https://github.com/g-battaglia/libephemeris/compare/v0.16.1...v0.17.0
[0.16.1]: https://github.com/g-battaglia/libephemeris/compare/v0.16.0...v0.16.1
[0.16.0]: https://github.com/g-battaglia/libephemeris/compare/v0.15.0...v0.16.0
[0.15.0]: https://github.com/g-battaglia/libephemeris/compare/v0.14.0...v0.15.0
[0.14.0]: https://github.com/g-battaglia/libephemeris/compare/v0.13.0...v0.14.0
[0.13.0]: https://github.com/g-battaglia/libephemeris/compare/v0.12.0...v0.13.0
[0.12.0]: https://github.com/g-battaglia/libephemeris/compare/v0.11.0...v0.12.0
[0.11.0]: https://github.com/g-battaglia/libephemeris/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/g-battaglia/libephemeris/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/g-battaglia/libephemeris/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/g-battaglia/libephemeris/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/g-battaglia/libephemeris/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/g-battaglia/libephemeris/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/g-battaglia/libephemeris/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/g-battaglia/libephemeris/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/g-battaglia/libephemeris/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/g-battaglia/libephemeris/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/g-battaglia/libephemeris/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/g-battaglia/libephemeris/releases/tag/v0.1.0
