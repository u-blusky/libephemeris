# Keplerian Fallback — Possible Improvements

This document catalogs every known improvement opportunity for the Keplerian
fallback pipeline in `libephemeris/minor_bodies.py`. Items are organized by
estimated cost and expected precision gain.

The Keplerian fallback is the last resort in the minor body calculation chain:

```
SPK kernel → Auto-download SPK → ASSIST n-body → Keplerian (this code)
```

It matters when SPK is unavailable AND ASSIST is not installed or its data
files are missing. The goal is to minimize error growth over time when this
fallback is the only option available.

---

## Current precision (measured)

From `tests/test_keplerian_precision_benchmark.py`, comparing Keplerian
positions against SPK truth for 5 bodies (Ceres, Pallas, Juno, Vesta, Chiron):

| Time from epoch | Max error | Dominant error source |
|-----------------|-----------|------------------------|
| At epoch        | 0.002"    | Floating point only |
| 1 month         | 7.5"     | Short-period perturbations |
| 6 months        | 49"      | Short-period perturbations |
| 1 year          | 1.9'     | Secular drift begins |
| 5 years         | 27'      | Secular + short-period |
| 10 years        | 45'      | Secular dominates |
| 25 years        | 2.6'     | Multi-epoch table helps |
| 50 years        | 3.6°     | Secular overwhelms |
| 100 years       | 3.5°     | Multi-epoch + secular |

---

## Current implementation summary

**37 bodies** with osculating elements at epoch JD 2461000.5 (2025-Sep-19):
13 main belt, 6 centaurs, 9 TNOs, 9 NEAs.

**Multi-epoch elements** for **6 bodies** only (Ceres, Pallas, Juno, Vesta,
Chiron, Pholus): 17 epochs each at 50-year intervals, 1650–2450 CE.

**Secular perturbations** from 4 planets (Jupiter, Saturn, Uranus, Neptune):
- omega (arg. perihelion): linear precession
- Omega (ascending node): linear regression
- e (eccentricity): oscillation via (h,k) vector formalism (Laplace-Lagrange)
- i (inclination): oscillation via (p,q) vector formalism (Laplace-Lagrange)
- n (mean motion): **NOT perturbed** (d_n = 0 always)
- a (semi-major axis): **NOT perturbed** (constant)

**Libration model** for 2 plutinos (Ixion, Orcus).

**Kepler equation solver**: Newton-Raphson, tolerance 1e-8, 30 iterations
(elliptic), 50 iterations (hyperbolic), Barker closed-form (parabolic).

---

## Low-cost improvements

### L1. Denser multi-epoch table (every 20 years instead of 50)

**Cost:** Low — only data generation, no algorithm changes.
**Impact:** High at 25-50 year scales.

Currently the multi-epoch table has 17 entries per body at 50-year intervals
(1650–2450). Maximum propagation distance from nearest epoch: ~25 years.
Reducing to 20-year intervals gives ~10-year max propagation, reducing the
error at 25-50 year timescales by roughly 2-3x.

**Changes required:**
- Regenerate `MINOR_BODY_ELEMENTS_MULTI` in `minor_bodies.py` with 20-year
  spacing: ~41 entries per body instead of 17
- Data source: SPK type 21 state vectors → Keplerian conversion (same method
  as current, script exists in `scripts/`)
- No code changes to `_get_closest_epoch_elements()` — it already does linear
  scan

**Estimated data size increase:** 6 bodies × 41 entries × ~8 floats × ~80
chars = ~15 KB (from current ~6 KB). Negligible.

**Expected improvement:**
- 25 years: 2.6' → < 1' (max propagation drops from 25 → 10 years)
- 50 years: 3.6° → < 1° (max propagation drops from 25 → 10 years)

### L2. Multi-epoch for all 37 bodies (not just 6)

**Cost:** Low-medium — data generation for 31 additional bodies.
**Impact:** High for TNOs, centaurs, NEAs at long timescales.

Currently only 6 bodies have multi-epoch data. The remaining 31 bodies
propagate from a single 2025 epoch. For a body like Sedna (a=549 AU,
P=~12,000 yr), a 100-year Keplerian propagation from a single epoch
accumulates significant error.

**Changes required:**
- Generate SPK files for all 37 bodies covering 1650-2450 (or wider for
  bodies with available SPK coverage)
- Convert state vectors to Keplerian elements at 20-year intervals
- Add entries to `MINOR_BODY_ELEMENTS_MULTI`

**Complications:**
- Some NEAs (Apophis, Bennu, Ryugu) have chaotic orbits — close planetary
  encounters make Keplerian propagation unreliable regardless of epoch density
- Some TNOs (Sedna, Gonggong) may have limited SPK coverage from Horizons
- NEAs with extreme eccentricity (Icarus e=0.827) may have poorly conditioned
  Keplerian elements at some epochs

**Expected improvement:** Variable by body; most benefit for centaurs
(Nessus, Asbolus, Chariklo) and non-resonant TNOs (Quaoar, Varuna, Haumea).

### L3. Semi-major axis secular perturbation (d_n ≠ 0)

**Cost:** Low — straightforward formula addition.
**Impact:** Low-moderate over decades.

Currently `d_n` is hardcoded to 0.0 in `calc_secular_perturbation_rates()`
(line 772). The semi-major axis (and thus mean motion) does receive secular
perturbations from planetary interactions, particularly near mean-motion
resonances.

The second-order secular theory gives a correction to `n`:

```
d_a/dt = 0  (first-order secular theory)
```

At first order, `a` is constant. But the mean longitude rate `n` receives
a correction from the interaction with the forced eccentricity:

```
d_n = -(3/2) * (n/a) * d_a
```

For most asteroids this is negligible (< 0.01"/year), but for near-resonant
bodies (e.g., Hildas near 3:2, plutinos near 2:3) it could matter at the
arcminute level over decades.

**Changes required:**
- Compute `d_n` from second-order secular theory in
  `calc_secular_perturbation_rates()`
- Use `n_pert = elements.n + d_n * dt` instead of `n_pert = elements.n`

**Expected improvement:** < 1' over 50 years for most bodies. More for
near-resonant bodies.

### L4. Evolving planet elements

**Cost:** Low — replace constants with linear functions of time.
**Impact:** Low-moderate over centuries.

The 4 perturbing planets (Jupiter, Saturn, Uranus, Neptune) use static
J2000.0 orbital elements. Over centuries, the planets' own elements drift
(Jupiter's eccentricity oscillates with a ~300,000 year period, etc.).
Using linear rates for the planet elements would improve the forced
eccentricity/inclination vectors.

**Changes required:**
- Replace static constants (JUPITER_E, JUPITER_I, etc.) with functions
  of time: `e_J(t) = e_J0 + de_J * (t - t_J2000)`
- Rates available from Simon et al. (1994) or Standish (1992)
- Update `_calc_forced_elements()` to accept `jd_tt` and compute
  planet elements at that date

**Expected improvement:** ~10-30" over 500-year propagations. Negligible
for < 100 years.

### L5. Kepler solver convergence warning

**Cost:** Very low — add a warning log.
**Impact:** Debugging/correctness only.

The Newton-Raphson solver returns the last iterate without warning if 30
(elliptic) or 50 (hyperbolic) iterations are exceeded. For highly eccentric
orbits (Icarus e=0.827, Sedna e=0.861), non-convergence is possible.

**Changes required:**
- Add `logger.warning()` if max iterations reached
- Optionally: implement Markley's (1995) or Raposo-Pulido & Pelaez (2017)
  starter for high eccentricity, which guarantees convergence in 2-3
  iterations for any eccentricity

**Expected improvement:** No precision change, but prevents silent failures.

### L6. Laplace coefficient integration accuracy

**Cost:** Very low — increase step count.
**Impact:** Low (sub-arcsecond) for near-resonant bodies.

The trapezoidal integration for Laplace coefficients uses 100 steps
(`_calc_laplace_coefficients()`, line 698). For large alpha (plutinos:
alpha ≈ 0.76), higher resolution may improve the secular rates.

**Changes required:**
- Increase `n_steps` from 100 to 500 or 1000
- Or: use Gauss-Legendre quadrature for exact integration with fewer points
- Or: use the closed-form recursion formula for Laplace coefficients
  (Brouwer & Clemence 1961, eq. 15.21)

**Expected improvement:** < 0.1" for most bodies. May matter for bodies
very close to resonance where alpha → 1.

---

## Medium-cost improvements

### M1. Short-period perturbation terms (analytical)

**Cost:** Medium-high — requires implementing perturbation series.
**Impact:** High at 1 month to 5 year timescales (currently 7-49" error).

This is the dominant error source at timescales < 5 years. Short-period
perturbations arise from conjunctions with Jupiter and Saturn and oscillate
with the synodic period (~400 days for Jupiter-Ceres).

**Theory:** Brouwer & Clemence (1961), Chapter 15. The first-order
short-period terms for the disturbing function give corrections to all
6 orbital elements as trigonometric series in the mean anomalies of the
asteroid and the perturbing planet.

**For the longitude of a main-belt asteroid, the dominant terms are:**

```
Δλ_short ≈ Σ_j A_j(a, e, i, a_J, e_J) × sin(k₁ M + k₂ M_J + k₃ ω + k₄ Ω)
```

where the amplitudes A_j depend on Laplace coefficients and eccentricity
functions. For Ceres, the dominant term has amplitude ~300" with period
~466 days (synodic period with Jupiter).

**What was tried and failed (Task 4.2):**
- Empirical Fourier fit over 800 years — secular drift dominated residuals
- Polynomial detrending + Fourier fit — correct amplitudes (~300" Ceres,
  ~780" Pallas) but WORSEN positions at short timescales because the
  amplitudes vary with time as elements drift. A static Fourier series
  captures the average amplitude which is wrong at any specific date.

**What would work (not yet attempted):**
- Analytical perturbation theory (Brouwer & Clemence): the amplitudes are
  functions of the current orbital elements, so they naturally evolve.
  Requires implementing the disturbing function expansion to second order
  in eccentricity and first order in inclination.
- Alternative: use the VSOP-like approach with precalculated coefficient
  tables per body from JPL data (Chapront-Touzé 1988 style).

**Estimated effort:** 2-4 days of implementation. ~200-400 lines of code.
Coefficient tables for Jupiter and Saturn perturbations on each body.

**Expected improvement:**
- 1 month: 7" → ~1"
- 6 months: 49" → ~5"
- 1 year: 1.9' → ~20"

### M2. Libration model for additional resonant TNOs

**Cost:** Medium — requires identifying resonances and fitting parameters.
**Impact:** High for specific bodies over decades.

Currently only Ixion (2:3) and Orcus (2:3) have libration corrections.
Other TNOs in mean-motion resonances with Neptune:

| Body | Resonance | Libration? | Currently modeled? |
|------|-----------|------------|-------------------|
| Ixion | 2:3 | Yes, ~78° amplitude | Yes |
| Orcus | 2:3 | Yes, ~68° amplitude | Yes |
| Haumea | 7:12 | Possibly | No |
| Makemake | Not resonant | N/A | N/A |
| Gonggong | 3:10 | Yes, ~30° amplitude | No |
| Eris | Not resonant | N/A | N/A |
| Varuna | Near 3:4? | Uncertain | No |
| Quaoar | Near 5:8? | Possible | No |
| Sedna | Not resonant | N/A | N/A |

**Changes required:**
- Identify which bodies are in confirmed resonances (from literature or
  numerical integration)
- Fit libration parameters (amplitude, period, phase) from SPK data
- Add entries to `PLUTINO_LIBRATION_PARAMS`
- The libration correction framework already exists (`calc_libration_correction()`)

**Expected improvement:** For resonant bodies, reduces error from degrees
to arcminutes over 100-year propagations.

### M3. Benchmark expansion to all 37 bodies

**Cost:** Medium — needs SPK truth values for all bodies.
**Impact:** Better understanding of where improvements are needed.

Currently the benchmark tests only 5 bodies. Expanding to all 37 would
reveal which bodies have the worst Keplerian precision and where to
prioritize improvements.

**Changes required:**
- Ensure SPK files are available for all 37 bodies at test time
- Extend `BENCHMARK_BODIES` in `test_keplerian_precision_benchmark.py`
- Group results by category (main belt, centaur, TNO, NEA)
- Identify outliers where Keplerian is catastrophically wrong

**Expected outcome:** A precision matrix showing which bodies and timescales
are the weakest, guiding further improvements.

### M4. Mean motion resonance corrections (beyond libration)

**Cost:** Medium — analytical or semi-analytical.
**Impact:** Moderate for near-resonant bodies.

Bodies near (but not in) mean-motion resonances experience amplified
secular perturbation rates. The current secular model treats all bodies
equally, but near-resonant bodies (e.g., Hildas at 3:2 with Jupiter)
have much larger perturbation effects.

**Changes required:**
- Detect near-resonance condition: `|p/q - n_body/n_planet| < threshold`
- Apply resonant secular perturbation formulae (Murray & Dermott 1999,
  Chapter 8)
- Different perturbation model for bodies within the resonance width

**Bodies affected:** Ixion, Orcus (2:3 Neptune), potentially Gonggong
(3:10 Neptune), some main belt bodies near 3:1, 5:2, 7:3 Jupiter
resonances.

**Expected improvement:** Factor of 2-5x for near-resonant bodies at
10-50 year timescales.

---

## High-cost improvements

### H1. Second-order secular perturbation theory

**Cost:** High — substantial mathematical complexity.
**Impact:** Moderate improvement over first-order secular theory.

The current implementation uses first-order Laplace-Lagrange theory. The
second-order theory (Hori 1966, Yuasa 1973) includes:

- Cross-terms between different perturbing planets (Jupiter × Saturn
  coupling)
- Second-order eccentricity/inclination effects
- Long-period terms with periods of tens of thousands of years

These become important for:
- Bodies with high eccentricity (Hidalgo e=0.66, Icarus e=0.83)
- Bodies with high inclination (Pallas i=35°, Eris i=44°)
- Bodies crossing multiple planetary orbits (centaurs)

**Changes required:**
- Implement second-order averaging of the disturbing function
- Compute cross-coupling coefficients between planet pairs
- Add higher-order eccentricity functions d₁, d₂ to the Laplace
  coefficient calculation

**Estimated effort:** 3-5 days. ~500-800 lines of code.

**Expected improvement:** Factor of 2-3x at 100-year timescales for
high-e/high-i bodies. Negligible for low-e/low-i main belt asteroids.

### H2. Analytical short-period + long-period perturbation theory (full Brouwer)

**Cost:** High — complete Brouwer artificial satellite theory adapted for
heliocentric orbits.
**Impact:** High — would bring Keplerian to ~1" precision at 1-year scales.

Full implementation of Brouwer & Clemence (1961) perturbation theory for
heliocentric orbits:

1. **Short-period terms** (~synodic period): δa, δe, δi, δω, δΩ, δM
   as Fourier series in mean anomalies
2. **Long-period terms** (~secular, but with period ∝ 1/e): δω, δΩ
   corrections to secular rates
3. **Mixed secular-periodic terms**: amplitude modulation of short-period
   terms by secular evolution

This is what analytical orbit propagators (SGP4 for satellites, VSOP for
planets) implement. For asteroids, the key perturbers are Jupiter and Saturn.

**Changes required:**
- Disturbing function expansion in Legendre polynomials or Hansen
  coefficients
- Coefficient tables for each perturber-asteroid pair
- At evaluation time: compute all periodic terms and sum corrections

**Estimated effort:** 5-10 days. ~1000-2000 lines of code. Extensive
validation required.

**Expected improvement:**
- 1 month: 7" → < 0.5"
- 1 year: 1.9' → < 5"
- 10 years: 45' → < 1'

### H3. Semi-analytical propagation (Encke's method)

**Cost:** High — requires numerical integration at evaluation time.
**Impact:** Very high — ephemeris-quality for short arcs.

Instead of pure Keplerian + perturbation corrections, use Encke's method:
propagate the osculating orbit analytically (Kepler), then integrate the
perturbation residuals numerically with a low-order integrator (e.g., RK4
with large step size).

This combines the speed of Keplerian propagation (for the dominant 2-body
motion) with the accuracy of numerical integration (for the perturbative
residuals, which are small and smooth).

**Changes required:**
- Implement planetary position lookup (could use Skyfield or precomputed
  tables)
- RK4 or Störmer-Cowell integrator for the perturbation equations of motion
- Step size: ~10 days for main belt, ~30 days for TNOs
- Caching of intermediate results for repeated queries

**Estimated effort:** 3-5 days. But adds a runtime dependency on planetary
positions (Skyfield or LEB).

**Expected improvement:** Sub-arcsecond for propagation up to decades.
Essentially equivalent to ASSIST but without the ASSIST data files.

**Trade-off:** This blurs the line between Keplerian fallback and numerical
integration. If Skyfield is available (which it always is — it's a required
dependency), this approach would effectively be "lightweight ASSIST" without
the external data files.

### H4. Non-gravitational forces for NEAs

**Cost:** High — requires Yarkovsky/YORP modeling.
**Impact:** Critical for specific bodies (Apophis, Bennu, Ryugu) at
decade timescales.

Near-Earth asteroids experience measurable non-gravitational accelerations:

| Force | Magnitude (AU/day²) | Timescale | Bodies affected |
|-------|---------------------|-----------|-----------------|
| Yarkovsky | ~1e-15 | Decades | All NEAs |
| Solar radiation pressure | ~1e-13 | Days | Small NEAs (< 100m) |
| Outgassing | Variable | Perihelion | Cometary bodies |

For Apophis, the Yarkovsky effect shifts the orbit by ~200 km/year, which
translates to ~0.1" per year in geocentric longitude — small, but it
accumulates.

**Changes required:**
- Store Yarkovsky parameter (da/dt in AU/My) per body
- Apply as `a(t) = a₀ + (da/dt) × dt` in the Keplerian propagation
- Data source: JPL SBDB (https://ssd.jpl.nasa.gov/) provides measured
  da/dt for several NEAs

**Expected improvement:** ~1" per decade for NEAs with known Yarkovsky.
Irrelevant for main belt and TNOs.

### H5. Cubic Hermite interpolation of multi-epoch elements

**Cost:** Medium-high.
**Impact:** Moderate — smoother transitions between epochs.

Currently `_get_closest_epoch_elements()` picks the single closest epoch.
This creates discontinuities in the computed position when the "closest
epoch" switches from one entry to the next (at midpoints between epochs).

A cubic Hermite interpolation of the orbital elements between adjacent
epochs would:
- Eliminate position discontinuities at epoch boundaries
- Use the velocity information implicit in the element rates
- Provide C¹-continuous positions across the full date range

**Changes required:**
- Compute element rates (from finite differences of adjacent epochs)
- Implement Hermite interpolation for each element (a, e, i, ω, Ω, M)
  with proper angle unwrapping
- Handle eccentricity bounds (0 < e < 1) and inclination bounds (0 < i < π)

**Complication:** Orbital elements are not smooth functions — they can have
discontinuities near resonances or close encounters. Hermite interpolation
assumes smoothness and would produce artifacts in these cases.

**Expected improvement:** Eliminates ~10" discontinuities at epoch
boundaries. Smooth position function across full date range.

---

## Structural issues (not precision improvements)

### S1. Missing convergence handling in Kepler solver

The Newton-Raphson solver silently returns a non-converged result after
30 (elliptic) or 50 (hyperbolic) iterations. Should log a warning and
optionally fall back to a bisection method.

### S2. No validation of orbital element sanity

No check that e < 1 for elliptic orbits, a > 0, 0 < i < 180°, etc.
Garbage elements would produce garbage positions silently.

### S3. Parabolic/hyperbolic orbits skip all perturbations

`calc_minor_body_position()` only applies secular perturbations for
elliptic orbits (e < 1). Hyperbolic comets and parabolic bodies get zero
perturbation corrections. This is correct for single-pass trajectories
but wrong for bodies with poorly determined orbits that happen to have
e ≈ 1.

### S4. Co-orbital detection too aggressive

`_calc_forced_elements()` skips a perturber if `|a - a_planet| < 0.1 AU`.
This threshold is arbitrary and could skip important perturbations for
Trojan asteroids or bodies in horseshoe orbits with Jupiter.

---

## Priority recommendations

### If the goal is maximum impact for minimum effort:

1. **L2** — Multi-epoch for all 37 bodies (biggest gap today)
2. **L1** — Denser multi-epoch (20-year intervals)
3. **M3** — Benchmark all 37 bodies (understand the problem first)
4. **L5** — Kepler solver convergence warning (correctness)

### If the goal is pushing Keplerian to its theoretical limit:

1. **M1** — Analytical short-period perturbations (eliminates 7-49" errors)
2. **H2** — Full Brouwer theory (sub-arcsecond at 1 year)
3. **H5** — Hermite interpolation (smooth transitions)
4. **L4** — Evolving planet elements (century-scale accuracy)

### If the goal is making the Keplerian fallback irrelevant:

1. **H3** — Semi-analytical Encke propagation (sub-arcsecond without ASSIST)

This would effectively replace the Keplerian fallback with a lightweight
numerical integrator that uses Skyfield's planetary positions directly.
The advantage over ASSIST: no extra 1.3 GB download. The disadvantage:
slower than pure Keplerian (but still much faster than Skyfield per-body
evaluation).

---

## References

- Brouwer, D. & Clemence, G.M. (1961). Methods of Celestial Mechanics.
  Academic Press. — Chapters 11-16 (perturbation theory).
- Murray, C.D. & Dermott, S.F. (1999). Solar System Dynamics. Cambridge
  University Press. — Chapters 6-8 (secular and resonant dynamics).
- Hori, G. (1966). Theory of general perturbation with unspecified
  canonical variable. Publ. Astron. Soc. Japan, 18, 287.
- Markley, F.L. (1995). Kepler equation solver. Celestial Mechanics, 63, 101.
- Raposo-Pulido, V. & Pelaez, J. (2017). An efficient code to solve the
  Kepler equation. MNRAS, 467, 1702.
- Simon, J.L. et al. (1994). Numerical expressions for precession formulae
  and mean elements for the Moon and planets. A&A, 282, 663.
