# Implementation Plan: Documentation Reorganization

> Date: March 2026
>
> Goal: Reorganize all 30 markdown files into a clean, thematic structure under `docs/`
> with uniform style, English language, and complete preservation of all technical content.

---

## Table of Contents

1. [Current State](#1-current-state)
2. [Problems Identified](#2-problems-identified)
3. [Target Structure](#3-target-structure)
4. [Writing Standards](#4-writing-standards)
5. [Phase 1: Create Directory Structure](#5-phase-1-create-directory-structure)
6. [Phase 2: Guides](#6-phase-2-guides)
7. [Phase 3: Reference](#7-phase-3-reference)
8. [Phase 4: Methodology](#8-phase-4-methodology)
9. [Phase 5: LEB](#9-phase-5-leb)
10. [Phase 6: Development](#10-phase-6-development)
11. [Phase 7: Marketing](#11-phase-7-marketing)
12. [Phase 8: Deletions](#12-phase-8-deletions)
13. [Phase 9: Cross-References and Final Checks](#13-phase-9-cross-references-and-final-checks)
14. [Content-Specific Fixes](#14-content-specific-fixes)
15. [Rewrite Classification Per Document](#15-rewrite-classification-per-document)
16. [Summary](#16-summary)

---

## 1. Current State

### File Inventory (30 files, ~14,800 lines)

**Root directory (8 files, ~4,100 lines):**

| File | Lines | Language | Purpose |
|------|-------|----------|---------|
| `README.md` | 501 | EN | Main project README |
| `CHANGELOG.md` | 1,303 | EN | Release history |
| `AGENTS.md` | 79 | EN | AI agent instructions |
| `TODO.md` | 532 | IT | Roadmap (~98% complete) |
| `KEPLERIAN_TODO.md` | 551 | EN | Keplerian improvement catalog |
| `PLAN.md` | 179 | IT | Moshier removal plan (completed) |
| `TODO_PRECISION.md` | 415 | IT | Precision fix spec (~85% completed) |
| `PRECISION_RECAP.md` | 543 | IT | Precision fix record (completed) |

**docs/ directory (18 files, ~9,800 lines):**

| File | Lines | Language | Purpose |
|------|-------|----------|---------|
| `LEB_GUIDE.md` | 1,753 | EN | LEB complete technical guide |
| `LEB_PLAN.md` | 1,612 | EN | LEB implementation plan |
| `AYANAMSHA.md` | 997 | EN | 43 ayanamsha modes |
| `HOUSE_SYSTEMS.md` | 710 | EN | 19 house system algorithms |
| `PRECISION.md` | 698 | EN | Scientific precision specs |
| `PERFORMANCE_AND_LEB.md` | 693 | EN | Codebase overview + LEB design proposal |
| `PRECISION_TUNING.md` | 510 | EN | Precision configuration guide |
| `INTERPOLATED_APOGEE.md` | 434 | EN | Interpolated apogee/perigee user guide |
| `migration-guide.md` | 431 | EN | pyswisseph migration guide |
| `testing.md` | 390 | EN | Testing guide |
| `LANDING.md` | 336 | EN | Landing page marketing copy |
| `TRUE_LILITH_METHODS.md` | 266 | EN | True Lilith calculation methods |
| `methodology.md` | 251 | EN | Core methodology vs Swiss Ephemeris |
| `PLANET_CENTERS_SPK.md` | 250 | EN | Planet center SPK generation |
| `interpolated_perigee_methodology.md` | 208 | EN | Perigee calibration v2.2 |
| `FULL_RANGE_COVERAGE_PLAN.md` | 179 | EN | Minor body range extension plan |
| `REBOUND_BENEFITS.md` | 174 | EN | N-body integration benefits |
| `PYERFA_BENEFITS.md` | 141 | EN | PyERFA integration benefits |
| `methodology_lunar_apsides.md` | 119 | EN | Lunar apsides methodology |
| `VERIFICATION.md` | 85 | EN | Branch verification record |

**Other (2 files):**

| File | Lines | Language | Purpose |
|------|-------|----------|---------|
| `implementation_plans/perigee_recalibration.md` | 431 | IT | Perigee recalibration plan |
| `.pytest_cache/README.md` | 8 | EN | Auto-generated (not documentation) |

---

## 2. Problems Identified

1. **8 markdown files in root** — TODO, PLAN, KEPLERIAN_TODO, TODO_PRECISION, PRECISION_RECAP should not be at top level
2. **18 files in docs/ with no subdirectories** — flat structure, no logical grouping
3. **1 orphaned file** in `implementation_plans/` (separate directory for a single file)
4. **Mixed languages** — Italian (TODO.md, PLAN.md, TODO_PRECISION.md, PRECISION_RECAP.md, perigee_recalibration.md) and English
5. **Mixed document types** — user guides, internal specs, completed plans, TODO lists, historical recaps all mixed together
6. **Inconsistent naming** — `SCREAMING_CASE.md`, `snake_case.md`, `kebab-case.md`
7. **Inconsistent formatting** — some have TOC, some don't; varied header styles, metadata, tone
8. **Superseded content** — `PERFORMANCE_AND_LEB.md` Sections 3-4 are 100% duplicated by `LEB_GUIDE.md` and `LEB_PLAN.md`
9. **Stale content** — `TRUE_LILITH_METHODS.md` says "interpolated apogee is not yet implemented" (it is)
10. **Contradictory content** — `PRECISION.md` Section 9 describes ELP2000 perturbations as applied ("Stage 2") to True Node; `PRECISION_RECAP.md` documents that they are NOT applied (produces tens-of-degrees errors)
11. **No documentation index** — no entry point for navigating the docs

---

## 3. Target Structure

### Root (unchanged)

```
README.md          # UNCHANGED — stays as-is
CHANGELOG.md       # UNCHANGED
AGENTS.md          # UNCHANGED (update doc paths in Lunar Calibration Workflow if needed)
```

### docs/

```
docs/
├── README.md                              # NEW — Documentation index
│
├── guides/
│   ├── migration-guide.md                 # ← docs/migration-guide.md
│   └── precision-tuning.md               # ← docs/PRECISION_TUNING.md
│
├── reference/
│   ├── ayanamsha.md                       # ← docs/AYANAMSHA.md
│   ├── house-systems.md                   # ← docs/HOUSE_SYSTEMS.md
│   └── precision.md                       # ← docs/PRECISION.md (+ fix Stage 2 contradiction)
│
├── methodology/
│   ├── overview.md                        # ← docs/methodology.md
│   ├── lunar-apsides.md                   # ← docs/methodology_lunar_apsides.md
│   ├── interpolated-apogee.md             # ← docs/INTERPOLATED_APOGEE.md
│   ├── interpolated-perigee.md            # ← docs/interpolated_perigee_methodology.md
│   │                                      #    (+ Gaussian smoothing failure from perigee_recalibration.md)
│   ├── true-lilith.md                     # ← docs/TRUE_LILITH_METHODS.md (fix stale note)
│   ├── planet-centers-spk.md              # ← docs/PLANET_CENTERS_SPK.md
│   ├── pyerfa-integration.md              # ← docs/PYERFA_BENEFITS.md
│   └── rebound-integration.md             # ← docs/REBOUND_BENEFITS.md
│
├── leb/
│   ├── guide.md                           # ← docs/LEB_GUIDE.md
│   └── design.md                          # ← docs/LEB_PLAN.md
│
└── development/
    ├── testing.md                         # ← docs/testing.md
    ├── architecture-overview.md           # ← docs/PERFORMANCE_AND_LEB.md (unique sections only)
    ├── roadmap.md                         # ← TODO.md (translated to EN)
    ├── keplerian-improvements.md          # ← KEPLERIAN_TODO.md
    ├── full-range-coverage.md             # ← docs/FULL_RANGE_COVERAGE_PLAN.md
    └── precision-history.md               # NEW — unique content from TODO_PRECISION.md
                                           #        and PRECISION_RECAP.md (translated to EN)
```

### marketing/ (root)

```
marketing/
└── LANDING.md                             # ← docs/LANDING.md (no modifications)
```

---

## 4. Writing Standards

All documents in the new structure must follow these rules:

| Aspect | Standard |
|--------|----------|
| Language | English |
| File naming | `kebab-case.md` |
| Title | Single H1 (`#`) as document title |
| Header | 1-2 sentence description immediately after H1 |
| TOC | Present if document exceeds ~150 lines |
| Tone | Technical reference (3rd person, no diary entries, no "I tried") |
| Cross-refs | Relative paths updated to new structure |
| "vs Swiss Ephemeris" sections | Preserved where present, neutral tone |
| Formulas and code | Preserved 100% — no content loss |
| Tables | Preserved 100% |
| References | At the end of each document |
| Emoji | None |

### Templates by Category

**Guides (`docs/guides/`):**
```
# Title
Brief description.
## Prerequisites (if any)
## [Content sections]
## Troubleshooting (if applicable)
```

**Reference (`docs/reference/`):**
```
# Title
Brief description.
## Overview
## [Topic sections organized logically]
## Technical Implementation (if applicable)
## References
```

**Methodology (`docs/methodology/`):**
```
# Title
Brief description of the computational approach.
## Background
## Method
  ### [Subsections per step/layer]
## Precision and Validation
## Comparison with Swiss Ephemeris (if applicable)
## References
```

**Development (`docs/development/`):**
```
# Title
Brief description and current status.
## Current State
## [Topic sections]
## Open Items / Future Work (if applicable)
```

---

## 5. Phase 1: Create Directory Structure

Create the following directories:

```bash
mkdir -p docs/guides
mkdir -p docs/reference
mkdir -p docs/methodology
mkdir -p docs/leb
mkdir -p docs/development
mkdir -p marketing
```

Create `docs/README.md` — the documentation index:

```markdown
# LibEphemeris Documentation

Comprehensive documentation for LibEphemeris, a high-precision astronomical
ephemeris library providing Swiss Ephemeris-compatible API using NASA JPL
DE440/DE441 via Skyfield.

## Guides

- **[Migration Guide](guides/migration-guide.md)** — Migrating from pyswisseph to libephemeris
- **[Precision Tuning](guides/precision-tuning.md)** — Configuring optional dependencies and ephemeris files for maximum precision

## Reference

- **[Precision](reference/precision.md)** — Scientific precision specifications, models, and measured accuracy
- **[House Systems](reference/house-systems.md)** — Mathematical documentation for all 19 house systems
- **[Ayanamsha](reference/ayanamsha.md)** — Complete reference for all 43 sidereal zodiac modes

## Methodology

- **[Overview](methodology/overview.md)** — Principal computational differences vs Swiss Ephemeris
- **[Lunar Apsides](methodology/lunar-apsides.md)** — Perigee and apogee computation methodology
- **[Interpolated Apogee](methodology/interpolated-apogee.md)** — SE_INTP_APOG and SE_INTP_PERG guide
- **[Interpolated Perigee](methodology/interpolated-perigee.md)** — Passage-interpolated harmonic fitting (v2.2)
- **[True Lilith](methodology/true-lilith.md)** — Osculating lunar apogee calculation methods
- **[Planet Centers](methodology/planet-centers-spk.md)** — Barycenter vs planet center corrections for outer planets
- **[PyERFA Integration](methodology/pyerfa-integration.md)** — IAU standard nutation, precession, and obliquity via PyERFA
- **[REBOUND Integration](methodology/rebound-integration.md)** — N-body minor body orbit propagation

## Binary Ephemeris (LEB)

- **[Technical Guide](leb/guide.md)** — Complete LEB reference: format specification, reader, pipelines, commands
- **[Design](leb/design.md)** — Original implementation plan and file format specification

## Development

- **[Testing](development/testing.md)** — Running tests, expected failures, downstream project fixtures
- **[Roadmap](development/roadmap.md)** — Current project status and open tasks
- **[Keplerian Improvements](development/keplerian-improvements.md)** — Catalog of Keplerian fallback improvement opportunities
- **[Full Range Coverage](development/full-range-coverage.md)** — Extending minor body coverage across the DE441 range
- **[Precision History](development/precision-history.md)** — Record of precision fixes, investigations, and open opportunities
- **[Architecture Overview](development/architecture-overview.md)** — Codebase metrics, performance bottleneck analysis, future Rust port vision
```

---

## 6. Phase 2: Guides

### 6.1 `docs/guides/migration-guide.md`

- **Source:** `docs/migration-guide.md`
- **Rewrite level:** Light
- **Actions:**
  - Add standard header (1-2 sentence description after H1)
  - Update cross-references: `TRUE_LILITH_METHODS.md` → `../methodology/true-lilith.md`, `LEB_GUIDE.md` → `../leb/guide.md`
  - Verify all code snippets and tables are preserved
  - No content changes needed — already well-written, English, tutorial tone

### 6.2 `docs/guides/precision-tuning.md`

- **Source:** `docs/PRECISION_TUNING.md`
- **Rewrite level:** Light
- **Actions:**
  - Rename from SCREAMING_CASE
  - Add standard header
  - No cross-references to update (self-contained)
  - Verify all 14 code snippets, 8 tables, and external references preserved

---

## 7. Phase 3: Reference

### 7.1 `docs/reference/ayanamsha.md`

- **Source:** `docs/AYANAMSHA.md`
- **Rewrite level:** Light
- **Actions:**
  - Rename from SCREAMING_CASE
  - Add standard header
  - No cross-references to update (self-contained)
  - Preserve all 43 ayanamsha entries with their property tables, historical basis, and reference points
  - Preserve the Technical Implementation section and all 13 references

### 7.2 `docs/reference/house-systems.md`

- **Source:** `docs/HOUSE_SYSTEMS.md`
- **Rewrite level:** Light
- **Actions:**
  - Rename from SCREAMING_CASE
  - Add standard header
  - Preserve all mathematical formulas, algorithms, and the summary comparison table
  - Preserve all 7 references

### 7.3 `docs/reference/precision.md`

- **Source:** `docs/PRECISION.md`
- **Rewrite level:** Light + content fix
- **Actions:**
  - Rename from SCREAMING_CASE
  - Add standard header
  - **FIX Stage 2 contradiction in Section 9 (Lunar Points → True Node):**
    The current text describes ELP2000-82B perturbation corrections as "Stage 2" that IS applied to the True Node.
    Investigation (documented in `PRECISION_RECAP.md` TODO 8) proved that the ELP2000 series is designed for the
    mean node, not the geometric node computed from `h = r × v`. Applying it to the geometric node produced errors
    of tens of degrees. Add a note clarifying:
    - The perturbation series exists in the codebase (~900 lines, 170+ terms)
    - It is NOT currently applied to the True Node calculation
    - The geometric two-body method is used instead
    - Residual vs Swiss Ephemeris: ~8.9 arcsec (~0.14 deg max)
  - Preserve all other content (18 sections, all measurement tables, all formulas, all 14 references)

---

## 8. Phase 4: Methodology

All 8 documents in this section receive **complete rewrite** for uniform tone, structure, and style.

### 8.1 `docs/methodology/overview.md`

- **Source:** `docs/methodology.md`
- **Rewrite level:** Complete
- **Current issues:** British English ("prioritises"), some sections thinner than others
- **Actions:**
  - Rewrite to Methodology template with consistent section structure
  - Unify to American English
  - Preserve all 5 methodology areas (Ephemeris Foundation, Outer Planet Centers, Lunar Apsides, Delta T, Minor Bodies)
  - Preserve the summary comparison table and all 8 references
  - Update cross-references: `interpolated_perigee_methodology.md` → `interpolated-perigee.md`, `INTERPOLATED_APOGEE.md` → `interpolated-apogee.md`

### 8.2 `docs/methodology/lunar-apsides.md`

- **Source:** `docs/methodology_lunar_apsides.md`
- **Rewrite level:** Complete
- **Current issues:** Academic tone slightly different from overview.md
- **Actions:**
  - Rewrite to Methodology template
  - Preserve the measured discrepancy table, the 12,000-passage methodology description, and all 4 references
  - Update cross-references to new paths

### 8.3 `docs/methodology/interpolated-apogee.md`

- **Source:** `docs/INTERPOLATED_APOGEE.md`
- **Rewrite level:** Complete
- **Current issues:** Mixed formal/tutorial tone, "libephemeris" lowercase, extensive but inconsistently structured
- **Actions:**
  - Rewrite to Methodology template with uniform tone
  - Capitalize "LibEphemeris" consistently
  - Preserve ALL technical content: body ID table, osculating problem analysis, both implementation methods (Moshier analytical + ELP2000 perturbation series), asymmetric perturbation analysis, precision tables, decision matrix, code examples, smoothness comparison
  - Preserve all 10 references (7 primary + 3 astrological)
  - Update cross-references: `TRUE_LILITH_METHODS.md` → `true-lilith.md`, `PRECISION.md` → `../reference/precision.md`

### 8.4 `docs/methodology/interpolated-perigee.md`

- **Source:** `docs/interpolated_perigee_methodology.md`
- **Rewrite level:** Complete
- **Current issues:** Lab-notebook style with version history narrative
- **Actions:**
  - Rewrite to Methodology template
  - Restructure version history (v2.0, v2.1, v2.2) into a "Failed Approaches" section alongside the successful v2.2 method
  - **ADD the Gaussian smoothing failure** from `implementation_plans/perigee_recalibration.md` Section 1.4:
    - 8-year Gaussian window (2922 days)
    - Failure mechanism: apsidal line precesses ~325 degrees in 8 years, causing `atan2(sum_sin, sum_cos)` to collapse near zero
    - Result: ~-0.16 degrees perturbation instead of the correct ~8-25 degrees
    - This becomes a fourth failed approach (v0 / Gaussian), preceding v2.0/v2.1/v2.2
  - Preserve ALL technical content: three-layer architecture, 61-term series table (12 categories), coefficient changes v1→v2.2, precision evolution, calibration commands
  - Preserve all 4 references

### 8.5 `docs/methodology/true-lilith.md`

- **Source:** `docs/TRUE_LILITH_METHODS.md`
- **Rewrite level:** Complete
- **Current issues:** Stale content (line 189), mixed tones, reference to non-existent `LIST.md`
- **Actions:**
  - Rewrite to Methodology template
  - **REMOVE stale note** "The interpolated apogee is not yet implemented" — it IS implemented (documented in `interpolated-apogee.md`)
  - **REMOVE reference** to `LIST.md` (does not exist)
  - Preserve ALL technical content: both algorithms (eccentricity vector + orbital elements), all 7 perturbation corrections with amplitudes, precision measurements (500-date sampling), sources of residual differences, the osculating apogee paradox analysis
  - Preserve all 8 references
  - Update cross-references to new paths

### 8.6 `docs/methodology/planet-centers-spk.md`

- **Source:** `docs/PLANET_CENTERS_SPK.md`
- **Rewrite level:** Complete
- **Current issues:** DevOps/runbook style, different tone from other methodology docs
- **Actions:**
  - Rewrite to Methodology template, keeping operational sections (SPK generation, troubleshooting) clearly separated
  - Preserve ALL content: barycenter offset table, three-tier correction strategy, SPK generation process, prerequisites, source files, technical details (SPK Type 2, Chebyshev polynomial specs, chaining), coverage/precision, troubleshooting, changelog
  - Update implicit cross-reference to `methodology.md` Section 2

### 8.7 `docs/methodology/pyerfa-integration.md`

- **Source:** `docs/PYERFA_BENEFITS.md`
- **Rewrite level:** Complete
- **Current issues:** Product-doc style, minimal scientific context
- **Actions:**
  - Rewrite to Methodology template with added scientific context on nutation/precession
  - Preserve ALL content: three nutation model comparison, error growth table, all 8 exposed functions, cached nutation, installation, when-to-use decision matrix

### 8.8 `docs/methodology/rebound-integration.md`

- **Source:** `docs/REBOUND_BENEFITS.md`
- **Rewrite level:** Complete
- **Current issues:** Mirrors PYERFA in product-doc style
- **Actions:**
  - Rewrite to Methodology template, matching `pyerfa-integration.md` in structure
  - Preserve ALL content: Keplerian vs N-body precision table, error growth table, 8 API items, 3 integrator descriptions, installation, when-to-use decision matrix
  - Preserve all 6 references (implicit — currently none listed, but concepts reference Brouwer & Clemence etc.)

---

## 9. Phase 5: LEB

### 9.1 `docs/leb/guide.md`

- **Source:** `docs/LEB_GUIDE.md`
- **Rewrite level:** Light
- **Actions:**
  - Add standard header
  - Update cross-reference: `LEB_PLAN.md` → `design.md` (in the same directory)
  - Update any references to `PERFORMANCE_AND_LEB.md` → `../development/architecture-overview.md`
  - Preserve ALL content (1,753 lines — the most comprehensive single document)

### 9.2 `docs/leb/design.md`

- **Source:** `docs/LEB_PLAN.md`
- **Rewrite level:** Light
- **Actions:**
  - Add standard header
  - Update cross-reference: `PERFORMANCE_AND_LEB.md` → `../development/architecture-overview.md`
  - Preserve ALL content (1,612 lines)

---

## 10. Phase 6: Development

### 10.1 `docs/development/testing.md`

- **Source:** `docs/testing.md`
- **Rewrite level:** Light
- **Actions:**
  - Add standard header
  - Preserve ALL content: test statuses, expected failures, skip categories, SPK fixture patterns, CI/CD configuration, troubleshooting

### 10.2 `docs/development/architecture-overview.md`

- **Source:** `docs/PERFORMANCE_AND_LEB.md`
- **Rewrite level:** Complete (restructure)
- **Current issues:** Marked "DESIGN PROPOSAL / NOT YET IMPLEMENTED" — misleading since LEB is implemented. Sections 3-4 are 100% duplicated by `leb/guide.md` and `leb/design.md`.
- **Actions:**
  - **PRESERVE Section 1 (Codebase Overview)** — unique content:
    - Project metrics table (37 source files, ~69,000 LOC, 519 public API entries, etc.)
    - Architecture summary
    - Runtime dependencies table
    - Module breakdown by LOC (24 modules)
  - **PRESERVE Section 2 (Performance Bottleneck Analysis)** — unique content:
    - Cost-per-operation table (swe_calc_ut ~1ms, houses ~500us, eclipse ~500ms-2s, full chart ~12ms)
    - CPU time breakdown (5 bottlenecks ranked: Skyfield .observe() 90%+, central difference 3x, eclipse combinatorial, no caching, erfa.nut06a)
    - Existing caching infrastructure table
    - numpy hot-path analysis
  - **REMOVE Section 3 (Precomputed Binary Format)** — 100% duplicated by `leb/guide.md` Section 3 and `leb/design.md` Section 3, with outdated values (48-byte BodyEntry vs actual 52-byte)
  - **REMOVE Section 4 (Python-First Implementation Plan)** — 100% duplicated by `leb/guide.md` Sections 2,4-7,11 and `leb/design.md` Sections 4-9,11
  - **PRESERVE Section 5 (Future Rust Port via PyO3)** — partially unique:
    - The full-library Rust port vision (~7,800 lines including houses, eclipses, crossing, lunar, stars, hypothetical) exists only here
    - The LEB-reader-only Rust port (~1,200 lines) is also in `leb/design.md` Section 10, but the ambitious scope is unique
    - Preserve the detailed crate structure, dependency list, speedup projections, 15-week timeline
    - Add note: "The LEB-reader-only Rust port is also documented in `../leb/design.md` Section 10"
  - **PRESERVE Appendix (Full Module Inventory)** — unique content:
    - A.1: Coordinate frames by body type (11 body categories)
    - A.2: Numerical algorithms used (14 algorithms with parameters)
    - A.3: House systems implemented (24 systems with properties)
    - A.4: Ayanamsha systems categorization (43 systems)
  - **UPDATE header:** Remove "STATUS: DESIGN PROPOSAL / NOT YET IMPLEMENTED", replace with contextual description explaining this document contains the codebase overview, performance analysis, and future Rust vision
  - **ADD cross-references** to `../leb/guide.md` and `../leb/design.md` for LEB-specific documentation

### 10.3 `docs/development/roadmap.md`

- **Source:** `TODO.md` (root)
- **Rewrite level:** Complete (translate IT→EN + restructure)
- **Current issues:** Italian language, ~98% completed tasks with only 3 minor sub-items open
- **Actions:**
  - Translate entire document to English
  - Restructure to Development template with clear completed/open status for each item
  - Mark completed tasks clearly (with brief summary of what was done)
  - Highlight the 3 open items prominently:
    1. ASSIST data files download + end-to-end verification (awaiting user approval)
    2. Extended tier LEB file not yet generated
    3. ASSIST performance verification for single evaluations
  - Preserve the architectural context, solution descriptions, and file references
  - Remove the "Ordine di implementazione raccomandato" checklist (all items checked off)

### 10.4 `docs/development/keplerian-improvements.md`

- **Source:** `KEPLERIAN_TODO.md` (root)
- **Rewrite level:** Light (structural uniformity)
- **Current issues:** Already in English, good structure, but naming and style differ from target
- **Actions:**
  - Rename from SCREAMING_CASE
  - Add standard header
  - Verify the document is entirely forward-looking (no completed items — confirmed)
  - Preserve ALL content: precision table, implementation summary, 6 low-cost items (L1-L6), 4 medium-cost items (M1-M4), 5 high-cost items (H1-H5), 4 structural issues (S1-S4), priority recommendations, 6 references

### 10.5 `docs/development/full-range-coverage.md`

- **Source:** `docs/FULL_RANGE_COVERAGE_PLAN.md`
- **Rewrite level:** Complete
- **Current issues:** Bug report / plan hybrid, pseudocode-heavy
- **Actions:**
  - Rewrite to Development template with clear problem/solution/status structure
  - Preserve ALL content: problem statement, current situation table, target coverage diagram, 3 bugs found (with file:line references), 5 implementation steps (with pseudocode), JPL Horizons SPK limits verification
  - Clarify status of each implementation step (which are done, which are pending)

### 10.6 `docs/development/precision-history.md`

- **Source:** NEW (content extracted from `TODO_PRECISION.md` and `PRECISION_RECAP.md`)
- **Rewrite level:** New document (translated IT→EN)
- **Purpose:** Preserve all unique technical content that would otherwise be lost when deleting the source files
- **Structure:**

```markdown
# Precision History

Record of precision improvements applied to LibEphemeris (February 2026),
including investigation results and architectural decisions that inform
future development.

## Completed Fixes

### PREC_RATE / PREC_RATE_QUAD Regression Bug (Critical)
- Mechanism: New _PREC_C1-_PREC_C5 constants defined but old PREC_RATE/PREC_RATE_QUAD
  names not updated → NameError crash for ALL formula-based ayanamshas
- Three fix locations with before/after (SE_SIDM_SURYASIDDHANTA_MSUN, general formula,
  SE_SIDM_J2000)
- Source: PRECISION_RECAP TODO 4

### Nutation Unification to IAU 2006/2000A
- Four code paths still using iau2000b_radians (77 terms):
  planets.py:_get_true_ayanamsa(), planets.py:_calc_ayanamsa_ex(),
  utils.py:azalt(), utils.py:azalt_rev()
- All replaced with erfa.nut06a() (1365 terms)
- Source: PRECISION_RECAP TODO 2

### Obliquity Unification to IAU 2006
- Three inconsistent formulas: cache.py (erfa.obl06), _calc_ayanamsa_ex (Laskar 1986,
  84381.448"), azalt (IAU 2006 constant 84381.406" but only 3 polynomial terms)
- 0.042" constant offset eliminated
- Source: PRECISION_RECAP TODO 3

### Precession Upgrade to IAU 2006 + Frame Bias
- Original code claimed "IAU 2006" in comments but actually used Lieske 1977 coefficients
- Intermediate fix: erfa.pmat06()
- Final fix: Skyfield pipeline rewrite (as part of aberration fix below)
- Source: PRECISION_RECAP TODO 5+6

### Annual Aberration in Star-Based Ayanamsha
- ~20.5" error in _get_star_position_ecliptic()
- Rewrite: ~130 lines of manual calculation → ~15 lines using Skyfield
  observe().apparent().frame_latlon()
- Automatically includes: aberration, light-time, gravitational deflection,
  IAU 2006 precession, IAU 2000A nutation
- Source: PRECISION_RECAP TODO 1

### Central Difference Velocity
- 7 forward-difference calculations converted to central-difference across 4 files
  (hypothetical.py, spk.py, planetary_moons.py, planets.py)
- ~100x precision improvement per timestep
- Source: PRECISION_RECAP TODO 7

## Investigations: ELP2000 Perturbations Not Applied

### True Node (SE_TRUE_NODE) — ~8.9" residual vs Swiss Ephemeris
- The ELP2000-82B perturbation series (~170 terms, ~900 lines) exists in the codebase
  but is designed for the MEAN node, not the geometric node computed from h = r × v
- Applying the series to the geometric node produced errors of tens of degrees
- Decision: NOT APPLIED. The geometric two-body method is used instead
- NOTE: docs/reference/precision.md Section 9 "Stage 2" describes these perturbations
  as applied — this is incorrect and should be fixed (see Phase 3 above)
- Future approaches: calibrate systematic offset vs SWE, or implement three-body
  tidal force model
- Source: PRECISION_RECAP TODO 8

### True Lilith (SE_OSCU_APOG) — ~54" residual vs Swiss Ephemeris
- Same conclusion: the ELP2000 apogee perturbation series is designed for the
  interpolated (mean) apogee, not the osculating eccentricity vector
- Decision: NOT APPLIED
- No clear improvement path identified
- Source: PRECISION_RECAP TODO 9

## Open Opportunities

### Analytical Chebyshev Velocities (TODO 7a — not implemented)
- Skyfield/jplephem provides analytical derivatives of JPL Chebyshev polynomials
- Would eliminate 2 recursive _calc_body() calls for velocity computation
- Expected: ~3x performance improvement for swe_calc_ut() with SEFLG_SPEED
- Expected: sub-arcsecond/day Moon velocity improvement
- Source: TODO_PRECISION TODO 7

### True Node Precision Improvement
- Current residual: ~8.9" (~0.14 deg max) vs Swiss Ephemeris
- Proposed approaches:
  1. Calibrate systematic offset against SWE reference values
  2. Implement osculating element integration (not ELP2000 term filtering)
- Source: PRECISION_RECAP TODO 8 future work

### True Lilith Precision Improvement
- Current residual: ~54" (~0.065 deg max) vs Swiss Ephemeris
- No clear path forward — the osculating apogee concept is inherently model-dependent
- Source: PRECISION_RECAP TODO 9

## Overall Precision Impact (February 2026 Fixes)

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Star-based ayanamsha | ~20.5" error | <0.01" | ~200x |
| Nutation consistency | ~1 mas inconsistency | <0.05 mas | ~20x |
| Obliquity offset | 0.042" systematic | 0 | eliminated |
| Stellar precession | Lieske 1977 | IAU 2006 | ~1000x |
| Velocity precision | O(h) forward diff | O(h²) central diff | ~100x |
| Formula-based ayanamshas | NameError crash | working | critical fix |
```

- **Content extracted from:**
  - `PRECISION_RECAP.md`: All fix details, investigation records (TODOs 8+9), impact table, accidentally deleted code incident
  - `TODO_PRECISION.md`: Analytical velocity proposal (TODO 7a), open improvement targets (TODOs 8+9)
- **Content NOT extracted (redundant):**
  - Problem descriptions that duplicate `docs/reference/precision.md` final-state documentation
  - Dependency graph between TODOs (all completed)
  - Test failure analysis (759 snapshot tests — transient issue)

---

## 11. Phase 7: Marketing

### 11.1 `marketing/LANDING.md`

- **Source:** `docs/LANDING.md`
- **Rewrite level:** None
- **Actions:**
  - Move to `marketing/LANDING.md`
  - No content modifications

---

## 12. Phase 8: Deletions

| File | Reason |
|------|--------|
| `PLAN.md` | 100% redundant — "SEFLG_MOSEPH silently ignored" documented in 7+ other files. All 12 implementation phases completed and recorded in CHANGELOG. The before/after impact table describes the current state, which is now the default. |
| `TODO_PRECISION.md` | All unique content extracted to `docs/development/precision-history.md`. Remaining content (problem descriptions, fix specs) is redundant with `docs/reference/precision.md` which documents the final state. |
| `PRECISION_RECAP.md` | All unique content extracted to `docs/development/precision-history.md`. Historical implementation details are preserved in git history. |
| `docs/VERIFICATION.md` | Branch-specific verification record (`dev/fix-tests`, commits `a59b110` and `c96fb46`). No unique technical content — changes described are covered by CHANGELOG and the code itself. |
| `implementation_plans/perigee_recalibration.md` | Three-layer architecture superseded by `docs/methodology/interpolated-perigee.md` (actual v2.2 with 61 terms vs planned 43). Gaussian smoothing failure integrated into interpolated-perigee.md. Calibration workflow documented in AGENTS.md. |
| `implementation_plans/` | Directory empty after above deletion. |
| `docs/PERFORMANCE_AND_LEB.md` | Content migrated: unique sections (1, 2, 5, Appendix) → `docs/development/architecture-overview.md`. Duplicate sections (3, 4) discarded. Original file no longer needed. |

After deletions, also remove:
- All original files in `docs/` that were moved to subdirectories (the old flat-structure files)

---

## 13. Phase 9: Cross-References and Final Checks

### 13.1 Cross-Reference Updates

All internal cross-references must be updated. Known references:

| Old Path | New Path | Referenced From |
|----------|----------|-----------------|
| `TRUE_LILITH_METHODS.md` | `methodology/true-lilith.md` | migration-guide.md, INTERPOLATED_APOGEE.md |
| `LEB_GUIDE.md` | `leb/guide.md` | migration-guide.md, LEB_PLAN.md |
| `LEB_PLAN.md` | `leb/design.md` | LEB_GUIDE.md |
| `PERFORMANCE_AND_LEB.md` | `development/architecture-overview.md` | LEB_PLAN.md |
| `PRECISION.md` | `reference/precision.md` | INTERPOLATED_APOGEE.md, TRUE_LILITH_METHODS.md |
| `interpolated_perigee_methodology.md` | `methodology/interpolated-perigee.md` | methodology.md, methodology_lunar_apsides.md |
| `INTERPOLATED_APOGEE.md` | `methodology/interpolated-apogee.md` | methodology.md, methodology_lunar_apsides.md |

### 13.2 AGENTS.md Update

Check if `AGENTS.md` references any doc paths that changed (specifically the "Lunar Calibration Workflow" section and "See `docs/LEB_PLAN.md` and `docs/LEB_GUIDE.md`" references). Update paths if needed.

### 13.3 Final Verification Checklist

For every document in the new structure, verify:

- [ ] Standard header present (1-2 sentence description after H1)
- [ ] TOC present if document exceeds ~150 lines
- [ ] Language is English throughout
- [ ] Tone is consistent with category (guides = tutorial, reference = formal, methodology = scientific, development = practical)
- [ ] All formulas, tables, code snippets, and measurements preserved
- [ ] All internal cross-references use correct relative paths
- [ ] References section present at end (where applicable)
- [ ] No stale content or contradictions
- [ ] File naming follows `kebab-case.md`
- [ ] "LibEphemeris" capitalized consistently (not "libephemeris")

---

## 14. Content-Specific Fixes

These are specific content corrections to be applied during rewriting:

| Document | Fix | Detail |
|----------|-----|--------|
| `reference/precision.md` | Fix Stage 2 contradiction | Section 9 True Node: clarify ELP2000 perturbations are NOT applied (see precision-history.md for investigation) |
| `methodology/true-lilith.md` | Remove stale note | Line 189 "interpolated apogee is not yet implemented" → remove (it IS implemented) |
| `methodology/true-lilith.md` | Remove dead reference | Remove reference to non-existent `LIST.md` |
| `methodology/interpolated-perigee.md` | Add Gaussian failure | Add failed Gaussian smoothing approach (8-year window, apsidal precession collapse) from perigee_recalibration.md Section 1.4 |
| `development/architecture-overview.md` | Remove misleading header | Replace "STATUS: DESIGN PROPOSAL / NOT YET IMPLEMENTED" with accurate description |
| `development/architecture-overview.md` | Remove duplicate sections | Remove Sections 3-4 (LEB format + Python plan) — 100% duplicated by leb/guide.md and leb/design.md with outdated values |

---

## 15. Rewrite Classification Per Document

### Light Rewrite (8 documents)

These documents are already well-written in English with appropriate tone. Changes limited to: header standardization, file renaming, cross-reference updates.

| Document | Source | Lines |
|----------|--------|-------|
| `guides/migration-guide.md` | `docs/migration-guide.md` | 431 |
| `guides/precision-tuning.md` | `docs/PRECISION_TUNING.md` | 510 |
| `reference/ayanamsha.md` | `docs/AYANAMSHA.md` | 997 |
| `reference/house-systems.md` | `docs/HOUSE_SYSTEMS.md` | 710 |
| `reference/precision.md` | `docs/PRECISION.md` | 698 |
| `leb/guide.md` | `docs/LEB_GUIDE.md` | 1,753 |
| `leb/design.md` | `docs/LEB_PLAN.md` | 1,612 |
| `development/testing.md` | `docs/testing.md` | 390 |

### Complete Rewrite (8 documents)

These documents need tone/style unification, structural restructuring, or language translation. All technical content preserved.

| Document | Source | Lines | Reason |
|----------|--------|-------|--------|
| `methodology/overview.md` | `docs/methodology.md` | 251 | Inconsistent British/American English, style unification |
| `methodology/lunar-apsides.md` | `docs/methodology_lunar_apsides.md` | 119 | Tone alignment with overview |
| `methodology/interpolated-apogee.md` | `docs/INTERPOLATED_APOGEE.md` | 434 | Mixed formal/tutorial tone, lowercase "libephemeris" |
| `methodology/interpolated-perigee.md` | `docs/interpolated_perigee_methodology.md` | 208 | Lab-notebook style, + Gaussian failure addition |
| `methodology/true-lilith.md` | `docs/TRUE_LILITH_METHODS.md` | 266 | Stale content, mixed tones, dead references |
| `methodology/planet-centers-spk.md` | `docs/PLANET_CENTERS_SPK.md` | 250 | DevOps/runbook style differs from methodology standard |
| `methodology/pyerfa-integration.md` | `docs/PYERFA_BENEFITS.md` | 141 | Product-doc style, needs scientific context |
| `methodology/rebound-integration.md` | `docs/REBOUND_BENEFITS.md` | 174 | Product-doc style, needs scientific context |

### Complete Rewrite + Translation (2 documents)

| Document | Source | Lines | Reason |
|----------|--------|-------|--------|
| `development/roadmap.md` | `TODO.md` | 532 | Italian → English, restructure completed/open items |
| `development/full-range-coverage.md` | `docs/FULL_RANGE_COVERAGE_PLAN.md` | 179 | Bug report hybrid style, clarify status |

### Complete Restructure (1 document)

| Document | Source | Lines | Reason |
|----------|--------|-------|--------|
| `development/architecture-overview.md` | `docs/PERFORMANCE_AND_LEB.md` | 693 | Remove duplicate sections 3-4, keep unique sections 1-2, 5, Appendix |

### New Documents (2)

| Document | Lines (est.) | Source |
|----------|-------------|--------|
| `docs/README.md` | ~60 | New (documentation index) |
| `development/precision-history.md` | ~150 | Extracted from TODO_PRECISION.md + PRECISION_RECAP.md |

### No Modification (2)

| Document | Source |
|----------|--------|
| `marketing/LANDING.md` | `docs/LANDING.md` (move only) |
| `development/keplerian-improvements.md` | `KEPLERIAN_TODO.md` (move + light header) |

---

## 16. Summary

### Quantitative Overview

| Operation | Count |
|-----------|-------|
| Directories created | 6 |
| Documents with light rewrite | 8 |
| Documents with complete rewrite | 8 |
| Documents translated IT→EN | 2 |
| Documents restructured (remove duplicates) | 1 |
| New documents created | 2 |
| Documents moved without modification | 2 |
| Documents deleted | 7 (+ 1 empty directory) |
| Root files unchanged | 3 (README.md, CHANGELOG.md, AGENTS.md) |
| Cross-references to update | ~15-20 |
| Content-specific fixes | 6 |
| **Total documents in new structure** | **21** |

### Execution Order

1. Create directory structure + `docs/README.md` index
2. Process `guides/` (2 files — low complexity)
3. Process `reference/` (3 files — low complexity, fix Stage 2)
4. Process `methodology/` (8 files — high complexity, complete rewrites)
5. Process `leb/` (2 files — low complexity)
6. Create `development/precision-history.md` (new, content extraction + translation)
7. Create `development/architecture-overview.md` (restructured from PERFORMANCE_AND_LEB.md)
8. Translate and rewrite `development/roadmap.md` (from TODO.md)
9. Process remaining `development/` files (3 files)
10. Move `LANDING.md` → `marketing/`
11. Delete obsolete files
12. Update all cross-references
13. Update `AGENTS.md` paths if needed
14. Final verification: every document passes the checklist

### Invariants

- **README.md** (root) is NOT modified
- **CHANGELOG.md** is NOT modified
- **100% of technical content** is preserved (formulas, tables, code, measurements, references)
- **No document is deleted** unless its unique content has been extracted to a new location
- **All cross-references** are updated to reflect new paths
