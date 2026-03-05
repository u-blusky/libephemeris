# LEB Precision Improvement — TODO

> Reference plan: `docs/leb/precision-improvement-plan.md`
> Branch: `leb-precision-improvement`

## Context

The LEB (LibEphemeris Binary) precision improvement work is split into 5 phases (0-4), executed in order 0 → 1 → 2 → 4 → 3.

The **base tier** (de440s, 1850-2150) is complete: all 5 phases executed, LEB file regenerated (94.3 MB), tests pass 404/404 with tight tolerances.

The **medium tier** (de440, 1550-2650) is complete: all 5 phases executed, LEB file regenerated (315 MB), tests pass 976/976 (+ 11 xfail, 12 skip) with tight tolerances.

### What was done (base tier)

- **Phase 0**: Removed xfails from asteroid tests, enabled SPK auto-download in CompareHelper
- **Phase 1**: Aggressive Chebyshev parameters (Moon/Earth 4d, Venus/Mars 16d, Jupiter/Saturn 32d, Asteroids 8d, Hypotheticals 32d, Nutation 16d)
- **Phase 2**: Analytical velocity via Chebyshev derivative transformed through rotation matrices (replaces central difference)
- **Phase 4**: Base LEB file regenerated with all 30 bodies verified
- **Phase 3**: Tolerances tightened to minimum, documentation updated (guide.md + design.md)
- **Bug fix**: Tau bug in verify_segment of the generator, arcsec reporting with geocentric distance

### What was done (medium tier)

- **Phase 4**: Medium LEB file regenerated (315 MB, 31 bodies, all verified)
- **Phase 3**: Tolerances tightened to minimum based on measured errors
- **Asteroid date filtering**: SPK range restricted to 1900-2100 CE (JPL Horizons SPK21 files only cover this range; outside it the LEB generator baked-in wrong Keplerian data into Chebyshev coefficients)
- **Asteroid velocity tolerances**: Separate tolerances for asteroids (`ASTEROID_SPEED_LON/LAT/DIST_DEG_DAY`) in `test_compare_leb_velocities.py` and `test_compare_leb_asteroids.py`
- **Crossing solver fixes**: Catch `RuntimeError` for Mars 180° and Saturn helio 180°/270° (pre-existing bug in `crossing.py`, not LEB-related)

### Measured errors (base tier, worst case)

- Planet position: from 0.0002" (Sun) to 4.85" (Saturn lat) — sub-arcsecond for all except Saturn lat
- Planet velocity: < 0.014 deg/day (Saturn lon is worst)
- Asteroid position: ~0.44"
- Asteroid lat velocity: 0.19-0.71 deg/day (architectural limit of ICRS→ecliptic pipeline)
- Ecliptic bodies: < 0.028"
- Hypotheticals: ~0 (1e-14)

### Measured errors (medium tier, worst case)

- Planet position: from 0.0003" (Sun) to 4.58" (Uranus lon ~1900 CE)
- Planet velocity: < 0.0014 deg/day (Moon lon)
- Lat velocity: < 0.0029 deg/day (OscuApogee/InterpApogee)
- Distance: < 1.92e-5 AU (Pluto)
- Asteroid position: ~0.29" (filtered to 1900-2100 CE)
- Asteroid lat velocity: 0.34 deg/day (Pallas, architectural limit)
- Ecliptic bodies: < 0.035" (OscuApogee)
- Hypotheticals: ~0
- Equatorial: < 0.37" (Uranus)

### Known architectural limits (not fixable without changing format)

1. **Uranus/Saturn position (~5")** — The ICRS→ecliptic pipeline amplifies errors by 1/geocentric_distance. Worst at extreme dates (1900 CE for Uranus in medium tier, entire range for Saturn in base tier).
2. **Asteroid latitude velocity (0.2-0.7 deg/day)** — Same mechanism, worse for nearby bodies
3. **Pluto distance velocity (~2.3e-5 AU/day)** — Combination of extreme distance and eccentric orbit
4. **Asteroids outside SPK range** — JPL Horizons SPK21 files only cover ~1900-2100 CE. The LEB generator produces wrong Keplerian data outside this range, baked-in to Chebyshev coefficients. Tests filtered with `filter_asteroid_dates()`.

---

## TODO

### Medium Tier (de440, 1550-2650)

- [x] Regenerate medium LEB file: `poe leb:generate:medium:groups`
- [x] Copy generated file to `/Volumes/Data/libephemeris/leb/ephemeris_medium.leb`
- [x] Run medium tests: `pytest tests/test_leb/compare/ -m leb_compare -v`
- [x] Tighten `TIER_DEFAULTS["medium"]` tolerances to minimum based on observed errors
- [x] Handle asteroid failures (SPK coverage restricted to 1900-2100 CE)
- [x] Verify all tests pass with tight tolerances (976 passed, 11 xfailed)

### Extended Tier (de441, -5000 to +5000) — Optional

- [ ] Regenerate extended LEB file: `poe leb:generate:extended:single`
- [ ] Run extended tests and tighten tolerances
- [ ] Handle asteroid SPK range (only covers ~1900-2100)

### Cleanup

- [x] Verify `poe typecheck` passes (mypy) — 0 errors
- [x] Translate `docs/leb/testing.md` to English
- [x] Add `poe test:leb:compare:medium` alias
- [x] Add `--single` mode for body-by-body LEB generation (lowest memory)
- [x] Add `--skip-aux` optimization (nutation/delta-T/stars generated only once)
- [x] Add nutation progress bar

### Release

- [ ] Upload updated LEB files to GitHub Releases: `poe release:leb <version>`
- [ ] Update hashes in `libephemeris/download.py`
- [ ] Commit updated hashes
- [ ] Merge branch `leb-precision-improvement` to main
