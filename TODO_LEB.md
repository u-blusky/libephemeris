# LEB Generator TODO

## Completed

- [x] Vectorized Skyfield ICRS evaluations (~150x speedup)
- [x] Vectorized nutation via direct erfa (~50x speedup)
- [x] Batched verification (fit nodes + verify in single array)
- [x] spktype21 integration for asteroids (~36x vs scalar fallback)
- [x] Parallelized analytical bodies via ProcessPoolExecutor
- [x] Linear extrapolation for SPK boundary overshoot
- [x] Default `--workers` to CPU count
- [x] Progress bars per-body (lightweight, zero dependencies)
- [x] Full-range asteroid SPK download (passes tier jd_start/jd_end)
- [x] Cached narrow-range SPK detection + force re-download
- [x] Removed Swiss Ephemeris naming from generator labels/comments

## Pending

### High priority

- [ ] **Regenerate `ephemeris_base.leb`**
  Run `poe leb:generate:base` with the updated code so asteroids use the
  spktype21 path instead of the Keplerian fallback. The current file was
  generated before the optimizations and has ~1500" errors on asteroids.

- [ ] **Validate asteroid precision after regeneration**
  Quick check: Chiron, Ceres, Pallas, Juno, Vesta should have position
  errors <1" against the Skyfield reference with SPK registered.
  ```bash
  python3 -c "
  from libephemeris.leb_reader import LEBReader
  from libephemeris.fast_calc import fast_calc_ut
  from libephemeris.constants import *
  import libephemeris as ephem, os
  os.environ['LIBEPHEMERIS_STRICT_PRECISION'] = '0'
  ephem.set_jpl_file('de440s.bsp')
  ephem.set_auto_spk_download(True)
  ephem.set_strict_precision(False)
  reader = LEBReader('data/leb/ephemeris_base.leb')
  for ipl, name in [(SE_CHIRON,'Chiron'),(SE_CERES,'Ceres'),(SE_VESTA,'Vesta')]:
      fast, _ = fast_calc_ut(reader, 2451545.0, ipl, SEFLG_SPEED)
      ref, _ = ephem.swe_calc_ut(2451545.0, ipl, SEFLG_SPEED)
      print(f'{name}: dLon={abs(fast[0]-ref[0])*3600:.4f} arcsec')
  reader.close()
  "
  ```

### Medium priority

- [ ] **Regenerate medium and extended tiers** (if needed)
  `poe leb:generate:medium` and `poe leb:generate:extended` with the same
  optimizations. Only needed if those tiers are in use.

- [ ] **Run full precision test suite**
  `poe test:leb:precision` on regenerated files. Currently slow (~hours).
  Consider reducing sample dates or adding `-n` parallelism.

### Low priority

- [ ] **Benchmark full base tier generation time**
  Document the final generation time for 300yr base tier with the
  spktype21 path active for all asteroids (expected ~3 min vs ~18 min
  with the old scalar fallback).

## Performance reference

| Configuration                        | Time       |
|--------------------------------------|------------|
| Old scalar (1 worker, 300yr)         | ~18 min    |
| New vectorized (10yr, 1 worker)      | 18.0s      |
| New vectorized (10yr, 4 workers)     | 6.4s       |
| New vectorized (300yr, 4w, no ast.)  | 170s       |
| New with spktype21 asteroids (300yr) | TBD        |
