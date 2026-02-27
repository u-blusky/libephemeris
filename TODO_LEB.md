# LEB Generator TODO

## Completed

- [x] Vectorized Skyfield ICRS evaluations (~150x speedup)
- [x] Vectorized nutation via direct erfa (~50x speedup)
- [x] Batched verification (fit nodes + verify in single array)
- [x] spktype21 integration for asteroids (~36x vs scalar fallback)
- [x] Linear extrapolation for SPK boundary overshoot
- [x] Default `--workers` to CPU count
- [x] Progress bars per-body (lightweight, zero dependencies)
- [x] Full-range asteroid SPK download (passes tier jd_start/jd_end)
- [x] Cached narrow-range SPK detection + force re-download
- [x] Removed Swiss Ephemeris naming from generator labels/comments
- [x] Fixed macOS fork deadlock — removed ProcessPoolExecutor, analytical
  bodies now run sequentially (deadlock-free, ~2-3 min for base tier)
- [x] `--group {planets,asteroids,analytical}` CLI option for independent
  generation of each body group
- [x] `--merge FILE [FILE ...]` CLI option to merge partial .leb files
- [x] `merge_leb_files()` function — validates JD range, checks duplicates,
  copies raw coefficient blobs (zero re-computation)
- [x] `BODY_GROUPS` dict mapping group names to body ID lists
- [x] New poe commands for per-group generation + merge (all 3 tiers)
- [x] `docs/LEB_GUIDE.md` — comprehensive technical guide (Section 6.15
  documents the group/merge workflow)
- [x] Fixed `verify_leb()` AU-to-arcsec conversion — now divides by body
  distance for correct angular error; threshold raised to 1" (was 10 mas)
- [x] Regenerated `ephemeris_base.leb` via group workflow — all 31 bodies
  pass verification (planets <0.35", asteroids <1", analytical 0.00")

## Pending

### Medium priority

- [ ] **Regenerate medium and extended tiers** (if needed)
  Use the group workflow:
  ```bash
  poe leb:generate:medium:groups
  poe leb:generate:extended:groups
  ```
  Only needed if those tiers are in use.

- [ ] **Run full precision test suite**
  `poe test:leb:precision` on regenerated files. Currently slow (~hours).
  Consider reducing sample dates or adding `-n` parallelism.

### Low priority

- [ ] **Benchmark full base tier generation time**
  Document the final generation time for 300yr base tier with the
  group workflow (planets ~20s, asteroids ~minutes, analytical ~minutes).

## Performance reference

| Configuration                        | Time       |
|--------------------------------------|------------|
| Old scalar (1 worker, 300yr)         | ~18 min    |
| New vectorized (10yr, 1 worker)      | 18.0s      |
| New vectorized (10yr, 4 workers)     | 6.4s       |
| New vectorized (300yr, 4w, no ast.)  | 170s       |
| Group workflow: planets (300yr)      | ~20s       |
| Group workflow: analytical (5yr)     | ~8s        |
| Group workflow: merge                | ~0.1s      |
