# TODO - House Systems Improvements

## High Priority
- [ ] **`gauquelin_sector()` API Alignment**: Align the signature of `gauquelin_sector` and `swe_gauquelin_sector` with the original Swiss Ephemeris API.
  - Current LE: `(jd, planet, lat, lon, altitude, pressure, temperature, flags, method)`
  - Original SE: `(tjdut, body, method, geopos, atpress, attemp, flags)`
- [ ] **Fix Skipped Tests**: Enable and fix the 5 skipped tests in `compare_scripts/tests/test_houses/test_gauquelin_sector.py` and `test_compare_houses_ext.py` once the API is aligned.

## Precision Improvements
- [ ] **Campanus `house_pos`**: Investigate the coordinate transformation discrepancy in `house_pos` for the Campanus system.
  - Current tolerance: 0.7° (vs 0.02° for other systems).
  - Goal: Sub-arcminute precision matching Swiss Ephemeris.

## Code Quality & Maintenance
- [ ] **Type Annotations**: Fix pre-existing LSP errors and type mismatches in `libephemeris/houses.py` (e.g., lines 561, 1976, 4037).
- [ ] **Methods 2-5 for Gauquelin**: Implement precise rise/set algorithms for Gauquelin sector calculation methods 2 through 5 (currently using approximations).
