# TODO - House Systems Improvements

## Completed

- [x] **`swe_gauquelin_sector()` API Alignment**: Aligned the signature with Swiss Ephemeris API.
  - Now uses: `(tjdut, body, method, geopos, atpress, attemp, flags)`
  - The `gauquelin_sector()` function retains the original libephemeris-style signature for convenience.

- [x] **Methods 2-5 for Gauquelin**: Implemented precise rise/set algorithms for Gauquelin sector calculation methods 2-5.
  - Uses actual rise/set times via `rise_trans()` function
  - Maximum difference vs Swiss Ephemeris: ~0.02 sectors (essentially exact match)
  - Methods 2: disc center, no refraction
  - Methods 3: disc center, with refraction
  - Methods 4: disc edge (upper limb), no refraction
  - Methods 5: disc edge (upper limb), with refraction

- [x] **Type Annotations**: Fixed critical type issues in `libephemeris/houses.py`:
  - Line 560: `float(t.gast)` to convert numpy.float64 to Python float
  - Line 1746: Fixed "possibly unbound" error in RA calculation
  - Lines 4043, 4340: Changed `Union[float, bytes, str] = None` to `Optional[...]`

- [x] **Campanus `house_pos`**: Fixed the coordinate transformation for Campanus house_pos.
  - Previous tolerance: 0.7° - now matches Swiss Ephemeris precisely (0.02° tolerance)
  - Fixed the meridian distance sign convention: uses `mdd - 90` (code's mdd is `ra - armc`)
  - Implemented SE's cotrans rotation formula correctly

- [x] **Methods 0-1 for Gauquelin**: Fixed the hour-angle based methods (0 and 1) to match Swiss Ephemeris.
  - Now uses `house_pos(armc, lat, eps, (lon, lat), "G")` which correctly implements Gauquelin sector calculation
  - Maximum difference vs Swiss Ephemeris: ~0.02 sectors (essentially exact match)
  - Method 0: with ecliptic latitude
  - Method 1: without ecliptic latitude (projected to ecliptic)

## Code Quality & Maintenance

- [ ] **Remaining LSP Errors**: Some pre-existing type errors remain in `houses.py` related to `Union` type handling (e.g., lines 4108, 4118, 4270-4311). These are related to the `lon_or_hsys` parameter type complexity.

