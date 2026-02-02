# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/g-battaglia/libephemeris/compare/v0.5.1...HEAD
[0.5.1]: https://github.com/g-battaglia/libephemeris/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/g-battaglia/libephemeris/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/g-battaglia/libephemeris/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/g-battaglia/libephemeris/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/g-battaglia/libephemeris/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/g-battaglia/libephemeris/releases/tag/v0.1.0
