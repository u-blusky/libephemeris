# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.8] - 2025-01-26

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

[Unreleased]: https://github.com/g-battaglia/libephemeris/compare/v0.1.8...HEAD
[0.1.8]: https://github.com/g-battaglia/libephemeris/compare/v0.1.0...v0.1.8
[0.1.0]: https://github.com/g-battaglia/libephemeris/releases/tag/v0.1.0
