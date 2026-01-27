# LibEphemeris Precision Improvement TODO List

This document contains detailed TODO items for improving libephemeris precision to match Swiss Ephemeris. Each TODO is written as a single comprehensive block with full context for autonomous execution.

---

## CRITICAL PRIORITY: True Lunar Node

- [x] IMPROVE TRUE LUNAR NODE PRECISION MAIN IMPLEMENTATION: Currently the file lunar.py in the function calc_true_lunar_node (lines 337-420) calculates the true lunar node using an osculating geometric method based on the Moon's angular momentum vector (h = r × v) plus some simplified perturbations from Sun and Jupiter defined in the function _calc_planetary_perturbations (lines 121-210), but this approach produces an error of approximately 1-2 degrees compared to Swiss Ephemeris as documented in docs/PRECISION.md lines 242-246 which states "True Node: ~1-2° Different oscillation model", so your task is to completely rewrite the True Node calculation by implementing the complete series of lunar perturbations according to the ELP2000-82B theory by Chapront-Touzé which includes over 60 periodic terms for node perturbations, where the main terms are already partially present in the _calc_planetary_perturbations function but missing are the terms for Venus, Mars, and the long-period terms like evection and variation, and you should also consider using the external library pyerfa (installable with pip install erfa) which provides the standard IAU SOFA functions for nutation and precession with sub-arcsecond precision, or alternatively you can use astropy.coordinates for lunar calculations, and the final result must reduce the error from 1-2 degrees to less than 0.01 degrees (approximately 36 arcseconds) compared to pyswisseph, and for validation you must create tests that compare calc_true_lunar_node with swe.calc_ut(jd, swe.TRUE_NODE, swe.FLG_SPEED) from pyswisseph on at least 100 random dates in the range 1900-2100, and the reference documentation is Swiss Ephemeris swisseph.htm section 2.2.2 "The True Node" which explains how Swiss Ephemeris calculates the true node using the instantaneous position and velocity of the Moon to derive the osculating orbital plane and then applies corrections for planetary perturbations.

- [x] ADD VENUS PERTURBATION TERMS TO TRUE NODE: Currently the function _calc_planetary_perturbations in lunar.py (lines 121-210) only includes perturbation terms from the Sun and Jupiter, but Venus also causes measurable perturbations to the lunar node with amplitudes of several arcminutes, so your task is to add Venus perturbation terms by first calculating Venus mean longitude using a polynomial similar to _calc_jupiter_mean_longitude (lines 96-118), then adding sinusoidal perturbation terms involving Venus longitude combined with the lunar fundamental arguments D, M, M_prime, F, the main Venus-related terms have amplitudes of approximately 0.001-0.005 degrees and should be added to the jupiter_perturbation section or a new venus_perturbation section, reference Chapront-Touzé lunar theory papers for the exact coefficients.

- [x] ADD MARS PERTURBATION TERMS TO TRUE NODE: Similar to Venus, Mars causes small but measurable perturbations to the lunar node that are currently not included in _calc_planetary_perturbations, so your task is to add Mars mean longitude calculation and corresponding perturbation terms with Mars-related arguments, Mars perturbation amplitudes are smaller than Venus (approximately 0.0005-0.002 degrees) but contribute to overall accuracy, add these terms following the same pattern as Jupiter perturbations.

- [x] ADD SATURN PERTURBATION TERMS TO TRUE NODE: Saturn also contributes to lunar node perturbations though with smaller amplitude than Jupiter, so your task is to add Saturn mean longitude and perturbation terms to _calc_planetary_perturbations, Saturn terms have amplitudes of approximately 0.001-0.003 degrees, this completes the set of major planetary perturbations affecting the lunar node.

- [x] ADD LONG PERIOD EVECTION TERM TO TRUE NODE: The evection is a major perturbation of the lunar orbit caused by the Sun's gravitational influence with a period of about 31.8 days and amplitude of 1.274 degrees in longitude, while the current implementation may partially account for this through the solar perturbation terms, verify that the evection effect on the node is fully captured, the evection affects the node through its influence on the lunar eccentricity which in turn affects the node calculation, add explicit evection correction if missing.

- [x] ADD LONG PERIOD VARIATION TERM TO TRUE NODE: The variation is another major lunar perturbation with period of 14.77 days and amplitude of 0.658 degrees, verify that this effect on the node is properly included in the perturbation series, the variation is caused by the difference in solar gravitational force between the Moon's position at quadrature versus conjunction/opposition, add explicit variation correction to the node calculation if not already present.

- [x] ADD ANNUAL EQUATION TERM TO TRUE NODE: The annual equation is a perturbation with period of one anomalistic year (365.26 days) and amplitude of 0.186 degrees caused by the varying Earth-Sun distance due to Earth's orbital eccentricity, verify this long-period term is included in the True Node perturbation series and add if missing.

- [x] ADD PARALLACTIC EQUATION TERM TO TRUE NODE: The parallactic equation has amplitude of about 0.035 degrees and is caused by the finite distance of the Sun, add this term to the perturbation series if not already included.

- [x] IMPLEMENT SECOND ORDER PERTURBATION TERMS FOR TRUE NODE: The current implementation only includes first-order perturbation terms (single sinusoids), but second-order terms involving products of perturbations can contribute at the arcsecond level, implement key second-order terms from the ELP theory to further improve precision.

- [x] IMPLEMENT PROPER FRAME ROTATION FOR TRUE NODE: Currently the True Node calculation transforms the angular momentum vector from ICRS to ecliptic using only mean obliquity (lines 398-406 in lunar.py), but for full precision this should include nutation correction to get the true ecliptic of date, add nutation terms using the IAU 2000A or 2000B model (already available via Skyfield's iau2000a_radians) to rotate to the true ecliptic frame.

- [x] CREATE TRUE NODE PERTURBATION TERM TABLE: Create a comprehensive table or data structure containing all the perturbation coefficients for the True Node based on ELP2000-82B theory, this table should include at least 30-50 terms with their amplitudes, arguments (combinations of D, M, M', F, and planetary longitudes), and phases, this will make the code more maintainable and allow easy verification against published values.

- [x] VALIDATE TRUE NODE AGAINST PYSWISSEPH WITH 1000 RANDOM DATES: Create a comprehensive test file tests/test_precision/test_true_node_precision.py that compares calc_true_lunar_node against pyswisseph swe.calc_ut(jd, swe.TRUE_NODE, flags) for 1000 random dates spanning 1900-2100, calculate the maximum difference, mean difference, and standard deviation, assert that maximum difference is less than 0.01 degrees (36 arcseconds), generate a report of dates with largest discrepancies for debugging.

- [x] VALIDATE TRUE NODE AT HISTORICAL DATES: Create tests comparing True Node precision at historically significant dates like eclipses where the node position is important, test dates from -1000 to 1900 to verify the polynomial validity warnings in lunar.py are appropriate and the precision degradation is as expected.

- [x] VALIDATE TRUE NODE VELOCITY CALCULATION: Currently the True Node returns only position (longitude, 0, 0) but Swiss Ephemeris also returns velocity, verify that when SEFLG_SPEED flag is used the velocity is calculated correctly via numerical differentiation, compare velocity values against pyswisseph and ensure they agree within 0.001 degrees/day.

- [x] DOCUMENT TRUE NODE CALCULATION METHOD: Add comprehensive docstring and comments documenting the True Node calculation method, the perturbation terms included, the expected precision, and references to the source algorithms.

---

## CRITICAL PRIORITY: True Lilith / Osculating Apogee

- [x] IMPROVE TRUE LILITH OSCULATING APOGEE PRECISION MAIN IMPLEMENTATION: Currently the file lunar.py in the function calc_true_lilith (lines 498-597) calculates the osculating lunar apogee using the eccentricity vector method where e = (v × h)/μ - r/|r| to find the direction of perigee and then adds 180 degrees to get apogee, but this simple two-body approach produces an error of approximately 5-7 degrees compared to Swiss Ephemeris as documented in docs/PRECISION.md lines 244 which states "True Lilith: ~5-7° Different orbital model", so your task is to implement a much more sophisticated orbital model that accounts for solar perturbations on the lunar eccentricity.

- [x] ADD EVECTION CORRECTION TO TRUE LILITH: The evection is the largest perturbation to the lunar eccentricity caused by the Sun with amplitude of about 1.274 degrees and period of 31.8 days, the current calc_true_lilith completely ignores this effect, implement the evection correction by adding the term 1.2739 * sin(2*D - l') where D is the mean elongation and l' is the Moon's mean anomaly, this single correction should reduce the error by about 1-2 degrees.

- [x] ADD VARIATION CORRECTION TO TRUE LILITH: The variation effect with amplitude 0.658 degrees and period 14.77 days also affects the apogee position, add the variation term 0.6583 * sin(2*D) to the True Lilith calculation.

- [x] ADD ANNUAL EQUATION TO TRUE LILITH: The annual equation with amplitude 0.186 degrees affects the lunar apogee, add this correction term 0.1856 * sin(M) based on the Sun's mean anomaly M.

- [x] ADD PARALLACTIC INEQUALITY TO TRUE LILITH: Add the parallactic inequality term with amplitude approximately 0.125 degrees.

- [x] ADD REDUCTION TO ECLIPTIC TERM TO TRUE LILITH: Add the reduction to ecliptic term which affects the longitude through the inclination of the lunar orbit.

- [x] ADD EVECTION-RELATED SECONDARY TERMS TO LILITH: Beyond the main evection term, there are several related terms that affect the lunar eccentricity and thus the apogee direction, implement the terms from Meeus "Astronomical Algorithms" Chapter 47 Table 47.B including terms with arguments l-2D, l+2D, 2l, 2l-2D.

- [x] IMPLEMENT PROPER MU VALUE FOR TRUE LILITH: Currently calc_true_lilith uses a fixed GM value for Earth (line 568: mu = 398600.435436 converted to AU³/day²), but for the osculating apogee calculation the effective mu should account for the Earth-Moon mass ratio, verify this is correct by comparing the effective gravitational parameter used in the calculation against the IAU value.

- [x] ADD SOLAR GRAVITATIONAL INFLUENCE ON ECCENTRICITY VECTOR: The eccentricity vector calculation in calc_true_lilith (lines 571-575) assumes pure two-body dynamics, but the Sun's gravity continuously perturbs the lunar eccentricity vector direction, implement a correction for this three-body effect by adding perturbation terms to the eccentricity vector before calculating the apogee direction.

- [x] IMPLEMENT ALTERNATIVE TRUE LILITH METHOD: Research and potentially implement an alternative method for True Lilith based on osculating orbital elements (computing semi-major axis, eccentricity, and argument of perigee from state vectors) rather than the eccentricity vector method, compare both approaches to determine which matches Swiss Ephemeris better.

- [x] VALIDATE TRUE LILITH AGAINST PYSWISSEPH: Create a comprehensive test file tests/test_precision/test_true_lilith_precision.py comparing calc_true_lilith against swe.calc_ut(jd, swe.OSCU_APOG, flags) for 500 random dates, the current error is 5-7 degrees, track progress as corrections are added, target final error less than 0.5 degrees.

- [x] VALIDATE TRUE LILITH LATITUDE: The True Lilith also has an ecliptic latitude component (unlike Mean Lilith which is always on the ecliptic), verify that the latitude calculation in calc_true_lilith (line 592-594) is correct by comparing against pyswisseph.

- [x] VALIDATE TRUE LILITH VELOCITY: When SEFLG_SPEED is set, verify that True Lilith velocity is calculated correctly, the osculating apogee can move rapidly so velocity is important.

- [x] COMPARE TRUE LILITH CALCULATION METHODS: Research and document the exact method Swiss Ephemeris uses for the osculating apogee, it may use a different approach than the eccentricity vector method, document findings and adjust implementation accordingly. COMPLETED: Swiss Ephemeris uses the same fundamental approach (computing osculating orbital elements from Moon's position and speed vectors). The 5-15 degree differences arise from different ephemeris sources, perturbation treatment, and the inherent ambiguity of the osculating apogee concept. Full documentation in docs/TRUE_LILITH_METHODS.md.

---

## CRITICAL PRIORITY: Interpolated Apogee and Perigee

- [x] IMPLEMENT INTERPOLATED LUNAR APOGEE SE_INTP_APOG MAIN FUNCTION: Currently the constants.py file defines SE_INTP_APOG as ID 21 (line 52) and SE_INTP_PERG as ID 22 (line 53) but these are NOT implemented in the library and will likely raise an error or return zero values when requested, so your task is to implement the interpolated apogee calculation as described in Swiss Ephemeris documentation section 2.2.4 "The Interpolated or Natural Apogee and Perigee" which explains that the osculating apogee has spurious oscillations because it is calculated from instantaneous orbital elements that change rapidly due to solar perturbations, and the interpolated apogee smooths out these oscillations by using a polynomial or spline interpolation through successive osculating apogee positions.

- [x] DESIGN INTERPOLATION ALGORITHM FOR APOGEE: Design the interpolation algorithm for the interpolated apogee, this should compute osculating apogee positions at multiple time points around the target date (for example at t-7d, t-3.5d, t, t+3.5d, t+7d spanning approximately half a synodic month) and then fit a smooth curve (polynomial, cubic spline, or similar) through these points to find the interpolated value at time t, research what method Swiss Ephemeris uses.

- [x] IMPLEMENT CALC_INTERPOLATED_APOGEE FUNCTION: Create a new function calc_interpolated_apogee(jd_tt) in lunar.py that computes the interpolated apogee following the designed algorithm, this function should call calc_true_lilith at multiple sample points, extract the longitude values, perform the interpolation, and return the smoothed longitude.

- [x] IMPLEMENT INTERPOLATED LUNAR PERIGEE SE_INTP_PERG: Using the same approach as the interpolated apogee, implement calc_interpolated_perigee(jd_tt) that computes the smoothed perigee position by interpolating through multiple osculating perigee positions, the perigee can be computed as the osculating apogee minus 180 degrees but the interpolation should be done on the perigee values directly for consistency.

- [x] INTEGRATE INTERPOLATED APOGEE INTO PLANETS.PY: Modify the _calc_body function or add a new branch in planets.py to recognize SE_INTP_APOG (21) and SE_INTP_PERG (22) body IDs and route them to the new calc_interpolated_apogee and calc_interpolated_perigee functions in lunar.py, ensure proper handling of flags including SEFLG_SPEED.

- [x] DETERMINE OPTIMAL INTERPOLATION WINDOW FOR APOGEE: Research what time window and how many sample points Swiss Ephemeris uses for the interpolated apogee, test different configurations (5 points over 14 days, 7 points over 21 days, 9 points over 28 days, etc.) and compare results against pyswisseph to find the optimal settings that minimize difference.

- [x] IMPLEMENT VELOCITY FOR INTERPOLATED APOGEE: The interpolated apogee should also return velocity when SEFLG_SPEED is set, implement velocity calculation by analytically differentiating the interpolation polynomial or by numerical differentiation using positions at t-dt and t+dt. COMPLETED: Implemented via numerical differentiation in planets.py:576-596 using forward difference with dt=1 second. Velocity is calculated for longitude, latitude, and distance when SEFLG_SPEED flag is set. Tests exist in tests/test_lunar/test_interpolated_apogee.py (TestInterpolatedApogeeVelocity) and tests/test_lunar/test_interpolated_perigee.py (TestInterpolatedPerigeeVelocity).

- [x] HANDLE EDGE CASES FOR INTERPOLATED APOGEE: Handle edge cases in the interpolated apogee calculation including wrapping around 360/0 degrees (the apogee might cross this boundary during the interpolation window), handle dates at the boundaries of the ephemeris range where future samples cannot be computed.

- [x] VALIDATE INTERPOLATED APOGEE AGAINST PYSWISSEPH: Create tests comparing calc_interpolated_apogee and calc_interpolated_perigee against pyswisseph swe.calc_ut(jd, swe.INTP_APOG, flags) and swe.calc_ut(jd, swe.INTP_PERG, flags) on at least 100 test dates spanning different lunar phases, target error less than 0.1 degrees.

- [x] DOCUMENT INTERPOLATED APOGEE IMPLEMENTATION: Add comprehensive documentation explaining the interpolated apogee/perigee, how it differs from the osculating (true) apogee, and when to use each variant. COMPLETED: Comprehensive documentation added in docs/INTERPOLATED_APOGEE.md covering the problem with osculating apsides, the interpolated solution, implementation details, when to use each variant, API usage examples, and precision limitations. Documentation tests added in tests/test_lunar/test_interpolated_apogee_documentation.py. PRECISION.md updated to include interpolated apogee/perigee entries.

---

## MEDIUM PRIORITY: Minor Bodies - Core Infrastructure

- [x] CREATE SPK_AUTO.PY MODULE FOR AUTOMATIC DOWNLOADS: Create a new module libephemeris/spk_auto.py that provides automatic SPK download and caching functionality for minor bodies, this module should work with the astroquery.jplhorizons library (pip install astroquery) to download SPK files from JPL Horizons on demand.

- [x] IMPLEMENT AUTO_GET_SPK FUNCTION: In spk_auto.py implement the main function auto_get_spk(body_id, jd_start, jd_end) that checks if an SPK file for the requested body and date range exists in a local cache directory (default ~/.libephemeris/spk/), and if not automatically downloads it from JPL Horizons using the Horizons class with ephemerides() method requesting SPK binary output.

- [x] IMPLEMENT SPK CACHE DIRECTORY MANAGEMENT: Create functions to manage the SPK cache directory including ensure_cache_dir() to create the cache directory if it doesn't exist, get_cache_path(body_id) to return the path where a body's SPK should be stored, and cache_info() to return information about cached files.

- [x] IMPLEMENT SPK CACHE LOOKUP: Implement function is_spk_cached(body_id, jd_start, jd_end) that checks if a cached SPK file exists and covers the requested date range, this requires parsing the SPK file header or maintaining a separate index.

- [x] IMPLEMENT JPL HORIZONS SPK DOWNLOAD: Implement function download_spk_from_horizons(body_id, jd_start, jd_end, output_path) that uses astroquery.jplhorizons to download an SPK file for the specified body and date range, handle errors like body not found, date range too large, or network failures gracefully.

- [x] IMPLEMENT SPK REGISTRATION AFTER DOWNLOAD: After downloading an SPK file, automatically register it with the existing spk.py infrastructure so it can be used by calc_ut, this should call the appropriate registration function from spk.py.

- [x] IMPLEMENT AUTOMATIC SPK FALLBACK IN PLANETS.PY: Modify planets.py to automatically try SPK lookup before falling back to Keplerian propagation for bodies in MINOR_BODY_ELEMENTS, add a configuration option (environment variable or function call) to enable/disable automatic SPK download.

- [x] IMPLEMENT SPK CACHE MANAGEMENT FUNCTIONS: Implement helper functions list_cached_spk() to show all cached files with their date ranges and sizes, clear_spk_cache() to delete all cached files, get_cache_size() to report total disk usage, and prune_old_cache(max_age_days) to remove files not accessed recently.

- [x] IMPLEMENT SPK BODY NAME MAPPING: Create a mapping from libephemeris body IDs (SE_CHIRON, SE_ERIS, etc.) to JPL Horizons target designations (2060, 136199, etc.) for use in SPK downloads.

- [x] ADD CONFIGURATION FOR SPK BEHAVIOR: Add configuration options for SPK behavior including auto_download_spk (bool) to enable/disable automatic downloads, spk_cache_dir (path) to override default cache location, and spk_date_padding (days) to add extra time buffer when downloading.

- [x] CREATE SPK DOWNLOAD SCRIPT: Create a command-line script scripts/download_spk.py that can be used to pre-download SPK files for common bodies, useful for users who want to populate the cache ahead of time rather than downloading on demand.

- [x] UPDATE MINOR BODY ORBITAL ELEMENTS TO LATEST EPOCH: Currently the file minor_bodies.py defines orbital elements at epoch 2025.0 (JD 2461000.5) as stated in line 334 and the docstring notes that accuracy degrades 10-50 arcsec per year from epoch, so your task is to update all orbital elements in the MINOR_BODY_ELEMENTS dictionary (lines 340-483) to the most recent epoch available from JPL Small Body Database (sbdb.jpl.nasa.gov).

- [x] CREATE ORBITAL ELEMENTS UPDATE SCRIPT: Create a maintenance script scripts/update_orbital_elements.py that fetches current orbital elements from JPL SBDB API for all bodies in MINOR_BODY_ELEMENTS, compares with current values, and either updates minor_bodies.py directly or generates a report of needed changes, include instructions to run this quarterly.

- [x] ADD URANUS PERTURBATIONS FOR TNO ACCURACY: Currently the file minor_bodies.py in the function calc_secular_perturbation_rates (lines 165-284) only includes perturbations from Jupiter and Saturn, but for Trans-Neptunian Objects the gravitational influence of Uranus is significant, add Uranus orbital elements (a=19.2184 AU, e=0.0457, i=0.772°, Omega=74.0°, omega=170.9°, n=0.01177°/day) and mass ratio (1/22902) to the perturbation calculation.

- [x] ADD NEPTUNE PERTURBATIONS FOR TNO ACCURACY: Add Neptune as a perturbation source with orbital elements (a=30.1104 AU, e=0.0086, i=1.769°, Omega=131.7°, omega=44.9°, n=0.006021°/day) and mass ratio (1/19412) for bodies with semi-major axis greater than 20 AU, this is especially important for plutinos like Ixion and Orcus which are in 2:3 mean motion resonance with Neptune.

- [x] IMPLEMENT MEAN MOTION RESONANCE DETECTION: For bodies in mean motion resonance with Neptune (plutinos 2:3, twotinos 1:2, etc.), the secular perturbation theory breaks down and can give incorrect results, implement detection of resonant orbits by checking if the mean motion is close to a resonant ratio, and flag these bodies or implement special handling.

- [x] VALIDATE MINOR BODY IMPROVEMENTS: After adding Uranus and Neptune perturbations, validate by comparing Eris, Makemake, Ixion, and Orcus positions against pyswisseph over a 50-year span (2000-2050), document the improvement in precision.

---

## MEDIUM PRIORITY: Minor Bodies - Adding New Bodies

- [x] ADD NESSUS CENTAUR TO CATALOG: Add asteroid 7066 Nessus to MINOR_BODY_ELEMENTS with orbital elements from JPL SBDB (a≈24.6 AU, e≈0.52, i≈15.6°, get current epoch values), Nessus is astrologically important as a centaur associated with psychological patterns and abuse cycles, add SE_NESSUS = 7066 + SE_AST_OFFSET to constants.py, add NAIF_NESSUS = 2007066 to constants.py for SPK support.

- [x] ADD ASBOLUS CENTAUR TO CATALOG: Add asteroid 8405 Asbolus to MINOR_BODY_ELEMENTS with current orbital elements from JPL SBDB, Asbolus is another significant centaur in astrological practice, add SE_ASBOLUS = 8405 + SE_AST_OFFSET and NAIF_ASBOLUS = 2008405 to constants.py.

- [x] ADD CHARIKLO CENTAUR TO CATALOG: Add asteroid 10199 Chariklo to MINOR_BODY_ELEMENTS, Chariklo is the largest known centaur (diameter ~250 km) and has a ring system discovered in 2014, add SE_CHARIKLO = 10199 + SE_AST_OFFSET to constants.py.

- [x] ADD GONGGONG TNO TO CATALOG: Add dwarf planet candidate 225088 Gonggong (formerly 2007 OR10) to MINOR_BODY_ELEMENTS with orbital elements from JPL SBDB (a≈67 AU, e≈0.50, i≈30.7°), this is one of the largest TNOs, add SE_GONGGONG = 225088 + SE_AST_OFFSET to constants.py.

- [x] ADD VARUNA TNO TO CATALOG: Add asteroid 20000 Varuna to MINOR_BODY_ELEMENTS, Varuna is a large classical Kuiper belt object with diameter ~670 km, the constant SE_VARUNA is already defined in constants.py line 58 but verify orbital elements are in MINOR_BODY_ELEMENTS.

- [x] ADD APOPHIS NEAR EARTH ASTEROID TO CATALOG: Add asteroid 99942 Apophis to MINOR_BODY_ELEMENTS with very current orbital elements as this asteroid's orbit is being actively refined due to close Earth approaches in 2029 and 2036, add SE_APOPHIS = 99942 + SE_AST_OFFSET to constants.py, note that Apophis orbital elements change measurably with each update.

- [x] ADD HYGIEA ASTEROID TO CATALOG: Add asteroid 10 Hygiea to MINOR_BODY_ELEMENTS, Hygiea is the fourth largest asteroid in the main belt and a dwarf planet candidate, add SE_HYGIEA = 10 + SE_AST_OFFSET to constants.py.

- [x] ADD INTERAMNIA ASTEROID TO CATALOG: Add asteroid 704 Interamnia to MINOR_BODY_ELEMENTS, one of the largest main belt asteroids, add SE_INTERAMNIA = 704 + SE_AST_OFFSET.

- [ ] ADD DAVIDA ASTEROID TO CATALOG: Add asteroid 511 Davida to MINOR_BODY_ELEMENTS, another large main belt asteroid, add SE_DAVIDA = 511 + SE_AST_OFFSET.

- [ ] ADD EUROPA ASTEROID TO CATALOG: Add asteroid 52 Europa (not to be confused with Jupiter's moon) to MINOR_BODY_ELEMENTS, add SE_EUROPA_AST = 52 + SE_AST_OFFSET.

- [ ] ADD SYLVIA ASTEROID TO CATALOG: Add asteroid 87 Sylvia to MINOR_BODY_ELEMENTS, Sylvia has two small moons (Romulus and Remus) making it a triple asteroid system, add SE_SYLVIA = 87 + SE_AST_OFFSET.

- [ ] ADD PSYCHE ASTEROID TO CATALOG: Add asteroid 16 Psyche to MINOR_BODY_ELEMENTS, Psyche is a metallic M-type asteroid and the target of NASA's Psyche mission launching 2023, add SE_PSYCHE = 16 + SE_AST_OFFSET.

- [ ] ADD EROS ASTEROID TO CATALOG: Add asteroid 433 Eros to MINOR_BODY_ELEMENTS, Eros is a near-Earth asteroid visited by the NEAR Shoemaker spacecraft and of astrological interest representing passion, add SE_EROS = 433 + SE_AST_OFFSET.

- [ ] ADD AMOR ASTEROID TO CATALOG: Add asteroid 1221 Amor to MINOR_BODY_ELEMENTS, the prototype of the Amor near-Earth asteroid class, add SE_AMOR = 1221 + SE_AST_OFFSET.

- [ ] ADD ICARUS ASTEROID TO CATALOG: Add asteroid 1566 Icarus to MINOR_BODY_ELEMENTS, an Apollo asteroid with very eccentric orbit that comes closer to the Sun than Mercury, add SE_ICARUS = 1566 + SE_AST_OFFSET.

- [ ] ADD TORO ASTEROID TO CATALOG: Add asteroid 1685 Toro to MINOR_BODY_ELEMENTS, add SE_TORO = 1685 + SE_AST_OFFSET.

- [ ] ADD SAPPHO ASTEROID TO CATALOG: Add asteroid 80 Sappho to MINOR_BODY_ELEMENTS, used in some astrological traditions to represent artistic expression and same-sex love, add SE_SAPPHO = 80 + SE_AST_OFFSET.

- [ ] ADD PANDORA ASTEROID TO CATALOG: Add asteroid 55 Pandora to MINOR_BODY_ELEMENTS, add SE_PANDORA_AST = 55 + SE_AST_OFFSET (distinct from Saturn moon Pandora).

- [ ] ADD LILITH ASTEROID TO CATALOG: Add asteroid 1181 Lilith (the asteroid, not to be confused with lunar apogee Lilith points) to MINOR_BODY_ELEMENTS, add SE_LILITH_AST = 1181 + SE_AST_OFFSET.

- [ ] ADD HIDALGO ASTEROID TO CATALOG: Add asteroid 944 Hidalgo to MINOR_BODY_ELEMENTS, Hidalgo has a comet-like orbit and is used in astrological research, add SE_HIDALGO = 944 + SE_AST_OFFSET.

- [ ] ADD TOUTATIS ASTEROID TO CATALOG: Add asteroid 4179 Toutatis to MINOR_BODY_ELEMENTS, a potentially hazardous asteroid studied by radar and spacecraft, add SE_TOUTATIS = 4179 + SE_AST_OFFSET.

- [ ] ADD ITOKAWA ASTEROID TO CATALOG: Add asteroid 25143 Itokawa to MINOR_BODY_ELEMENTS, target of the Hayabusa sample return mission, add SE_ITOKAWA = 25143 + SE_AST_OFFSET.

- [ ] ADD BENNU ASTEROID TO CATALOG: Add asteroid 101955 Bennu to MINOR_BODY_ELEMENTS, target of OSIRIS-REx sample return mission, add SE_BENNU = 101955 + SE_AST_OFFSET.

- [ ] ADD RYUGU ASTEROID TO CATALOG: Add asteroid 162173 Ryugu to MINOR_BODY_ELEMENTS, target of Hayabusa2 sample return mission, add SE_RYUGU = 162173 + SE_AST_OFFSET.

- [ ] IMPLEMENT GENERIC ASTEROID LOOKUP BY NUMBER: Implement a function calc_asteroid_by_number(asteroid_number, jd_tt) that can calculate position for any numbered asteroid by fetching orbital elements from JPL SBDB API on demand or by checking if an SPK file is available, this allows users to request any of the 1+ million known asteroids without needing each one hardcoded.

- [ ] ADD ASTEROID NAME LOOKUP: Implement function get_asteroid_number(name) that looks up an asteroid's catalog number by name using a local database or JPL SBDB query, useful for users who know the name but not the number.

---

## MEDIUM PRIORITY: Hypothetical Bodies - Uranian Planets

- [ ] CREATE HYPOTHETICAL.PY MODULE: Create a new module libephemeris/hypothetical.py to contain calculation functions for all hypothetical bodies including Uranian planets, Transpluto, and other fictitious bodies, this separates hypothetical calculations from real astronomical bodies.

- [ ] IMPLEMENT CUPIDO URANIAN PLANET: Implement the first Hamburg School Uranian planet Cupido with orbital elements from Swiss Ephemeris seorbel.txt or documentation (a≈40.99837 AU, e≈0.00, nearly circular orbit), add SE_CUPIDO = 40 to constants.py using SE_FICT_OFFSET, implement calc_cupido(jd_tt) in hypothetical.py using simple Keplerian propagation.

- [ ] IMPLEMENT HADES URANIAN PLANET: Implement the second Uranian planet Hades with its orbital elements from seorbel.txt, add SE_HADES = 41 to constants.py, implement calc_hades(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT ZEUS URANIAN PLANET: Implement the third Uranian planet Zeus, add SE_ZEUS = 42 to constants.py, implement calc_zeus(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT KRONOS URANIAN PLANET: Implement the fourth Uranian planet Kronos, add SE_KRONOS = 43 to constants.py, implement calc_kronos(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT APOLLON URANIAN PLANET: Implement the fifth Uranian planet Apollon, add SE_APOLLON = 44 to constants.py, implement calc_apollon(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT ADMETOS URANIAN PLANET: Implement the sixth Uranian planet Admetos, add SE_ADMETOS = 45 to constants.py, implement calc_admetos(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT VULKANUS URANIAN PLANET: Implement the seventh Uranian planet Vulkanus, add SE_VULKANUS = 46 to constants.py, implement calc_vulkanus(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT POSEIDON URANIAN PLANET: Implement the eighth Uranian planet Poseidon, add SE_POSEIDON = 47 to constants.py, implement calc_poseidon(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT GENERIC URANIAN PLANET CALCULATOR: Create a generic function calc_uranian_planet(body_id, jd_tt) that handles all Uranian planets by looking up their orbital elements from a data structure and performing the Keplerian propagation.

- [ ] INTEGRATE URANIAN PLANETS INTO PLANETS.PY: Modify planets.py to recognize Uranian planet body IDs (40-47) and route them to hypothetical.py for calculation.

- [ ] VALIDATE ALL URANIAN PLANETS AGAINST PYSWISSEPH: Create tests comparing all eight Uranian planet positions against pyswisseph swe.calc_ut(jd, swe.CUPIDO, flags) etc. for 50 random dates each, these should match exactly since they use simple orbital propagation from the same elements.

---

## MEDIUM PRIORITY: Hypothetical Bodies - Other Fictitious Bodies

- [ ] IMPLEMENT TRANSPLUTO ISIS: Add support for Transpluto also known as Isis which is a hypothetical planet beyond Pluto proposed by astrologer Ram that Swiss Ephemeris supports as documented in section 2.7.2, get orbital elements from seorbel.txt, add SE_TRANSPLUTO = 48 (or appropriate offset) to constants.py, implement calc_transpluto(jd_tt) in hypothetical.py using Keplerian propagation, validate against pyswisseph.

- [ ] IMPLEMENT VULCAN HYPOTHETICAL PLANET: Add support for the hypothetical intramercurial planet Vulcan as documented in Swiss Ephemeris section 2.7.5, various astrologers have proposed different orbital elements for Vulcan, implement the version Swiss Ephemeris uses from seorbel.txt, add SE_VULCAN to constants.py, implement in hypothetical.py.

- [ ] IMPLEMENT WALDEMATH BLACK MOON: Add support for Dr. Waldemath's hypothetical second moon of Earth also called the Waldemath Moon or Dark Moon as documented in Swiss Ephemeris section 2.7.7, this is different from Mean Lilith and True Lilith which are lunar apogee points, get orbital elements from seorbel.txt, add SE_WALDEMATH to constants.py, implement calc_waldemath(jd_tt) in hypothetical.py.

- [ ] IMPLEMENT WHITE MOON SELENA: Add support for White Moon also called Selena which in most systems is defined as the lunar perigee opposite to Black Moon Lilith (Mean Lilith + 180° or True Lilith + 180°), add SE_WHITE_MOON to constants.py, implement in hypothetical.py, verify which definition Swiss Ephemeris uses (mean or true based).

- [ ] IMPLEMENT PROSERPINA HYPOTHETICAL: Add the hypothetical planet Proserpina as used by some astrologers, get elements from seorbel.txt if present, add SE_PROSERPINA to constants.py, implement in hypothetical.py.

- [ ] IMPLEMENT PLANET X LEVERRIER: Add the hypothetical Planet X as calculated by Leverrier (the one that led to Neptune's discovery), as documented in Swiss Ephemeris section 2.7.8, add SE_PLANET_X_LEVERRIER to constants.py.

- [ ] IMPLEMENT PLANET X ADAMS: Add the hypothetical Planet X as calculated by Adams (similar to Leverrier's but independently derived), add SE_PLANET_X_ADAMS to constants.py.

- [ ] IMPLEMENT PLANET X LOWELL: Add the hypothetical Planet X as calculated by Percival Lowell (the one that led to Pluto's discovery though Pluto was too small to be Lowell's Planet X), add SE_PLANET_X_LOWELL to constants.py.

- [ ] IMPLEMENT PLANET X PICKERING: Add the hypothetical Planet X as calculated by Pickering, add SE_PLANET_X_PICKERING to constants.py.

- [ ] CREATE SEORBEL.TXT PARSER: Create a function parse_seorbel(filepath) that can parse the Swiss Ephemeris seorbel.txt file format to extract orbital elements for any fictitious body, this allows users to add custom hypothetical bodies by providing a seorbel.txt format file.

- [ ] DOWNLOAD AND INCLUDE SEORBEL.TXT: Download the seorbel.txt file from Swiss Ephemeris GitHub repository and include it in the libephemeris package data, or extract the needed orbital elements into Python constants.

---

## MEDIUM PRIORITY: Eclipses - Besselian Elements

- [ ] IMPLEMENT BESSELIAN ELEMENT X CALCULATION: Implement calculation of the Besselian x coordinate which is the x-component of the Moon's shadow axis in the fundamental plane (the plane through Earth's center perpendicular to the Moon-Sun line), this requires precise Moon and Sun positions and the shadow cone geometry.

- [ ] IMPLEMENT BESSELIAN ELEMENT Y CALCULATION: Implement calculation of the Besselian y coordinate which is the y-component of the Moon's shadow axis in the fundamental plane.

- [ ] IMPLEMENT BESSELIAN ELEMENT D CALCULATION: Implement calculation of d, the declination of the Moon's shadow axis relative to the fundamental plane.

- [ ] IMPLEMENT BESSELIAN ELEMENT L1 CALCULATION: Implement calculation of l1, the radius of the penumbral shadow cone where it intersects the fundamental plane.

- [ ] IMPLEMENT BESSELIAN ELEMENT L2 CALCULATION: Implement calculation of l2, the radius of the umbral (or antumbral for annular eclipses) shadow cone where it intersects the fundamental plane.

- [ ] IMPLEMENT BESSELIAN ELEMENT MU CALCULATION: Implement calculation of mu, the Greenwich hour angle of the shadow axis.

- [ ] IMPLEMENT BESSELIAN ELEMENT TIME DERIVATIVES: Calculate the hourly changes in all Besselian elements (dx/dt, dy/dt, dd/dt, dl1/dt, dl2/dt, dmu/dt) for interpolation during the eclipse.

- [ ] CREATE BESSELIAN ELEMENTS DATA STRUCTURE: Create a dataclass or named tuple BesselianElements to hold all elements and their derivatives for a given reference time.

- [ ] IMPLEMENT BESSELIAN ELEMENTS INTERPOLATION: Implement function to interpolate Besselian elements to any time during the eclipse using the elements and derivatives at the reference time.

---

## MEDIUM PRIORITY: Eclipses - Timing and Circumstances

- [ ] IMPROVE ECLIPSE TIMING PRECISION OVERALL: Currently the file eclipse.py calculates solar and lunar eclipses with timing precision of approximately 1-2 minutes as noted in docs/PRECISION.md lines 305-307, but Swiss Ephemeris achieves precision of seconds, the task is to reduce timing error to less than 10 seconds by implementing proper Besselian element calculations.

- [ ] IMPLEMENT ECLIPSE FIRST CONTACT TIME C1: Using Besselian elements, implement precise calculation of first external contact time C1 (when the Moon's disk first touches the Sun's disk externally).

- [ ] IMPLEMENT ECLIPSE SECOND CONTACT TIME C2: Implement precise calculation of second contact time C2 (when totality or annularity begins, the Moon is fully inside/outside the Sun's disk).

- [ ] IMPLEMENT ECLIPSE MAXIMUM TIME: Implement precise calculation of the time of maximum eclipse when the separation between Moon and Sun centers is minimum.

- [ ] IMPLEMENT ECLIPSE THIRD CONTACT TIME C3: Implement precise calculation of third contact time C3 (when totality or annularity ends).

- [ ] IMPLEMENT ECLIPSE FOURTH CONTACT TIME C4: Implement precise calculation of fourth external contact time C4 (when the Moon's disk completely separates from the Sun's disk).

- [ ] IMPLEMENT LUNAR ECLIPSE PENUMBRAL CONTACTS P1 P4: For lunar eclipses, implement calculation of penumbral contact times P1 (Moon enters penumbra) and P4 (Moon exits penumbra).

- [ ] IMPLEMENT LUNAR ECLIPSE UMBRAL CONTACTS U1 U2 U3 U4: For lunar eclipses, implement calculation of umbral contact times U1 (Moon enters umbra), U2 (totality begins), U3 (totality ends), U4 (Moon exits umbra).

- [ ] IMPLEMENT ECLIPSE DURATION CALCULATION: Calculate the duration of totality or annularity for solar eclipses, and duration of umbral/total phase for lunar eclipses.

- [ ] IMPLEMENT SAROS SERIES NUMBER CALCULATION: Implement function get_saros_number(jd_eclipse) that determines which Saros series (numbered approximately 1-180 for solar, 1-180 for lunar) an eclipse belongs to, the Saros cycle is 6585.32 days, implement by relating to known reference eclipses.

- [ ] IMPLEMENT INEX SERIES NUMBER CALCULATION: Implement function get_inex_number(jd_eclipse) that determines which Inex series an eclipse belongs to, the Inex cycle is 10571.95 days (358 synodic months, 28.945 years).

- [ ] VALIDATE ECLIPSE TIMING AGAINST PYSWISSEPH: Create tests comparing eclipse times against pyswisseph for a list of known eclipses including the total solar eclipses of 2017-08-21, 2024-04-08, and the lunar eclipses of 2018-01-31, 2022-11-08, verify timing agrees within 10 seconds.

---

## MEDIUM PRIORITY: Eclipses - Geography and Paths

- [ ] IMPLEMENT ECLIPSE PATH WIDTH CALCULATION: For central solar eclipses, calculate the width of the path of totality or annularity in kilometers at any point along the central line.

- [ ] IMPLEMENT ECLIPSE CENTRAL LINE COORDINATES: Implement function to calculate the geographic coordinates (latitude, longitude) of points along the central line of a solar eclipse.

- [ ] IMPLEMENT ECLIPSE NORTHERN LIMIT: Calculate the northern limit of the umbral/antumbral shadow path.

- [ ] IMPLEMENT ECLIPSE SOUTHERN LIMIT: Calculate the southern limit of the umbral/antumbral shadow path.

- [ ] IMPLEMENT COMPLETE SOL_ECLIPSE_WHERE: Implement or complete swe_sol_eclipse_where(jd, flags) that returns comprehensive geographic information about where an eclipse is visible including central line coordinates, northern and southern limits, and path width.

- [ ] IMPLEMENT SOL_ECLIPSE_HOW WITH FULL DETAILS: Enhance swe_sol_eclipse_how to return all eclipse circumstances at a given location including times of all contacts, maximum obscuration percentage, position angle of contacts, altitude and azimuth of Sun during eclipse, etc.

- [ ] IMPLEMENT ECLIPSE MAGNITUDE AT LOCATION: Calculate the eclipse magnitude (fraction of solar diameter covered) at any specified geographic location.

- [ ] IMPLEMENT ECLIPSE OBSCURATION AT LOCATION: Calculate the obscuration (fraction of solar disk area covered) at any specified geographic location.

---

## MEDIUM PRIORITY: Eclipses - Lunar and Occultations

- [ ] IMPLEMENT LUNAR ECLIPSE UMBRAL MAGNITUDE: Calculate the umbral magnitude (fraction of Moon's diameter within the umbral shadow) for lunar eclipses.

- [ ] IMPLEMENT LUNAR ECLIPSE PENUMBRAL MAGNITUDE: Calculate the penumbral magnitude for lunar eclipses.

- [ ] IMPLEMENT LUNAR ECLIPSE GAMMA: Calculate the gamma parameter for lunar eclipses (distance of Moon center from shadow axis in Earth radii).

- [ ] IMPLEMENT LUN_OCCULT_WHEN_LOC: Implement swe_lun_occult_when_loc(jd_start, planet, starname, flags, geopos, direction) that finds when an occultation of a planet or star by the Moon is visible from a given location, the calculation requires finding when the Moon's disk overlaps the target as seen from the observer accounting for parallax.

- [ ] IMPLEMENT LUN_OCCULT_WHEN_GLOB: Implement swe_lun_occult_when_glob(jd_start, planet, starname, flags, direction) that finds the next occultation visible anywhere on Earth.

- [ ] IMPLEMENT LUN_OCCULT_WHERE: Implement swe_lun_occult_where(jd, planet, starname, flags) that determines where on Earth an occultation is visible at a given time.

- [ ] IMPLEMENT PLANETARY OCCULTATION SEARCH: Extend occultation functions to handle mutual planetary occultations (e.g., when Jupiter occults a star).

- [ ] IMPLEMENT GRAZING OCCULTATION DETECTION: Detect and flag grazing occultations where the star passes near the lunar limb.

---

## LOW PRIORITY: Heliacal Events

- [ ] CREATE HELIACAL.PY MODULE: Create a new module libephemeris/heliacal.py to contain all heliacal event calculation functions including swe_heliacal_ut, swe_vis_limit_mag, and swe_heliacal_pheno_ut.

- [ ] IMPLEMENT ATMOSPHERIC EXTINCTION MODEL: Implement the atmospheric extinction model needed for heliacal calculations, extinction increases the apparent magnitude of objects near the horizon, use the formula from Schaefer or Green where extinction in magnitudes is approximately 0.28 * sec(z) where z is zenith angle, with adjustments for wavelength and atmospheric conditions.

- [ ] IMPLEMENT TWILIGHT SKY BRIGHTNESS MODEL: Implement model for sky brightness during civil twilight (Sun 0° to -6°), nautical twilight (-6° to -12°), and astronomical twilight (-12° to -18°) as a function of Sun altitude, azimuth relative to target, and atmospheric conditions.

- [ ] IMPLEMENT VISIBILITY THRESHOLD MODEL: Implement the Schaefer visibility model that determines whether an object of given magnitude is visible against a sky of given brightness, accounting for the observer's eye adaptation and experience.

- [ ] IMPLEMENT HELIACAL_UT MAIN FUNCTION: Implement swe_heliacal_ut(jd_start, geopos, datm, dobs, object_name, event_type, hel_flag) as described in Swiss Ephemeris documentation section 5.1, where event_type specifies heliacal rising, heliacal setting, evening first, or morning last.

- [ ] IMPLEMENT VIS_LIMIT_MAG FUNCTION: Implement swe_vis_limit_mag(jd, geopos, datm, dobs, object_name, flags) that calculates the limiting magnitude for visibility of a celestial object given atmospheric conditions and observer parameters.

- [ ] IMPLEMENT HELIACAL_PHENO_UT FUNCTION: Implement swe_heliacal_pheno_ut for computing detailed heliacal phenomena information.

- [ ] IMPLEMENT HELIACAL EVENTS FOR PLANETS: Ensure heliacal calculations work correctly for Mercury, Venus, Mars, Jupiter, Saturn including proper handling of inferior/superior conjunction geometry for inner planets.

- [ ] IMPLEMENT HELIACAL EVENTS FOR STARS: Ensure heliacal calculations work for bright fixed stars, determine visibility based on star magnitude and atmospheric conditions.

- [ ] IMPLEMENT LUNAR CRESCENT VISIBILITY: Implement calculation of lunar crescent visibility for Islamic calendar applications, determining when the new crescent moon first becomes visible after conjunction.

- [ ] IMPLEMENT AKHET RISING: Implement the ancient Egyptian concept of akhet rising (heliacal rising of a star).

- [ ] VALIDATE HELIACAL EVENTS AGAINST PYSWISSEPH: Create tests comparing heliacal event times against pyswisseph for several planets and bright stars.

---

## LOW PRIORITY: Fixed Stars - Catalog Expansion

- [ ] ADD ALL 15 BEHENIAN FIXED STARS: Verify that all 15 Behenian fixed stars from medieval astrology are included in fixed_stars.py: Algol, Alcyone (Pleiades), Aldebaran, Regulus, Alkaid (Polaris), Algorab/Gienah, Spica, Arcturus, Alphecca, Antares, Vega, Deneb Algedi, Fomalhaut, Deneb, Markab, add any missing with proper motion data.

- [ ] ADD ALL FOUR ROYAL STARS: Verify that the four Royal Stars of Persia are included with accurate data: Aldebaran (watcher of the East), Regulus (watcher of the North), Antares (watcher of the West), Fomalhaut (watcher of the South).

- [ ] ADD PLEIADES CLUSTER STARS: Add the visible Pleiades stars: Alcyone (brightest), Asterope, Celaeno, Electra, Maia, Merope, Taygeta, Atlas, Pleione.

- [ ] ADD HYADES CLUSTER STARS: Add significant Hyades stars including Prima Hyadum (Gamma Tauri), Secunda Hyadum (Delta Tauri), Theta Tauri, Epsilon Tauri.

- [ ] ADD ORION CONSTELLATION STARS: Ensure all major Orion stars are included: Betelgeuse, Rigel, Bellatrix, Alnilam, Alnitak, Mintaka, Saiph, Meissa.

- [ ] ADD URSA MAJOR STARS: Add all Big Dipper stars: Dubhe, Merak, Phecda, Megrez, Alioth, Mizar (with Alcor), Alkaid.

- [ ] ADD SOUTHERN CROSS STARS: Add Crux constellation stars: Acrux (Alpha Crucis), Mimosa/Becrux (Beta Crucis), Gacrux (Gamma Crucis), Delta Crucis.

- [ ] ADD CENTAURUS BRIGHT STARS: Add Alpha Centauri (Rigil Kentaurus), Beta Centauri (Hadar), and other bright Centaurus stars.

- [ ] ADD SCORPIUS CONSTELLATION STARS: Ensure Scorpius stars are complete: Antares, Shaula, Sargas, Dschubba, Graffias, Lesath.

- [ ] ADD LEO CONSTELLATION STARS: Ensure Leo stars are complete: Regulus, Denebola, Algieba, Zosma.

- [ ] ADD ZODIACAL CONSTELLATION BRIGHT STARS: Add key bright stars from each zodiacal constellation used in astrological interpretation.

- [ ] EXPAND TO 100 BRIGHTEST STARS: Expand the star catalog to include the 100 brightest stars visible from Earth, using Hipparcos catalog data.

- [ ] IMPORT FULL SEFSTARS.TXT CATALOG: Create a script to parse and import the complete Swiss Ephemeris sefstars.txt file containing 800+ stars, convert to StarCatalogEntry format.

- [ ] ADD RADIAL VELOCITY DATA: Add radial velocity field to StarData dataclass and populate for nearby high proper motion stars where this significantly affects position over centuries.

- [ ] ADD PARALLAX DATA: Add parallax field to StarData for nearby stars where annual parallax is significant (>10 mas).

- [ ] ADD SPECTRAL TYPE DATA: Add spectral type (O/B/A/F/G/K/M class) to StarCatalogEntry for informational purposes.

- [ ] ADD VARIABLE STAR FLAG: Add flag indicating if a star is a known variable (like Algol, Mira) with variability type and period if applicable.

- [ ] IMPLEMENT STAR SEARCH BY BAYER DESIGNATION: Implement search for stars by Bayer designation like "Alpha Leonis", "Beta Persei", "Gamma Virginis" parsing Greek letter names.

- [ ] IMPLEMENT STAR SEARCH BY FLAMSTEED NUMBER: Implement search for stars by Flamsteed number like "32 Leonis", "87 Virginis".

- [ ] IMPLEMENT STAR SEARCH BY HIPPARCOS NUMBER: Implement search for stars by HIP number like "HIP 49669", "HIP 65474".

- [ ] IMPLEMENT STAR SEARCH BY HD NUMBER: Implement search for stars by Henry Draper catalog number like "HD 87901".

- [ ] IMPLEMENT STAR SEARCH BY COMMON NAME: Implement fuzzy search for stars by common name handling alternate spellings (Betelgeuse/Betelgeux, Fomalhaut/Formalhaut) and transliterations.

---

## LOW PRIORITY: Planetary Phenomena

- [ ] VERIFY SWE_PHENO_UT IMPLEMENTATION: Verify the current implementation status of swe_pheno_ut in libephemeris, document what is implemented and what is missing.

- [ ] IMPLEMENT PHASE ANGLE CALCULATION: Ensure phase angle (Sun-planet-Earth angle) is calculated correctly for all planets.

- [ ] IMPLEMENT ILLUMINATED FRACTION CALCULATION: Implement phase (illuminated fraction of disk visible from Earth) calculation for all planets using the phase angle.

- [ ] IMPLEMENT MERCURY MAGNITUDE FORMULA: Implement apparent magnitude calculation for Mercury using the formula from Meeus Chapter 41 which accounts for Mercury's large phase angle range.

- [ ] IMPLEMENT VENUS MAGNITUDE FORMULA: Implement apparent magnitude for Venus with its unique magnitude curve that peaks near dichotomy due to atmospheric effects.

- [ ] IMPLEMENT MARS MAGNITUDE FORMULA: Implement apparent magnitude for Mars including the opposition surge effect.

- [ ] IMPLEMENT JUPITER MAGNITUDE FORMULA: Implement apparent magnitude for Jupiter based on phase angle and distance.

- [ ] IMPLEMENT SATURN MAGNITUDE FORMULA: Implement apparent magnitude for Saturn including the effect of ring tilt angle on total brightness.

- [ ] IMPLEMENT SATURN RING TILT CALCULATION: Calculate the tilt angle of Saturn's rings as seen from Earth, needed for magnitude calculation.

- [ ] IMPLEMENT URANUS NEPTUNE MAGNITUDE FORMULAS: Implement simplified magnitude formulas for Uranus and Neptune.

- [ ] IMPLEMENT PLANET APPARENT DIAMETER CALCULATION: Calculate apparent angular diameter in arcseconds for each planet based on physical radius and geocentric distance.

- [ ] IMPLEMENT DETAILED MOON PHASE FUNCTION: Implement swe_lun_phase(jd) returning phase angle (0-360°), illumination fraction (0.0-1.0), phase name, and lunar age in days.

- [ ] IMPLEMENT MOON PHASE NAME DETERMINATION: Return appropriate phase name (New, Waxing Crescent, First Quarter, Waxing Gibbous, Full, Waning Gibbous, Last Quarter, Waning Crescent) based on phase angle.

- [ ] IMPLEMENT ELONGATION CALCULATION: Ensure elongation from Sun is calculated correctly, properly distinguish between morning star (western elongation) and evening star (eastern elongation).

- [ ] IMPLEMENT RETROGRADE STATION FINDER: Implement function to find exact times when a planet stations retrograde or direct, using root-finding on the velocity.

- [ ] IMPLEMENT CONJUNCTION FINDER: Implement function to find planetary conjunctions with the Sun (superior and inferior) and with other planets.

- [ ] IMPLEMENT OPPOSITION FINDER: Implement function to find when outer planets are at opposition (opposite the Sun in the sky).

- [ ] IMPLEMENT GREATEST ELONGATION FINDER: Implement function to find times of greatest eastern and western elongation for Mercury and Venus.

---

## LOW PRIORITY: House System Improvements

- [ ] IMPLEMENT CO-ASCENDANT CALCULATION: Currently docs/PRECISION.md line 121 states "Co-Ascendant: Not implemented (returns 0.0)", implement the Co-Ascendant in houses.py as defined by Walter Koch (Ascendant calculated at the same sidereal time but at latitude 0° equator).

- [ ] IMPLEMENT POLAR ASCENDANT CALCULATION: Currently docs/PRECISION.md line 122 states "Polar Ascendant: Not implemented (returns 0.0)", implement the Polar Ascendant in houses.py as the Ascendant calculated at the same sidereal time at latitude 90°.

- [ ] VERIFY GAUQUELIN SECTOR IMPLEMENTATION: Currently docs/PRECISION.md line 106 notes that Gauquelin sectors use "Placidus approximation (not true 36-sector)", verify the implementation correctly divides diurnal and nocturnal arcs into 18 sectors each.

- [ ] IMPROVE PLACIDUS POLAR LATITUDE HANDLING: Currently Placidus houses fall back to Porphyry above 66.5° latitude (Arctic Circle), implement better handling with clearer warnings.

- [ ] IMPROVE KOCH POLAR LATITUDE HANDLING: Similar to Placidus, improve Koch house calculation for polar latitudes with appropriate fallbacks.

- [ ] VERIFY ALL 19 HOUSE SYSTEMS: Create comprehensive tests comparing all 19 house systems (Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign, Porphyry, Alcabitius, Topocentric, Morinus, Meridian, Vehlow, Horizontal, Carter, Krusinski, Natural, Gauquelin, APC, Sripati) against pyswisseph at 100+ random locations and times.

- [ ] DOCUMENT HOUSE SYSTEM ALGORITHMS: Document the mathematical algorithm and formula used for each house system in code comments or separate documentation.

- [ ] IMPLEMENT HOUSE CUSP VELOCITY: When SEFLG_SPEED is set, calculate the rate of change of house cusps for progressed chart applications.

---

## LOW PRIORITY: Ayanamsha Improvements

- [ ] VERIFY ALL 43 AYANAMSHA MODES: Create comprehensive tests comparing all 43 ayanamsha modes against pyswisseph at multiple dates (J2000.0, 1900, 2000, 2100) to verify accuracy.

- [ ] IMPROVE TRUE CITRA AYANAMSHA PRECISION: Currently docs/PRECISION.md notes star-based ayanamshas have ±0.06° precision, improve True Citra by using more precise Spica coordinates with full proper motion correction.

- [ ] IMPROVE TRUE REVATI AYANAMSHA PRECISION: Similar improvement for True Revati ayanamsha using precise Zeta Piscium coordinates.

- [ ] IMPROVE TRUE PUSHYA AYANAMSHA PRECISION: Improve True Pushya ayanamsha precision.

- [ ] IMPROVE TRUE MULA AYANAMSHA PRECISION: Improve True Mula ayanamsha precision.

- [ ] IMPROVE GALACTIC CENTER AYANAMSHAS: Improve Galactic Center based ayanamshas (0 Sag, Rgilbrand, Cochrane, Mula Wilhelm) using current best coordinates for Sgr A* (the radio source at the Galactic Center).

- [ ] IMPLEMENT CUSTOM AYANAMSHA: Implement swe_set_sid_mode with SE_SIDM_USER option allowing users to define custom ayanamsha with their own initial value at a reference epoch and annual precession rate.

- [ ] DOCUMENT AYANAMSHA DEFINITIONS: Create documentation explaining the astronomical and historical basis for each ayanamsha mode, what reference point it uses, and when it was zero.

- [ ] VALIDATE LAHIRI AGAINST IAE: Validate the Lahiri ayanamsha against the Indian Astronomical Ephemeris official values as described in Swiss Ephemeris Appendix E.

---

## OPTIONAL: Performance and Architecture

- [ ] BENCHMARK AGAINST PYSWISSEPH: Create comprehensive benchmarks comparing calculation speed for planets, houses, and other functions against pyswisseph, quantify the performance difference.

- [ ] PROFILE HOT PATHS: Use Python profiling tools to identify the most time-consuming functions and optimize them.

- [ ] IMPLEMENT BATCH CALCULATION: Implement batch calculation mode for computing positions at many dates efficiently, leveraging Skyfield's vectorized calculation capability.

- [ ] ADD POSITION CACHING: Implement optional caching of recently computed positions using an LRU cache to speed up repeated queries.

- [ ] EVALUATE DE440 DE441 UPGRADE: Currently libephemeris uses DE421 as default, evaluate upgrading to DE440 (2020) or DE441 which have improved accuracy especially for outer planets and extended time range.

- [ ] ADD DE440 SUPPORT: Add support for DE440 ephemeris with documentation on how to download and configure it.

- [ ] ADD DE441 SUPPORT: Add support for DE441 ephemeris for extended time range calculations.

- [ ] VERIFY THREAD SAFETY: Verify that EphemerisContext provides true thread safety for multi-threaded applications.

---

## OPTIONAL: Planetary Moons

- [ ] IMPLEMENT JUPITER GALILEAN MOON POSITIONS: Implement positions for Io, Europa, Ganymede, Callisto using JPL satellite ephemeris file jup365.bsp.

- [ ] IMPLEMENT SATURN MAJOR MOON POSITIONS: Implement positions for Titan, Rhea, Iapetus, Dione, Tethys, Enceladus using sat441.bsp.

- [ ] IMPLEMENT MARS MOON POSITIONS: Implement positions for Phobos and Deimos using mar097.bsp.

- [ ] IMPLEMENT URANUS MOON POSITIONS: Implement positions for major Uranian moons (Miranda, Ariel, Umbriel, Titania, Oberon) using ura111.bsp.

- [ ] IMPLEMENT NEPTUNE MOON POSITIONS: Implement positions for Triton and Nereid using nep095.bsp.

- [ ] IMPLEMENT PLUTO MOON POSITIONS: Implement positions for Charon, Nix, Hydra, Kerberos, Styx using plu058.bsp.

---

## OPTIONAL: Time and Coordinates

- [ ] VERIFY DELTA T IMPLEMENTATION: Verify Delta T implementation in time_utils.py matches the Stephenson/Morrison/Hohenkerk 2016 model.

- [ ] ADD IERS DELTA T DOWNLOAD: Implement automatic download of observed Delta T values from IERS for recent dates.

- [ ] VERIFY TDB TT HANDLING: Verify that the library correctly handles the distinction between Barycentric Dynamical Time (TDB) and Terrestrial Time (TT) where relevant.

- [ ] IMPLEMENT TAI TIME SCALE: Add support for International Atomic Time if not present.

- [ ] VERIFY UTC LEAP SECOND HANDLING: Verify that UTC to UT1 conversion correctly handles leap seconds using current IERS data.

- [ ] VERIFY COORDINATE TRANSFORMATIONS: Verify all coordinate transformation flags (ICRS, J2000, equinox of date, ecliptic, equatorial) work correctly.

---

## DOCUMENTATION AND TESTING

- [ ] CREATE PRECISION VALIDATION TEST SUITE: Create a new test directory tests/test_precision_validation/ with comprehensive tests comparing every calculation type against pyswisseph.

- [ ] TEST ALL PLANETS 1000 DATES: Create test comparing Sun, Moon, Mercury through Pluto positions at 1000 random dates.

- [ ] TEST ALL HOUSE SYSTEMS 100 LOCATIONS: Create test comparing all 19 house systems at 100 random lat/lon/time combinations.

- [ ] TEST ALL 43 AYANAMSHAS: Create test comparing all ayanamsha modes at multiple dates.

- [ ] TEST MINOR BODIES WITH WITHOUT SPK: Create tests verifying minor body precision with Keplerian fallback vs SPK files.

- [ ] GENERATE PRECISION REPORT: Create script that generates a comprehensive precision report showing max/mean/std deviation for each calculation type.

- [ ] UPDATE PRECISION.MD: After implementing improvements, update docs/PRECISION.md to reflect new accuracy levels.

- [ ] CREATE PRECISION TUNING GUIDE: Create docs/PRECISION_TUNING.md explaining how to achieve maximum precision with optional dependencies.

- [ ] DOCUMENT OPTIONAL DEPENDENCIES: Document what each optional dependency (pyerfa, astropy, astroquery, rebound) provides.

- [ ] CREATE MIGRATION GUIDE: Create documentation helping users migrate from pyswisseph to libephemeris.

- [ ] ADD COMPREHENSIVE DOCSTRINGS: Verify all public API functions have complete docstrings with parameters, returns, and examples.

- [ ] CREATE API REFERENCE: Generate or write complete API reference documentation.

- [ ] ADD USAGE EXAMPLES: Add example scripts demonstrating common use cases.

---

## ERROR HANDLING AND EDGE CASES

- [ ] IMPROVE DATE RANGE ERROR MESSAGES: When calculations fail due to date outside ephemeris range, provide clear error message with supported range.

- [ ] HANDLE MISSING SPK GRACEFULLY: When SPK file is requested but not available, provide helpful error message explaining how to obtain it.

- [ ] HANDLE POLAR EDGE CASES: Improve error handling and fallback behavior for house calculations at extreme latitudes (>80°).

- [ ] VALIDATE INPUT COORDINATES: Add validation for geographic coordinates (lat -90 to 90, lon -180 to 180) with clear error messages.

- [ ] VALIDATE JULIAN DAY RANGE: Add validation that Julian Day is within supported ephemeris range before calculation.

- [ ] HANDLE UNKNOWN BODY IDS: Provide clear error message when unknown body ID is requested.

- [ ] HANDLE RETROGRADE STATIONS: Ensure calculations remain stable when planet is near stationary point with near-zero velocity.

- [ ] HANDLE ECLIPSE EDGE CASES: Handle edge cases in eclipse calculations like very shallow partial eclipses.

- [ ] IMPROVE EXCEPTION HIERARCHY: Review and improve the exception classes in exceptions.py for better error categorization.

---

## DEPENDENCIES TO EVALUATE

- [ ] EVALUATE PYERFA INTEGRATION: Research the pyerfa library (pip install erfa) which provides Python bindings to IAU SOFA/ERFA routines, evaluate precision improvements from using erfa.nut00a(), erfa.pnm06a().

- [ ] DOCUMENT PYERFA BENEFITS: Document what precision improvements pyerfa would provide if integrated.

- [ ] EVALUATE ASTROPY INTEGRATION: Research how astropy.coordinates and astropy.time could supplement calculations.

- [ ] DOCUMENT ASTROPY BENEFITS: Document what astropy provides that could improve libephemeris.

- [ ] EVALUATE REBOUND FOR N-BODY: Research the REBOUND library (pip install rebound) for n-body gravitational integration, evaluate for asteroid orbit propagation.

- [ ] DOCUMENT REBOUND BENEFITS: Document how REBOUND could improve asteroid precision.

- [ ] EVALUATE JPLEPHEM VERSION: Verify using latest jplephem version for optimal SPK reading performance.

- [ ] EVALUATE SKYFIELD VERSION: Verify using latest Skyfield version and document any improvements from updates.

---

This list contains approximately 300 detailed TODO items covering all aspects of improving libephemeris precision to match Swiss Ephemeris.
