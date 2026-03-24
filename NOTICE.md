# NOTICE — Provenance and Intellectual Property

## Independent Implementation

LibEphemeris is an **independent implementation** of an astronomical ephemeris
library for Python. It does not contain, and has never shipped in any release,
code derived from the Swiss Ephemeris (SE) source code by Astrodienst AG.

All astronomical computations are based on:

- **JPL DE440/DE441** planetary and lunar ephemerides, accessed via the
  [Skyfield](https://rhodesmill.org/skyfield/) library (Brandon Rhodes, MIT license)
- **IAU 2006/2000A** precession-nutation model, via
  [pyerfa](https://github.com/liberfa/pyerfa) (BSD-3-Clause license)
- **Peer-reviewed academic sources**, including but not limited to:
  - Meeus, J. — *Astronomical Algorithms*, 2nd ed. (1998)
  - Simon, J.L. et al. — "Numerical expressions for precession formulae
    and mean elements for the Moon and the planets", A&A 282, 663 (1994)
  - Chapront, J. et al. — "A new determination of lunar orbital parameters,
    precession constant and tidal acceleration from LLR", A&A 387, 700 (2002)
  - Chapront-Touze, M. & Chapront, J. — ELP 2000-85 lunar theory
  - Park, R.S. et al. — "The JPL Planetary and Lunar Ephemerides
    DE440 and DE441", AJ 161, 105 (2021)
- **Primary historical sources** for hypothetical bodies:
  - Witte, A. & Sieggrun, F. — *Regelwerk fur Planetenbilder* (1928)
  - Neely, J. — refined orbital elements (1988)
  - Makransky, B. — *Primary Directions* (1988), for the Sunshine house system

## API Compatibility

LibEphemeris provides an API that is **signature-compatible** with
[pyswisseph](https://github.com/astrorigin/pyswisseph) (the Python binding
for Swiss Ephemeris). Function names use the `swe_` prefix and accept the
same parameters and flags to allow drop-in migration.

API compatibility does not imply code derivation. The underlying algorithms,
data sources, and implementation are entirely independent. API signatures
and interface conventions are not copyrightable subject matter
(see *Google LLC v. Oracle America, Inc.*, 593 U.S. 1 (2021)).

## Development History

During early development, some experimental branches temporarily included
data from Swiss Ephemeris sources (e.g., Moshier trigonometric tables,
`seorbel.txt` orbital element file). These were identified, removed, and
replaced with independently sourced alternatives before any stable release:

- Moshier analytical backend: removed entirely in favor of JPL DE440/DE441
  via Skyfield (no analytical approximations are used in production)
- `seorbel.txt`: replaced with `data/fictitious_orbits.csv`, compiled from
  primary published sources (Witte/Sieggrun 1928, Neely 1988)
- All algorithms reference peer-reviewed publications or JPL data products
  as their primary sources

The git history of this repository reflects this progression transparently.

## License

LibEphemeris is distributed under the **GNU Affero General Public License v3**
(AGPL-3.0-only). See the [LICENSE](LICENSE) file for the full text.

This project has no license dependency on, and no license obligation toward,
the Swiss Ephemeris or Astrodienst AG.
