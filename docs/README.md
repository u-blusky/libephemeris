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
