"""
Moshier mode comparison tests for ICRS, J2000, and NONUT frame flags.

Tests compare pyswisseph (C library) vs libephemeris (Python) with SEFLG_MOSEPH
combined with SEFLG_ICRS, SEFLG_J2000, and SEFLG_NONUT flags.

Frame flag handling in Moshier mode:
- The C library (pyswisseph) with Moshier computes VSOP87 positions and then
  applies the requested frame transformations (precession, nutation, ICRS tie).
- libephemeris Moshier path currently outputs mean ecliptic of date from VSOP87
  series, and does NOT handle SEFLG_ICRS, SEFLG_J2000, or SEFLG_NONUT.

Expected behavior differences:
- SEFLG_ICRS: C library applies ICRS frame tie rotation (~50 mas offset from
  dynamical frame). libephemeris ignores the flag entirely -- it is unimplemented
  in all calculation paths (both Moshier and SPK). Differences grow with distance
  from J2000 epoch due to precession effects.
- SEFLG_J2000: C library returns J2000 ecliptic coordinates (no precession,
  no nutation). libephemeris Moshier path returns mean ecliptic of date
  (VSOP87 native frame) regardless of the flag.
- SEFLG_NONUT: C library suppresses nutation correction. libephemeris Moshier
  path never applies nutation, so this flag is accidentally a no-op.

These tests are marked with pytest.mark.xfail where the flags are known to be
unimplemented, documenting the gaps for future implementation.

Total: 85 test cases:
- 3 frame flag classes x 5 planets x 3 dates = 45 tests
- 1 edge case class x 5 planets x 8 extreme dates = 40 tests
"""

from __future__ import annotations

import pytest

swe = pytest.importorskip("swisseph", reason="pyswisseph required for comparison tests")

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_ICRS,
    SEFLG_J2000,
    SEFLG_NONUT,
)

# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# 5 planets covering inner, outer, and the Sun
PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MERCURY, "Mercury"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# 3 dates: J2000 epoch, a modern date, and a historical date
TEST_DATES = [
    (2451545.0, "J2000.0 (2000-01-01)"),
    (2460000.5, "2023-02-25"),
    (2440000.5, "1968-05-24"),
]

# Tolerance for Moshier comparison: 0.01 degrees = 36 arcseconds
# This accounts for differences between C and Python VSOP87 implementations
MOSHIER_TOLERANCE_DEG = 0.01

# Distance tolerance as percentage
DISTANCE_TOLERANCE_PCT = 1.0


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMoshierICRSPositions:
    """Compare SEFLG_MOSEPH | SEFLG_ICRS between pyswisseph and libephemeris.

    ICRS (International Celestial Reference System) is the IAU standard
    reference frame aligned with distant quasars. It differs from the
    dynamical J2000 frame by a small rotation (~50 mas frame tie).

    In the C library, SEFLG_ICRS with Moshier applies the ICRS frame
    tie rotation to the computed positions. In libephemeris, the
    SEFLG_ICRS flag is completely unimplemented -- it is silently
    ignored in all calculation paths (both Moshier and SPK).

    Expected differences: The ICRS frame tie is ~50 mas, but since
    libephemeris also lacks proper precession/nutation handling in
    Moshier mode, the total difference depends on the epoch. At J2000,
    the difference is dominated by the VSOP87 implementation mismatch.
    At other epochs, the frame tie + precession effects accumulate.
    """

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason=(
            "SEFLG_ICRS is not implemented in libephemeris Moshier path. "
            "The flag is silently ignored, causing ~50 mas frame tie "
            "offset plus potential precession differences at non-J2000 epochs."
        ),
        strict=False,
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    @pytest.mark.parametrize("body_id,body_name", PLANETS)
    def test_moshier_icrs_position(self, jd, date_desc, body_id, body_name):
        """MOSEPH|ICRS position should match between C and Python.

        Compares longitude, latitude, and distance for each planet/date
        combination using SEFLG_MOSEPH | SEFLG_ICRS | SEFLG_SPEED.
        """
        flag = SEFLG_MOSEPH | SEFLG_ICRS | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, body_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        if pos_swe[2] > 0:
            diff_dist_pct = abs(pos_swe[2] - pos_py[2]) / pos_swe[2] * 100.0
        else:
            diff_dist_pct = 0.0

        assert diff_lon < MOSHIER_TOLERANCE_DEG, (
            f"{body_name} MOSEPH|ICRS at {date_desc}: "
            f"longitude diff {diff_lon:.6f} deg exceeds "
            f"{MOSHIER_TOLERANCE_DEG} deg tolerance\n"
            f"  pyswisseph: {pos_swe[0]:.6f} deg\n"
            f"  libephemeris: {pos_py[0]:.6f} deg"
        )
        assert diff_lat < MOSHIER_TOLERANCE_DEG, (
            f"{body_name} MOSEPH|ICRS at {date_desc}: "
            f"latitude diff {diff_lat:.6f} deg exceeds "
            f"{MOSHIER_TOLERANCE_DEG} deg tolerance\n"
            f"  pyswisseph: {pos_swe[1]:.6f} deg\n"
            f"  libephemeris: {pos_py[1]:.6f} deg"
        )
        assert diff_dist_pct < DISTANCE_TOLERANCE_PCT, (
            f"{body_name} MOSEPH|ICRS at {date_desc}: "
            f"distance diff {diff_dist_pct:.4f}% exceeds "
            f"{DISTANCE_TOLERANCE_PCT}% tolerance\n"
            f"  pyswisseph: {pos_swe[2]:.8f} AU\n"
            f"  libephemeris: {pos_py[2]:.8f} AU"
        )


class TestMoshierJ2000Positions:
    """Compare SEFLG_MOSEPH | SEFLG_J2000 between pyswisseph and libephemeris.

    SEFLG_J2000 requests positions in the J2000.0 ecliptic reference frame,
    suppressing precession from J2000 to equinox of date. In the C library,
    Moshier mode with SEFLG_J2000 computes VSOP87 positions and skips the
    precession step (since VSOP87 natively produces mean ecliptic coordinates,
    the C library's internal handling is equivalent to skipping the
    precession-to-date that would normally be applied).

    In libephemeris, the Moshier path returns mean ecliptic of date coordinates
    from the VSOP87 series. SEFLG_J2000 is not checked in the Moshier path.
    Since the VSOP87 series used by Meeus produce coordinates in the dynamical
    ecliptic of the date (not J2000.0), the libephemeris output at dates far
    from J2000 will differ from the C library's J2000-frame output by the
    accumulated precession.

    At J2000 epoch itself, precession is zero, so the results may coincidentally
    agree. At distant epochs (e.g., 1968, 2023), differences grow at
    ~0.014 deg/year due to precession.
    """

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason=(
            "SEFLG_J2000 is not handled in libephemeris Moshier path. "
            "The Moshier path returns mean ecliptic of date (VSOP87 native), "
            "while the C library with SEFLG_J2000 returns J2000 ecliptic. "
            "At non-J2000 epochs, the difference equals accumulated precession."
        ),
        strict=False,
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    @pytest.mark.parametrize("body_id,body_name", PLANETS)
    def test_moshier_j2000_position(self, jd, date_desc, body_id, body_name):
        """MOSEPH|J2000 position should match between C and Python.

        Compares longitude, latitude, and distance for each planet/date
        combination using SEFLG_MOSEPH | SEFLG_J2000 | SEFLG_SPEED.
        """
        flag = SEFLG_MOSEPH | SEFLG_J2000 | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, body_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        if pos_swe[2] > 0:
            diff_dist_pct = abs(pos_swe[2] - pos_py[2]) / pos_swe[2] * 100.0
        else:
            diff_dist_pct = 0.0

        assert diff_lon < MOSHIER_TOLERANCE_DEG, (
            f"{body_name} MOSEPH|J2000 at {date_desc}: "
            f"longitude diff {diff_lon:.6f} deg exceeds "
            f"{MOSHIER_TOLERANCE_DEG} deg tolerance\n"
            f"  pyswisseph: {pos_swe[0]:.6f} deg\n"
            f"  libephemeris: {pos_py[0]:.6f} deg"
        )
        assert diff_lat < MOSHIER_TOLERANCE_DEG, (
            f"{body_name} MOSEPH|J2000 at {date_desc}: "
            f"latitude diff {diff_lat:.6f} deg exceeds "
            f"{MOSHIER_TOLERANCE_DEG} deg tolerance\n"
            f"  pyswisseph: {pos_swe[1]:.6f} deg\n"
            f"  libephemeris: {pos_py[1]:.6f} deg"
        )
        assert diff_dist_pct < DISTANCE_TOLERANCE_PCT, (
            f"{body_name} MOSEPH|J2000 at {date_desc}: "
            f"distance diff {diff_dist_pct:.4f}% exceeds "
            f"{DISTANCE_TOLERANCE_PCT}% tolerance\n"
            f"  pyswisseph: {pos_swe[2]:.8f} AU\n"
            f"  libephemeris: {pos_py[2]:.8f} AU"
        )


class TestMoshierNoNutation:
    """Compare SEFLG_MOSEPH | SEFLG_NONUT between pyswisseph and libephemeris.

    SEFLG_NONUT suppresses the nutation correction, yielding mean ecliptic
    positions rather than true (apparent) ecliptic positions. The difference
    between mean and true positions is the nutation in longitude (delta_psi),
    which oscillates with a period of ~18.6 years and amplitude ~17 arcseconds.

    In the C library, Moshier mode with SEFLG_NONUT computes VSOP87 positions,
    applies precession to equinox of date, but skips nutation.

    In libephemeris, the Moshier path returns mean ecliptic of date (VSOP87
    native output) WITHOUT applying either precession or nutation. Since
    SEFLG_NONUT is not explicitly checked, it is silently ignored. However,
    because nutation is already absent from the Moshier output, the NONUT
    flag is effectively a no-op.

    The key difference is that the C library still applies precession when
    SEFLG_NONUT is set (only nutation is skipped), while libephemeris does
    not apply precession in the Moshier path. So the comparison results
    depend on how much the VSOP87 native frame differs from precessed-to-date
    without nutation. Since VSOP87 Meeus series produce ecliptic-of-date
    coordinates natively, the precession component may coincidentally cancel.
    """

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason=(
            "SEFLG_NONUT is not explicitly handled in libephemeris Moshier path. "
            "The C library applies precession but skips nutation; libephemeris "
            "applies neither. The VSOP87 native frame (mean ecliptic of date) "
            "may coincidentally approximate the C library's NONUT output, but "
            "differences depend on epoch and VSOP87 implementation details."
        ),
        strict=False,
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    @pytest.mark.parametrize("body_id,body_name", PLANETS)
    def test_moshier_nonut_position(self, jd, date_desc, body_id, body_name):
        """MOSEPH|NONUT position should match between C and Python.

        Compares longitude, latitude, and distance for each planet/date
        combination using SEFLG_MOSEPH | SEFLG_NONUT | SEFLG_SPEED.
        """
        flag = SEFLG_MOSEPH | SEFLG_NONUT | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, body_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        if pos_swe[2] > 0:
            diff_dist_pct = abs(pos_swe[2] - pos_py[2]) / pos_swe[2] * 100.0
        else:
            diff_dist_pct = 0.0

        assert diff_lon < MOSHIER_TOLERANCE_DEG, (
            f"{body_name} MOSEPH|NONUT at {date_desc}: "
            f"longitude diff {diff_lon:.6f} deg exceeds "
            f"{MOSHIER_TOLERANCE_DEG} deg tolerance\n"
            f"  pyswisseph: {pos_swe[0]:.6f} deg\n"
            f"  libephemeris: {pos_py[0]:.6f} deg"
        )
        assert diff_lat < MOSHIER_TOLERANCE_DEG, (
            f"{body_name} MOSEPH|NONUT at {date_desc}: "
            f"latitude diff {diff_lat:.6f} deg exceeds "
            f"{MOSHIER_TOLERANCE_DEG} deg tolerance\n"
            f"  pyswisseph: {pos_swe[1]:.6f} deg\n"
            f"  libephemeris: {pos_py[1]:.6f} deg"
        )
        assert diff_dist_pct < DISTANCE_TOLERANCE_PCT, (
            f"{body_name} MOSEPH|NONUT at {date_desc}: "
            f"distance diff {diff_dist_pct:.4f}% exceeds "
            f"{DISTANCE_TOLERANCE_PCT}% tolerance\n"
            f"  pyswisseph: {pos_swe[2]:.8f} AU\n"
            f"  libephemeris: {pos_py[2]:.8f} AU"
        )


# ============================================================================
# EDGE CASE CONFIGURATIONS
# ============================================================================

# Dates spanning the full Moshier range for extreme date testing.
# Format: (year, month, day, hour, description)
# Year uses astronomical year numbering: year 0 = 1 BCE, -1 = 2 BCE, etc.
EDGE_CASE_DATES = [
    (-2000, 1, 1, 12.0, "-2000 CE (deep ancient)"),
    (-1000, 6, 15, 12.0, "-1000 CE (ancient)"),
    (0, 1, 1, 12.0, "0 CE / 1 BCE"),
    (500, 3, 21, 12.0, "500 CE"),
    (1000, 12, 25, 12.0, "1000 CE"),
    (1200, 7, 4, 12.0, "1200 CE"),
    (2700, 1, 1, 12.0, "2700 CE (post-DE440)"),
    (3000, 6, 15, 12.0, "3000 CE (Moshier edge)"),
]

# 5 planets for edge case testing: Sun, Moon, and outer planets
EDGE_CASE_PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]


def _tolerance_for_year(year: int, body_id: int) -> float:
    """Return graduated longitude tolerance in degrees based on era and body.

    VSOP87 (planets) and ELP2000-82B (Moon) are series expansions about
    J2000.0 that degrade with distance from that epoch. The C library
    (pyswisseph) and Python (libephemeris) use different truncations of
    these series, so differences between them grow with distance from J2000.

    The Moon uses ELP2000-82B which diverges more than VSOP87 between the
    C and Python implementations, especially at dates far from J2000.
    Empirically, Moon differences reach ~0.12 deg at 3000 CE while planets
    stay below ~0.02 deg at the same epoch.

    Tolerance bands for planets (Sun, Mars, Jupiter, Saturn):
        >= 1000 CE:     0.02 deg (72 arcsec)  - close to J2000
        0 - 1000 CE:    0.1  deg (6 arcmin)   - moderate distance
        -1000 - 0 CE:   0.5  deg (30 arcmin)  - far from J2000
        < -1000 CE:     1.0  deg (1 degree)   - very far from J2000

    Tolerance bands for Moon (ELP2000-82B):
        >= 1000 CE:     0.15 deg (9 arcmin)   - ELP2000-82B divergence
        0 - 1000 CE:    0.5  deg (30 arcmin)  - moderate distance
        -1000 - 0 CE:   1.0  deg (1 degree)   - far from J2000
        < -1000 CE:     2.0  deg (2 degrees)  - very far from J2000
    """
    if body_id == SE_MOON:
        if year >= 1000:
            return 0.15
        elif year >= 0:
            return 0.5
        elif year >= -1000:
            return 1.0
        else:
            return 2.0

    # Planets: VSOP87 series
    if year >= 1000:
        return 0.02
    elif year >= 0:
        return 0.1
    elif year >= -1000:
        return 0.5
    else:
        return 1.0


class TestMoshierEdgeCases:
    """Test Moshier mode at extreme dates spanning -2000 to +3000 CE.

    This class verifies positional continuity and accuracy across the full
    Moshier range, comparing pyswisseph (C) vs libephemeris (Python) with
    SEFLG_MOSEPH | SEFLG_SPEED.

    The primary use case for Moshier is dates outside DE440 range (1550-2650 CE):
    archaeological astronomy, astro-chronology, and ancient civilization studies.
    VSOP87/ELP2000-82B accuracy degrades with distance from J2000.0, so tolerances
    are graduated by era to detect regressions without false positives.

    Total: 40 test cases (8 dates x 5 planets).
    """

    @pytest.mark.comparison
    @pytest.mark.edge_case
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc",
        EDGE_CASE_DATES,
        ids=[d[4] for d in EDGE_CASE_DATES],
    )
    @pytest.mark.parametrize(
        "body_id,body_name",
        EDGE_CASE_PLANETS,
        ids=[p[1] for p in EDGE_CASE_PLANETS],
    )
    def test_moshier_extreme_date(
        self, year, month, day, hour, date_desc, body_id, body_name
    ):
        """MOSEPH position at extreme date should match between C and Python.

        Compares longitude and latitude for each planet/date combination
        using SEFLG_MOSEPH | SEFLG_SPEED with graduated tolerances based
        on distance from J2000. Prints degradation report for analysis
        with ``pytest -s``.
        """
        jd = swe.julday(year, month, day, hour)
        flag = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, body_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        tol = _tolerance_for_year(year, body_id)
        dist_from_j2000 = abs(year - 2000)

        # Degradation report visible with pytest -s
        print(
            f"\n  [DEGRADATION] {body_name} @ {date_desc} "
            f"(JD {jd:.1f}, {dist_from_j2000}y from J2000): "
            f"dlon={diff_lon:.6f} deg, dlat={diff_lat:.6f} deg, "
            f"tol={tol:.3f} deg"
        )

        assert diff_lon < tol, (
            f"{body_name} MOSEPH at {date_desc}: "
            f"longitude diff {diff_lon:.6f} deg exceeds "
            f"{tol} deg tolerance\n"
            f"  pyswisseph:   {pos_swe[0]:.6f} deg\n"
            f"  libephemeris: {pos_py[0]:.6f} deg\n"
            f"  Distance from J2000: {dist_from_j2000} years"
        )
        assert diff_lat < tol, (
            f"{body_name} MOSEPH at {date_desc}: "
            f"latitude diff {diff_lat:.6f} deg exceeds "
            f"{tol} deg tolerance\n"
            f"  pyswisseph:   {pos_swe[1]:.6f} deg\n"
            f"  libephemeris: {pos_py[1]:.6f} deg\n"
            f"  Distance from J2000: {dist_from_j2000} years"
        )
