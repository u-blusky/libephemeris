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

Total: 45 test cases (3 classes x 5 planets x 3 dates).
"""

from __future__ import annotations

import pytest

swe = pytest.importorskip("swisseph", reason="pyswisseph required for comparison tests")

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
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
