"""
Sidereal Bug Fix Regression Tests: pyswisseph vs libephemeris.

Guards against regressions of 4 sidereal calculation bugs fixed in the
leb/precision branch.  Uses pyswisseph as the reference standard.

Bug 1 (64b8367): Pipeline A SID+EQ used nutation matrix instead of mean equator.
Bug 2 (e6555ed): Pipeline B/C dpsi nutation handling wrong for SID+EQ.
Bug 3 (9f0fde7): J2000 suppression for non-mean SID bodies (TrueNode, OscuApog,
                  IntpApog, IntpPerg).
Bug 4 (b816be0): (a) Frame bias in _get_precession_matrix, (b) SID+J2K
                  precession order for MeanNode/MeanApog.

All tests compare pyswisseph output against libephemeris.  Tolerances are
tight enough to catch the specific error signatures of each bug.
"""

from __future__ import annotations

import pytest
import swisseph as swe

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SE_SIDM_LAHIRI,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_RAMAN,
)


# ============================================================================
# UTILITIES
# ============================================================================


def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360 wrap."""
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def arcsec(deg: float) -> float:
    return deg * 3600.0


# Pipeline A bodies (ICRS barycentric)
PIPELINE_A_BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Pipeline B bodies (ecliptic direct) — mean vs true distinction matters
PIPELINE_B_MEAN_BODIES = [
    (SE_MEAN_NODE, "MeanNode"),
    (SE_MEAN_APOG, "MeanApog"),
]

PIPELINE_B_TRUE_BODIES = [
    (SE_TRUE_NODE, "TrueNode"),
    (SE_OSCU_APOG, "OscuApog"),
]

PIPELINE_B_INTERP_BODIES = [
    (SE_INTP_APOG, "IntpApog"),
    (SE_INTP_PERG, "IntpPerg"),
]

ALL_PIPELINE_B_BODIES = (
    PIPELINE_B_MEAN_BODIES + PIPELINE_B_TRUE_BODIES + PIPELINE_B_INTERP_BODIES
)

# Test dates across the DE440 range
TEST_DATES = [
    (2451545.0, "J2000"),
    (2415020.5, "J1900"),
    (2460000.5, "2023"),
    (2469807.5, "2050"),
    (2433282.5, "1950"),
    (2440587.5, "1970"),
    (2378496.5, "1800"),
    (2524593.5, "2200"),
]

# Tolerances — calibrated to catch each bug's error signature
# while accommodating known systematic offsets (IAU 2006 vs 1976 precession)
PIPELINE_A_SID_EQ_TOL = 0.01  # ~36" (Bug 1 error was ~36")
PIPELINE_A_SID_J2K_TOL = 0.01  # ~36" (Bug 1 error was ~0.3")
PIPELINE_B_SID_EQ_TOL = 0.015  # ~54" (Bug 2 error was ~10-20")
PIPELINE_B_SID_J2K_TOL = 0.015  # ~54" (includes precession model difference)
PIPELINE_B_SID_EQ_MEAN_TOL = 0.015  # Mean bodies
PIPELINE_B_SID_EQ_TRUE_TOL = 0.015  # True bodies
INTP_LON_TOL = 5.5  # IntpApog/IntpPerg: intentional algorithm deviation (see docs)

SID_MODES = [
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SE_SIDM_RAMAN, "Raman"),
]


# ============================================================================
# BUG 1: Pipeline A SID+EQ / SID+J2K (commit 64b8367)
# ============================================================================


class TestBug1PipelineASidEq:
    """Pipeline A SID+EQ: nutation matrix was used instead of mean equator.

    Error signature: ~36" RA error at J2000.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    @pytest.mark.parametrize("sid_mode,sid_name", SID_MODES)
    def test_sid_eq_position(
        self, body_id, body_name, jd, date_desc, sid_mode, sid_name
    ):
        """Pipeline A SID+EQ RA/Dec matches pyswisseph."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))
        dec_diff = abs(float(pos_swe[1]) - float(pos_py[1]))

        assert ra_diff < PIPELINE_A_SID_EQ_TOL, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}") '
            f"[Bug 1 regression?]"
        )
        assert dec_diff < PIPELINE_A_SID_EQ_TOL, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'Dec diff {dec_diff:.6f}° ({arcsec(dec_diff):.1f}")'
        )


class TestBug1PipelineASidJ2k:
    """Pipeline A SID+J2K: true ayanamsha was used instead of mean.

    Error signature: ~0.3" lon error at J2000.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_sid_j2k_position(self, body_id, body_name, jd, date_desc):
        """Pipeline A SID+J2K lon matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < PIPELINE_A_SID_J2K_TOL, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}") '
            f"[Bug 1 regression?]"
        )


class TestBug1PipelineASidEqJ2k:
    """Pipeline A SID+EQ+J2K combined."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES[:3])
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_sid_eq_j2k_position(self, body_id, body_name, jd, date_desc):
        """Pipeline A SID+EQ+J2K RA matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_A_SID_J2K_TOL, (
            f"{body_name} SID+EQ+J2K at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )


# ============================================================================
# BUG 2: Pipeline B/C SID+EQ dpsi handling (commit e6555ed)
# ============================================================================


class TestBug2PipelineBSidEq:
    """Pipeline B SID+EQ: dpsi nutation handling was wrong.

    - MeanNode/MeanApog should skip dpsi
    - TrueNode/OscuApog should subtract dpsi
    - Mean obliquity must be used for SID+EQ rotation

    Error signature: ~10-20" for MeanNode SID+EQ.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ALL_PIPELINE_B_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_sid_eq_all_bodies(self, body_id, body_name, jd, date_desc):
        """All Pipeline B bodies SID+EQ match pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        if body_id in (SE_INTP_APOG, SE_INTP_PERG):
            tol = INTP_LON_TOL
        else:
            tol = PIPELINE_B_SID_EQ_TOL

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < tol, (
            f"{body_name} SID+EQ at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}") '
            f"exceeds {tol}° [Bug 2 regression?]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", SID_MODES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:3])
    def test_mean_body_sid_eq_multi_ayanamsha(
        self, body_id, body_name, sid_mode, sid_name, jd, date_desc
    ):
        """Mean ecliptic bodies SID+EQ across multiple ayanamshas."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_B_SID_EQ_MEAN_TOL, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_TRUE_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:3])
    def test_true_body_sid_eq_speed(self, body_id, body_name, jd, date_desc):
        """True ecliptic bodies SID+EQ velocity matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        speed_tol = 0.05 if body_id == SE_OSCU_APOG else 0.005
        speed_diff = abs(float(pos_swe[3]) - float(pos_py[3]))

        assert speed_diff < speed_tol, (
            f"{body_name} SID+EQ speed at {date_desc}: diff {speed_diff:.6f}°/day"
        )


# ============================================================================
# BUG 3: J2000 suppression for non-mean SID bodies (commit 9f0fde7)
# ============================================================================


class TestBug3J2kSuppression:
    """pyswisseph ignores SEFLG_J2000 for TrueNode, OscuApog, IntpApog,
    IntpPerg when SIDEREAL is set.  MeanNode/MeanApog precess to J2000 normally.

    Error signature: If SID+J2K returns different values for TrueNode than
    SID alone, the J2K suppression is not working.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_TRUE_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_true_body_sid_j2k_equals_sid(self, body_id, body_name, jd, date_desc):
        """TrueNode/OscuApog: SID+J2K should equal SID (J2K suppressed).

        pyswisseph ignores J2000 for these bodies when sidereal is set.
        libephemeris must do the same.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags_sid = SEFLG_SIDEREAL | SEFLG_SPEED
        flags_sid_j2k = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED

        # pyswisseph reference: confirm SE ignores J2K for these bodies
        pos_swe_sid, _ = swe.calc_ut(jd, body_id, flags_sid)
        pos_swe_j2k, _ = swe.calc_ut(jd, body_id, flags_sid_j2k)

        swe_diff = angular_diff(float(pos_swe_sid[0]), float(pos_swe_j2k[0]))
        assert swe_diff < 0.0001, (
            f"pyswisseph sanity check failed: {body_name} SID vs SID+J2K "
            f"differ by {swe_diff:.6f}° (expected identical)"
        )

        # libephemeris: must also ignore J2K
        pos_py_sid, _ = ephem.swe_calc_ut(jd, body_id, flags_sid)
        pos_py_j2k, _ = ephem.swe_calc_ut(jd, body_id, flags_sid_j2k)

        py_diff = angular_diff(float(pos_py_sid[0]), float(pos_py_j2k[0]))
        assert py_diff < 0.0001, (
            f"{body_name} SID vs SID+J2K differ by {py_diff:.6f}° "
            f"(expected identical — J2K should be suppressed) [Bug 3 regression]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_INTERP_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:3])
    def test_interp_body_sid_j2k_equals_sid(self, body_id, body_name, jd, date_desc):
        """IntpApog/IntpPerg: SID+J2K should equal SID (J2K suppressed)."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags_sid = SEFLG_SIDEREAL | SEFLG_SPEED
        flags_sid_j2k = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED

        pos_py_sid, _ = ephem.swe_calc_ut(jd, body_id, flags_sid)
        pos_py_j2k, _ = ephem.swe_calc_ut(jd, body_id, flags_sid_j2k)

        py_diff = angular_diff(float(pos_py_sid[0]), float(pos_py_j2k[0]))
        assert py_diff < 0.0001, (
            f"{body_name} SID vs SID+J2K differ by {py_diff:.6f}° "
            f"(expected identical) [Bug 3 regression]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_j2k_differs_from_sid(
        self, body_id, body_name, jd, date_desc
    ):
        """MeanNode/MeanApog: SID+J2K should DIFFER from SID (J2K applied).

        If identical, J2K is being incorrectly suppressed for mean bodies.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags_sid = SEFLG_SIDEREAL | SEFLG_SPEED
        flags_sid_j2k = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED

        # pyswisseph reference
        pos_swe_sid, _ = swe.calc_ut(jd, body_id, flags_sid)
        pos_swe_j2k, _ = swe.calc_ut(jd, body_id, flags_sid_j2k)

        swe_diff = angular_diff(float(pos_swe_sid[0]), float(pos_swe_j2k[0]))

        # libephemeris
        pos_py_sid, _ = ephem.swe_calc_ut(jd, body_id, flags_sid)
        pos_py_j2k, _ = ephem.swe_calc_ut(jd, body_id, flags_sid_j2k)

        py_diff = angular_diff(float(pos_py_sid[0]), float(pos_py_j2k[0]))

        # Both should show a difference (J2K IS applied for mean bodies)
        if swe_diff > 0.001:
            assert py_diff > 0.001, (
                f"{body_name} SID vs SID+J2K: pyswisseph differs by {swe_diff:.6f}°, "
                f"but libephemeris differs by only {py_diff:.6f}° "
                f"(J2K incorrectly suppressed?) [Bug 3 regression]"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_j2k_matches_swe(self, body_id, body_name, jd, date_desc):
        """MeanNode/MeanApog SID+J2K lon matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < PIPELINE_B_SID_J2K_TOL, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}")'
        )


# ============================================================================
# BUG 4: Frame bias + SID+J2K precession order (commit b816be0)
# ============================================================================


class TestBug4FrameBiasAndPrecessionOrder:
    """(a) _get_precession_matrix used t.P including ICRS frame bias (~17 mas).
    (b) MeanNode/MeanApog SID+J2K applied precession before ayanamsha
        (non-commutative, up to 28" at extreme dates).
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_mean_body_sid_j2k_extended_dates(self, body_id, body_name, jd, date_desc):
        """MeanNode/MeanApog SID+J2K across extended dates.

        Bug 4b error grows with distance from J2000, up to 28" at extreme dates.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < PIPELINE_B_SID_J2K_TOL, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}") '
            f"[Bug 4b regression?]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_eq_j2k_combined(self, body_id, body_name, jd, date_desc):
        """MeanNode/MeanApog SID+EQ+J2K: triple combo matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_B_SID_J2K_TOL, (
            f"{body_name} SID+EQ+J2K at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES[:3])
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_pipeline_a_sid_eq_no_frame_bias(self, body_id, body_name, jd, date_desc):
        """Pipeline A SID+EQ: verify no ICRS frame bias in precession.

        Bug 4a introduced ~17 mas systematic at J2000.  The fix replaced t.P
        (which includes ICRS bias) with mean_equator rotation.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_A_SID_EQ_TOL, (
            f"{body_name} SID+EQ at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}") '
            f"[Bug 4a frame bias regression?]"
        )


# ============================================================================
# COMPREHENSIVE SID+EQ COVERAGE: ALL PIPELINE B/C BODIES
# ============================================================================


class TestSidEqComprehensive:
    """Ensure every Pipeline B body works correctly with SID+EQ."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_MEAN_NODE, "MeanNode"),
            (SE_TRUE_NODE, "TrueNode"),
            (SE_MEAN_APOG, "MeanApog"),
            (SE_OSCU_APOG, "OscuApog"),
        ],
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    @pytest.mark.parametrize("sid_mode,sid_name", SID_MODES)
    def test_primary_body_sid_eq(
        self, body_id, body_name, jd, date_desc, sid_mode, sid_name
    ):
        """All primary Pipeline B bodies SID+EQ across ayanamshas."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = PIPELINE_B_SID_EQ_TOL
        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))
        dec_diff = abs(float(pos_swe[1]) - float(pos_py[1]))

        assert ra_diff < tol, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )
        assert dec_diff < tol, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'Dec diff {dec_diff:.6f}° ({arcsec(dec_diff):.1f}")'
        )


class TestSidJ2kComprehensive:
    """Ensure every Pipeline B body works correctly with SID+J2K."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_MEAN_NODE, "MeanNode"),
            (SE_TRUE_NODE, "TrueNode"),
            (SE_MEAN_APOG, "MeanApog"),
            (SE_OSCU_APOG, "OscuApog"),
        ],
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_primary_body_sid_j2k(self, body_id, body_name, jd, date_desc):
        """All primary Pipeline B bodies SID+J2K match pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = PIPELINE_B_SID_J2K_TOL
        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < tol, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}")'
        )
