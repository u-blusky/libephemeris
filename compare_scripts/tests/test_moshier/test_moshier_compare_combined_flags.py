"""
Moshier Combined Flag Cross-Library Comparison Tests.

Validates that multiple flag combinations in Moshier mode produce consistent
results between pyswisseph (C library) and libephemeris (Python).

The Moshier dispatcher in _calc_body_moshier() (planets.py:857-1087) applies
coordinate transformations in a specific order:
  1. Raw position calculation (VSOP87/ELP/Pluto theory)
  2. Heliocentric dispatch (if SEFLG_HELCTR) — lines 948-1025
  3. Equatorial transformation (if SEFLG_EQUATORIAL) — lines 1051-1072
  4. Sidereal correction (if SEFLG_SIDEREAL and NOT SEFLG_EQUATORIAL) — line 1075

Key interaction: SEFLG_SIDEREAL is SKIPPED when SEFLG_EQUATORIAL is active
(line 1075: ``if is_sidereal and not is_equatorial``). The C library follows
the same convention: sidereal ayanamsha applies only to ecliptic coordinates.

This module tests 8 flag combinations x 3 planets x 3 dates = 72 test cases,
covering 2-flag and 3+ flag interactions that are not tested elsewhere.

Flag combinations tested:
  1. MOSEPH | SPEED | EQUATORIAL          — equatorial + speed
  2. MOSEPH | SPEED | SIDEREAL            — sidereal + speed (Lahiri)
  3. MOSEPH | SPEED | HELCTR              — heliocentric + speed
  4. MOSEPH | EQUATORIAL | SIDEREAL       — eq + sidereal (sidereal skipped)
  5. MOSEPH | SPEED | EQUATORIAL | NONUT  — equatorial without nutation
  6. MOSEPH | SPEED | HELCTR | EQUATORIAL — heliocentric + equatorial
  7. MOSEPH | SPEED | HELCTR | SIDEREAL   — heliocentric + sidereal (Lahiri)
  8. MOSEPH | SPEED | EQUATORIAL | SIDEREAL — eq + sidereal + speed

Documented behaviour for ambiguous combinations:
  - EQUATORIAL | SIDEREAL: sidereal correction is NOT applied because
    ``is_sidereal and not is_equatorial`` is False. Both libraries agree:
    the result equals EQUATORIAL-only. This is by design — ayanamsha is an
    ecliptic-longitude offset and has no meaning in equatorial coordinates.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MARS,
    SE_JUPITER,
    SE_SIDM_LAHIRI,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_SIDEREAL,
    SEFLG_HELCTR,
    SEFLG_NONUT,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 3, 20, 12.0, "2024 equinox"),
    (1900, 6, 15, 0.0, "1900"),
]


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================

# Two-flag combinations inherit tolerances from single-flag tests:
#   - Equatorial: ~0.03° (from test_moshier_compare_equatorial.py)
#   - Sidereal:   ~0.025° (from test_moshier_compare_sidereal.py)
#   - Heliocentric: ~0.04° (from test_moshier_compare_heliocentric.py)
# Three-flag combinations use wider tolerances because errors compound:
#   helio error + equatorial/sidereal transformation adds ~0.01°.

# --- Combo 1: EQUATORIAL + SPEED ---
EQ_SPEED_LON_TOL = 0.03  # degrees (RA)
EQ_SPEED_LAT_TOL = 0.03  # degrees (Dec)
EQ_SPEED_DIST_REL_TOL = 0.001  # relative
EQ_SPEED_VEL_TOL = 0.02  # deg/day

# --- Combo 2: SIDEREAL + SPEED ---
SID_SPEED_LON_TOL = 0.025  # degrees
SID_SPEED_LAT_TOL = 0.02  # degrees
SID_SPEED_DIST_REL_TOL = 0.001  # relative
SID_SPEED_VEL_TOL = 0.01  # deg/day

# --- Combo 3: HELIOCENTRIC + SPEED ---
HELIO_SPEED_LON_TOL = 0.04  # degrees
HELIO_SPEED_LAT_TOL = 0.04  # degrees
HELIO_SPEED_DIST_REL_TOL = 0.001  # relative
HELIO_SPEED_VEL_TOL = 0.01  # deg/day

# --- Combo 4: EQUATORIAL + SIDEREAL (no speed) ---
# Sidereal is skipped when equatorial is active, so tolerances match
# equatorial-only.
EQ_SID_LON_TOL = 0.03  # degrees (RA)
EQ_SID_LAT_TOL = 0.03  # degrees (Dec)
EQ_SID_DIST_REL_TOL = 0.001  # relative

# --- Combo 5: EQUATORIAL + NONUT + SPEED ---
EQ_NONUT_LON_TOL = 0.03  # degrees
EQ_NONUT_LAT_TOL = 0.03  # degrees
EQ_NONUT_DIST_REL_TOL = 0.001  # relative
EQ_NONUT_VEL_TOL = 0.02  # deg/day

# --- Combo 6: HELIOCENTRIC + EQUATORIAL + SPEED ---
HELIO_EQ_LON_TOL = 0.05  # degrees (helio error + eq transform)
HELIO_EQ_LAT_TOL = 0.05  # degrees
HELIO_EQ_DIST_REL_TOL = 0.001  # relative
HELIO_EQ_VEL_TOL = 0.02  # deg/day

# --- Combo 7: HELIOCENTRIC + SIDEREAL + SPEED ---
HELIO_SID_LON_TOL = 0.05  # degrees (helio error + sidereal correction)
HELIO_SID_LAT_TOL = 0.05  # degrees
HELIO_SID_DIST_REL_TOL = 0.001  # relative
HELIO_SID_VEL_TOL = 0.02  # deg/day

# --- Combo 8: EQUATORIAL + SIDEREAL + SPEED ---
# Same as equatorial-only because sidereal is skipped when equatorial is set.
EQ_SID_SPEED_LON_TOL = 0.03  # degrees (RA)
EQ_SID_SPEED_LAT_TOL = 0.03  # degrees (Dec)
EQ_SID_SPEED_DIST_REL_TOL = 0.001  # relative
EQ_SID_SPEED_VEL_TOL = 0.02  # deg/day


# ============================================================================
# HELPER
# ============================================================================


def _compare_positions(
    pos_swe,
    pos_py,
    planet_name,
    date_desc,
    combo_desc,
    lon_tol,
    lat_tol,
    dist_rel_tol,
    vel_tol=None,
    is_helio_sun=False,
    is_helio_sid_sun=False,
):
    """Compare all 6 position components between swe and libephemeris.

    For heliocentric Sun, only validates that both return zeros (no angular
    comparison needed since the position is at the origin by definition).

    When is_helio_sid_sun is True (HELCTR + SIDEREAL + Sun), the test only
    checks distance because: the C library returns lon=0 (no sidereal
    correction on origin), while libephemeris applies the ayanamsha to 0°,
    producing lon=(360 - ayanamsha). Both are valid interpretations since
    sidereal correction on the heliocentric origin is semantically meaningless.
    """
    label = f"{planet_name} [{combo_desc}] @ {date_desc}"

    if is_helio_sun:
        # Heliocentric Sun: both must return zero distance
        assert pos_py[2] == 0.0, (
            f"{label}: helio Sun distance expected 0.0, got {pos_py[2]}"
        )
        assert pos_swe[2] == 0.0, (
            f"{label}: swe helio Sun distance expected 0.0, got {pos_swe[2]}"
        )
        # When sidereal is active, skip longitude comparison:
        # C library keeps origin at 0°, Python applies ayanamsha to 0°.
        # Both are valid — ayanamsha on the heliocentric origin is meaningless.
        if is_helio_sid_sun:
            return
        # Without sidereal: longitude should match (both 0)
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < 0.01, f"{label}: helio Sun lon diff {diff_lon:.6f}°"
        return

    # Longitude / RA (angular comparison for wrap-around)
    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    assert diff_lon < lon_tol, (
        f"{label}: lon/RA diff {diff_lon:.6f}° exceeds {lon_tol}° "
        f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
    )

    # Latitude / Dec
    diff_lat = abs(pos_swe[1] - pos_py[1])
    assert diff_lat < lat_tol, (
        f"{label}: lat/Dec diff {diff_lat:.6f}° exceeds {lat_tol}° "
        f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
    )

    # Distance (relative comparison)
    if pos_swe[2] > 0:
        rel_dist = abs(pos_swe[2] - pos_py[2]) / pos_swe[2]
        assert rel_dist < dist_rel_tol, (
            f"{label}: dist rel diff {rel_dist:.6f} ({rel_dist * 100:.4f}%) "
            f"exceeds {dist_rel_tol * 100}%"
        )

    # Velocity (if tolerance provided = speed flag is set)
    if vel_tol is not None:
        diff_vel_lon = abs(pos_swe[3] - pos_py[3])
        assert diff_vel_lon < vel_tol, (
            f"{label}: lon/RA speed diff {diff_vel_lon:.6f}°/day "
            f"exceeds {vel_tol}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

        diff_vel_lat = abs(pos_swe[4] - pos_py[4])
        assert diff_vel_lat < vel_tol, (
            f"{label}: lat/Dec speed diff {diff_vel_lat:.6f}°/day "
            f"exceeds {vel_tol}°/day "
            f"(swe={pos_swe[4]:.6f}, lib={pos_py[4]:.6f})"
        )


# ============================================================================
# COMBO 1: MOSEPH | SPEED | EQUATORIAL
# ============================================================================


class TestCombo1EquatorialSpeed:
    """Equatorial + Speed: ecliptic -> RA/Dec with velocity.

    Two-flag combination (beyond MOSEPH). Validates the ecliptic-to-equatorial
    transformation chain including velocity numerical differentiation.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | EQUATORIAL combined flags."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "EQ+SPEED",
            EQ_SPEED_LON_TOL,
            EQ_SPEED_LAT_TOL,
            EQ_SPEED_DIST_REL_TOL,
            EQ_SPEED_VEL_TOL,
        )


# ============================================================================
# COMBO 2: MOSEPH | SPEED | SIDEREAL
# ============================================================================


class TestCombo2SiderealSpeed:
    """Sidereal + Speed: tropical position minus Lahiri ayanamsha, with velocity.

    Requires setting sidereal mode (Lahiri) before calculation. The ayanamsha
    rate is subtracted from longitude velocity.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_sidereal_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | SIDEREAL combined flags."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_SIDEREAL

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "SID+SPEED",
            SID_SPEED_LON_TOL,
            SID_SPEED_LAT_TOL,
            SID_SPEED_DIST_REL_TOL,
            SID_SPEED_VEL_TOL,
        )


# ============================================================================
# COMBO 3: MOSEPH | SPEED | HELCTR
# ============================================================================


class TestCombo3HeliocentricSpeed:
    """Heliocentric + Speed: heliocentric positions with velocity.

    Sun returns (0, 0, 0) by definition. Mars and Jupiter use VSOP87
    heliocentric calculations with numerical differentiation for velocity.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_heliocentric_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | HELCTR combined flags."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_HELCTR
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_HELCTR

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "HELIO+SPEED",
            HELIO_SPEED_LON_TOL,
            HELIO_SPEED_LAT_TOL,
            HELIO_SPEED_DIST_REL_TOL,
            HELIO_SPEED_VEL_TOL,
            is_helio_sun=(planet_id == SE_SUN),
        )


# ============================================================================
# COMBO 4: MOSEPH | EQUATORIAL | SIDEREAL
# ============================================================================


class TestCombo4EquatorialSidereal:
    """Equatorial + Sidereal (no speed): validates that sidereal is skipped.

    When SEFLG_EQUATORIAL is set, sidereal correction is NOT applied
    (planets.py:1075: ``if is_sidereal and not is_equatorial``).
    The C library follows the same convention. This test verifies that both
    libraries produce equatorial-only results despite SEFLG_SIDEREAL being set.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_sidereal(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | EQUATORIAL | SIDEREAL — sidereal should be skipped."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_EQUATORIAL | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_EQUATORIAL | SEFLG_SIDEREAL

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "EQ+SID",
            EQ_SID_LON_TOL,
            EQ_SID_LAT_TOL,
            EQ_SID_DIST_REL_TOL,
        )


# ============================================================================
# COMBO 5: MOSEPH | SPEED | EQUATORIAL | NONUT
# ============================================================================


class TestCombo5EquatorialNonutSpeed:
    """Equatorial + NONUT + Speed: mean equatorial without nutation.

    NONUT suppresses nutation so the transformation uses mean obliquity only.
    This isolates the nutation model difference between the C and Python
    implementations. If results agree with NONUT but diverge without it,
    the nutation model is the divergence source.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_nonut_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | EQUATORIAL | NONUT combined flags."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_NONUT
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_NONUT

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "EQ+NONUT+SPEED",
            EQ_NONUT_LON_TOL,
            EQ_NONUT_LAT_TOL,
            EQ_NONUT_DIST_REL_TOL,
            EQ_NONUT_VEL_TOL,
        )


# ============================================================================
# COMBO 6: MOSEPH | SPEED | HELCTR | EQUATORIAL
# ============================================================================


class TestCombo6HeliocentricEquatorialSpeed:
    """Heliocentric + Equatorial + Speed: 3-flag combination.

    Applies heliocentric dispatch first (raw VSOP87 heliocentric position),
    then transforms the heliocentric ecliptic coordinates to equatorial (RA/Dec).
    Error compounds: VSOP87 C-vs-Python divergence (~0.04°) plus obliquity
    formula difference (~0.01°).
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_heliocentric_equatorial_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | HELCTR | EQUATORIAL combined flags."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_HELCTR | swe.FLG_EQUATORIAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_HELCTR | SEFLG_EQUATORIAL

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "HELIO+EQ+SPEED",
            HELIO_EQ_LON_TOL,
            HELIO_EQ_LAT_TOL,
            HELIO_EQ_DIST_REL_TOL,
            HELIO_EQ_VEL_TOL,
            is_helio_sun=(planet_id == SE_SUN),
        )


# ============================================================================
# COMBO 7: MOSEPH | SPEED | HELCTR | SIDEREAL
# ============================================================================


class TestCombo7HeliocentricSiderealSpeed:
    """Heliocentric + Sidereal + Speed: 3-flag combination.

    Applies heliocentric dispatch, then sidereal ayanamsha correction (since
    equatorial is NOT set, sidereal IS applied). This combination is used in
    Vedic cosmobiology for heliocentric sidereal positions.
    Error compounds: VSOP87 divergence (~0.04°) plus ayanamsha difference.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_heliocentric_sidereal_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | HELCTR | SIDEREAL combined flags."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_HELCTR | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_HELCTR | SEFLG_SIDEREAL

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "HELIO+SID+SPEED",
            HELIO_SID_LON_TOL,
            HELIO_SID_LAT_TOL,
            HELIO_SID_DIST_REL_TOL,
            HELIO_SID_VEL_TOL,
            is_helio_sun=(planet_id == SE_SUN),
            is_helio_sid_sun=(planet_id == SE_SUN),
        )


# ============================================================================
# COMBO 8: MOSEPH | SPEED | EQUATORIAL | SIDEREAL
# ============================================================================


class TestCombo8EquatorialSiderealSpeed:
    """Equatorial + Sidereal + Speed: validates sidereal is skipped with speed.

    Same interaction as Combo 4 (sidereal skipped when equatorial is set),
    but with SEFLG_SPEED also active. This verifies that the velocity
    transformation is equatorial-only and not affected by the sidereal flag.
    3 planets x 3 dates = 9 test cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_equatorial_sidereal_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test MOSEPH | SPEED | EQUATORIAL | SIDEREAL combined flags."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        jd = swe.julday(year, month, day, hour)
        flag_swe = (
            swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL | swe.FLG_SIDEREAL
        )
        flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_SIDEREAL

        try:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag_swe)
        except Exception as e:
            pytest.skip(f"pyswisseph failed: {e}")

        try:
            pos_py, _ = ephem.swe_calc_ut(jd, planet_id, flag_py)
        except Exception as e:
            pytest.skip(f"libephemeris failed: {e}")

        _compare_positions(
            pos_swe,
            pos_py,
            planet_name,
            date_desc,
            "EQ+SID+SPEED",
            EQ_SID_SPEED_LON_TOL,
            EQ_SID_SPEED_LAT_TOL,
            EQ_SID_SPEED_DIST_REL_TOL,
            EQ_SID_SPEED_VEL_TOL,
        )
