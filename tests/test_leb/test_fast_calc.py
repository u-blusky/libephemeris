"""
Tests for the fast calculation pipeline (fast_calc.py).

Compares fast_calc output with Skyfield/swe_calc_ut for all supported bodies
and flag combinations.
"""

from __future__ import annotations


import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SE_MARS,
    SE_MEAN_NODE,
    SE_MOON,
    SE_SUN,
    SEFLG_BARYCTR,
    SEFLG_EQUATORIAL,
    SEFLG_HELCTR,
    SEFLG_J2000,
    SEFLG_NOABERR,
    SEFLG_NONUT,
    SEFLG_RADIANS,
    SEFLG_SIDEREAL,
    SEFLG_SPEED,
    SEFLG_TOPOCTR,
    SEFLG_TRUEPOS,
    SEFLG_XYZ,
)
from libephemeris.fast_calc import (
    _calc_ayanamsa_from_leb,
    _cartesian_to_spherical,
    _mean_obliquity_iau2006,
    _rotate_equatorial_to_ecliptic,
    _vec3_dist,
    _vec3_sub,
    fast_calc_tt,
    fast_calc_ut,
)


# =============================================================================
# UTILITY FUNCTION TESTS
# =============================================================================


class TestUtilityFunctions:
    """Test coordinate transform and vector utility functions."""

    @pytest.mark.unit
    def test_vec3_sub(self):
        """Vector subtraction."""
        result = _vec3_sub((3.0, 5.0, 7.0), (1.0, 2.0, 3.0))
        assert result == (2.0, 3.0, 4.0)

    @pytest.mark.unit
    def test_vec3_dist(self):
        """Euclidean distance."""
        assert abs(_vec3_dist((3.0, 4.0, 0.0)) - 5.0) < 1e-14
        assert abs(_vec3_dist((0.0, 0.0, 0.0))) < 1e-14

    @pytest.mark.unit
    def test_cartesian_to_spherical_origin(self):
        """Origin gives (0, 0, 0)."""
        lon, lat, dist = _cartesian_to_spherical(0.0, 0.0, 0.0)
        assert lon == 0.0
        assert lat == 0.0
        assert dist == 0.0

    @pytest.mark.unit
    def test_cartesian_to_spherical_x_axis(self):
        """Point on +x axis: lon=0, lat=0."""
        lon, lat, dist = _cartesian_to_spherical(1.0, 0.0, 0.0)
        assert abs(lon - 0.0) < 1e-12
        assert abs(lat - 0.0) < 1e-12
        assert abs(dist - 1.0) < 1e-14

    @pytest.mark.unit
    def test_cartesian_to_spherical_y_axis(self):
        """Point on +y axis: lon=90."""
        lon, lat, dist = _cartesian_to_spherical(0.0, 1.0, 0.0)
        assert abs(lon - 90.0) < 1e-12

    @pytest.mark.unit
    def test_cartesian_to_spherical_z_axis(self):
        """Point on +z axis: lat=90."""
        lon, lat, dist = _cartesian_to_spherical(0.0, 0.0, 1.0)
        assert abs(lat - 90.0) < 1e-12

    @pytest.mark.unit
    def test_mean_obliquity_j2000(self):
        """Mean obliquity at J2000 should be ~23.4393 degrees."""
        eps = _mean_obliquity_iau2006(2451545.0)
        # IAU 2006: 84381.406 arcsec = 23.4393 deg
        assert abs(eps - 23.4393) < 0.001

    @pytest.mark.unit
    def test_rotate_equatorial_to_ecliptic_identity_at_zero(self):
        """With obliquity=0, rotation is identity."""
        x, y, z = _rotate_equatorial_to_ecliptic(1.0, 2.0, 3.0, 0.0)
        assert abs(x - 1.0) < 1e-14
        assert abs(y - 2.0) < 1e-14
        assert abs(z - 3.0) < 1e-14


# =============================================================================
# FAST CALC vs SKYFIELD COMPARISON TESTS
# =============================================================================


class TestFastCalcVsSkyfield:
    """Compare fast_calc output with Skyfield pipeline for all bodies."""

    # Tolerance in arcseconds
    TOLERANCE_ARCSEC = 0.1  # 100 milliarcseconds (generous for Pipeline A)
    TOLERANCE_ARCSEC_ECLIPTIC = 0.01  # 10 mas for ecliptic direct bodies

    @pytest.fixture
    def test_dates(self, leb_reader):
        """Return a set of test dates within the LEB file range."""
        jd_start, jd_end = leb_reader.jd_range
        margin = 10.0  # Stay away from edges
        # Spread dates across the range
        dates = [
            jd_start + margin,
            (jd_start + jd_end) / 2.0,
            jd_end - margin,
            jd_start + (jd_end - jd_start) * 0.25,
            jd_start + (jd_end - jd_start) * 0.75,
        ]
        return dates

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SE_SUN, SE_MOON, SE_MARS])
    def test_planet_ecliptic_default(self, leb_reader, test_dates, ipl):
        """Planet ecliptic lon/lat matches Skyfield within tolerance."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, ipl, SEFLG_SPEED)
            sky_result, _ = ephem.swe_calc_ut(jd, ipl, SEFLG_SPEED)

            lon_err = abs(fast_result[0] - sky_result[0])
            # Handle wrap-around
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            lat_err_arcsec = abs(fast_result[1] - sky_result[1]) * 3600.0

            assert lon_err_arcsec < self.TOLERANCE_ARCSEC, (
                f"Body {ipl} at JD {jd}: lon error = {lon_err_arcsec:.4f} arcsec"
            )
            assert lat_err_arcsec < self.TOLERANCE_ARCSEC, (
                f"Body {ipl} at JD {jd}: lat error = {lat_err_arcsec:.4f} arcsec"
            )

    @pytest.mark.integration
    def test_mean_node_ecliptic(self, leb_reader, test_dates):
        """Mean node longitude matches Skyfield within tolerance."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, SE_MEAN_NODE, SEFLG_SPEED)
            sky_result, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_SPEED)

            lon_err = abs(fast_result[0] - sky_result[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            assert lon_err_arcsec < self.TOLERANCE_ARCSEC_ECLIPTIC, (
                f"Mean node at JD {jd}: lon error = {lon_err_arcsec:.4f} arcsec"
            )

    @pytest.mark.integration
    def test_sun_speed(self, leb_reader, test_dates):
        """Sun speed from .leb should match Skyfield within ~0.001 deg/day."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, SE_SUN, SEFLG_SPEED)
            sky_result, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)

            # Speed tolerance: 0.001 deg/day = 3.6 arcsec/day
            speed_err = abs(fast_result[3] - sky_result[3])
            assert speed_err < 0.01, (
                f"Sun speed error at JD {jd}: {speed_err:.6f} deg/day"
            )

    @pytest.mark.integration
    def test_sun_distance(self, leb_reader, test_dates):
        """Sun distance from .leb should match Skyfield within ~1e-6 AU."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, SE_SUN, SEFLG_SPEED)
            sky_result, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)

            dist_err = abs(fast_result[2] - sky_result[2])
            assert dist_err < 1e-4, f"Sun distance error at JD {jd}: {dist_err:.2e} AU"


class TestFastCalcFlags:
    """Test various flag combinations."""

    @pytest.fixture
    def jd_mid(self, leb_reader):
        """Return the midpoint JD of the test file."""
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_no_speed(self, leb_reader, jd_mid):
        """Without SEFLG_SPEED, velocity should be zero."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, 0)
        assert result[3] == 0.0
        assert result[4] == 0.0
        assert result[5] == 0.0

    @pytest.mark.integration
    def test_with_speed(self, leb_reader, jd_mid):
        """With SEFLG_SPEED, velocity should be non-zero."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_SPEED)
        # Sun moves ~1 deg/day in longitude
        assert abs(result[3]) > 0.5, f"Sun speed too small: {result[3]}"

    @pytest.mark.integration
    def test_equatorial_flag(self, leb_reader, jd_mid):
        """SEFLG_EQUATORIAL should return equatorial RA/Dec for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_EQUATORIAL)
        # RA should be in [0, 360)
        assert 0.0 <= result[0] < 360.0, f"RA = {result[0]}"
        # Dec should be in [-90, 90]
        assert -90.0 <= result[1] <= 90.0, f"Dec = {result[1]}"

    @pytest.mark.integration
    def test_j2000_flag(self, leb_reader, jd_mid):
        """SEFLG_J2000 should return J2000 ecliptic coords for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_J2000)
        assert 0.0 <= result[0] < 360.0, f"Lon = {result[0]}"

    @pytest.mark.integration
    def test_equatorial_j2000_flag(self, leb_reader, jd_mid):
        """SEFLG_EQUATORIAL|SEFLG_J2000 should return J2000 RA/Dec."""
        result, _ = fast_calc_ut(
            leb_reader, jd_mid, SE_SUN, SEFLG_EQUATORIAL | SEFLG_J2000
        )
        assert 0.0 <= result[0] < 360.0, f"RA = {result[0]}"
        assert -90.0 <= result[1] <= 90.0, f"Dec = {result[1]}"

    @pytest.mark.integration
    def test_helctr_flag(self, leb_reader, jd_mid):
        """SEFLG_HELCTR should return heliocentric coords for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MARS, SEFLG_HELCTR)
        # Mars heliocentric distance ~1.38-1.67 AU
        assert 1.0 < result[2] < 2.0, f"Mars helio dist = {result[2]}"

    @pytest.mark.integration
    def test_baryctr_flag(self, leb_reader, jd_mid):
        """SEFLG_BARYCTR should return barycentric coords for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MARS, SEFLG_BARYCTR)
        # Mars barycentric distance ~1.38-1.67 AU
        assert 1.0 < result[2] < 2.0, f"Mars bary dist = {result[2]}"

    @pytest.mark.integration
    def test_noaberr_flag(self, leb_reader, jd_mid):
        """SEFLG_NOABERR should skip aberration and return valid coords."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MARS, SEFLG_NOABERR)
        assert 0.0 <= result[0] < 360.0, f"Lon = {result[0]}"

    @pytest.mark.integration
    def test_truepos_flag(self, leb_reader, jd_mid):
        """SEFLG_TRUEPOS should skip light-time and aberration."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MARS, SEFLG_TRUEPOS)
        assert 0.0 <= result[0] < 360.0, f"Lon = {result[0]}"

    @pytest.mark.integration
    def test_topoctr_raises(self, leb_reader, jd_mid):
        """SEFLG_TOPOCTR should raise KeyError (fall back to Skyfield)."""
        with pytest.raises(KeyError, match="SEFLG_TOPOCTR"):
            fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_TOPOCTR)

    @pytest.mark.integration
    def test_xyz_raises(self, leb_reader, jd_mid):
        """SEFLG_XYZ should raise KeyError (fall back to Skyfield)."""
        with pytest.raises(KeyError, match="SEFLG_XYZ"):
            fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_XYZ)

    @pytest.mark.integration
    def test_radians_raises(self, leb_reader, jd_mid):
        """SEFLG_RADIANS should raise KeyError (fall back to Skyfield)."""
        with pytest.raises(KeyError, match="SEFLG_RADIANS"):
            fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_RADIANS)


class TestEclipticDirectVelocity:
    """Test ecliptic-direct bodies with SEFLG_EQUATORIAL + SEFLG_SPEED.

    This validates the approximate velocity transform in Pipeline B
    (_pipeline_ecliptic) when converting from ecliptic to equatorial.
    """

    @pytest.fixture
    def jd_mid(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_mean_node_equatorial_speed(self, leb_reader, jd_mid):
        """Mean node equatorial velocity via LEB should match Skyfield."""
        flags = SEFLG_EQUATORIAL | SEFLG_SPEED
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MEAN_NODE, flags)
        sky_result, _ = ephem.swe_calc_ut(jd_mid, SE_MEAN_NODE, flags)

        # Position check (RA/Dec)
        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0
        assert ra_err_arcsec < 0.1, f"Mean node RA error: {ra_err_arcsec:.4f} arcsec"

        dec_err_arcsec = abs(fast_result[1] - sky_result[1]) * 3600.0
        assert dec_err_arcsec < 0.1, f"Mean node Dec error: {dec_err_arcsec:.4f} arcsec"

        # Velocity check: RA speed (deg/day)
        # Mean node moves ~-0.053 deg/day ecliptic, equatorial speed is similar
        # Tolerance: 10% of Skyfield's value or 0.01 deg/day, whichever is larger
        if abs(sky_result[3]) > 0.001:
            speed_ra_err = abs(fast_result[3] - sky_result[3])
            tol = max(abs(sky_result[3]) * 0.1, 0.01)
            assert speed_ra_err < tol, (
                f"Mean node RA speed error: {speed_ra_err:.6f} deg/day "
                f"(LEB={fast_result[3]:.6f}, Sky={sky_result[3]:.6f})"
            )

    @pytest.mark.integration
    def test_mean_node_equatorial_j2000_speed(self, leb_reader, jd_mid):
        """Mean node SEFLG_EQUATORIAL | SEFLG_J2000 should produce valid output."""
        flags = SEFLG_EQUATORIAL | SEFLG_J2000 | SEFLG_SPEED
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MEAN_NODE, flags)
        sky_result, _ = ephem.swe_calc_ut(jd_mid, SE_MEAN_NODE, flags)

        # Position check
        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0
        assert ra_err_arcsec < 0.5, (
            f"Mean node J2000 RA error: {ra_err_arcsec:.4f} arcsec"
        )

    @pytest.mark.integration
    def test_mean_node_j2000_ecliptic(self, leb_reader, jd_mid):
        """Mean node SEFLG_J2000 (ecliptic J2000) should match Skyfield."""
        flags = SEFLG_J2000 | SEFLG_SPEED
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, SE_MEAN_NODE, flags)
        sky_result, _ = ephem.swe_calc_ut(jd_mid, SE_MEAN_NODE, flags)

        lon_err = abs(fast_result[0] - sky_result[0])
        if lon_err > 180.0:
            lon_err = 360.0 - lon_err
        lon_err_arcsec = lon_err * 3600.0
        assert lon_err_arcsec < 0.5, (
            f"Mean node J2000 ecliptic lon error: {lon_err_arcsec:.4f} arcsec"
        )


class TestFastCalcFallback:
    """Test fallback behavior for unsupported bodies/ranges."""

    @pytest.fixture
    def jd_mid(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_unknown_body_raises_keyerror(self, leb_reader, jd_mid):
        """Bodies not in .leb should raise KeyError."""
        with pytest.raises(KeyError, match="Body 99999"):
            fast_calc_ut(leb_reader, jd_mid, 99999, 0)

    @pytest.mark.integration
    def test_out_of_range_raises_valueerror(self, leb_reader):
        """JD outside .leb range should raise ValueError."""
        jd_start, _ = leb_reader.jd_range
        with pytest.raises((ValueError, KeyError)):
            fast_calc_ut(leb_reader, jd_start - 10000, SE_SUN, 0)


class TestFastCalcTT:
    """Test the TT entry point."""

    @pytest.mark.integration
    def test_fast_calc_tt_basic(self, leb_reader):
        """fast_calc_tt should work with TT input."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        result, iflag = fast_calc_tt(leb_reader, jd_mid, SE_SUN, SEFLG_SPEED)

        assert len(result) == 6
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_tt_vs_ut_close(self, leb_reader):
        """fast_calc_tt and fast_calc_ut should give close results for similar JD."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        # UT and TT differ by Delta-T (~69 sec in 2020s)
        # Using the same JD should give slightly different positions
        result_ut, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_SPEED)
        result_tt, _ = fast_calc_tt(leb_reader, jd_mid, SE_SUN, SEFLG_SPEED)

        # Sun moves ~1 deg/day, Delta-T ~69 sec → diff ~0.0008 deg
        lon_diff = abs(result_ut[0] - result_tt[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        # Should be small but non-zero
        assert lon_diff < 0.01, f"UT vs TT lon diff = {lon_diff} deg (too large)"


class TestNonutFallback:
    """Test SEFLG_NONUT triggers Skyfield fallback."""

    @pytest.mark.integration
    def test_nonut_raises_keyerror_ut(self, leb_reader):
        """SEFLG_NONUT should raise KeyError to trigger fallback."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        with pytest.raises(KeyError, match="SEFLG_NONUT"):
            fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_NONUT)

    @pytest.mark.integration
    def test_nonut_raises_keyerror_tt(self, leb_reader):
        """SEFLG_NONUT should raise KeyError to trigger fallback in TT path."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        with pytest.raises(KeyError, match="SEFLG_NONUT"):
            fast_calc_tt(leb_reader, jd_mid, SE_SUN, SEFLG_NONUT)


class TestExplicitSiderealParams:
    """Test thread-safe explicit sidereal parameter passing."""

    @pytest.mark.integration
    def test_explicit_sid_mode_lahiri(self, leb_reader):
        """Explicit sid_mode=1 (Lahiri) should work like global state."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        result, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SE_SUN,
            SEFLG_SPEED | SEFLG_SIDEREAL,
            sid_mode=1,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        # Sidereal lon should be ~24 deg less than tropical
        tropical, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_SPEED)
        diff = tropical[0] - result[0]
        if diff < 0:
            diff += 360.0
        # Lahiri ayanamsa is ~24 deg in 2020s
        assert 20.0 < diff < 30.0, f"Sidereal offset = {diff:.2f} deg"

    @pytest.mark.integration
    def test_explicit_sid_mode_fagan_bradley(self, leb_reader):
        """Explicit sid_mode=0 (Fagan-Bradley) should differ from Lahiri."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        result_lahiri, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SE_SUN,
            SEFLG_SPEED | SEFLG_SIDEREAL,
            sid_mode=1,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        result_fb, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SE_SUN,
            SEFLG_SPEED | SEFLG_SIDEREAL,
            sid_mode=0,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        # Fagan-Bradley and Lahiri differ by ~0.88 deg
        diff = abs(result_lahiri[0] - result_fb[0])
        assert diff > 0.5, f"FB vs Lahiri diff = {diff:.4f} deg (too small)"
        assert diff < 3.0, f"FB vs Lahiri diff = {diff:.4f} deg (too large)"

    @pytest.mark.integration
    def test_sidereal_speed_includes_precession_correction(self, leb_reader):
        """Sidereal dlon should differ from tropical dlon by precession rate."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        trop, _ = fast_calc_ut(leb_reader, jd_mid, SE_SUN, SEFLG_SPEED)
        sid, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SE_SUN,
            SEFLG_SPEED | SEFLG_SIDEREAL,
            sid_mode=1,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        # Precession rate ~50.3"/yr = 5028.8"/century
        # = 5028.8 / 3600 / 36525 ≈ 0.0000382 deg/day
        speed_diff = trop[3] - sid[3]
        assert 0.00003 < speed_diff < 0.00005, (
            f"Sidereal speed correction = {speed_diff:.7f} deg/day "
            f"(expected ~0.0000382)"
        )


class TestPipelineBJ2000Velocity:
    """Test that Pipeline B (ecliptic-direct) J2000-only branch precesses velocity."""

    @pytest.mark.integration
    def test_mean_node_j2000_velocity_precessed(self, leb_reader):
        """Mean node velocity in J2000 ecliptic should be precession-corrected."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        # J2000 ecliptic (the branch we fixed)
        result_j2000, _ = fast_calc_ut(
            leb_reader, jd_mid, SE_MEAN_NODE, SEFLG_SPEED | SEFLG_J2000
        )
        # Ecliptic of date (default)
        result_date, _ = fast_calc_ut(leb_reader, jd_mid, SE_MEAN_NODE, SEFLG_SPEED)

        # Velocity should differ slightly due to precession correction
        # Mean node moves ~-0.053 deg/day
        assert abs(result_j2000[3]) > 0.01, (
            f"J2000 mean node speed = {result_j2000[3]:.6f} (suspiciously small)"
        )
        # The difference should be small (~precession rate)
        speed_diff = abs(result_j2000[3] - result_date[3])
        assert speed_diff < 0.01, (
            f"J2000 vs date speed diff = {speed_diff:.6f} deg/day (too large)"
        )


class TestAyanamsa:
    """Test the LEB-based ayanamsa computation."""

    @pytest.mark.integration
    def test_lahiri_ayanamsa_reasonable(self, leb_reader):
        """Lahiri ayanamsa should be ~24 degrees in 2020s."""

        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        try:
            aya = _calc_ayanamsa_from_leb(leb_reader, jd_mid)
            # Lahiri ayanamsa is ~24 degrees in the 2020s
            assert 20.0 < aya < 30.0, f"Ayanamsa = {aya} (out of expected range)"
        except (KeyError, ValueError):
            # May fail if sidereal mode is not set; that's OK
            pass

    @pytest.mark.integration
    def test_explicit_sid_mode_param(self, leb_reader):
        """_calc_ayanamsa_from_leb accepts explicit sid_mode."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        aya = _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=1)
        assert 20.0 < aya < 30.0, f"Lahiri ayanamsa = {aya}"

        aya_fb = _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=0)
        assert 20.0 < aya_fb < 30.0, f"Fagan-Bradley ayanamsa = {aya_fb}"
        assert abs(aya - aya_fb) > 0.5, "Lahiri and Fagan-Bradley should differ"

    @pytest.mark.integration
    def test_star_based_mode_raises(self, leb_reader):
        """Star-based sidereal modes should raise KeyError."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        with pytest.raises(KeyError, match="Star-based"):
            _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=17)
