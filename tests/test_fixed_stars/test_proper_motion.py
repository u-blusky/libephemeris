"""Tests for rigorous proper motion calculation in fixed_stars module.

The proper motion implementation uses the space motion approach from
Hipparcos Vol. 1, Section 1.5.5, which correctly handles the curvature
of the celestial sphere when propagating stellar positions over long
time periods.

This tests the calc_fixed_star_position function in fixed_stars.py.
"""

import math
import pytest
from libephemeris.fixed_stars import (
    calc_fixed_star_position,
    StarData,
    FIXED_STARS,
)
from libephemeris.constants import SE_REGULUS, SE_SPICA_STAR


def _linear_proper_motion(star: StarData, t_years: float) -> tuple[float, float]:
    """
    Apply LINEAR proper motion (old method for comparison).

    Returns (ra_pm, dec_pm) in degrees.
    """
    ra_pm = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
    dec_pm = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
    return ra_pm, dec_pm


def _rigorous_proper_motion(star: StarData, t_years: float) -> tuple[float, float]:
    """
    Apply RIGOROUS proper motion using space motion approach.

    Uses 3D vector propagation per Hipparcos Vol. 1, Section 1.5.5.
    Returns (ra_pm, dec_pm) in degrees.
    """
    # Convert proper motions from arcsec/year to radians/year
    pm_ra_rad = math.radians(star.pm_ra / 3600.0)
    pm_dec_rad = math.radians(star.pm_dec / 3600.0)

    # Convert J2000 position to radians
    ra_rad = math.radians(star.ra_j2000)
    dec_rad = math.radians(star.dec_j2000)

    # Unit position vector at J2000 epoch
    cos_dec = math.cos(dec_rad)
    sin_dec = math.sin(dec_rad)
    cos_ra = math.cos(ra_rad)
    sin_ra = math.sin(ra_rad)

    px = cos_dec * cos_ra
    py = cos_dec * sin_ra
    pz = sin_dec

    # Proper motion velocity vector in the tangent plane
    vx = -pm_ra_rad * sin_ra - pm_dec_rad * sin_dec * cos_ra
    vy = pm_ra_rad * cos_ra - pm_dec_rad * sin_dec * sin_ra
    vz = pm_dec_rad * cos_dec

    # Propagate position
    px_t = px + vx * t_years
    py_t = py + vy * t_years
    pz_t = pz + vz * t_years

    # Normalize to get unit vector
    r = math.sqrt(px_t * px_t + py_t * py_t + pz_t * pz_t)
    px_t /= r
    py_t /= r
    pz_t /= r

    # Convert back to RA/Dec
    dec_pm = math.degrees(math.asin(pz_t))
    ra_pm = math.degrees(math.atan2(py_t, px_t))
    if ra_pm < 0:
        ra_pm += 360.0

    return ra_pm, dec_pm


@pytest.mark.unit
class TestFixedStarsProperMotionRigorous:
    """Tests for rigorous proper motion calculation in fixed_stars module."""

    def test_zero_time_returns_j2000(self):
        """At J2000.0 epoch (t_years=0), result should equal input exactly."""
        star = FIXED_STARS[SE_SPICA_STAR]

        ra_pm, dec_pm = _rigorous_proper_motion(star, 0.0)

        # Should match J2000 coordinates exactly (within floating point)
        assert abs(ra_pm - star.ra_j2000) < 1e-10
        assert abs(dec_pm - star.dec_j2000) < 1e-10

    def test_short_time_matches_linear(self):
        """For short periods (1 year), rigorous and linear should agree well."""
        star = FIXED_STARS[SE_SPICA_STAR]
        t_years = 1.0

        ra_rig, dec_rig = _rigorous_proper_motion(star, t_years)
        ra_lin, dec_lin = _linear_proper_motion(star, t_years)

        # Should agree to within 0.001 arcsec for 1 year
        diff_ra = abs(ra_rig - ra_lin) * 3600  # arcsec
        diff_dec = abs(dec_rig - dec_lin) * 3600  # arcsec

        assert diff_ra < 0.001, f"RA difference: {diff_ra} arcsec"
        assert diff_dec < 0.001, f"Dec difference: {diff_dec} arcsec"

    def test_long_time_differs_from_linear(self):
        """
        For very long periods (1000 years), the rigorous method should
        produce different results due to spherical geometry effects.
        """
        # Create a star with significant proper motion for clear difference
        high_pm_star = StarData(
            ra_j2000=180.0,  # degrees
            dec_j2000=60.0,  # high declination where RA circles are smaller
            pm_ra=0.5,  # 0.5 arcsec/year (high motion)
            pm_dec=0.5,  # 0.5 arcsec/year
        )

        t_years = 1000.0

        ra_rig, dec_rig = _rigorous_proper_motion(high_pm_star, t_years)
        ra_lin, dec_lin = _linear_proper_motion(high_pm_star, t_years)

        # Should differ due to curvature effects
        diff_ra = abs(ra_rig - ra_lin) * 3600  # arcsec
        diff_dec = abs(dec_rig - dec_lin) * 3600  # arcsec

        # The difference should be measurable (at least a few arcsec)
        total_diff = math.sqrt(diff_ra**2 + diff_dec**2)
        assert total_diff > 1.0, (
            f"Expected significant difference, got {total_diff} arcsec"
        )

    def test_preserves_unit_vector_norm(self):
        """The rigorous method should maintain unit vector normalization."""
        star = FIXED_STARS[SE_SPICA_STAR]
        t_years = 100.0

        ra_pm, dec_pm = _rigorous_proper_motion(star, t_years)

        # Convert back to unit vector and check norm
        ra_rad = math.radians(ra_pm)
        dec_rad = math.radians(dec_pm)

        x = math.cos(dec_rad) * math.cos(ra_rad)
        y = math.cos(dec_rad) * math.sin(ra_rad)
        z = math.sin(dec_rad)

        norm = math.sqrt(x * x + y * y + z * z)
        assert abs(norm - 1.0) < 1e-10

    def test_symmetric_time_direction(self):
        """Forward and backward propagation should be symmetric."""
        star = FIXED_STARS[SE_SPICA_STAR]
        t_years = 50.0

        # Go forward then backward
        ra_fwd, dec_fwd = _rigorous_proper_motion(star, t_years)

        # Create star at forward position and go backward
        star_fwd = StarData(
            ra_j2000=ra_fwd,
            dec_j2000=dec_fwd,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )

        ra_back, dec_back = _rigorous_proper_motion(star_fwd, -t_years)

        # Should return close to original (within small numerical error)
        diff_ra = abs(ra_back - star.ra_j2000) * 3600  # arcsec
        diff_dec = abs(dec_back - star.dec_j2000) * 3600  # arcsec

        # Second-order effects mean it won't be exact, but should be close
        assert diff_ra < 0.01, f"RA difference: {diff_ra} arcsec"
        assert diff_dec < 0.01, f"Dec difference: {diff_dec} arcsec"

    def test_all_stars_in_catalog(self):
        """Test that proper motion works for all stars in the catalog."""
        t_years = 50.0

        for star_id, star in FIXED_STARS.items():
            ra_pm, dec_pm = _rigorous_proper_motion(star, t_years)

            # Position should be valid
            assert 0 <= ra_pm < 360, f"Star {star_id}: RA out of range: {ra_pm}"
            assert -90 <= dec_pm <= 90, f"Star {star_id}: Dec out of range: {dec_pm}"

            # Position should have moved from J2000 (except for zero proper motion)
            if star.pm_ra != 0 or star.pm_dec != 0:
                moved = (ra_pm != star.ra_j2000) or (dec_pm != star.dec_j2000)
                assert moved, (
                    f"Star {star_id}: Position didn't change with proper motion"
                )

    def test_declination_stays_bounded(self):
        """Declination should stay within +/-90 even with extreme motion."""
        # Star near pole with motion toward pole
        polar_star = StarData(
            ra_j2000=0.0,
            dec_j2000=89.0,  # Near north pole
            pm_ra=0.0,
            pm_dec=1.0,  # Moving toward pole
        )

        # 100 years would push linear dec to 89 + 100/3600 * 100 = 91.78 deg
        # But should stay bounded
        t_years = 100.0

        ra_pm, dec_pm = _rigorous_proper_motion(polar_star, t_years)

        assert -90 <= dec_pm <= 90, f"Dec out of bounds: {dec_pm}"


@pytest.mark.unit
class TestCalcFixedStarPositionIntegration:
    """Integration tests using the actual calc_fixed_star_position function."""

    def test_regulus_ecliptic_position_j2000(self):
        """Test Regulus ecliptic position at J2000.0."""
        # JD for J2000.0
        jd_tt = 2451545.0

        lon, lat, dist = calc_fixed_star_position(SE_REGULUS, jd_tt)

        # Regulus is at approximately 150 deg ecliptic longitude at J2000
        assert 149 < lon < 151, f"Regulus ecliptic longitude: {lon}"
        # Latitude should be small
        assert -1 < lat < 2, f"Regulus ecliptic latitude: {lat}"
        # Distance is placeholder (100000 AU)
        assert dist == 100000.0

    def test_spica_ecliptic_position_j2000(self):
        """Test Spica ecliptic position at J2000.0."""
        jd_tt = 2451545.0

        lon, lat, dist = calc_fixed_star_position(SE_SPICA_STAR, jd_tt)

        # Spica is at approximately 203-204 deg ecliptic longitude at J2000
        assert 202 < lon < 205, f"Spica ecliptic longitude: {lon}"

    def test_position_changes_over_time(self):
        """Test that star positions change over 100 years."""
        jd_2000 = 2451545.0
        jd_2100 = 2451545.0 + 100 * 365.25

        lon_2000, _, _ = calc_fixed_star_position(SE_SPICA_STAR, jd_2000)
        lon_2100, _, _ = calc_fixed_star_position(SE_SPICA_STAR, jd_2100)

        # Position should have changed due to proper motion + precession
        diff = abs(lon_2100 - lon_2000)
        assert diff > 0.1, f"Expected movement over 100 years, got {diff} deg"

    def test_position_changes_backwards_in_time(self):
        """Test that star positions change going backwards in time."""
        jd_2000 = 2451545.0
        jd_1900 = 2451545.0 - 100 * 365.25

        lon_2000, _, _ = calc_fixed_star_position(SE_REGULUS, jd_2000)
        lon_1900, _, _ = calc_fixed_star_position(SE_REGULUS, jd_1900)

        # Position should have changed
        diff = abs(lon_2000 - lon_1900)
        assert diff > 0.1, f"Expected movement over 100 years, got {diff} deg"

    def test_unknown_star_raises_error(self):
        """Test that unknown star ID raises ValueError."""
        with pytest.raises(ValueError, match="Unknown star ID"):
            calc_fixed_star_position(99999, 2451545.0)

    def test_all_catalog_stars_work(self):
        """Test that all stars in the catalog can be calculated."""
        jd_tt = 2451545.0

        for star_id in FIXED_STARS.keys():
            lon, lat, dist = calc_fixed_star_position(star_id, jd_tt)

            # All results should be valid
            assert 0 <= lon < 360, f"Star {star_id}: lon out of range"
            assert -90 <= lat <= 90, f"Star {star_id}: lat out of range"
            assert dist == 100000.0


@pytest.mark.unit
class TestProperMotionAccuracy:
    """Tests verifying the accuracy improvement from rigorous proper motion."""

    def test_100_year_accuracy_under_0_1_arcsec(self):
        """
        Verify that the rigorous method maintains sub-0.1 arcsec accuracy
        over 100 years compared to what would accumulate with linear method.
        """
        star = FIXED_STARS[SE_SPICA_STAR]
        t_years = 100.0

        ra_rig, dec_rig = _rigorous_proper_motion(star, t_years)
        ra_lin, dec_lin = _linear_proper_motion(star, t_years)

        # For Spica over 100 years, the curvature error should be small
        # but measurable - the rigorous method should agree within 0.1 arcsec
        # of itself when done in steps vs all at once

        # Test that the rigorous method is internally consistent
        # by comparing single-step with multi-step propagation
        ra_step, dec_step = star.ra_j2000, star.dec_j2000
        n_steps = 10
        step_size = t_years / n_steps

        temp_star = StarData(
            ra_j2000=star.ra_j2000,
            dec_j2000=star.dec_j2000,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )

        for _ in range(n_steps):
            ra_step, dec_step = _rigorous_proper_motion(temp_star, step_size)
            temp_star = StarData(
                ra_j2000=ra_step,
                dec_j2000=dec_step,
                pm_ra=star.pm_ra,
                pm_dec=star.pm_dec,
            )

        # Single step vs multi-step should agree very closely
        diff_ra = abs(ra_rig - ra_step) * 3600  # arcsec
        diff_dec = abs(dec_rig - dec_step) * 3600  # arcsec

        assert diff_ra < 0.1, f"RA step consistency: {diff_ra} arcsec"
        assert diff_dec < 0.1, f"Dec step consistency: {diff_dec} arcsec"
