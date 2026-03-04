"""
LEB vs Skyfield Comparison: Velocity Precision.

Dedicated, exhaustive velocity validation for all 30 LEB bodies.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from .conftest import (
    TOLS,
    ICRS_PLANETS,
    ECLIPTIC_BODIES,
    ASTEROID_BODIES,
    HYPOTHETICAL_BODIES,
    CompareHelper,
)

ALL_LEB_BODIES = ICRS_PLANETS + ECLIPTIC_BODIES + ASTEROID_BODIES + HYPOTHETICAL_BODIES


class TestLongitudeVelocity:
    """Longitude velocity (result[3]) precision for all bodies."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        test_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        """Longitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in test_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.SPEED_LON_DEG_DAY, (
            f"{body_name}: max lon speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestLatitudeVelocity:
    """Latitude velocity (result[4]) precision for all bodies."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_speed_latitude(
        self,
        compare: CompareHelper,
        test_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        """Latitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in test_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[4] - leb[4])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.SPEED_LAT_DEG_DAY, (
            f"{body_name}: max lat speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestDistanceVelocity:
    """Distance velocity (result[5]) precision for ICRS bodies."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS + ASTEROID_BODIES)
    def test_speed_distance(
        self,
        compare: CompareHelper,
        test_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        """Distance speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in test_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day at JD {worst_jd:.1f}"
        )


class TestVelocityEquatorial:
    """Velocity in equatorial coordinates for ICRS planets."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_equatorial_velocity(
        self,
        compare: CompareHelper,
        test_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Equatorial velocity matches Skyfield within tolerance."""
        from libephemeris.constants import SEFLG_EQUATORIAL

        flags = SEFLG_SPEED | SEFLG_EQUATORIAL
        max_err = 0.0

        for jd in test_dates_50:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS.SPEED_LON_DEG_DAY, (
            f"{body_name}: max equatorial speed error = {max_err:.6f} deg/day"
        )
