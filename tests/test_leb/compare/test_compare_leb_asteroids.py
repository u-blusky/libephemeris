"""
LEB vs Skyfield Comparison: Main-Belt Asteroids.

Validates the 5 asteroids in LEB (all Pipeline A).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from .conftest import (
    TOLS,
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    lon_error_arcsec,
    _ASTEROID_BODY_IDS,
)


class TestAsteroidPosition:
    """Asteroid position precision."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_position(
        self,
        compare: CompareHelper,
        test_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid longitude and latitude match Skyfield within tolerance."""
        max_lon_err = 0.0
        max_lat_err = 0.0
        worst_jd = 0.0

        dates = filter_asteroid_dates(test_dates_200, body_id)
        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_jd = jd
            if lat_err > max_lat_err:
                max_lat_err = lat_err

        assert max_lon_err < TOLS.POSITION_ARCSEC, (
            f'{body_name}: max lon error = {max_lon_err:.4f}" at JD {worst_jd:.1f}'
        )
        assert max_lat_err < TOLS.POSITION_ARCSEC, (
            f'{body_name}: max lat error = {max_lat_err:.4f}"'
        )


class TestAsteroidSpeed:
    """Asteroid velocity precision."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_speed(
        self,
        compare: CompareHelper,
        test_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        dates = filter_asteroid_dates(test_dates_200, body_id)
        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        tol = (
            TOLS.ASTEROID_SPEED_LON_DEG_DAY
            if body_id in _ASTEROID_BODY_IDS
            else TOLS.SPEED_LON_DEG_DAY
        )
        assert max_err < tol, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestAsteroidDistance:
    """Asteroid distance precision."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_distance(
        self,
        compare: CompareHelper,
        test_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        dates = filter_asteroid_dates(test_dates_200, body_id)
        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
