"""
LEB vs Skyfield Comparison: Main-Belt Asteroids (Medium Tier).

Validates the 5 asteroids in LEB across the medium tier range (1560-2640).
Asteroid dates are filtered to SPK coverage range (~1920-2080) to avoid
Keplerian-fallback data baked into Chebyshev coefficients.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from tests.test_leb.compare.conftest import (
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    lon_error_arcsec,
)

from .conftest import TOLS_MEDIUM


class TestMediumAsteroidPosition:
    """Asteroid position precision."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_position(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid longitude and latitude match Skyfield within tolerance."""
        dates = filter_asteroid_dates(medium_dates_300, body_id)
        assert len(dates) > 0, f"{body_name}: no dates in SPK coverage range"

        max_lon_err = 0.0
        max_lat_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_jd = jd
            if lat_err > max_lat_err:
                max_lat_err = lat_err

        assert max_lon_err < TOLS_MEDIUM.ASTEROID_ARCSEC, (
            f'{body_name}: max lon error = {max_lon_err:.4f}" at JD {worst_jd:.1f}'
        )
        assert max_lat_err < TOLS_MEDIUM.ASTEROID_ARCSEC, (
            f'{body_name}: max lat error = {max_lat_err:.4f}"'
        )


class TestMediumAsteroidSpeed:
    """Asteroid velocity precision."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_speed(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid speed matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(medium_dates_300, body_id)
        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestMediumAsteroidDistance:
    """Asteroid distance precision."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_distance(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid distance matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(medium_dates_300, body_id)
        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
