"""
LEB vs Skyfield Comparison: Main-Belt Asteroids (Base Tier).

Validates the 5 asteroids in LEB across the base tier range (1860-2140).
Asteroid SPK coverage may limit effective date range.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from tests.test_leb.compare.conftest import (
    ASTEROID_BODIES,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_BASE

# The base tier LEB file was generated without asteroid SPK files.
# All asteroid data uses Keplerian orbital element fallback, producing
# catastrophically wrong positions (7,000-14,000" errors). These tests
# are xfailed until the base tier LEB is regenerated with proper SPK
# data (Phase 4 of the precision improvement plan).
_XFAIL_ASTEROID_REASON = (
    "Base tier LEB generated without asteroid SPK data; "
    "Keplerian fallback produces catastrophic errors. "
    "Regenerate with `poe leb:generate:base:groups` after downloading SPK files."
)


class TestBaseAsteroidPosition:
    """Asteroid position precision."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.xfail(reason=_XFAIL_ASTEROID_REASON, strict=False)
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_position(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid longitude and latitude match Skyfield within tolerance."""
        max_lon_err = 0.0
        max_lat_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                # Asteroid may not have SPK coverage for all dates
                continue

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_jd = jd
            if lat_err > max_lat_err:
                max_lat_err = lat_err

        assert max_lon_err < TOLS_BASE.ASTEROID_ARCSEC, (
            f'{body_name}: max lon error = {max_lon_err:.4f}" at JD {worst_jd:.1f}'
        )
        assert max_lat_err < TOLS_BASE.ASTEROID_ARCSEC, (
            f'{body_name}: max lat error = {max_lat_err:.4f}"'
        )


class TestBaseAsteroidSpeed:
    """Asteroid velocity precision."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.xfail(reason=_XFAIL_ASTEROID_REASON, strict=False)
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_speed(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestBaseAsteroidDistance:
    """Asteroid distance precision."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.xfail(reason=_XFAIL_ASTEROID_REASON, strict=False)
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_distance(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
