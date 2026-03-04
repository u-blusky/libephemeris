"""
LEB vs Skyfield Comparison: Lunar/Ecliptic Bodies (Base Tier).

Validates all 6 Pipeline B (ecliptic-direct) bodies with per-body tolerances
across the base tier range (1860-2140).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from tests.test_leb.compare.conftest import (
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_BASE


class TestBaseLunarLongitude:
    """Lunar body longitude precision with per-body tolerance."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_longitude(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body longitude matches Skyfield within per-body tolerance."""
        tol_arcsec = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "lon", TOLS_BASE.ECLIPTIC_ARCSEC
        )
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < tol_arcsec, (
            f'{body_name}: max lon error = {max_err:.4f}" at JD {worst_jd:.1f} '
            f'(tol={tol_arcsec}")'
        )


class TestBaseLunarLatitude:
    """Lunar body latitude precision."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_latitude(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body latitude matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[1] - leb[1]) * 3600.0
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.ECLIPTIC_ARCSEC, (
            f'{body_name}: max lat error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestBaseLunarSpeed:
    """Lunar body velocity precision with per-body tolerance."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body longitude speed matches Skyfield within per-body tolerance."""
        tol_speed = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "speed", TOLS_BASE.SPEED_LON_DEG_DAY
        )
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < tol_speed, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f} "
            f"(tol={tol_speed})"
        )


class TestBaseLunarDistance:
    """Lunar body distance precision."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_distance(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
