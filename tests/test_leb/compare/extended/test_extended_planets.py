"""
LEB vs Skyfield Comparison: ICRS Planet Positions (Extended Tier).

Validates ecliptic longitude, latitude, distance, and speed for all 11 ICRS
pipeline planets across the extended tier range (-5000 to 5000).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_EXT


class TestExtPlanetLongitude:
    """Planet longitude precision across extended tier dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_longitude(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet ecliptic longitude matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.POSITION_ARCSEC, (
            f'{body_name}: max lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestExtPlanetLatitude:
    """Planet latitude precision across extended tier dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_latitude(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet ecliptic latitude matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[1] - leb[1]) * 3600.0
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.POSITION_ARCSEC, (
            f'{body_name}: max lat error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestExtPlanetDistance:
    """Planet distance precision across extended tier dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_distance(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet geocentric distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestExtPlanetSpeed:
    """Planet velocity precision across extended tier dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet longitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )
