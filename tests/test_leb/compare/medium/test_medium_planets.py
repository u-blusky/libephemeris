"""
LEB vs Skyfield Comparison: ICRS Planet Positions (Medium Tier).

Validates ecliptic longitude, latitude, distance, and speed for all 11 ICRS
pipeline planets across the medium tier range (1560-2640).
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

from .conftest import TOLS_MEDIUM


class TestMediumPlanetLongitude:
    """Planet longitude precision across medium tier dates."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_longitude(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet ecliptic longitude matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.POSITION_ARCSEC, (
            f'{body_name}: max lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestMediumPlanetLatitude:
    """Planet latitude precision across medium tier dates."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_latitude(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet ecliptic latitude matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[1] - leb[1]) * 3600.0
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.POSITION_ARCSEC, (
            f'{body_name}: max lat error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestMediumPlanetDistance:
    """Planet distance precision across medium tier dates."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_distance(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet geocentric distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestMediumPlanetSpeed:
    """Planet velocity precision across medium tier dates."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet longitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_dates_300:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestMediumPlanetAllComponents:
    """Combined check of all 6 output components."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_all_components(
        self,
        compare: CompareHelper,
        medium_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        """All 6 output components match within tolerance."""
        max_lon_err = 0.0
        max_lat_err = 0.0
        max_dist_err = 0.0
        max_speed_lon_err = 0.0

        for jd in medium_dates_150:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            dist_err = abs(ref[2] - leb[2])
            speed_lon_err = abs(ref[3] - leb[3])

            max_lon_err = max(max_lon_err, lon_err)
            max_lat_err = max(max_lat_err, lat_err)
            max_dist_err = max(max_dist_err, dist_err)
            max_speed_lon_err = max(max_speed_lon_err, speed_lon_err)

        assert max_lon_err < TOLS_MEDIUM.POSITION_ARCSEC, (
            f'{body_name}: lon error {max_lon_err:.4f}"'
        )
        assert max_lat_err < TOLS_MEDIUM.POSITION_ARCSEC, (
            f'{body_name}: lat error {max_lat_err:.4f}"'
        )
        assert max_dist_err < TOLS_MEDIUM.DISTANCE_AU, (
            f"{body_name}: dist error {max_dist_err:.2e} AU"
        )
        assert max_speed_lon_err < TOLS_MEDIUM.SPEED_LON_DEG_DAY, (
            f"{body_name}: speed error {max_speed_lon_err:.6f} deg/day"
        )
