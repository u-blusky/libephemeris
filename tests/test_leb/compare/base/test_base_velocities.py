"""
LEB vs Skyfield Comparison: Velocity Precision (Base Tier).

Exhaustive velocity validation for all 30 LEB bodies across base tier (1860-2140).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED, SEFLG_EQUATORIAL

from tests.test_leb.compare.conftest import (
    ALL_LEB_BODIES,
    ASTEROID_BODIES,
    ECLIPTIC_TOLERANCES,
    ICRS_PLANETS,
    CompareHelper,
)

from .conftest import TOLS_BASE


class TestBaseLongitudeVelocity:
    """Longitude velocity (result[3]) precision for all bodies."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        """Longitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_150:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        tol = ecl_tol if ecl_tol is not None else TOLS_BASE.SPEED_LON_DEG_DAY
        assert max_err < tol, (
            f"{body_name}: max lon speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestBaseLatitudeVelocity:
    """Latitude velocity (result[4]) precision for all bodies."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_speed_latitude(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        """Latitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_150:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[4] - leb[4])
            if err > max_err:
                max_err = err
                worst_jd = jd

        # Asteroid latitude velocity is architecturally limited: the
        # ICRS→ecliptic pipeline amplifies errors by 1/geocentric_distance,
        # which is severe for nearby asteroids.  Use a separate tolerance.
        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        if ecl_tol is not None:
            tol = ecl_tol
        else:
            asteroid_ids = {b[0] for b in ASTEROID_BODIES}
            tol = (
                TOLS_BASE.ASTEROID_SPEED_LAT_DEG_DAY
                if body_id in asteroid_ids
                else TOLS_BASE.SPEED_LAT_DEG_DAY
            )
        assert max_err < tol, (
            f"{body_name}: max lat speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestBaseDistanceVelocity:
    """Distance velocity (result[5]) precision for ICRS bodies."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS + ASTEROID_BODIES)
    def test_speed_distance(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        """Distance speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_150:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day at JD {worst_jd:.1f}"
        )


class TestBaseVelocityEquatorial:
    """Velocity in equatorial coordinates for ICRS planets."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_equatorial_velocity(
        self,
        compare: CompareHelper,
        base_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        """Equatorial velocity matches Skyfield within tolerance."""
        flags = SEFLG_SPEED | SEFLG_EQUATORIAL
        max_err = 0.0

        for jd in base_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name}: max equatorial speed error = {max_err:.6f} deg/day"
        )
