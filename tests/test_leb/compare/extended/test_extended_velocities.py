"""
LEB vs Skyfield Comparison: Velocity Precision (Extended Tier).

Exhaustive velocity validation for all LEB bodies (excluding asteroids
whose SPK coverage does not span the extended range).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED, SEFLG_EQUATORIAL

from tests.test_leb.compare.conftest import (
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    HYPOTHETICAL_BODIES,
    ICRS_PLANETS,
    CompareHelper,
)

from .conftest import TOLS_EXT

# Asteroids excluded: SPK coverage does not span -5000 to 5000
EXT_BODIES = ICRS_PLANETS + ECLIPTIC_BODIES + HYPOTHETICAL_BODIES


class TestExtLongitudeVelocity:
    """Longitude velocity (result[3]) for planets + ecliptic + hypothetical."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", EXT_BODIES)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Longitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        tol = ecl_tol if ecl_tol is not None else TOLS_EXT.SPEED_LON_DEG_DAY
        assert max_err < tol, (
            f"{body_name}: max lon speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestExtLatitudeVelocity:
    """Latitude velocity (result[4]) for planets + ecliptic + hypothetical."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", EXT_BODIES)
    def test_speed_latitude(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Latitude speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[4] - leb[4])
            if err > max_err:
                max_err = err
                worst_jd = jd

        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        tol = ecl_tol if ecl_tol is not None else TOLS_EXT.SPEED_LAT_DEG_DAY
        assert max_err < tol, (
            f"{body_name}: max lat speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )


class TestExtDistanceVelocity:
    """Distance velocity (result[5]) for ICRS planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_speed_distance(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Distance speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day at JD {worst_jd:.1f}"
        )


class TestExtVelocityEquatorial:
    """Velocity in equatorial coordinates for ICRS planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_equatorial_velocity(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        """Equatorial velocity matches Skyfield within tolerance."""
        flags = SEFLG_SPEED | SEFLG_EQUATORIAL
        max_err = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name}: max equatorial speed error = {max_err:.6f} deg/day"
        )
