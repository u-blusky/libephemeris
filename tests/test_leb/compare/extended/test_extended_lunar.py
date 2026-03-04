"""
LEB vs Skyfield Comparison: Lunar/Ecliptic Bodies (Extended Tier).

Validates all 6 Pipeline B (ecliptic-direct) bodies with per-body tolerances
across the extended tier range (-5000 to 5000).
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

from .conftest import TOLS_EXT


class TestExtLunarLongitude:
    """Lunar body longitude precision with per-body tolerance."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_longitude(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body longitude matches Skyfield within per-body tolerance."""
        tol_arcsec = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "lon", TOLS_EXT.ECLIPTIC_ARCSEC
        )
        # Extended tier may need slightly looser ecliptic tolerances
        tol_arcsec = max(tol_arcsec, TOLS_EXT.ECLIPTIC_ARCSEC)
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
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


class TestExtLunarLatitude:
    """Lunar body latitude precision."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_latitude(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body latitude matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[1] - leb[1]) * 3600.0
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.ECLIPTIC_ARCSEC, (
            f'{body_name}: max lat error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestExtLunarSpeed:
    """Lunar body velocity precision with per-body tolerance."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        ext_dates_500: list[float],
        body_id: int,
        body_name: str,
    ):
        """Lunar body longitude speed matches Skyfield within tolerance."""
        tol_speed = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "speed", TOLS_EXT.SPEED_LON_DEG_DAY
        )
        # Extended tier may need slightly looser speed tolerances
        tol_speed = max(tol_speed, TOLS_EXT.SPEED_LON_DEG_DAY)
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_500:
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
