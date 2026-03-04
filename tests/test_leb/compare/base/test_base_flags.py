"""
LEB vs Skyfield Comparison: Flag Combinations (Base Tier).

Validates all flag combinations supported natively by LEB
across the base tier range (1860-2140).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SEFLG_HELCTR,
    SEFLG_BARYCTR,
    SEFLG_TRUEPOS,
    SEFLG_NOABERR,
)

from tests.test_leb.compare.conftest import (
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_BASE

FLAG_BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]

FLAG_COMBINATIONS = [
    (SEFLG_SPEED, "default"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL, "equatorial"),
    (SEFLG_SPEED | SEFLG_J2000, "J2000"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_J2000, "equatorial_J2000"),
    (SEFLG_SPEED | SEFLG_HELCTR, "heliocentric"),
    (SEFLG_SPEED | SEFLG_BARYCTR, "barycentric"),
    (SEFLG_SPEED | SEFLG_TRUEPOS, "truepos"),
    (SEFLG_SPEED | SEFLG_NOABERR, "noaberr"),
]


class TestBaseFlagCombinations:
    """All 8 flag combinations for 4 planets."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize("flags,flag_name", FLAG_COMBINATIONS)
    def test_flags(
        self,
        compare: CompareHelper,
        base_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """All 6 components match within tolerance for each flag combo."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            err = max(lon_err, lat_err)

            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.EQUATORIAL_ARCSEC, (
            f'{body_name} {flag_name}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestBaseFlagVelocity:
    """Velocity under different flag combinations."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize(
        "flags,flag_name",
        [
            (SEFLG_SPEED | SEFLG_EQUATORIAL, "equatorial"),
            (SEFLG_SPEED | SEFLG_J2000, "J2000"),
            (SEFLG_SPEED | SEFLG_HELCTR, "heliocentric"),
            (SEFLG_SPEED | SEFLG_BARYCTR, "barycentric"),
        ],
    )
    def test_flag_speed(
        self,
        compare: CompareHelper,
        base_dates_50: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """Speed matches within tolerance for each flag combo."""
        max_err = 0.0

        for jd in base_dates_50:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name} {flag_name}: max speed error = {max_err:.6f} deg/day"
        )
