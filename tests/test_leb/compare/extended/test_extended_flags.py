"""
LEB vs Skyfield Comparison: Flag Combinations (Extended Tier).

Validates all flag combinations supported natively by LEB
across the extended tier range (-5000 to 5000).
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

from .conftest import TOLS_EXT

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


class TestExtFlagCombinations:
    """All 8 flag combinations for 4 planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize("flags,flag_name", FLAG_COMBINATIONS)
    def test_flags(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """All 6 components match within tolerance for each flag combo."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            err = max(lon_err, lat_err)

            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.EQUATORIAL_ARCSEC, (
            f'{body_name} {flag_name}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )
